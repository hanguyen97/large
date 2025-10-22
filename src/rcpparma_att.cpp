// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "stdlib.h"

using namespace Rcpp;
using namespace arma;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// Function to update support_ss if necessary
void updateSupport(std::vector<int>& support_ss, const std::vector<int>& sel_b) {
  
  std::unordered_set<int> support_set(support_ss.begin(), support_ss.end());
  bool all_found = true;
  for (int num : sel_b) {
    if (support_set.find(num) == support_set.end()) { // Element not found
      all_found = false;
      break;
    }
  }
  
  // If all elements are already present, return early
  if (all_found) {
    return;
  }
  
  // Otherwise, add missing elements to support_ss
  for (int num : sel_b) {
    if (support_set.find(num) == support_set.end()) {
      support_ss.push_back(num);
      support_set.insert(num);
    }
  }
}

bool haveSameElements(std::vector<int> sel_b, std::vector<int> support_ss) {
  return std::multiset<int>(sel_b.begin(), sel_b.end()) == std::multiset<int>(support_ss.begin(), support_ss.end());
}

// Function to compute LS estimator of sigma^2
double get_LSsigma2(const arma::colvec& y, const arma::mat& X) {
  int n = X.n_rows; 
  int p = X.n_cols; 
  
  arma::mat XtX_inv = arma::solve(X.t() * X, arma::eye(X.n_cols, X.n_cols));
  arma::mat H = X * XtX_inv * X.t();
  
  arma::mat I = eye(n, n);
  arma::mat M = I - H;
  
  double sigma2_hat = as_scalar(y.t() * M * y) / (n - p);
  
  return sigma2_hat;
}

std::vector<arma::uvec> get_sorted_indices(const arma::mat& R) {
  std::vector<arma::uvec> sorted_indices_per_col;
  
  for (arma::uword j = 0; j < R.n_cols; ++j) {
    arma::vec col = arma::abs(R.col(j));
    col.shed_row(j);  // Exclude diagonal
    arma::uvec indices = arma::sort_index(col, "descend");
    sorted_indices_per_col.push_back(indices);
  }
  
  return sorted_indices_per_col;
}

double avg_offd_abs(const arma::mat& W) {
  if (W.n_rows != W.n_cols) {
    throw std::invalid_argument("Matrix must be square.");
  }
  arma::mat abs_W = arma::abs(W);
  arma::mat diag_abs_W = arma::diagmat(abs_W);
  
  return arma::accu(abs_W - diag_abs_W) / (W.n_elem - W.n_rows);
}

//' Inner Loop using Autotune Lasso  
// [[Rcpp::export]]
List lasso_autotune(const arma::mat& X_X, const arma::colvec& X_Y, const arma::uvec& r_XY, 
                    const arma::vec& lambdas, 
                     double sigma2, int n, double s_22, 
                     const arma::colvec& y, const arma::mat& Z, 
                     int node, int outer_iter, 
                     double alpha, const arma::vec& F_crit_values, 
                     double lambda0 = -1, 
                     bool verbose_i = false,
                     bool penalize_diag = false) {
   
   int p = X_X.n_rows;
   int d = std::min(n-2, p);
   arma::colvec X_r_old = X_Y;
   arma::vec b_old = arma::zeros<arma::colvec>(p);
   arma::vec b = b_old;
   double e_old = 1e6;
   double sel_sigma2 = var(y);
   double sigma2_old = 1e6;
   bool F_test = true;
   bool sis = true;
   double thresh = -1;
   std::vector<int> support_ss;
   std::vector<int> support_ss_old;
   std::vector<int> sel_b;
   arma::vec lambdas_sub = arma::join_vert(lambdas.head(node), lambdas.tail(lambdas.n_elem-node-1));
   
   if (lambda0 == -1) {
     lambda0 = max(abs(X_Y)) * (1 / s_22);
   }
   
   for (int iter = 0; iter <= 1000; iter++) {
     
     b = b_old;
     vec X_r = X_r_old;
     vec sd_r = zeros<vec>(p);
     
     thresh = lambda0 * sigma2;
     for (int j = 0; j < p; j++) {
       X_r += X_X.col(j) * b[j];
       double X_rj = X_r[j];
       
       if (abs(X_rj) < 0.5 * thresh) {
         b[j] = 0;
       } else {
         if (penalize_diag) {
           b[j] = (abs(X_rj) - 0.5 * thresh) * ((X_rj > 0) ? 1 : -1) / (X_X(j, j) + lambdas_sub(j));
         } else {
           b[j] = (abs(X_rj) - 0.5 * thresh) * ((X_rj > 0) ? 1 : -1) / X_X(j, j);
         }
       }
       
       X_r -= X_X.col(j) * b[j];
       if (F_test) {
         arma::colvec pr = y - Z.cols(find(linspace<uvec>(0, p - 1, p) != j)) * b(find(linspace<uvec>(0, p - 1, p) != j));
         sd_r[j] = sqrt(as_scalar(pr.t() * pr) / n);
       }
     }
     
     
     if (F_test) {
       arma::uvec sorted_sd_idx;
       if (sis == true and iter == 0) {
         sorted_sd_idx = r_XY;
       } 
       
       std::vector<int> sel_b;
       std::vector<int> new_b = sel_b;
       sel_sigma2 = var(y);
       
       arma::uvec sel_idx;
       // Sequential F test for variable selection
       for (size_t j = 0; j < d; j++) {
         // Set sigma2 to sigma2 ols 
         sigma2 = sel_sigma2;
         
         // Find index of the maximum element
         auto max_it = std::max_element(sd_r.begin(), sd_r.end());
         int max_idx = std::distance(sd_r.begin(), max_it);
         // Set the maximum element to -1
         sd_r[max_idx] = -1;
         
         arma::uword j_idx;
         if (sis == true and iter == 0) {
           j_idx = sorted_sd_idx[j];
         } else {
           j_idx = max_idx;
         }
         new_b.push_back(static_cast<int>(j_idx));
         sel_idx = arma::join_vert(sel_idx, arma::uvec{j_idx});
         double new_sigma2 = get_LSsigma2(y, Z.cols(sel_idx));
         double F_stat = (sel_sigma2 - new_sigma2) / (new_sigma2 / (n-(j+1)));
         
         if (F_stat > F_crit_values[j]) {
           sel_b = new_b;
           sel_sigma2 = new_sigma2;
         } else {
           break;
         }
       }
       
       updateSupport(support_ss, sel_b);
       if (iter > 0) {
         // Check if support supper set converges
         if (haveSameElements(support_ss, support_ss_old)) {
           F_test = false;
         } 
         
         // Check if support supper set converges
         if (abs(sel_sigma2 - sigma2_old) < 0.0001) {
           F_test = false;
         }
       }
     }
     
     b_old = b;
     X_r_old = X_r;
     sigma2_old = sel_sigma2;
     support_ss_old = support_ss;
     double e = mean(square(X_r_old));
     
     if (abs(e - e_old) > 0.0001) {
       e_old = e;
     } else {
       break;
     }
   }
   
   return List::create(Named("b")=b,
                       Named("sigma2")=sigma2,
                       Named("thresh")=thresh);
 }


//' Locally Adaptive Regularization for Graph Estimation
//' 
//' Estimates a sparse inverse covariance matrix using a lasso (L1) penalty
//' with locally adaptive regularization
//' 
//' @param X A numeric data matrix
//' @param alpha Significance level of the F-test used to determine nodewise regularization. 
//'   Default is 0.02.
//' @param penalize_diag Logical; whether to penalize diagonal entries of the precision matrix. 
//'   Default is \code{FALSE}.
//' @param thr Convergence threshold. Default is 0.05.
//' @param maxit Maximum number of iterations for the outer loop. Default is 20.
//' @param verbose Logical; if \code{TRUE}, print overall iteration progress. Default is \code{TRUE}.
//' 
//' @return
//' \itemize{
//'   \item \code{Theta}: Estimated precision matrix.
//'   \item \code{sigma2.hat}: Estimated residual variances from nodewise regressions.
//'   \item \code{lambdas}: Vector of adaptively selected regularization parameters for each node.
//'   \item \code{niter}: Number of outer iterations performed.
//'   \item \code{converged}: Logical indicating whether the algorithm converged.
//' }
//' 
// [[Rcpp::export]]
List fit_large(const arma::mat& X, double alpha = 0.02, 
                      double penalize_diag = false,
                      double thr = 0.05, int maxit = 20, 
                      bool verbose = true) {
   
   int n = X.n_rows;
   int p = X.n_cols;
   bool verbose_i = false;
     
   Rcpp::Function qf("qf");
   // Prevent non-positive df2 of F-test
   int d = std::min(n-2, p-1);
   arma::vec F_crit_values = arma::linspace<arma::vec>(1, d, d);  
   F_crit_values.transform([&](double j) { 
     return Rcpp::as<double>(qf(1 - alpha, 1, n-(j+1)));
   });
   
   arma::mat S = (X.t() * X) / n;  
   arma::mat R = cor(X);
   arma::vec sigma2_hat = S.diag();
   
   arma::mat W = S;
   arma::mat Theta = zeros<mat>(p, p);
   
   double e_old = 1e6;
   double e = 0;
   bool final_cycle = false;
   bool valid_diag = true;
   int niter = -1;
   
   double thresh = -1;
   arma::vec lambdas(p);
   lambdas.fill(-1);
   arma::vec lambda0s(p);
   lambda0s.fill(-1);
   arma::vec b_hat(p-1, fill::zeros);
   
   arma::mat W_old = W;
   
   if (penalize_diag) {
     for (int j = 0; j < p; j++) {
       arma::uvec idx = regspace<uvec>(0, p - 1);
       idx.shed_row(j);
       arma::colvec s_12 = S.submat(idx, uvec{(unsigned int)j});
       lambdas(j) = 0.5 * max(abs(s_12)) ;
     }
   }
   
   // Sort indices of cor matrix
   std::vector<arma::uvec> R_sorted_idx = get_sorted_indices(R);
   
   for (int iter = 0; iter < maxit; iter++) {
     if (verbose) {
       Rcout << "iter = " << iter+1 << "; relative F-norm change = " << round(e_old * 1000) / 1000 << std::endl;
     }
     
     for (int j = 0; j < p; j++) {
       arma::uvec idx = regspace<uvec>(0, p - 1);
       idx.shed_row(j);
       
       arma::mat W_11 = W(idx, idx);
       arma::colvec s_12 = S.submat(idx, uvec{(unsigned int)j});
       double s_22 = S(j, j);
       
       List fitted = lasso_autotune(W_11, s_12, R_sorted_idx[j],
                                    lambdas, sigma2_hat(j), 
                                    n, s_22, 
                                    X.col(j), X.cols(idx), 
                                    j, iter, alpha, F_crit_values,
                                    lambda0s(j), 
                                    verbose_i, penalize_diag);
       
       b_hat = as<arma::vec>(fitted["b"]);
       sigma2_hat(j) = fitted["sigma2"];
       thresh = fitted["thresh"];
       lambda0s(j) = thresh / sigma2_hat(j);
       lambdas(j) = 0.5 * thresh;
       
       arma::vec lambdas_sub = arma::join_vert(lambdas.head(j), lambdas.tail(lambdas.n_elem-j-1));
       arma::mat Wsub = (W_11 + arma::diagmat(lambdas_sub)) * b_hat;
       W.submat(idx, uvec{(unsigned int)j}) = Wsub;
       W.submat(uvec{(unsigned int)j}, idx) = trans(Wsub);
       
       if (final_cycle) {
         if (penalize_diag) {
           Theta(j, j) = 1.0 / (sigma2_hat(j)+lambdas(j));
         } else {
           Theta(j, j) = 1.0 / sigma2_hat(j);
         }
         arma::mat Thetasub = -Theta(j, j) * b_hat;
         Theta.submat(idx, uvec{(unsigned int)j}) = Thetasub;
         Theta.submat(uvec{(unsigned int)j}, idx) = trans(Thetasub);
         
         if (Theta(j, j) < 0) {
           Rcout << "Diagonal of node " << j + 1 << " is negative. Consider a smaller alpha threshold!" << std::endl;
           valid_diag = false;
         }
         
       }
       
     }
     
     arma::mat W_diff = W - W_old;
     e = norm(W_diff, "fro") / norm(W_old, "fro");

     if (final_cycle) { 
       if (iter < maxit-1) {
         niter = iter + 1;
       }
       break;
     }
     
     if (e < thr) {
       final_cycle = true;
     } else {
       W_old = W;
       e_old = e;
     }
     if (iter == maxit-2) {
       final_cycle = true;
     }
   }
   
   bool converged = (niter != -1);
   if (converged == 0) {
     Rcout << "did not converged; final iter = " << maxit << std::endl;
   } else {
     Rcout << "final iter = " << niter << std::endl;
   }
   if (!valid_diag) {
     converged = 0;
   }
   return List::create(Named("Theta") = Theta, 
                       Named("sigma2.hat") = sigma2_hat,
                       Named("lambdas") = lambdas,
                       Named("niter") = niter,
                       Named("converged") = converged);
 }
