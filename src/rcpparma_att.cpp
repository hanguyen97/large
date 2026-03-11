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

bool is_subset(std::vector<int> a, std::vector<int> b) {
  // sort both vectors
  std::sort(a.begin(), a.end());
  std::sort(b.begin(), b.end());
  
  // check if 'a' is a subset of 'b'
  return std::includes(b.begin(), b.end(), a.begin(), a.end());
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
   std::vector<int> support_set;
   std::vector<int> support_set_old;
   // std::vector<int> sel_b;
   arma::vec lambdas_sub = arma::join_vert(lambdas.head(node), lambdas.tail(lambdas.n_elem-node-1));
   
   if (lambda0 == -1) {
     lambda0 = 0.5 * max(abs(X_Y)) * (1 / s_22);
   }
   
   for (int iter = 0; iter <= 1000; iter++) {
     
     b = b_old;
     vec X_r = X_r_old;
     vec sd_r = zeros<vec>(p);
     
     thresh = lambda0 * sigma2;
     for (int j = 0; j < p; j++) {
       X_r += X_X.col(j) * b[j];
       double X_rj = X_r[j];
       
       if (abs(X_rj) < thresh) {
         b[j] = 0;
       } else {
         if (penalize_diag) {
           b[j] = (abs(X_rj) - thresh) * ((X_rj > 0) ? 1 : -1) / (X_X(j, j) + lambdas_sub(j));
         } else {
           b[j] = (abs(X_rj) - thresh) * ((X_rj > 0) ? 1 : -1) / X_X(j, j);
         }
       }
       
       X_r -= X_X.col(j) * b[j];
     }
     
     
     if (F_test) {
       // ---------------------------------------- // 
       // ------------- Sigma update ------------- //
       // ---------------------------------------- // 
       
       // get SD of partial residual and sort 
       for (int j = 0; j < p; j++) {
         arma::colvec pr = y - Z.cols(find(linspace<uvec>(0, p - 1, p) != j)) * b(find(linspace<uvec>(0, p - 1, p) != j));
         sd_r[j] = sqrt(as_scalar(pr.t() * pr) / n);
       }
       arma::uvec pred_ranking;
       if (sis == true and iter == 0) {
         pred_ranking = r_XY;
       } else {
         pred_ranking = arma::sort_index(sd_r, "descend");
       }
       
       // initialize sigma update variables
       // std::vector<int> sel_b;
       // std::vector<int> new_b = sel_b;
       // sel_sigma2 = var(y);
       // arma::uvec sel_idx;
       arma::colvec y_temp = y;
       arma::mat u = Z.cols(pred_ranking);
       support_set = std::vector<int>{};
       
       // sequential F test for variable selection
       for (size_t j = 0; j < d; j++) {
         // set sigma2 to sigma2 ols 
         sigma2 = sel_sigma2;
         
         // Gram-Schmidt process
         arma::vec uj = Z.col(pred_ranking[j]); 
         if (j > 0) {
           for (arma::uword k = 0; k < j; ++k) {
             arma::vec uk = u.col(k);
             double num = arma::dot(uj, uk);        
             double den = arma::dot(uk, uk);      
             if (den > 1e-12) {
               uj = u.col(j) - (num / den) * uk;
             }
           }
         }
         u.col(j) = uj;
         double den = arma::dot(uj, uj);
         arma::vec y_hat = arma::zeros<arma::vec>(uj.n_elem);
         if (den > 1e-12) {
           y_hat = (arma::dot(y_temp, uj) / den) * uj;
         }
         
         // calculate F-statistics
         double new_sigma2 = as_scalar(arma::dot(y_temp - y_hat, y_temp - y_hat)) / (y.n_elem - (j+1)) ;
         double new_RSS = as_scalar(arma::dot(y_temp - y_hat, y_temp - y_hat));
         double old_RSS = as_scalar(arma::dot(y_temp, y_temp));
         double F_stat = (old_RSS - new_RSS) / (new_RSS / (n-(j+1)));
         
         if (F_stat > F_crit_values[j]) {
           y_temp = y_temp - y_hat;
           sel_sigma2 = new_sigma2;
           support_set.push_back(pred_ranking[j]);
         } else {
           break;
         }
       }
       
       // ----------------------------------------------- // 
       // ------------- End of Sigma update ------------- //
       // ----------------------------------------------- // 
       
       if (iter > 0) {
         
         // // check if support supper set converges
         if (is_subset(support_set, support_set_old)) {
           // cout << "support set converges" << endl;
           //   for (int x : support_set) {
           //     std::cout << x << " ";
           //   }
           //   std::cout << std::endl; // New line after printing
           
           F_test = false;
         }
         
         // check if sigma2 converges
         if (abs(sel_sigma2 - sigma2_old) < 1e-6) {
           // cout << "sigma2 converges" << endl;
           F_test = false;
         }
       }
     }
     
     b_old = b;
     X_r_old = X_r;
     sigma2_old = sel_sigma2;
     support_set_old = support_set;
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
                bool penalize_diag = false,
                bool penalize_test = false,
                double thr = 0.05, int maxit = 20, 
                bool verbose = true) {
   
   std::vector<double> runtimes;
   auto start = std::chrono::high_resolution_clock::now();
   
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
   
   if (penalize_diag) {
     for (int j = 0; j < p; j++) {
       arma::uvec idx = regspace<uvec>(0, p - 1);
       idx.shed_row(j);
       arma::colvec s_12 = S.submat(idx, uvec{(unsigned int)j});
       lambdas(j) = 0.5 * max(abs(s_12));
     }
   }
   
   if (penalize_test) {
     // for (int j = 0; j < p; j++) {
     //   arma::uvec idx = regspace<uvec>(0, p - 1);
     //   idx.shed_row(j);
     //   arma::colvec s_12 = S.submat(idx, uvec{(unsigned int)j});
     //   lambdas(j) = 0.5 * max(abs(s_12));
     // }
     // W.diag() += lambdas;
     W.diag() += mean(S.diag());
   }
   
   arma::mat W_old = W;
   
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
       lambdas(j) = thresh;
       
       arma::vec lambdas_sub = arma::join_vert(lambdas.head(j), lambdas.tail(lambdas.n_elem-j-1));
       arma::mat Wsub;
       if (penalize_diag) {
         Wsub = (W_11 + arma::diagmat(lambdas_sub)) * b_hat;
       } else {
         Wsub = W_11 * b_hat;
       }
       W.submat(idx, uvec{(unsigned int)j}) = Wsub;
       W.submat(uvec{(unsigned int)j}, idx) = trans(Wsub);
       
       if (final_cycle) {
         if (penalize_diag) {
           Theta(j, j) = 1.0 / sigma2_hat(j);
           // Theta(j, j) = 1.0 / (sigma2_hat(j)+lambdas(j));
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
       auto end = std::chrono::high_resolution_clock::now();
       std::chrono::duration<double> elapsed = end - start;
       runtimes.push_back(elapsed.count());
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
     
     auto end = std::chrono::high_resolution_clock::now();
     std::chrono::duration<double> elapsed = end - start;
     start = end;
     runtimes.push_back(elapsed.count());
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
                       Named("converged") = converged,
                       Named("runtimes") = runtimes);
 }