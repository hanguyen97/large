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


//' Inner Loop using Autotune Lasso  
// [[Rcpp::export]]
List lasso_autotune(const arma::mat& X_X, const arma::colvec& X_Y, 
                     double sigma2, int n, double s_22, 
                     const arma::colvec& y, const arma::mat& Z, 
                     int node, int outer_iter, 
                     double alpha, const arma::vec& F_crit_values, 
                     double lambda0 = -1, 
                     bool verbose_i = false) {
   
   int p = X_X.n_rows;
   arma::colvec X_r_old = X_Y;
   arma::vec b_old = arma::zeros<arma::colvec>(p);
   double e_old = 1e6;
   bool F_test = true;
   double thresh = -1;
   
   if (lambda0 == -1) {
     lambda0 = max(abs(X_Y)) * (1 / s_22);
   }
   
   for (int iter = 0; iter <= 1000; iter++) {
     
     if (verbose_i) {
       Rcout << "node " << node << " inner iter " << iter << " sigma2 " << sigma2 << std::endl;
     }
     
     if (e_old > 0.0001) {
       vec b = b_old;
       vec X_r = X_r_old;
       vec sd_r = zeros<vec>(p);
       
       thresh = lambda0 * sigma2;
       for (int j = 0; j < p; j++) {
         X_r += X_X.col(j) * b[j];
         double X_rj = X_r[j];
         
         if (abs(X_rj) < 0.5 * thresh) {
           b[j] = 0;
         } else {
           b[j] = (abs(X_rj) - 0.5 * thresh) * ((X_rj > 0) ? 1 : -1) / X_X(j, j);
         }
         
         X_r -= X_X.col(j) * b[j];
         sd_r[j] = stddev(y - Z.cols(find(linspace<uvec>(0, p - 1, p) != j)) * b(find(linspace<uvec>(0, p - 1, p) != j)));
       }
       
       arma::uvec sorted_sd_idx = sort_index(sd_r, "descend");
       
       std::vector<int> support_set;
       if (F_test) {
         std::vector<int> sel_b;
         double sel_sigma2 = var(y);
         std::vector<int> new_b = sel_b;
         
         // Sequential F test for variable selection
         for (size_t j = 0; j < sorted_sd_idx.size(); j++) {
           int j_idx = sorted_sd_idx[j];
           new_b.push_back(j_idx);
           double new_sigma2 = get_LSsigma2(y, Z.cols(sorted_sd_idx.subvec(0, j)));
           double F_stat = (sel_sigma2 - new_sigma2) / (new_sigma2 / (n-(j+1)));
           
           if (F_stat > F_crit_values[j]) {
             sel_b = new_b;
             sel_sigma2 = new_sigma2;
           } else {
             if (verbose_i) {
               Rcout << "sel_b.size " << sel_b.size() << std::endl;
             }
             break;
           }
         }
         
         // Set sigma2 to sigma2 ols 
         sigma2 = sel_sigma2;
         
         // Check if support set converges
         if (support_set == sel_b) {
           F_test = false;
         } else {
           double e1 = 0.0;
           double e2 = 0.0;
           if (sel_b.size() > 0) {
             e1 = sum(abs(b(sorted_sd_idx.subvec(0, sel_b.size()-1)))); 
           }
           if (support_set.size() > 0) {
             e2 = sum(abs(b_old(sorted_sd_idx.subvec(0, support_set.size()-1))));
           }
           if (abs(e1 - e2) < 0.0001) {
             F_test = false;
           } else {
             support_set = sel_b;
           }
         }
       }
       
       b_old = b;
       X_r_old = X_r;
       double e = mean(square(X_r_old));
       
       if (abs(e - e_old) > 0.0001) {
         e_old = e;
       } else {
         break;
       }
     }
   }
   
   return List::create(Named("b")=b_old,
                       Named("sigma2")=sigma2,
                       Named("lambda")=thresh);
 }


//' Autotune Graphical Lasso 
//'
//' @param X Data matrix 
//' @param alpha alpha value of F test
//' @param thr Threshold for convergence. Default value is 1e-4.
//' @param maxit Maximum number of iterations of outer loop. Default 10,000
//' @return Estimated precision matrix
// [[Rcpp::export]]
List glasso_autotune(const arma::mat& X, double alpha = 0.1, double thr = 1e-4, 
                      int maxit = 1e4, bool verbose = true) {
   
   int n = X.n_rows;
   int p = X.n_cols;
   
   Rcpp::Function qf("qf");
   arma::vec F_crit_values = arma::linspace<arma::vec>(1, p-1, p-1);  // Replace loop
   F_crit_values.transform([&](double j) { 
     return Rcpp::as<double>(qf(1 - alpha, 1, n-(j+1)));
     });
   
   arma::mat S = (X.t() * X) / n;  
   arma::vec sigma2_hat = S.diag();
   
   arma::mat W_old = S;
   arma::mat W = S;
   arma::mat Theta = zeros<mat>(p, p);
   
   double e_old = 1e6;
   double e = 0;
   bool final_cycle = false;
   int niter = -1;
   
   double lambda_tmp = -1;
   arma::vec lambdav(p);
   lambdav.fill(-1);
   arma::vec b_hat(p-1, fill::zeros);
   
   for (int iter = 0; iter < maxit; iter++) {
     if (verbose) {
       Rcout << "glasso iter = " << iter+1 << "; error = " << round(e_old * 1000) / 1000 << std::endl;
     }
     
     for (int j = 0; j < p; j++) {
       arma::uvec idx = regspace<uvec>(0, p - 1);
       idx.shed_row(j);
       
       arma::mat W_11 = W(idx, idx);
       arma::colvec s_12 = S.submat(idx, uvec{(unsigned int)j});
       double s_22 = S(j, j);
       List fitted = lasso_autotune(W_11, s_12, 
                                    sigma2_hat(j), n, s_22, 
                                    X.col(j), X.cols(idx), 
                                    j, iter, alpha, F_crit_values,
                                    lambdav(j));
       
       b_hat = as<arma::vec>(fitted["b"]);
       sigma2_hat(j) = fitted["sigma2"];
       lambda_tmp = fitted["lambda"];
       lambdav(j) = lambda_tmp / sigma2_hat(j);
       //lambda(j) = fitted["lambda"]/fitted["sigma2"];
       
       arma::mat Wsub = W_11 * b_hat;
       W.submat(idx, uvec{(unsigned int)j}) = Wsub;
       W.submat(uvec{(unsigned int)j}, idx) = trans(Wsub);
       
       if (final_cycle) {
         Theta(j, j) = 1.0 / (W(j, j) - dot(W.submat(idx, uvec{(unsigned int)j}), b_hat));
         arma::mat Thetasub = -Theta(j, j) * b_hat;
         Theta.submat(idx, uvec{(unsigned int)j}) = Thetasub;
         Theta.submat(uvec{(unsigned int)j}, idx) = trans(Thetasub);
       }
       
     }
     
     e = norm(W - W_old, "fro");
     if (verbose) {
       Rcout << "change in W.err = " << std::abs(e - e_old) << std::endl;
     }
     
     if (final_cycle) { 
       niter = iter + 1;
       break;
     }
     if (std::abs(e - e_old) < thr) {
       final_cycle = true;
     } else {
       W_old = W;
       e_old = e;
     }
   }
   
   bool converged = (niter != -1);
   if (converged == 0) {
     Rcout << "did not converged; final glasso iter = " << maxit << std::endl;
   } else {
     Rcout << "final glasso iter = " << niter << std::endl;
   }
   return List::create(Named("Theta") = Theta, 
                       Named("sigma2.hat") = sigma2_hat,
                       Named("niter") = niter,
                       Named("converged") = converged);
 }
