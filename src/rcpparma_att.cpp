// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <vector>
#include <algorithm>
#include <chrono>
#include <cmath>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

bool is_subset(std::vector<int> a, std::vector<int> b) {
  std::sort(a.begin(), a.end());
  std::sort(b.begin(), b.end());
  return std::includes(b.begin(), b.end(), a.begin(), a.end());
}

std::vector<arma::uvec> get_sorted_indices(const arma::mat& R) {
  std::vector<arma::uvec> sorted_indices_per_col;
  for (arma::uword j = 0; j < R.n_cols; ++j) {
    arma::vec col = arma::abs(R.col(j));
    col.shed_row(j);  // Exclude diagonal
    sorted_indices_per_col.push_back(arma::sort_index(col, "descend"));
  }
  return sorted_indices_per_col;
}

//' Inner Loop using Autotune Lasso
 // [[Rcpp::export]]
 List lasso_autotune(const arma::mat& X_X, const arma::colvec& X_Y, const arma::uvec& r_XY,
                     const arma::vec& lambdas,
                     double sigma2, int n, double s_22,
                     const arma::colvec& y, const arma::mat& Z,
                     int node, double alpha, const arma::vec& F_crit_values,
                     double lambda0 = -1,
                     bool penalize_diag = false,
                     bool update_sigma = true,     // false = skip F-test, keep sigma2 fixed (used when frozen)
                     double inner_tol = -1) {      // > 0: stop on relative beta change
   
   int p = X_X.n_rows;
   int d = std::min(n - 2, p);
   arma::colvec X_r_old = X_Y;
   arma::vec b_old = arma::zeros<arma::colvec>(p);
   arma::vec b = b_old;
   double e_old = 1e6;
   double sel_sigma2 = var(y);
   double sigma2_old = 1e6;
   bool F_test = update_sigma;
   double thresh = -1;
   std::vector<int> support_set;
   std::vector<int> support_set_old;
   arma::vec lambdas_sub = arma::join_vert(lambdas.head(node), lambdas.tail(lambdas.n_elem - node - 1));
   
   if (lambda0 == -1) {
     lambda0 = 0.5 * max(abs(X_Y)) * (1 / s_22);
   }
   
   for (int iter = 0; iter <= 1000; iter++) {
     
     b = b_old;
     vec X_r = X_r_old;
     
     // ------------- Coordinate descent sweep ------------- //
     thresh = lambda0 * sigma2;
     for (int j = 0; j < p; j++) {
       X_r += X_X.col(j) * b[j];
       double X_rj = X_r[j];
       
       if (std::abs(X_rj) < thresh) {
         b[j] = 0;
       } else {
         double denom = penalize_diag ? (X_X(j, j) + lambdas_sub(j)) : X_X(j, j);
         b[j] = (std::abs(X_rj) - thresh) * ((X_rj > 0) ? 1 : -1) / denom;
       }
       
       X_r -= X_X.col(j) * b[j];
     }
     
     // ------------- Sigma update (sequential F-test) ------------- //
     if (F_test) {
       sel_sigma2 = arma::dot(y, y) / y.n_elem;
       
       // rank predictors: SIS ranking on the first pass, otherwise by SD of partial residual
       arma::uvec pred_ranking;
       if (iter == 0) {
         pred_ranking = r_XY;
       } else {
         // partial residual for j = full residual + Z_j * b_j  (avoids rebuilding Z without column j)
         arma::colvec r_full = y - Z * b;
         arma::vec sd_r(p);
         for (int j = 0; j < p; j++) {
           arma::colvec pr = r_full + Z.col(j) * b(j);
           sd_r[j] = std::sqrt(arma::dot(pr, pr) / n);
         }
         pred_ranking = arma::sort_index(sd_r, "descend");
       }
       
       arma::colvec y_temp = y;
       arma::mat u(y.n_elem, std::max(d, 1));   // orthogonalized columns, filled as we go
       support_set.clear();
       
       for (int j = 0; j < d; j++) {
         sigma2 = sel_sigma2;
         
         // Gram-Schmidt
         arma::vec uj = Z.col(pred_ranking[j]);
         for (int k = 0; k < j; ++k) {
           double den_k = arma::dot(u.col(k), u.col(k));
           if (den_k > 1e-12) uj -= (arma::dot(uj, u.col(k)) / den_k) * u.col(k);
         }
         u.col(j) = uj;
         double den = arma::dot(uj, uj);
         arma::vec y_hat = arma::zeros<arma::vec>(uj.n_elem);
         if (den > 1e-12) {
           y_hat = (arma::dot(y_temp, uj) / den) * uj;
         }
         
         // F-statistic
         arma::vec res = y_temp - y_hat;
         double new_RSS = arma::dot(res, res);
         double old_RSS = arma::dot(y_temp, y_temp);
         double new_sigma2 = new_RSS / (y.n_elem - (j + 1));
         double F_stat = (old_RSS - new_RSS) / (new_RSS / (n - (j + 1)));
         
         if (F_stat > F_crit_values[j]) {
           y_temp = res;
           sel_sigma2 = new_sigma2;
           support_set.push_back(pred_ranking[j]);
         } else {
           break;
         }
       }
       
       if (iter > 0) {
         if (is_subset(support_set, support_set_old)) F_test = false;   // support converged
         if (std::abs(sel_sigma2 - sigma2_old) < 1e-6) F_test = false;  // sigma2 converged
       }
     }
     
     // relative beta change (only needed for the inner_tol stopping rule)
     double db1 = 0;
     if (inner_tol > 0) {
       db1 = arma::norm(b - b_old, 2) / std::max(arma::norm(b_old, 2), 1e-12);
     }
     
     b_old = b;
     X_r_old = X_r;
     sigma2_old = sel_sigma2;
     support_set_old = support_set;
     double e = mean(square(X_r_old));
     
     bool stop;
     if (inner_tol > 0) {
       stop = (iter > 0) && !F_test && (db1 < inner_tol);
     } else {
       stop = !(std::abs(e - e_old) > 0.0001);
     }
     if (stop) break;
     e_old = e;
   }
   
   return List::create(Named("b") = b,
                       Named("sigma2") = sigma2,
                       Named("thresh") = thresh,
                       Named("lambda0") = lambda0);
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
 //' @param freeze_after Freeze the nodewise variances (and hence the penalties) after this many
 //'   outer iterations (must be >= 2), at the average of the last two iterations.
 //'   Default \code{-1} (never).
 //' @param auto_freeze Logical; freeze automatically when a 2-cycle is detected
 //'   (\code{e(t,t-2) < 0.2 * e(t,t-1)} while not converged). Default \code{FALSE}.
 //' @param inner_tol If positive, stop each inner loop when the relative change in the
 //'   regression coefficients falls below this value. Default \code{-1} (original rule).
 //'
 //' @return
 //' \itemize{
 //'   \item \code{Theta}: Estimated precision matrix.
 //'   \item \code{sigma2.hat}: Estimated residual variances from nodewise regressions.
 //'   \item \code{lambdas}: Vector of adaptively selected regularization parameters for each node.
 //'   \item \code{niter}: Number of outer iterations performed.
 //'   \item \code{converged}: Logical indicating whether the convergence threshold was reached.
 //'   \item \code{frozen_at}: Outer iteration from which the variances were frozen (NA if never).
 //'   \item \code{runtimes}: Runtime of each outer iteration (seconds).
 //'   \item \code{e_trace}: Relative F-norm change of W vs. the previous outer iteration.
 //'   \item \code{e2_trace}: Relative F-norm change of W vs. two outer iterations back
 //'     (only computed when \code{auto_freeze = TRUE}, otherwise NA).
 //' }
 //'
 // [[Rcpp::export]]
 List fit_large(const arma::mat& X, double alpha = 0.02,
                bool penalize_diag = false,
                bool penalize_test = false,
                double thr = 0.05, int maxit = 20,
                bool verbose = true,
                int freeze_after = -1,
                bool auto_freeze = false,
                double inner_tol = -1) {
   
   std::vector<double> runtimes;
   auto start = std::chrono::high_resolution_clock::now();
   
   int n = X.n_rows;
   int p = X.n_cols;
   
   Rcpp::Function qf("qf");
   // Prevent non-positive df2 of F-test
   int d = std::min(n - 2, p - 1);
   arma::vec F_crit_values = arma::linspace<arma::vec>(1, d, d);
   F_crit_values.transform([&](double j) {
     return Rcpp::as<double>(qf(1 - alpha, 1, n - (j + 1)));
   });
   
   arma::mat S = (X.t() * X) / n;
   arma::mat R = cor(X);
   arma::vec sigma2_hat = S.diag();
   
   arma::mat W = S;
   arma::mat Theta = zeros<mat>(p, p);
   
   double e = 0;
   bool final_cycle = false;
   bool valid_diag = true;
   bool hit_thr = false;
   int niter = -1;
   
   double thresh = -1;
   arma::vec lambdas(p);
   lambdas.fill(-1);
   arma::vec lambda0s(p);
   lambda0s.fill(-1);
   arma::vec b_hat(p - 1, fill::zeros);
   
   // freezing state
   bool frozen = false;
   bool freeze_pending = false;
   int frozen_at = NA_INTEGER;
   arma::vec sigma2_prev = sigma2_hat;   // sigma2 from the sweep before the most recent one
   
   if (penalize_diag) {
     for (int j = 0; j < p; j++) {
       arma::uvec idx = regspace<uvec>(0, p - 1);
       idx.shed_row(j);
       arma::colvec s_12 = S.submat(idx, uvec{(unsigned int)j});
       lambdas(j) = 0.5 * max(abs(s_12));
     }
   }
   
   if (penalize_test) {
     for (int j = 0; j < p; j++) {
       arma::uvec idx = regspace<uvec>(0, p - 1);
       idx.shed_row(j);
       arma::colvec s_12 = S.submat(idx, uvec{(unsigned int)j});
       lambdas(j) = 0.5 * mean(abs(s_12));
     }
     W.diag() += lambdas;
   }
   
   arma::mat W_old  = W;   // W one outer iteration back
   arma::mat W_old2;       // W two outer iterations back (only kept for auto_freeze)
   if (auto_freeze) W_old2 = W;
   std::vector<double> e_trace, e2_trace;
   
   // Sort indices of cor matrix
   std::vector<arma::uvec> R_sorted_idx = get_sorted_indices(R);
   
   for (int iter = 0; iter < maxit; iter++) {
     
     // apply a pending freeze: average sigma2 over the last two sweeps (phase-independent)
     if (freeze_pending && !frozen) {
       sigma2_hat = 0.5 * (sigma2_hat + sigma2_prev);
       lambdas = lambda0s % sigma2_hat;
       frozen = true;
       frozen_at = iter + 1;
       if (verbose) {
         Rcout << "  sigma2 frozen from sweep " << iter + 1 << " (average of last two sweeps)" << std::endl;
       }
     }
     sigma2_prev = sigma2_hat;
     
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
                                    j, alpha, F_crit_values,
                                    lambda0s(j),
                                    penalize_diag,
                                    !frozen, inner_tol);
       
       b_hat = as<arma::vec>(fitted["b"]);
       sigma2_hat(j) = fitted["sigma2"];
       thresh = fitted["thresh"];
       lambda0s(j) = thresh / sigma2_hat(j);         
       lambdas(j) = thresh;
       
       arma::vec lambdas_sub = arma::join_vert(lambdas.head(j), lambdas.tail(lambdas.n_elem - j - 1));
       arma::mat Wsub;
       if (penalize_diag) {
         Wsub = (W_11 + arma::diagmat(lambdas_sub)) * b_hat;
       } else {
         Wsub = W_11 * b_hat;
       }
       W.submat(idx, uvec{(unsigned int)j}) = Wsub;
       W.submat(uvec{(unsigned int)j}, idx) = trans(Wsub);
       
       if (final_cycle) {
         Theta(j, j) = 1.0 / sigma2_hat(j);
         arma::mat Thetasub = -Theta(j, j) * b_hat;
         Theta.submat(idx, uvec{(unsigned int)j}) = Thetasub;
         Theta.submat(uvec{(unsigned int)j}, idx) = trans(Thetasub);
         
         if (Theta(j, j) < 0) {
           Rcout << "Diagonal of node " << j + 1 << " is negative. Consider a smaller alpha threshold!" << std::endl;
           valid_diag = false;
         }
       }
     }
     
     // change vs. one sweep back, and (for 2-cycle detection) vs. two sweeps back
     e = norm(W - W_old, "fro") / norm(W_old, "fro");
     double e2 = NA_REAL;
     if (auto_freeze && iter > 0) {
       e2 = norm(W - W_old2, "fro") / norm(W_old2, "fro");
     }
     e_trace.push_back(e);
     e2_trace.push_back(e2);
     
     if (verbose) {
       Rcout << "sweep " << iter + 1 << ": e(t,t-1) = " << e;
       if (auto_freeze && iter > 0) Rcout << ", e(t,t-2) = " << e2;
       Rcout << std::endl;
     }
     
     if (auto_freeze) W_old2 = W_old;
     W_old = W;
     
     if (final_cycle) {
       niter = iter + 1;
       auto end = std::chrono::high_resolution_clock::now();
       std::chrono::duration<double> elapsed = end - start;
       runtimes.push_back(elapsed.count());
       break;
     }
     
     if (e < thr) {
       final_cycle = true;
       hit_thr = true;
     } else if (!frozen && !freeze_pending) {
       // decide whether to freeze sigma2 from the next sweep on
       if (freeze_after >= 2 && iter + 1 >= freeze_after) {
         freeze_pending = true;
       }
       if (auto_freeze && iter >= 2 && e2 < 0.2 * e) {
         freeze_pending = true;
         if (verbose) {
           Rcout << "  2-cycle detected at sweep " << iter + 1
                 << " (e(t,t-2) / e(t,t-1) = " << e2 / e << ")" << std::endl;
         }
       }
     }
     if (iter == maxit - 2) {
       final_cycle = true;
     }
     
     auto end = std::chrono::high_resolution_clock::now();
     std::chrono::duration<double> elapsed = end - start;
     start = end;
     runtimes.push_back(elapsed.count());
   }
   
   bool converged = hit_thr && valid_diag;
   if (!hit_thr) {
     Rcout << "did not converge; final iter = " << niter << std::endl;
   } else {
     Rcout << "converged; final iter = " << niter << std::endl;
   }
   
   return List::create(Named("Theta") = Theta,
                       Named("sigma2.hat") = sigma2_hat,
                       Named("lambdas") = lambdas,
                       Named("niter") = niter,
                       Named("converged") = converged,
                       Named("frozen_at") = frozen_at,
                       Named("runtimes") = runtimes,
                       Named("e_trace") = wrap(e_trace),
                       Named("e2_trace") = wrap(e2_trace));
 }