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

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::export]]
arma::mat rcpparma_hello_world() {
  arma::mat m1 = arma::eye<arma::mat>(3, 3);
  arma::mat m2 = arma::eye<arma::mat>(3, 3);
  
  return m1 + 3 * (m1 + m2);
}


// another simple example: outer product of a vector, 
// returning a matrix
//
// [[Rcpp::export]]
arma::mat rcpparma_outerproduct(const arma::colvec & x) {
  arma::mat m = x * x.t();
  return m;
}

// and the inner product returns a scalar
//
// [[Rcpp::export]]
double rcpparma_innerproduct(const arma::colvec & x) {
  double v = arma::as_scalar(x.t() * x);
  return v;
}


// and we can use Rcpp::List to return both at the same time
//
// [[Rcpp::export]]
Rcpp::List rcpparma_bothproducts(const arma::colvec & x) {
  arma::mat op = x * x.t();
  double    ip = arma::as_scalar(x.t() * x);
  return Rcpp::List::create(Rcpp::Named("outer")=op,
                            Rcpp::Named("inner")=ip);
}

// Function to compute LS estimator of sigma^2
double get_LSsigma2(arma::colvec y, arma::mat X) {
  int n = X.n_rows; 
  int p = X.n_cols; 
  
  arma::mat XtX_inv = inv(X.t() * X);
  arma::mat H = X * XtX_inv * X.t();
  
  arma::mat I = eye(n, n);
  arma::mat M = I - H;
  
  double sigma2_hat = as_scalar(y.t() * M * y) / (n - p);
  
  return sigma2_hat;
}


//' Get an initial set of putative variables for the GAM algorithm
//'
//' @param X_X matrix of putative variables
//' @param X_Y the (possibly transformed) response
//' @param sigma2 a threshold of (absolute) correlation above which a pair is considered highly correlated
//' @return a list containing  variables to ignore because they are highly correlated with other, and SLR coefficients
// [[Rcpp::export]]
List lasso_autotune(arma::mat X_X, arma::colvec X_Y, 
                    double sigma2, int n, double s_22, 
                    arma::colvec y, arma::mat Z, 
                    int node, int outer_iter, 
                    bool verbose,
                    double lambda0 = -1, double alpha = 0.1) {
  
  int p = X_X.n_rows;
  arma::colvec X_r_old = X_Y;
  vec b_old = arma::zeros<arma::colvec>(p);
  double e_old = 1e6;
  bool F_test = true;
  double thresh = -1;

  if (lambda0 == -1) {
    lambda0 = max(abs(X_Y)) * (1 / s_22);
  }
  
  for (int iter = 0; iter <= 1000; iter++) {
    
    if (verbose) {
      Rcout << "inner iter " << iter << std::endl;
    }
    
    if (e_old > 0.0001) {
      vec b = b_old;
      vec X_r = X_r_old;
      vec sd_r = zeros<vec>(p);
      
      for (int j = 0; j < p; j++) {
        X_r += X_X.col(j) * b[j];
        double X_rj = X_r[j];
        thresh = lambda0 * sigma2;
        
        if (abs(X_rj) < 0.5 * thresh) {
          b[j] = 0;
        } else {
          b[j] = (abs(X_rj) - 0.5 * thresh) * ((X_rj > 0) ? 1 : -1) / X_X(j, j);
        }
        
        X_r -= X_X.col(j) * b[j];
        sd_r[j] = stddev(y - Z.cols(find(linspace<uvec>(0, p - 1, p) != j)) * b(find(linspace<uvec>(0, p - 1, p) != j)));
      }
      
      arma::uvec sorted_sd_idx = sort_index(sd_r, "descend");
      Rcout << "sorted_sd_idx" << sorted_sd_idx << std::endl;
      
      std::vector<int> support_set;
      if (F_test) {
        std::vector<int> sel_b;
        double sel_sigma2 = var(y);
        std::vector<int> new_b = sel_b;
        
        for (size_t j = 0; j < sorted_sd_idx.size(); j++) {
          int j_idx = sorted_sd_idx[j];
          new_b.push_back(j_idx);
          double new_sigma2 = get_LSsigma2(y, Z.cols(sorted_sd_idx.subvec(0, j)));
          double F_stat = (sel_sigma2 - new_sigma2) / (new_sigma2 / (j+1));
            
          Rcout << "new_sigma2" << new_sigma2 << std::endl;
          
          Rcpp::Function qf("qf");
          double F_crit = Rcpp::as<double>(qf(1-alpha, 1, j+1));
            
          if (F_stat < F_crit) {
            Rcout << "F continue" << std::endl;
            sel_b = new_b;
            sel_sigma2 = new_sigma2;
          } else {
            Rcout << "F ends" << std::endl;
            break;
          }
          
          if (support_set == sel_b) {
            F_test = false;
          } else if (j > 1) {
            double e1 = sum(abs(b(sorted_sd_idx.subvec(0, j)))); 
            double e2 = sum(abs(b_old(sorted_sd_idx.subvec(0, j-1))));
            if (abs(e1 - e2) < 0.0001) {
              F_test = false;
            } else {
              support_set = sel_b;
            }
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
  
  Rcout << "This is a fixed text message." << F_test << std::endl;
  
  return List::create(Named("b")=b_old,
                      Named("sigma2")=sigma2,
                      Named("lambda")=thresh);
}

