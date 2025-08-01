AutotuneGLASSO: An Automatic Approach to Variable-Specific Tuning for
Gaussian Graphical Models
================

## Introduction

We provide the `ATTglasso` package for automatic tuning of
regularization parameters in graphical Lasso (GLASSO), enhancing both
estimation accuracy and graph recovery by leveraging penalties. GLASSO
estimates the precision matrix $\Theta$ of a Gaussian graphical model
(GGM) by maximizing the $\ell_1$-penalized log-likelihood over the space
of positive semi-definite matrices: $$
    \hat{\Theta} \in \underset{\Theta \succeq 0}{\arg\max} \left\{ \log \det(\Theta) - \mathrm{trace}(S\Theta) - \lambda \|\Theta\|_1 \right\},
$$ where $S = \frac{1}{n} \sum_{i=1}^n x_i x_i^\top$ is the sample
covariance matrix and $\|\Theta\|_1$ denotes the elementwise $\ell_1$
norm. The tuning parameter $\lambda \geq 0$ controls the sparsity of the
estimate.

Unlike standard GLASSO, which relies on a single global penalty,
`AutotuneGLASSO` adaptively learns a set of node-specific penalties
$\lambda_j$. It does so by augmenting the nodewise Lasso regression step
to jointly estimate both regression coefficients and error variances,
allowing more flexible and data-driven regularization across nodes.

## Installation

You can download the `ATTglasso` package from Github.

``` r
# You can install the development version from GitHub:
# install.packages("devtools")
devtools::install_github("hanguyen97/ATTglasso")
```

## Quick Start: ATTglasso

We create a data set for illustration:

``` r
set.seed(1)

# Generate data from AR(1) model
p <- 10
n <- 200

Theta <- matrix(data=0, nrow=p, ncol=p)
diag(Theta) <- 1
offd1 <- 0.3
diag(Theta[1:(p-1), 2:(p)]) <- offd1
diag(Theta[2:(p), 1:(p-1)]) <- offd1

Sigma <- solve(Theta)
diag(Sigma)
```

    ##  [1] 1.111111 1.234568 1.248285 1.249809 1.249976 1.249976 1.249809 1.248285
    ##  [9] 1.234568 1.111111

``` r
library(MASS)
X <- mvrnorm(n=n, mu=rep(0,p), Sigma=Sigma)
```

We can estimate the precision matrix using `glasso_autotune`.

``` r
library(ATTglasso)
start.T <- Sys.time()
out.att.glasso <- glasso_autotune(X=X, alpha=0.02, thr=1e-4)
```

    ## glasso iter = 1; error = 1e+06
    ## glasso iter = 2; error = 0.268
    ## glasso iter = 3; error = 0.003
    ## glasso iter = 4; error = 0.003
    ## final glasso iter = 4

``` r
(Sys.time()-start.T )
```

    ## Time difference of 0.01940703 secs

``` r
round(out.att.glasso$Theta,4)
```

    ##         [,1]   [,2]   [,3]   [,4]   [,5]   [,6]   [,7]   [,8]   [,9]  [,10]
    ##  [1,] 0.7329 0.1153 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
    ##  [2,] 0.1153 0.6474 0.0933 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
    ##  [3,] 0.0000 0.0933 0.6348 0.0892 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
    ##  [4,] 0.0000 0.0000 0.0892 0.7712 0.1050 0.0000 0.0000 0.0000 0.0000 0.0000
    ##  [5,] 0.0000 0.0000 0.0000 0.1050 0.7860 0.0724 0.0000 0.0000 0.0000 0.0000
    ##  [6,] 0.0000 0.0000 0.0000 0.0000 0.0724 0.6977 0.0834 0.0000 0.0000 0.0000
    ##  [7,] 0.0000 0.0000 0.0000 0.0000 0.0000 0.0834 0.6626 0.1159 0.0000 0.0000
    ##  [8,] 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.1159 0.6550 0.0345 0.0000
    ##  [9,] 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0345 0.6537 0.1196
    ## [10,] 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.1196 0.7490
