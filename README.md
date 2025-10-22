LARGE: Locally Adaptive Regularization for Graph Estimation
================

## Introduction

We provide the `large` package for learning adaptive tuning parameters
in graph estimation. Given $n$ i.i.d. observations of
$X \sim N_p(0, \Theta^{-1})$, LARGE estimates the precision matrix
$\Theta$ of a Gaussian graphical model (GGM) by maximizing the
$\ell_1$-penalized log-likelihood with column-wise penalties:
<p>
$$
    \hat{\Theta}_{\mathrm{LARGE}} = \underset{\Theta \succeq 0}{\arg\max} \left\{ \log \det(\Theta) - \mathrm{trace}(S\Theta) - \sum_{j=1}^p \lambda_j \|\Theta_{\cdot j}\|_{1, \mathrm{off}} \right\},
$$
</p>

where $S = \frac{1}{n} \sum_{i=1}^n X_i X_i^\top$ is the sample
covariance matrix, and
$\|\Theta_{\cdot j}\|_{1, \mathrm{off}} = \sum_{k \ne j} |\Theta_{kj}|$.

Unlike standard GLASSO, which relies on a single global penalty, LARGE
adaptively learns a set of nodewise penalties
$\lambda_j, j = 1, \ldots, p$. It does so by augmenting the nodewise
Lasso regression step to jointly estimate both regression coefficients
and error variances, allowing more flexible and data-driven
regularization across nodes.

## Installation

You can download the `large` package from Github.

``` r
# You can install the development version from GitHub:
# install.packages("devtools")
devtools::install_github("hanguyen97/large")
```

## Quick Start:

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

We can estimate the precision matrix using `large`. $\alpha = 0.02$
denotes the significance level of the sequential F-test procedure used
for edge selection at each nodewise Lasso regression step.

``` r
library(large)
start.T <- Sys.time()
out <- fit_large(X=X, alpha=0.02, thr=1e-4)
```

    ## iter = 1; relative F-norm change = 1e+06
    ## iter = 2; relative F-norm change = 0.201
    ## iter = 3; relative F-norm change = 0.005
    ## iter = 4; relative F-norm change = 0.005
    ## final iter = 4

``` r
(Sys.time()-start.T )
```

    ## Time difference of 0.01915598 secs

``` r
round(out$Theta,4)
```

    ##         [,1]   [,2]   [,3]   [,4]   [,5]   [,6]   [,7]   [,8]   [,9]  [,10]
    ##  [1,] 0.9501 0.1784 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
    ##  [2,] 0.1784 0.8542 0.1386 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
    ##  [3,] 0.0000 0.1386 0.8217 0.1294 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
    ##  [4,] 0.0000 0.0000 0.1294 0.9799 0.1507 0.0000 0.0000 0.0000 0.0000 0.0000
    ##  [5,] 0.0000 0.0000 0.0000 0.1507 0.9906 0.1051 0.0000 0.0000 0.0000 0.0000
    ##  [6,] 0.0000 0.0000 0.0000 0.0000 0.1051 0.8989 0.1324 0.0000 0.0000 0.0000
    ##  [7,] 0.0000 0.0000 0.0000 0.0000 0.0000 0.1324 0.9346 0.1860 0.0000 0.0000
    ##  [8,] 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.1860 0.8923 0.0534 0.0000
    ##  [9,] 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0534 0.8687 0.1847
    ## [10,] 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.1847 0.9928
