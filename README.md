ATTglasso
================

## Overview

`ATTglasso` is an R package that performs Autotune Graphical LASSO.

## Installation

``` r
# You can install the development version from GitHub:
# install.packages("devtools")
devtools::install_github("hanguyen97/ATTglasso")
```

## Usage

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

# Run autotune GLASSO
library(ATTglasso)
start.T <- Sys.time()
out.att.glasso <- glasso_autotune(X=X, alpha=0.02, thr=1e-4)
```

    ## glasso iter = 1; error = 1e+06
    ## change in W.err = 999999
    ## glasso iter = 2; error = 0.978
    ## change in W.err = 0.959196
    ## glasso iter = 3; error = 0.019
    ## change in W.err = 0.0189771
    ## glasso iter = 4; error = 0
    ## change in W.err = 0.000213065
    ## glasso iter = 5; error = 0
    ## change in W.err = 4.99287e-06
    ## glasso iter = 6; error = 0
    ## change in W.err = 4.99101e-06
    ## final glasso iter = 6

``` r
(Sys.time()-start.T )
```

    ## Time difference of 0.02983999 secs

``` r
round(out.att.glasso$Theta,4)
```

    ##         [,1]   [,2]   [,3]   [,4]   [,5]   [,6]   [,7]   [,8]   [,9]  [,10]
    ##  [1,] 0.8766 0.1595 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
    ##  [2,] 0.1595 0.7612 0.1248 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
    ##  [3,] 0.0000 0.1248 0.7378 0.1181 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
    ##  [4,] 0.0000 0.0000 0.1181 0.8918 0.1379 0.0000 0.0000 0.0000 0.0000 0.0000
    ##  [5,] 0.0000 0.0000 0.0000 0.1379 0.9050 0.0951 0.0000 0.0000 0.0000 0.0000
    ##  [6,] 0.0000 0.0000 0.0000 0.0000 0.0951 0.8108 0.1137 0.0000 0.0000 0.0000
    ##  [7,] 0.0000 0.0000 0.0000 0.0000 0.0000 0.1137 0.7970 0.1635 0.0000 0.0000
    ##  [8,] 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.1635 0.7842 0.0477 0.0000
    ##  [9,] 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0477 0.7723 0.1685
    ## [10,] 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.1685 0.9059
