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
    ## glasso iter = 2; error = 1.188
    ## change in W.err = 1.17184
    ## glasso iter = 3; error = 0.016
    ## change in W.err = 0.0163259
    ## glasso iter = 4; error = 0
    ## change in W.err = 0.000113608
    ## glasso iter = 5; error = 0
    ## change in W.err = 1.95857e-06
    ## glasso iter = 6; error = 0
    ## change in W.err = 1.95824e-06
    ## final glasso iter = 6

``` r
(Sys.time()-start.T )
```

    ## Time difference of 0.03056479 secs

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
