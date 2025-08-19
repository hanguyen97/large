library(pROC)
library(plot3D)
library(gdata)
library(clime)
library(Matrix)

col_palette <- gray((100:0/100)^1.5)

rel.diff <- function(Est, Tru){
  return(sum((Est - Tru)^2)/sum(Tru^2))
}

ROC.Theta <- function(Theta.hat, Theta){
  truth <- (!!upperTriangle(round(Theta, 5), diag = FALSE)) + 0
  pred <- upperTriangle(round((Theta.hat), 5), diag = FALSE)
  roc.obj <- pROC::roc(truth, pred, direction = "<", levels = c(0, 1))
  return(roc.obj)
}

performance <- function(Theta, Theta.hat){
  truth <- as.factor((!!upperTriangle(round(Theta, 5), diag = FALSE)) + 0)
  preds <- as.factor((!!upperTriangle(round(Mod(Theta.hat), 5), diag = FALSE)) + 0)
  t <- table(truth, preds)
  if(ncol(t) == 1)
    t <- cbind(t, c(0, 0))
  a <- t[1, 1]
  b <- t[1, 2]
  c <- t[2, 1]
  d <- t[2, 2]

  pre <- ifelse(b+d, d/(b+d), NA)
  recall <- ifelse(c+d, d/(c+d), NA)
  acc <- ifelse(a+b+c+d, (a+d)/(a+b+c+d), NA)
  f1 <- ifelse(b+c+2*d, d/(d + (b+c)/2), NA)
  fdr <- ifelse(b+d, b/(b+d), NA) 

  return(c(pre, recall, acc, f1, fdr, tp=d, fp=b, fn=c, tn=a))
}

get_band_graph <- function(p, band = 1, rho1=0.3, rho2=0.3, same.diag = 1) {
  if (p < 3) stop("Matrix size must be at least 3.")
  
  M <- matrix(0, p, p)  # Initialize an p x p matrix with zeros
  
  if (same.diag > 0) {
    diag(M) <- same.diag  
    for (i in 1:(p-1)) {
      M[i, i+1] <- rho1  # Upper diagonal
      M[i+1, i] <- rho1  # Lower diagonal (symmetry)
    }
    if (band == 2) {
      for (i in 1:(p-2)) {
        M[i, i+2] <- rho2  # Second upper diagonal
        M[i+2, i] <- rho2  # Second lower diagonal (symmetry)
      }
    }
  } else {
    blocks <- c(20, 10, 1)  # Different diagonal values for blocks
    block_sizes <- sort(rep(blocks, length.out = p), decreasing=TRUE)  # Repeat to match p
    diag(M) <- block_sizes  # Assign diagonal values
    for (i in 1:(p-1)) {
      if (M[i, i+1] == 0) {
        M[i, i+1] <- rho1 * min(block_sizes[i], block_sizes[i+1])  # Normalize with min()
        M[i+1, i] <- M[i, i+1]  # Lower diagonal (symmetry)
      }
    }
    if (band == 2) {
      for (i in 1:(p-2)) {
        if (M[i, i+2] == 0) {
          M[i, i+2] <- rho2 * min(block_sizes[i], block_sizes[i+2])  # Normalize with min()
          M[i+2, i] <- M[i, i+2]  # Lower diagonal (symmetry)
        }
      }
    }
  }
  
  return(M)
}

get_scalefree_graph <- function(p, same.diag = 1) {
  Omega_0 <- huge.generator(d=p, graph="scale-free")$omega
  
  # Standardize to have unit diagonal
  if (same.diag > 0) {
    D <- sqrt(diag(Omega_0))
    Omega_0 <- diag(sqrt(same.diag)/D) %*% Omega_0 %*% diag(sqrt(same.diag)/D)
  } else {
    blocks <- c(1, 2, 4)  # Different diagonal values for blocks
    block_sizes <- sort(rep(blocks, length.out = p))  # Repeat to match p
    D <- sqrt(diag(Omega_0))
    Omega_0 <- diag(sqrt(block_sizes)/D) %*% Omega_0 %*% diag(sqrt(block_sizes)/D)
  }
  
  return(Omega_0)
}


get_random_graph <- function(p, prob=0.02, rho=0.3, same.diag = 0) {
  if (same.diag == 0) {
    sc <- c(10, 1, 0.5) 
    p1 <- round(p/3)
    p2 <- round(p/3)
    p3 <- p - (p1+p2)
    p_blocks <- c(p1, p2, p3)
    
    Theta_blocks <- list()
    for (i in 1:3) {
      p <- p_blocks[i]
      B <- matrix(0, p, p)
      off_diag_indices <- which(upper.tri(B), arr.ind = TRUE)
      random_vals <- runif(nrow(off_diag_indices)) < prob
      B[off_diag_indices] <- rho * random_vals
      B <- B + t(B)  # Make symmetric
      diag(B) <- 1
      Theta_blocks[[i]] <- sc[i] * B
    }
    
    Theta <- as.matrix(bdiag(Theta_blocks[[1]], Theta_blocks[[2]], Theta_blocks[[3]]))
    return(Theta)
    
  }
}


