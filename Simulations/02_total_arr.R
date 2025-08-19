library(foreach)
library(doMC)
registerDoMC(40)
options(cores = 40)

# Load data
setwd("~/ATTglasso/")
load("out/DGP.RData")
R <- 10

for (i in 10:nrow(dgp.df)) {
  graph.type <- dgp.df[i, ]$graph.type
  graph.param <- dgp.df[i, ]$graph.param
  p <- dgp.df[i, ]$p
  n <- dgp.df[i, ]$n
  same.diag <- dgp.df[i, ]$same.diag
  
  print(paste(i, graph.type, graph.param, p, n, same.diag, R))
  data.name <- paste0(graph.type, graph.param, 
                      "_p", p, "_n", n, "_diag", same.diag, "_R", R )

  load(paste0("out/", data.name, "_data.RData"))
  
  R <- 20
  lower = -13
  upper = 3
  lambda.grid <- 2^seq(upper, lower, by = -0.1); l <- length(lambda.grid)
  
  # Run glasso at each lambda values
  startT <- Sys.time()
  total.arr = foreach(r = 1:R, .combine = cbind) %dopar% {
    
    source("library_test.R")
    library(glasso)
    
    Theta <- Theta_list[[r]]
    X <- X_list[[r]]
    
    rmse.arr.glasso = bic.arr.glasso = aic.arr.glasso = array(NA, l) 
    auroc.arr.glasso = array(NA, l)
    
    flag <- 0
    
    for(i in 1:l){
      # Estimate
      lambda <- lambda.grid[i]
      # print(paste("start of", i, "-th exp for", r,"-th rep, ", log2(lambda)))
      
      Theta.hat.glasso <- glasso(s=cov(X),rho=lambda)$wi
      
      if(any(is.na(Theta.hat.glasso))){
        rmse.arr.glasso[i] = auroc.arr.glasso[i] = NA
      } else {
        rmse.arr.glasso[i]<- rel.diff(Theta.hat.glasso, Theta)
        auroc.arr.glasso[i] = ROC.Theta(Theta.hat.glasso, Theta)$auc
      }
      
      i.best <- 0
      if(i==1){
        Theta.hat.best <- Theta.hat.glasso
        start.rmse <- rmse.arr.glasso[i]
        min.rmse <- NULL
        i.best <- i
      }
      
      if(!is.na(rmse.arr.glasso[i])){
        # now.bic <- bic.arr.glasso[i]
        now.rmse <- rmse.arr.glasso[i]
        min.rmse <- min(min.rmse, now.rmse)
        if(min.rmse >= now.rmse){
          i.best <- i
          Theta.hat.best <- Theta.hat.glasso
        }
        
        if(now.rmse > min.rmse && now.rmse - min.rmse >= 0.4 * (start.rmse - min.rmse))
          flag <- 1
      }
      if(flag == 1)
        break
      
      # print(paste("end of", i, "-th exp for", r,"-th rep,", log2(lambda)))
    }
    temp <- data.frame(c(rmse.arr.glasso, 
                         auroc.arr.glasso
    )) # , rmse.arr.clime
    
    # save(temp, file = paste0("total", r, ".RData"))
    temp
  }
  print(Sys.time() - startT)
  
  # Save output
  metrics <- c(rep("rmse", l), rep("auroc", l))
  total.arr <- cbind(metrics, log2(lambda.grid), total.arr)
  colnames(total.arr) <- c("metric", "log_lambda", paste0("rep", 1:R))
  save(total.arr, lambda.grid, l, file=paste0("out/", data.name, "_total.arr.RData"))
  
}



