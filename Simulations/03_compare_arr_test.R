rm(list = ls())

setwd("~/ATTglasso/")
library(foreach)
library(doMC)
registerDoMC(10)

options(cores = 10)

# Load data
load("out/DGP.RData")
R <- 10
for (i in 7) {
  graph.type <- dgp.df[i, ]$graph.type
  graph.param <- dgp.df[i, ]$graph.param
  p <- dgp.df[i, ]$p
  n <- dgp.df[i, ]$n
  same.diag <- dgp.df[i, ]$same.diag
  
  print(paste(i, graph.type, graph.param, p, n, same.diag, R))
  data.name <- paste0(graph.type, graph.param, 
                      "_p", p, "_n", n, "_diag", same.diag, "_R", R)
  
  load(paste0("out/", data.name, "_data.RData"))
  
  # Run autotune glasso and benchmark
  compare.arr = foreach(r = 1:R, .combine = rbind) %dopar% {
    
    set.seed(r)
    source("library_test.R")
    
    library(huge)
    library(glasso)
    library(CVglasso)
    library(ATTglasso)
    
    # Get ground truth
    Theta <- Theta_list[[r]]
    Theta.off <- Theta
    diag(Theta.off) <- 0
    X <- X_list[[r]]
    
    # Run autotune glasso 1%
    startT <- Sys.time()
    out.att.glasso <- glasso_autotune(X=X, alpha=0.01, thr=1e-3, maxit=100, penalize_diag=TRUE)
    endT <- Sys.time()-startT
    Theta.hat.att.glasso <- out.att.glasso$Theta
    rmse.att.glasso <- rel.diff(Est=Theta.hat.att.glasso, Tru=Theta)

    Theta.hat.off.att.glasso <- Theta.hat.att.glasso
    diag(Theta.hat.off.att.glasso) <- 0
    rmse.off.att.glasso <- rel.diff(Est=Theta.hat.off.att.glasso, Tru=Theta.off)

    auroc.att.glasso <- ROC.Theta(Theta.hat=Theta.hat.att.glasso, Theta=Theta)$auc
    perf.att.glasso <- performance(Theta=Theta, Theta.hat=Theta.hat.att.glasso)
    temp.att.glasso_a1 <- data.frame(method="att.glasso_a1", 
                                  rep=r, 
                                  # sel.lambda.max=max(out.att.glasso$lambda), 
                                  # sel.lambda.min=min(out.att.glasso$lambda), 
                                  sel.lambda.max = NA,
                                  sel.lambda.min = NA,
                                  rmse=rmse.att.glasso,
                                  rmse.off=rmse.off.att.glasso,
                                  auroc=auroc.att.glasso,
                                  pre=perf.att.glasso[1],
                                  recall=perf.att.glasso[2],
                                  acc=perf.att.glasso[3],
                                  f1=perf.att.glasso[4],
                                  fdr=perf.att.glasso[5],
                                  converged=out.att.glasso$converged,
                                  run.time=endT)       
    
    # Run autotune glasso 2%
    startT <- Sys.time()
    out.att.glasso <- glasso_autotune(X=X, alpha=0.02, thr=1e-3, maxit=100, penalize_diag=TRUE)
    endT <- Sys.time()-startT
    Theta.hat.att.glasso <- out.att.glasso$Theta
    rmse.att.glasso <- rel.diff(Est=Theta.hat.att.glasso, Tru=Theta)

    Theta.hat.off.att.glasso <- Theta.hat.att.glasso
    diag(Theta.hat.off.att.glasso) <- 0
    rmse.off.att.glasso <- rel.diff(Est=Theta.hat.off.att.glasso, Tru=Theta.off)

    auroc.att.glasso <- ROC.Theta(Theta.hat=Theta.hat.att.glasso, Theta=Theta)$auc
    perf.att.glasso <- performance(Theta=Theta, Theta.hat=Theta.hat.att.glasso)
    temp.att.glasso_a2 <- data.frame(method="att.glasso_a2", 
                                  rep=r, 
                                  # sel.lambda.max=max(out.att.glasso$lambda), 
                                  # sel.lambda.min=min(out.att.glasso$lambda), 
                                  sel.lambda.max = NA,
                                  sel.lambda.min = NA,
                                  rmse=rmse.att.glasso,
                                  rmse.off=rmse.off.att.glasso,
                                  auroc=auroc.att.glasso,
                                  pre=perf.att.glasso[1],
                                  recall=perf.att.glasso[2],
                                  acc=perf.att.glasso[3],
                                  f1=perf.att.glasso[4],
                                  fdr=perf.att.glasso[5],
                                  converged=out.att.glasso$converged,
                                  run.time=endT)       
    
    # Run autotune glasso 0.05
    startT <- Sys.time()
    out.att.glasso <- glasso_autotune(X=X, alpha=0.05, thr=1e-3, maxit=100, penalize_diag=TRUE)
    endT <- Sys.time()-startT
    Theta.hat.att.glasso <- out.att.glasso$Theta
    rmse.att.glasso <- rel.diff(Est=Theta.hat.att.glasso, Tru=Theta)

    Theta.hat.off.att.glasso <- Theta.hat.att.glasso
    diag(Theta.hat.off.att.glasso) <- 0
    rmse.off.att.glasso <- rel.diff(Est=Theta.hat.off.att.glasso, Tru=Theta.off)

    auroc.att.glasso <- ROC.Theta(Theta.hat=Theta.hat.att.glasso, Theta=Theta)$auc
    perf.att.glasso <- performance(Theta=Theta, Theta.hat=Theta.hat.att.glasso)
    temp.att.glasso_a5 <- data.frame(method="att.glasso_a5", 
                                  rep=r, 
                                  # sel.lambda.max=max(out.att.glasso$lambda), 
                                  # sel.lambda.min=min(out.att.glasso$lambda), 
                                  sel.lambda.max = NA,
                                  sel.lambda.min = NA,
                                  rmse=rmse.att.glasso,
                                  rmse.off=rmse.off.att.glasso,
                                  auroc=auroc.att.glasso,
                                  pre=perf.att.glasso[1],
                                  recall=perf.att.glasso[2],
                                  acc=perf.att.glasso[3],
                                  f1=perf.att.glasso[4],
                                  fdr=perf.att.glasso[5],
                                  converged=out.att.glasso$converged,
                                  run.time=endT)       
    
    # # Run huge ebic
    # startT <- Sys.time()
    # out.huge <- huge(X, method = "glasso")
    # out.huge.ebic <- huge.select(out.huge, criterion = "ebic")
    # # rmse.huge.ebic <- rel.diff(out.huge.ebic$opt.icov, Theta)
    # Theta.hat.huge.ebic <- glasso(s=cov(X),rho=out.huge.ebic$opt.lambda)$wi
    # rmse.huge.ebic <- rel.diff(Est=Theta.hat.huge.ebic, Tru=Theta)

    # Theta.hat.off.huge.ebic <- Theta.hat.huge.ebic
    # diag(Theta.hat.off.huge.ebic) <- 0
    # rmse.off.huge.ebic <- rel.diff(Est=Theta.hat.off.huge.ebic, Tru=Theta.off)

    # auroc.huge.ebic <- ROC.Theta(Theta.hat=Theta.hat.huge.ebic, Theta=Theta)$auc
    # perf.huge.ebic <- performance(Theta=Theta, Theta.hat=Theta.hat.huge.ebic)
    # temp.huge.ebic <- data.frame(method="huge.ebic", 
    #                              rep=r, 
    #                              sel.lambda.max=out.huge.ebic$opt.lambda, 
    #                              sel.lambda.min=out.huge.ebic$opt.lambda, 
    #                              rmse=rmse.huge.ebic,
    #                              rmse.off=rmse.off.huge.ebic,
    #                              auroc=auroc.huge.ebic,
    #                              pre=perf.huge.ebic[1],
    #                              recall=perf.huge.ebic[2],
    #                              acc=perf.huge.ebic[3],
    #                              f1=perf.huge.ebic[4],
    #                              converged=TRUE,
    #                              run.time=Sys.time()-startT) 
    
    # # Run huge ric
    # startT <- Sys.time()
    # out.huge <- huge(X, method = "glasso")
    # out.huge.ric <- huge.select(out.huge, criterion = "ric")
    # if (is.null(out.huge.ric$opt.icov)) {
    #   rmse.huge.ric <- NA
    #   sel.lambda.huge.ric <- NA
    #   converged.huge.ric <- FALSE
    # } else {
    #   # rmse.huge.ric <- rel.diff(out.huge.ric$opt.icov, Theta)
    #   Theta.hat.huge.ric <- glasso(s=cov(X),rho=out.huge.ric$opt.lambda)$wi
    #   rmse.huge.ric <- rel.diff(Est=Theta.hat.huge.ric, Tru=Theta)

    #   Theta.hat.off.huge.ric <- Theta.hat.huge.ric
    #   diag(Theta.hat.off.huge.ric) <- 0
    #   rmse.off.huge.ric <- rel.diff(Est=Theta.hat.off.huge.ric, Tru=Theta.off)

    #   auroc.huge.ric <- ROC.Theta(Theta.hat=Theta.hat.huge.ric, Theta=Theta)$auc
    #   sel.lambda.huge.ric <- out.huge.ric$opt.lambda
    #   perf.huge.ric <- performance(Theta=Theta, Theta.hat=Theta.hat.huge.ric)
    #   converged.huge.ric <- TRUE
    # }
    # temp.huge.ric <- data.frame(method="huge.ric",
    #                             rep=r,
    #                             sel.lambda.max=sel.lambda.huge.ric,
    #                             sel.lambda.min=sel.lambda.huge.ric,
    #                             rmse=rmse.huge.ric,
    #                             rmse.off=rmse.off.huge.ric,
    #                             auroc=auroc.huge.ric,
    #                             pre=perf.huge.ric[1],
    #                             recall=perf.huge.ric[2],
    #                             acc=perf.huge.ric[3],
    #                             f1=perf.huge.ric[4],
    #                             converged=converged.huge.ric,
    #                             run.time=Sys.time()-startT)
    
    # Run CVglasso loglik - use original glasso function
    startT <- Sys.time()
    out.cv.loglik = CVglasso(X=X, S=cov(X), crit.cv="loglik", K=5, diagonal=TRUE, trace="none")
    endT <- Sys.time()-startT
    # Theta.hat.cv.loglik <- out.cv.loglik$Omega
    # rmse.cv.loglik <- rel.diff(out.cv.loglik$Omega, Theta)
    Theta.hat.cv.loglik <- glasso(s=cov(X),rho=out.cv.loglik$Tuning[2])$wi
    rmse.cv.loglik <- rel.diff(Theta.hat.cv.loglik, Theta)

    Theta.hat.off.cv.loglik <- Theta.hat.cv.loglik
    diag(Theta.hat.off.cv.loglik) <- 0
    rmse.off.cv.loglik <- rel.diff(Est=Theta.hat.off.cv.loglik, Tru=Theta.off)

    auroc.cv.loglik <- ROC.Theta(Theta.hat.cv.loglik, Theta)$auc
    perf.cv.loglik <- performance(Theta=Theta, Theta.hat=Theta.hat.cv.loglik)
    temp.cv.loglik <- data.frame(method="cv.loglik", 
                                 rep=r, 
                                 sel.lambda.max=out.cv.loglik$Tuning[2], 
                                 sel.lambda.min=out.cv.loglik$Tuning[2], 
                                 rmse=rmse.cv.loglik,
                                 rmse.off=rmse.off.cv.loglik,
                                 auroc=auroc.cv.loglik,
                                 pre=perf.cv.loglik[1],
                                 recall=perf.cv.loglik[2],
                                 acc=perf.cv.loglik[3],
                                 f1=perf.cv.loglik[4],
                                 fdr=perf.cv.loglik[5],
                                 converged=TRUE,
                                 run.time=endT) 

    # # Run huge Stars
    # startT <- Sys.time()
    # out.huge <- huge(X, method = "glasso")
    # out.huge.stars <- huge.select(out.huge, criterion = "stars")
    # if (is.null(out.huge.stars$opt.icov)) {
    #   rmse.huge.stars <- NA
    #   sel.lambda.huge.stars <- NA
    # } else {
    #   # rmse.huge.stars <- rel.diff(out.huge.stars$opt.icov, Theta)
    #   Theta.hat.huge.stars <- glasso(s=cov(X),rho=out.huge.stars$opt.lambda)$wi
    #   rmse.huge.stars <- rel.diff(Est=Theta.hat.huge.stars, Tru=Theta)

    #   Theta.hat.off.huge.stars <- Theta.hat.huge.stars
    #   diag(Theta.hat.off.huge.stars) <- 0
    #   rmse.off.huge.stars <- rel.diff(Est=Theta.hat.off.huge.stars, Tru=Theta.off)

    #   auroc.huge.stars <- ROC.Theta(Theta.hat=Theta.hat.huge.stars, Theta=Theta)$auc
    #   sel.lambda.huge.stars <- out.huge.stars$opt.lambda
    #   perf.huge.stars <- performance(Theta=Theta, Theta.hat=Theta.hat.huge.stars)
    # }
    # temp.huge.stars <- data.frame(method="huge.stars",
    #                             rep=r,
    #                             sel.lambda.max=sel.lambda.huge.stars,
    #                             sel.lambda.min=sel.lambda.huge.stars,
    #                             rmse=rmse.huge.stars,
    #                             rmse.off=rmse.off.huge.stars,
    #                             auroc=auroc.huge.stars,
    #                             pre=perf.huge.stars[1],
    #                             recall=perf.huge.stars[2],
    #                             acc=perf.huge.stars[3],
    #                             f1=perf.huge.stars[4],
    #                             converged=TRUE,
    #                             run.time=Sys.time()-startT)

    temp <- rbind(temp.att.glasso_a1, temp.att.glasso_a2, 
                  temp.att.glasso_a5, 
                  # temp.huge.ebic, 
                  # temp.huge.ric, 
                  # temp.huge.stars,
                  temp.cv.loglik)
    
    temp
  }
  
  save(compare.arr, file=paste0("out/", data.name, "_compare.arr.RData"))
  
}
