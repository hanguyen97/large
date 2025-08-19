rm(list = ls())

setwd("/Users/hanguyen/Documents/Cornell/Projects/FC_SC/Code/Simulations/")

library(huge)
library(MASS)
library(corrplot)
source("library_test.R")

# ------------- DGPs -------------# 
# simulate data - 100 replications
# p = {100}, n = {200, 500, 1000}
# graph type: random, band, cluster

set.seed(1)

# replication
R <- 2
p <- c(100, 300)
dgp <- data.frame(p = rep(p,2), n=c(300, 300, 500, 500))

# graph settings
# graph.type <- c("band", "toeplitz", "random", "random", "scale-free") 
graph.type <- c("band", "band", "random") 
graph.param <- c(1, 2, 0.01)

graph <- cbind(graph.type=graph.type, 
               graph.param=graph.param)
dgp.df <- merge(dgp, graph)

dgp.df$same.diag <- 0

dgp.df$graph.param <- as.numeric(dgp.df$graph.param)
dgp.df <- dgp.df[order(dgp.df$graph.type, dgp.df$graph.param, 
                      dgp.df$same.diag, dgp.df$p, dgp.df$n), ]  



# ------------- simulate all data -------------# 
diag_sc <- c(10, 1, 0.5) 

for (i in 1:nrow(dgp.df)) {
  graph.type <- dgp.df[i, ]$graph.type
  graph.param <- dgp.df[i, ]$graph.param
  p <- dgp.df[i, ]$p
  n <- dgp.df[i, ]$n
  same.diag <- dgp.df[i, ]$same.diag
  
  print(paste(i, graph.type, graph.param, p, n, same.diag))
  data.name <- paste0(graph.type, graph.param, 
                      "_p", p, "_n", n, "_diag", same.diag, "_R", R)
  Theta_list <- list()
  X_list <- list()
  
  for (r in 1:R) {
    if (graph.type=="random") {
      Theta <- get_random_graph(p=p, same.diag=same.diag, prob=graph.param) 
    } else if (graph.type=="band" && graph.param==1) {
      Theta <- get_band1_graph(p=p, same.diag=same.diag)
    } else if (graph.type=="band" && graph.param==2) {
      Theta <- get_band2_graph(p=p, same.diag=same.diag)
    } else if (graph.type=="scale-free") {
      Theta <- get_scalefree_graph(p=p, same.diag=same.diag)
    } else if (graph.type=="toeplitz") {
      Theta <- get_toeplitz_graph(p=p, same.diag=same.diag)
    }
    
    if (min(eigen(Theta)$values)<=1e-3) {
      stop("Theta is not positive definite")
    }
    Theta_list[[r]] <- Theta
    
    Sigma <- solve(Theta)
    if (min(eigen(Sigma)$values)<=1e-3) {
      print("Sigma is not positive definite")
    }
    X_list[[r]] <- mvrnorm(n=n, mu=rep(0,p), Sigma=Sigma)
    
  }
  
  save(X_list, Theta_list, file=paste0("out/", data.name, "_data.RData"))
}

save(dgp.df, file=paste0("out/DGP.RData"))



