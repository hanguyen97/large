rm(list = ls())
library(corrplot)

# Load data
setwd("~/ATTglasso/")
load("out/DGP.RData")

pdf(file = paste0("out/rmse_V7.pdf"), width = 8, height = 6)

dgp.df <- dgp.df[order(dgp.df$graph.type, dgp.df$p, dgp.df$n),]

R <- 20
for (i in 1:nrow(dgp.df)) {
  graph.type <- dgp.df[i, ]$graph.type
  graph.param <- dgp.df[i, ]$graph.param
  p <- dgp.df[i, ]$p
  n <- dgp.df[i, ]$n
  same.diag <- dgp.df[i, ]$same.diag
  if (same.diag != 0) { next }
  
  print(paste(i, graph.type, graph.param, p, n, same.diag, R))
  data.name <- paste0(graph.type, graph.param, 
                      "_p", p, "_n", n, "_diag", same.diag, "_R", R)
    
  load(paste0("out/", data.name, "_data.RData"))
  load(paste0("out/", data.name, "_total.arr.RData"))
  load(paste0("out/", data.name, "_compare.arr.RData"))
  
  # Calculate mean performance 
  total.arr$mean.performance <- rowMeans(total.arr[,3:(R+2)])
  
  # Subset data
  att.glasso <- subset(compare.arr, method=="att.glassoV4")
  huge.ebic <- subset(compare.arr, method=="huge.ebic")
  huge.ric <- subset(compare.arr, method=="huge.ric")
  cv.loglik <- subset(compare.arr, method=="cv.loglik")
  huge.stars <- subset(compare.arr, method=="huge.stars")
  # if we reduce n, check CV performance
  
  # Plot RMSE
  # pdf(file = paste0("out/rmse_", data.name, ".pdf"), width = 8, height = 6)
  
  corrplot(cov2cor(Theta_list[[1]]), method = "color", tl.pos="n")
  
  plot(log2(lambda.grid), total.arr$mean.performance[1:l], type="l",
       xlab="log2(lambda)", ylab="RMSE", ylim=c(0,1.1),
       main=gsub("_", " ", data.name))
  
  points(log2(huge.ebic$sel.lambda.max), huge.ebic$rmse, col="red", pch=19)
  abline(h = median(huge.ebic$rmse), col = "red", lwd = 2, lty = 2)  
  
  points(log2(huge.ric$sel.lambda.max), huge.ric$rmse, col="orange", pch=19)
  abline(h = median(huge.ric$rmse), col = "orange", lwd = 2, lty = 2)  
  
  points(log2(cv.loglik$sel.lambda.max), cv.loglik$rmse, col="green", pch=19)
  abline(h = median(cv.loglik$rmse), col = "green", lwd = 2, lty = 2)  
  
  points(log2(huge.stars$sel.lambda.max), huge.stars$rmse, col="purple", pch=19)
  abline(h = median(huge.stars$rmse), col = "purple", lwd = 2, lty = 2)  
  
  points(rep(-13,R), att.glasso$rmse, col="blue", pch=3)
  # points(log2(att.glasso$sel.lambda.max), att.glasso$rmse, col="blue", pch=3)
  # points(log2(att.glasso$sel.lambda.min), att.glasso$rmse, col="blue", pch=4)
  abline(h = median(att.glasso$rmse), col = "blue", lwd = 2, lty = 2) 
      
  # Add legend
  legend("topright", title="median RMSE",
         legend = c("huge.ebic", "huge.ric", "cv.loglik", "huge.stars", "att.glassoV4"),
         col = c("red", "orange", "green", "purple",  "blue"),
         lwd = 1,
         lty = 1)
  
}

dev.off()





