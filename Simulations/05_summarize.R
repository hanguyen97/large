rm(list = ls())

setwd("~/ATTglasso/")

library(data.table)

load("out/DGP.RData")
dgp.df <- data.table(dgp.df)

R <- 20

perf <- data.table()

for (i in c(1:12)) {
  graph.type <- dgp.df[i, ]$graph.type
  graph.param <- dgp.df[i, ]$graph.param
  p <- dgp.df[i, ]$p
  n <- dgp.df[i, ]$n
  same.diag <- dgp.df[i, ]$same.diag
  
  print(paste(i, graph.type, graph.param, p, n, same.diag, R))
  data.name <- paste0(graph.type, graph.param, 
                      "_p", p, "_n", n, "_diag", same.diag, "_R", R)

  load(paste0("out/", data.name, "_compare.arr.RData"))
  compare.arr$dgp <- data.name
  perf <- rbind(perf, compare.arr)
  
}

# Output all metrics
write.csv(perf, "out/all_perf.csv")

# Summarize 
metrics <- c("rmse", "auroc", "pre", "recall", "acc", "f1", "run.time", "converged")
sum_tbl <- perf[, lapply(.SD, function(x) paste0(round( mean(x),2), " (", round(sd( x, na.rm=TRUE),2), ")")), 
          by=c("dgp", "method"), .SDcols = metrics]
write.csv(sum_tbl, "out/sum_perf.csv")

# metrics <- c("rmse", "auroc", "pre", "recall", "acc", "f1", "run.time", "converged")
# sum_tbl <- perf[, lapply(.SD, function(x) paste0(round( mean(x),3), " (", round(sd( x, na.rm=TRUE),3), ")")), 
#           by=c("dgp", "method"), .SDcols = metrics]

# # Output all metrics
# write.csv(sum_tbl, "out/sum_tbl_test.csv")

# # AUROC only
# short_tbl <- sum_tbl[, c("dgp", "method","auroc")]
# short_tbl[, n := as.numeric(sub(".*_n(\\d+)_.*", "\\1", dgp))]
# short_tbl[, p := as.numeric(sub(".*_p(\\d+)_.*", "\\1", dgp))]
# short_tbl[, dgp := gsub("_p[0-9]+_n[0-9]+", "", dgp)]

# auroc_tbl <- dcast(short_tbl, dgp+method ~p+n, value.var = "auroc")
# write.csv(auroc_tbl, "out/auroc_test.csv")
