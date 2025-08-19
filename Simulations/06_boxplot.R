rm(list = ls())

setwd("~/ATTglasso/")

library(data.table)
library(ggplot2)

load("out/DGP.RData")
dgp.df <- data.table(dgp.df)

df <- data.table(read.csv("out/all_perf.csv", header=TRUE))
# rename
df[dgp=="band1_p100_n50_diag0_R10", dgp:="band1_p100_n050_diag0_R10"]
df[dgp=="band2_p100_n50_diag0_R20", dgp:="band2_p100_n050_diag0_R20"]
df[dgp=="toeplitz0_p100_n50_diag0_R20", dgp:="toeplitz0_p100_n050_diag0_R20"]

# Convert converged to factor for coloring
df$converged <- as.factor(df$converged)

# Define the plotting function
plot_box <- function(metric, diag=0, y_scale="fixed") {
  tmp_df <- df[grepl(paste0("diag", diag), df$dgp),]
  ggplot(tmp_df, aes(x = method, y = as.numeric(get(metric)))) +
    geom_boxplot(outlier.shape = NA) +  # Boxplot without outliers
    geom_jitter(width = 0.2, alpha = 0.7, aes(color = converged)) +  # Add jittered points
    facet_wrap(~ dgp, scales = y_scale) +  # Separate plots by dgp
    labs(y = metric, title = paste("Boxplot of", metric, "by Method")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x labels
    scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "red"))  # Colors for convergence
}

box_plot_rmse_no_outlier <- function(metric, diag=0) {
  tmp_df <- df[grepl(paste0("diag", diag), df$dgp) & get(metric) < 2,]
  ggplot(tmp_df, aes(x = method, y = get(metric))) +
    geom_boxplot(outlier.shape = NA) +  # Boxplot without outliers
    # geom_jitter(width = 0.2, alpha = 0.7, aes(color = converged)) +  # Add jittered points
    facet_wrap(~ dgp, scales = "fixed") +  # Separate plots by dgp
    labs(y = metric, title = paste("Boxplot of", metric, "by Method - no outlier")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x labels
    scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "red"))  # Colors for convergence
}
# Generate the boxplots

pdf(file = paste0("out/all_boxplots_pen.pdf"), width = 8, height = 6)
plot_box("auroc", diag=0, y_scale="fixed")
# plot_box("auroc", diag=1, y_scale="fixed")

plot_box("fdr", diag=0, y_scale="fixed")
# plot_box("auroc", diag=1, y_scale="fixed")

plot_box("rmse", diag=0, y_scale="free_y")
# plot_box("rmse", diag=1, y_scale="free_y")

box_plot_rmse_no_outlier("rmse", diag=0)

plot_box("rmse.off", diag=0, y_scale="free_y")
# plot_box("rmse.off", diag=1, y_scale="free_y")

box_plot_rmse_no_outlier("rmse.off", diag=0)

plot_box("run.time", diag=0, y_scale="free_y")
# plot_box("run.time", diag=1, y_scale="fixed")

dev.off()