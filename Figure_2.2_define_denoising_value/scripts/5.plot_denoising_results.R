#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

###
# write.table(
#    mean_CNV_df,
#    paste0(table_dir, "mean_CNV_score_per_cell.txt"),
#    sep = "\t",
#    quote = F,
#    col.names = T,
#    row.names = F
#  )
#accuracy_metrics <- readRDS(
#      paste0(Robject_dir, "2.accuracy_metrics.Rdata")
#    )
###

project_name <- "thesis"
subproject_name <- "Figure_2.2_define_denoising_value"
args = commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
analysis_mode <- args[2]

project_name <- "thesis"
subproject_name <- "Figure_2.2_define_denoising_value"
args = commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
analysis_mode <- args[2]

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))

library(ggplot2)
theme_set(theme_minimal())
library(dplyr)
library(cowplot)
  
lib_loc <- "/share/ScratchGeneral/jamtor/R/3.5dev/"
home_dir <- "/share/ScratchGeneral/jamtor/"

project_dir <- paste0(home_dir, "projects/", 
  project_name, "/", subproject_name, "/")
results_dir <- seurat_path <- paste0(project_dir, "results/")
func_dir <- paste0(project_dir, "scripts/functions/")

in_path <- paste0(results_dir, "infercnv/", sample_name, "/")

out_dir <- paste0(in_path, "denoising_results/")
system(paste0("mkdir -p ", out_dir))

# convert downsample_proportions to list:
downsample_list <- as.list(downsample_proportions)

for (i in 1:length(simulation_numbers)) {

  metrics_list <- lapply(downsample_list, function(x) {

    print(paste0("Fetching metrics for simulation ", 
      simulation_numbers[i], " downsampled to ", x, "..."))

    if (downsample_type == "UMI") {
      in_dir <- paste0(in_path, "/", simulation_numbers[i], "/", 
        x, "_downsampling/", analysis_mode, "_mode/Rdata/")
    } else if (downsample_type == "gene") {
      if (x != "no") {
        in_dir <- paste0(in_path, "/", simulation_numbers[i], "/", 
          x, "_gene_downsampling/", analysis_mode, "_mode/Rdata/")
      } else if (x == "no") {
        in_dir <- paste0(in_path, "/", simulation_numbers[i], "/", 
          x, "_downsampling/", analysis_mode, "_mode/Rdata/")
      }
    }
    
    accuracy_metrics <- readRDS(
      paste0(in_dir, "3.accuracy_metrics_with_correlation.Rdata")
    )

    if ("pearson_p_val" %in% rownames(accuracy_metrics)) {
      print(paste0("Correlation p-value = ", 
      accuracy_metrics["pearson_p_val",]))
    }
    
    if (x=="no") {

      metric_df <- data.frame(
        Proportion = rep(100, 5),
        Metric = c("Sensitivity", "Specificity", "Precision", "F1 score", "Correlation"),
        Score = c(
          accuracy_metrics["sensitivity",],
          accuracy_metrics["specificity",],
          accuracy_metrics["precision",],
          accuracy_metrics["F1",],
          1
        ),
        Sim_no = rep(simulation_numbers[i], 5)
      )

      metric_df$Score[nrow(metric_df)] <- 1

    } else {

      metric_df <- data.frame(
        Proportion = rep(as.numeric(x)*100, 5),
        Metric = c("Sensitivity", "Specificity", "Precision", "F1 score", "Correlation"),
        Score = c(
          accuracy_metrics["sensitivity",],
          accuracy_metrics["specificity",],
          accuracy_metrics["precision",],
          accuracy_metrics["F1",],
          accuracy_metrics["pearson_R_squared",]
        ),
        Sim_no = rep(simulation_numbers[i], 5)
      )

    }

    return(metric_df)

  })

  if (i==1) {
    plot_metrics <- do.call("rbind", metrics_list)
  } else {
    plot_metrics <- rbind(
      plot_metrics,
      do.call("rbind", metrics_list)
    )
  }

}


# split df based upon upon downsampling proportion:
split_metrics1 <- split(plot_metrics, plot_metrics$Proportion)
split_metrics2 <- lapply(split_metrics1, function(x) {
  return(split(x, x$Metric))
})
split_metrics3 <- do.call("c", split_metrics2)

final_metrics <- lapply(split_metrics3, function(x) {
  
  res_df <- data.frame(
    Proportion = x$Proportion[1],
    Metric = x$Metric[1],
    Score = mean(x$Score),
    SD = sd(x$Score),
    SE = sd(x$Score)/sqrt(nrow(x))
  )

})
final_metric_df <- do.call("rbind", final_metrics)
levels(final_metric_df$Metric) <- c("Specificity", "Precision", "Correlation",
  "Sensitivity", "F1 score")

cols <- c("#C02456", "#58B9DB", "#7CBA61", "#F4D30B", "#B066B2")

p <- ggplot(final_metric_df, aes(x = Proportion, y = Score)) 
p <- p + geom_line(aes(color = Metric))
p <- p + geom_errorbar(aes(ymin=Score-SE, ymax=Score+SE))
p <- p + scale_color_manual(values = cols[1:5])
p <- p + scale_x_continuous(
  breaks = c(0, 5, seq(10, 100, 10))
)
p <- p + scale_y_continuous(
  breaks = seq(
    round(
      min(
        plot_metrics$Score
      ), 1
    ), 1, 0.1
  )
)
p <- p + theme(
  panel.grid.minor = element_blank(),
  axis.title.x = element_text(size = 9),
  axis.title.y = element_blank(),
  legend.title = element_blank()
)
if (downsample_type == "UMI") {
  p <- p + xlab("Proportion UMIs (%)")
} else {
  p <- p + xlab("Proportion genes (%)")
}


pdf(
  paste0(out_dir, "downsampling_results.pdf"), 
  width = 7, height = 5
)
  print(p)
dev.off()

png(
  paste0(out_dir, "downsampling_results.png"), 
  width = 7, height = 5, units = "in", res = 300
)
  print(p)
dev.off()

# plot precision and specificity only:
subset_df <- final_metric_df[
  final_metric_df$Metric %in% c("Specificity", "Precision", "Correlation"),
]

p <- ggplot(subset_df, aes(x = Proportion, y = Score)) 
p <- p + geom_line(aes(color = Metric))
p <- p + geom_errorbar(aes(ymin=Score-SE, ymax=Score+SE))
p <- p + scale_color_manual(values = cols[1:3])
p <- p + scale_x_continuous(
  breaks = c(0, 5, seq(10, 100, 10))
)
p <- p + scale_y_continuous(
  breaks = seq(
    round(
      min(
        plot_metrics$Score
      ), 1
    ), 1, 0.1
  )
)
p <- p + theme(
  panel.grid.minor = element_blank(),
  axis.title.x = element_text(size = 9),
  axis.title.y = element_blank(),
  legend.title = element_blank()
)
if (downsample_type == "UMI") {
  p <- p + xlab("Proportion UMIs (%)")
} else {
  p <- p + xlab("Proportion genes (%)")
}

pdf(
  paste0(out_dir, "downsampling_specificity_precision_correlation.pdf"), 
  width = 7, height = 5
)
  print(p)
dev.off()

png(
  paste0(out_dir, "downsampling_specificity_precision_correlation.png"), 
  width = 7, height = 5, units = "in", res = 300
)
  print(p)
dev.off()

print(paste0("All output plots in", out_dir))