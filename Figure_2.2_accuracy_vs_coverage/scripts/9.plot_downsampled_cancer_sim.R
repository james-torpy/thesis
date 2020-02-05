#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

RStudio <- TRUE

project_name <- "thesis"
subproject_name <- "Figure_2.3_simulations"
args = commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
include_t_cells <- as.logical(args[2])
samples_to_include <- strsplit(args[3], split = "_")

sample_name <- "CID4520N_cancer_sim"
include_t_cells <- TRUE
downsample_proportions <- unlist(
  strsplit(
    "0.05_0.1_0.15_0.2_0.3_0.4_0.5_0.6_0.7_0.8_0.9_1.0", 
    split = "_"
  )
)

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("T-cells included? ", as.character(include_t_cells)))

library(ggplot2)
theme_set(theme_minimal())
library(dplyr)
library(cowplot)

if (RStudio) {
  
  home_dir <- "/Users/jamestorpy/clusterHome/"
  
} else {
  
  lib_loc <- "/share/ScratchGeneral/jamtor/R/3.5dev/"
  
  home_dir <- "/share/ScratchGeneral/jamtor/"
  
}

project_dir <- paste0(home_dir, "projects/", 
  project_name, "/", subproject_name, "/")
results_dir <- seurat_path <- paste0(project_dir, "results/")

if (include_t_cells) {
  in_path <- paste0(results_dir, "infercnv/t_cells_included/", 
    sample_name, "/downsampling/")
} else {
  in_path <- paste0(results_dir, "infercnv/t_cells_excluded/", 
    sample_name, "/downsampling/")
}

for (i in 1:length(downsample_proportions)) {
  print(i)

  if (downsample_proportions[i] == "1.0") {
    
    in_dir <- paste0(results_dir, "infercnv/t_cells_included/", 
      sample_name, "/downsampling/tables/")
    
    if (i==1) {
      accuracy_metrics <- read.table(
        paste0(in_dir, "accuracy_metrics.txt"),
        header = F, sep = "\t", as.is = T
      )
      print(paste0("Correlation p-value = ", 
        accuracy_metrics$V2[accuracy_metrics$V1 == "pearson_p_val"]))
      plot_metrics <- data.frame(
        Proportion = rep(as.numeric(downsample_proportions[i])*100, 3),
        Metric = c("Sensitivity", "Specificity", "Correlation"),
        Score = c(
          accuracy_metrics$V2[accuracy_metrics$V1 == "sensitivity"],
          accuracy_metrics$V2[accuracy_metrics$V1 == "specificity"],
          1
        )
      )
    } else {
      accuracy_metrics <- read.table(
        paste0(in_dir, "accuracy_metrics.txt"),
        header = F, sep = "\t", as.is = T
      )
      print(paste0("Correlation p-value = ", 
                   accuracy_metrics$V2[accuracy_metrics$V1 == "pearson_p_val"]))
      plot_metrics <- rbind(
        plot_metrics,
        data.frame(
          Proportion = rep(as.numeric(downsample_proportions[i])*100, 3),
          Metric = c("Sensitivity", "Specificity", "Correlation"),
          Score = c(
            accuracy_metrics$V2[accuracy_metrics$V1 == "sensitivity"],
            accuracy_metrics$V2[accuracy_metrics$V1 == "specificity"],
            1
          )
        )
      )
    }
    
  } else {
    in_dir <- paste0(results_dir, "infercnv/t_cells_included/", 
      sample_name, "/downsampling/", downsample_proportions[i], 
      "_downsampling/tables/")

  
    if (i==1) {
      accuracy_metrics <- read.table(
        paste0(in_dir, "accuracy_metrics.txt"),
        header = F, sep = "\t", as.is = T
      )
      print(paste0("Correlation p-value = ", 
                   accuracy_metrics$V2[accuracy_metrics$V1 == "pearson_p_val"]))
      plot_metrics <- data.frame(
        Proportion = rep(as.numeric(downsample_proportions[i])*100, 3),
        Metric = c("Sensitivity", "Specificity", "Correlation"),
        Score = c(
          accuracy_metrics$V2[accuracy_metrics$V1 == "sensitivity"],
          accuracy_metrics$V2[accuracy_metrics$V1 == "specificity"],
          accuracy_metrics$V2[accuracy_metrics$V1 == "pearson_R_squared"]
        )
      )
      
    } else {
      accuracy_metrics <- read.table(
        paste0(in_dir, "accuracy_metrics.txt"),
        header = F, sep = "\t", as.is = T
      )
      print(paste0("Correlation p-value = ", 
                   accuracy_metrics$V2[accuracy_metrics$V1 == "pearson_p_val"]))
      plot_metrics <- rbind(
        plot_metrics,
        data.frame(
          Proportion = rep(as.numeric(downsample_proportions[i])*100, 3),
          Metric = c("Sensitivity", "Specificity", "Correlation"),
          Score = c(
            accuracy_metrics$V2[accuracy_metrics$V1 == "sensitivity"],
            accuracy_metrics$V2[accuracy_metrics$V1 == "specificity"],
            accuracy_metrics$V2[accuracy_metrics$V1 == "pearson_R_squared"]
          )
        )
      )
    }
  }
}

cols <- c("#B488B4", "#7CBA61", "#9DC9DC")

p <- ggplot(plot_metrics, aes(x = Proportion, y = Score)) 
p <- p + geom_line(aes(color = Metric))
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
p <- p + xlab("Proportion UMIs (%)")
p

pdf(
  paste0(in_path, "UMI_downsampling_results.pdf"), 
  width = 7, height = 5
)
  print(p)
dev.off()

png(
  paste0(in_path, "UMI_downsampling_results.png"), 
  width = 7, height = 5, units = "in", res = 300
)
  print(p)
dev.off()


