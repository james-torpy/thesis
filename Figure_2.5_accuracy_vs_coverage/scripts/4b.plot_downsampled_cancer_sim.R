#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

### add F1 scores to plots ###

RStudio <- FALSE

project_name <- "thesis"
subproject_name <- "Figure_2.5_accuracy_vs_coverage"
args = commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
print(sample_name)
include_t_cells <- as.logical(args[2])
print(include_t_cells)
downsample_proportions <- unlist(strsplit(args[3], split = "_"))
print(downsample_proportions)
downsample_type <- args[4]
print(downsample_type)
CNV_type <- args[5]
print(CNV_type)
simulation_number_start <- as.numeric(args[6])
print(simulation_number_start)
simulation_number_end <- as.numeric(args[7])
print(simulation_number_end)
analysis_mode <- args[8]
print(analysis_mode)

sample_name <- "CID4520N_cancer_sim"
include_t_cells <- TRUE
downsample_proportions <- unlist(
  strsplit(
    "0.05_0.1_0.15_0.2_0.3_0.4_0.5_0.6_0.7_0.8_0.9_no", 
    split = "_"
  )
)
downsample_type <- "gene"
CNV_type <- "both"
simulation_number_start <- as.numeric("1")
simulation_number_end <- as.numeric("30")
analysis_mode <- "samples"

simulation_numbers <- simulation_number_start:simulation_number_end

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("T-cells included? ", as.character(include_t_cells)))

lib_loc = "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(ggplot2)
theme_set(theme_minimal())
library(dplyr)
library(cowplot)
library(inflection, lib.loc = lib_loc)

if (RStudio) {
  
  home_dir <- "/Users/jamestorpy/clusterHome/"
  
} else {
  
  lib_loc <- "/share/ScratchGeneral/jamtor/R/3.5dev/"
  
  home_dir <- "/share/ScratchGeneral/jamtor/"
  
}

project_dir <- paste0(home_dir, "projects/", 
  project_name, "/", subproject_name, "/")
results_dir <- seurat_path <- paste0(project_dir, "results/")
func_dir <- paste0(project_dir, "scripts/functions/")

if (include_t_cells) {
  in_path <- paste0(results_dir, "infercnv/t_cells_included/", 
    sample_name, "/", CNV_type, "/")
} else {
  in_path <- paste0(results_dir, "infercnv/t_cells_excluded/", 
    sample_name, "/", CNV_type, "/")
}

sim_path <- paste0(results_dir, "cancer_simulation/", sample_name,
  "/", CNV_type, "/")

out_dir <- paste0(in_path, downsample_type, "_downsampling_final_results/")
table_dir <- paste0(out_dir, "tables/")
system(paste0("mkdir -p ", table_dir))
plot_dir <- paste0(out_dir, "plots/")
system(paste0("mkdir -p ", plot_dir))


################################################################################
### 1. Fetch data ###
################################################################################

if (!file.exists(paste0(table_dir, "final_metric_df.txt"))) {

  # convert downsample_proportions to list:
  downsample_list <- as.list(downsample_proportions)
  
  for (i in 1:length(simulation_numbers)) {
  
    metrics_list <- lapply(downsample_list, function(x) {
  
      print(paste0("Fetching metrics for simulation ", 
        simulation_numbers[i], " downsampled to ", x, "..."))
  
      if (downsample_type == "UMI") {
        in_dir <- paste0(in_path, "/", simulation_numbers[i], "/", 
          x, "_downsampling/", analysis_mode, "_mode/tables/")
      } else if (downsample_type == "gene") {
        if (x != "no") {
          in_dir <- paste0(in_path, "/", simulation_numbers[i], "/", 
            x, "_gene_downsampling/", analysis_mode, "_mode/tables/")
        } else if (x == "no") {
          in_dir <- paste0(in_path, "/", simulation_numbers[i], "/", 
            x, "_downsampling/", analysis_mode, "_mode/tables/")
        }
      }
      
      accuracy_metrics <- read.table(
        paste0(in_dir, "accuracy_metrics.txt"),
        sep = "\t",
        header = F,
        stringsAsFactors = F,
        fill = T
      )
      accuracy_metrics <- data.frame(
        row.names = accuracy_metrics$V1,
        score = accuracy_metrics$V2
      )
      
      if (x=="no") {
  
        metric_df <- data.frame(
          Proportion = rep(100, 4),
          Metric = c("Sensitivity", "Specificity", "Precision", "F1 score"),
          Score = c(
            accuracy_metrics["sensitivity",],
            accuracy_metrics["specificity",],
            accuracy_metrics["precision",],
            accuracy_metrics["F1",]
          ),
          Sim_no = rep(simulation_numbers[i], 4)
        )
  
        metric_df$Score[nrow(metric_df)] <- 1
  
      } else {
  
        metric_df <- data.frame(
          Proportion = rep(as.numeric(x)*100, 4),
          Metric = c("Sensitivity", "Specificity", "Precision", "F1 score"),
          Score = c(
            accuracy_metrics["sensitivity",],
            accuracy_metrics["specificity",],
            accuracy_metrics["precision",],
            accuracy_metrics["F1",]
          ),
          Sim_no = rep(simulation_numbers[i], 4)
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

  write.table(
  	final_metric_df,
  	paste0(table_dir, "final_metric_df.txt"),
  	sep = "\t",
  	quote = F,
  	col.names = T
  )

} else {

  final_metric_df <- read.table(
  	paste0(table_dir, "final_metric_df.txt"),
    sep = "\t",
    header = T
  )

}


################################################################################
### 2. Fetch mean nUMI and nGene per proportion ###
################################################################################

for (k in 1:length(simulation_numbers)) {

  temp_coverage <- read.table(
    paste0(sim_path, simulation_numbers[k], "/tables/mean_coverage_values.txt")
  )

  if (k==1) {
    coverage_df <- data.frame(
      simulation_number = k,
      downsample_proportion = temp_coverage$downsample_proportion,
      mean_nUMI = temp_coverage$nUMI,
      mean_nGene = temp_coverage$nGene
    )
  } else {
    coverage_df <- rbind(
      coverage_df,
      data.frame(
        simulation_number = k,
        downsample_proportion = temp_coverage$downsample_proportion,
        mean_nUMI = temp_coverage$nUMI,
        mean_nGene = temp_coverage$nGene
      )
    )
  }

}

# split coverage_df by downsample proportion:
split_coverage <- split(coverage_df, coverage_df$downsample_proportion)
proportion_coverages <- lapply(split_coverage, function(x) {
  return(
    data.frame(
      nUMI = round(mean(x$mean_nUMI), 0),
      nGene = round(mean(x$mean_nGene), 0)
    )
  )
})
final_coverages <- do.call("rbind", proportion_coverages)

if (downsample_type == "UMI") {
  final_coverages <- final_coverages[
    grep("gene", rownames(final_coverages), invert=T),
  ]
  final_coverages$Proportion <- as.numeric(
    gsub("no", "1", rownames(final_coverages))
  )*100
} else {
  final_coverages <- final_coverages[
    grep("gene|no", rownames(final_coverages)),
  ]
  final_coverages$Proportion <- as.numeric(
    gsub(
      "_gene", "", gsub(
        "no", "1", rownames(final_coverages)
      )
    )
  )*100
}
final_metric_df <- merge(final_metric_df, final_coverages, by="Proportion")


################################################################################
### 3. Plot sensitivity and specificity only ###
################################################################################

subset_df <- final_metric_df[
  final_metric_df$Metric %in% c("Sensitivity", "Specificity"),
]
subset_df$Metric <- factor(
  subset_df$Metric, 
  levels = c("Sensitivity", "Specificity")
)
cols <- c("#58B9DB", "#D95F02")

if (downsample_type == "UMI") {
  p <- ggplot(
    subset_df, 
    aes(x = nUMI, y = Score, group = Metric, color = Metric)
  ) 
  p <- p + xlab("Mean UMI/cell")
  p <- p + geom_errorbar(
    aes(ymin=Score-SE, ymax=Score+SE), 
    width=floor(max(subset_df$nUMI)/50)
  )
  p <- p + scale_x_continuous(
    breaks = seq(0, 14000, 2000),
    labels = seq(0, 14000, 2000)
  )
} else {
  p <- ggplot(
    subset_df, 
    aes(x = nGene, y = Score, group = Metric, color = Metric)
  )
  p <- p + xlab("Mean gene number/cell")
  p <- p + geom_errorbar(
    aes(ymin=Score-SE, ymax=Score+SE), 
    width=floor(max(subset_df$nGene)/50)
  )
  p <- p + scale_x_continuous(
    breaks = seq(0, 2500, 500),
    labels = seq(0, 2500, 500)
  )
}
p <- p + ylab("Accuracy score")
p <- p + geom_line()
p <- p + scale_color_manual(values = cols)
p <- p + theme_cowplot(12)
p <- p + theme(
  axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
  axis.title.y.right = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
  
)
pdf(
  paste0(plot_dir, "downsampling_sensitivity_specificity.pdf"), 
  width = 7, height = 5
)
  print(p)
dev.off()

png(
  paste0(plot_dir, "downsampling_sensitivity_specificity.png"), 
  width = 7, height = 5, units = "in", res = 300
)
  print(p)
dev.off()

print(paste0("All output plots in", out_dir))


################################################################################
### 4. Plot all metrics ###
################################################################################

cols <- c("#C02456", "#58B9DB", "#F4D30B", "#B066B2")

if (downsample_type == "UMI") {
  p <- ggplot(final_metric_df, aes(x = nUMI, y = Score)) 
  p <- p + xlab("Mean UMI/cell")
} else {
  p <- ggplot(final_metric_df, aes(x = nGene, y = Score)) 
  p <- p + xlab("Mean gene number/cell")
}
p <- p + geom_line(aes(color = Metric))
p <- p + geom_errorbar(aes(ymin=Score-SE, ymax=Score+SE))
p <- p + scale_color_manual(values = cols[1:5])
p <- p + scale_x_continuous(
  breaks = seq(0, 14000, 2000)
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
p <- p + theme_cowplot(12)
p <- p + theme(
  panel.grid.minor = element_blank(),
  axis.title.x = element_text(size = 9),
  axis.title.y = element_blank(),
  legend.title = element_blank()
)

pdf(
  paste0(plot_dir, "all_downsampling_results.pdf"), 
  width = 7, height = 5
)
  print(p)
dev.off()

png(
  paste0(plot_dir, "all_downsampling_results.png"), 
  width = 7, height = 5, units = "in", res = 300
)
  print(p)
dev.off()



################################################################################
### 5. Find derivative and inflection point of senstivity curve ###
################################################################################

# isolate sensitivity curve:
sensitivity_df <- final_metric_df[final_metric_df$Metric == "Sensitivity",]

# find first inflection point by calculating the bisection extremum distance 
# estimator:
if (downsample_type == "UMI") {
  temp_bede=bede(sensitivity_df$nUMI, sensitivity_df$Score, 0)
} else {
  temp_bede=bede(sensitivity_df$nGene, sensitivity_df$Score, 0)
}

# plot sensitivity with inflection point:
cols <- c("#58B9DB")

if (downsample_type == "UMI") {
  p <- ggplot(
    sensitivity_df, 
    aes(x = nUMI, y = Score, group = Metric, color = Metric)
  ) 
  p <- p + xlab("Mean UMI/cell")
  p <- p + geom_errorbar(
    aes(ymin=Score-SE, ymax=Score+SE), 
    width=floor(max(sensitivity_df$nUMI)/50)
  )
  p <- p + scale_x_continuous(
    breaks = seq(0, 14000, 2000),
    labels = seq(0, 14000, 2000)
  )
} else {
  p <- ggplot(
    sensitivity_df, 
    aes(x = nGene, y = Score, group = Metric, color = Metric)
  )
  p <- p + xlab("Mean gene number/cell")
  p <- p + geom_errorbar(
    aes(ymin=Score-SE, ymax=Score+SE), 
    width=floor(max(sensitivity_df$nGene)/50)
  )
  p <- p + scale_x_continuous(
    breaks = seq(0, 2500, 500),
    labels = seq(0, 2500, 500)
  )
}
p <- p + ylab("Accuracy score")
p <- p + geom_line()
p <- p + geom_vline(xintercept=temp_bede$iplast, colour = "blue")
p <- p + scale_color_manual(values = cols)
p <- p + theme_cowplot(12)
p <- p + theme(
  axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
  axis.title.y.right = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
  
)
pdf(
  paste0(plot_dir, "sensitivity_inflection_point.pdf"), 
  width = 7, height = 5
)
  print(p)
dev.off()

png(
  paste0(plot_dir, "specificity_inflection_point.png"), 
  width = 7, height = 5, units = "in", res = 300
)
  print(p)
dev.off()

print(paste0("All output plots in", plot_dir))
