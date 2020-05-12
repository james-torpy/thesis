#! /share/ClusterShare/software/contrib/briglo/octR/src/R-3.6.0/builddir/bin/RRscript

project_name <- "thesis"
subproject_name <- "Figure_2.2_and_2.3_define_denoising_value"
args = commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
analysis_mode <- args[2]
sim_names <- args[3]

project_name <- "thesis"
subproject_name <- "Figure_2.2_and_2.3_define_denoising_value"
args = commandArgs(trailingOnly=TRUE)
sample_name <- "CID4520N"
analysis_mode <- "samples"
sim_names <- paste0("sim", 1:30, collapse = ".")

# split multi-element variable vectors:
sim_names <- unlist(
  strsplit(
    sim_names,
    split = "\\."
  )
)
sim_names <- sim_names[grep("normal", sim_names, invert = T)]

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("Analysis mode = ", analysis_mode))
print("Simulation names = ")
print(sim_names)

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library("rlang", lib.loc = lib_loc)
library("dplyr", lib.loc = lib_loc)
library("ggpubr", lib.loc = lib_loc)
#library(ggplot2)
theme_set(theme_minimal())
library(dplyr)
library(cowplot)
library(scales)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", 
  project_name, "/", subproject_name, "/")
results_dir <- seurat_path <- paste0(project_dir, "results/")
in_path <- paste0(results_dir, "infercnv/", sample_name, "/")

out_path <- paste0(in_path, "distinguish_between_CNV_results/")
plot_dir <- paste0(out_path, "plots/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(out_path, "tables/")
system(paste0("mkdir -p ", table_dir))


#################################################################################
#### 1. Load data for each simulation ###
#################################################################################

for (i in 1:length(sim_names)) {

  if (i==1) {

    CNV_data <- list(
      readRDS(
        paste0(in_path, sim_names[i], "/no_denoising/", analysis_mode, 
          "_mode/Rdata/CNV_data.Rdata")
      )
    )

  } else {

    CNV_data[[i]] <- readRDS(
      paste0(in_path, sim_names[i], "/no_denoising/", analysis_mode, 
        "_mode/Rdata/CNV_data.Rdata")
    )

  }

}


#################################################################################
#### 2. Format and combine data ###
#################################################################################

#i=1
signal_dfs <- lapply(CNV_data, function(x) {

  # split CNV_indices by multiplier:
  split_CNV_indices <- split(x$CNV_indices, x$CNV_indices$multiplier)
  
  # boxplot all average signal values for each copy number:
  # split CNV_indices by multiplier:
  split_CNV_indices <- split(x$CNV_indices, x$CNV_indices$multiplier)
  split_CNV_indices <- split_CNV_indices[names(split_CNV_indices) != 1]
  # calculate average infercnv signal for each multiplier:
  signal_per_multiplier <- lapply(split_CNV_indices, function(y) {
    # fetch true positive signal values for all regions:
    for (r in 1:nrow(y)) {
      if (r==1) {
        average_signal <- x$average_signal[(y$start[r]):(y$end[r])]
        acc_call <- x$accuracy_annotation_vector[(y$start[r]):(y$end[r])]
        true_positive_signal <- mean(average_signal[acc_call == "true_positive"])
        
#        if (any(is.na(true_positive_signal))) {
#          print(paste0("NAs present in sim ", i, ", copy number ",
#          	unique(y$multiplier), ", row ", r))
#        }

      } else {
        average_signal <- x$average_signal[(y$start[r]):(y$end[r])]
        acc_call <- x$accuracy_annotation_vector[(y$start[r]):(y$end[r])]
        if (any(acc_call != "wrong_call")) {
          true_positive_signal <- c(
            true_positive_signal,
            mean(average_signal[acc_call == "true_positive"])
          )
        }
#        if (any(is.na(true_positive_signal))) {
#          print(paste0("NAs present in sim ", i, ", copy number ",
#          	unique(y$multiplier), ", row ", r))
#        }

      }
    }
    
    return(data.frame(signal = true_positive_signal))
  })
  # label copy number in each df:
  for (j in 1:length(signal_per_multiplier)) {
    signal_per_multiplier[[j]]$copy_number <- names(signal_per_multiplier)[j]
  }
  
  # convert to one df:
  signal_df <- do.call("rbind", signal_per_multiplier)

  i <<- i+1
  return(signal_df)

})

# label each df with sim number:
for (l in 1:length(signal_dfs)) {
  signal_dfs[[l]]$sim = l
}

# bind together dfs:
all_signal <- do.call("rbind", signal_dfs)
all_signal <- all_signal[order(all_signal$copy_number),]

# remove NAs:
all_signal <- all_signal[!is.na(all_signal$signal),]

# split into signal per multiplier:
signal_per_multiplier <- split(all_signal, all_signal$copy_number)

# plot distibutions of each copy number:
for (k in 1:length(signal_per_multiplier)) {
  print(k)
  png(
    paste0(
    plot_dir, "copy_number_", names(signal_per_multiplier)[k], 
    "_distribution.png"
    )
  )
    plot(density(signal_per_multiplier[[k]]$signal))
  dev.off()
}

# check n for each group:
group_lengths <- unlist(lapply(signal_per_multiplier, nrow))
write.table(
  data.frame(group_lengths), 
  paste0(table_dir, "/copy_number_group_lengths.txt"),
  quote = F,
  row.names = T,
  col.names = F
)

# separate gains and losses:
gain_df <- all_signal[as.numeric(all_signal$copy_number) > 1,]
loss_df <- all_signal[as.numeric(all_signal$copy_number) < 1,]
loss_df$copy_number <- factor(
  loss_df$copy_number,
  levels = c("0.5", "0")
)


#################################################################################
#### 3. Check for significant differences ###
#################################################################################

# use pairwise wilcox to test whether loss groups are significantly different,
# using Benjamini & Hochberg method to correct for multiple testing (FDR guys):
gain_df$copy_number <- factor(gain_df$copy_number, levels = c("1.5", "2", "3"))
gain_wilcox <- pairwise.wilcox.test(gain_df$signal, gain_df$copy_number,
  alternative = "greater", p.adjust.method = "BH")

# create gain boxplot:
my_comparisons = list( c("1.5", "2"), c("2", "3") )

p <- ggboxplot(gain_df, x = "copy_number", y = "signal",
  fill = "copy_number", palette = c("#BF889F", "#BD5D89", "#C03667"))
p <- p + xlab("Copy number fold change")
p <- p + ylab("CNV signal")
p <- p + theme(
  legend.position = "none"
)
p <- p + stat_compare_means(
  comparisons = my_comparisons,
  method = "wilcox.test",
  label = "p.signif",
  label.y = c(
    max(
      c(signal_per_multiplier[[3]]$signal, signal_per_multiplier[[4]]$signal)
    ) + 0.01,
    max(
      c(signal_per_multiplier[[4]]$signal, signal_per_multiplier[[5]]$signal)
    ) + 0.012
  )
)
p <- p + theme(
  axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
  axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
)
png(
  paste0(plot_dir, "signal_vs_gain_copy_number_boxplot.png"),
  height = 5,
  width = 7,
  res = 300,
  units = "in"
)   
  p 
dev.off()

# use wilcox to test whether loss groups are significantly different:
loss_df$copy_number <- factor(gain_df$copy_number, levels = c("0.5", "0"))
loss_wilcox <- wilcox.test(
  signal_per_multiplier[[2]]$signal, 
  signal_per_multiplier[[1]]$signal, 
  alternative = "greater"
)
loss_p_val <- loss_wilcox$p.value

# create loss boxplot:
my_comparisons = list( c("0.5", "0") )

p <- ggboxplot(loss_df, x = "copy_number", y = "signal",
  fill = "copy_number", palette = c("#82BFCE", "#618AC7"))
p <- p + xlab("Copy number fold change")
p <- p + ylab("CNV signal")
p <- p + theme(
  legend.position = "none"
)
p <- p + stat_compare_means(
  comparisons = my_comparisons,
  method = "wilcox.test",
  label = "p.signif",
  label.y = max(loss_df$signal)+0.01
)
p <- p + theme(
  axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
  axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
)
png(
  paste0(plot_dir, "signal_vs_loss_copy_number_boxplot.png"),
  height = 5,
  width = 7,
  res = 300,
  units = "in"
)   
  p 
dev.off()

