#! /share/ClusterShare/software/contrib/briglo/octR/src/R-3.6.0/builddir/bin/RRscript

project_name <- "thesis"
subproject_name <- "Figure_1.2_and_1.3_simulate_cancer_and_define_denoising_value"
args = commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
analysis_mode <- args[2]
sim_names <- args[3]

#project_name <- "thesis"
#subproject_name <- "Figure_1.2_and_1.3_simulate_cancer_and_define_denoising_value"
#args = commandArgs(trailingOnly=TRUE)
#sample_name <- "CID4520N"
#analysis_mode <- "samples"
#sim_names <- paste0("sim", 1:30, collapse = ".")

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
Robject_dir <- paste0(out_path, "Rdata/")
system(paste0("mkdir -p ", Robject_dir))


#################################################################################
#### 1. Load data for each simulation ###
#################################################################################

if (!file.exists(paste0(Robject_dir, "/1.simulation_CNV_data.Rdata"))) {
  for (i in 1:length(sim_names)) {

    if (i==1) {
  
      CNV_data <- list(
        readRDS(
          paste0(in_path, sim_names[i], "/1.3_denoising/", analysis_mode, 
            "_mode/Rdata/CNV_data.Rdata")
        )
      )
  
    } else {
  
      CNV_data[[i]] <- readRDS(
        paste0(in_path, sim_names[i], "/1.3_denoising/", analysis_mode, 
          "_mode/Rdata/CNV_data.Rdata")
      )
  
    }
  
  }

  saveRDS(CNV_data, paste0(Robject_dir, "/1.simulation_CNV_data.Rdata"))

} else {
  CNV_data <- readRDS(paste0(Robject_dir, "/1.simulation_CNV_data.Rdata"))
}


#################################################################################
#### 2. Format and combine data ###
#################################################################################

signal_dfs <- lapply(CNV_data, function(x) {

  # split CNV_indices by multiplier:
  split_CNV_indices <- split(x$CNV_indices, x$CNV_indices$multiplier)
  
  # boxplot all average signal values for each copy number:
  # split CNV_indices by multiplier:
  split_CNV_indices <- split(x$CNV_indices, x$CNV_indices$multiplier)
  split_CNV_indices <- split_CNV_indices[names(split_CNV_indices) != 1]
  # calculate average infercnv signal for each multiplier:
  multiplier_signal <- lapply(split_CNV_indices, function(y) {
    # fetch true positive signal values for all regions and save both 
    # individual values and mean:
    for (r in 1:nrow(y)) {
      if (r==1) {

        average_signal <- x$average_signal[(y$start[r]):(y$end[r])]
        acc_call <- x$accuracy_annotation_vector[(y$start[r]):(y$end[r])]
        true_positive_signal <- average_signal[acc_call == "true_positive"]

      } else {

        average_signal <- x$average_signal[(y$start[r]):(y$end[r])]
        acc_call <- x$accuracy_annotation_vector[(y$start[r]):(y$end[r])]
        if (any(acc_call != "wrong_call")) {
          true_positive_signal <- c(
            true_positive_signal,
            average_signal[acc_call == "true_positive"]
          )
        }

      }
    }
    return(data.frame(signal = true_positive_signal))
  })
  
  # label copy number in each df:
  for (j in 1:length(multiplier_signal)) {
    multiplier_signal[[j]]$copy_number <- names(multiplier_signal)[j]
  }
  
  # convert to one df:
  signal_df <- do.call("rbind", multiplier_signal)

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

# plot distributions of signal values for each multiplier:
for (t in 1:length(signal_per_multiplier)) {
  distr_plot <- density(signal_per_multiplier[[t]]$signal)
  png(
    paste0(
      plot_dir, "distribution_of_", 
      names(signal_per_multiplier)[t], "_copy_number_signal.png"
    )
  )
    plot(distr_plot)
  dev.off()
}

# organise signal values into gain and loss dfs:
gain_df <- all_signal[all_signal$signal > 0,]
loss_df <- all_signal[all_signal$signal < 0,]

# save all signal scores:
saveRDS(signal_per_multiplier, paste0(
  Robject_dir, "2.all_true_positive_CNV_signal_per_multiplier.Rdata")
)


#################################################################################
#### 3. Check for significant differences ###
#################################################################################

# use pairwise wilcox to test whether loss groups are significantly different,
# using Benjamini & Hochberg method to correct for multiple testing (FDR guys):
gain_df$copy_number <- factor(gain_df$copy_number, levels = c("1.5", "2", "3"))
gain_wilcox <- pairwise.wilcox.test(gain_df$signal, gain_df$copy_number,
  alternative = "greater", p.adjust.method = "BH")

# write gain signifcance values as table:
write.table(
  as.data.frame(gain_wilcox$p.value),
  paste0(table_dir, "gain_wilcox_p_vals.txt"),
  sep = "\t",
  quote = F,
  col.names = T,
  row.names = T
)

# use wilcox to test whether loss groups are significantly different:
loss_df$copy_number <- factor(loss_df$copy_number, levels = c("0.5", "0"))
loss_wilcox <- wilcox.test(
  loss_df$signal[loss_df$copy_number == 0.5],
  loss_df$signal[loss_df$copy_number == 0], 
  alternative = "greater"
)
loss_p_val <- loss_wilcox$p.value

# write loss signifcance values as table:
write.table(
  loss_wilcox$p.value,
  paste0(table_dir, "loss_wilcox_p_vals.txt"),
  sep = "\t",
  quote = F,
  col.names = F,
  row.names = F
)

# save n for each copy number:
write.table(
  as.data.frame(
    unlist(
      lapply(signal_per_multiplier, nrow)
    )
  ),
  paste0(table_dir, "signal_value_number_per_CNV_type.txt"),
  sep = "\t",
  quote = F,
  col.names = F,
  row.names = T
)

# bind loss and gain dfs together:
both_df <- rbind(loss_df, gain_df)

# create copy number distribution boxplot:
my_comparisons = list( c("0", "0.5"), c("1.5", "2"), c("2", "3") )

p <- ggboxplot(
  both_df, 
  x = "copy_number", 
  y = "signal",
  fill = "copy_number", 
  palette = c("#618AC7", "#82BFCE", "#BF889F", "#BD5D89", "#C03667"),
  bxp.errorbar = TRUE
)
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
  	max(loss_df$signal)+0.01,
    max(
      c(
        gain_df$signal[gain_df$copy_number == 1.5], 
        gain_df$signal[gain_df$copy_number == 2])
      ) + 0.01,
    max(
      c(
        gain_df$signal[gain_df$copy_number == 2], 
        gain_df$signal[gain_df$copy_number == 3])
      ) + 0.014
  )
)
p <- p + theme(
  axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
  axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
)
png(
  paste0(plot_dir, "signal_vs_copy_number_boxplot.png"),
  height = 5,
  width = 7,
  res = 300,
  units = "in"
)   
  p 
dev.off()


#################################################################################
#### 4. Fetch counts of correctly estimated CNV peaks across all simulations ###
#################################################################################

for (i in 1:length(sim_names)) {

  estimated_vs_known_data <- readRDS(
    paste0(in_path, sim_names[i], "/1.3_denoising/", analysis_mode, 
      "_mode/Rdata/estimated_vs_known_counts.Rdata")
  )
  correct_counts <- lapply(estimated_vs_known_data, function(x) {
    split_df <- split(x$estimated_vs_known, x$estimated_vs_known$copy_no)
    temp_counts <- lapply(split_df, function(y) {
      return(
        data.frame(
          copy_no = y$copy_no[1],
          correct = length(which(y$correct)),
          total = nrow(y)
        )
      )
    })
    return(do.call("rbind", temp_counts))
  })
  correct_counts <- do.call("rbind", correct_counts)
  correct_counts$sim <- sim_names[i]
  if (i==1) {
    all_counts <- correct_counts
  } else {
    all_counts <- rbind(all_counts, correct_counts)
  }
 
}

split_counts <- split(all_counts, all_counts$copy_no)
proportion_correct <- lapply(split_counts, function(x) {

  x$proportion <- round(x$correct/x$total, 2)*100

  return(
    data.frame(
      mean = round(mean(x$proportion), 1),
      SE = round(sd(x$proportion)/sqrt(nrow(x)), 1)
    )
  )

})
plot_df <- do.call("rbind", proportion_correct)

plot_df$type <- "gain"
plot_df$type[as.numeric(rownames(plot_df)) < 1] <- "loss"

plot_df$copy_no <- factor(
  rownames(plot_df),
  levels = sort(rownames(plot_df))
)

# write as table:
write.table(
  plot_df, 
  paste0(table_dir, "percentage_peaks_called_per_copy_number.txt"),
  quote = F,
  sep = "\t",
  row.names = T,
  col.names = F
)

# plot as barplot:
p <- ggplot(plot_df, aes(x=copy_no, y = mean, fill = type))
p <- p + geom_bar(stat="identity", width = 0.75)
p <- p + geom_errorbar(
  aes(ymin=mean-SE, ymax=mean+SE), 
  width=.1
)
p <- p + scale_y_continuous(
  expand = c(0, 0),
  limits = c(0, 100),
  breaks = seq(0, 100, 25),
  labels = seq(0, 100, 25)
)
p <- p + theme_cowplot(12)
p <- p + xlab("Copy number")
p <- p + ylab("Proportion correctly called (%)")
p <- p + scale_fill_manual(values = c("#BF3667", "#58B9DB"))
p <- p + theme(
  axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
  axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)),
  legend.position = "none"
)
png(
  paste0(plot_dir, "correctly_called_peaks.png"), 
  width = 5, height = 5, unit = "in", res = 300
)
  print(p)
dev.off()







