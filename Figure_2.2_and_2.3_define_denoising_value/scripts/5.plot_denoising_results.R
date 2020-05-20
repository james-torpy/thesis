#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

project_name <- "thesis"
subproject_name <- "Figure_2.2_and_2.3_define_denoising_value"
args = commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
analysis_mode <- args[2]
sim_names <- args[3]
denoise_values <- args[4]
ignore_denoising <- as.numeric(args[5])

#project_name <- "thesis"
#subproject_name <- "Figure_2.2_and_2.3_define_denoising_value"
#args = commandArgs(trailingOnly=TRUE)
#sample_name <- "CID4520N"
#analysis_mode <- "samples"
#sim_names <- paste0(
#  "filtered_normal.",
#  paste0("sim", 1:30, collapse = ".")
#)
##sim_names <- "filtered_normal.sim1"
#denoise_values <- "no_0.5_1_1.1_1.2_1.3_1.4_1.5"
#ignore_denoising <- as.numeric("2")

# split multi-element variable vectors:
sim_names <- unlist(
  strsplit(
    sim_names,
    split = "\\."
  )
)
sim_names <- sim_names[grep("normal", sim_names, invert = T)]

denoise_values <- unlist(
  strsplit(
    denoise_values,
    split = "_"
  )
)
denoise_values <- denoise_values[!(denoise_values %in% ignore_denoising)]

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("Analysis mode = ", analysis_mode))
print("Simulation names = ")
print(sim_names)
print("Denoising values = ")
print(denoise_values)
print(paste0("Ignore desnoising value = ", ignore_denoising))

library(ggplot2)
theme_set(theme_minimal())
library(dplyr)
library(cowplot)
library(scales)
  
lib_loc <- "/share/ScratchGeneral/jamtor/R/3.5dev/"
home_dir <- "/share/ScratchGeneral/jamtor/"

project_dir <- paste0(home_dir, "projects/", 
  project_name, "/", subproject_name, "/")
results_dir <- seurat_path <- paste0(project_dir, "results/")
func_dir <- paste0(project_dir, "scripts/functions/")

in_path <- paste0(results_dir, "infercnv/", sample_name, "/")

out_dir <- paste0(in_path, "denoising_results/")
table_dir <- paste0(out_dir, "tables/")
system(paste0("mkdir -p ", table_dir))


#######################################################################################
### 1. Load, calculate stats for and plot normal data ###
#######################################################################################ż

# load each file of fd from mean CNV score and calculate scaled denoise score:
for (i in 1:length(denoise_values)) {
  if (i==1) {

    normal_input <- paste0(denoise_values[i], "_denoising")
    score_df <- read.table(
      paste0(in_path, "/normal/", normal_input, "/", analysis_mode, 
        "_mode/tables/mean_CNV_score_per_cell.txt"),
      sep = "\t",
      header = T,
      as.is = T
    )

    # plot distribution of fd from mean score:
    system(paste0("mkdir -p ", out_dir, normal_input))
    pdf(paste0(out_dir, normal_input, "/fd_from_mean_distribution.pdf"))
      hist(score_df$score)
    dev.off()

    normal_CNV_score <- data.frame(
      denoise_value = gsub("_denoising", "", normal_input),
      mean = round(mean(score_df$score), 6),
      SE = round(sd(score_df$score)/sqrt(nrow(score_df)), 6)
    )

  } else {

    normal_input <- paste0(denoise_values[i], "_denoising")
    score_df <- read.table(
      paste0(in_path, "/normal/", normal_input, "/", analysis_mode, 
        "_mode/tables/mean_CNV_score_per_cell.txt"),
      sep = "\t",
      header = T,
      as.is = T
    )

    # plot distribution of fd from mean score:
    system(paste0("mkdir -p ", out_dir, normal_input))
    pdf(paste0(out_dir, normal_input, "/fd_from_mean_distribution.pdf"))
      hist(score_df$score)
    dev.off()

    normal_CNV_score <- rbind(
      normal_CNV_score,
      data.frame(
        denoise_value = gsub("_denoising", "", normal_input),
        mean = round(mean(score_df$score), 6),
        SE = round(sd(score_df$score)/sqrt(nrow(score_df)), 6)
      )
    )

  }
}
normal_CNV_score$denoise_value <- factor(
  normal_CNV_score$denoise_value,
  levels = c(
    "no", 
    levels(normal_CNV_score$denoise_value)[levels(normal_CNV_score$denoise_value) != "no"]
  )
)

# plot fd from mean InferCNV score:
p <- ggplot(normal_CNV_score, aes(x = denoise_value, y = mean)) 
p <- p + geom_line(group = 1, color="#C02456")
p <- p + geom_errorbar(
  aes(ymin=mean-SE, ymax=mean+SE), 
  color="#C02456",
  width=.1
)
p <- p + theme(
  panel.grid.minor = element_blank(),
  axis.title.x = element_text(size = 9),
  axis.title.y = element_text(size = 9),
  legend.title = element_blank()
)
p <- p + xlab("Denoising range (std dev from mean)")
p <- p + ylab("InferCNV signal")

pdf(
  paste0(out_dir, "normal_infercnv_score_denoising_results.pdf"), 
  width = 7, height = 5
)
  print(p)
dev.off()

png(
  paste0(out_dir, "normal_infercnv_score_denoising_results.png"), 
  width = 7, height = 5, units = "in", res = 300
)
  print(p)
dev.off()


#######################################################################################
### 2. Load, calculate stats for and plot simulation accuracy data ###
#######################################################################################

for (k in 1:length(sim_names)) {
  if (k==1) {
    sim_input <- paste0(
      in_path, sim_names[k], "/", denoise_values, "_denoising"
    )
  } else {
    sim_input <- c(
      sim_input, 
      paste0(
        in_path, sim_names[k], "/", denoise_values, "_denoising"
      )
    )
  }
}

metrics_list <- list("sensitivity", "specificity", "precision", "F1")
for (i in 1:length(sim_input)) {
  if (i==1) {

    score_df <- read.table(
      paste0(sim_input[i], "/", analysis_mode, 
        "_mode/tables/accuracy_metrics.txt"),
      sep = "\t",
      header = T,
      as.is = T
    )
    denoise_value <- gsub(
      "_denoising", "", gsub(
      	"^.*/", "", sim_input[i]
      )
    )

    accuracy_score <- lapply(metrics_list, function(x) {
      return(
      	data.frame(
          denoise_value = gsub("_denoising", "", denoise_value),
          score = score_df[x,]
        )
      )
    })

  } else {

    score_df <- read.table(
      paste0(sim_input[i], "/", analysis_mode, 
        "_mode/tables/accuracy_metrics.txt"),
      sep = "\t",
      header = T,
      as.is = T
    )
    denoise_value <- gsub(
      "_denoising", "", gsub(
      	"^.*/", "", sim_input[i]
      )
    )

    temp_accuracy_score <- lapply(metrics_list, function(x) {
      return(
      	data.frame(
          denoise_value = gsub("_denoising", "", denoise_value),
          score = score_df[x,]
        )
      )
    })
   accuracy_score <- Map(rbind, accuracy_score, temp_accuracy_score) 
  }
}
names(accuracy_score) <- unlist(metrics_list)

# calculate mean and SEs for all metrics:
accuracy_means <- lapply(accuracy_score, function(x) {
  split_values <- split(x, x$denoise_value)
  for (v in 1:length(split_values)) {
  	if (v==1) {
  	  mean_df <- data.frame(
  	  	denoise_value = split_values[[v]]$denoise_value[1],
        mean = mean(split_values[[v]]$score),
        SE = sd(split_values[[v]]$score)/sqrt(nrow(split_values[[v]]))
      )
  	} else {
  	  mean_df <- rbind(
  	  	mean_df,
  	  	data.frame(
  	  	  denoise_value = split_values[[v]]$denoise_value[1],
          mean = mean(split_values[[v]]$score),
          SE = sd(split_values[[v]]$score)/sqrt(nrow(split_values[[v]]))
        )
      )
  	}
  }

  return(mean_df)
})

# add metric name to each df and bind them together:
for (d in 1:length(accuracy_means)) {
  accuracy_means[[d]]$metric = names(accuracy_means)[d]
}
accuracy_mean_df <- do.call("rbind", accuracy_means)
accuracy_mean_df$denoise_value <- factor(
  accuracy_mean_df$denoise_value,
  levels = c(
    "no", 
    levels(accuracy_mean_df$denoise_value)[levels(accuracy_mean_df$denoise_value) != "no"]
  )
)

accuracy_mean_df$metric <- factor(
  accuracy_mean_df$metric,
  levels = c("sensitivity", "specificity", "precision", "F1")
)

write.table(
  accuracy_mean_df,
  paste0(table_dir, "mean_accuracy_values.txt"),
  sep = "\t",
  quote = F,
  col.names = T,
  row.names = F
)

sim_cols <- c("#58B9DB", "#D95F02", "#F4D30B", "#B066B2")
p <- ggplot(
  accuracy_mean_df, 
  aes(
  	x = denoise_value, 
  	y = mean, 
  	group = metric, 
  	color = metric
  )
)
p <- p + geom_line()
p <- p + geom_errorbar(
  aes(ymin=mean-SE, ymax=mean+SE), 
  width=.1
)
p <- p + scale_color_manual(values=sim_cols)
pdf(
  paste0(out_dir, "sim_denoising_results.pdf"), 
  width = 7, height = 5
)
  print(p)
dev.off()

png(
  paste0(out_dir, "sim_denoising_results.png"), 
  width = 7, height = 5, unit = "in", res = 300
)
  print(p)
dev.off()


#######################################################################################
### 3. Plot normal CNV score and simulation accuracy scores together ###
#######################################################################################

# convert normal CNV score to fit the same scale as accuracy metrics:
normal_CNV_score$scaled_mean <- normal_CNV_score$mean*20
normal_CNV_score$scaled_SE <- normal_CNV_score$SE*20

# bind scaled normal CNV score and accuracy score together:
plot_norm_score <- subset(
  normal_CNV_score, select = c(denoise_value, scaled_mean, scaled_SE)
)
plot_norm_score$metric =  "CNV_score"
colnames(plot_norm_score) <- gsub("scaled_", "", colnames(plot_norm_score))

all_score <- rbind(plot_norm_score, accuracy_mean_df)
all_score$metric <- factor(
  all_score$metric,
  levels = c("CNV_score", "sensitivity", "specificity", "precision", "F1")
)

all_cols <- c("#C02456", "#58B9DB", "#D95F02", "#F4D30B", "#B066B2")
p <- ggplot(
  all_score, 
  aes(
  	x = denoise_value, 
  	y = mean, 
  	group = metric, 
  	color = metric
  )
)
p <- p + geom_line()
p <- p + geom_errorbar(
  aes(ymin=mean-SE, ymax=mean+SE), 
  width=.1
)
p <- p + scale_y_continuous(sec.axis = sec_axis(~ . / 20, name = "Normal breast noise"))
p <- p + scale_color_manual(values=all_cols)
p <- p + xlab("Denoising range (SD from mean)")
p <- p + ylab("Cancer simulation accuracy score")

pdf(
  paste0(out_dir, "all_denoising_results.pdf"), 
  width = 7, height = 5
)
  print(p)
dev.off()

png(
  paste0(out_dir, "all_denoising_results.png"), 
  width = 7, height = 5, unit = "in", res = 300
)
  print(p)
dev.off()


#######################################################################################
### 4. Plot normal CNV score, sensitivity and specificity ###
#######################################################################################

plot_scores <- all_score[
  all_score$metric %in% c("CNV_score", "sensitivity", "specificity"),
]

cols_2 <- c("#C02456", "#58B9DB", "#D95F02")
p <- ggplot(
  plot_scores, 
  aes(
  	x = denoise_value, 
  	y = mean, 
  	group = metric, 
  	color = metric
  )
)
p <- p + geom_line()
p <- p + geom_errorbar(
  aes(ymin=mean-SE, ymax=mean+SE), 
  width=.1
)
p <- p + scale_y_continuous(sec.axis = sec_axis(~ . / 20, name = "Normal breast noise"))
p <- p + scale_color_manual(values=all_cols, 
  labels=c("CNV score", "Sensitivity", "Specificity"))
p <- p + xlab("Denoising range (SD from mean)")
p <- p + ylab("Cancer simulation accuracy score")

pdf(
  paste0(out_dir, "CNV_sensitivity_specificity_denoising_results.pdf"), 
  width = 7, height = 5
)
  print(p)
dev.off()

png(
  paste0(out_dir, "CNV_sensitivity_specificity_denoising_results.png"), 
  width = 7, height = 5, unit = "in", res = 300
)
  print(p)
dev.off()


#######################################################################################
### 5. Plot normal CNV score, recall and precision ###
#######################################################################################

plot_scores <- all_score[
  all_score$metric %in% c("CNV_score", "sensitivity", "precision"),
]

cols_2 <- c("#C02456", "#58B9DB", "#F4D30B")
p <- ggplot(
  plot_scores, 
  aes(
  	x = denoise_value, 
  	y = mean, 
  	group = metric, 
  	color = metric
  )
)
p <- p + geom_line()
p <- p + geom_errorbar(
  aes(ymin=mean-SE, ymax=mean+SE), 
  width=.1
)
p <- p + scale_y_continuous(sec.axis = sec_axis(~ . / 20, name = "Normal breast noise"))
p <- p + scale_color_manual(values=cols_2, labels=c("CNV score", "Recall", "Precision"))
p <- p + xlab("Denoising range (SD from mean)")
p <- p + ylab("Cancer simulation accuracy score")

pdf(
  paste0(out_dir, "CNV_recall_precision_denoising_results.pdf"), 
  width = 7, height = 5
)
  print(p)
dev.off()

png(
  paste0(out_dir, "CNV_recall_precision_denoising_results.png"), 
  width = 7, height = 5, unit = "in", res = 300
)
  print(p)
dev.off()

