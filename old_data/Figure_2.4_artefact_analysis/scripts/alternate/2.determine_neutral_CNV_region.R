
args = commandArgs(trailingOnly=TRUE)

project_name <- "thesis"
subproject_name <- "Figure_2.4_artefact_analysis"

project_name <- "thesis"
subproject_name <- "Figure_2.4_artefact_analysis"

print(paste0("Project name = ", project_name))
print(paste0("Subproject name = ", subproject_name))

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", 
  project_name, "/", subproject_name, "/")
results_dir <- paste0(project_dir, "results/")

in_dir <- paste0(project_dir, "raw_files/infercnv/")
out_dir <- paste0(results_dir, "/determine_neutral_regions/") 

system(paste0("mkdir -p ", in_dir))
system(paste0("mkdir -p ", out_dir))

print(paste0("In directory = ", in_dir))
print(paste0("Out directory = ", out_dir))

print("Determining CNV neutral regions...")

library(ggplot2)
library(cowplot)


################################################################################
### 1. Load infercnv output for each cancer sample and plot distribution of
# all CNV scores  ###
################################################################################

# fetch cancer sample names:
sample_files <- list.files(
  	in_dir, 
  	pattern = "infercnv.observations.txt", 
  	full.names = T,
  	recursive = T
  )
# load each infercnv output:
for (i in 1:length(sample_files)) {
  print(paste0("Loading ", sample_files[i]))
  if (i==1) {
    all_infercnv <- unlist(
      read.table(
        x, 
        as.is=T, 
        header = T
      )
    )
  } else {
    c(
      all_infercnv, 
      unlist(
        read.table(
          x, 
          as.is=T, 
          header = T
        )
      )
    )
  }
}


################################################################################
### 2. Define neutral region as 2 standard deviations either side of mean of
# main peak ###
################################################################################

# create density plot and determine rough netural window:
all_df <- data.frame(score = all_infercnv)
p <- ggplot(all_df, aes(score)) + 
  geom_density() + 
  geom_vline(xintercept = 0.95, colour = "red") + 
  geom_vline(xintercept = 1.05, colour = "red")
pdf(paste0(out_dir, "all_cancer_CNV_score_density_plot.pdf"))
  p
dev.off()
png(paste0(out_dir, "all_cancer_CNV_score_density_plot.png"))
  p
dev.off()

# isolate neutral peak and surrounds:
rough_neutral_df <- data.frame(
  score = all_df[all_df$score > 0.95 & all_df$score < 1.05,]
)

# plot ad mark 1 and 2 std devs from mean:
p <- ggplot(rough_neutral_df, aes(score)) + 
  geom_density() 
#  + 
#  geom_vline(xintercept = 0.955, colour = "red") + 
#  geom_vline(xintercept = 1.045, colour = "red")
pdf(paste0(out_dir, "all_cancer_CNV_neutral_score_density_plot.pdf"))
  p
dev.off()
png(paste0(out_dir, "all_cancer_CNV_neutral_score_density_plot.png"))
  p
dev.off()

######
# check normal infercnv neutral distribution as cancer is weird:
normal <- unlist(
  read.table(
    "/share/ScratchGeneral/jamtor/projects/thesis/Figure_2.4_artefact_analysis/results/infercnv/CID4520N/No_removed/samples_mode/infercnv.observations.txt", 
    as.is=T, header=T
  )
)
######

#density_plot <- density(infercnv_bulk)
## divide values into 5 clusters:
#infercnv_clusters <- kmeans(infercnv_bulk, 5)
#cluster_signal <- data.frame(
#  score = infercnv_bulk,
#  cluster = infercnv_clusters$cluster
#)
#split_signal <- split(cluster_signal, cluster_signal$cluster)
#
## determine means of signal groups:
#group_means <- lapply(split_signal, function(x) mean(x$score))
#group_ranges <- lapply(split_signal, function(x) range(x$score))



