#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

project_name <- "single_cell"
subproject_name <- "brca_mini_atlas_131019"

include_t_cells <- TRUE
min_CNV_proportion <- 0.5
min_CNV_length <- 10

print(paste0("Subproject name = ", subproject_name))
print(paste0("T-cells included?", as.character(include_t_cells)))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(data.table)
library(ComplexHeatmap, lib.loc = lib_loc)
library(grid)
library(ggplot2)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/",
  subproject_name, "/")
ref_dir <- paste0(project_dir, "/refs/")
func_dir <- paste0(project_dir, "/scripts/functions/")
results_dir <- paste0(project_dir, "results/")

if (include_t_cells) {
  in_dir <- paste0(results_dir, "infercnv/t_cells_included/")
} else {
  in_dir <- paste0(results_dir, "infercnv/t_cells_excluded/")
}

out_dir <- paste0(in_dir, "/CNV_distributions/")

Robject_dir <- paste0(out_dir, "Rdata/")
system(paste0("mkdir -p ", Robject_dir))
plot_dir <- paste0(out_dir, "plots/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(out_dir, "tables/")
system(paste0("mkdir -p ", table_dir))

print(paste0("In directory = ", in_dir))
print(paste0("Reference directory = ", ref_dir))
print(paste0("R function directory = ", func_dir))
print(paste0("R object directory = ", Robject_dir))
print(paste0("Table directory = ", table_dir))
print(paste0("Plot directory = ", plot_dir))


################################################################################
### 0. Load functions ###
################################################################################

detect_CNVs <- dget(paste0(func_dir, "detect_CNVs.R"))
plot_signal <- dget(paste0(func_dir, "plot_signal.R"))
print_signal <- dget(paste0(func_dir, "print_signal.R"))


################################################################################
### 1. Load inferred CNVs ###
################################################################################

sample_names <- list.files(in_dir, pattern = "CID")

# remove any matched datasets:
sample_names <- sample_names[!substring(sample_names, 8, 8) == 2]

# load inferred CNVs:
for (i in 1:length(sample_names)) {

  print(paste0("Loading inferred CNVs from ", sample_names[i], "..."))

  if (i==1) {
  	epithelial_heatmaps <- list(
      as.data.frame(
        t(
       	  read.table(
       		paste0(
       		  in_dir, sample_names[i], "/infercnv.12_denoised.observations.txt"
       		)
       	  )
        )
      )
    )
  } else {
  	epithelial_heatmaps[[i]] <- as.data.frame(
      t(
       	read.table(
       	  paste0(
       	  	in_dir, sample_names[i], "/infercnv.12_denoised.observations.txt"
       	  )
       	)
      )
    )
  }

}


################################################################################
### 2. Detect and plot CNVs in each dataset ###
################################################################################

# plot CNV signal:
signal_plots <- lapply(
  epithelial_heatmaps, 
  plot_signal,
  ref_dir
)

# detect CNVs:
detected_CNV_data <- lapply(
  epithelial_heatmaps, 
  detect_CNVs, 
  min_CNV_proportion,
  min_CNV_length,
  ref_dir
)
detected_CNVs <- lapply(detected_CNV_data, function(x) x$detected_CNVs)
CNV_by_gene <- lapply(detected_CNV_data, function(x) x$CNV_by_gene)

# annotate signal plot with detected CNVs:
CNV_annots <- lapply(CNV_by_gene, function(x) {

  CNV_annot_vector <- factor(
    x,
    levels = c("neutral", "gain", "loss")
  )
  CNV_cols <- structure(
    c("#E5E4DF", "#BF3667", "#58B9DB"),
    names = levels(CNV_annot_vector)
  )
  CNV_annot <- Heatmap(
    t(matrix(CNV_annot_vector)),
    name = "CNV_annot",
    col = CNV_cols,
    show_heatmap_legend = FALSE
  )
  CNV_annot@name <- "CNV_annot"
  CNV_annot_obj <- grid.grabExpr(
    draw(CNV_annot, heatmap_legend_side = "left")
  )
  dev.off()

  return(CNV_annot_obj)

})

# print signal plot:
for (i in 1:length(sample_names)) {
  print_signal(
  	sample_names[i],
  	signal_plots[[i]],
  	CNV_annots[[i]],
  	epithelial_heatmaps[[i]],
  	plot_dir
  )
}


################################################################################
### 3. Plot CNV counts and lengths ###
################################################################################

# count CNVs in each dataset:
CNV_counts <- lapply(detected_CNVs, function(x) {
  return(nrow(x[x$type != "neutral",]))
})

# record CNV lengths in each dataset:
CNV_lengths <- lapply(detected_CNVs, function(x) {
  CNV_only <- x[x$type != "neutral",]
  return(CNV_only$end-CNV_only$start)
})

# collate data and save:
CNV_dist_data <- data.frame(
  sample_name = sample_names,
  CNV_count = unlist(CNV_counts),
  mean_CNV_length = round(unlist(lapply(CNV_lengths, mean)), 0)
)

write.table(
  CNV_dist_data,
  paste0(table_dir, "CNV_metadata.txt"),
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

# plot distribution of CNV numbers across all:
count_density_plot <- density(CNV_dist_data$count)
pdf(paste0(plot_dir, "count_density_plot.pdf"))
  plot(count_density_plot, main=NA, xlab = "CNV count")
dev.off()

# plot distribution of CNV lengths across all:
length_density_plot <- density(unlist(CNV_lengths))
pdf(paste0(plot_dir, "length_density_plot.pdf"))
  plot(length_density_plot, main=NA, xlab = "CNV length")
dev.off()




#################################################################################################


#png(
#  paste0(plot_dir, "signal_plot_test.png"),
#  height = 8, 
#  width = 22,
#  res = 300,
#  units = "in"
#)
#
#  grid.newpage()
#
#    # plot signal:
#    pushViewport(viewport(x = 0.52, y = 0.5, width = 0.93, height = 0.9))
#      grid.draw(signal_plots[[2]])
#    popViewport()
#
#dev.off()










