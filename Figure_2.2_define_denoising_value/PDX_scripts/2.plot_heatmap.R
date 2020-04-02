#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

project_name <- "thesis"
subproject_name <- "Figure_2.2_define_denoising_value"
args = commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
analysis_mode <- args[2]
denoise_value <- args[3]

sample_name <- "PDX4386"
analysis_mode <- "samples"
denoise_value <- "1.5"

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("analysis_mode = ", nUMI_threshold))
print(paste0("denoise_value = ", as.character(nGene_threshold)))

library(reshape2)
library(dplyr)
library(cowplot)
library(viridis)
library(GenomicRanges)
library(ggplot2)

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(Polychrome, lib.loc=lib_loc)
library(ComplexHeatmap, lib.loc=lib_loc)
library(circlize, lib.loc = lib_loc)
#library(scales, lib.loc = lib_loc)
library(fpc, lib.loc = lib_loc)
library(naturalsort, lib.loc = lib_loc)
  
home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", 
  project_name, "/", subproject_name, "/")
ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/functions/")
results_dir <- seurat_path <- paste0(project_dir, "results/")

in_dir <- paste0(results_dir, "infercnv/", sample_name, "/",  
  analysis_mode, "_mode/", denoise_value, "_denoising/")

Robject_dir <- paste0(in_dir, "Rdata/")
system(paste0("mkdir -p ", Robject_dir))
plot_dir <- paste0(in_dir, "plots/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(in_dir, "tables/")
system(paste0("mkdir -p ", table_dir))

input_dir <- paste0(in_path, "input_files/")

print(paste0("In directory = ", in_dir))
print(paste0("Reference directory = ", ref_dir))
print(paste0("R function directory = ", func_dir))
print(paste0("R object directory = ", Robject_dir))
print(paste0("Table directory = ", table_dir))
print(paste0("Plot directory = ", plot_dir))


################################################################################
### 0. Define functions ###
################################################################################

fetch_chromosome_boundaries <- dget(paste0(func_dir, 
  "fetch_chromosome_boundaries.R"))


################################################################################
### 1. Load InferCNV output and create heatmap and metadata dfs ###
################################################################################

# load InferCNV output:
print("Loading InferCNV output files...")
if (denoise_value == "no") {
  epithelial_heatmap <- as.data.frame(t(read.table(paste0(in_dir, 
    "infercnv.observations.txt"))))
} else {
  epithelial_heatmap <- as.data.frame(t(read.table(paste0(in_dir, 
    "infercnv.12_denoised.observations.txt"))))
}



######

all_array_CNVs <- read.table(paste0(ref_dir, "all_array_CNVs.txt"))
colnames(all_array_CNVs) <- gsub("CID4499_1", "CID44991", colnames(all_array_CNVs))
if (sample_type == "cancer") {
  if (!file.exists(paste0(Robject_dir, "array_CNV_annotation.Rdata"))) {
    array_CNV_annotation <- create_array_CNV_annotation(epithelial_heatmap, all_array_CNVs)
    saveRDS(array_CNV_annotation, paste0(Robject_dir, "array_CNV_annotation.Rdata"))
    grid_array_heatmap <- grid.grabExpr(draw(array_CNV_annotation$array_CNV_heatmap, 
      heatmap_legend_side = "left"))
  } else {
    array_CNV_annotation <- readRDS(paste0(Robject_dir, "array_CNV_annotation.Rdata"))
    grid_array_heatmap <- grid.grabExpr(draw(array_CNV_annotation$array_CNV_heatmap, 
      heatmap_legend_side = "left"))
  }

  # remove normals and unassigned from epithelial_heatmap:
  cancer_only_heatmap <- epithelial_heatmap[
    epithelial_metadata$cell_ids[epithelial_metadata$normal_cell_call == "cancer"],
  ]

  # correlate infercnv CNV profiles with array profiles:
  average_heatmap <- apply(cancer_only_heatmap, 2, mean)

  cor_pearson <- cor.test(
    as.numeric(array_CNV_annotation$array_CNVs), 
    as.numeric(average_heatmap), 
    method = "pearson",
    alternative = "greater"
  )
  CNV_correlation_result <- data.frame(cor_pearson$estimate, cor_pearson$p.value)

}


################################################################################
### 10. Generate heatmap ###
################################################################################

# prepare df for plotting:
plot_object <- epithelial_heatmap
colnames(plot_object) <- rep("la", ncol(plot_object))
plot_object <- as.matrix(plot_object)

# define heatmap colours:
na_less_vector <- unlist(plot_object)
na_less_vector <- na_less_vector[!is.na(na_less_vector)]
heatmap_cols <- colorRamp2(c(min(na_less_vector), 1, max(na_less_vector)), 
      c("#00106B", "white", "#680700"), space = "sRGB")
print("Generating final heatmap...")

final_heatmap <- Heatmap(
  plot_object, name = paste0("hm"), 
  col = heatmap_cols,
  cluster_columns = F, cluster_rows = F,
  show_row_names = F, show_column_names = F,
  #column_names_gp = gpar(col = "white"),
  show_row_dend = F,
  bottom_annotation = accuracy_annotation,
  heatmap_legend_param = list(title = "CNV\nscore", 
  at = c(round(min(na_less_vector), 1), 1, round(max(na_less_vector), 1)),
  color_bar = "continuous", grid_height = unit(1.5, "cm"), 
  grid_width = unit(1.5, "cm"), legend_direction = "horizontal",
  title_gp = gpar(fontsize = 18, fontface = "bold"), 
  labels_gp = gpar(fontsize = 12)),
  use_raster = T, raster_device = c("png")
)

annotated_heatmap <- grid.grabExpr(
  draw(final_heatmap, gap = unit(6, "mm"), heatmap_legend_side = "left",
  annotation_legend_side = "left")
)
dev.off()

# determine where starting co-ordinates for heatmap are based upon longest cluster name
# (0.00604 units per character):
longest_cluster_name <- max(nchar(unique(as.character(epithelial_metadata$cell_type))))
x_coord <- longest_cluster_name*0.0037


################################################################################
### 11. Plot heatmap ###
################################################################################

# plot final annotated heatmap:
pdf(paste0(plot_dir, "infercnv_plot_with_accuracy_metrics.pdf"), height = 13, width = 20)   
  
  grid.newpage()
  pushViewport(viewport(x = 0.01, y = 0.16, width = 0.99, height = 0.78, 
    just = c("left", "bottom")))
    grid.draw(annotated_heatmap)
    decorate_heatmap_body("hm", {
      for ( e in 1:length(chr_data$end_pos) ) {
        grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), 
          gp = gpar(lwd = 1, col = "#383838"))
        grid.text(gsub("chr", "", names(chr_data$lab_pos)[e]), chr_data$lab_pos[e], 
          unit(0, "npc") + unit(-3.5, "mm"), gp=gpar(fontsize=18))
      }
    })
  popViewport()

  pushViewport(viewport(x = x_coord + 0.965, y = 0.05, 
    width = 0.88, height = 0.088, just = c("right", "bottom")))
    grid.draw(grid_sim_plot)
  popViewport()

  pushViewport(viewport(x=x_coord + 0.0765, y=0.89, width = 0.1, height = 0.1, 
    just = "right"))
    grid.text(paste0("No. genes = ", ncol(epithelial_heatmap)), 
      gp=gpar(fontsize=16))
  popViewport()
  pushViewport(viewport(x=x_coord + 0.071, y=0.86, width = 0.1, height = 0.1, 
    just = "right"))
    grid.text(
      paste0(
      "Mean UMI = ", final_accuracy_metrics["mean_UMI",]
      ), gp=gpar(fontsize=16)
    )
  popViewport()
  pushViewport(viewport(x=x_coord + 0.075, y=0.83, width = 0.1, height = 0.1, 
    just = "right"))
    grid.text(
      paste0(
      "Sensitivity = ", 
      final_accuracy_metrics["sensitivity",]
      ), gp=gpar(fontsize=16)
    )
  popViewport()
  pushViewport(viewport(x=x_coord + 0.07, y=0.8, width = 0.1, height = 0.1, 
    just = "right"))
    grid.text(
      paste0(
      "Specificity = ", 
      final_accuracy_metrics["specificity",]
      ), gp=gpar(fontsize=16)
    )
  popViewport() 
  pushViewport(viewport(x=x_coord + 0.073, y=0.77, width = 0.1, height = 0.1, 
    just = "right"))
    grid.text(
      paste0(
      "Precision = ", 
      final_accuracy_metrics["precision",]
      ), gp=gpar(fontsize=16)
    )
  popViewport() 
  pushViewport(viewport(x=x_coord + 0.072, y=0.74, width = 0.1, height = 0.1, 
    just = "right"))
    grid.text(
      paste0(
      "F1 score = ", 
      final_accuracy_metrics["F1",]
      ), gp=gpar(fontsize=16)
    )
  popViewport() 
  pushViewport(viewport(x=x_coord + 0.073, y=0.71, width = 0.1, height = 0.1, 
    just = "right"))
    grid.text(
      paste0(
      "No. wrong calls = ", 
      final_accuracy_metrics["wrong_call",]
      ), gp=gpar(fontsize=16)
    )
  popViewport() 

  if (exists("cor_result")) {
    pushViewport(viewport(x=x_coord + 0.06, y=0.67, width = 0.1, height = 0.1, 
      just = "right"))
      grid.text(paste0("Correlation with\noriginal = ", round(cor_result$R_squared, 3)), 
        gp=gpar(fontsize=16))
    popViewport()
    pushViewport(viewport(x=x_coord + 0.05, y=0.635, width = 0.1, height = 0.1, 
      just = "right"))
      grid.text(paste0("p.val = ", round(cor_result$p_val, 3)), 
        gp=gpar(fontsize=16))
    popViewport()
  }
  
dev.off()

print(paste0("Heatmap created, output in ", plot_dir))

# print epithelial_metadata as table:
write.table(epithelial_metadata, paste0(table_dir, "epithelial_metadata.txt"), 
  sep="\t", quote=F, row.names=F, col.name=T)

print(paste0("Heatmap created, output in ", plot_dir))

##convert pdf to png:
#system(paste0("for p in ", plot_dir, "*.pdf; do echo $p; f=$(basename $p); echo $f; ",
#              "new=$(echo $f | sed 's/.pdf/.png/'); echo $new; ", 
#              "convert -density 150 ", plot_dir, "$f -quality 90 ", plot_dir, "$new; done"))
#


################################################################################
### development ###
################################################################################

#false_positive_segments <- list(
#  gain = false_positive_segments$gain,
#  loss = data.frame(
#    start = c(5400, 5560, 5800),
#    end = c(5420, 5580, 5820),
#    type = rep("neutral", 3),
#    CNV_call = c("loss", "neutral", "loss"),
#    accuracy = c("false_positive", "true_negative", "false_positive"),
#    slice_type = rep("full", 3)
#  )
#)
#
#false_positive_grs <- list(
#  gain = false_positive_grs$gain,
#  loss = GRanges(
#    seqnames = rep("all", 3),
#    ranges = IRanges(
#      start = c(5400, 5560, 5800),
#      end = c(5420, 5580, 5820)
#    ),
#    strand = Rle(rep("*", nrow(temp_coords))),
#    accuracy = c("false_positive", "true_negative", "false_positive"),
#    CNV_call = c("loss", "neutral", "loss"),
#    slice_type = rep("full", 3) 
#  )
#)
#
#o=1
#lapply(neutral_regions_list, function(neutral_region) {
#  print(o)
#  # set default accuracy and CNV_calls to be changed if necessary:
#  neutral_region$accuracy <- "true negative"
#  neutral_region$CNV_call <- "neutral"
#  neutral_region$slice_type <- "full"
#  neutral_region$slice_type[
#    neutral_region$end-neutral_region$start+1 < 21
#  ] <- "partial"
#  
#  # subset neutral_region:
#  neutral_region <- subset(neutral_region, select = c(start, end, type, 
#    CNV_call, accuracy, slice_type))
#
#  # fetch the mean signal at each gene:
#  mean_segment_signal <- average_epithelial[
#    neutral_region$start:neutral_region$end
#  ]
#  print(paste0(length(mean_segment_signal), " genes in segment"))
#
#  index_types <- c("gain", "loss")
#
#  for (i in 1:length(index_types)) {
#
#    # determine genes with a CNV signal:
#    if (index_types[i] == "gain") {
#      signal_indices <- which(mean_segment_signal > neutral_signal_range[2])
#      print(paste0("Length of gain genes = ", length(signal_indices)))
#    } else if (index_types[i] == "loss") {
#      signal_indices <- which(mean_segment_signal < neutral_signal-range[1])
#      print(paste0("Length of loss genes = ", length(signal_indices)))
#    }
#  }
#  o <<- o+1
#  return(neutral_region)
#})


