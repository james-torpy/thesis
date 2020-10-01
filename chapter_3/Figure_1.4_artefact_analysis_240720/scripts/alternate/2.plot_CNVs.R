#! /share/ClusterShare/software/contrib/briglo/octR/src/R-3.6.0/builddir/bin/Rscript

###############################################################################
### Simulates cancer datasets taking as inputs number of CNVs, CNV lengths, ###
### gain/loss multipliers and downsamples required                          ###
###############################################################################

args = commandArgs(trailingOnly=TRUE)

project_name <- "thesis"
subproject_name <- "Figure_2.4_artefact_analysis"
sample_name <- args[1]
print(paste0("sample name = ", sample_name))
analysis_mode <- args[2]
print(paste0("Analysis mode of normal InferCNV run = ", 
  analysis_mode, "_mode"))
remove_cell_types <- args[3]
print("Cell types removed:")
print(remove_cell_types)
<<<<<<< HEAD:Figure_2.4_artefact_analysis/scripts/2.plot_normal.R
min_artefact_proportion <- as.numeric(args[4])
print("Minimum artefact proportion = ")
print(min_artefact_proportion)
min_artefact_length <- as.numeric(args[5])
print("Minimum artefact length = ")
print(min_artefact_length)

#project_name <- "thesis"
#subproject_name <- "Figure_2.4_artefact_analysis"
#sample_name <- "CID4520N"
#print(paste0("sample name = ", sample_name))
#analysis_mode <- "samples"
#print(paste0("Analysis mode of normal InferCNV run = ", 
#  analysis_mode, "_mode"))
#remove_cell_types <- "T_cells"
#print("Cell types removed:")
#min_artefact_proportion <- 0.5
#print("Minimum artefact proportion = ")
#print(min_artefact_proportion)
#min_artefact_length <- 10
#print("Minimum artefact length = ")
#print(min_artefact_length)
#neutral_signal_range <- c(0.97, 1.03)
=======
neutral_signal_range <- unlist(strsplit(args[4], split = "_"))
print("Neutral signal range = ")
print(neutral_signal_range)
min_artefact_proportion <- as.numeric(args[5])
print("Minimum artefact proportion = ")
print(min_artefact_proportion)
min_artefact_length <- as.numeric(args[6])
print("Minimum artefact length = ")
print(min_artefact_length)

project_name <- "thesis"
subproject_name <- "Figure_2.4_artefact_analysis"
sample_name <- "CID4520N"
print(paste0("sample name = ", sample_name))
analysis_mode <- "samples"
print(paste0("Analysis mode of normal InferCNV run = ", 
  analysis_mode, "_mode"))
remove_cell_types <- "Endothelial"
print("Cell types removed:")
neutral_signal_range <- unlist(strsplit("0.97_1.03", split = "_"))
print("Neutral signal range = ")
print(neutral_signal_range)
min_artefact_proportion <- 0.5
print("Minimum artefact proportion = ")
print(min_artefact_proportion)
min_artefact_length <- 10
print("Minimum artefact length = ")
print(min_artefact_length)
remove_cell_types <- "No"
>>>>>>> starting_fig_2.7:Figure_2.6_artefact_analysis/scripts/alternate/2.plot_CNVs.R

print(paste0("Project name = ", project_name))
print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(Seurat)
library(cowplot)
library(rlist, lib.loc = lib_loc)
library(GenomicRanges)
library(naturalsort, lib.loc = lib_loc)
library(splatter, lib.loc = lib_loc)
library(ggplot2)
library(org.Hs.eg.db)
library(data.table)
library(ComplexHeatmap, lib.loc = lib_loc)
library(circlize, lib.loc = lib_loc)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", 
  project_name, "/", subproject_name, "/")
ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/functions/")
results_dir <- paste0(project_dir, "results/")

in_dir <- paste0(results_dir, "infercnv/", sample_name, "/", 
<<<<<<< HEAD:Figure_2.4_artefact_analysis/scripts/2.plot_normal.R
  remove_cell_types, "_removed/", analysis_mode, "_mode/")
=======
  analysis_mode, "_mode/")
>>>>>>> starting_fig_2.7:Figure_2.6_artefact_analysis/scripts/alternate/2.plot_CNVs.R
input_dir <- paste0(results_dir, "infercnv/", sample_name, 
  "/input_files/")

Robject_dir <- paste0(in_dir, "Rdata/")
system(paste0("mkdir -p ", Robject_dir))
plot_dir <- paste0(in_dir, "plots/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(in_dir, "tables/")
system(paste0("mkdir -p ", table_dir))

print(paste0("Sample directory = ", in_dir))
print(paste0("Robject directory = ", Robject_dir))
print(paste0("Plot dir = ", plot_dir))
print(paste0("Table directory = ", table_dir))

print(paste0("Plotting InferCNV heatmap of ", sample_name))


################################################################################
### 0. Define functions and choose seeds ###
################################################################################

fetch_chromosome_boundaries <- dget(paste0(func_dir, 
  "fetch_chromosome_boundaries.R"))

create_extended_vector <- function(df, column) {
  for (j in 1:nrow(df)) {
    vector_length <- length(seq(df$start[j], df$end[j]))
    if (j==1) {
      result_vector <- as.character(
        rep(
          eval(parse(text=paste0("df$", column, "[j]"))), vector_length
        )
      )
    } else {
      result_vector <- c(
        result_vector,
        as.character(
          rep(
            eval(parse(text=paste0("df$", column, "[j]"))), vector_length
          )
        )
      )
    }
  }
  return(result_vector)
}


###################################################################################
### 1. Load InferCNV output and screen for artefacts ###
###################################################################################

if (
  !file.exists(
    paste0(Robject_dir, "1a.artefact_indices.Rdata")
  ) | 
  !file.exists(
    paste0(Robject_dir, "1b.artefact_genes.Rdata")
  ) |
  !file.exists(
    paste0(Robject_dir, "1c.artefact_annotation.Rdata")
  )
) {

  # load normal infercnv output:
  print("Loading normal InferCNV output files...")
  epithelial_heatmap <- as.data.frame(t(read.table(paste0(in_dir, 
    "infercnv.12_denoised.observations.txt"))))
  
  # determine average signal across all cells for each gene:
  heatmap_averages <- apply(epithelial_heatmap, 2, mean)
  
  for ( i in 1:(ncol(epithelial_heatmap)-(min_artefact_length-1)) ) {
<<<<<<< HEAD:Figure_2.4_artefact_analysis/scripts/2.plot_normal.R

   averages_per_cell <- apply(epithelial_heatmap[i:(i+(min_artefact_length-1))], 1, mean)
    
    # if at least 90% of cells have average signal in loss range, label as loss
    # and vice versa for gains:
    if ( length(which(averages_per_cell < neutral_signal_range[1])) > 
      floor(min_artefact_proportion*length(averages_per_cell)) ) {

      segment_record <- data.frame(
        start = i,
        end = i+(min_artefact_length-1),
        average_signal = mean(heatmap_averages[i:(i+(min_artefact_length-1))]),
        type = "loss"
      )

=======

   averages_per_cell <- apply(epithelial_heatmap[i:(i+(min_artefact_length-1))], 1, mean)
    
    # if at least 90% of cells have average signal in loss range, label as loss
    # and vice versa for gains:
    if ( length(which(averages_per_cell < neutral_signal_range[1])) > 
      floor(min_artefact_proportion*length(averages_per_cell)) ) {

      segment_record <- data.frame(
        start = i,
        end = i+(min_artefact_length-1),
        average_signal = mean(heatmap_averages[i:(i+(min_artefact_length-1))]),
        type = "loss"
      )

>>>>>>> starting_fig_2.7:Figure_2.6_artefact_analysis/scripts/alternate/2.plot_CNVs.R
    } else if ( length(which(averages_per_cell > neutral_signal_range[2])) > 
      floor(min_artefact_proportion*length(averages_per_cell)) ) {

      segment_record <- data.frame(
        start = i,
        end = i+(min_artefact_length-1),
        average_signal = mean(heatmap_averages[i:(i+(min_artefact_length-1))]),
        type = "gain"
      )

    } else {

      segment_record <- data.frame(
        start = i,
        end = i+(min_artefact_length-1),
        average_signal = mean(heatmap_averages[i:(i+(min_artefact_length-1))]),
        type = "no_artefact"
      )

    }

    if (i==1) {
      artefact_record <- segment_record
    } else {
      artefact_record <- rbind(artefact_record, segment_record)
    }

    artefact_record$type <- as.character(artefact_record$type)
  
  }
  
  split_record <- split(artefact_record, rleid(artefact_record$type))

  final_artefact_record <- lapply(split_record, function(x) {
    return(
      data.frame(
        start = x$start[1],
        end = x$end[nrow(x)],
        type = x$type[1]
      )
    )
  })
  final_artefact_record <- do.call("rbind", final_artefact_record)
  final_artefact_record$type <- as.character(final_artefact_record$type)

  # fix overlapping indices by prioritising artefact annotation over 
  # no_artefact:
  no_artefact_indices <- which(final_artefact_record$type == "no_artefact")
  for (ind in no_artefact_indices) {
    print(ind)
    # if the first non-artefact region is not at the start of the genome,
    #  make start of non-artefact region one after the end of the previous 
    # artefact region
    if (final_artefact_record$start[ind] != 1) {
      final_artefact_record$start[ind] <- final_artefact_record$end[ind-1] + 1
    }
    # if the last non-artefact region is not at the end of the genome,
    # make end of non-artefact region one before the start of the next 
    # artefact region: 
    if (final_artefact_record$end[ind] != ncol(epithelial_heatmap)) {
      if (final_artefact_record$start[ind] < final_artefact_record$start[ind+1]) {
        final_artefact_record$end[ind] <- final_artefact_record$start[ind+1] - 1        
      } else {
        final_artefact_record$type[ind] <- "remove"
        final_artefact_record$start[ind+1] <- final_artefact_record$end[ind-1] + 1
      }
    }
  }
  final_artefact_record <- final_artefact_record[
    final_artefact_record$type != "remove",
  ]

  artefact_indices <- final_artefact_record[
    final_artefact_record$type != "no_artefact",
  ]

  # fetch chromosome info:
  chr_data <- fetch_chromosome_boundaries(epithelial_heatmap, ref_dir)

  if (nrow(artefact_indices) > 0) {
    for (m in 1:nrow(artefact_indices)) {
      if (m==1) {
        heatmap_artefact <- artefact_indices$start[m]:artefact_indices$end[m]
      } else {
        heatmap_artefact <- c(heatmap_artefact,
          artefact_indices$start[m]:artefact_indices$end[m]
        )
      }   
    }  

    for (k in 1:length(chr_data$ends)) {
      if (k==1) {
  
        artefact_indices$start_chr[
          artefact_indices$start <= chr_data$ends[k]
        ] <- names(chr_data$ends)[k]
  
        artefact_indices$end_chr[
          artefact_indices$end <= chr_data$ends[k]
        ] <- names(chr_data$ends)[k]
  
      } else {
  
        artefact_indices$start_chr[
          artefact_indices$start <= chr_data$ends[k] & 
          artefact_indices$start > chr_data$ends[k-1]
        ] <- names(chr_data$ends)[k]
  
        artefact_indices$end_chr[
          artefact_indices$end <= chr_data$ends[k] & 
          artefact_indices$end > chr_data$ends[k-1]
        ] <- names(chr_data$ends)[k]
  
      }
    }
  
    # identify artefact genes and split by artefact:
    for (l in 1:nrow(artefact_indices)) {

      gene_df <- data.frame(
        indices = artefact_indices$start[l]:artefact_indices$end[l],
        genes = colnames(epithelial_heatmap)[
          artefact_indices$start[l]:artefact_indices$end[l]
        ],
        start_chr = artefact_indices$start_chr[l],
        end_chr = artefact_indices$end_chr[l]
      )

      if (l==1) {
        artefact_genes <- list(gene_df)
      } else {
        artefact_genes[[l]] <- gene_df
      }

    }

    # create artefact annotation:
    artefact_vector <- factor(
      create_extended_vector(final_artefact_record, "type"),
      levels = c("no_artefact", "gain", "loss")
    )
    cols <- c("#E5E4DF", "#BF3667", "#58B9DB")
    names(cols) <- levels(artefact_vector)
    
    artefact_annotation <- HeatmapAnnotation(
      artefact_call = artefact_vector,
      col = list(artefact_call = cols),
      annotation_name_side = "left",
      annotation_legend_param = list(title = "", 
        labels_gp = gpar(fontsize = 12))
    )
    artefact_annotation@name <- "artefact"

    saveRDS(
      artefact_indices, 
      paste0(Robject_dir, "1a.artefact_indices.Rdata")
    )
    saveRDS(
      artefact_genes, 
      paste0(Robject_dir, "1b.artefact_genes.Rdata")
    )
    saveRDS(
      artefact_annotation, 
      paste0(Robject_dir, "1c.artefact_annotation.Rdata")
    )

  }

} else {

  artefact_record <- readRDS(
    paste0(Robject_dir, "1a.artefact_indices.Rdata")
  )
  artefact_genes <- readRDS(
    paste0(Robject_dir, "1b.artefact_genes.Rdata")
  )
  artefact_annotation <- readRDS(
    paste0(Robject_dir, "1c.artefact_annotation.Rdata")
  )

  # load normal infercnv output:
  print("Loading normal InferCNV output files...")
  epithelial_heatmap <- as.data.frame(t(read.table(paste0(in_dir, 
    "infercnv.12_denoised.observations.txt"))))

  # fetch chromosome info:
  chr_data <- fetch_chromosome_boundaries(epithelial_heatmap, ref_dir)

}


###################################################################################
### 2. Plot CNV heatmap ###
###################################################################################

# prepare df for plotting:
plot_object <- epithelial_heatmap
colnames(plot_object) <- rep("la", ncol(plot_object))
plot_object <- as.matrix(plot_object)

# define heatmap colours:
na_less_vector <- unlist(plot_object)
na_less_vector <- na_less_vector[!is.na(na_less_vector)]
heatmap_cols <- colorRamp2(c(min(na_less_vector), 1, max(na_less_vector)), 
      c("#00106B", "white", "#680700"), space = "sRGB")
print("Generating filtered normal heatmap...")
final_heatmap <- Heatmap(
  plot_object, name = paste0("hm"), 
  col = heatmap_cols,
  cluster_columns = F, cluster_rows = F,
  show_row_names = F, show_column_names = F,
  #column_names_gp = gpar(col = "white"),
  show_row_dend = F,
  bottom_annotation = artefact_annotation,
  heatmap_legend_param = list(title = "CNV\nscore", 
  at = c(round(min(na_less_vector), 2), 1, round(max(na_less_vector), 2)),
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
# plot final annotated heatmap:
pdf(
  paste0(plot_dir, "infercnv_artefact_plot.pdf"), 
  height = 13, width = 20
)   
  
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

dev.off()


