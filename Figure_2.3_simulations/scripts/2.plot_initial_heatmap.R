#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

#### Generate CNV heatmaps with following annotations: ###
# nUMI
# nGene
# CNA values
# SNP array CNVs

project_name <- "thesis"
subproject_name <- "Figure_2.3_simulations"
args = commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
include_t_cells <- as.logical(args[2])
cancer_x_threshold_sd_multiplier <- args[3]
cancer_y_threshold_sd_multiplier <- args[4]
normal_x_threshold_sd_multiplier <- args[5]
normal_y_threshold_sd_multiplier <- args[6]
remove_normals <- as.logical(args[7])
reclustered_group_annotation <- as.logical(args[8])

#sample_name <- "CID4520N"
#include_t_cells <- TRUE
#cancer_x_threshold_sd_multiplier <- 2
#cancer_y_threshold_sd_multiplier <- 1.5
#normal_x_threshold_sd_multiplier <- 1
#normal_y_threshold_sd_multiplier <- 1.25
#include_normals <- FALSE
#reclustered_group_annotation <- TRUE

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("T-cells included?", as.character(include_t_cells)))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.5dev/"
library(cluster, lib.loc = lib_loc)
library(Seurat)
library(RColorBrewer)
library(ComplexHeatmap, lib.loc=lib_loc)
library(circlize, lib.loc = lib_loc)
library(reshape2)
library(ggplot2)
library(scales, lib.loc = lib_loc)
library(fpc, lib.loc = lib_loc)
library(dplyr)
library(naturalsort, lib.loc = lib_loc)
library(cowplot)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", 
  project_name, "/", subproject_name, "/")
ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/functions/")
results_dir <- seurat_path <- paste0(project_dir, "results/")
seurat_dir <- paste0(project_dir, "raw_files/seurat_objects/", 
  sample_name, "/")

if (include_t_cells) {
  in_dir <- paste0(results_dir, "infercnv/t_cells_included/", sample_name, "/")
} else {
  in_dir <- paste0(results_dir, "infercnv/t_cells_excluded/", sample_name, "/")
}

Robject_dir <- paste0(in_dir, "Rdata/")
system(paste0("mkdir -p ", Robject_dir))
plot_dir <- paste0(in_dir, "plots/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(in_dir, "tables/")
system(paste0("mkdir -p ", table_dir))

print(paste0("In directory = ", in_dir))
print(paste0("Reference directory = ", ref_dir))
print(paste0("R function directory = ", func_dir))
print(paste0("R object directory = ", Robject_dir))
print(paste0("Table directory = ", table_dir))
print(paste0("Plot directory = ", plot_dir))


################################################################################
### 0. Define functions ###
################################################################################

gg_color_hue <- dget(paste0(func_dir, "gg_color_hue.R"))
prepare_infercnv_metadata <- dget(paste0(func_dir, 
  "prepare_infercnv_metadata.R"))
fetch_chromosome_boundaries <- dget(paste0(func_dir, 
  "fetch_chromosome_boundaries.R"))


################################################################################
### 1. Load InferCNV output and create heatmap and metadata dfs ###
################################################################################

if (!file.exists(paste0(Robject_dir, "/1b.initial_epithelial_metadata.Rdata"))) {
  # load InferCNV output:
  print("Loading InferCNV output files...")
  infercnv_output <- as.data.frame(t(read.table(paste0(in_dir, 
  	"infercnv.12_denoised.observations.txt"))))

  # create cluster metadata df:
  seurat_10X <- readRDS(paste0(seurat_dir, "03_seurat_object_processed.Rdata"))
  Idents(seurat_10X) <- seurat_10X@meta.data$PC_A_res.1
  metadata <- prepare_infercnv_metadata(seurat_10X, subset_data = F, 
    as.data.frame(t(infercnv_output)), for_infercnv=F)
  metadata_df <- metadata$metadata
  # determine the epithelial cells and only include these in heatmap:
  print(paste0("Number of heatmap rows before non-epithelial thrown: ", 
  	nrow(infercnv_output)))
  epithelial_ids <- metadata_df$cell_ids[grep("pithelial", metadata_df$cell_type)]
  epithelial_heatmap <- infercnv_output[rownames(infercnv_output) %in% epithelial_ids,]
  print(paste0("Number of heatmap rows after non-epithelial thrown: ", 
  	nrow(epithelial_heatmap)))

  # create epithelial_metadata df and only include epithelial cells in epithelial_heatmap:
  print("Creating epithelial metadata df...")
  epithelial_metadata <- metadata_df[rownames(epithelial_heatmap),]

  saveRDS(epithelial_heatmap, paste0(Robject_dir, "/1a.initial_epithelial_heatmap.Rdata"))
  saveRDS(epithelial_metadata, paste0(Robject_dir, "/1b.initial_epithelial_metadata.Rdata"))

} else {

	print("Loading heatmap and metadata dfs...")
	epithelial_heatmap <- readRDS(paste0(Robject_dir, "/1a.initial_epithelial_heatmap.Rdata"))
	epithelial_metadata <- readRDS(paste0(Robject_dir, "/1b.initial_epithelial_metadata.Rdata"))

}


################################################################################
### 2. Add QC metadata ###
################################################################################

if (!file.exists(paste0(Robject_dir, 
  "2b.epithelial_metadata_with_cell_type_and_QC.Rdata"))) {
  # order epithelial metadata cell type cluster levels:
  epithelial_metadata$cell_type <- factor(
    epithelial_metadata$cell_type,
    levels = naturalsort(unique(epithelial_metadata$cell_type))
  )
  print(paste0(
    "Are epithelial_metadata rownames in the same order as epithelial_heatmap? ",
    identical(rownames(epithelial_heatmap), rownames(epithelial_metadata))
  ))

  # add nUMI and nGene data to epithelial_metadata:
  print("Adding QC metrics to epithelial metadata df...")
  if (!exists("seurat_10X")) {
    seurat_10X <- readRDS(
      paste0(seurat_dir, "03_seurat_object_processed.Rdata")
    )
    Idents(seurat_10X) <- seurat_10X@meta.data$PC_A_res.1
  }
  QC <- data.frame(
    row.names = names(Idents(seurat_10X)),
    nUMI = seurat_10X@meta.data$nCount_RNA,
    nGene = seurat_10X@meta.data$nFeature_RNA
  )
  QC <- QC[rownames(epithelial_metadata),]
  epithelial_metadata <- cbind(epithelial_metadata, QC)
  print(paste0(
    "Are epithelial_metadata rownames still in the same order as epithelial_heatmap?? ",
    identical(rownames(epithelial_heatmap), rownames(epithelial_metadata))
  ))
  saveRDS(
    epithelial_heatmap, paste0(Robject_dir, 
    "2a.epithelial_heatmap_with_cell_type_and_QC.Rdata")
  )
  saveRDS(
    epithelial_metadata, paste0(Robject_dir, 
    "2b.epithelial_metadata_with_cell_type_and_QC.Rdata")
  )
} else {
   epithelial_heatmap <- readRDS(
    paste0(Robject_dir, 
    "2a.epithelial_heatmap_with_cell_type_and_QC.Rdata")
  )
  epithelial_metadata <- readRDS(
    paste0(Robject_dir, 
    "2b.epithelial_metadata_with_cell_type_and_QC.Rdata")
  )
  print(paste0(
    "Are epithelial_metadata rownames still in the same order as epithelial_heatmap?? ",
    identical(rownames(epithelial_heatmap), rownames(epithelial_metadata))
  ))
}


################################################################################
### 3. Add CNA metadata ###
################################################################################

if (!file.exists(paste0(Robject_dir, 
  "3b.epithelial_heatmap_with_cell_type_QC_and_CNA_values.Rdata"))) {
  
  # determine CNA values and add to epithelial_metadata:
  print("Determining CNA values and adding to epithelial metadata df...")
  # scale infercnv values to -1:1, square values and take the mean:
  scaled_df <- as.data.frame(rescale(as.matrix(epithelial_heatmap), c(-1,1)))
  CNA_values <- apply(scaled_df, 1, function(y) {
    #y[is.na(y)] <- 0
    #scaled_y <- rescale(y, c(-1, 1))
    return(mean(y^2))
  })
  CNA_value_df <- data.frame(
    row.names = names(CNA_values),
    CNA_value = CNA_values
  )
  epithelial_metadata <- cbind(epithelial_metadata, CNA_value_df)
  print(paste0(
    "Are epithelial_metadata rownames still in the same order as epithelial_heatmap?? ",
    identical(rownames(epithelial_heatmap), rownames(epithelial_metadata))
  ))
  # determine correlation with top 5% cancer values and add to epithelial_metadata:
  print(paste0(
    "Determining correlation with top 5% cancer values and adding to epithelial ", 
    "metadata df..."
  ))
  
  saveRDS(
    epithelial_heatmap, paste0(Robject_dir, 
    "3a.epithelial_heatmap_with_cell_type_QC_and_CNA_values.Rdata")
  )
  saveRDS(
    epithelial_metadata, paste0(Robject_dir, 
    "3b.epithelial_heatmap_with_cell_type_QC_and_CNA_values.Rdata")
  )

} else {
  epithelial_heatmap <- readRDS(
    paste0(Robject_dir, 
    "3a.epithelial_heatmap_with_cell_type_QC_and_CNA_values.Rdata")
  )
  epithelial_metadata <- readRDS(
    paste0(Robject_dir, 
    "3b.epithelial_heatmap_with_cell_type_QC_and_CNA_values.Rdata")
  )
  print(paste0(
    "Are epithelial_metadata rownames still in the same order as epithelial_heatmap?? ",
    identical(rownames(epithelial_heatmap), rownames(epithelial_metadata))
  ))
}


################################################################################
### 4. Order heatmap, metadata and create heatmap annotations ###
################################################################################

# create CNA annotation:
CNA_value_annotation <- rowAnnotation(
  correlation_annotation = anno_barplot(
    epithelial_metadata$CNA_value,
    gp = gpar(
      col = "#D95F02", 
      width = unit(4, "cm")
    ), 
    border = FALSE, 
    which = "row", 
    axis = F
  )
)
CNA_value_annotation@name <- "CNA_value"
# create QC annotations:
nUMI_annotation <- rowAnnotation(
  correlation_annotation = anno_barplot(
    epithelial_metadata$nUMI,
    gp = gpar(
      col = "#D8B72E", 
      width = unit(4, "cm")
    ), 
    border = FALSE, 
    which = "row", 
    axis = F
  )
)
nUMI_annotation@name <- "nUMI"
nGene_annotation <- rowAnnotation(
  correlation_annotation = anno_barplot(
    epithelial_metadata$nGene, name = "nGene",
    gp = gpar(
      col = "#9ECAE1", 
      width = unit(4, "cm")
    ), 
    border = FALSE, 
    which = "row", 
    axis = F
  )
)
nGene_annotation@name <- "nGene"


################################################################################
### 5. Generate heatmap ###
################################################################################

# fetch chromosome boundary co-ordinates:
if (!file.exists(paste0(Robject_dir, "chromosome_data.Rdata"))) {
  chr_data <- fetch_chromosome_boundaries(epithelial_heatmap, ref_dir)
  saveRDS(chr_data, paste0(Robject_dir, "chromosome_data.Rdata"))
} else {
  chr_data <- readRDS(paste0(Robject_dir, "chromosome_data.Rdata"))
}
# prepare df for plotting:
plot_object <- epithelial_heatmap
colnames(plot_object) <- rep("la", ncol(plot_object))

# define heatmap colours:
na_less_vector <- unlist(plot_object)
na_less_vector <- na_less_vector[!is.na(na_less_vector)]
heatmap_cols <- colorRamp2(c(min(na_less_vector), 1, max(na_less_vector)), 
      c("#00106B", "white", "#680700"), space = "sRGB")
print("Generating final heatmap...")

# create main CNV heatmap:
final_heatmap <- Heatmap(
  plot_object, name = paste0("hm"), 
  col = heatmap_cols,
  cluster_columns = F, cluster_rows = F,
  show_row_names = F, show_column_names = T,
  column_names_gp = gpar(col = "white"),
  show_row_dend = F,
  heatmap_legend_param = list(title = "CNV\nvalue", color_bar = "continuous", 
    grid_height = unit(1.5, "cm"), grid_width = unit(1.5, "cm"), 
    legend_direction = "horizontal", title_gp = gpar(fontsize = 18, fontface = "bold"), 
    labels_gp = gpar(fontsize = 12)),
  use_raster = T, raster_device = c("png")
)

ht_list <- final_heatmap + CNA_value_annotation + nUMI_annotation + nGene_annotation

annotated_heatmap <- grid.grabExpr(
  draw(ht_list, gap = unit(6, "mm"), heatmap_legend_side = "left")
)

# determine where starting co-ordinates for heatmap are based upon longest cluster name
# (0.00604 units per character):
longest_cluster_name <- max(nchar(unique(as.character(epithelial_metadata$cell_type))))
x_coord <- longest_cluster_name*0.0037

# plot final annotated heatmap:
pdf(paste0(plot_dir, "infercnv_plot.pdf"), height = 13, width = 18)   
  
  grid.newpage()
    pushViewport(viewport(x = 0.005, y = 0.065, width = 0.99, height = 0.78, 
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

    pushViewport(viewport(x=x_coord + 0.875, y=0.025, width = 0.1, height = 0.1, 
      just = "bottom"))
      grid.text("CNA", rot=65)
    popViewport()
    pushViewport(viewport(x=x_coord + 0.892, y=0.025, width = 0.1, height = 0.1, 
      just = "bottom"))
      grid.text("nUMI", rot=65)
    popViewport()
    pushViewport(viewport(x=x_coord + 0.917, y=0.025, width = 0.1, height = 0.1, 
      just = "bottom"))
      grid.text("nGene", rot=65)
    popViewport()
    
dev.off()

print(paste0("Heatmap created, output in ", plot_dir))

# print epithelial_metadata as table:
write.table(epithelial_metadata, paste0(table_dir, "epithelial_metadata.txt"), 
  sep="\t", quote=F, row.names=F, col.name=T)

#convert pdf to png:
system(paste0("for p in ", plot_dir, "*.pdf; do echo $p; f=$(basename $p); echo $f; ",
"new=$(echo $f | sed 's/.pdf/.png/'); echo $new; ", 
"convert -density 150 ", plot_dir, "$f -quality 90 ", plot_dir, "$new; done"))

