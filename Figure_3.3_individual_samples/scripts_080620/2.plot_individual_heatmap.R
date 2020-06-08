#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

#### Generate CNV heatmaps with following annotations: ###
# nUMI
# nGene
# Expression clusters
# CNV subclusters

#### Generate following t-SNEs/UMAPS: ###
# Best matching resolution clusters annotated by cell type
# Broad cell marker feature plot
# Epithelial vs myoepithelial feature plot
# Luminal vs basal feature plot

project_name <- "thesis"
subproject_name <- "Figure_3.3_individual_samples"
args = commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
all_cell_PC <- args[2]
all_cell_res <- args[3]
malignant_PC <- args[4]
malignant_res <- args[5]
broad_markers <- unlist(
  strsplit(
    args[6],
    "_"
  )
)
epithelial_markers <- unlist(
  strsplit(
    args[7],
    "_"
  )
)
nUMI_threshold <- as.numeric(args[8])
nGene_threshold <- as.numeric(args[9])
min_proportion_cells_for_subcluster <- as.numeric(args[10])
random_tree_p <- as.numeric(args[11])
subcluster_by <- args[12]

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("All cell PC = ", all_cell_PC))
print(paste0("All cell res = ", all_cell_res))
print(paste0("Malignant PC = ", malignant_PC))
print(paste0("Malignant res = ", malignant_res))
print("Broad markers: ")
print(broad_markers)
print("Epithelial markers: ")
print(epithelial_markers)
print(paste0("nUMI threshold = ", nUMI_threshold))
print(paste0("nGene threshold = ", nGene_threshold))
print(paste0("Min proportion of cells for subcluster to be defined = ", 
  min_proportion_cells_for_subcluster))
print(paste0("Subclustering by pval ", random_tree_p))
print(paste0("Subclustering by ", subcluster_by))

#sample_name <- "CID4495"
#all_cell_PC <- "C"
#all_cell_res <- "PC_C_res.0.6"
#malignant_PC <- "D"
#malignant_res <- "SUBSET_D_res.0.8"
#broad_markers <- unlist(
#  strsplit(
#    c(
#      "ACTB_PTPRC_CD19_CD3D_CD68_PDGFRB_PECAM1_EPCAM"
#    ),
#    "_"
#  )
#)
#epithelial_markers <- unlist(
#  strsplit(
#    c(
#      "ACTB_EPCAM_KRT18_KRT5_MKI67"
#    ),
#    "_"
#  )
#)
#nUMI_threshold <- as.numeric("8000")
#nGene_threshold <- as.numeric("1300")
#min_proportion_cells_for_subcluster <- as.numeric("0.005")
#random_tree_p <- 0.05
#subcluster_by <- "random_trees"

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(MAST, lib.loc = lib_loc)
library(scales, lib.loc = lib_loc)
library(ggplot2)
library(Seurat)
library(infercnv, lib.loc = lib_loc)
library(RColorBrewer)
library(ComplexHeatmap, lib.loc=lib_loc)
library(circlize, lib.loc = lib_loc)
library(reshape2)
library(fpc, lib.loc = lib_loc)
library(dplyr)
library(naturalsort, lib.loc = lib_loc)
library(cowplot)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/",
  subproject_name, "/")
ref_dir <- paste0(project_dir, "/refs/")
func_dir <- paste0(project_dir, "/scripts/functions/")
results_dir <- paste0(project_dir, "results/")
raw_dir <- paste0(project_dir, "raw_files/")
seurat_dir <- paste0(raw_dir, "seurat_objects/", sample_name, "/")
in_dir <- paste0(results_dir, "infercnv/", sample_name, "/p_", 
  random_tree_p, "/")
input_dir <- paste0(in_dir, "/input_files/")

out_dir <- paste0(in_dir, "clustered_by_", subcluster_by, "/")
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
### 0. Define functions ###
################################################################################

fetch_chromosome_boundaries <- dget(paste0(func_dir, 
  "fetch_chromosome_boundaries.R"))
extra_colours <- c(brewer.pal(12, "Set3"), brewer.pal(12, "Paired"))
col_palette <- c(brewer.pal(8, "Dark2")[c(3:6,8)], brewer.pal(12, "Set3"), brewer.pal(8, "Accent"),
  "#660200", "#918940", "black", "#9ECAE1", "#3F007D", "#67000D", "#FDD0A2", "#08519C",
  "#DBB335", "#5AA050", "#807DBA", "#1A6000", "#F16913", "#FD8D3C", "#DEEBF7", 
  "#7F2704", "#DADAEB", "#FC9272", "#BCBDDC", extra_colours, "#b2182b", "#85929E", 
  "#9B59B6", "#74add1","#1b7837", "#b8e186", "#fed976","#e7298a", "#18ffff", "#ef6c00",
  "#A93226", "#E611ED","orange", "#b8bc53", "#5628ce", "#fa909c", "#8ff331","#270e26")
col_palette <- col_palette[-7]


################################################################################
### 1. Load InferCNV output and create heatmap and metadata dfs ###
################################################################################

if (!file.exists(paste0(Robject_dir, "/1b.initial_epithelial_metadata.Rdata"))) {

  # load InferCNV output:
  print("Loading InferCNV output files...")
  epithelial_heatmap <- as.data.frame(t(read.table(paste0(in_dir, 
  	"infercnv.12_denoised.observations.txt"))))

  seurat_malignant <- readRDS(paste0(seurat_dir, "05_seurat_object_malignant.Rdata"))
  Idents(seurat_malignant) <- paste0(
    seurat_malignant@meta.data$garnett_call_ext_major,
    " ",
    eval(parse(text = paste0("seurat_malignant@meta.data$", malignant_res)))
  )
  levels(Idents(seurat_malignant)) <- naturalsort(levels(Idents(seurat_malignant)))

  # define cluster annotation colours:
  cluster_number <- length(unique(Idents(seurat_malignant)))
  cluster_cols <- col_palette[1:cluster_number]
  names(cluster_cols) <- unique(Idents(seurat_malignant))

  saveRDS(cluster_cols, paste0(Robject_dir, "cluster_colours.Rdata"))

  # create cluster metadata df:
  epithelial_metadata <- data.frame(
  	cell_ids = names(Idents(seurat_malignant)),
  	cell_type = Idents(seurat_malignant),
  	nUMI = seurat_malignant@meta.data$nCount_RNA,
  	nGene = seurat_malignant@meta.data$nFeature_RNA
  )

  # remove non-malignant cells from heatmap:
  epithelial_heatmap <- epithelial_heatmap[
    rownames(epithelial_heatmap) %in% rownames(epithelial_metadata),
  ]
  # make metadata same order as heatmap:
  epithelial_metadata <- epithelial_metadata[
  	rownames(epithelial_heatmap),
  ]

  print(paste0(
    "Are epithelial_metadata rownames in the same order as epithelial_heatmap?? ",
    identical(rownames(epithelial_heatmap), rownames(epithelial_metadata))
  ))
  
  saveRDS(epithelial_heatmap, paste0(Robject_dir, 
    "/1a.initial_epithelial_heatmap.Rdata"))
  saveRDS(epithelial_metadata, paste0(Robject_dir, 
    "/1b.initial_epithelial_metadata.Rdata"))

} else {

  print("Loading heatmap and metadata dfs...")
  epithelial_heatmap <- readRDS(paste0(Robject_dir, 
    "/1a.initial_epithelial_heatmap.Rdata"))
  epithelial_metadata <- readRDS(paste0(Robject_dir, 
    "/1b.initial_epithelial_metadata.Rdata"))
  print(paste0(
    "Are epithelial_metadata rownames in the same order as epithelial_heatmap?? ",
    identical(rownames(epithelial_heatmap), rownames(epithelial_metadata))
  ))
  cluster_cols <- readRDS(paste0(Robject_dir, "cluster_colours.Rdata"))

}


################################################################################
### 2. Add annotations and create heatmap ###
################################################################################

if (!file.exists(paste0(
  plot_dir, "prefiltered_infercnv_plot_with_expression_clusters.png"
))) {

  # create group annotation:
  group_annotation_df <- subset(epithelial_metadata, select = cell_type)
  group_annotation <- Heatmap(
    as.matrix(group_annotation_df), 
    col = cluster_cols, 
    name = "Expression\nclusters", 
    width = unit(4, "mm"), 
    show_row_names = F, show_column_names = F,
    show_heatmap_legend = T
  )
  # create QC annotations:
  nUMI_annotation <- rowAnnotation(
    nUMI = anno_barplot(
      epithelial_metadata$nUMI,
      gp = gpar(
        col = "#D8B72E", 
        width = unit(4, "cm")
      ), 
      border = FALSE, 
      which = "row", 
      axis = F
    ), show_annotation_name = FALSE
  )
  nUMI_annotation@name <- "nUMI"
  nGene_annotation <- rowAnnotation(
    nGene = anno_barplot(
      epithelial_metadata$nGene, name = "nGene",
      gp = gpar(
        col = "#9ECAE1", 
        width = unit(4, "cm")
      ), 
      border = FALSE, 
      which = "row", 
      axis = F
    ), show_annotation_name = FALSE
  )
  nGene_annotation@name <- "nGene"
  
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
    as.matrix(plot_object), name = paste0("hm"), 
    col = heatmap_cols,
    cluster_columns = F, cluster_rows = F,
    show_row_names = F, show_column_names = T,
    column_names_gp = gpar(col = "white"),
    show_row_dend = F,
    show_heatmap_legend = F,
    use_raster = T, raster_device = c("png")
  )
  
  # determine where starting co-ordinates for heatmap are based upon longest cluster name
  # (0.00604 units per character):
  longest_cluster_name <- max(nchar(unique(as.character(epithelial_metadata$cell_type))))
  x_coord <- longest_cluster_name*0.0037
  
  # generate heatmap legend:
  signal_ranges <- round(range(unlist(plot_object)), 2)
  lgd <- Legend(
    at = c(signal_ranges[1], 1, signal_ranges[2]),
    col_fun = heatmap_cols, 
    title = "CNV signal", 
    direction = "horizontal",
    grid_height = unit(2.5, "cm"),
    grid_width = unit(0.1, "cm"),
    labels_gp = gpar(fontsize = 16),
    title_gp = gpar(fontsize = 22, fontface = "plain")
  )
  
  
  ################################################################################
  ### 3. Plot pre-filtered heatmap with expression clusters ###
  ################################################################################
  
  ht_list <- group_annotation + final_heatmap + 
    nUMI_annotation + nGene_annotation
  
  annotated_heatmap <- grid.grabExpr(
    draw(ht_list, gap = unit(6, "mm"), heatmap_legend_side = "left")
  )
  dev.off()
  
  # plot final annotated heatmap:
  png(paste0(plot_dir, "prefiltered_infercnv_plot_with_expression_clusters.png"), 
    height = 13, width = 20, res = 300, units = "in") 
  
    grid.newpage()
      pushViewport(viewport(x = 0.155, y = 0.065, width = 0.752, height = 0.78, 
        just = c("left", "bottom")))
        grid.draw(annotated_heatmap)
        decorate_heatmap_body("hm", {
          for ( e in 1:length(chr_data$end_pos) ) {
          grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), 
            gp = gpar(lwd = 1, col = "#383838"))
          }
        })
      popViewport()
      # plot legend:
      pushViewport(viewport(x = unit(2, "cm"), y = unit(15, "cm"), width = unit(0.1, "cm"), 
        height = unit(0.4, "cm"), just = c("right", "bottom")))
        draw(lgd, x = unit(0.1, "cm"), y = unit(0.1, "cm"), just = c("left", "bottom"))
      popViewport()
      # label annotations:
      pushViewport(viewport(x=x_coord + 0.807, y=0.0001, width = 0.1, height = 0.1, 
        just = "bottom"))
        grid.text("nUMI", rot=65)
      popViewport()
      pushViewport(viewport(x=x_coord + 0.825, y=0.0001, width = 0.1, height = 0.1, 
        just = "bottom"))
        grid.text("nGene", rot=65)
      popViewport()
      
  dev.off()
  
  print(paste0("Heatmap created, output in ", plot_dir))
  
  
  ################################################################################
  ### 4. Plot pre-filtered heatmap without expression clusters ###
  ################################################################################
  
  ht_list <- final_heatmap + nUMI_annotation + nGene_annotation
  
  annotated_heatmap <- grid.grabExpr(
    draw(ht_list, gap = unit(6, "mm"), heatmap_legend_side = "left")
  )
  dev.off()
  
  # plot final annotated heatmap:
  png(paste0(plot_dir, "prefiltered_infercnv_plot_without_expression_clusters.png"), 
    height = 13, width = 20, res = 300, units = "in") 
  
    grid.newpage()
      pushViewport(viewport(x = 0.155, y = 0.065, width = 0.752, height = 0.78, 
        just = c("left", "bottom")))
        grid.draw(annotated_heatmap)
        decorate_heatmap_body("hm", {
          for ( e in 1:length(chr_data$end_pos) ) {
          grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), 
            gp = gpar(lwd = 1, col = "#383838"))
          }
        })
      popViewport()
      # plot legend:
      pushViewport(viewport(x = unit(2, "cm"), y = unit(15, "cm"), width = unit(0.1, "cm"), 
        height = unit(0.4, "cm"), just = c("right", "bottom")))
        draw(lgd, x = unit(0.1, "cm"), y = unit(0.1, "cm"), just = c("left", "bottom"))
      popViewport()
      # label annotations:
      pushViewport(viewport(x=x_coord + 0.807, y=0.0001, width = 0.1, height = 0.1, 
        just = "bottom"))
        grid.text("nUMI", rot=65)
      popViewport()
      pushViewport(viewport(x=x_coord + 0.825, y=0.0001, width = 0.1, height = 0.1, 
        just = "bottom"))
        grid.text("nGene", rot=65)
      popViewport()
      
  dev.off()
  
  print(paste0("Heatmap created, output in ", plot_dir))


}


if (
  !file.exists(
    paste0(
      Robject_dir, "/2b.epithelial_metadata_with_CNV_subclusters.Rdata"
    )
  )
) {

  ################################################################################
  ### 5a. Add tumour subclusters by random tree ###
  ################################################################################

  if (subcluster_by == "random_trees") {

    # load final output object and fetch subcluster data:
    final_infercnv_obj <- readRDS(paste0(in_dir, "run.final.infercnv_obj"))
    subcluster_labels <- final_infercnv_obj@tumor_subclusters$subclusters$all_observations
  
    # rename each cluster:
    names(subcluster_labels) <- paste0("CNV ", 1:length(subcluster_labels))
    
    # match cell ids with subcluster labels and bind in df:
    for (s in 1:length(subcluster_labels)) {
      if (s==1) {
        subcluster_df <- data.frame(
          row.names = names(subcluster_labels[[s]]),
          subcluster_id = rep(
            names(subcluster_labels)[s], 
            length(subcluster_labels[[s]])
          ) 
        )
      } else {
        subcluster_df <- rbind(
          subcluster_df,
          data.frame(
            row.names = names(subcluster_labels[[s]]),
            subcluster_id = rep(
              names(subcluster_labels)[s], 
              length(subcluster_labels[[s]])
            ) 
          )
        )  
      }
    }

  }


#  ######
#  else if (subcluster_by == "Leiden") {
#
#    library(copynumber)
#
#    library(leiden)
#
#    leiden_res <- leiden(epithelial_heatmap)
#
#  }
#  ######
  

  # add to epithelial_metadata:
  epithelial_metadata <- merge(epithelial_metadata, subcluster_df, by="row.names")
  rownames(epithelial_metadata) <- epithelial_metadata$Row.names
  epithelial_metadata <- subset(epithelial_metadata, select = -Row.names)

  # order cells by subcluster:
  epithelial_metadata <- epithelial_metadata[
    order(epithelial_metadata$subcluster_id),
  ]
  epithelial_heatmap <- epithelial_heatmap[rownames(epithelial_metadata),]

  print(paste0(
    "Are epithelial_metadata rownames in the same order as epithelial_heatmap?? ",
    identical(rownames(epithelial_heatmap), rownames(epithelial_metadata))
  ))
  
  saveRDS(epithelial_heatmap, paste0(Robject_dir, 
    "/2a.epithelial_heatmap_with_CNV_subclusters.Rdata"))
  saveRDS(epithelial_metadata, paste0(Robject_dir, 
    "/2b.epithelial_metadata_with_CNV_subclusters.Rdata"))

} else {

  print("Loading heatmap and metadata dfs...")
  epithelial_heatmap <- readRDS(paste0(Robject_dir, 
    "/2a.epithelial_heatmap_with_CNV_subclusters.Rdata"))
  epithelial_metadata <- readRDS(paste0(Robject_dir, 
    "/2b.epithelial_metadata_with_CNV_subclusters.Rdata"))
  print(paste0(
    "Are epithelial_metadata rownames in the same order as epithelial_heatmap?? ",
    identical(rownames(epithelial_heatmap), rownames(epithelial_metadata))
  ))

}


################################################################################
### 6. Add annotations and create heatmap ###
################################################################################

if(!exists("final_infercnv_obj")) {
  final_infercnv_obj <- readRDS(paste0(in_dir, "run.final.infercnv_obj"))
}

if (.hasSlot(final_infercnv_obj, "tumor_subclusters")) {

  if (!file.exists(paste0(
    plot_dir, "prefiltered_infercnv_plot_with_subclone_and_expression_clusters.png"
  ))) {
  
    # create group annotation:
    group_annotation_df <- subset(epithelial_metadata, select = cell_type)
    group_annotation <- Heatmap(
      as.matrix(group_annotation_df), 
      col = cluster_cols, 
      name = "Expression\nclusters", 
      width = unit(4, "mm"), 
      show_row_names = F, show_column_names = F,
      show_heatmap_legend = T
    )
    # create QC annotations:
    nUMI_annotation <- rowAnnotation(
      nUMI = anno_barplot(
        epithelial_metadata$nUMI,
        gp = gpar(
          col = "#D8B72E", 
          width = unit(4, "cm")
        ), 
        border = FALSE, 
        which = "row", 
        axis = F
      ), show_annotation_name = FALSE
    )
    nUMI_annotation@name <- "nUMI"
    nGene_annotation <- rowAnnotation(
      nGene = anno_barplot(
        epithelial_metadata$nGene, name = "nGene",
        gp = gpar(
          col = "#9ECAE1", 
          width = unit(4, "cm")
        ), 
        border = FALSE, 
        which = "row", 
        axis = F
      ), show_annotation_name = FALSE
    )
    nGene_annotation@name <- "nGene"
    
    # create subcluster annotation:
    subcluster_annot_df <- subset(epithelial_metadata, select = subcluster_id)
    subcluster_cols <- rev(col_palette)[
      1:length(unique(subcluster_annot_df$subcluster_id))
    ]
    names(subcluster_cols) <- unique(subcluster_annot_df$subcluster_id)
    subcluster_annotation <- Heatmap(as.matrix(subcluster_annot_df), 
      col = subcluster_cols, 
      name = "subcluster_annotation", width = unit(6, "mm"), 
      show_row_names = F, show_column_names = F, 
      show_heatmap_legend = T,
      heatmap_legend_param = list(
        title = "",
        labels_gp = gpar(fontsize = 16)
      )
    )
  
  #  # write legend function:
  #  subcluster_legend <- function(subclusters, subcluster_cols) {
  #
  #    pushViewport(viewport())
  #
  #    for (i in 1:length(subclusters)) {
  #      # plot legend text:
  #      pushViewport(viewport(x = 0.48, y = 0.05 + (i*.05), 
  #        width = unit(1, "cm"), height = unit(0.5, "cm"), 
  #        just = c("left")))
  #        grid.text(subclusters[i], gp=gpar(fontsize=18))
  #      popViewport()
  #    }
  #
  #  }
  
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
      as.matrix(plot_object), name = paste0("hm"), 
      col = heatmap_cols,
      cluster_columns = F, cluster_rows = F,
      show_row_names = F, show_column_names = T,
      column_names_gp = gpar(col = "white"),
      show_row_dend = F,
      show_heatmap_legend = F,
      use_raster = T, raster_device = c("png")
    )
  
    # determine where starting co-ordinates for heatmap are based upon longest cluster name
    # (0.00604 units per character):
    longest_cluster_name <- max(nchar(unique(as.character(epithelial_metadata$cell_type))))
    x_coord <- longest_cluster_name*0.0037
    
    # generate heatmap legend:
    signal_ranges <- round(range(unlist(plot_object)), 2)
    lgd <- Legend(
      at = c(signal_ranges[1], 1, signal_ranges[2]),
      col_fun = heatmap_cols, 
      title = "CNV signal", 
      direction = "horizontal",
      grid_height = unit(2.5, "cm"),
      grid_width = unit(0.1, "cm"),
      labels_gp = gpar(fontsize = 16),
      title_gp = gpar(fontsize = 22, fontface = "plain")
    )
  
  
    ################################################################################
    ### 7. Generate heatmap ordered by subcluster without expression clusters ###
    ################################################################################
    
    ht_list <- subcluster_annotation + final_heatmap + 
      nUMI_annotation + nGene_annotation
  
    annotated_heatmap <- grid.grabExpr(
      draw(ht_list, gap = unit(6, "mm"), heatmap_legend_side = "left")
    )
    dev.off()
  
    # determine where starting co-ordinates for heatmap are based upon longest cluster name
    # (0.00604 units per character):
    longest_cluster_name <- max(nchar(unique(as.character(epithelial_metadata$cell_type))))
    x_coord <- longest_cluster_name*0.0037
    
    # generate heatmap legend:
    signal_ranges <- round(range(unlist(plot_object)), 2)
    lgd <- Legend(
      at = c(signal_ranges[1], 1, signal_ranges[2]),
      col_fun = heatmap_cols, 
      title = "CNV signal", 
      direction = "horizontal",
      grid_height = unit(2.5, "cm"),
      grid_width = unit(0.1, "cm"),
      labels_gp = gpar(fontsize = 16),
      title_gp = gpar(fontsize = 22, fontface = "plain")
    )
  
    # plot final annotated heatmap:
    png(paste0(
      plot_dir, 
      "prefiltered_infercnv_plot_with_subclone_without_expression_clusters.png"), 
      height = 13, 
      width = 20, 
      res = 300, 
      units = "in"
    ) 
    
      grid.newpage()
        pushViewport(viewport(x = 0.155, y = 0.065, width = 0.752, height = 0.78, 
          just = c("left", "bottom")))
          grid.draw(annotated_heatmap)
          decorate_heatmap_body("hm", {
            for ( e in 1:length(chr_data$end_pos) ) {
            grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), 
              gp = gpar(lwd = 1, col = "#383838"))
            }
          })
        popViewport()
        # plot legend:
        pushViewport(viewport(x = unit(2, "cm"), y = unit(15, "cm"), width = unit(0.1, "cm"), 
          height = unit(0.4, "cm"), just = c("right", "bottom")))
          draw(lgd, x = unit(0.1, "cm"), y = unit(0.1, "cm"), just = c("left", "bottom"))
        popViewport()
        # label annotations:
        pushViewport(viewport(x=x_coord + 0.807, y=0.0001, width = 0.1, height = 0.1, 
          just = "bottom"))
          grid.text("nUMI", rot=65)
        popViewport()
        pushViewport(viewport(x=x_coord + 0.825, y=0.0001, width = 0.1, height = 0.1, 
          just = "bottom"))
          grid.text("nGene", rot=65)
        popViewport()
        
    dev.off()
    
    print(paste0("Heatmap created, output in ", plot_dir))
  
  
    ################################################################################
    ### 8. Generate heatmap ordered by subcluster with expression clusters ###
    ################################################################################
    
    ht_list <- subcluster_annotation + group_annotation + final_heatmap + 
      nUMI_annotation + nGene_annotation
  
    annotated_heatmap <- grid.grabExpr(
      draw(ht_list, gap = unit(6, "mm"), heatmap_legend_side = "left")
    )
    dev.off()
    
    # plot final annotated heatmap:
    png(paste0(
      plot_dir, 
      "prefiltered_infercnv_plot_with_subclone_and_expression_clusters.png"), 
      height = 13, 
      width = 20, 
      res = 300, 
      units = "in"
    ) 
    
      grid.newpage()
        pushViewport(viewport(x = 0.155, y = 0.065, width = 0.752, height = 0.78, 
          just = c("left", "bottom")))
          grid.draw(annotated_heatmap)
          decorate_heatmap_body("hm", {
            for ( e in 1:length(chr_data$end_pos) ) {
            grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), 
              gp = gpar(lwd = 1, col = "#383838"))
            }
          })
        popViewport()
        # plot legend:
        pushViewport(viewport(x = unit(2, "cm"), y = unit(15, "cm"), width = unit(0.1, "cm"), 
          height = unit(0.4, "cm"), just = c("right", "bottom")))
          draw(lgd, x = unit(0.1, "cm"), y = unit(0.1, "cm"), just = c("left", "bottom"))
        popViewport()
        # label annotations:
        pushViewport(viewport(x=x_coord + 0.807, y=0.0001, width = 0.1, height = 0.1, 
          just = "bottom"))
          grid.text("nUMI", rot=65)
        popViewport()
        pushViewport(viewport(x=x_coord + 0.825, y=0.0001, width = 0.1, height = 0.1, 
          just = "bottom"))
          grid.text("nGene", rot=65)
        popViewport()
        
    dev.off()
    
    print(paste0("Heatmap created, output in ", plot_dir))
  
  }

}


################################################################################
### 9. Remove cells below nUMI/nGene threshold ###
################################################################################

if (!file.exists(paste0(Robject_dir, "/3a.epithelial_heatmap_good_coverage_only.Rdata"))) {

  print(paste0(
    "Number of cells before filtering for good UMI/gene coverage = ", 
    nrow(epithelial_metadata)
  ))

  filtered_metadata <- epithelial_metadata[
    epithelial_metadata$nUMI > nUMI_threshold & epithelial_metadata$nUMI > nGene_threshold,
  ]

  if (.hasSlot(final_infercnv_obj, "tumor_subclusters")) {
    # throw subclusters with number of cells < (min_proportion_cells_for_subcluster*total cell number):
    split_metadata <- split(
      filtered_metadata, filtered_metadata$subcluster_id
    )
    subcluster_cell_no <- lapply(split_metadata, nrow)
    to_remove <- split_metadata[
      subcluster_cell_no < (min_proportion_cells_for_subcluster*sum(unlist(subcluster_cell_no)))
    ]
    cells_to_remove <- as.character(do.call("rbind", to_remove)$cell_ids)
    filtered_metadata <- filtered_metadata[!(filtered_metadata$cell_ids %in% cells_to_remove),]
  }

  filtered_heatmap <- epithelial_heatmap[
    rownames(epithelial_heatmap) %in% filtered_metadata$cell_ids,
  ]

  print(paste0(
    "Number of cells after filtering for good UMI/gene coverage = ", 
    nrow(filtered_metadata)
  ))

  print(paste0(
    "Are filtered_metadata rownames in the same order as filtered_heatmap?? ",
    identical(rownames(filtered_heatmap), rownames(filtered_metadata))
  ))
  
  saveRDS(filtered_heatmap, paste0(Robject_dir, 
    "/3a.filtered_heatmap_good_coverage_only.Rdata"))
  saveRDS(filtered_metadata, paste0(Robject_dir, 
    "/3b.filtered_metadata_good_coverage_only.Rdata"))

} else {

  print("Loading heatmap and metadata dfs...")
  filtered_heatmap <- readRDS(paste0(Robject_dir, 
    "/3a.filtered_heatmap_good_coverage_only.Rdata"))
  filtered_metadata <- readRDS(paste0(Robject_dir, 
    "/3b.filtered_metadata_good_coverage_only.Rdata"))

  print(paste0(
    "Are filtered_metadata rownames in the same order as filtered_heatmap? ",
    identical(rownames(filtered_heatmap), rownames(filtered_metadata))
  ))

}


################################################################################
### 8. Add annotations and prepare heatmap ###
################################################################################

if (!file.exists(paste0(
    plot_dir, "infercnv_plot_with_expression_clusters.png"
  ))) {

  # create group annotation:
  group_annotation_df <- subset(filtered_metadata, select = cell_type)
  group_annotation <- Heatmap(
    as.matrix(group_annotation_df), 
    col = cluster_cols, 
    name = "Expression\nclusters", 
    width = unit(4, "mm"), 
    show_row_names = F, show_column_names = F,
    show_heatmap_legend = T
  )
  # create QC annotations:
  nUMI_annotation <- rowAnnotation(
    nUMI = anno_barplot(
      filtered_metadata$nUMI,
      gp = gpar(
        col = "#D8B72E", 
        width = unit(4, "cm")
      ), 
      border = FALSE, 
      which = "row", 
      axis = F
    ), show_annotation_name = FALSE
  )
  nUMI_annotation@name <- "nUMI"
  nGene_annotation <- rowAnnotation(
    nGene = anno_barplot(
      filtered_metadata$nGene, name = "nGene",
      gp = gpar(
        col = "#9ECAE1", 
        width = unit(4, "cm")
      ), 
      border = FALSE, 
      which = "row", 
      axis = F
    ), show_annotation_name = FALSE
  )
  nGene_annotation@name <- "nGene"
  
  if (.hasSlot(final_infercnv_obj, "tumor_subclusters")) {
    # create subcluster annotation:
    subcluster_annot_df <- subset(filtered_metadata, select = subcluster_id)
    subcluster_cols <- rev(col_palette)[
      1:length(unique(subcluster_annot_df$subcluster_id))
    ]
    names(subcluster_cols) <- unique(subcluster_annot_df$subcluster_id)
    subcluster_annotation <- Heatmap(as.matrix(subcluster_annot_df), 
      col = subcluster_cols, 
      name = "subcluster_annotation", width = unit(6, "mm"), 
      show_row_names = F, show_column_names = F, 
      show_heatmap_legend = T,
      heatmap_legend_param = list(
        title = "",
        labels_gp = gpar(fontsize = 16)
      )
    )
  }
  
  # fetch chromosome boundary co-ordinates:
  if (!file.exists(paste0(Robject_dir, "chromosome_data.Rdata"))) {
    chr_data <- fetch_chromosome_boundaries(filtered_heatmap, ref_dir)
    saveRDS(chr_data, paste0(Robject_dir, "chromosome_data.Rdata"))
  } else {
    chr_data <- readRDS(paste0(Robject_dir, "chromosome_data.Rdata"))
  }
  
  # prepare df for plotting:
  plot_object <- filtered_heatmap
  colnames(plot_object) <- rep("la", ncol(plot_object))
  # define heatmap colours:
  na_less_vector <- unlist(plot_object)
  na_less_vector <- na_less_vector[!is.na(na_less_vector)]
  heatmap_cols <- colorRamp2(c(min(na_less_vector), 1, max(na_less_vector)), 
        c("#00106B", "white", "#680700"), space = "sRGB")
  
  print("Generating final heatmap...")
  # create main CNV heatmap:
  final_heatmap <- Heatmap(
    as.matrix(plot_object), name = paste0("hm"), 
    col = heatmap_cols,
    cluster_columns = F, cluster_rows = F,
    show_row_names = F, show_column_names = T,
    column_names_gp = gpar(col = "white"),
    show_row_dend = F,
    show_heatmap_legend = F,
    heatmap_legend_param = list(labels_gp = gpar(col = "red", fontsize = 12)),
    use_raster = T, raster_device = c("png")
  )
  
  
  ################################################################################
  ### 8. Plot heatmap with expression clusters ###
  ################################################################################
  
  if (.hasSlot(final_infercnv_obj, "tumor_subclusters")) {
    ht_list <- subcluster_annotation + group_annotation + final_heatmap + 
      nUMI_annotation + nGene_annotation
  } else {
    ht_list <- group_annotation + final_heatmap + 
      nUMI_annotation + nGene_annotation
  }
  annotated_heatmap <- grid.grabExpr(
    draw(ht_list, gap = unit(6, "mm"), heatmap_legend_side = "left")
  )
  dev.off()
  
  # determine where starting co-ordinates for heatmap are based upon longest cluster name
  # (0.00604 units per character):
  longest_cluster_name <- max(nchar(unique(as.character(filtered_metadata$cell_type))))
  x_coord <- longest_cluster_name*0.0037
  
  # generate heatmap legend:
  signal_ranges <- round(range(unlist(plot_object)), 2)
  lgd <- Legend(
    at = c(signal_ranges[1], 1, signal_ranges[2]),
    col_fun = heatmap_cols, 
    title = "CNV signal", 
    direction = "horizontal",
    grid_height = unit(2.5, "cm"),
    grid_width = unit(0.1, "cm"),
    labels_gp = gpar(fontsize = 16),
    title_gp = gpar(fontsize = 22, fontface = "plain")
  )
  
  # plot final annotated heatmap:
  png(paste0(plot_dir, "infercnv_plot_with_expression_clusters.png"), 
    height = 13, width = 20, res = 300, units = "in") 
  
    grid.newpage()
      pushViewport(viewport(x = 0.155, y = 0.065, width = 0.752, height = 0.78, 
        just = c("left", "bottom")))
        grid.draw(annotated_heatmap)
        decorate_heatmap_body("hm", {
          for ( e in 1:length(chr_data$end_pos) ) {
            grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), 
              gp = gpar(lwd = 3, col = "#383838"))
            if (e==1) {
              grid.text(gsub("chr", "", names(chr_data$lab_pos)[e]), chr_data$lab_pos[e], 
              unit(0, "npc") + unit(-3.5, "mm"), gp=gpar(fontsize=20))
            } else if (e==21) {
              grid.text(paste0("\n", gsub("chr", "", names(chr_data$lab_pos)[e])), chr_data$lab_pos[e], 
              unit(0, "npc") + unit(-3.5, "mm"), gp=gpar(fontsize=20))
            } else {
              grid.text(gsub("chr", "", names(chr_data$lab_pos)[e]), chr_data$lab_pos[e], 
              unit(0, "npc") + unit(-3.5, "mm"), gp=gpar(fontsize=20))        
            }
          }
        })
      popViewport()
      # plot legend:
      pushViewport(viewport(x = unit(2, "cm"), y = unit(15, "cm"), width = unit(0.1, "cm"), 
        height = unit(0.4, "cm"), just = c("right", "bottom")))
        draw(lgd, x = unit(0.1, "cm"), y = unit(0.1, "cm"), just = c("left", "bottom"))
      popViewport()
      # label annotations:
      pushViewport(viewport(x=x_coord + 0.807, y=0.0001, width = 0.1, height = 0.1, 
        just = "bottom"))
        grid.text("nUMI", rot=65)
      popViewport()
      pushViewport(viewport(x=x_coord + 0.825, y=0.0001, width = 0.1, height = 0.1, 
        just = "bottom"))
        grid.text("nGene", rot=65)
      popViewport()
      
  dev.off()
  
  print(paste0("Heatmap created, output in ", plot_dir))
  
  
  ################################################################################
  ### 8. Plot heatmap without expression clusters ###
  ################################################################################
  
  if (.hasSlot(final_infercnv_obj, "tumor_subclusters")) {
    ht_list <- subcluster_annotation + final_heatmap + 
      nUMI_annotation + nGene_annotation
  } else {
    ht_list <- group_annotation + final_heatmap + 
      nUMI_annotation + nGene_annotation
  }
  annotated_heatmap <- grid.grabExpr(
    draw(ht_list, gap = unit(6, "mm"), heatmap_legend_side = "left")
  )
  dev.off()
  
  # determine where starting co-ordinates for heatmap are based upon longest cluster name
  # (0.00604 units per character):
  longest_cluster_name <- max(nchar(unique(as.character(filtered_metadata$cell_type))))
  x_coord <- longest_cluster_name*0.0037
  
  # generate heatmap legend:
  signal_ranges <- round(range(unlist(plot_object)), 2)
  lgd <- Legend(
    at = c(signal_ranges[1], 1, signal_ranges[2]),
    col_fun = heatmap_cols, 
    title = "CNV signal", 
    direction = "horizontal",
    grid_height = unit(2.5, "cm"),
    grid_width = unit(0.1, "cm"),
    labels_gp = gpar(fontsize = 16),
    title_gp = gpar(fontsize = 22, fontface = "plain")
  )
  
  # plot final annotated heatmap:
  png(paste0(plot_dir, "infercnv_plot_without_expression_clusters.png"), 
    height = 13, width = 20, res = 300, units = "in") 
  
    grid.newpage()
      pushViewport(viewport(x = 0.155, y = 0.065, width = 0.752, height = 0.78, 
        just = c("left", "bottom")))
        grid.draw(annotated_heatmap)
        decorate_heatmap_body("hm", {
          for ( e in 1:length(chr_data$end_pos) ) {
            grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), 
              gp = gpar(lwd = 3, col = "#383838"))
            if (e==1) {
              grid.text(gsub("chr", "", names(chr_data$lab_pos)[e]), chr_data$lab_pos[e], 
              unit(0, "npc") + unit(-3.5, "mm"), gp=gpar(fontsize=20))
            } else if (e==21) {
              grid.text(paste0("\n", gsub("chr", "", names(chr_data$lab_pos)[e])), chr_data$lab_pos[e], 
              unit(0, "npc") + unit(-3.5, "mm"), gp=gpar(fontsize=20))
            } else {
              grid.text(gsub("chr", "", names(chr_data$lab_pos)[e]), chr_data$lab_pos[e], 
              unit(0, "npc") + unit(-3.5, "mm"), gp=gpar(fontsize=20))        
            }
          }
        })
      popViewport()
      # plot legend:
      pushViewport(viewport(x = unit(2, "cm"), y = unit(15, "cm"), width = unit(0.1, "cm"), 
        height = unit(0.4, "cm"), just = c("right", "bottom")))
        draw(lgd, x = unit(0.1, "cm"), y = unit(0.1, "cm"), just = c("left", "bottom"))
      popViewport()
      # label annotations:
      pushViewport(viewport(x=x_coord + 0.807, y=0.0001, width = 0.1, height = 0.1, 
        just = "bottom"))
        grid.text("nUMI", rot=65)
      popViewport()
      pushViewport(viewport(x=x_coord + 0.825, y=0.0001, width = 0.1, height = 0.1, 
        just = "bottom"))
        grid.text("nGene", rot=65)
      popViewport()
      
  dev.off()
  
  print(paste0("Heatmap created, output in ", plot_dir))

}
  
  