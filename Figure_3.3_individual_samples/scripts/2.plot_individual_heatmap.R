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
#sample_name <- "CID4515"
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

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
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
#integrated_dir <- paste0("/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/",
#	"Jun2019/04_reclustering_analysis/run06_v1.2.1/output/Epithelial/02_Rdata/")
in_dir <- paste0(results_dir, "infercnv/CID4515/")
input_dir <- paste0(results_dir, "infercnv/CID4515/input_files/")

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


#################################################################################
### 1. Load all cells seurat object and create tSNEs/UMAPs ###
#################################################################################

if (!file.exists(paste0(plot_dir, "epithelial_feature.png"))) {

  seurat_10X <- readRDS(paste0(seurat_dir, "03_seurat_object_processed.Rdata"))
  Idents(seurat_10X) <- paste0(
    seurat_10X@meta.data$garnett_seurat_cluster_call_major_PC_C_res.0.6,
    " ",
    seurat_10X@meta.data$PC_C_res.0.6
  )
  
  tSNE <- DimPlot(
    object = seurat_10X,
    label.size = 70,
    pt.size = 2,
    reduction = paste0("TSNE", all_cell_PC),
     order = T
  )
  png(
    paste0(plot_dir, all_cell_res, "_tsne.png"),
     width = 10, 
     height = 8, 
     res = 300, 
     units = 'in'
  )
    print(tSNE)
  dev.off()
  
  UMAP <- DimPlot(
    object = seurat_10X,
    label.size = 70,
    pt.size = 2,
    reduction = paste0("UMAP", all_cell_PC),
    order = T
  )
  
  png(
    paste0(plot_dir, all_cell_res, "_umap.png"),
     width = 10, 
     height = 8, 
     res = 300, 
     units = 'in'
  )
    print(UMAP)
  dev.off()
  
  broad_feature <- FeaturePlot(
    object = seurat_10X,
    features = broad_markers,
    pt.size = 1.5,
    reduction = paste0("TSNE", all_cell_PC),
    order = T
  )
  dev.off()
  png(
    paste0(plot_dir, "broad_feature.png"),
     width = 10, 
     height = 8, 
     res = 300, 
     units = 'in'
  )
    print(broad_feature)
  dev.off()
  
  epithelial_feature <- FeaturePlot(
    object = seurat_10X,
    features = epithelial_markers,
    pt.size = 1.5,
    reduction = paste0("TSNE", all_cell_PC),
    order = T
  )
  dev.off()
  png(
    paste0(plot_dir, "epithelial_feature.png"),
     width = 10, 
     height = 8, 
     res = 300, 
     units = 'in'
  )
    print(epithelial_feature)
  dev.off()

}


#################################################################################
### 1. Load malignant epithelial seurat object and create tSNEs/UMAPs to choose 
# best PC/resolution ###
#################################################################################

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


if (!file.exists(paste0(plot_dir, malignant_res, "_epithelial_tsne.png"))) {

  # isolate and plot epithelial cells:
  epi_tSNE <- DimPlot(
    object = seurat_malignant,
    label.size = 40,
    pt.size = 2,
    reduction = paste0("TSNE", malignant_PC),
    cols = cluster_cols
  )
  
  png(
    paste0(plot_dir, malignant_res, "_epithelial_tsne.png"),
     width = 10, 
     height = 8, 
     res = 300, 
     units = 'in'
  )
    print(epi_tSNE)
  dev.off()
  
  epi_UMAP <- DimPlot(
    object = seurat_malignant,
    label.size = 40,
    pt.size = 2,
    reduction = paste0("UMAP", malignant_PC),
    cols = cluster_cols
  )
  
  png(
    paste0(plot_dir, malignant_res, "_epithelial_umap.png"),
     width = 10, 
     height = 8, 
     res = 300, 
     units = 'in'
  )
    print(epi_UMAP)
  dev.off()

}

if (!file.exists(paste0(plot_dir, "expression_cluster_DGE_heatmap.png"))) {

  epi_DE <- FindAllMarkers(
    only.pos = T,
    object = seurat_malignant,
    min.pct = 0.5, 
    logfc.threshold = 0.7, 
    test.use = 'MAST'
  )
  
  epi_DE_sorted <- arrange(epi_DE, cluster, desc(avg_logFC))
  
  write.table(epi_DE_sorted, paste0(table_dir, "expression_cluster_DGE.txt"))
  
  heatmap_genes <- epi_DE_sorted %>% 
    group_by(cluster) %>% 
    top_n(10, avg_logFC)
  
  png(
    paste0(plot_dir, "expression_cluster_DGE_heatmap.png"),
    width = 12,
    height = 9,
    res = 300,
    units = "in"
  )
    print(DoHeatmap(
      seurat_malignant,
      features = heatmap_genes$gene,
      group.by = "ident",
      group.colors = cluster_cols,
      size = 3
    ))
  dev.off()
  
}


################################################################################
### 2. Load InferCNV output and create heatmap and metadata dfs ###
################################################################################

if (!file.exists(paste0(Robject_dir, "/1b.initial_epithelial_metadata.Rdata"))) {

  # load InferCNV output:
  print("Loading InferCNV output files...")
  epithelial_heatmap <- as.data.frame(t(read.table(paste0(in_dir, 
  	"infercnv.12_denoised.observations.txt"))))

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

}


################################################################################
### 3. Add tumour subclusters ###
################################################################################

if (
  !file.exists(
  	paste0(
  	  Robject_dir, "/2b.epithelial_metadata_with_CNV_subclusters.Rdata"
  	)
  )
) {

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
### 4. Add annotations ###
################################################################################

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
  show_heatmap_legend = F
)


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

ht_list <- subcluster_annotation + group_annotation + final_heatmap + 
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
#png(paste0(plot_dir, "legend_test.png"), 
#  height = 13, width = 20, res = 300, units = "in") 
#  grid.newpage()
#    pushViewport(viewport(x = unit(2, "cm"), y = unit(15, "cm"), width = unit(0.5, "cm"), height = unit(1, "cm"), 
#      just = c("right", "bottom")))
#      draw(lgd, x = unit(0.1, "cm"), y = unit(0.1, "cm"), just = c("left", "bottom"))
#    popViewport()
#dev.off()


################################################################################
### 6. Plot heatmap ###
################################################################################

# plot final annotated heatmap:
png(paste0(plot_dir, "infercnv_plot.png"), 
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

pdf(paste0(plot_dir, "infercnv_plot.pdf"), 
  height = 13, width = 20) 

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

# print epithelial_metadata as table:
write.table(epithelial_metadata, paste0(table_dir, "epithelial_metadata.txt"), 
  sep="\t", quote=F, row.names=F, col.name=T)

#convert pdf to png:
system(paste0("for p in ", plot_dir, "*.pdf; do echo $p; f=$(basename $p); echo $f; ",
"new=$(echo $f | sed 's/.pdf/.png/'); echo $new; ", 
"convert -density 150 ", plot_dir, "$f -quality 90 ", plot_dir, "$new; done"))

