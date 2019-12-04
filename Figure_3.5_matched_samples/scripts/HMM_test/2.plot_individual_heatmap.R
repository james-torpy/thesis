#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

#### Generate CNV heatmaps with following annotations: ###
# nUMI
# nGene
# CNA values
# SNP array CNVs

project_name <- "thesis"
subproject_name <- "Figure_3"
args = commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
include_t_cells <- as.logical(args[2])

#sample_name <- "CID45171_epithelial_clustering_supervised"
#include_t_cells <- TRUE

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
  in_dir <- paste0(results_dir, "infercnv/t_cells_included/", 
    sample_name, "/")
} else {
  in_dir <- paste0(results_dir, "infercnv/t_cells_excluded/", 
    sample_name, "/")
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
create_array_CNV_annotation <- dget(paste0(func_dir, 
  "create_array_CNV_annotation.R"))


################################################################################
### 1. Load InferCNV output and create heatmap and metadata dfs ###
################################################################################

if (!file.exists(paste0(Robject_dir, "/epithelial_HMM_df.Rdata")) | 
	!file.exists(paste0(Robject_dir, "/orig_metadata_HMM_df.Rdata"))) {
  # load InferCNV output:
  print("Loading InferCNV output files...")
  infercnv_output <- as.data.frame(t(read.table(paste0(in_dir, 
  	"infercnv.13_HMM_pred.Bayes_Net.Pnorm_0.5.observations.txt"))))

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
  heatmap_df <- infercnv_output[rownames(infercnv_output) %in% epithelial_ids,]
  print(paste0("Number of heatmap rows after non-epithelial thrown: ", 
  	nrow(heatmap_df)))

  saveRDS(heatmap_df, paste0(Robject_dir, "/epithelial_HMM_df.Rdata"))
  saveRDS(metadata_df, paste0(Robject_dir, "/orig_metadata_HMM_df.Rdata"))
} else {
	print("Loading heatmap and metadata dfs...")
	heatmap_df <- readRDS(paste0(Robject_dir, "/epithelial_HMM_df.Rdata"))
	metadata_df <- readRDS(paste0(Robject_dir, "/orig_metadata_HMM_df.Rdata"))
}


################################################################################
### 2. Add annotation metadata ###
################################################################################

if (!file.exists(paste0(Robject_dir, 
  "epithelial_HMM_metadata_with_cell_type_and_QC.Rdata"))) {
  # create epithelial_metadata df and only include epithelial cells in heatmap_df:
  print("Creating epithelial metadata df...")
  epithelial_metadata <- metadata_df[rownames(heatmap_df),]
  # order epithelial metadata cell type cluster levels:
  epithelial_metadata$cell_type <- factor(
    epithelial_metadata$cell_type,
    levels = naturalsort(unique(epithelial_metadata$cell_type))
  )
  print(paste0(
    "Are epithelial_metadata rownames in the same order as heatmap_df? ",
    identical(rownames(heatmap_df), rownames(epithelial_metadata))
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
    "Are epithelial_metadata rownames still in the same order as heatmap_df?? ",
    identical(rownames(heatmap_df), rownames(epithelial_metadata))
  ))
  saveRDS(
    epithelial_metadata, paste0(Robject_dir, 
    "epithelial_HMM_metadata_with_cell_type_and_QC.Rdata")
  )
} else {
  epithelial_metadata <- readRDS(
    paste0(Robject_dir, 
    "epithelial_HMM_metadata_with_cell_type_and_QC.Rdata")
  )
}
if (!file.exists(paste0(Robject_dir, 
  "epithelial_HMM_metadata_with_cell_type_QC_GIN_and_correlation_values.Rdata"))) {
  
  # determine GIN and add to epithelial_metadata:
  print("Determining CNA values and adding to epithelial metadata df...")
  GIN <- apply(heatmap_df, 1, function(y), mean(y^2))
  GIN_df <- data.frame(
    row.names = names(GIN),
    GIN = GIN
  )
  epithelial_metadata <- cbind(epithelial_metadata, GIN_df)
  print(paste0(
    "Are epithelial_metadata rownames still in the same order as heatmap_df?? ",
    identical(rownames(heatmap_df), rownames(epithelial_metadata))
  ))
    # determine correlation with top 5% cancer values and add to epithelial_metadata:
  print(paste0(
    "Determining correlation with top 5% cancer values and adding to epithelial ", 
    "metadata df..."
  ))
  # determine top 5% cancer cells:
  GIN_order <- order(GIN_df$GIN, decreasing=T)
  ordered_GIN  <- data.frame(
    row.names = rownames(GIN_df)[GIN_order],
    GIN = GIN_df[GIN_order,]
  )
  top_cancer <- head(ordered_GIN, nrow(ordered_GIN)*0.05)
  
  # find average genome-wide CNV predictions across genome:
  top_cancer_GIN_average <- apply(heatmap_df[rownames(top_cancer),], 2, mean)
  # find correlations of each cell's CNVs with top_GIN_CNV_average:
  cancer_correlations <- apply(heatmap_df, 1, function(x) {
    if (length(unique(as.numeric(x))) == 1) {
      cor_result <- data.frame(cor.estimate="no_CNVs_recorded", 
        cor.p.value="no_CNVs_recorded")
    } else {
      cor <- cor.test(as.numeric(x), top_cancer_GIN_average, method = "kendall")
      cor_result <- data.frame(cor$estimate, cor$p.value)
    }
    return(cor_result)
  })
  correlation_df <- do.call("rbind", cancer_correlations)

  # add to epithelial_metadata:
  epithelial_metadata <- cbind(epithelial_metadata, correlation_df)
  print(paste0(
    "Are epithelial_metadata rownames still in the same order as heatmap_df?? ",
    identical(rownames(heatmap_df), rownames(epithelial_metadata))
  ))
  saveRDS(
    epithelial_metadata, paste0(Robject_dir, 
    "epithelial_HMM_metadata_with_cell_type_QC_GIN_and_correlation_values.Rdata")
  )
} else {
  epithelial_metadata <- readRDS(
    paste0(Robject_dir, 
    "epithelial_HMM_metadata_with_cell_type_QC_GIN_and_correlation_values.Rdata")
  )
}


################################################################################
### 3. Remove clusters without enough epithelial cells ###
################################################################################

if (!file.exists(paste0(Robject_dir, 
    "epithelial_HMM_metadata_filtered_with_cell_type_QC_CNA_and_correlation_values.Rdata"))) {
  print("Removing non-epithelial cluster cells from heatmap and metadata dfs...")
  
  if (!exists("seurat_10X")) {
    seurat_10X <- readRDS(
      paste0(seurat_dir, "03_seurat_object_processed.Rdata")
    )
  }
  Idents(seurat_10X) <- seurat_10X@meta.data$PC_A_res.1
  cluster_ids <- as.character(unique(epithelial_metadata$cell_type))
  for (l in 1:length(cluster_ids)) {
    # fetch garnett calls for cluster id no:
    simple_id <- gsub("^.*_", "", cluster_ids[l])
    garnett_calls <- table(
      seurat_10X@meta.data$garnett_call_ext_major[
        Idents(seurat_10X) == simple_id
      ]
    )
    # if 1.5x or more non-epithelial than epithelial cells, add to remove list:
    non_epithelial <- garnett_calls[grep("pithelial", names(garnett_calls), 
      invert=T)]
    epithelial <- garnett_calls[grep("pithelial", names(garnett_calls))]
    if (any(non_epithelial > (1.5*epithelial))) {
      if (!exists("remove_clusters")) {
        remove_clusters <- c(cluster_ids[l])
      } else {
        remove_clusters <- c(remove_clusters, cluster_ids[l])
      }
    }
  }

  # remove cells attributed to clusters to be removed:
  if (exists("remove_clusters")) {
    remove_cells <- epithelial_metadata$cell_ids[
      epithelial_metadata$cell_type %in% remove_clusters
    ]
    print(paste0("No. epithelial cells to be removed = ", length(remove_cells)))
    print(paste0("No. rows in epithelial metadata df before removing ",
      "non-epithelial cluster cells = ", nrow(epithelial_metadata)))
    epithelial_metadata <- epithelial_metadata[
      !(epithelial_metadata$cell_ids %in% remove_cells),
    ]
    print(paste0("No. rows in epithelial metadata df after removing ",
    "non-epithelial cluster cells = ", nrow(epithelial_metadata)))
  
    heatmap_df <- heatmap_df[!(rownames(heatmap_df) %in% remove_cells),]
    print(paste0(
      "Are epithelial_metadata rownames in the same order as heatmap_df?? ",
      identical(rownames(heatmap_df), rownames(epithelial_metadata))
    ))
  } else {
    print(paste0(
      "Are epithelial_metadata rownames in the same order as heatmap_df?? ",
      identical(rownames(heatmap_df), rownames(epithelial_metadata))
    ))
  }
  saveRDS(epithelial_metadata,  paste0(Robject_dir, 
    "epithelial_metadata_filtered_with_cell_type_QC_CNA_and_correlation_values.Rdata"))
} else {
  epithelial_metadata <- readRDS(
    paste0(Robject_dir, 
    "epithelial_HMM_metadata_filtered_with_cell_type_QC_CNA_and_correlation_values.Rdata")
  )
}


################################################################################
### 4. Add sc50 calls to metadata ###
################################################################################

if (!file.exists(paste0(Robject_dir, "epithelial_HMM_metadata_final.Rdata"))) {
  print(paste0(
    "Adding sc50 calls to epithelial metadata df..."
  ))
  # read in data and select data for current sample:
  sc50 <- read.table(paste0(ref_dir, "/sc50_calls.txt"), header = T, sep = "\t",
    stringsAsFactors = F)
  sc50$SampleID <- gsub("4290", "4290A", sc50$SampleID)
  sample_sc50 <- sc50[grep(sample_name, sc50$SampleID),]

  if (length(sample_sc50$Sept.SC50_Calls) > 0) {

    # select only cells contained in epithelial_metadata:
    print(paste0("Original sc50 call no. for ", sample_name, " = ", 
      nrow(sample_sc50)))
    rownames(sample_sc50) <- sample_sc50$SampleID
    sample_sc50 <- sample_sc50[rownames(epithelial_metadata),]
     print(paste0("sc50 call no. for ", sample_name, 
      " after selecting only cells in epithelial metadata df = ", nrow(sample_sc50)))

    # cbind to epithelial_metadata:
    epithelial_metadata <- cbind(epithelial_metadata, sample_sc50$Sept.SC50_Calls)
    colnames(epithelial_metadata)[ncol(epithelial_metadata)] <- "sc50"

  } else {

    print(paste0(
      "No sc50 calls found for sample ", sample_name, ", adding NAs for this column..."
    ))
    epithelial_metadata$sc50 <- NA
  }

  print(paste0(
    "Are epithelial_metadata rownames still in the same order as heatmap_df?? ",
    identical(rownames(heatmap_df), rownames(epithelial_metadata))
  ))

  saveRDS(epithelial_metadata, paste0(Robject_dir, "epithelial_HMM_metadata_final.Rdata"))
} else {
  epithelial_metadata <- readRDS(paste0(Robject_dir, "epithelial_HMM_metadata_final.Rdata"))
}

heatmap_metadata <- epithelial_metadata

saveRDS(heatmap_df, paste0(Robject_dir, "heatmap_HMM_df.Rdata"))
saveRDS(heatmap_metadata, paste0(Robject_dir, "heatmap_HMM_metadata.Rdata"))


################################################################################
### 5. Order heatmap, metadata and create heatmap annotations ###
################################################################################

# define group annotation colours:
extra_colours <- c(brewer.pal(12, "Set3"), brewer.pal(12, "Paired"))
col_palette <- c(brewer.pal(8, "Dark2")[c(3:6,8)], brewer.pal(12, "Set3"), brewer.pal(8, "Accent"),
  "#660200", "#918940", "black", "#9ECAE1", "#3F007D", "#67000D", "#FDD0A2", "#08519C",
  "#DBB335", "#5AA050", "#807DBA", "#1A6000", "#F16913", "#FD8D3C", "#DEEBF7", 
  "#7F2704", "#DADAEB", "#FC9272", "#BCBDDC", extra_colours, "#b2182b", "#85929E", 
  "#9B59B6", "#74add1","#1b7837", "#b8e186", "#fed976","#e7298a", "#18ffff", "#ef6c00",
  "#A93226", "black","orange", "#b8bc53", "#5628ce", "#fa909c", "#8ff331","#270e26")
cluster_number <- length(unique(heatmap_metadata$cell_type))
cluster_cols <- col_palette[1:cluster_number]
names(cluster_cols) <- unique(heatmap_metadata$cell_type)
cluster_cols <- cluster_cols[levels(heatmap_metadata$cell_type)]

# create group annotation:
group_annotation_df <- subset(heatmap_metadata, select = cell_type)
group_annotation <- Heatmap(
  group_annotation_df, 
  col = cluster_cols, 
  name = "group_annotation", 
  width = unit(4, "mm"), 
  show_row_names = F, show_column_names = F,
  show_heatmap_legend = F
)
# create GIN annotation:
GIN_annotation <- rowAnnotation(
  correlation_annotation = anno_barplot(
    heatmap_metadata$GIN,
    gp = gpar(
      col = "#D95F02", 
      width = unit(4, "cm")
    ), 
    border = FALSE, 
    which = "row", 
    axis = F
  )
)
GIN_annotation@name <- "GIN"
# create QC annotations:
nUMI_annotation <- rowAnnotation(
  correlation_annotation = anno_barplot(
    heatmap_metadata$nUMI,
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
    heatmap_metadata$nGene, name = "nGene",
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

# create array CNV annotation:
all_array_CNVs <- read.table(paste0(ref_dir, "all_array_CNVs.txt"))
colnames(all_array_CNVs) <- gsub("CID4499_1", "CID44991", colnames(all_array_CNVs))
if (any(colnames(all_array_CNVs) %in% sample_name)) {
  if (!file.exists(paste0(Robject_dir, "array_CNV_annotation.Rdata"))) {
    array_CNV_annotation <- create_array_CNV_annotation(heatmap_df, all_array_CNVs)
    saveRDS(array_CNV_annotation, paste0(Robject_dir, "array_CNV_annotation.Rdata"))
    grid_array_heatmap <- grid.grabExpr(draw(array_CNV_annotation$array_CNV_heatmap, 
      heatmap_legend_side = "left"))
  } else {
    array_CNV_annotation <- readRDS(paste0(Robject_dir, "array_CNV_annotation.Rdata"))
    grid_array_heatmap <- grid.grabExpr(draw(array_CNV_annotation$array_CNV_heatmap, 
      heatmap_legend_side = "left"))
  }
}


################################################################################
### 6. Generate heatmap ###
################################################################################

# fetch chromosome boundary co-ordinates:
if (!file.exists(paste0(Robject_dir, "chromosome_data.Rdata"))) {
  chr_data <- fetch_chromosome_boundaries(heatmap_df, ref_dir)
  saveRDS(chr_data, paste0(Robject_dir, "chromosome_data.Rdata"))
} else {
  chr_data <- readRDS(paste0(Robject_dir, "chromosome_data.Rdata"))
}
# prepare df for plotting:
plot_object <- heatmap_df
old_scores <- c(0, 0.5, 1, 1.5, 2, 3)
new_scores <- c(-2, -1, 0, 1, 2, 3)
for (s in 1:length(new_scores)) {
  plot_object[plot_object == old_scores[s]] <- new_scores[s]
}

colnames(plot_object) <- rep("la", ncol(plot_object))

# define heatmap colours:
na_less_vector <- unlist(plot_object)
na_less_vector <- na_less_vector[!is.na(na_less_vector)]
heatmap_cols <- colorRamp2(c(-2, -1, 0, 1, 2, 3), 
      c("#00106B", "#9191CC", "white", "#DDB6B6", "#AB4848", "#930707"), 
      space = "sRGB")
print("Generating final heatmap...")

# create main CNV heatmap:
final_heatmap <- Heatmap(
  plot_object, name = paste0("hm"), 
  col = heatmap_cols,
  cluster_columns = F, cluster_rows = F,
  show_row_names = F, show_column_names = T,
  column_names_gp = gpar(col = "white"),
  show_row_dend = F,
  heatmap_legend_param = list(title = "Modified\nexpression", color_bar = "continuous", 
    grid_height = unit(1.5, "cm"), grid_width = unit(1.5, "cm"), legend_direction = "horizontal",
    title_gp = gpar(fontsize = 12, fontface = "bold"), labels_gp = gpar(fontsize = 12)),
  use_raster = T, raster_device = c("png")
)

ht_list <- group_annotation + final_heatmap + 
  GIN_annotation + nUMI_annotation + nGene_annotation

annotated_heatmap <- grid.grabExpr(
  draw(ht_list, gap = unit(6, "mm"), heatmap_legend_side = "left")
)

# determine where starting co-ordinates for heatmap are based upon longest cluster name
# (0.00604 units per character):
longest_cluster_name <- max(nchar(unique(as.character(heatmap_metadata$cell_type))))
x_coord <- longest_cluster_name*0.0037

# plot final annotated heatmap:
pdf(paste0(plot_dir, "HMM_plot.pdf"), height = 13, width = 18)   
  grid.newpage()
  pushViewport(viewport(x = 0.005, y = 0.065, width = 0.99, height = 0.78, 
  	just = c("left", "bottom")))
      grid.draw(annotated_heatmap)
      decorate_heatmap_body("hm", {
        for ( e in 1:length(chr_data$end_pos) ) {
        grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), 
        	gp = gpar(lwd = 1, col = "#383838"))
        }
      })
    popViewport()
    pushViewport(viewport(x=x_coord-0.05, y=0.98, width = 0.1, height = 0.1, 
      just = "bottom"))
      grid.text(as.character(sample_name), gp = gpar(fontsize = 14))
    popViewport()
    if (exists("grid_array_heatmap")) {
    	pushViewport(viewport(x = x_coord + 0.845, y = 0.86, 
     	  width = 0.8657, height = 0.13, just = c("right", "bottom")))
      grid.draw(grid_array_heatmap)
      popViewport()
    }
    
dev.off()

print(paste0("Heatmap created, output in ", plot_dir))

# print epithelial_metadata as table:
write.table(heatmap_metadata, paste0(table_dir, "epithelial_HMM_metadata.txt"), 
  sep="\t", quote=F, row.names=F, col.name=T)

#convert pdf to png:
system(paste0("for p in ", plot_dir, "*.pdf; do echo $p; f=$(basename $p); echo $f; ",
"new=$(echo $f | sed 's/.pdf/.png/'); echo $new; ", 
"convert -density 150 ", plot_dir, "$f -quality 90 ", plot_dir, "$new; done"))

