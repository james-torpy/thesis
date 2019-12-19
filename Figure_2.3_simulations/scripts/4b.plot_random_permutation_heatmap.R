#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

RStudio <- FALSE

project_name <- "thesis"
subproject_name <- "Figure_2.3_simulations"
args = commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
include_t_cells <- as.logical(args[2])
random_proportion <- args[3]

#sample_name <- "CID4520N"
#include_t_cells <- TRUE
#random_proportion <- "0.01"

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("T-cells included? ", as.character(include_t_cells)))
print(paste0("Proportion of genes randomly permutated = ", 
  as.character(random_proportion)))

library(Seurat)
library(ggplot2)
library(cowplot)


if (RStudio) {
  
  library(ComplexHeatmap)
  library(circlize)
  library(scales)
  library(fpc)
  library(naturalsort)
  
  home_dir <- "/Users/jamestorpy/clusterHome/"
  
} else {
  
  lib_loc <- "/share/ScratchGeneral/jamtor/R/3.5dev/"
  library(ComplexHeatmap, lib.loc=lib_loc)
  library(circlize, lib.loc = lib_loc)
  library(scales, lib.loc = lib_loc)
  library(fpc, lib.loc = lib_loc)
  library(naturalsort, lib.loc = lib_loc)
  
  home_dir <- "/share/ScratchGeneral/jamtor/"
  
}

project_dir <- paste0(home_dir, "projects/", 
  project_name, "/", subproject_name, "/")
func_dir <- paste0(project_dir, "scripts/functions/")
ref_dir <- paste0(project_dir, "refs/")
results_dir <- seurat_path <- paste0(project_dir, "results/")
seurat_dir <- paste0(project_dir, "raw_files/seurat_objects/", 
  sample_name, "/")

paste0(results_dir, "infercnv/t_cells_included/",
    sample_name, "/", random_proportion, "_genes_randomly_permutated/")

if (include_t_cells) {
  in_dir <- paste0(results_dir, "infercnv/t_cells_included/",
    sample_name, "/", random_proportion, "_genes_randomly_permutated/")
  non_permutated_dir <- paste0(results_dir, "infercnv/t_cells_included/",
    sample_name, "/0_genes_randomly_permutated/")
} else {
  in_dir <- paste0(results_dir, "infercnv/t_cells_excluded/",
    sample_name, "/", random_proportion, "_genes_randomly_permutated/")
  non_permutated_dir <- paste0(results_dir, "infercnv/t_cells_excluded/",
    sample_name, "/0_genes_randomly_permutated/")
}

input_dir <- paste0(in_dir, "input_files/")
Robject_dir <- paste0(in_dir, "Rdata/")
system(paste0("mkdir -p ", Robject_dir))
plot_dir <- paste0(in_dir, "plots/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(in_dir, "tables/")
system(paste0("mkdir -p ", table_dir))

print(paste0("In directory = ", in_dir))
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

  # load metadata df:
  metadata_df <- read.table(paste0(input_dir, "metadata.txt"), header = F,
    sep = "\t", as.is = TRUE)
  colnames(metadata_df) <- c("cell_ids", "cell_type")
  row.names(metadata_df) <- metadata_df$cell_ids

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
### 2. Add annotation metadata ###
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
  count_df <- read.table(paste0(input_dir, "input_matrix.txt"), header = TRUE,
    sep = "\t", as.is = TRUE)
  count_df_rownames <- rownames(count_df)
  nUMI <- apply(count_df, 2, sum)
  nGene <- apply(count_df, 2, function(x) length(x[x!=0]))

  QC <- data.frame(
    row.names = colnames(count_df),
    nUMI = nUMI,
    nGene = nGene
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
### 3. Create QC annotations  ###
################################################################################

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
### 4. Calculate correlations with non-permutated InferCNV output  ###
################################################################################

if (random_proportion != "0") {

  # load InferCNV output:
  print("Loading non-permutated InferCNV heatmap...")
  non_permutated_output <- as.data.frame(t(read.table(paste0(non_permutated_dir, 
    "infercnv.12_denoised.observations.txt"))))

  # load metadata df:
  non_permutated_metadata <- read.table(paste0(non_permutated_dir, 
    "input_files/metadata.txt"), header = F, sep = "\t", as.is = TRUE)
  colnames(non_permutated_metadata) <- c("cell_ids", "cell_type")
  row.names(non_permutated_metadata) <- non_permutated_metadata$cell_ids

  # ensure only epithelial cells in non-permutated data:
  non_permutated_ids <- non_permutated_metadata$cell_ids[
    grep("pithelial", non_permutated_metadata$cell_type)
  ]
  non_permutated_heatmap <- non_permutated_output[
    rownames(non_permutated_output) %in% non_permutated_ids,
  ]

  # calculate mean of permutated and non-permutated data and make them
  # the same length:
  non_permutated_mean_CNV <- apply(non_permutated_heatmap, 2, mean)
  permutated_mean_CNV <- apply(epithelial_heatmap, 2, mean)
  non_permutated_mean_CNV <- non_permutated_mean_CNV[
    names(permutated_mean_CNV)
  ]

  # calculate correlation between non-permutated and permutated CNVs:
  correlation_with_original <- cor.test(
    as.numeric(permutated_mean_CNV), 
    as.numeric(non_permutated_mean_CNV), 
    method = "pearson"
  )
  cor_result <- data.frame(
    R_squared = correlation_with_original$estimate, 
    p_val = correlation_with_original$p.value
  )

  write.table(
    cor_result,
    paste0(table_dir, "correlation_result.txt"),
    sep = "\t",
    row.names = F,
    col.names = F,
    quote = F
  )
}


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

ht_list <- final_heatmap + nUMI_annotation + nGene_annotation

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

    pushViewport(viewport(x=x_coord + 0.917, y=0.025, width = 0.1, height = 0.1, 
      just = "bottom"))
      grid.text("nUMI", rot=65)
    popViewport()
    pushViewport(viewport(x=x_coord + 0.942, y=0.025, width = 0.1, height = 0.1, 
      just = "bottom"))
      grid.text("nGene", rot=65)
    popViewport()

    if (exists("cor_result")) {
      pushViewport(viewport(x=x_coord + 0.072, y=0.71, width = 0.1, height = 0.1, 
        just = "right"))
        grid.text(paste0("Correlation with\noriginal = ", round(cor_result$R_squared, 3)), 
          gp=gpar(fontsize=16))
      popViewport()
      pushViewport(viewport(x=x_coord + 0.055, y=0.66, width = 0.1, height = 0.1, 
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

#convert pdf to png:
system(paste0("for p in ", plot_dir, "*.pdf; do echo $p; f=$(basename $p); echo $f; ",
"new=$(echo $f | sed 's/.pdf/.png/'); echo $new; ", 
"convert -density 150 ", plot_dir, "$f -quality 90 ", plot_dir, "$new; done"))


