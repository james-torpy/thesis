#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

subproject_name <- "Figure_3.5"

#args = commandArgs(trailingOnly=TRUE)
#include_t_cells <- as.logical(args[1])
#subset_data <- as.logical(args[2])
#subset_samples <- as.logical(args[3])
#na_colour <- args[4]
#include_normals <- as.logical(args[5])

include_t_cells <- T
subset_data <- F
subset_samples <- F
na_colour <- "white"
include_normals <- F
reclustered_group_annotation <- TRUE

#sample_names <- c("CID4386", "CID43862", "CID43863")
sample_names <- c("CID45171", "CID45172")
samples_string <- paste(sample_names, collapse="_")
sample_types <- c("primary", "metastasis")

heatmap_prefix <- paste0(samples_string,
  "_combined_infercnv_", na_colour, "_missing_values_heatmap"
)

if (include_normals) {
  heatmap_prefix <- paste0(heatmap_prefix, "_normals_included")
}


print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample names = ", sample_names))
print(paste0("Subset data? ", as.character(subset_data)))
print(paste0("Subset samples? ", as.character(subset_samples)))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.5dev/"
library(cluster, lib.loc = lib_loc)
library(ComplexHeatmap, lib.loc=lib_loc)
library(circlize, lib.loc = lib_loc)
library(scales, lib.loc = lib_loc)
library(fpc, lib.loc = lib_loc)
library(Seurat)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(dplyr)

project_name <- "thesis"
home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/",
  subproject_name, "/")
ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/functions/")
results_dir <- paste0(project_dir, "results/")

if (include_t_cells) {
  in_path <- paste0(results_dir, "infercnv/t_cells_included/")
  out_dir <- paste0(results_dir, "infercnv/t_cells_included/",
    samples_string, "_combined/")
} else {
  in_path <- paste0(results_dir, "infercnv/t_cells_excluded/")
  out_dir <- paste0(results_dir, "infercnv/t_cells_excluded/",
    samples_string, "_combined/")
}

if (subset_data & subset_samples) {
  Robject_dir <- paste0(out_dir, "Rdata_double_sub/")
  plot_dir <- paste0(out_dir, "plots_double_sub/")
} else if (subset_data) {
  Robject_dir <- paste0(out_dir, "Rdata_sub/")
  plot_dir <- paste0(out_dir, "plots_sub/")
} else if (subset_samples) {
  Robject_dir <- paste0(out_dir, "Rdata_sample_sub/")
  plot_dir <- paste0(out_dir, "plots_sample_sub/")
} else {
  Robject_dir <- paste0(out_dir, "Rdata/")
  plot_dir <- paste0(out_dir, "plots/")
}

system(paste0("mkdir -p ", Robject_dir))
system(paste0("mkdir -p ", plot_dir))

print(paste0("Reference directory = ", ref_dir))
print(paste0("R function directory = ", func_dir))
print(paste0("Output directory = ", out_dir))
print(paste0("R object directory = ", Robject_dir))
print(paste0("Plot directory = ", plot_dir))


################################################################################
### 0. Define functions ###
################################################################################

fetch_chromosome_boundaries <- dget(paste0(func_dir, 
  "fetch_chromosome_boundaries.R"))
create_PAM50_CNV_annotation <- dget(paste0(func_dir, 
  "create_PAM50_CNV_annotation.R"))


################################################################################
### 1. Load heatmap data and combine ###
################################################################################

print("Loading heatmap and metadata dfs for each sample...")

if (!file.exists(paste0(Robject_dir, "initial_combined_heatmap.Rdata")) | 
  !file.exists(paste0(Robject_dir, "initial_combined_heatmap_metadata.Rdata"))) {
  for (s in 1:length(sample_names)) {
    print(sample_names[s])

    sample_Robject_dir <- paste0(in_path, sample_names[[s]], 
      "/Rdata/")
    heatmap_df <- readRDS(
      paste0(sample_Robject_dir, "7a.epithelial_heatmap_final.Rdata")
    )
    heatmap_metadata <- readRDS(
      paste0(sample_Robject_dir, "7b.epithelial_metadata_final.Rdata")
    )
    if (reclustered_group_annotation & 
      "reclustered_cell_type" %in% colnames(heatmap_metadata)) {
      heatmap_metadata <- subset(heatmap_metadata, 
        select = c(cell_ids, reclustered_cell_type, nUMI, nGene, CNA_value, cor.estimate, 
        cor.p.value, sc50))
      colnames(heatmap_metadata) <- gsub(
        "reclustered_cell_type", "cell_type",
        colnames(heatmap_metadata)
      )
    } else {
      heatmap_metadata <- subset(heatmap_metadata, 
        select = c(cell_ids, cell_type, nUMI, nGene, CNA_value, cor.estimate, 
        cor.p.value, sc50))
    }
    

    print(paste0(
      "Are heatmap_metadata rownames in the same order as heatmap_df?? ",
      identical(rownames(heatmap_df), rownames(heatmap_metadata))
    ))
  
    # if > 2 cells present, add heatmap and metadata to list:
    if (nrow(heatmap_df) > 2) {
      if (s==1) {
        heatmap_dfs <- list(heatmap_df)
        names(heatmap_dfs) <- sample_names[s]
        gene_list <- colnames(heatmap_df)
        heatmap_metadatas <- list(heatmap_metadata)
        names(heatmap_metadatas) <- sample_names[s]
      } else if (exists("heatmap_dfs")) {
        heatmap_dfs[[s]] <- heatmap_df
        names(heatmap_dfs)[s] <- sample_names[s]
        gene_list <- c(gene_list, colnames(heatmap_df))
        heatmap_metadatas[[s]] <- heatmap_metadata
        names(heatmap_metadatas)[s] <- sample_names[s]
      } else {
        heatmap_dfs <- list(NULL)
        heatmap_dfs[[s]] <- heatmap_df
        names(heatmap_dfs)[s] <- sample_names[s]
        gene_list <- c(gene_list, colnames(heatmap_df))
        heatmap_metadatas <- list(NULL)
        heatmap_metadatas[[s]] <- heatmap_metadata
        names(heatmap_metadatas)[s] <- sample_names[s]
      }
    } else {
      sample_names <- sample_names[grep(sample_name, sample_names, invert=T)]
    }
  }
  # remove NULL values from lists:
  heatmap_dfs[sapply(heatmap_dfs, is.null)] <- NULL
  heatmap_metadatas[sapply(heatmap_metadatas, is.null)] <- NULL
  
  # get complete gene list from samples and order:
  print("Adding missing genes and collating heatmap and metadata dfs...")
  gene_order <- read.table(paste0(ref_dir, "infercnv_gene_order.txt"))
  gene_list <- unique(gene_list)
  gene_list <- as.character(gene_order$V1[gene_order$V1 %in% gene_list])
  
  # add all missing genes as columns to all heatmap dfs:
  complete_heatmap_dfs <- lapply(heatmap_dfs, function(x) {
    
    missing_genes <- gene_list[!(gene_list %in% colnames(x))]
    missing_genes_df <- data.frame(matrix(NA, nrow = nrow(x), 
      ncol = length(missing_genes)))
    colnames(missing_genes_df) <- missing_genes
  
    complete_df <- cbind(x, missing_genes_df)
    m <- match(gene_list, colnames(complete_df))
    complete_df <- complete_df[,m]
  
    return(complete_df)
  
  })
  
  # collate all dfs and check rows line up:
  group_heatmap <- do.call(rbind, complete_heatmap_dfs)
  rownames(group_heatmap) <- gsub("^.*\\.C", "C", rownames(group_heatmap))
  print("Are all rows present in group heatmap df?")
  identical(nrow(group_heatmap), sum(unlist(lapply(complete_heatmap_dfs, nrow))))
  print("Are all genes present in group heatmap df?")
  identical(colnames(group_heatmap), gene_list)
  
  # subset group_heatmap if needed:
  if (subset_data) {
    group_heatmap <- group_heatmap[,1:300]
  }
  
  # collate heatmap metadata:
  group_heatmap_metadata <- do.call("rbind", heatmap_metadatas)
  rownames(group_heatmap_metadata) <- gsub(
    "^.*\\.", "", rownames(group_heatmap_metadata)
  )

  # select only group_heatmap_metadata rows in group_heatmap and order:
  group_heatmap_metadata <- group_heatmap_metadata[rownames(group_heatmap),]
  m <- match(rownames(group_heatmap), rownames(group_heatmap_metadata))
  group_heatmap_metadata <- group_heatmap_metadata[m,]

  group_heatmap <- group_heatmap[rownames(group_heatmap_metadata),]


  ################################################################################
  ### 2. Add type and sample columns to group_heatmap_metadata ###
  ################################################################################
  
  # add sample column to heatmap_metadata:
  group_heatmap_metadata$sample <- gsub(
    "_.*$", "", group_heatmap_metadata$cell_id
  )
  group_heatmap_metadata$type <- "primary"
  for (t in 1:length(sample_types)) {
    group_heatmap_metadata$type[
      grep(sample_names[t], group_heatmap_metadata$sample)
    ] <- sample_types[t]
  }

  group_heatmap_metadata$cell_type <- factor(
    group_heatmap_metadata$cell_type,
    levels = naturalsort(unique(group_heatmap_metadata$cell_type))
  )
  
  print(paste0(
      "Are group_heatmap_metadata rownames still in the same order as group_heatmap?? ",
      identical(rownames(group_heatmap), rownames(group_heatmap_metadata))
  ))
  
  saveRDS(group_heatmap, paste0(Robject_dir, "initial_group_heatmap.Rdata"))
  saveRDS(
    group_heatmap_metadata, 
    paste0(Robject_dir, "initial_group_heatmap_metadata.Rdata")
  )

} else {

  group_heatmap <- readRDS(paste0(Robject_dir, "initial_group_heatmap.Rdata"))
  group_heatmap_metadata <- readRDS(
    paste0(Robject_dir, "initial_group_heatmap_metadata.Rdata")
  )
}


################################################################################
### 3a. Create heatmap row annotations ###
################################################################################

# create grouping colours:
extra_colours <- c(brewer.pal(12, "Set3"), brewer.pal(12, "Paired"))
col_palette <- c(brewer.pal(8, "Dark2"), brewer.pal(12, "Set3"), brewer.pal(8, "Accent"),
  "#660200", "#918940", "black", "#9ECAE1", "#3F007D", "#67000D", "#FDD0A2", "#08519C",
  "#DBB335", "#5AA050", "#807DBA", "#1A6000", "#F16913", "#FD8D3C", "#DEEBF7", 
  "#7F2704", "#DADAEB", "#FC9272", "#BCBDDC", extra_colours, "#b2182b", "#85929E", 
  "#9B59B6", "#74add1","#1b7837", "#b8e186", "#fed976","#e7298a", "#18ffff", "#ef6c00",
  "#A93226", "black","orange", "#b8bc53", "#5628ce", "#fa909c", "#8ff331","#270e26")
col_palette2 <- c(extra_colours, "#b2182b", "#85929E", 
  "#9B59B6", "#74add1","#1b7837", "#b8e186", "#fed976","#e7298a", "#18ffff", "#ef6c00",
  "#A93226", "black","orange", "#b8bc53", "#5628ce", "#fa909c", "#8ff331","#270e26",
  brewer.pal(8, "Dark2"), brewer.pal(12, "Set3"), brewer.pal(8, "Accent"),
  "#660200", "#918940", "black", "#9ECAE1", "#3F007D", "#67000D", "#FDD0A2", "#08519C",
  "#DBB335", "#5AA050", "#807DBA", "#1A6000", "#F16913", "#FD8D3C", "#DEEBF7", 
  "#7F2704", "#DADAEB", "#FC9272", "#BCBDDC")

# create sample column for group_heatmap_metadata:
group_heatmap_metadata$sample <- gsub(
  "_.*$", "",
  group_heatmap_metadata$cell_ids,
)

# order group_heatmap the same way as group_heatmap_metadata:
group_heatmap <- group_heatmap[rownames(group_heatmap_metadata),]
print(paste0("Are group_heatmap and group_metadata_df rownames identical? ",
  identical(rownames(group_heatmap), rownames(group_heatmap_metadata))))

# create group annotation:
group_annotation_df <- subset(group_heatmap_metadata, select = cell_type)
print(paste0("Are sample_annotation_df and group_heatmap rownames identical? ",
  identical(rownames(group_heatmap), rownames(group_annotation_df))))
# determine colours and create annotation:
group_number <- length(unique(group_annotation_df$cell_type))
group_cols <- col_palette[1:group_number]
group_annotation <- Heatmap(
  group_annotation_df, 
  col = group_cols, 
  name = "group_annotation", 
  width = unit(4, "mm"), 
  show_row_names = F, show_column_names = F,
  heatmap_legend_param = list(
    title = "Cluster", title_gp = gpar(fontsize = 20, fontface = "bold",
    lineheight = 1.5), 
    labels_gp = gpar(fontsize = 20, lineheight = 1.5), 
    at = as.character(levels(group_annotation_df$cell_type))
  )
)

# create sample annotation:
sample_annotation_df <- subset(group_heatmap_metadata, select = sample)
sample_annotation_df$sample <- factor(sample_annotation_df$sample,
  levels=unique(sample_annotation_df$sample))
print(paste0("Are sample_annotation_df and group_heatmap rownames identical? ",
  identical(rownames(group_heatmap), rownames(sample_annotation_df))))
# determine colours and create annotation:
sample_number <- length(unique(sample_annotation_df$sample))
sample_cols <- rev(col_palette)[1:sample_number]
sample_annotation <- Heatmap(
  sample_annotation_df, 
  col = sample_cols, 
  name = "sample_annotation", 
  width = unit(4, "mm"), 
  show_row_names = F, show_column_names = F,
  heatmap_legend_param = list(
    title = "Sample", title_gp = gpar(fontsize = 20, fontface = "bold",
    lineheight = 1.5), 
    labels_gp = gpar(fontsize = 20, lineheight = 1.5), 
    at = as.character(levels(sample_annotation_df$sample))
  )
)

# create type annotation:
type_annotation_df <- subset(group_heatmap_metadata, select = type)
type_annotation_df$type <- factor(type_annotation_df$type,
  levels=unique(type_annotation_df$type))
print(paste0("Are type_annotation_df and group_heatmap rownames identical? ",
  identical(rownames(group_heatmap), rownames(type_annotation_df))))
# determine colours and create annotation:
type_number <- length(unique(type_annotation_df$type))
type_cols <- col_palette[1:type_number]
type_annotation <- Heatmap(
  type_annotation_df, 
  col = type_cols, 
  name = "type_annotation", 
  width = unit(4, "mm"), 
  show_row_names = F, show_column_names = F,
  heatmap_legend_param = list(
      title = "Type", title_gp = gpar(fontsize = 20, fontface = "bold"), 
      labels_gp = gpar(fontsize = 20, lineheight = 1.5), 
      at = as.character(levels(type_annotation_df$type))
    )
)


################################################################################
### 3b. Create heatmap row annotations ###
################################################################################

# create QC annotations:
nUMI_annotation <- rowAnnotation(
  correlation_annotation = anno_barplot(
    group_heatmap_metadata$nUMI, name = "nUMI",
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
    group_heatmap_metadata$nGene, name = "nGene",
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
### 4. Create heatmap column annotations ###
################################################################################

# create CNA annotation:
CNA_value_annotation <- rowAnnotation(
  correlation_annotation = anno_barplot(
    group_heatmap_metadata$CNA_value,
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

# determine co-ordinates of vertical lines at chromosome borders:
chr_data <- fetch_chromosome_boundaries(group_heatmap, ref_dir)
saveRDS(chr_data, paste0(Robject_dir, "chromosome_data.Rdata"))


#############################################################################
### 5. Prepare heatmap and annotations ###
################################################################################

if (!file.exists(paste0(Robject_dir, "final_group_heatmap.Rdata"))) {
  if (!exists("split_heatmap")) {
    sample_rownames <- split(
      rownames(group_heatmap), gsub("_.*$", "", rownames(group_heatmap))
    )
    split_heatmap <- lapply(sample_rownames, function(x) group_heatmap[x,])
  }
  
  # scale all samples the same:
  print(paste0("Summaries of split heatmap df before rescaling to min/max values: "))
  print(lapply(split_heatmap, function(x) summary(unlist(x))))
  
  rescaled_split_heatmap <- lapply(split_heatmap, function(x) {
    return(as.data.frame(rescale(as.matrix(x), c(-1, 1))))
  })
  
  print(paste0("Summaries of split heatmap df after rescaling: "))
  print(lapply(rescaled_split_heatmap, function(x) summary(unlist(x))))
  
  
  rescaled_heatmap_df <- do.call("rbind", rescaled_split_heatmap)
  rownames(rescaled_heatmap_df) <- gsub("^.*\\.", "", rownames(rescaled_heatmap_df))
  group_heatmap <- rescaled_heatmap_df[rownames(group_heatmap_metadata),]
  print(paste0("Are group_heatmap and group_heatmap_metadata rownames identical? ",
    identical(rownames(group_heatmap), rownames(group_heatmap_metadata))))
  
  saveRDS(group_heatmap, paste0(Robject_dir, "final_group_heatmap.Rdata"))

  saveRDS(
    group_heatmap_metadata, 
    paste0(Robject_dir, "final_group_heatmap_metadata.Rdata")
  )

} else {
  group_heatmap <- readRDS(paste0(Robject_dir, "final_group_heatmap.Rdata"))
  group_heatmap_metadata <- readRDS(paste0(Robject_dir, "final_group_heatmap_metadata.Rdata"))
}

# prepare df for plotting:
plot_object <- group_heatmap
colnames(plot_object) <- rep("la", ncol(plot_object))

# define heatmap colours:
na_less_vector <- unlist(plot_object)
na_less_vector <- na_less_vector[!is.na(na_less_vector)]
heatmap_cols <- colorRamp2(c(min(na_less_vector), 0, max(na_less_vector)), 
      c("#00106B", "white", "#680700"), space = "sRGB")

print("Generating final heatmap...")

# create main CNV heatmap:
final_heatmap <- Heatmap(
  plot_object, name = "hm",
  na_col = na_colour,
  col = heatmap_cols,
  cluster_columns = F, cluster_rows = F,
  show_row_names = F, show_column_names = T,
  column_names_gp = gpar(col = "white"),
  show_row_dend = FALSE,
#  bottom_annotation = CNV_genes_annotation, bottom_annotation_height = unit(2, "cm"),
#    gap = unit(1, "cm"),
  heatmap_legend_param = list(title = "Modified\nexpression", color_bar = "continuous", 
  grid_height = unit(2, "cm"), grid_width = unit(1, "cm"), legend_direction = "horizontal",
  title_gp = gpar(fontsize = 20, fontface = "bold"), labels_gp = gpar(fontsize = 18)),
  use_raster = T, raster_device = c("png")
)

# determine co-ordinates of horizontal lines at group borders:
spl_groups <- split(group_heatmap_metadata$sample, 
  group_heatmap_metadata$sample)
spl_groups <- spl_groups[unique(group_heatmap_metadata$sample)]
if (length(spl_groups) > 1) {
  for ( n in 1:(length(spl_groups)-1) ) {
    if (n==1) {
      hlines <- c(length(spl_groups[[n]])/length(group_heatmap_metadata$sample))
    } else {
      hlines[n] <- hlines[n-1] + length(spl_groups[[n]])/length(group_heatmap_metadata$sample)
    }
  }
  hlines <- 1-hlines
} else {
  hlines <- c(length(spl_groups[[1]])/length(group_heatmap_metadata$sample))
}


ht_list <- type_annotation + sample_annotation + group_annotation +
  final_heatmap + CNA_value_annotation + nUMI_annotation + nGene_annotation
#ht_list <- type_annotation + final_heatmap

annotated_heatmap <- grid.grabExpr(
  draw(ht_list, gap = unit(3, "mm"), heatmap_legend_side = "left")
)

# determine where starting co-ordinates for heatmap are based upon longest cluster name
# (0.00604 units per character):
sample_names_and_types <- c(unique(as.character(group_heatmap_metadata$sample)), 
  unique(as.character(group_heatmap_metadata$type))
)
longest_label <- max(nchar(sample_names_and_types))
x_coord <- longest_label*0.0037


#############################################################################
### 5. Plot heatmap and annotations ###
################################################################################

#if (include_normals) {
#  pdf(paste0(plot_dir, heatmap_prefix, "_rescaled_with_normals.pdf"), height = 21, 
#    width = 20)
#} else {
  pdf(paste0(plot_dir, heatmap_prefix, "_rescaled.pdf"), height = 21, width = 20)
#}
grid.newpage()

  pushViewport(viewport(x = 0, y = 0.3, 
                      width = 0.99, height = 0.68, just = c("left", "bottom")))
    grid.draw(annotated_heatmap)
    decorate_heatmap_body("hm", {
      for ( e in 1:length(chr_data$end_pos) ) {
        grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), gp = gpar(lwd = 2, 
          col = "#383838"))
        grid.text(gsub("chr", "", names(chr_data$lab_pos)[e]), chr_data$lab_pos[e], 
          unit(0, "npc") + unit(-2.1, "mm"), gp=gpar(fontsize=20))
      }
      for ( m in 1:length(hlines) ) {
        grid.lines(c(0, 1), c(hlines[m], hlines[m]), gp = gpar(lwd = 2, col = "#383838"))
      }
    })
  popViewport()

  pushViewport(viewport(x=x_coord + 0.91, y=0.3, width = 0.1, height = 0.1, just = "bottom"))
    grid.text("CNA", rot=65, gp=gpar(fontsize=20))
  popViewport()
  pushViewport(viewport(x=x_coord + 0.925, y=0.3, width = 0.1, height = 0.1, just = "bottom"))
    grid.text("nUMI", rot=65, gp=gpar(fontsize=20))
  popViewport()
  pushViewport(viewport(x=x_coord + 0.94, y=0.3, width = 0.1, height = 0.1, just = "bottom"))
    grid.text("nGene", rot=65, gp=gpar(fontsize=20))
  popViewport()
    
dev.off()

print(paste0("Heatmap created, output in ", plot_dir))

#convert pdf to png:
system(paste0("for p in ", plot_dir, "*.pdf; do echo $p; f=$(basename $p); echo $f; ",
"new=$(echo $f | sed 's/.pdf/.png/'); echo $new; ", 
"convert -density 150 ", plot_dir, "$f -quality 90 ", plot_dir, "$new; done"))


