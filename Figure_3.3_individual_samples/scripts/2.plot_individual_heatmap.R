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
subcluster_method <- args[2]
subcluster_p <- args[3]
if (subcluster_p != "none") {
  subcluster_p <- as.numeric(subcluster_p)
}
nUMI_threshold <- as.numeric(args[4])
nGene_threshold <- as.numeric(args[5])
min_proportion_cells_for_subcluster <- as.numeric(args[6])
order_by <- args[7]
subcluster_annot <- as.logical(args[8])
subcluster_legend <- as.logical(args[9])
expression_annot <- as.logical(args[10])
expression_legend <- as.logical(args[11])
QC_annot <- as.logical(args[12])
normal_annotation <- as.logical(args[13])
remove_normals <- as.logical(args[14])

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("Subclustered by ", subcluster_method))
print(paste0("Subclustered by pval ", subcluster_p))
print(paste0("nUMI threshold = ", nUMI_threshold))
print(paste0("nGene threshold = ", nGene_threshold))
print(paste0("Min proportion of cells for subcluster to be defined = ", 
  min_proportion_cells_for_subcluster))
print(paste0("Order cells by = ", order_by))
print(paste0("Print subcluster annotation? ", subcluster_annot))
print(paste0("Print subcluster legend? ", subcluster_legend))
print(paste0("Print expression annotation? ", expression_annot))
print(paste0("Print expression legend?", expression_legend))
print(paste0("Print QC annotations? ", QC_annot))

project_name <- "thesis"
subproject_name <- "Figure_3.3_individual_samples"
sample_name <- "prca_filt"
subcluster_method <- "none"
subcluster_p <- "none"
if (subcluster_p != "none") {
  subcluster_p <- as.numeric(subcluster_p)
}
nUMI_threshold <- as.numeric("1100")
nGene_threshold <- as.numeric("300")
coverage_filter <- FALSE
normal_epithelial <- "none"
min_proportion_cells_for_subcluster <- as.numeric("0.005")
order_by <- "expression"
subcluster_annot <- FALSE
subcluster_legend <- FALSE
expression_annot <- TRUE
expression_legend <- TRUE
QC_annot <- TRUE
normal_annot <- FALSE
normal_legend <- FALSE
plot_references <- TRUE
seurat_filename <- "Final_FILT.Rds"
cancer_x_threshold_sd_multiplier <- 2
cancer_y_threshold_sd_multiplier <- 1.5
normal_x_threshold_sd_multiplier <- 1
normal_y_threshold_sd_multiplier <- 1.25
remove_normals <- FALSE

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(scales, lib.loc = lib_loc)
library(cluster)
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
in_dir <- paste0(results_dir, "infercnv/", sample_name, "/", 
  subcluster_method, "/p_", subcluster_p, "/")
input_dir <- paste0(in_dir, "/input_files/")

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
get_subpops <- dget(paste0(func_dir, "get_subpops.R"))
define_normals <- dget(paste0(func_dir, "define_normals.R"))

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

  ######
  # subset:
  # epithelial_heatmap <- epithelial_heatmap[1:200,1:500]
  ######

  # load and isolate epithelial metadata:
  initial_metadata <- readRDS(paste0(Robject_dir, "initial_metadata.Rdata"))
  epithelial_metadata <- initial_metadata[
    grep("Epithelial", initial_metadata$cell_type),
  ]
  # replace dashes with full stops for cell ids:
  epithelial_metadata$cell_ids <- gsub("-", ".", epithelial_metadata$cell_ids)

  # remove cells not present in heatmap:
  print(paste0(
    "No cells in metadata before filtering for those in heatmap only = ", 
    nrow(epithelial_metadata)
  ))
  epithelial_metadata <- epithelial_metadata %>%
    filter(cell_ids %in% rownames(epithelial_heatmap))
  print(paste0(
    "No cells in metadata after filtering for those in heatmap only = ", 
    nrow(epithelial_metadata)
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
 
}


################################################################################
### 2. Add subcluster and normal annotations ###
################################################################################

# remove cells below coverage thresholds:
if (coverage_filter) {
  epithelial_metadata <- epithelial_metadata[
    epithelial_metadata$nUMI >= nUMI_threshold & 
      epithelial_metadata$nGene >= nGene_threshold 
  ]
  epithelial_heatmap <- epithelial_heatmap[epithelial_metadata$cell_ids,]
}

# load final output object and fetch subcluster data:  
if (subcluster_annot) {
  final_infercnv_obj <- readRDS(paste0(in_dir, "run.final.infercnv_obj"))
  epithelial_metadata <- get_subpops(final_infercnv_obj, epithelial_metadata)
}

if (normal_annot | remove_normals) {
  if (!file.exists(
    paste0(Robject_dir, "/2.epithelial_metadata_with_normal_annot.Rdata")
  )) {

    epithelial_metadata <- define_normals(
      epithelial_heatmap, 
      epithelial_metadata,
      plot_dir,
      Robject_dir
    )

    saveRDS(
      epithelial_metadata,
      paste0(Robject_dir, "/2.epithelial_metadata_with_normal_annot.Rdata")
    )

  } else {

    epithelial_metadata <- readRDS(
      paste0(Robject_dir, "/2.epithelial_metadata_with_normal_annot.Rdata")
    )

  }
}

if (remove_normals) {
  epithelial_metadata <- epithelial_metadata[
    epithelial_metadata$normal_cell_call != "normal|unassigned"
  ]
  epithelial_heatmap <- epithelial_heatmap[
    rownames(epithelial_heatmap %in% epithelial_metadata$cell_ids)
  ]

  saveRDS(
      epithelial_metadata,
      paste0(Robject_dir, "/3.epithelial_metadata_without_normals.Rdata")
    )
}

 
################################################################################
### 3. Order metadata ###
################################################################################

if (order_by == "CNV") {

  # if normals annotated, change the order of clusters to plot most normal first:
  if (normal_annot) {

    # split metadata by cluster for ordering:
    metadata_split <- split(
      epithelial_metadata, epithelial_metadata$subcluster_id
    )

    # determine which clusters have most normals:
    normal_nos <- lapply(metadata_split, function(x) {
      res <- data.frame(
        cluster = x$subcluster_id[1],
        normal_no = length(which(x$normal_cell_call == "normal"))
      )
      return(res)
    })
    normal_no_df <- do.call("rbind", normal_nos)

    # order according to normal number:
    normal_no_df <- normal_no_df[
      naturalorder(normal_no_df$normal_no, decreasing = TRUE),
    ]
    metadata_split <- metadata_split[normal_no_df$cluster]

    # rbind metadata together in that order:
    epithelial_metadata <- do.call("rbind", metadata_split)

  } else {

    epithelial_metadata <- epithelial_metadata[
      naturalorder(epithelial_metadata$subcluster_id),
    ]

  }

} else if (order_by == "expression") {

  # if normals annotated, change the order of clusters to plot most normal first:
  if (normal_annot) {

    # split metadata by cluster for ordering:
    metadata_split <- split(
      epithelial_metadata, epithelial_metadata$cell_type
    )

    # determine which clusters have most normals:
    normal_nos <- lapply(metadata_split, function(x) {
      res <- data.frame(
        cluster = x$cell_type[1],
        normal_no = length(which(x$normal_cell_call == "normal"))
      )
      return(res)
    })
    normal_no_df <- do.call("rbind", normal_nos)

    # order according to normal number:
    normal_no_df <- normal_no_df[
      naturalorder(normal_no_df$normal_no, decreasing = TRUE),
    ]
    metadata_split <- metadata_split[normal_no_df$cluster]

    # rbind metadata together in that order:
    epithelial_metadata <- do.call("rbind", metadata_split)

  } else {

    epithelial_metadata <- epithelial_metadata[
      naturalorder(epithelial_metadata$cell_type),
    ]

  }

} else if (order_by == "normals_then_CNVs") {

  epithelial_metadata <- epithelial_metadata[
    naturalorder(epithelial_metadata$subcluster_id),
  ]

  # split metadata by normal status for ordering:
  metadata_split <- split(
    epithelial_metadata, epithelial_metadata$normal_cell_call
  )

  epithelial_metadata <- do.call(
    "rbind",
    list(
      metadata_split$normal,
      metadata_split$unassigned,
      metadata_split$cancer
    )
  )

} else if (order_by == "normals_then_expression") {

   epithelial_metadata <- epithelial_metadata[
    naturalorder(epithelial_metadata$cell_type),
  ]

  # split metadata by normal status for ordering:
  metadata_split <- split(
    epithelial_metadata, epithelial_metadata$normal_cell_call
  )

  epithelial_metadata <- do.call(
    "rbind",
    list(
      metadata_split$normal,
      metadata_split$unassigned,
      metadata_split$cancer
    )
  )

}

# adjust order of heatmap:
epithelial_heatmap <- epithelial_heatmap[epithelial_metadata$cell_ids,]

print(paste0(
  "Is heatmap order the same as metadata? ", 
  identical(rownames(epithelial_heatmap), epithelial_metadata$cell_ids)
))


################################################################################
### 4. Create heatmap and annotations ###
################################################################################

# create expression cluster annotation:
if (expression_annot) {
  # define cluster annotation colours:
  expr_number <- length(unique(epithelial_metadata$cell_type))
  expr_cols <- col_palette[1:expr_number]
  names(expr_cols) <- unique(epithelial_metadata$cell_type)

  expr_annot_df <- subset(epithelial_metadata, select = cell_type)
  expr_annot <- Heatmap(
    as.matrix(expr_annot_df), 
    col = expr_cols,
    name = "Expression\nclusters", 
    width = unit(4, "mm"), 
    show_row_names = F, show_column_names = F,
    show_heatmap_legend = FALSE
  )
}

# create CNV subcluster annotation:
if (subcluster_annot) {
  # create CNV subcluster annotation:
  subcluster_annot_df <- subset(epithelial_metadata, select = subcluster_id)
  subcluster_cols <- rev(col_palette)[
    1:length(unique(subcluster_annot_df$subcluster_id))
  ]
  names(subcluster_cols) <- unique(subcluster_annot_df$subcluster_id)
  CNV_subcluster_annot <- Heatmap(as.matrix(subcluster_annot_df), 
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

# create expression cluster annotation:
if (normal_annot) {
  # define cluster annotation colours:
  normal_number <- length(unique(epithelial_metadata$cell_type))
  normal_cols <- c(
    "normal" = "#1B7837", 
    "unassigned" = "#E7E4D3", 
    "cancer" = "#E7298A"
  )
 
  normal_annot_df <- subset(epithelial_metadata, select = cell_type)
  normal_call_annot <- Heatmap(
    as.matrix(normal_annot_df), 
    col = normal_cols,
    name = "Expression\nclusters", 
    width = unit(4, "mm"), 
    show_row_names = F, show_column_names = F,
    show_heatmap_legend = normal_legend
  )
}

# create QC annotations:
if (QC_annot) {
  nUMI_annot <- rowAnnotation(
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
  nUMI_annot@name <- "nUMI"
  nGene_annot <- rowAnnotation(
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
  nGene_annot@name <- "nGene"
}


################################################################################
### 5. Create epithelial and reference heatmap ###
################################################################################

# determine heatmap colours
if (plot_references) {

  # load InferCNV output:
  print("Loading InferCNV reference output files...")
  ref_heatmap <- as.data.frame(t(read.table(paste0(in_dir, 
  	"infercnv.12_denoised.references.txt"))))

  both_heatmap <- rbind(epithelial_heatmap, ref_heatmap)

  # prepare df for plotting:
  col_heatmap <- both_heatmap
  # define heatmap colours:
  na_less_vector <- unlist(col_heatmap)
  na_less_vector <- na_less_vector[!is.na(na_less_vector)]
  heatmap_cols <- colorRamp2(c(min(na_less_vector), 1, max(na_less_vector)), 
        c("#00106B", "white", "#680700"), space = "sRGB")

  # prepare df for plotting:
  ref_object <- ref_heatmap
  colnames(ref_object) <- rep("la", ncol(ref_object))
  
  print("Generating reference heatmap...")
  # create main CNV heatmap:
  final_ref_heatmap <- Heatmap(
    as.matrix(ref_object), name = paste0("ref_hm"), 
    col = heatmap_cols,
    cluster_columns = F, cluster_rows = F,
    show_row_names = F, show_column_names = T,
    column_names_gp = gpar(col = "white"),
    show_row_dend = F,
    show_heatmap_legend = F,
    use_raster = T, raster_device = c("png")
  )

  # converrt to grid object:
  grid_ref_heatmap <- grid.grabExpr(
    draw(final_ref_heatmap, gap = unit(6, "mm"))
  )
  dev.off()

} else {
 
  # define heatmap colours:
  na_less_vector <- unlist(col_heatmap)
  na_less_vector <- na_less_vector[!is.na(na_less_vector)]
  heatmap_cols <- colorRamp2(c(min(na_less_vector), 1, max(na_less_vector)), 
        c("#00106B", "white", "#680700"), space = "sRGB")

}

# prepare df for plotting:
plot_object <- epithelial_heatmap
colnames(plot_object) <- rep("la", ncol(plot_object))

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

chr_data <- fetch_chromosome_boundaries(epithelial_heatmap, ref_dir)
 
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

# determine filename:
filename <- paste0(plot_dir, "infercnv_plot")
if (any(c(expression_annot, subcluster_annot, QC_annot, normal_annot))) {
  if (expression_annot) {
    filename <- paste0(filename, "_exp_clusters")
  }
  if (subcluster_annot) {
    filename <- paste0(filename, "_CNV_subclusters")
  }
  if (normal_annot) {
    filename <- paste0(filename, "_normals")
  }
  if (QC_annot) {
    filename <- paste0(filename, "_QC")
  }
   filename <- paste0(filename, "_annotated")
}
if (plot_references) {
  filename <- paste0(filename, "_references_plotted")
}

# organise chromosome labels:
chr_labels <- chr_data$lab_pos
names(chr_labels) <- gsub("chr", "", names(chr_labels))
names(chr_labels) <- gsub("^1$", "chr1", names(chr_labels))
names(chr_labels) <- gsub("21", "\n21", names(chr_labels))


################################################################################
### 4. Create and plot annotated heatmap ###
################################################################################

if (subcluster_annot) {

  ht_list <- CNV_subcluster_annot
  if (expression_annot) {
  	ht_list <- ht_list + expr_annot + final_heatmap
  } else {
  	ht_list <- ht_list + final_heatmap
  }
  if (normal_annot) {
    	ht_list <- ht_list + normal_call_annot
    }
  if (QC_annot) {
  	ht_list <- ht_list + nUMI_annot + nGene_annot
  }

} else {

  if (expression_annot) {

    ht_list <- expr_annot + final_heatmap
    if (normal_annot) {
    	ht_list <- ht_list + normal_call_annot
    }
    if (QC_annot) {
    	ht_list <- ht_list + nUMI_annot + nGene_annot
    }

  } else {
  
    ht_list <- final_heatmap
    if (normal_annot) {
    	ht_list <- ht_list + normal_call_annot
    }
    if (QC_annot) {
    	ht_list <- ht_list + nUMI_annot + nGene_annot
    }
  
  }

}

annotated_heatmap <- grid.grabExpr(
  draw(ht_list, gap = unit(6, "mm"), heatmap_legend_side = "left")
)
dev.off()


#save.image(paste0(Robject_dir, "temp.Rdata"))
#load(paste0(Robject_dir, "temp.Rdata"))


if (plot_references) {

  # plot final annotated heatmap:
  png(paste0(filename, "_test.png"), 
    height = 13, width = 20, res = 300, units = "in") 
  
    grid.newpage()

      # create reference heatmap viewport:
      pushViewport(viewport(x = 0.555, y = 0.68, width = 0.76, height = 0.3,
      	just = "bottom"))
      	grid.rect()
        # plot heatmap:
#        grid.draw(grid_ref_heatmap)
#        decorate_heatmap_body("ref_hm", {
#          for ( e in 1:length(chr_data$end_pos) ) {
#          grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), 
#            gp = gpar(lwd = 3, col = "#383838"))
#        })
      popViewport()

      # create reference heatmap label:
      pushViewport(viewport(x = 0.075, y = 0.95, 
      	width = unit(2, "cm"), height = unit(0.5, "cm")))
        #grid.rect()
        grid.text("Reference cells:", gp=gpar(fontsize=20))
      popViewport()

      # create epithelial heatmap viewport:
      pushViewport(viewport(x = 0.5, y = 0, width = 1, height = 0.75,
      	just = "bottom"))

         # plot heatmap:
        pushViewport(viewport(x = 0.155, y = 0.065, width = 0.85, height = 0.85, 
          just = c("left", "bottom")))
          #grid.rect()
          grid.draw(annotated_heatmap)
          decorate_heatmap_body("hm", {
            for ( e in 1:length(chr_data$end_pos) ) {
            grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), 
              gp = gpar(lwd = 1, col = "#383838"))
            grid.text(names(chr_labels)[e], chr_labels[e], 
                unit(0, "npc") + unit(-3.5, "mm"), gp=gpar(fontsize=18))
            }
          })
        popViewport()
   
        # plot heatmap legend:
        pushViewport(viewport(x = unit(2, "cm"), y = unit(19, "cm"), width = unit(0.1, "cm"), 
          height = unit(0.4, "cm"), just = c("right", "bottom")))
          draw(lgd, x = unit(0.1, "cm"), y = unit(0.1, "cm"), just = c("left", "bottom"))
        popViewport()
    
        # plot expression legend:
        if (expression_legend) {
    
          l_labels <- gsub("_", " ", unique(epithelial_metadata$cell_type))
    
          pushViewport(viewport(x = unit(4, "cm"), y = unit(10, "cm"), width = unit(6, "cm"), 
            height = unit(10, "cm")))
            #grid.rect()
    
            # add title:
            pushViewport(viewport(x = 0.06, y = 1, width = unit(2, "cm"), height = unit(0.5, "cm")))
              grid.text("Expression\nclusters", gp=gpar(fontsize=20), just = "left")
              #grid.rect()
            popViewport()
    
            for (l in 1:length(l_labels)) {
              # add labels:
              pushViewport(viewport(
              	x = 0.25, 
              	y = 0.95-(0.12*l), 
              	width = unit(2, "cm"), 
              	height = unit(0.5, "cm")
              ))
                #grid.rect()
                grid.text(l_labels[l], gp=gpar(fontsize=18), just = "left")
              popViewport()
              # add boxes:
              pushViewport(viewport(
              	x = 0.14, 
              	y = 0.95-(0.12*l), 
              	width = unit(0.7, "cm"), 
              	height = unit(0.7, "cm")
              ))
                grid.rect(gp=gpar(col = expr_cols[l], fill = expr_cols[l]))
              popViewport()
            }
          popViewport()
        }
     
        # label annotations:
        pushViewport(viewport(x=0.95, y=0.1, width = 0.1, height = 0.1, 
          just = "top"))
          grid.text("nUMI", rot=65, gp=gpar(fontsize=18))
        popViewport()
        pushViewport(viewport(x=0.98, y=0.1, width = 0.1, height = 0.1, 
          just = "top"))
          grid.text("nGene", rot=65, gp=gpar(fontsize=18))
        popViewport()

      popViewport()
        
  dev.off()

} else {

  # plot final annotated heatmap:
  png(paste0(filename, ".png"), 
    height = 13, width = 20, res = 300, units = "in") 
  
    grid.newpage()
  
      # plot heatmap:
      pushViewport(viewport(x = 0.155, y = 0.065, width = 0.8, height = 0.85, 
        just = c("left", "bottom")))
        #grid.rect()
        grid.draw(annotated_heatmap)
        decorate_heatmap_body("hm", {
          for ( e in 1:length(chr_data$end_pos) ) {
          grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), 
            gp = gpar(lwd = 1, col = "#383838"))
          grid.text(names(chr_labels)[e], chr_labels[e], 
              unit(0, "npc") + unit(-3.5, "mm"), gp=gpar(fontsize=18))
          }
        })
      popViewport()
  
      # plot heatmap legend:
      pushViewport(viewport(x = unit(2, "cm"), y = unit(26, "cm"), width = unit(0.1, "cm"), 
        height = unit(0.4, "cm"), just = c("right", "bottom")))
        draw(lgd, x = unit(0.1, "cm"), y = unit(0.1, "cm"), just = c("left", "bottom"))
      popViewport()
  
      # plot expression legend:
      if (expression_legend) {
  
        l_labels <- unique(epithelial_metadata$cell_type)
  
        pushViewport(viewport(x = unit(4, "cm"), y = unit(18.5, "cm"), width = unit(6, "cm"), 
          height = unit(10, "cm")))
          #grid.rect()
  
          # add title:
          pushViewport(viewport(x = 0.06, y = 1, width = unit(2, "cm"), height = unit(0.5, "cm")))
            grid.text("Expression\nclusters", gp=gpar(fontsize=20), just = "left")
            #grid.rect()
          popViewport()
  
          for (l in 1:length(l_labels)) {
            # add labels:
            pushViewport(viewport(
            	x = 0.25, 
            	y = 0.95-(0.12*l), 
            	width = unit(2, "cm"), 
            	height = unit(0.5, "cm")
            ))
              #grid.rect()
              grid.text(l_labels[l], gp=gpar(fontsize=18), just = "left")
            popViewport()
            # add boxes:
            pushViewport(viewport(
            	x = 0.14, 
            	y = 0.95-(0.12*l), 
            	width = unit(0.7, "cm"), 
            	height = unit(0.7, "cm")
            ))
              grid.rect(gp=gpar(col = expr_cols[l], fill = expr_cols[l]))
            popViewport()
          }
        popViewport()
      }
   
      # label annotations:
      pushViewport(viewport(x=0.9, y=0.1, width = 0.1, height = 0.1, 
        just = "top"))
        grid.text("nUMI", rot=65, gp=gpar(fontsize=18))
      popViewport()
      pushViewport(viewport(x=0.93, y=0.1, width = 0.1, height = 0.1, 
        just = "top"))
        grid.text("nGene", rot=65, gp=gpar(fontsize=18))
      popViewport()
      
  dev.off()
  
  pdf(paste0(filename, ".pdf"), 
    height = 13, width = 20) 
  
    grid.newpage()
  
      # plot heatmap:
      pushViewport(viewport(x = 0.155, y = 0.065, width = 0.8, height = 0.85, 
        just = c("left", "bottom")))
        #grid.rect()
        grid.draw(annotated_heatmap)
        decorate_heatmap_body("hm", {
          for ( e in 1:length(chr_data$end_pos) ) {
          grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), 
            gp = gpar(lwd = 1, col = "#383838"))
          grid.text(names(chr_labels)[e], chr_labels[e], 
              unit(0, "npc") + unit(-3.5, "mm"), gp=gpar(fontsize=18))
          }
        })
      popViewport()
  
      # plot heatmap legend:
      pushViewport(viewport(x = unit(2, "cm"), y = unit(26, "cm"), width = unit(0.1, "cm"), 
        height = unit(0.4, "cm"), just = c("right", "bottom")))
        draw(lgd, x = unit(0.1, "cm"), y = unit(0.1, "cm"), just = c("left", "bottom"))
      popViewport()
  
      # plot expression legend:
      if (expression_legend) {
  
        l_labels <- unique(epithelial_metadata$cell_type)
  
        pushViewport(viewport(x = unit(4, "cm"), y = unit(18.5, "cm"), width = unit(6, "cm"), 
          height = unit(10, "cm")))
          #grid.rect()
  
          # add title:
          pushViewport(viewport(x = 0.06, y = 1, width = unit(2, "cm"), height = unit(0.5, "cm")))
            grid.text("Expression\nclusters", gp=gpar(fontsize=20), just = "left")
            #grid.rect()
          popViewport()
  
          for (l in 1:length(l_labels)) {
            # add labels:
            pushViewport(viewport(
              x = 0.25, 
              y = 0.95-(0.12*l), 
              width = unit(2, "cm"), 
              height = unit(0.5, "cm")
            ))
              #grid.rect()
              grid.text(l_labels[l], gp=gpar(fontsize=18), just = "left")
            popViewport()
            # add boxes:
            pushViewport(viewport(
              x = 0.14, 
              y = 0.95-(0.12*l), 
              width = unit(0.7, "cm"), 
              height = unit(0.7, "cm")
            ))
              grid.rect(gp=gpar(col = expr_cols[l], fill = expr_cols[l]))
            popViewport()
          }
        popViewport()
      }
   
      # label annotations:
      pushViewport(viewport(x=0.9, y=0.1, width = 0.1, height = 0.1, 
        just = "top"))
        grid.text("nUMI", rot=65, gp=gpar(fontsize=18))
      popViewport()
      pushViewport(viewport(x=0.93, y=0.1, width = 0.1, height = 0.1, 
        just = "top"))
        grid.text("nGene", rot=65, gp=gpar(fontsize=18))
      popViewport()
      
  dev.off()

}

print(paste0("Heatmap created, output in ", plot_dir))
  
  
  