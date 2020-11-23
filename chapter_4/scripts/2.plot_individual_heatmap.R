#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

project_name <- "thesis"
subproject_name <- "chapter_4"
args = commandArgs(trailingOnly=TRUE)

sample_name <- args[1]
subcluster_method <- args[2]
subcluster_p <- args[3]
if (subcluster_p != "none") {
  subcluster_p <- as.numeric(subcluster_p)
}
coverage_filter <- args[4]
remove_artefacts <- args[5]
epi_res <- args[6]
epi_PC <- args[7]
x_outlier_multiplier <- as.numeric(args[8])
x_thresh_multiplier <- as.numeric(args[9])
y_outlier_multiplier <- as.numeric(args[10])
y_thresh_multiplier <- as.numeric(args[11])
min_cluster_cells <- as.numeric(args[12])
subcluster_merge <- as.logical(args[13])
merge_thresh <- as.numeric(args[14])
merge_diff_prop <- as.numeric(args[15])
order_by <- args[16]
QC_annot <- as.logical(args[17])
plot_references <- as.logical(args[18])
array_CNVs <- as.logical(args[19])
plot_type <- args[20]

project_name <- "thesis"
subproject_name <- "chapter_4"
sample_name <- "CID3948"
subcluster_method <- "random_trees"
subcluster_p <- "0.05"
if (subcluster_p != "none") {
  subcluster_p <- as.numeric(subcluster_p)
}
coverage_filter <- "filtered"
remove_artefacts <- "artefacts_not_removed"
epi_res <- "PC_C_res.1"
epi_PC <- "C"
x_outlier_multiplier <- 1.5
x_thresh_multiplier <- 3
y_outlier_multiplier <- 1.5
y_thresh_multiplier <- 3
min_cluster_cells <- 5
subcluster_merge <- TRUE
merge_thresh <- 0.95
merge_diff_prop <- 0.75
order_by <- "CNV"
QC_annot <- TRUE
plot_references <- FALSE
array_CNVs <- FALSE
plot_type <- "normals_annotated"

if (plot_type == "normals_annotated") {
  remove_normals <- FALSE
  subcluster_annot <- FALSE
  subcluster_legend <- FALSE
  expression_annot <- FALSE
  expression_legend <- FALSE
  normal_annot <- TRUE
  normal_legend <- TRUE
} else if (plot_type == "clusters_annotated") {
  remove_normals <- TRUE
  subcluster_annot <- TRUE
  subcluster_legend <- TRUE
  expression_annot <- TRUE
  expression_legend <- TRUE
  normal_annot <- FALSE
  normal_legend <- FALSE
}


print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("Subclustered by ", subcluster_method))
print(paste0("Subclustered by pval ", subcluster_p))
print(paste0("Filter by coverage? ", coverage_filter))
print(paste0("Remove artefacts? ", remove_artefacts))
print(paste0("X-axis outlier multiplier = ", x_outlier_multiplier))
print(paste0("X-axis threshold multiplier = ", x_thresh_multiplier))
print(paste0("Y-axis outlier threshold multiplier = ", y_outlier_multiplier))
print(paste0("Y-axis threshold multiplier = ", y_thresh_multiplier))
print(paste0("Remove normal cells? ", remove_normals))
print(paste0("Minimum cells required to call CNV subcluster = ", min_cluster_cells))
print(paste0("Merge identical subclusters? ", subcluster_merge))
print(paste0("Minimum correlation required to merge subclusters ", merge_thresh))
print(paste0("Proportion of less coverage subcluster to more coverage required to merge subclusters ", merge_diff_prop))
print(paste0("Order cells by = ", order_by))
print(paste0("Print subcluster annotation? ", subcluster_annot))
print(paste0("Print subcluster legend? ", subcluster_legend))
print(paste0("Print expression annotation? ", expression_annot))
print(paste0("Print expression legend? ", expression_legend))
print(paste0("Print QC annotations? ", QC_annot))
print(paste0("Print normal annotation? ", normal_annot))
print(paste0("Print normal legend? ", normal_legend))
print(paste0("Plot reference cells? ", plot_references))
print(paste0("Array CNVs? ", array_CNVs))

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
library(tibble)
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
in_dir <- paste0(results_dir, "infercnv/", sample_name, "/", coverage_filter, "/", 
  subcluster_method, "/p_", subcluster_p, "/")
input_dir <- paste0(in_dir, "/input_files/")
common_Robject_dir <- Robject_dir <- paste0(in_dir, "Rdata/")

out_dir <- paste0(in_dir, remove_artefacts, "/")

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


################################################################################
### 0. Define functions ###
################################################################################

fetch_chromosome_boundaries <- dget(paste0(func_dir, 
  "fetch_chromosome_boundaries.R"))
get_subpops <- dget(paste0(func_dir, "get_subpops.R"))
define_normals <- dget(paste0(func_dir, "define_normals.R"))
merge_subclusters <- dget(paste0(func_dir, "merge_subclusters.R"))
create_array_CNV_annotation <- dget(paste0(func_dir, 
  "create_array_CNV_annotation.R"))
create_legend <- dget(paste0(func_dir, "create_legend.R"))

extra_colours <- c(brewer.pal(12, "Set3"), brewer.pal(12, "Paired"))
col_palette <- c(brewer.pal(8, "Dark2")[c(3:6,8)], brewer.pal(12, "Set3"), brewer.pal(8, "Accent"),
  "#660200", "#918940", "#9ECAE1", "#3F007D", "#67000D", "#FDD0A2", "#08519C",
  "#DBB335", "#5AA050", "#807DBA", "#1A6000", "#F16913", "#FD8D3C", "#DEEBF7", 
  "#7F2704", "#DADAEB", "#FC9272", "#BCBDDC", extra_colours, "#b2182b", "#85929E", 
  "#9B59B6", "#74add1","#1b7837", "#b8e186", "#fed976","#e7298a", "#18ffff", "#ef6c00",
  "#A93226", "#E611ED","orange", "#b8bc53", "#5628ce", "#fa909c", "#8ff331")
col_palette <- col_palette[-7]

expr_cols <- read.table(
  paste0(ref_dir, "expression_colour_palette.txt"),
  header = F,
  stringsAsFactors = F,
  comment.char = ""
)[,1]

subcluster_cols <- read.table(
  paste0(ref_dir, "CNV_colour_palette.txt"),
  header = F,
  stringsAsFactors = F,
  comment.char = ""
)[,1]


################################################################################
### 1. Load InferCNV output and create heatmap and metadata dfs ###
################################################################################

if (!file.exists(paste0(Robject_dir, "/1b.initial_epithelial_metadata.Rdata"))) {

  # load InferCNV output:
  print("Loading InferCNV output files...")
  epithelial_heatmap <- as.data.frame(t(read.table(paste0(in_dir, 
  	"infercnv.12_denoised.observations.txt"))))

  # load and isolate epithelial metadata:
  initial_metadata <- readRDS(paste0(common_Robject_dir, "initial_metadata.Rdata"))
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

  if (file.exists(
    paste0(seurat_dir, "05_seurat_object_epithelial_reclustered.Rdata")
  )) {
    # load reclustered epithelial seurat object:
    seurat_epi <- readRDS(
      paste0(seurat_dir, "05_seurat_object_epithelial_reclustered.Rdata")
    )
    # choose resolution if needed:
    if (epi_res != "none") {
      Idents(seurat_epi) <- eval(
        parse(
          text = paste0("seurat_epi@meta.data$", epi_res)
        )
      )
    }
    # add expression cluster column to metadata:
    epithelial_metadata$expression_id <- paste0(
      "Expression_",
      Idents(seurat_epi)[
        match(epithelial_metadata$cell_ids, names(Idents(seurat_epi)))
      ]
    )
  } else {

    epithelial_metadata$expression_id <- gsub(
      "Expression",
      "Epithlelial",
      epithelial_metadata$cell_type
    )

  }
  
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
### 2. Remove artefacts ###
################################################################################

if (remove_artefacts == "artefacts_removed") {

  if (!file.exists(paste0(Robject_dir, 
    "/2.epithelial_heatmap_artefacts_removed.Rdata"))) {

    # load known artefact genes:
    artefact_genes <- readRDS(paste0(ref_dir, "artefact_gene_record.Rdata"))
    
    # determine diploid signal
    signal_table <- table(round(unlist(epithelial_heatmap), 6))
    diploid_signal <- as.numeric(names(signal_table)[which.max(signal_table)])
    
    # locate artefact within heatmap and make all signal within diploid:
    for (a in 1:length(artefact_genes)) {
    	print(a)
      # identify first and last gene in artefact present in epithelial heatmap:
      first_genes <- as.character(
        artefact_genes[[a]]$first_ten[
          artefact_genes[[a]]$first_ten %in% colnames(epithelial_heatmap)
        ]
      )
      first_gene <- first_genes[1]
    
      last_genes <- as.character(
        artefact_genes[[a]]$last_ten[
          artefact_genes[[a]]$last_ten %in% colnames(epithelial_heatmap)
        ]
      )
      last_gene <- last_genes[length(last_genes)]
    
      heatmap_coords <- c(
        grep(paste0("^", first_gene, "$"), colnames(epithelial_heatmap)),
        grep(paste0("^", last_gene, "$"), colnames(epithelial_heatmap))
      )
  
      # add gene tails equal to 1/3 of the artefact length to start and end of 
      # artefact coordinates to account for residue signal:
      tail_length <- round((heatmap_coords[2]-heatmap_coords[1])/3, 0)
      heatmap_coords[1] <- heatmap_coords[1] - tail_length
      heatmap_coords[2] <- heatmap_coords[2] + tail_length
      coord_vector <- heatmap_coords[1]:heatmap_coords[2]
  
      epithelial_heatmap[,coord_vector] <- diploid_signal
    
    }

    saveRDS(epithelial_heatmap, paste0(Robject_dir, 
      "/2.epithelial_heatmap_artefacts_removed.Rdata"))
    
  } else {

  	epithelial_heatmap <- readRDS(paste0(Robject_dir, 
      "/2.epithelial_heatmap_artefacts_removed.Rdata"))

  }
}


################################################################################
### 3. Add subcluster and normal annotations ###
################################################################################

# load final output object and fetch subcluster data:  
if (subcluster_method == "random_trees") {
  final_infercnv_obj <- readRDS(paste0(in_dir, "run.final.infercnv_obj"))
  epithelial_metadata <- get_subpops(final_infercnv_obj, epithelial_metadata)
}

# define cells as normal, unassigned or cancer:
if (!file.exists(
  paste0(Robject_dir, "/3.epithelial_metadata_with_normals.Rdata")
)) {

  epithelial_metadata <- define_normals(
    epithelial_heatmap, 
    epithelial_metadata,
    x_outlier_multiplier,
    x_thresh_multiplier,
    y_outlier_multiplier,
    y_thresh_multiplier,
    plot_dir,
    Robject_dir
  )

  # if needed, label normals and unassigned in subcluster annotation:
  if (subcluster_method == "random_trees") {
  	epithelial_metadata$subcluster_id <- as.character(
  	  epithelial_metadata$subcluster_id
  	)
  	epithelial_metadata$subcluster_id[
  	  epithelial_metadata$normal_cell_call == "normal"
  	] <- "normal"
  	epithelial_metadata$subcluster_id[
  	  epithelial_metadata$normal_cell_call == "unassigned"
  	] <- "unassigned"
  }

  # do the same for expression clusters:
  epithelial_metadata$expression_id <- as.character(
    epithelial_metadata$expression_id
  )
  epithelial_metadata$expression_id[
    epithelial_metadata$normal_cell_call == "normal"
  ] <- "normal"
  epithelial_metadata$expression_id[
    epithelial_metadata$normal_cell_call == "unassigned"
  ] <- "unassigned"

  saveRDS(
    epithelial_metadata,
    paste0(Robject_dir, "/3.epithelial_metadata_with_normals.Rdata")
  )

} else {

  epithelial_metadata <- readRDS(
    paste0(Robject_dir, "/3.epithelial_metadata_with_normals.Rdata")
  )

}

if (remove_normals) {
  epithelial_metadata <- epithelial_metadata[
    epithelial_metadata$normal_cell_call != "normal" & 
    epithelial_metadata$normal_cell_call != "unassigned",
  ]
  epithelial_heatmap <- epithelial_heatmap[
    rownames(epithelial_heatmap) %in% epithelial_metadata$cell_ids,
  ]

  saveRDS(
    epithelial_metadata,
    paste0(Robject_dir, "/3.epithelial_metadata_without_normals.Rdata")
  )
}


################################################################################
### 4. Filter and merge subclusters ###
################################################################################

# remove expression clusters with less than n cells:
remove_expression_clusters <- names(
  which(
    lapply(
      split(epithelial_metadata, epithelial_metadata$expression_id),
      nrow
    ) < min_cluster_cells
  )
)
epithelial_metadata <- epithelial_metadata[
  !(epithelial_metadata$expression_id %in% remove_expression_clusters),
]
# update heatmap:
temp_heatmap <- epithelial_heatmap %>%
  rownames_to_column() %>%
  filter(rownames(epithelial_heatmap) %in% epithelial_metadata$cell_ids) %>%
  column_to_rownames()

# remove subclusters with less than n cells:
remove_subclusters <- names(
  which(
    lapply(
      split(epithelial_metadata, epithelial_metadata$subcluster_id),
      nrow
    ) < min_cluster_cells
  )
)
epithelial_metadata <- epithelial_metadata[
  !(epithelial_metadata$subcluster_id %in% remove_subclusters),
]
# update heatmap:
temp_heatmap <- epithelial_heatmap %>%
  rownames_to_column() %>%
  filter(rownames(epithelial_heatmap) %in% epithelial_metadata$cell_ids) %>%
  column_to_rownames()

# merge subclusters with greater than r correlation:
if (subcluster_merge) {

  epithelial_metadata <- merge_subclusters(
    epithelial_heatmap,
    epithelial_metadata,
    merge_thresh,
    merge_diff_prop
  )
 
}

if (subcluster_method == "random_trees") {
 
  # update subcluster ids:
  no_subclusters <- length(
    unique(
      epithelial_metadata$subcluster_id[
        grep("CNV", epithelial_metadata$subcluster_id)
      ]
    )
  )
  epithelial_metadata$subcluster_id <- factor(
    epithelial_metadata$subcluster_id
  )
  levels(epithelial_metadata$subcluster_id) <- c(
    paste0("CNV_", 1:no_subclusters), "normal", "unassigned"
  )

} else {

  # update subcluster ids:
  no_subclusters <- length(
    unique(
      epithelial_metadata$subcluster_id[
        grep("CNV", epithelial_metadata$subcluster_id)
      ]
    )
  )
  epithelial_metadata$subcluster_id <- factor(
    epithelial_metadata$subcluster_id
  )
  levels(epithelial_metadata$subcluster_id) <- paste0(
    "CNV_", 1:no_subclusters
  )

}

# reassign expression cluster names in numerical order:
epithelial_metadata$expression_id <- factor(
  epithelial_metadata$expression_id,
  levels = naturalsort(unique(epithelial_metadata$expression_id))
)

if (!remove_normals) {

  levels(epithelial_metadata$expression_id) <- c(
    paste0(
      "Expression_",
      1:(length(levels(epithelial_metadata$expression_id))-2)
    ),
    "normal",
    "unassigned"
  )

} else {

  levels(epithelial_metadata$expression_id) <- paste0(
    "Expression_",
    1:length(levels(epithelial_metadata$expression_id))
  )

}
 
epithelial_metadata$subcluster_id <- as.factor(epithelial_metadata$subcluster_id)
if (!remove_normals) {
  levels(epithelial_metadata$subcluster_id)[
    grep("CNV", levels(epithelial_metadata$subcluster_id))
  ] <- paste0(
    "CNV_",
    1:length(levels(epithelial_metadata$subcluster_id)[
      grep("CNV", levels(epithelial_metadata$subcluster_id))
    ])
  )
} else {
  levels(epithelial_metadata$subcluster_id) <- paste0(
    "CNV_",
    1:length(levels(epithelial_metadata$subcluster_id))
  )
}


################################################################################
### 5. Order metadata ###
################################################################################

if (order_by == "CNV") {

  # if normals annotated, change the order of clusters to plot most normal first:
  if (normal_annot & !remove_normals) {

    # split metadata by normal call:
    normal_split <- split(
      epithelial_metadata, epithelial_metadata$normal_cell_call
    )
    # order by normal, unassigned then cancer:
    normal_split <- list(
      normal_split$normal,
      normal_split$unassigned,
      normal_split$cancer
    )

    # split metadata by cluster for ordering:
    subcluster_split <- lapply(normal_split, function(x) {      
      return(x[naturalorder(x$subcluster_id),])
    })

    # rbind metadata together in that order:
    epithelial_metadata <- do.call("rbind", subcluster_split)

  } else {

    epithelial_metadata <- epithelial_metadata[
      naturalorder(epithelial_metadata$subcluster_id),
    ]

  }

  epithelial_metadata$subcluster_id <- factor(
    epithelial_metadata$subcluster_id,
    levels = naturalsort(unique(epithelial_metadata$subcluster_id))
  )

} else if (order_by == "expression") {

  # if normals annotated, change the order of clusters to plot most normal first:
  if (!remove_normals) {

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

}

# adjust order of heatmap:
epithelial_heatmap <- epithelial_heatmap[epithelial_metadata$cell_ids,]

print(paste0(
  "Is heatmap order the same as metadata? ", 
  identical(rownames(epithelial_heatmap), epithelial_metadata$cell_ids)
))

# save final objects:
if (remove_normals) {
  if (
    !file.exists(
      paste0(
        Robject_dir, "/5b.final_epithelial_heatmap_without_normals.Rdata"
      )
    )
  ) {
    saveRDS(
      epithelial_heatmap,
      paste0(Robject_dir, "/5a.final_epithelial_heatmap_without_normals.Rdata")
    )
    saveRDS(
      epithelial_metadata,
      paste0(
        Robject_dir, "/5b.final_epithelial_metadata_without_normals.Rdata"
      )
    )
  }
} else {
  if (
    !file.exists(
      paste0(
        Robject_dir, "/4a.final_epithelial_heatmap_with_normals.Rdata"
        )
      )
    ) {
    saveRDS(
      epithelial_heatmap,
      paste0(Robject_dir, "/4a.final_epithelial_heatmap_with_normals.Rdata")
    )
    saveRDS(
      epithelial_metadata,
      paste0(Robject_dir, "/4b.final_epithelial_metadata_with_normals.Rdata")
    )
  }
}


################################################################################
### 6. Create heatmap and annotations ###
################################################################################

# create expression cluster annotation:
if (expression_annot) {

  # define cluster annotation colours:
  expr_number <- length(unique(epithelial_metadata$expression_id))
  expr_cols <- expr_cols[1:expr_number]
  names(expr_cols) <- levels(epithelial_metadata$expression_id)

  expr_annot_df <- subset(epithelial_metadata, select = expression_id)

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
  subcluster_cols <- subcluster_cols[
    1:length(unique(subcluster_annot_df$subcluster_id))
  ]
  names(subcluster_cols) <- unique(subcluster_annot_df$subcluster_id)
  CNV_subcluster_annot <- Heatmap(as.matrix(subcluster_annot_df), 
    col = subcluster_cols, 
    name = "CNV\nsubclusters", width = unit(6, "mm"), 
    show_row_names = F, show_column_names = F, 
    show_heatmap_legend = F,
    heatmap_legend_param = list(
      title = "",
      labels_gp = gpar(fontsize = 16)
    )
  )
}

# create expression cluster annotation:
if (normal_annot) {
  # define cluster annotation colours:
  normal_cols <- c(
    "normal" = "#1B7837", 
    "unassigned" = "#E7E4D3", 
    "cancer" = "#D95F02"
  )
 
  normal_annot_df <- subset(epithelial_metadata, select = normal_cell_call)
  normal_call_annot <- Heatmap(
    as.matrix(normal_annot_df), 
    col = normal_cols,
    name = "Normal\nannotation", 
    width = unit(4, "mm"), 
    show_row_names = F, show_column_names = F,
    show_heatmap_legend = FALSE
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

  CNA_annot <- rowAnnotation(
    CNA = anno_barplot(
      epithelial_metadata$CNA_value, name = "nGene",
      gp = gpar(
        col = "#D95F02", 
        width = unit(4, "cm")
      ), 
      border = FALSE, 
      which = "row", 
      axis = F
    ), show_annotation_name = FALSE
  )
  CNA_annot@name <- "CNA"

  if (plot_type == "normals_annotated") {
    corr_annot <- rowAnnotation(
      corr = anno_barplot(
        epithelial_metadata$cor.estimate, name = "nGene",
        gp = gpar(
          col = "#430F82", 
          width = unit(4, "cm")
        ), 
        border = FALSE, 
        which = "row", 
        axis = F
      ), show_annotation_name = FALSE
    )
    corr_annot@name <- "corr"
  }
  
}

if (array_CNVs) {
  # create array CNV annotation:
  all_array_CNVs <- read.table(paste0(ref_dir, "all_array_CNVs.txt"))
  colnames(all_array_CNVs) <- gsub("CID4499_1", "CID44991", colnames(all_array_CNVs))
  if (any(colnames(all_array_CNVs) %in% sample_name)) {
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
  }
}


################################################################################
### 7. Create epithelial and reference heatmap ###
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
  na_less_vector <- unlist(epithelial_heatmap)
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
  labels = c("loss", "", "gain"),
  col_fun = heatmap_cols, 
  title = "CNV signal", 
  direction = "horizontal",
  grid_height = unit(2.5, "cm"),
  grid_width = unit(0.1, "cm"),
  labels_gp = gpar(fontsize = 20),
  title_gp = gpar(fontsize = 26, fontface = "plain")
)

# organise chromosome labels:
chr_labels <- chr_data$lab_pos
names(chr_labels) <- gsub("chr", "", names(chr_labels))
names(chr_labels) <- gsub("^1$", "chr1", names(chr_labels))
names(chr_labels) <- gsub("21", "\n21", names(chr_labels))


################################################################################
### 8. Create and plot annotated heatmap ###
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
      ht_list <- normal_call_annot + ht_list
    }
    if (QC_annot) {
      ht_list <- ht_list + nUMI_annot + nGene_annot + CNA_annot + corr_annot
    }
  
  }

}

annotated_heatmap <- grid.grabExpr(
  draw(ht_list, gap = unit(6, "mm"), heatmap_legend_side = "left")
)
dev.off()

if (plot_references) {

  pdf(paste0(filename, ".pdf"), 
    height = 13, width = 20) 
  
    grid.newpage()

    # print reference heatmap:
    pushViewport(viewport(x = 0.56, y = 0.68, width = 0.765, height = 0.3,
      just = "bottom"))
      #grid.rect()
      # plot heatmap:
      grid.draw(grid_ref_heatmap)
      decorate_heatmap_body("ref_hm", {
        for ( e in 1:length(chr_data$end_pos) ) {
          grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), 
            gp = gpar(lwd = 1, col = "#383838"))
        }
      })
    popViewport()

    # label reference heatmap:
    pushViewport(viewport(x = 0.1, y = 0.95, 
      width = unit(2, "cm"), height = unit(0.5, "cm")))
      #grid.rect()
      grid.text("Reference cells:", gp=gpar(fontsize=20))
    popViewport()

    # print epithelial heatmap:
    pushViewport(viewport(x = 0.5, y = 0, width = 1, height = 0.75,
      just = "bottom"))
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
 
      # print heatmap legend:
      pushViewport(viewport(x = unit(2, "cm"), y = unit(19, "cm"), width = unit(0.1, "cm"), 
        height = unit(0.4, "cm"), just = c("right", "bottom")))
        draw(lgd, x = unit(0.1, "cm"), y = unit(0.1, "cm"), just = c("left", "bottom"))
      popViewport()

      # print subcluster legend:
      if (subcluster_legend) {
        create_legend("CNV subcluster", epithelial_metadata$subcluster_id, subcluster_cols) 
      }
  
      # print expression legend:
      if (expression_legend) {
        create_legend("Expression", epithelial_metadata$expression_id, expr_cols) 
      }
   
      # label QC annotations:
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
  
  pdf(paste0(filename, ".pdf"), 
    height = 13, width = 20) 
  
    grid.newpage()
  
      # plot heatmap:
      pushViewport(viewport(x = 0.155, y = 0.04, width = 0.82, height = 0.95, 
        just = c("left", "bottom")))
        #grid.rect()
        grid.draw(annotated_heatmap)
        decorate_heatmap_body("hm", {
          for ( e in 1:length(chr_data$end_pos) ) {
          grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), 
            gp = gpar(lwd = 1, col = "#383838"))
          grid.text(names(chr_labels)[e], chr_labels[e], 
              unit(0, "npc") + unit(-3.5, "mm"), gp=gpar(fontsize=22))
          }
        })
      popViewport()
  
      # plot heatmap legend:
      pushViewport(viewport(x = unit(1.7, "cm"), y = unit(28.3, "cm"), width = unit(0.1, "cm"), 
        height = unit(0.4, "cm"), just = c("right", "bottom")))
        draw(lgd, x = unit(0.1, "cm"), y = unit(0.1, "cm"), just = c("left", "bottom"))
      popViewport()
      
      # print subcluster legend:
      if (subcluster_legend) {
        pushViewport(viewport(x = 0.075, y = 0.476, width = unit(6, "cm"), 
          height = unit(10, "cm"), just="bottom"))
          create_legend(
            "CNV subcluster", 
            epithelial_metadata$subcluster_id,
            sort_labs = FALSE,
            subcluster_cols,
            lib_loc
          ) 
        popViewport()
      }
  
      # print expression legend:
      if (expression_legend) {
        pushViewport(viewport(x = 0.075, y = 0.15, width = unit(6, "cm"), 
          height = unit(10, "cm"), just="bottom"))
        #grid.rect()
          if (file.exists(
            paste0(seurat_dir, "05_seurat_object_epithelial_reclustered.Rdata")
          )) {
            create_legend(
              "Expression",
              epithelial_metadata$expression_id,
              sort_labs = TRUE,
              expr_cols,
              lib_loc
            ) 
          } else {
            create_legend(
              "Orig. expression",
              epithelial_metadata$expression_id,
              sort_labs = TRUE,
              expr_cols,
              lib_loc
            ) 
          }
        popViewport()
      }

      # print normal annot legend:
      if (normal_legend) {
        pushViewport(viewport(x = 0.075, y = 0.16, width = unit(6, "cm"), 
          height = unit(10, "cm"), just="bottom"))
          create_legend(
            "Normal vs cancer", 
            epithelial_metadata$normal_cell_call, 
            sort_labs = FALSE,
            normal_cols,
            lib_loc
          ) 
        popViewport()
      }
   
      # label annotations:
      pushViewport(viewport(x=0.86, y=0.08, width = 0.1, height = 0.1, 
        just = "top"))
        grid.text("nUMI", rot=65, gp=gpar(fontsize=18))
      popViewport()
      pushViewport(viewport(x=0.89, y=0.08, width = 0.1, height = 0.1, 
        just = "top"))
        grid.text("nGene", rot=65, gp=gpar(fontsize=18))
      popViewport()
      pushViewport(viewport(x=0.92, y=0.08, width = 0.1, height = 0.1, 
        just = "top"))
        grid.text("CNA", rot=65, gp=gpar(fontsize=18))
      popViewport()
      pushViewport(viewport(x=0.95, y=0.08, width = 0.1, height = 0.1, 
        just = "top"))
        grid.text("Corr.", rot=65, gp=gpar(fontsize=18))
      popViewport()

      # draw array CNV heatmap:
      if (array_CNVs) {
        pushViewport(viewport(x = 0.155, y = 0.86, 
          width = 0.8657, height = 0.13, just = c("right", "bottom")))
        grid.draw(grid_array_heatmap)
        popViewport()
      }
      
  dev.off()

}

print(paste0("Heatmap created, output in ", plot_dir))
  
  