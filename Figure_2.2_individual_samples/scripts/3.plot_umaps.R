#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

project_name <- "thesis"
subproject_name <- "Figure_2.2_individual_samples"
args = commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
subcluster_method <- args[2]
subcluster_p <- args[3]
if (subcluster_p != "none") {
  subcluster_p <- as.numeric(subcluster_p)
}
coverage_filter <- args[4]
remove_artefacts <- args[5]
normals_removed <- as.logical(args[6])
epi_res <- args[7]
epi_PC <- args[8]
garnett_slot <- args[9]
remove_outliers <- as.logical(args[10])

project_name <- "thesis"
subproject_name <- "Figure_2.2_individual_samples"
sample_name <- "CID4463"
subcluster_method <- "random_trees"
subcluster_p <- "0.05"
if (subcluster_p != "none") {
  subcluster_p <- as.numeric(subcluster_p)
}
coverage_filter <- "filtered"
remove_artefacts <- "artefacts_not_removed"
normals_removed <- TRUE
epi_res <- "PC_C_res.1"
epi_PC <-"C"
remove_outliers <- TRUE

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(naturalsort, lib.loc = lib_loc)
library(Seurat)
library(ggplot2)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/",
  subproject_name, "/")
ref_dir <- paste0(project_dir, "/refs/")
func_dir <- paste0(project_dir, "/scripts/functions/")
results_dir <- paste0(project_dir, "results/")
seurat_dir <- paste0(project_dir, "raw_files/seurat_objects/", 
  sample_name, "/")
in_path <- paste0(
  results_dir, "infercnv/", sample_name, "/",
  coverage_filter, "/", subcluster_method, "/p_", subcluster_p, "/", 
  remove_artefacts, "/"
)

Robject_dir <- paste0(in_path, "Rdata/")
system(paste0("mkdir -p ", Robject_dir))
plot_dir <- paste0(in_path, "plots/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(in_path, "tables/")
system(paste0("mkdir -p ", table_dir))

print(paste0("In path = ", in_path))
print(paste0("Reference directory = ", ref_dir))
print(paste0("R function directory = ", func_dir))
print(paste0("R object directory = ", Robject_dir))
print(paste0("Table directory = ", table_dir))
print(paste0("Plot directory = ", plot_dir))


################################################################################
### 0. Define colours ###
################################################################################

expr_cols <- read.table(
  paste0(ref_dir, "expression_colour_palette.txt"),
  header = F,
  stringsAsFactors = F,
  comment.char = ""
)[,1]

normal_col_df <- read.table(
  paste0(ref_dir, "normal_colour_palette.txt"),
  header = F,
  stringsAsFactors = F,
  comment.char = ""
)
normal_cols <- normal_col_df[,2]
names(normal_cols) <- normal_col_df[,1]
normal_cols["unassigned"] <- "#3752D3"

subcluster_cols <- read.table(
  paste0(ref_dir, "CNV_colour_palette.txt"),
  header = F,
  stringsAsFactors = F,
  comment.char = ""
)[,1]


################################################################################
### 1. Load and format data ###
################################################################################

print("Loading metadata df...")


epithelial_metadata <- readRDS(paste0(Robject_dir, 
  "/4b.final_epithelial_metadata_with_normals.Rdata"))

subcluster_meta <- subset(
  epithelial_metadata, 
  select = c(cell_ids, normal_cell_call, subcluster_id)
)

print("Loading seurat objects...")

# load all cell object:
seurat_10X <- readRDS(paste0(in_dir, "03_seurat_object_processed.Rdata"))
  
# choose resolution if needed:
if (res != "none") {
  Idents(seurat_10X) <- eval(
    parse(
      text = paste0("seurat_10X@meta.data$", res)
    )
  )
}

# update idents will cell type from Garnett:
annotated_idents <- gsub(
  " ", 
  "_", 
  paste0(
    eval(parse(text=paste0("seurat_10X@meta.data$", garnett_slot))), " ", 
    Idents(seurat_10X)
  )
)
Idents(seurat_10X) <- factor(
  as.character(annotated_idents),
  levels = naturalsort(unique(as.character(annotated_idents)))
)

# load reclustered epithelial object:
seurat_epi <- readRDS(
  paste0(
    seurat_dir, "05_seurat_object_epithelial_reclustered.Rdata"
  )
)

# choose resolution if needed:
if (epi_res != "none") {
  Idents(seurat_epi) <- eval(
    parse(
      text = paste0("seurat_epi@meta.data$", epi_res)
    )
  )
}


################################################################################
### 2. Plot UMAP of all cells ###
################################################################################

# order idents so epithelial cells are printed on top:
ident_order <- levels(Idents(seurat_10X))
ident_order <- c(
  ident_order[grep("Epithelial", ident_order)],
  ident_order[grep("Epithelial", ident_order, invert = T)]
)

all_umap <- DimPlot(
  seurat_10X,
  cols = expr_cols,
  pt.size = 1.5,
  reduction = paste0("UMAP", PC),
  label = F,
  order = ident_order
)

png(
  file = paste0(plot_dir, "all_cells_UMAP.png"), 
  width = 15, 
  height = 8, 
  res = 300, 
  units = 'in'
)
  print(all_umap)
dev.off()


################################################################################
### 3. Plot UMAP of reclustered epithelial cells without normals ###
################################################################################

# order idents so epithelial cells are printed on top:
ident_order <- levels(Idents(seurat_epi))

if (remove_outliers) {
  # remove outliers - defined as those > 5 standard devs from mean:
  embeddings <- eval(parse(
    text = paste0("seurat_epi@reductions$UMAP", epi_PC, "@cell.embeddings")
  ))
  mean_embeddings <- apply(embeddings, 2, mean)
  std_dev_embeddings <- apply(embeddings, 2,sd)
  
  outliers1 <- rownames(embeddings)[
    embeddings[,1] > (mean_embeddings[1] + (5*std_dev_embeddings[1])) | 
    embeddings[,1] < (mean_embeddings[1] - (5*std_dev_embeddings[1]))
  ]
  if (length(outliers1) < 5) {
    include_cells <- rownames(embeddings)[
      !(rownames(embeddings) %in% outliers1)
    ]
  }
  
  outliers2 <- rownames(embeddings)[
    embeddings[,2] > (mean_embeddings[2] + (5*std_dev_embeddings[2])) | 
    embeddings[,2] < (mean_embeddings[2] - (5*std_dev_embeddings[2]))
  ]
  if (length(outliers2) < 5) {
    include_cells <- include_cells[
      !(include_cells %in% outliers2)
    ]
  }
  
  # remvoe epithelial clusters with <10 cells:
  split_idents <- split(
    Idents(seurat_epi),
    Idents(seurat_epi)
  )
  small_clusters <- levels(Idents(seurat_epi))[
    unlist(
      lapply(split_idents, function(x) length(x) < 10)
    )
  ]
  include_cells <- include_cells[
    include_cells %in% names(Idents(seurat_epi))[
      !(Idents(seurat_epi) %in% small_clusters)
    ]
  ]
} else {
  include_cells <- rownames(seurat_epi@meta.data)
}

include_cells <- rownames(seurat_epi@meta.data)[
  rownames(seurat_epi@meta.data) %in% subcluster_meta$cell_ids
]

no_normals <- include_cells[
  include_cells %in% subcluster_meta$cell_ids[
    subcluster_meta$normal_cell_call == "cancer"
  ]
]

epi_umap <- DimPlot(
  seurat_epi,
  cells = no_normals,
  pt.size = 1.5,
  reduction = paste0("UMAP", epi_PC),
  label = F
) + 
scale_color_manual(
  labels = paste0("Epithelial ", 1:length(ident_order)), 
  values = expr_cols[1:length(ident_order)]
)

png(
  file = paste0(plot_dir, "epithelial_cell_expr_clusters_UMAP.png"), 
  width = 15, 
  height = 8, 
  res = 300, 
  units = 'in'
)
  print(epi_umap)
dev.off()


################################################################################
### 4. Plot UMAP with normals labelled ###
################################################################################

# update idents with subcluster ids:
seurat_epi@meta.data$normal_cell_call <- NA
seurat_epi@meta.data[subcluster_meta$cell_ids,]$normal_cell_call <- 
  as.character(subcluster_meta$normal_cell_call)

Idents(seurat_epi) <- factor(
  as.character(seurat_epi@meta.data$normal_cell_call),
  levels = naturalsort(
    unique(as.character(seurat_epi@meta.data$normal_cell_call))
  )
)

labs <- c("normal", "unassigned", "cancer")
labs <- labs[
  labs %in% gsub(
    "_", " ", unique(
      Idents(seurat_epi)[
        !is.na(Idents(seurat_epi))
      ]
    )
  )
]

epi_umap <- DimPlot(
  seurat_epi,
  cells = include_cells,
  pt.size = 1.5,
  reduction = paste0("UMAP", epi_PC),
  label = F
) + scale_color_manual(values = normal_cols[1:length(labs)])


png(
  file = paste0(plot_dir, "epithelial_cell_normal_call_UMAP.png"), 
  width = 15, 
  height = 8, 
  res = 300, 
  units = 'in'
)
  print(epi_umap)
dev.off()


################################################################################
### 6. Plot UMAP with CNV subclusters labelled ###
################################################################################

# update idents with subcluster ids:
seurat_epi@meta.data$subcluster_ids <- NA
seurat_epi@meta.data[subcluster_meta$cell_ids,]$subcluster_ids <- 
  as.character(subcluster_meta$subcluster_ids)

Idents(seurat_epi) <- factor(
  as.character(seurat_epi@meta.data$subcluster_ids),
  levels = naturalsort(unique(as.character(seurat_epi@meta.data$subcluster_ids)))
)

labs <- gsub(
  "_", " ", unique(
    Idents(seurat_epi)[
      !is.na(Idents(seurat_epi))
    ]
  )
)

epi_umap <- DimPlot(
  seurat_epi,
  cells = include_cells,
  pt.size = 1.5,
  reduction = paste0("UMAP", epi_PC),
  label = F
) + 
scale_color_manual(
  labels = labs, 
  values = subcluster_cols[1:length(labs)]
)

png(
  file = paste0(plot_dir, "epithelial_cell_CNV_subclusters_UMAP.png"), 
  width = 15, 
  height = 8, 
  res = 300, 
  units = 'in'
)
  print(epi_umap)
dev.off()




