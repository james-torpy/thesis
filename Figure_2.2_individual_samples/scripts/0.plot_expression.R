
args = commandArgs(trailingOnly=TRUE)

project_name <- "thesis"
subproject_name <- "Figure_2.2_individual_samples"
sample_name <- args[1]
res <- args[2]
epi_res <- args[3]
garnett_slot <- args[4]
nUMI_threshold <- as.numeric(args[5])
nGene_threshold <- as.numeric(args[6])

project_name <- "thesis"
subproject_name <- "Figure_2.2_individual_samples"
sample_name <- "CID4463"
res <- "PC_C_res.1"
PC <- "C"
epi_res <- "PC_C_res.1"
epi_PC <- "C"
garnett_slot <- "garnett_call_ext_major"
nUMI_threshold <- as.numeric("3200")
nGene_threshold <- as.numeric("700")

print(paste0("Project name = ", project_name))
print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("Cluster resolution = ", res))
print(paste0("Epithelial reclustered resolution = ", res))
print(paste0("Garnett slot specified = ", garnett_slot))
print(paste0("nUMI filter = ", manual_epithelial))
print(paste0("nGene filter = ", exclude_clusters))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(Seurat)
library(naturalsort, lib.loc = lib_loc)
library(RColorBrewer)
library(ggplot2)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/", 
  subproject_name, "/")
results_dir <- paste0(project_dir, "results/")
ref_dir <- paste0(project_dir, "refs/")
in_dir <- paste0(project_dir, "raw_files/seurat_objects/", 
  sample_name, "/")

out_dir <- paste0(results_dir, "seurat/", sample_name, "/")
plot_dir <- paste0(out_dir, "plots/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(out_dir, "tables/")
system(paste0("mkdir -p ", table_dir))
Robject_dir <- paste0(out_dir, "Rdata/")
system(paste0("mkdir -p ", Robject_dir))

print(paste0("Sample directory = ", in_dir))
print(paste0("Output directory = ", out_dir))
print(paste0("Plot directory = ", plot_dir))
print(paste0("Table directory = ", table_dir))


################################################################################
### 0. Load colours ###
################################################################################

expr_cols <- read.table(
  paste0(ref_dir, "expression_colour_palette.txt"),
  header = F,
  stringsAsFactors = F,
  comment.char = ""
)[,1]


################################################################################
### 1. Load and format all cell seurat object ###
################################################################################

# load seurat object:
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
### 3. Load and format epithelial cell seurat object ###
################################################################################

# load seurat object:
seurat_epi <- readRDS(
  paste0(
    in_dir, "05_seurat_object_epithelial_reclustered.Rdata"
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

# update idents will cell type from Garnett:
annotated_idents <- gsub(
  " ", 
  "_", 
  paste0(
    eval(parse(text=paste0("seurat_epi@meta.data$", garnett_slot))), " ", 
    Idents(seurat_epi)
  )
)
Idents(seurat_epi) <- factor(
  as.character(annotated_idents),
  levels = naturalsort(unique(as.character(annotated_idents)))
)


################################################################################
### 2. Plot UMAP of reclustered epithelial cells ###
################################################################################

# order idents so epithelial cells are printed on top:
ident_order <- levels(Idents(seurat_epi))

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

epi_umap <- DimPlot(
  seurat_epi,
  cells = include_cells,
  pt.size = 1.5,
  reduction = paste0("UMAP", epi_PC),
  label = F
) + 
scale_color_manual(
  labels = paste0("Epithelial ", 1:length(ident_order)), 
  values = expr_cols[1:length(ident_order)]
)

png(
  file = paste0(plot_dir, "reclustered_epithelial_cells_UMAP.png"), 
  width = 15, 
  height = 8, 
  res = 300, 
  units = 'in'
)
  print(epi_umap)
dev.off()









