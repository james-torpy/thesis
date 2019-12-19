#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

project_name <- "thesis"
subproject_name <- "Figure_2.3_simulations"
sample_name <- args[1]
numcores <- as.numeric(args[2])
subset_data <- as.logical(args[3])
include_t_cells <- as.logical(args[4])
nUMI_threshold <- as.numeric(args[5])
nGene_threshold <- as.numeric(args[6])
random_proportion <- as.numeric(args[7])

#numcores <- 10
#sample_name <- "CID4520N"
#subset_data <- FALSE
#include_t_cells <- TRUE
#analysis_mode <- "samples"
#nUMI_threshold <- 25000
#nGene_threshold <- 5000
#random_proportion <- as.numeric("0.01")


print(paste0("Project name = ", project_name))
print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("Subset data? ", as.character(subset_data)))
print(paste0("Number cores = ", numcores))
print(paste0("Include T cells? ", as.character(include_t_cells)))
print(paste0("nUMI threshold =  ", nUMI_threshold))
print(paste0("nGene threshold =  ", nGene_threshold))
print(paste0("Proportions of genes to randomise =  ", random_proportion))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.5dev/"
library(cluster, lib.loc = lib_loc)
library(Seurat)
library(infercnv, lib.loc = lib_loc)


home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", 
  project_name, "/", subproject_name, "/")
results_dir <- seurat_path <- paste0(project_dir, "results/")
func_dir <- paste0(project_dir, "scripts/functions/")
ref_dir <- paste0(project_dir, "refs/")
in_dir <- paste0(project_dir, "raw_files/seurat_objects/", 
  sample_name, "/")

if (include_t_cells) {
  infercnv_dir <- paste0(results_dir, "infercnv/t_cells_included/",
    sample_name, "/filtered_by_coverage/")
  out_dir <- paste0(results_dir, "infercnv/t_cells_included/",
    sample_name, "/", random_proportion, "_genes_randomly_permutated/")
} else {
  infercnv_dir <- paste0(results_dir, "infercnv/t_cells_excluded/",
    sample_name, "/filtered_by_coverage/")
  out_dir <- paste0(results_dir, "infercnv/t_cells_excluded/",
    sample_name, "/", random_proportion, "_genes_randomly_permutated/")
}

input_dir <- paste0(out_dir, "input_files/")
system(paste0("mkdir -p ", input_dir))

setwd(out_dir)

print(paste0("Sample directory = ", in_dir))
print(paste0("Output directory = ", out_dir))

print(
  paste0(
    "Randomly permutating ", random_proportion, " of genes for ", sample_name
  )
)


################################################################################
### 0. Define functions ###
################################################################################

prepare_infercnv_metadata <- dget(paste0(func_dir, "prepare_infercnv_metadata.R"))


################################################################################
### 1. Generate initial input matrix and metadata files ###
################################################################################

# load seurat object:
seurat_10X <- readRDS(paste0(in_dir, "03_seurat_object_processed.Rdata"))
Idents(seurat_10X) <- seurat_10X@meta.data$PC_A_res.1

# create raw matrix input file and subset if necessary:
count_df <- as.matrix(GetAssayData(seurat_10X , slot = "counts"))
if (subset_data) {
  count_df <- count_df[1:500, 1:500]
}
print(
  paste0(
    "Dimensions of count df = ", paste(as.character(dim(count_df)), collapse=",")
  )
)

# create metadata df:
print("Creating inferCNV metadata file...")
infercnv_metadata <- prepare_infercnv_metadata(seurat_10X, subset_data=subset_data, 
  count_df, for_infercnv=T)
seurat_10X <- infercnv_metadata$seurat
print(paste0("Cell types are: ", unique(infercnv_metadata$metadata$cell_type)))

# only keep cells in metadata df:
print(paste0("No cells in count df before filtering for those in metadata df = ", 
    ncol(count_df)))
count_df <- count_df[,colnames(count_df) %in% infercnv_metadata$metadata$cell_ids]
print(paste0("No cells in count df after filtering for those in metadata df = ", 
    ncol(count_df)))

# if necessary, remove T cells from analysis:
if (!include_t_cells) {
  infercnv_metadata$metadata <- infercnv_metadata$metadata[
    grep("[t,T][-_][c,C]ell", infercnv_metadata$metadata$cell_type, invert=T),
  ]
}
# collapse all stromal cells into 'stromal' cell type:
infercnv_metadata$metadata$cell_type[
  grep("pithelial", infercnv_metadata$metadata$cell_type, invert=T)
] <- "Stromal"
# collapse all epithelial cells into 'epithelial' cell type:
infercnv_metadata$metadata$cell_type[
  grep("pithelial", infercnv_metadata$metadata$cell_type)
] <- "Epithelial"

# filter out cells with nUMI < nUMI_threshold and nGene < nGene_threshold
print(paste0("Number of cells before filtering out low coverage: ",
  ncol(count_df)))
QC <- data.frame(
  row.names = names(Idents(seurat_10X)),
  nUMI = seurat_10X@meta.data$nCount_RNA,
  nGene = seurat_10X@meta.data$nFeature_RNA
)
QC <- QC[colnames(count_df),]

cells_to_remove <- colnames(count_df)[QC$nUMI < nUMI_threshold & QC$nGene < nGene_threshold]
print(paste0("Number cells initially in to remove list = ", length(cells_to_remove)))
# only filter out epithelial cells:
stromal_cells <- infercnv_metadata$metadata$cell_id[
  infercnv_metadata$metadata$cell_type == "Stromal"
]
epithelial_cells <- infercnv_metadata$metadata$cell_id[
  infercnv_metadata$metadata$cell_type == "Epithelial"
]
print(paste0("Number stromal cells = ", length(stromal_cells)))
print(paste0("Number epithelial cells = ", length(epithelial_cells)))
print(paste0("Number stromal cells in cells_to_remove = ",
  length(cells_to_remove[cells_to_remove %in% stromal_cells])))
print(paste0("Number epithelial cells in cells_to_remove = ",
  length(cells_to_remove[cells_to_remove %in% epithelial_cells])))
cells_to_remove <- cells_to_remove[
  !(cells_to_remove %in% stromal_cells)
]
print(paste0("Number cells to remove after removing stromal = ", length(cells_to_remove)))
count_df <- count_df[
  ,!(colnames(count_df) %in% cells_to_remove)
]
print(paste0("Number of cells after filtering out low coverage: ",
  ncol(count_df)))
print(paste0("Number of cells in metadata before removing low coverage cells = ",
  nrow(infercnv_metadata$metadata)))
infercnv_metadata$metadata <- infercnv_metadata$metadata[colnames(count_df),]
print(paste0("Number of cells in metadata after removing low coverage cells = ",
  nrow(infercnv_metadata$metadata)))

epithelial_clusters <- grep("pithelial", unique(infercnv_metadata$metadata$cell_type), value=T)


################################################################################
### 2. Randomly permutate epithelial count_df ###
################################################################################

# load previous InferCNV results to determine which genes weren't filtered out:
infercnv_output <- read.table(paste0(infercnv_dir, 
  "infercnv.12_denoised.observations.txt"))

non_filtered_genes <- rownames(count_df)[
  rownames(count_df) %in% rownames(infercnv_output)
]

# isolate epithelial counts in data frame and save stromal separately:
epithelial_count_df <- count_df[, colnames(count_df) %in% epithelial_cells]
stromal_count_df <- count_df[, colnames(count_df) %in% stromal_cells]

# fetch random selection of genes of sepcified proportion and multiply this
# gene by 3 for all epithelial cells:
total_permutations <- floor(
  as.numeric(random_proportion)*length(non_filtered_genes)
)
permutated_gene_indices <- sample(
  1:length(non_filtered_genes), 
  total_permutations, 
  replace=FALSE
)
permutated_genes <- non_filtered_genes[permutated_gene_indices]
permutated_df <- epithelial_count_df
permutated_df[
  rownames(permutated_df) %in% permutated_genes,
] <- permutated_df[
  rownames(permutated_df) %in% permutated_genes,
]*3

# add back to stromal_count_df:
print(
  paste0(
    "Number of cells before adding stromal_count_df = ", ncol(permutated_df)
  )
)
permutated_df <- cbind(permutated_df, stromal_count_df)
print(
  paste0(
    "Number of cells after adding stromal_count_df = ", ncol(permutated_df)
  )
)

# check gene counts have been changed:
print(
  paste0(
    "Are permutated gene counts different from original? ", 
    !(identical(count_df[permutated_genes,], permutated_df[permutated_genes,]))
  )
)  


################################################################################
### 3. Save InferCNV input files ###
################################################################################

# if no epithelial clusters present, abort:
if (length(epithelial_clusters) < 1) {
  print(paste0("No epithelial/myoepithelial clusters detected, aborting..."))
} else {
  # write count, metadata files and new seurat object:
  if (!file.exists(paste0(input_dir, "input_matrix.txt"))) {
    print("Creating inferCNV raw counts file...")
    write.table(permutated_df, paste0(input_dir, "input_matrix.txt"), quote=F,
    sep="\t", col.names=T, row.names=T)
  }

  if (!file.exists(paste0(input_dir, "metadata.txt"))) {
    write.table(infercnv_metadata$metadata, paste0(input_dir, "metadata.txt"), 
      quote=F, sep="\t", col.names=F, row.names=F)
    write.table(infercnv_metadata$number_per_group, paste0(input_dir, 
      "number_per_group.txt"), quote=F, col.names=F, row.names=F, sep="\t")
    saveRDS(seurat_10X, paste0(in_dir, "04_seurat_object_annotated.Rdata"))
  }
  
  # define normals which will act as InferCNV reference cells:
  normals <- grep(
    "[e,E]pithelial|[m,M]yoepithelial|CAF|[u,U]nassigned|[u,U]nknown|[t,T]umour|[t,T]umor", 
    unique(infercnv_metadata$metadata$cell_type[
      infercnv_metadata$metadata$cell_ids %in% colnames(permutated_df)
    ]), value=T, 
    invert=T
  )


  ################################################################################
  ### 4. Run InferCNV ###
  ################################################################################
  
  print(paste0("Normal is: ", normals))

  print("Creating inferCNV object...")
  raw_path <- paste0(input_dir, "input_matrix.txt")
  annotation_path <- paste0(input_dir, "metadata.txt")
  gene_path <- paste0(ref_dir, "infercnv_gene_order.txt")
  initial_infercnv_object <- CreateInfercnvObject(
    raw_counts_matrix=raw_path,
    annotations_file=annotation_path,
    delim="\t",
    gene_order_file=gene_path,
    ref_group_names=normals
  )

  print("InferCNV object created, running inferCNV...")
  system.time(
    infercnv_output <- try(
      infercnv::run(
        initial_infercnv_object,
        num_threads=numcores-1,
        out_dir=".",
        cutoff=0.1,
        window_length=101,
        max_centered_threshold=3,
        plot_steps=F,
        denoise=T,
        sd_amplifier=1.3,
        analysis_mode = "samples"
      )
    )
  )
}

# remove temp files:
system(paste0("rm ", out_dir, "/*dat"))
system(paste0("rm ", out_dir, "/*preliminary*"))
system(paste0("rm ", out_dir, "/01_*"))
system(paste0("rm ", out_dir, "/02_*"))
system(paste0("rm ", out_dir, "/03_*"))
system(paste0("rm ", out_dir, "/04_*"))
system(paste0("rm ", out_dir, "/05_*"))
system(paste0("rm ", out_dir, "/06_*"))
system(paste0("rm ", out_dir, "/07_*"))
system(paste0("rm ", out_dir, "/08_*"))
system(paste0("rm ", out_dir, "/09_*"))
system(paste0("rm ", out_dir, "/10_*"))
system(paste0("rm ", out_dir, "/11_*"))
system(paste0("rm ", out_dir, "/12_*"))

