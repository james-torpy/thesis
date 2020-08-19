#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

project_name <- "thesis"
subproject_name <- "Figure_2.2_and_2.3_define_denoising_value"
sample_name <- args[1]
print(paste0("Sample name = ", sample_name))

#sample_name <- "CID4520N"

print(paste0("Project name = ", project_name))
print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(Seurat)
library(ggplot2)
library(infercnv, lib.loc=lib_loc)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", 
  project_name, "/", subproject_name, "/")
ref_dir <- paste0(project_dir, "refs/")
results_dir <- paste0(project_dir, "results/")
func_dir <- paste0(project_dir, "scripts/functions/")

in_dir <- paste0(project_dir, "raw_files/seurat_objects/", sample_name, "/")

out_path <- paste0(results_dir, "infercnv/", sample_name, "/normal/")
Robject_dir <- paste0(out_path, "/Rdata/")
plot_dir <- paste0(out_path, "/plots/")
system(paste0("mkdir -p ", Robject_dir))
system(paste0("mkdir -p ", plot_dir))

input_dir <- paste0(out_path, "input_files/")
system(paste0("mkdir -p ", input_dir))

print(paste0("Sample directory = ", input_dir))
print(paste0("Reference directory = ", ref_dir))
print(paste0("InferCNV input file directory = ", input_dir))

print(paste0("Preparing normal breast InferCNV input files for ", sample_name,
  "..."))


################################################################################
### 0. Define functions and choose seeds ###
################################################################################

prepare_infercnv_metadata <- dget(paste0(func_dir, "prepare_infercnv_metadata.R"))


################################################################################
### 1. Fetch counts matrix from normal sample, create original infercnv input 
# files and filter for high coverage cells  ###
################################################################################

if (!file.exists(paste0(Robject_dir, "/1a.original_epithelial_df.Rdata")) | 
  !file.exists(paste0(Robject_dir, "/1b.original_non_epithelial_df.Rdata")) |
  !file.exists(paste0(Robject_dir, "/1c.original_infercnv_metadata.Rdata"))
) {

  # load seurat object:
  seurat_10X <- readRDS(paste0(in_dir, "03_seurat_object_processed.Rdata"))
  
  # create raw matrix input file and subset if necessary:
  count_df <- as.matrix(GetAssayData(seurat_10X , slot = "counts"))
  print(
    paste0(
      "Dimensions of count df = ", paste(as.character(dim(count_df)), collapse=",")
    )
  )

  ######
  #count_df <- count_df[1:6000,]
  ######

  # create metadata df:
  print("Creating inferCNV metadata file...")
  infercnv_metadata <- prepare_infercnv_metadata(seurat_10X, count_df, for_infercnv=T,
    garnett="garnett_call_ext_major")
  seurat_10X <- infercnv_metadata$seurat

  cell_types <- paste(unique(infercnv_metadata$metadata$cell_type), collapse = " ")
  print(paste0("Cell types are: ", cell_types))
  
  saveRDS(infercnv_metadata, paste0(Robject_dir, "/1c.original_infercnv_metadata.Rdata"))
  
  # only keep cells in metadata df:
  print(paste0("No cells in count df before filtering for those in metadata df = ", 
      ncol(count_df)))
  count_df <- count_df[,colnames(count_df) %in% infercnv_metadata$metadata$cell_ids]
  print(paste0("No cells in count df after filtering for those in metadata df = ", 
      ncol(count_df)))
  
  epithelial_df <- count_df[
    ,as.character(
      infercnv_metadata$metadata$cell_ids[
        grep("pithelial", infercnv_metadata$metadata$cell_type)
      ]
    )
  ]

  print(
    paste0(
      "Dimensions of epithelial df = ", paste(
        as.character(dim(epithelial_df)), collapse=","
      )
    )
  )

  non_epithelial_df <- count_df[
    ,as.character(
      infercnv_metadata$metadata$cell_ids[
        grep("pithelial", infercnv_metadata$metadata$cell_type, invert = T)
      ]
    )
  ]
  saveRDS(non_epithelial_df, paste0(Robject_dir, "/1b.original_non_epithelial_df.Rdata"))
  
} else {

  epithelial_df <- readRDS(paste0(Robject_dir, "/1a.original_epithelial_df.Rdata"))
  non_epithelial_df <- readRDS(
    paste0(Robject_dir, "/1b.original_non_epithelial_df.Rdata")
  )
  infercnv_metadata <- readRDS(
    paste0(Robject_dir, "/1c.original_infercnv_metadata.Rdata")
  )

}


################################################################################
### 2. Format annotation and counts matrices ###
################################################################################

if (
  !file.exists(paste0(Robject_dir, "/2a.epithelial_df.Rdata")) | 
  !file.exists(paste0(Robject_dir, "/2b.non_epithelial_df.Rdata")) |
  !file.exists(paste0(Robject_dir, "/2c.gene_annotation.Rdata")) | 
  !file.exists(paste0(input_dir, "input_matrix.txt")) |
  !file.exists(paste0(input_dir, "metadata.txt"))
) {

  # load in gene annotation and determine chromosome lengths:
  gene_annotation <- read.table(
    paste0(ref_dir, "infercnv_gene_order.txt"),
    header = F,
    sep = "\t",
    as.is = T
  )
  colnames(gene_annotation) <- c("symbol", "chromosome", "start", "end")
  
  # subset gene annotation and epithelial_df so they contain the same genes:
  rownames(gene_annotation) <- gene_annotation$symbol

  genes_not_in_gene_annotation <- rownames(epithelial_df)[
    !(rownames(epithelial_df) %in% rownames(gene_annotation))
  ]
  print(paste0(length(genes_not_in_gene_annotation), 
    " genes not present in InferCNV gene annotation but present in ", sample_name, 
    " counts matrix"))
  
  genes_not_in_epithelial_df <- rownames(gene_annotation)[
    !(rownames(gene_annotation) %in% rownames(epithelial_df))
  ]
  print(paste0(length(genes_not_in_epithelial_df), " genes not present in ", sample_name, 
    " counts matrix but present in InferCNV gene annotation"))
  
  print("Removing genes not present in both InferCNV gene annotation and counts matrix...")
  gene_annotation <- gene_annotation[
    rownames(gene_annotation) %in% rownames(epithelial_df),
  ]
  epithelial_df <- epithelial_df[rownames(epithelial_df) %in% rownames(gene_annotation),]

  print(paste0("Gene numbers of InferCNV gene annotation and counts matrix are now ",
    nrow(gene_annotation), " and ", nrow(epithelial_df), " respectively"))

  # remove genes in chromosomes M, X, Y:
  print("Removing genes from chromosomes M, X, Y...")
  gene_annotation <- gene_annotation[
    !(gene_annotation$chromosome %in% c("chrM", "chrX", "chrY")),
  ]
  # number genes in gene annotation:
  gene_annotation$number <- seq(1, nrow(gene_annotation))

  epithelial_df <- epithelial_df[
    gene_annotation$symbol[!(gene_annotation$chromosome %in% c("chrM", "chrX", "chrY"))],
  ]

  non_epithelial_df <- non_epithelial_df[rownames(epithelial_df),]

  combined_df <- cbind(
    epithelial_df,
    non_epithelial_df
  )

  # keep only cells in combined_df for metadata:
  final_metadata <- infercnv_metadata$metadata[
    rownames(infercnv_metadata$metadata) %in% colnames(combined_df),
  ]
  # change non-epithelial cell type to "Non_epithelial":
  final_metadata$cell_type[grep("pithelial", final_metadata$cell_type)] <- "Epithelial"
  final_metadata$cell_type[
    grep("pithelial", final_metadata$cell_type, invert=T)
  ] <- "Non_epithelial"

  print(
    paste0(
      "Gene numbers of InferCNV epithelial, non_epithelial, combined and gene_annotation",
      " dfs are now ", nrow(epithelial_df), ", ", nrow(non_epithelial_df), ", ", 
      nrow(combined_df), " and ", nrow(gene_annotation), " respectively"
    )
  )

  print(
    paste0(
      "Cell numbers of InferCNV epithelial, non_epithelial, combined and metadata dfs are ", 
      ncol(epithelial_df), ", ", ncol(non_epithelial_df), ", ", ncol(combined_df), " and ",
      nrow(final_metadata), " respectively"
    )
  )
  
  # save all dfs:
  saveRDS(epithelial_df, paste0(Robject_dir, "/2a.epithelial_df.Rdata"))
  saveRDS(non_epithelial_df, paste0(Robject_dir, "/2b.non_epithelial_df.Rdata"))
  saveRDS(gene_annotation, paste0(Robject_dir, "/2c.gene_annotation.Rdata"))

  # save normal dataset infercnv input files:
  write.table(combined_df, paste0(input_dir, "input_matrix.txt"), quote=F,
    sep="\t", col.names=T, row.names=T)
  write.table(final_metadata, paste0(input_dir, "metadata.txt"), 
    sep = "\t", quote = F, col.names = F, row.names = F)

} else {

  epithelial_df <- readRDS(paste0(Robject_dir, "/2a.epithelial_df.Rdata"))
  non_epithelial_df <- readRDS(paste0(Robject_dir, "/2b.non_epithelial_df.Rdata"))
  gene_annotation <- readRDS(paste0(Robject_dir, "/2c.gene_annotation.Rdata"))

}


