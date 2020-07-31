#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

project_name <- "thesis"
subproject_name <- "Figure_2.4_min_length_and_gap"
sample_name <- args[1]
print(paste0("Sample name = ", sample_name))
numcores <- as.numeric(args[2])
print(paste0("Numcores = ", numcores))
nUMI_threshold <- as.numeric(args[3])
print(paste0("nUMI_threshold = ", nUMI_threshold))
nGene_threshold <- as.numeric(args[4])
print(paste0("nGene_threshold = ", nGene_threshold))
include_t_cells <- as.logical(args[5])
print(paste0("Include T-cells? ", include_t_cells))
analysis_mode <- args[6]
print(paste0("Analysis mode = ", analysis_mode))

#sample_name <- "CID4520N"
#numcores <- 10
#nUMI_threshold <- 25000
#nGene_threshold <- 5000
#include_t_cells <- TRUE
#analysis_mode <- "samples"

print(paste0("Project name = ", project_name))
print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("Number cores = ", numcores))
print(paste0("Include T cells? ", as.character(include_t_cells)))

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

if (include_t_cells) {
  out_path <- paste0(results_dir, "infercnv/t_cells_included/")
  Robject_dir <- paste0(out_path, sample_name, "/Rdata/")
  plot_dir <- paste0(out_path, sample_name, "/plots/")
  out_dir <- paste0(out_path, sample_name, "/", analysis_mode, "_mode/")
} else {
  out_path <- paste0(results_dir, "infercnv/t_cells_included/")
  Robject_dir <- paste0(out_path, sample_name, "/Rdata/")
  plot_dir <- paste0(out_path, sample_name, "/plots/")
  out_dir <- paste0(out_path, sample_name, "/", analysis_mode, "_mode/")
}
input_dir <- paste0(out_path, sample_name, "/input_files/")

system(paste0("mkdir -p ", out_dir))
system(paste0("mkdir -p ", Robject_dir))
system(paste0("mkdir -p ", plot_dir))
system(paste0("mkdir -p ", input_dir))

setwd(out_dir)

print(paste0("Sample directory = ", input_dir))
print(paste0("Reference directory = ", ref_dir))
print(paste0("Output directory = ", out_dir))

print(paste0("Running InferCNV identify normals pipeline on ", sample_name ,
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
  # create metadata df:
  print("Creating inferCNV metadata file...")
  infercnv_metadata <- prepare_infercnv_metadata(seurat_10X, subset_data=subset_data, 
    count_df, for_infercnv=T)
  seurat_10X <- infercnv_metadata$seurat
  print(paste0("Cell types are: ", unique(infercnv_metadata$metadata$cell_type)))
  saveRDS(infercnv_metadata, paste0(Robject_dir, "/1c.original_infercnv_metadata.Rdata"))
  saveRDS(seurat_10X, paste0(in_dir, "04_seurat_object_annotated.Rdata"))
  
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

  # create density plots of nUMI and nGene:
  QC <- data.frame(
    row.names = colnames(epithelial_df),
    nUMI = apply(epithelial_df, 2, sum),
    nGene = apply(epithelial_df, 2, function(x) length(x[x!=0]))
  )
  QC <- QC[colnames(epithelial_df),]
  nUMI_density_plot <- density(QC$nUMI)
  pdf(paste0(plot_dir, "nUMI_density_plot.pdf"))
    plot(nUMI_density_plot, main=NA, xlab = "nUMI")
  dev.off()
  png(paste0(plot_dir, "nUMI_density_plot.png"))
    plot(nUMI_density_plot, main=NA, xlab = "nUMI")
  dev.off()
  nGene_density_plot <- density(QC$nGene)
  pdf(paste0(plot_dir, "nGene_density_plot.pdf"))
    plot(nGene_density_plot, main=NA, xlab = "nGene")
  dev.off()
  png(paste0(plot_dir, "nGene_density_plot.png"))
    plot(nGene_density_plot, main=NA, xlab = "nGene")
  dev.off()
  log_nUMI_density_plot <- density(log10(QC$nUMI))
  pdf(paste0(plot_dir, "log10_nUMI_density_plot.pdf"))
    plot(log_nUMI_density_plot, main=NA, xlab = "log10 nUMI")
  dev.off()
  png(paste0(plot_dir, "log10_nUMI_density_plot.png"))
    plot(log_nUMI_density_plot, main=NA, xlab = "log10 nUMI")
  dev.off()
  log_nGene_density_plot <- density(log10(QC$nGene))
  pdf(paste0(plot_dir, "log10_nGene_density_plot.pdf"))
    plot(log_nGene_density_plot, main=NA, xlab = "log10 nGene")
  dev.off()
  png(paste0(plot_dir, "log10_nGene_density_plot.png"))
    plot(log_nGene_density_plot, main=NA, xlab = "log10 nGene")
  dev.off()
  # filter out cells with nUMI < nUMI_threshold and nGene < nGene_threshold
  print(paste0("Number of cells before filtering out low coverage: ",
    nrow(QC)))
  cells_to_keep <- rownames(QC)[QC$nUMI > nUMI_threshold & QC$nGene > nGene_threshold]
  print(paste0("Number of cells after filtering out low coverage: ",
    length(cells_to_keep)))
  epithelial_df <- epithelial_df[
    ,colnames(epithelial_df) %in% cells_to_keep
  ]
  saveRDS(epithelial_df, paste0(Robject_dir, "/1a.original_epithelial_df.Rdata"))
  
  p <- ggplot(QC, aes(x=nUMI, y=nGene))
  p <- p + geom_point()
  p <- p + xlab("nUMI")
  p <- p + ylab("nGene")
  p <- p + theme(legend.title = element_blank())
  pdf(paste0(plot_dir, "QC_quad_plot.pdf"), 
    width = 10, height = 6)
    print(p)
  dev.off()
  png(paste0(plot_dir, "QC_quad_plot.png"), 
    width = 450, height = 270)
    print(p)
  dev.off()
  # create total count density quad plot:
  print("Determining total count...")
  total_counts <- apply(epithelial_df, 2, sum)
  total_count_density_plot <- density(total_counts, bw="SJ")
  pdf(paste0(plot_dir, "total_count_density_plot.pdf"))
    plot(total_count_density_plot, main=NA, xlab = "Total counts")
  dev.off()
  png(paste0(plot_dir, "total_count_density_plot.png"))
    plot(total_count_density_plot, main=NA, xlab = "Total counts")
  dev.off()
  # create log10 total count density quad plot:
  log_total_count_density_plot <- density(log10(total_counts), bw="SJ")
  pdf(paste0(plot_dir, "log10_total_count_density_plot.pdf"))
    plot(log_total_count_density_plot, main=NA, xlab = "Total counts")
  dev.off()
  png(paste0(plot_dir, "log10_total_count_density_plot.png"))
    plot(log_total_count_density_plot, main=NA, xlab = "Total counts")
  dev.off()
  
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


################################################################################
### 2. Run InferCNV ###
################################################################################

# define normals which will act as InferCNV reference cells:
normals <- "Non_epithelial"
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
if (exists("initial_infercnv_object")) {
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
        analysis_mode = analysis_mode
      )
    )
  )
} else {
  print("InferCNV object not created properly, check and rerun...")
}


## remove temp files:
#system(paste0("rm ", out_dir, "/*dat"))
#system(paste0("rm ", out_dir, "/*preliminary*"))
#system(paste0("rm ", out_dir, "/01_*"))
#system(paste0("rm ", out_dir, "/02_*"))
#system(paste0("rm ", out_dir, "/03_*"))
#system(paste0("rm ", out_dir, "/04_*"))
#system(paste0("rm ", out_dir, "/05_*"))
#system(paste0("rm ", out_dir, "/06_*"))
#system(paste0("rm ", out_dir, "/07_*"))
#system(paste0("rm ", out_dir, "/08_*"))
#system(paste0("rm ", out_dir, "/09_*"))
#system(paste0("rm ", out_dir, "/10_*"))
#system(paste0("rm ", out_dir, "/11_*"))
#system(paste0("rm ", out_dir, "/12_*"))

