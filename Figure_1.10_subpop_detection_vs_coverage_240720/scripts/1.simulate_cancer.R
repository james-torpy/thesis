#! /share/ClusterShare/software/contrib/briglo/octR/src/R-3.6.0/builddir/bin/Rscript


###############################################################################
### Simulates cancer datasets taking as inputs number of CNVs, CNV lengths, ###
### gain/loss multipliers and downsamples required                          ###
###############################################################################

args = commandArgs(trailingOnly=TRUE)

project_name <- "thesis"
subproject_name <- "Figure_1.9_subpop_detection_vs_coverage"
sample_name <- args[1]
print(paste0("sample name = ", sample_name))
subset_data <- as.logical(args[2])
print(paste0("subset_data = ", subset_data))
nUMI_threshold <- as.numeric(args[3])
print(paste0("nUMI_threshold = ", nUMI_threshold))
nGene_threshold <- as.numeric(args[4])
print(paste0("nGene_threshold = ", nGene_threshold))
# range of number of CNVs per simulation, chosen randomly
# defines number of CNV loop iterations with one CNV added per iteration:
CNV_no_range <- as.numeric(
  unlist(
    strsplit(
      args[5],
      split = "_"
    )
  )
)
print(paste0("Range of CNVs possible = ", CNV_no_range))

# range of lengths of CNVs, chosen randomly:
CNV_lengths <- as.numeric(
  unlist(
    strsplit(
      args[6],
      split = "_"
    )
  )
)
print(paste0("Possible CNV lengths = ", CNV_lengths))

# vector of possible gain/loss multipliers, chosen randomly:
CNV_multipliers <- as.numeric(
  unlist(
    strsplit(
      args[7],
      split = "_"
    )
  )
)
print(paste0("Possible CNV multipliers = ", CNV_multipliers))

downsample <- as.logical(args[8])
print(paste0("Downsample = ", downsample))

downsample_proportions <- as.numeric(
  unlist(
    strsplit(
      args[9],
      split = "_"
    )
  )
)
print(paste0("Proportions to downsample to = ", downsample_proportions))

simulation_number <- as.numeric(args[10])
print(paste0("Simulation number = ", simulation_number))

noise_cell_no <- as.numeric(args[11])
print(paste0("Noise input cell number = ", noise_cell_no))

project_name <- "thesis"
subproject_name <- "Figure_1.9_subpop_detection_vs_coverage"
#sample_name <- "CID4520N"
#subset_data <- FALSE
#nUMI_threshold <- 25000
#nGene_threshold <- 5000
## range of number of CNVs per simulation, chosen randomly
## defines number of CNV loop iterations with one CNV added per iteration:
#CNV_no_range <- as.numeric(
#  unlist(
#    strsplit(
#      "10_50",
#      split = "_"
#    )
#  )
#)
## range of lengths of CNVs, chosen randomly:
#CNV_lengths <- as.numeric(
#  unlist(
#    strsplit(
#      "20_50_75_100_150_150_200_200_250_250_300_300_350_400_450_500_550_600_650_700_750_800",
#      split = "_"
#    )
#  )
#)
## vector of possible gain/loss multipliers, chosen randomly:
#CNV_multipliers <- as.numeric(
#  unlist(
#    strsplit(
#      "3_2_1.5_0.5_0",
#      split = "_"
#    )
#  )
#)
#downsample <- TRUE
#downsample_proportions <- as.numeric(
#  unlist(
#    strsplit(
##      "0.5",
#      "0.9_0.8_0.7_0.6_0.5_0.4_0.3_0.2_0.15_0.1_0.05",
#      split = "_"
#    )
#  )
#)
#simulation_number <- 1
#noise_cell_no <- 5000

print(paste0("Project name = ", project_name))
print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("Subset data? ", as.character(subset_data)))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(Seurat)
library(cowplot)
library(rlist, lib.loc = lib_loc)
library(GenomicRanges)
library(naturalsort, lib.loc = lib_loc)
library(splatter, lib.loc = lib_loc)
library(ggplot2)
library(org.Hs.eg.db)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", 
  project_name, "/", subproject_name, "/")
ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/functions/")
results_dir <- paste0(project_dir, "results/")
in_dir <- paste0(project_dir, "raw_files/seurat_objects/", 
  sample_name, "/")
annot_dir <- paste0(results_dir, "infercnv/", 
  sample_name, "/input_files/")
noise_dir <- paste0(results_dir, "cancer_simulation/", sample_name, 
  "_cancer_sim/noise_generation/")
system(paste0("mkdir -p ", noise_dir))

emptydrops_dir <- paste0(in_dir, "/emptydrops/")

sim_out_path <- paste0(results_dir, "cancer_simulation/", sample_name, 
	"_cancer_sim/")

common_Robject_dir <- paste0(sim_out_path, "Rdata/")
system(paste0("mkdir -p ", common_Robject_dir))
common_plot_dir <- paste0(sim_out_path, "plots/")
system(paste0("mkdir -p ", common_plot_dir))
common_table_dir <- paste0(sim_out_path, "tables/")
system(paste0("mkdir -p ", common_table_dir))

Robject_dir <- paste0(sim_out_path, "/", simulation_number, "/Rdata/")
system(paste0("mkdir -p ", Robject_dir))
plot_dir <- paste0(sim_out_path, "/", simulation_number, "/plots/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(sim_out_path, "/", simulation_number, "/tables/")
system(paste0("mkdir -p ", table_dir))

out_path <- paste0(results_dir, "infercnv/", sample_name, 
  "_cancer_sim/")
out_dir <- paste0(out_path, "/", simulation_number, "/")
no_downsample_dir <- paste0(out_dir, "no_downsampling/input_files/")
system(paste0("mkdir -p ", no_downsample_dir))

original_dir <- paste0(out_path, "original/input_files/")
system(paste0("mkdir -p ", original_dir))

print(paste0("Sample directory = ", in_dir))
print(paste0("Reference directory = ", ref_dir))
print(paste0("Output directory = ", out_dir))
print(paste0("Plot directory = ", plot_dir))
print(paste0("Table directory = ", table_dir))
print(paste0("Original directory = ", original_dir))

print(paste0("Generating simulated cancer data set from ", sample_name))
print(paste0("Filtering out cells with less than ", nUMI_threshold, " UMIs and ",
  nGene_threshold, " genes"))
print(paste0("Each sample will be downsampled to the following proportions:"))
paste(downsample_proportions, collapse = ", ")


################################################################################
### 0. Define functions and choose seeds ###
################################################################################

prepare_infercnv_metadata <- dget(paste0(func_dir, "prepare_infercnv_metadata.R"))
plot_expression <- dget(paste0(func_dir, "plot_expression.R"))
plot_median_FC <- dget(paste0(func_dir, "plot_median_FC.R"))
add_CNV <- dget(paste0(func_dir, "add_CNV.R"))

# choose and record seeds:
if (file.exists(paste0(table_dir, "random_seed_record.txt"))) {
  seed_record <- read.table(
    paste0(table_dir, "random_seed_record.txt"),
    sep = "\t",
    header = T,
    as.is = T
  )
} else {
  seed_record <- data.frame(
    row.names = c("CNV_no", "start_position", "length", "multiplier", "UMI_downsample",
      "gene_downsample"),
    seed = sample(1:999, 6)
  )
  write.table(
    seed_record, 
    paste0(table_dir, "random_seed_record.txt"),
    sep = "\t",
    row.names = T,
    col.names = T,
    quote = F
  )
}

set.seed(seed_record["CNV_no",])
CNV_no <- sample(seq(CNV_no_range[1], CNV_no_range[2]), 1)
print(paste0(
  "Adding ", CNV_no, " CNVs ",
  "with lengths between ", CNV_lengths[1], " and ", CNV_lengths[length(CNV_lengths)]
))


###################################################################################
### 1. Fetch counts matrix from normal sample, create original infercnv input files 
# and plot nUMI and nGene ###
###################################################################################

if (!file.exists(paste0(common_Robject_dir, "/1.original_epithelial_df.Rdata"))) {
  
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

  saveRDS(infercnv_metadata, paste0(in_dir, "metadata.Rdata"))
  saveRDS(seurat_10X, paste0(in_dir, "04_seurat_object_annotated.Rdata"))
  
  # only keep cells in metadata df:
  print(paste0("No cells in count df before filtering for those in metadata df = ", 
      ncol(count_df)))
  count_df <- count_df[,colnames(count_df) %in% infercnv_metadata$metadata$cell_ids]
  print(paste0("No cells in count df after filtering for those in metadata df = ", 
      ncol(count_df)))

  # save original infercnv input files:
  if (!file.exists(paste0(original_dir, "input_matrix.txt"))) {
    write.table(count_df, paste0(original_dir, "input_matrix.txt"), quote=F,
      sep="\t", col.names=T, row.names=T)
    write.table(infercnv_metadata$metadata, paste0(original_dir, "metadata.txt"), sep = "\t",
      quote = F, col.names = F, row.names = F)
  }

  # create density plots of nUMI and nGene:
  QC <- data.frame(
    row.names = colnames(count_df),
    nUMI = apply(count_df, 2, sum),
    nGene = apply(count_df, 2, function(x) length(x[x!=0]))
  )
  QC <- QC[colnames(count_df),]
  nUMI_density_plot <- density(QC$nUMI)
  pdf(paste0(common_plot_dir, "nUMI_density_plot.pdf"))
    plot(nUMI_density_plot, main=NA, xlab = "nUMI")
  dev.off()
  png(paste0(common_plot_dir, "nUMI_density_plot.png"))
    plot(nUMI_density_plot, main=NA, xlab = "nUMI")
  dev.off()
  nGene_density_plot <- density(QC$nGene)
  pdf(paste0(common_plot_dir, "nGene_density_plot.pdf"))
    plot(nGene_density_plot, main=NA, xlab = "nGene")
  dev.off()
  png(paste0(common_plot_dir, "nGene_density_plot.png"))
    plot(nGene_density_plot, main=NA, xlab = "nGene")
  dev.off()
  log_nUMI_density_plot <- density(log10(QC$nUMI))
  pdf(paste0(common_plot_dir, "log10_nUMI_density_plot.pdf"))
    plot(log_nUMI_density_plot, main=NA, xlab = "log10 nUMI")
  dev.off()
  png(paste0(common_plot_dir, "log10_nUMI_density_plot.png"))
    plot(log_nUMI_density_plot, main=NA, xlab = "log10 nUMI")
  dev.off()
  log_nGene_density_plot <- density(log10(QC$nGene))
  pdf(paste0(common_plot_dir, "log10_nGene_density_plot.pdf"))
    plot(log_nGene_density_plot, main=NA, xlab = "log10 nGene")
  dev.off()
  png(paste0(common_plot_dir, "log10_nGene_density_plot.png"))
    plot(log_nGene_density_plot, main=NA, xlab = "log10 nGene")
  dev.off()
  
  # filter out cells with nUMI < nUMI_threshold and nGene < nGene_threshold
  print(paste0("Number of cells before filtering out low coverage: ",
    ncol(count_df)))
  cells_to_keep <- rownames(QC)[QC$nUMI > nUMI_threshold & QC$nGene > nGene_threshold]
  count_df <- count_df[
    ,colnames(count_df) %in% cells_to_keep
  ]
  print(paste0("Number of cells after filtering out low coverage: ",
    ncol(count_df)))

  # remove outfiltered cells from metadata:
  infercnv_metadata$metadata <- infercnv_metadata$metadata[
    cells_to_keep,
  ]
  epithelial_df <- count_df[
    ,as.character(
      infercnv_metadata$metadata$cell_ids[
        grep("pithelial", infercnv_metadata$metadata$cell_type)
      ]
    )
  ]

  if (!file.exists(paste0(common_plot_dir, "nUMI_density_plot.pdf"))) {
    
    p <- ggplot(QC, aes(x=nUMI, y=nGene))
    p <- p + geom_point()
    p <- p + xlab("nUMI")
    p <- p + ylab("nGene")
    p <- p + theme(legend.title = element_blank())
    pdf(paste0(common_plot_dir, "QC_quad_plot.pdf"), 
      width = 10, height = 6)
      print(p)
    dev.off()
    png(paste0(common_plot_dir, "QC_quad_plot.png"), 
      width = 450, height = 270)
      print(p)
    dev.off()
    # create total count density quad plot:
    print("Determining total count...")
    total_counts <- apply(epithelial_df, 2, sum)
    total_count_density_plot <- density(total_counts, bw="SJ")
    pdf(paste0(common_plot_dir, "total_count_density_plot.pdf"))
      plot(total_count_density_plot, main=NA, xlab = "Total counts")
    dev.off()
    png(paste0(common_plot_dir, "total_count_density_plot.png"))
      plot(total_count_density_plot, main=NA, xlab = "Total counts")
    dev.off()
    # create log10 total count density quad plot:
    log_total_count_density_plot <- density(log10(total_counts), bw="SJ")
    pdf(paste0(common_plot_dir, "log10_total_count_density_plot.pdf"))
      plot(log_total_count_density_plot, main=NA, xlab = "Total counts")
    dev.off()
    png(paste0(common_plot_dir, "log10_total_count_density_plot.png"))
      plot(log_total_count_density_plot, main=NA, xlab = "Total counts")
    dev.off()
  }

  saveRDS(epithelial_df, paste0(common_Robject_dir, "/1.original_epithelial_df.Rdata"))

} else {
  epithelial_df <- readRDS(paste0(common_Robject_dir, "/1.original_epithelial_df.Rdata"))
}


################################################################################
### 2. Format annotation and counts matrices ###
################################################################################

if (!file.exists(paste0(common_Robject_dir, 
  "/1b.epithelial_df.Rdata")) | !file.exists(paste0(common_Robject_dir, 
  "/1c.gene_annotation.Rdata"))) {

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
  epithelial_df <- epithelial_df[
    gene_annotation$symbol[!(gene_annotation$chromosome %in% c("chrM", "chrX", "chrY"))],
  ]

  print(paste0("Gene numbers of InferCNV gene annotation and counts matrix are now ",
    nrow(gene_annotation), " and ", nrow(epithelial_df), " respectively"))
  
  # number genes in gene annotation:
  gene_annotation$number <- seq(1, nrow(gene_annotation))

  saveRDS(epithelial_df, paste0(common_Robject_dir, "/1b.epithelial_df.Rdata"))
  saveRDS(gene_annotation, paste0(common_Robject_dir, 
    "/1c.gene_annotation.Rdata"))

} else {

  epithelial_df <- readRDS(paste0(common_Robject_dir, "/1b.epithelial_df.Rdata"))
  gene_annotation <- readRDS(paste0(common_Robject_dir, 
    "/1c.gene_annotation.Rdata"))

}


################################################################################
### 3. Prepare chromosome information ###
################################################################################

# determine chromosome information:
chromosome_lengths <- lapply(split(gene_annotation, gene_annotation$chromosome), 
  nrow)
chromosome_lengths <- chromosome_lengths[
  naturalsort(names(chromosome_lengths))
]

for (n in 1:length(chromosome_lengths)) {
  print(n)
  if (n==1) {
    chromosome_coords <- list(1:chromosome_lengths[[n]])
  } else {
    # define starting coordinate as last coordinate of previous chromosome + 1:
    start_coord <- chromosome_coords[[n-1]][length(chromosome_coords[[n-1]])]+1
    chromosome_coords[[n]] <- start_coord:(start_coord+chromosome_lengths[[n]]-1)
  }
}
names(chromosome_coords) <- names(chromosome_lengths)

# determine chromosome ends and midpoints:
chromosome_ends <- lapply(chromosome_coords, max)
chromosome_midpoints <- lapply(chromosome_coords, function(x) x[floor(length(x)/2)])


################################################################################
### 4. Plot average fold difference from median ###
################################################################################

# generate average original counts vector with each value representing a gene:
average_original_counts <- apply(epithelial_df, 1, mean)
# add 0.1 to zero values:
average_original_counts <- average_original_counts + 3e-4
# determine median:
median_average_original_counts <- median(average_original_counts)

# divide by median to get fold change from median and add to df for plotting:
original_fold_change <- average_original_counts/median_average_original_counts
original_fold_change_df <- data.frame(
  number = seq(1, length(original_fold_change)),
  count = original_fold_change
)

# take mean of original fold change:
median_original_fold_change <- median(original_fold_change)

# take the log10:
log_original_fold_change <- log10(original_fold_change)
# tabulate for plotting:
log_original_fold_change_df <- data.frame(
  number = seq(1, length(log_original_fold_change)),
  count = log_original_fold_change
)

# calculate the median of the counts vector:
log_median_original_fold_change <- log10(median_original_fold_change)

# plot counts:
if (!file.exists(
  paste0(common_plot_dir, "1.log10_original_fold_change_from_median.png")
)) {
  
  p <- plot_expression(
    log_original_fold_change_df,
    chromosome_midpoints,
    chromosome_ends,
    y_lims = c(-4, 2),
    CNV_record = "none"
  )

  png(
    paste0(common_plot_dir, "1.log10_original_fold_change_from_median.png"), 
    width = 20,
    height = 10,
    res = 300,
    units = "in"
  )
    print(p)
  dev.off()

}

#save.image(paste0(Robject_dir, "temp.Rdata"))
#rm(list=ls())
#simulation_number <- 1
#Robject_dir <- paste0("/share/ScratchGeneral/jamtor/projects/thesis/Figure_1.8_accuracy_vs_coverage/results/cancer_simulation/CID4520N_cancer_sim/",
#  simulation_number, "/Rdata/")
#load(paste0(Robject_dir, "temp.Rdata"))
#add_CNV <- dget(paste0(func_dir, "add_CNV.R"))


################################################################################
### 5. Add CNVs to df, determine fold change from original median of each 
# gene and median fold change for each CNV region ###
################################################################################

######
#rm(CNV_record)
#rm(modified_df)
#rm(log_modified_fold_change_df)
######

# determine total gene length of genome to ensure total CNV length < 0.4 x 
# genome length to comply with InferCNV assumptions:
genome_length <- nrow(epithelial_df)

if ( !file.exists(paste0(Robject_dir, "/2a.pre_noise_simulated_epithelial_dfs.Rdata")) | 
  !file.exists(paste0(Robject_dir, "/2b.pre_noise_log_modified_fold_change_dfs.Rdata")) |
  !file.exists(paste0(Robject_dir, "/2c.CNV_records.Rdata")) ) {

  print("Adding CNVs to dataset...")
  writeLines("\n")

  both_CNV_data <- add_CNV(
    epithelial_df,
    log_original_fold_change_df,
    CNV_no, 
    seed_record,
    subpops = "subpops",
    genome_length,
    chromosome_ends
  )
  modified_dfs <- both_CNV_data$modified_dfs
  log_modified_fold_change_dfs <- both_CNV_data$log_mod_FC_dfs
  CNV_records <- both_CNV_data$CNV_record
  
  saveRDS(modified_dfs, paste0(Robject_dir, "/2a.pre_noise_simulated_epithelial_dfs.Rdata"))
  saveRDS(
    log_modified_fold_change_dfs, 
    paste0(Robject_dir, "/2b.pre_noise_log_modified_fold_change_dfs.Rdata")
  )
  saveRDS(CNV_records, paste0(Robject_dir, "/2c.CNV_records.Rdata"))

} else {

  modified_dfs <- readRDS(paste0(Robject_dir, 
    "/2a.pre_noise_simulated_epithelial_dfs.Rdata"))
  log_modified_fold_change_dfs <- readRDS(paste0(Robject_dir, 
    "/2b.pre_noise_log_modified_fold_change_dfs.Rdata"))
  CNV_records <- readRDS(paste0(Robject_dir, "/2c.CNV_records.Rdata"))

}


################################################################################
### 6. Create modified counts and line CNV plots ###
################################################################################

for (i in 1:length(CNV_records)) {
  # plot median fold change from original median for modified data:
  if (!file.exists(paste0(
    plot_dir, 
    "2a.pre_noise_log_modified_fold_change_from_median_line_only_", 
    names(CNV_records)[i], 
    ".png"
  ))) {

    p <- suppressWarnings(
      plot_expression(
        log_modified_fold_change_dfs[[i]],
        chromosome_midpoints,
        chromosome_ends,
        y_lims = c(-4, 4), 
        CNV_record = CNV_records[[i]],
        line_only = TRUE
      )
    )

    png(
      paste0(
        plot_dir, 
        "2a.pre_noise_log_modified_fold_change_from_median_line_only_", 
        names(CNV_records)[i], 
        ".png"
      ), 
      width = 20,
      height = 10,
      res = 300,
      units = "in"
    )
      print(p)
    dev.off()
    
  }
  
  # plot counts:
  if (!file.exists(paste0(
    plot_dir, 
    "2b.pre_noise_log_modified_fold_change_from_median_",
    names(CNV_records)[i],
    ".png"
  ))) {
   
    p <- suppressWarnings(
      plot_expression(
        log_modified_fold_change_dfs[[i]],
        chromosome_midpoints,
        chromosome_ends,
        y_lims = c(-4, 4), 
        CNV_record = CNV_records[[i]],
        line_only = FALSE
      )
    )

    png(
      paste0(
    plot_dir, 
    "2b.pre_noise_log_modified_fold_change_from_median_",
    names(CNV_records)[i],
    ".png"
  ), 
      width = 20,
      height = 10,
      res = 300,
      units = "in"
    )
      print(p)
    dev.off()

  }
}
#save.image(paste0(Robject_dir, "temp2.Rdata"))
#load(paste0(Robject_dir, "temp2.Rdata"))


################################################################################
### 7. Bind new epithelial and original stromal counts and check against 
# metadata ###
################################################################################

if (!exists("count_df")) {
  # load seurat object:
  seurat_10X <- readRDS(paste0(in_dir, "04_seurat_object_annotated.Rdata"))
  # create raw matrix input file and subset if necessary:
  count_df <- as.matrix(GetAssayData(seurat_10X , slot = "counts"))
}

# only keep genes in both count and modified_df:
print(paste0(
  "Number of genes in count_df before removing genes not in modified_df = ", 
  nrow(count_df))
)
count_df <- count_df[rownames(modified_dfs[[1]]),]
print(paste0(
  "Dimensions of count_df after removing genes not in modified_df = ", 
  nrow(count_df))
)

if (!exists("infercnv_metadata")) {
  infercnv_metadata <- readRDS(paste0(in_dir, "metadata.Rdata"))
}

print(
  paste0(
    "Columns in stromal part of count_df = ", 
    ncol(
      count_df[
        ,as.character(infercnv_metadata$metadata$cell_ids[
          grep("pithelial", infercnv_metadata$metadata$cell_type, invert = T)
          ]
        )
      ]
    )
  )
)
print(paste0("Columns in modified_df = ", ncol(modified_dfs[[1]])))

# bind non-epithelial count_df and modified_df together as new_counts:
temp_non_epithelial_df <- count_df[
  ,as.character(infercnv_metadata$metadata$cell_ids[
    grep("pithelial", infercnv_metadata$metadata$cell_type, invert = T)
    ]
  )
]
new_counts <- lapply(modified_dfs, function(x) {

  nc <- cbind(x, temp_non_epithelial_df)
  print(paste0("Columns in new_counts = ", ncol(nc)))

  # only keep cells in cell_annotation:
  nc <- nc[
    ,colnames(nc) %in% as.character(infercnv_metadata$metadata$cell_ids)
  ]

  print(paste0(
    "Columns in new_counts after removing cells not in cell_annotation = ", 
    ncol(nc))
  )

  return(nc)
})


###################################################################################
### 8. Simulate noise dataset from noise counts ###
###################################################################################

if (!file.exists(paste0(noise_dir, "/noise_df.Rdata"))) {

  all_cells <- read.csv(paste0(emptydrops_dir, "01_Emptydrops_out.csv"))
  filtered_cells <- read.csv(
    paste0(emptydrops_dir, "02_emptydrops_filtered_cell_ids.csv")
  )
  raw_10X <- readRDS(paste0(in_dir, "01_Read10X_raw_data.Rdata"))
  seurat_10X <- readRDS(paste0(in_dir, "03_seurat_object_processed.Rdata"))
  
  # isolate outfiltered cells from raw counts matrix:
  outfiltered_cells <- all_cells[!(all_cells$X %in% filtered_cells$cell_ids),]
  cells_to_select <- colnames(raw_10X)[
    which(colnames(raw_10X) %in% outfiltered_cells$X)
  ]
  outfiltered_counts <- raw_10X[,colnames(raw_10X) %in% cells_to_select]
  
  # sanity checks:
  print(paste0(
    "Total number of cells in EmptyDrops out object = ", nrow(all_cells)
  ))
  print(paste0(
    "Total number of cells in raw matrix object = ", ncol(raw_10X)
  ))
  print(paste0(
    "Number of cells in filtered EmptyDrops object = ", nrow(filtered_cells)
  ))
  print(paste0(
    "Number of cells in seurat object = ", ncol(GetAssayData(seurat_10X , slot = "counts"))
  ))
  print(paste0(
    "Does total cell number - number of outfiltered cells = number of filtered cells? ", 
    nrow(all_cells) - nrow(outfiltered_cells) == nrow(filtered_cells)
  ))
  print(paste0(
    "Number of cells in outfiltered_counts = ", ncol(outfiltered_counts)
  ))

  # order by total counts:
  column_sums <- Matrix::colSums(outfiltered_counts)
  outfiltered_counts <- outfiltered_counts[
    ,order(column_sums, decreasing = T)
  ]

  # plot distribution of outfiltered counts:
  count_sums <- colSums(outfiltered_counts)
  outfiltered_counts <- outfiltered_counts[,count_sums != 0]
  outfiltered_counts <- outfiltered_counts[,1:noise_cell_no]
  outfiltered_means <- apply(outfiltered_counts, 2, mean)
  if (!file.exists(paste0(noise_dir, "outfiltered_count_density_plot.pdf"))) {
    outfiltered_count_density_plot <- density(as.vector(outfiltered_means))
    pdf(paste0(noise_dir, "outfiltered_count_density_plot.pdf"))
      plot(outfiltered_count_density_plot, main=NA, xlab = "nUMI")
    dev.off()
  }

  # select for cells to use for simulated noise dataset:
  noise_input <- outfiltered_counts[,1:noise_cell_no]
 
  if (!file.exists(paste0(noise_dir, "/noise_params.Rdata"))) {

    # calculate params from input to use for simulation:
    splat_params <- splatEstimate(as.matrix(noise_input))
    # reduce simulation to number of cells in simulation:
    splat_params <- setParam(splat_params, "nGenes", nrow(new_counts[[1]]))

    saveRDS(splat_params, paste0(noise_dir, "/noise_params.Rdata"))

  } else {

    splat_params <- readRDS(paste0(noise_dir, "/noise_params.Rdata"))

  }

  # simulate the data:
  sim <- splatSimulate(
  	splat_params, 
  	batchCells = ncol(new_counts[[1]]) + ncol(new_counts[[2]])
  )
  noise_counts <- counts(sim)
 
  print("Table of simulated noise counts:")
  table(unlist(noise_counts))

  # find chromosome ends to mark on plot:
  # load co-ordinates of genes/chromosomes:
  library(TxDb.Hsapiens.UCSC.hg38.knownGene, lib.loc=lib_loc)
  gene_coords <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
  gene_coords <- gene_coords[grep("chr[0-9]", seqnames(gene_coords))]
  gene_coords <- gene_coords[grep("_.*$", seqnames(gene_coords), invert=T)]

  # convert gene number ids to symbols:
  egSYMBOL <- toTable(org.Hs.egSYMBOL)
  m <- match(gene_coords$gene_id, egSYMBOL$gene_id)
  gene_coords$symbol <- egSYMBOL$symbol[m]

  # isolate last end position/gene symbol of each chromosome: 
  gene_chr <- data.frame(
  	as.character(seqnames(gene_coords)),
  	end(ranges(gene_coords)),
  	as.character(gene_coords$symbol)
  )
  colnames(gene_chr) <- c("seqnames", "ends", "symbol")

  # keep only genes in noise_input:
  input_gene_chr <- gene_chr[gene_chr$symbol %in% rownames(noise_input),]
  input_gene_chr <- input_gene_chr[order(input_gene_chr$end),]
  input_gene_chr <- input_gene_chr[order(input_gene_chr$seqnames),]
  input_gene_chr$symbol <- as.character(input_gene_chr$symbol)

  # find genes at end of each chromosome to mark on count plots:
  input_gene_chr_split <- split(input_gene_chr, input_gene_chr$seqnames)
  input_chr_ends <- lapply(input_gene_chr_split, function(x) return(x$symbol[nrow(x)]))
  input_end_genes <- unlist(input_chr_ends)
  input_end_genes <- input_end_genes[naturalsort(names(input_end_genes))]

  input_chr_mids <- lapply(input_gene_chr_split, function(x) return(x$symbol[floor(nrow(x)/2)]))
  input_mid_genes <- unlist(input_chr_mids)
  input_mid_genes <- input_mid_genes[naturalsort(names(input_mid_genes))]
  input_mid_ind <- match(input_mid_genes, rownames(outfiltered_counts))

  input_means <- apply(noise_input, 1, mean)
  input_means <- input_means[input_means < 0.9]
  input_df <- data.frame(
    number = 1:length(input_means),
    input = input_means
  )

  # log10 the count values:
  input_df$input <- log10(input_df$input+1e-4)
  
  p1 <- ggplot(input_df, aes(x=number, y=input))
  p1 <- p1 + geom_point()
  p1 <- p1 + xlab("Genomic location")
  p1 <- p1 + ylab("log10 mean count")
  for (g in 1:length(input_end_genes)) {
    #print(input_end_genes[g])

    input_df_rownames <- gsub("-|_", "", rownames(input_df))
    input_end_gene <- gsub("-|_", "", input_end_genes[g])
    ind <- grep(paste0("\\b", input_end_gene, "\\b"), input_df_rownames)
    
    # create chromosome end vertical line:
    p1 <- p1 + geom_segment(
      x=ind,
      xend=ind,
      y=min(input_df$input), 
      yend=max(input_df$input), 
      size=0.3, color="red"
    )

  }

  # create chromosome labels:
  p1 <- p1 + scale_x_continuous(
    breaks = input_mid_ind,
    labels = 1:length(input_mid_ind)
  )

  pdf(paste0(noise_dir, "log10_noise_input_scatterplot.pdf"))
    p1
  dev.off()

  # keep only genes in noise_counts:
  rownames(noise_counts) <- rownames(new_counts[[1]])
  sim_gene_chr <- gene_chr[gene_chr$symbol %in% rownames(noise_counts),]
  sim_gene_chr <- sim_gene_chr[order(sim_gene_chr$end),]
  sim_gene_chr <- sim_gene_chr[order(sim_gene_chr$seqnames),]
  sim_gene_chr$symbol <- as.character(sim_gene_chr$symbol)

  # find genes at end of each chromosome to mark on count plots:
  sim_gene_chr_split <- split(sim_gene_chr, sim_gene_chr$seqnames)
  sim_chr_ends <- lapply(sim_gene_chr_split, function(x) return(x$symbol[nrow(x)]))
  sim_end_genes <- unlist(sim_chr_ends)
  sim_end_genes <- sim_end_genes[naturalsort(names(sim_end_genes))]

  sim_chr_mids <- lapply(sim_gene_chr_split, function(x) return(x$symbol[floor(nrow(x)/2)]))
  sim_mid_genes <- unlist(sim_chr_mids)
  sim_mid_genes <- sim_mid_genes[naturalsort(names(sim_mid_genes))]

  sim_means <- apply(noise_counts, 1, mean)
  #sim_means <- sim_means[sim_means < 0.9]
  sim_df <- data.frame(
    number = 1:length(sim_means),
    sim = sim_means
  )

  sim_mid_ind <- match(sim_mid_genes, rownames(sim_df))

  # log10 the count values:
  sim_df$sim <- log10(sim_df$sim+1e-4)

  p2 <- ggplot(sim_df, aes(x=number, y=sim))
  p2 <- p2 + geom_point()
  p2 <- p2 + xlab("Genomic location")
  p2 <- p2 + ylab("log10 mean count")
  for (g in 1:length(sim_end_genes)) {
    #print(sim_end_genes[g])

    sim_df_rownames <- gsub("-|_", "", rownames(sim_df))
    sim_end_gene <- gsub("-|_", "", sim_end_genes[g])
    ind <- grep(paste0("\\b", sim_end_gene, "\\b"), sim_df_rownames)
#    names(ind) <- rownames(outfiltered_counts)[ind]
#    ind <- ind[grep("-|_", names(ind), invert=T)]
    
    # create chromosome end vertical line:
    p2 <- p2 + geom_segment(
      x=ind,
      xend=ind,
      y=min(sim_df$sim), 
      yend=max(sim_df$sim), 
      size=0.3, color="red"
    )
  }
  # create chromosome labels:
  p2 <- p2 + scale_x_continuous(
    breaks = sim_mid_ind,
    labels = 1:length(sim_mid_ind)
  )

  pdf(paste0(noise_dir, "log10_noise_sim_scatterplot.pdf"))
    p2
  dev.off()
  
  saveRDS(noise_counts, paste0(noise_dir, "/noise_df.Rdata"))

} else {
  noise_counts <- readRDS(paste0(noise_dir, "/noise_df.Rdata"))
}

#save.image(paste0(Robject_dir, "temp3.Rdata"))
#load(paste0(Robject_dir, "temp3.Rdata"))


###################################################################################
### 9. Add noise to simulated counts ###
###################################################################################

# ensure metadata has same cells as modified_df and label stromal cells:
epithelial_metadata <- infercnv_metadata$metadata[
  grep("pithelial", infercnv_metadata$metadata$cell_type),
]
epithelial_metadata <- epithelial_metadata[
  c(
    colnames(modified_dfs[[1]]), 
    colnames(modified_dfs[[2]])
  ),
]

epithelial_metadata$cell_type <- "Epithelial"
non_epithelial_metadata <- infercnv_metadata$metadata[
  grep("pithelial", infercnv_metadata$metadata$cell_type, invert = T),
]
non_epithelial_metadata$cell_type <- "Non_epithelial"
metadata_df <- rbind(
  non_epithelial_metadata,
  epithelial_metadata
)

# cbind together both subpopulations:
nondownsampled_counts_without_noise <- do.call("cbind", new_counts)

print("Summary of new counts before noise added:")
print(summary(as.vector(nondownsampled_counts_without_noise)))
nondownsampled_counts_with_noise <- nondownsampled_counts_without_noise+noise_counts
print("Summary of new counts after noise added:") 
print(summary(as.vector(nondownsampled_counts_with_noise)))

# create updated CNV scatterplots post-noise addition:
# generate mean original segment fold change vector with each value representing 
# a gene:
final_epithelial_df <- nondownsampled_counts_with_noise[
  ,colnames(nondownsampled_counts_with_noise) %in% epithelial_metadata$cell_ids[
    epithelial_metadata$cell_type == "Epithelial"
  ]
]

# split final epithelial and original epithelial dfs according to new counts:
orig_epithelial_dfs <- list(
  sub1 = epithelial_df[, colnames(epithelial_df) %in% colnames(new_counts$sub1)],
  sub2 = epithelial_df[, colnames(epithelial_df) %in% colnames(new_counts$sub2)]
)

final_epithelial_dfs <- list(
  sub1 = final_epithelial_df[
    , colnames(final_epithelial_df) %in% colnames(new_counts$sub1)
  ],
  sub2 = final_epithelial_df[
    , colnames(final_epithelial_df) %in% colnames(new_counts$sub2)
  ]
)

for (k in 1:length(final_epithelial_dfs)) {
  
  plot_data <- plot_median_FC(
    orig_epithelial_dfs[[k]],
    final_epithelial_dfs[[k]],
    log_original_fold_change_df,
    CNV_records[[k]]
  )

  if (k==1) {
  	final_plot_data <- list(sub1 = plot_data)
  } else {
  	final_plot_data$sub2 <- plot_data
  }

}


###################################################################################
### 10. Plot post-noise gene expression profiles ###
###################################################################################

for (i in 1:length(final_plot_data)) {

#  # plot median fold change from original median for modified data:
#  if (!file.exists(paste0(plot_dir, 
#        "3a.log_modified_fold_change_from_median_with_noise_line_only_",
#        names(final_plot_data)[i], ".png"))) {
#    p <- suppressWarnings(
#      plot_expression(
#        final_plot_data[[i]]$log_mod_FC_df,
#        chromosome_midpoints,
#        chromosome_ends,
#        y_lims = c(-4, 4), 
#        CNV_record = final_plot_data[[i]]$CNV_record,
#        line_only = TRUE
#      )
#    )
#  
#    png(
#      paste0(plot_dir, 
#        "3a.log_modified_fold_change_from_median_with_noise_line_only_",
#        names(final_plot_data)[i], ".png"), 
#      width = 20,
#      height = 10,
#      res = 300,
#      units = "in"
#    )
#      print(p)
#    dev.off()
#  }
  if (!file.exists(paste0(plot_dir, "3b.log_modified_fold_change_from_median_with_noise_",
    names(final_plot_data)[i], ".png"))) {
   
    p <- suppressWarnings(
        plot_expression(
          final_plot_data[[i]]$log_mod_FC_df,
          chromosome_midpoints,
          chromosome_ends,
          y_lims = c(-4, 4), 
          CNV_record = final_plot_data[[i]]$CNV_record,
          line_only = FALSE
        )
      )
    png(
      paste0(plot_dir, "3b.log_modified_fold_change_from_median_with_noise_",
        names(final_plot_data)[i], ".png"), 
      width = 20,
      height = 10,
      res = 300,
      units = "in")
      print(p)
    dev.off()
  }

}

# fetch indices of differentiating CNV:
diff_CNV <- CNV_records[[1]][
  CNV_records[[1]]$differentiating,
]$start:CNV_records[[1]][
  CNV_records[[1]]$differentiating,
]$end
# summary of differentiating CNV before adding:
print("Summary of differentiating CNV range in subclone 1 before adding:")
summary(unlist(as.data.frame(orig_epithelial_dfs[[1]][diff_CNV,])))
print("Summary of differentiating CNV range in subclone 2 before adding:")
summary(unlist(as.data.frame(orig_epithelial_dfs[[2]][diff_CNV,])))

# summary of differentiating CNV after adding:
print("Summary of differentiating CNV range in subclone 1 after adding:")
summary(unlist(as.data.frame(final_epithelial_dfs[[1]][diff_CNV,])))
print("Summary of differentiating CNV range in subclone 2 after adding:")
summary(unlist(as.data.frame(final_epithelial_dfs[[2]][diff_CNV,])))

# save cell ids of subpops:
subpop_ids <- lapply(modified_dfs, colnames)
saveRDS(subpop_ids, paste0(Robject_dir, "subpop_ids.Rdata"))

# save other data:
final_CNV_records <- list(
  final_plot_data[[1]]$CNV_record,
  final_plot_data[[2]]$CNV_record
)
final_log_mod_FC_dfs <- list(
  final_plot_data[[1]]$log_mod_FC_df,
  final_plot_data[[2]]$log_mod_FC_df
)

saveRDS(final_CNV_records, paste0(Robject_dir, "CNV_records.Rdata"))
saveRDS(
  final_log_mod_FC_dfs, 
  paste0(Robject_dir, "log_modified_fold_change_dfs.Rdata")
)

#######
#summary(unlist(as.data.frame(
#  nondownsampled_counts_with_noise[
#    ,colnames(nondownsampled_counts_with_noise) %in% colnames(final_epithelial_dfs[[1]])
#  ][diff_CNV,]
#)))
#
#summary(unlist(as.data.frame(
#  nondownsampled_counts_with_noise[
#    ,colnames(nondownsampled_counts_with_noise) %in% colnames(final_epithelial_dfs[[2]])
#  ][diff_CNV,]
#)))
#######

# save non-downsampled counts:
if (!file.exists(paste0(no_downsample_dir, "input_matrix.txt"))) {
  write.table(
  	nondownsampled_counts_with_noise, 
  	paste0(no_downsample_dir, "input_matrix.txt"), 
  	quote=F,
    sep="\t", 
    col.names=T, 
    row.names=T
  )
  write.table(
  	metadata_df, 
  	paste0(no_downsample_dir, "metadata.txt"), 
  	sep = "\t",
    quote = F, 
    col.names = F, 
    row.names = F
  )
}

# determine nUMI and nGene per cell:
mean_coverage <- data.frame(
  downsample_proportion = "no",
  nUMI = round(mean(apply(nondownsampled_counts_with_noise, 2, sum)), 0),
  nGene = round(
    mean(
      apply(
        nondownsampled_counts_with_noise, 2, function(x) length(which(x > 0))
      )
    ), 0
  )
)


################################################################################
### 11. Downsample new counts and save ###
################################################################################

if (downsample) {

  # downsample by UMI:
  library(DropletUtils, lib.loc = lib_loc)

  print("Total counts before downsampling:")
  print(sum(as.vector(nondownsampled_counts_with_noise)))

  # choose seeds for downsampling:
  set.seed(seed_record["UMI_downsample",])
  random_UMI_seeds <- sample(1:999, length(downsample_proportions))

  m=1
  for (proportion in downsample_proportions) {

    # downsample simulated dataset with noise added:
    print(paste0("Downsampling counts to ", proportion, " total UMIs..."))

    set.seed(random_UMI_seeds[m])
    downsampled_counts <- downsampleMatrix(
      nondownsampled_counts_with_noise, proportion, bycol=T
    )
    print("Total counts after downsampling:")
    print(sum(as.vector(downsampled_counts)))

    # determine nUMI and nGene per cell:
    mean_coverage <- rbind(
      mean_coverage,
      data.frame(
        downsample_proportion = as.character(proportion),
        nUMI = round(mean(apply(downsampled_counts, 2, sum)), 0),
        nGene = round(
          mean(
            apply(
              downsampled_counts, 2, function(x) length(which(x > 0))
            )
          ), 0
        )
      )
    )

    downsample_dir <- paste0(out_dir, proportion, 
      "_downsampling/input_files/")
    system(paste0("mkdir -p ", downsample_dir))
    write.table(
      downsampled_counts, paste0(downsample_dir, "input_matrix.txt"), 
      quote=F, sep="\t", col.names=T, row.names=T)
    write.table(metadata_df, paste0(downsample_dir, "metadata.txt"), sep = "\t",
      quote = F, col.names = F, row.names = F)

    m <<- m+1
  
  }

  # choose seeds for downsampling:
  set.seed(seed_record["gene_downsample",])
  random_gene_seeds <- sample(1:999, length(downsample_proportions))

  s=1
  for (proportion in downsample_proportions) {

    original_gene_no <- nrow(nondownsampled_counts_with_noise)
    downsampled_gene_no <- original_gene_no*proportion
    to_remove_no <- original_gene_no - downsampled_gene_no

    print("Total genes before downsampling:")
    print(original_gene_no)

    # downsample simulated dataset with noise added:
    print(paste0("Downsampling counts to ", proportion, " total genes..."))

    set.seed(random_gene_seeds[s])
    gene_downsampled_counts <- nondownsampled_counts_with_noise[
      -sample(1:original_gene_no, to_remove_no),
    ]

    print("Total genes after downsampling:")
    print(nrow(gene_downsampled_counts))

    # determine nUMI and nGene per cell:
    mean_coverage <- rbind(
      mean_coverage,
      data.frame(
        downsample_proportion = paste0(proportion, "_gene"),
        nUMI = round(mean(apply(gene_downsampled_counts, 2, sum)), 0),
        nGene = round(
          mean(
            apply(
              gene_downsampled_counts, 2, function(x) length(which(x > 0))
            )
          ), 0
        )
      )
    )

    gene_downsample_dir <- paste0(out_dir, proportion, 
      "_gene_downsampling/input_files/")
    system(paste0("mkdir -p ", gene_downsample_dir))
    write.table(
      gene_downsampled_counts, paste0(gene_downsample_dir, "input_matrix.txt"), 
      quote=F, sep="\t", col.names=T, row.names=T)
    write.table(metadata_df, paste0(gene_downsample_dir, "metadata.txt"), sep = "\t",
      quote = F, col.names = F, row.names = F)

    s <<- s+1
  
  }

}

# save mean coverage values:
write.table(
  mean_coverage, 
  paste0(table_dir, "/mean_coverage_values.txt"),
  sep = "\t",
  quote = F,
  col.names = T
)

