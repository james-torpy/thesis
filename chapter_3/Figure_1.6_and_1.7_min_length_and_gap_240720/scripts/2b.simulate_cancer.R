#! /share/ClusterShare/software/contrib/briglo/octR/src/R-3.6.0/builddir/bin/Rscript


###############################################################################
### Simulates cancer datasets taking as inputs number of CNVs, CNV lengths, ###
### gain/loss multipliers and downsamples required                          ###
###############################################################################

args = commandArgs(trailingOnly=TRUE)

project_name <- "thesis"
subproject_name <- "Figure_2.4_min_length_and_gap"
sample_name <- args[1]
print(paste0("sample name = ", sample_name))
nUMI_threshold <- as.numeric(args[2])
print(paste0("nUMI_threshold = ", nUMI_threshold))
nGene_threshold <- as.numeric(args[3])
print(paste0("nGene_threshold = ", nGene_threshold))
gap_or_CNV <- args[4]
print(paste0("Gap or CNV = ", gap_or_CNV))
CNV_type <- args[5]
print(paste0("CNV_type = ", CNV_type))
# range of lengths of CNVs/gaps:
feature_lengths <- as.numeric(
  unlist(
    strsplit(
      args[6],
      split = "_"
    )
  )
)
print(paste0("Possible CNV/gap lengths = ", feature_lengths))
gap_CNV_length <- as.numeric(args[7])
simulation_number <- as.numeric(args[8])
print(paste0("Simulation number = ", simulation_number))
noise_cell_no <- as.numeric(args[9])
print(paste0("Noise input cell number = ", noise_cell_no))
t_cells_included <- as.logical(args[10])
print(paste0("T-cells included = ", t_cells_included))
analysis_mode <- args[11]
print(paste0("Analysis mode = ", analysis_mode))

project_name <- "thesis"
subproject_name <- "Figure_2.4_min_length_and_gap"
sample_name <- "CID4520N"
nUMI_threshold <- 25000
nGene_threshold <- 5000
gap_or_CNV <- "CNV"
CNV_type <- "both"
# range of lengths of CNVs/gaps:
feature_lengths <- as.numeric(
  unlist(
    strsplit(
      "400_300_200_150_125_100_75_50_40_30_20_15_10_5",
      split = "_"
    )
  )
)
gap_CNV_length <- 100
simulation_number <- 1
noise_cell_no <- 5000
t_cells_included <- TRUE
analysis_mode <- "samples"

RStudio <- FALSE

print(paste0("Project name = ", project_name))
print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(Seurat)
library(cowplot)
library(rlist, lib.loc = lib_loc)
library(GenomicRanges)
library(naturalsort, lib.loc = lib_loc)
library(splatter, lib.loc = lib_loc)
library(ggplot2)
library(org.Hs.eg.db)
library(grid)

if (RStudio) {
  home_dir <- "/Users/jamestorpy/clusterHome/"
} else {
  home_dir <- "/share/ScratchGeneral/jamtor/"
}
project_dir <- paste0(home_dir, "projects/", 
  project_name, "/", subproject_name, "/")
ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/functions/")
results_dir <- paste0(project_dir, "results/")

if (t_cells_included) {
  in_path <- paste0(results_dir, "infercnv/t_cells_included/", sample_name, "/")
  in_dir <- paste0(in_path, "/input_files/")
  normal_infercnv_dir <- paste0(in_path, analysis_mode, "_mode/")
  in_Robject_dir <- paste0(results_dir, "infercnv/t_cells_included/", sample_name,
    "/Rdata/")
} else {
  in_path <- paste0(results_dir, "infercnv/t_cells_excluded/", sample_name, "/")
  in_dir <- paste0(in_path, "/input_files/")
  normal_infercnv_dir <- paste0(in_path, analysis_mode, "_mode/")
  in_Robject_dir <- paste0(in_path, "/Rdata/")
}

seurat_dir <- paste0(project_dir, "raw_files/seurat_objects/", 
  sample_name, "/")
emptydrops_dir <- paste0(seurat_dir, "/emptydrops/")

noise_dir <- paste0(results_dir, "cancer_simulation/", sample_name, 
  "_cancer_sim/noise_generation/")
system(paste0("mkdir -p ", noise_dir))

sim_out_path <- paste0(results_dir, "cancer_simulation/", sample_name, 
	"_cancer_sim/")

common_Robject_dir <- paste0(sim_out_path, "Rdata/")
system(paste0("mkdir -p ", common_Robject_dir))
common_plot_dir <- paste0(sim_out_path, "plots/")
system(paste0("mkdir -p ", common_plot_dir))
common_table_dir <- paste0(sim_out_path, "tables/")
system(paste0("mkdir -p ", common_table_dir))

if (gap_or_CNV == "CNV") {

  Robject_dir <- paste0(sim_out_path, gap_or_CNV, "/",
   simulation_number, "/Rdata/")
  system(paste0("mkdir -p ", Robject_dir))
  plot_dir <- paste0(sim_out_path, gap_or_CNV, "/", 
    simulation_number, "/plots/")
  system(paste0("mkdir -p ", plot_dir))
  table_dir <- paste0(sim_out_path, gap_or_CNV, "/", 
    simulation_number, "/tables/")
  system(paste0("mkdir -p ", table_dir))

} else if (gap_or_CNV == "gap") {

  Robject_dir <- paste0(sim_out_path, gap_or_CNV, "/",
   simulation_number, "/", CNV_type, "/Rdata/")
  system(paste0("mkdir -p ", Robject_dir))
  plot_dir <- paste0(sim_out_path, gap_or_CNV, "/", 
    simulation_number, "/", CNV_type, "/plots/")
  system(paste0("mkdir -p ", plot_dir))
  table_dir <- paste0(sim_out_path, gap_or_CNV, "/", 
    simulation_number, "/", CNV_type, "/tables/")
  system(paste0("mkdir -p ", table_dir))

}

out_path <- paste0(results_dir, "infercnv/t_cells_included/", sample_name, 
  "_cancer_sim/")

if (gap_or_CNV == "CNV") {
  out_dir <- paste0(out_path, gap_or_CNV, "/", simulation_number, 
    "/input_files/")
    system(paste0("mkdir -p ", out_dir))
  print(paste0("Output directory = ", out_dir))
} else if (gap_or_CNV == "gap") {
  out_dir <- paste0(out_path, gap_or_CNV, "/", simulation_number, 
    "/", CNV_type, "/input_files/")
  system(paste0("mkdir -p ", out_dir))
  print(paste0("Output directory = ", out_dir))
}

print(paste0("Sample directory = ", in_dir))
print(paste0("Reference directory = ", ref_dir))

print(paste0("Plot dir = ", plot_dir))
print(paste0("Table directory = ", table_dir))

print(paste0("Generating simulated cancer data set from ", sample_name))
print(paste0("Filtering out cells with less than ", nUMI_threshold, " UMIs and ",
  nGene_threshold, " genes"))


################################################################################
### 0. Define functions and choose seeds ###
################################################################################

introduce_CNVs <- dget(paste0(func_dir, "introduce_CNVs.R"))
introduce_gaps <- dget(paste0(func_dir, "introduce_gaps.R"))

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
    row.names = c("gain_start", "loss_start", "gain_gap_start", "loss_gap_start"),
    seed = sample(1:999, 4)
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


#################################################################################
#### 1. Filter genes filtered out in normal datset infercnv run from input 
# matrix and create new gene_annotation df ###
#################################################################################

if (
  !file.exists(paste0(common_Robject_dir, "1a.gene_annotation.Rdata")) | 
  !file.exists(paste0(common_Robject_dir, "1b.original_epithelial_df.Rdata")) |
  !file.exists(paste0(common_Robject_dir, "1c.original_metadata.Rdata"))
) {

  # load unfiltered matrix from last run:
  raw_matrix <- read.table(paste0(in_dir, "input_matrix.txt"), header = T, 
    sep = "\t")
  
  # load metadata:
  infercnv_metadata <- read.table(paste0(in_dir, "metadata.txt"), header = F, 
    sep = "\t")
  colnames(infercnv_metadata) <- c("cell_ids", "cell_type")
  rownames(infercnv_metadata) <- infercnv_metadata$cell_ids
  
  # isolate epithelial cells only:
  epithelial_ids <- as.character(
    infercnv_metadata$cell_ids[
      infercnv_metadata$cell_type == "Epithelial"
    ]
  )
  epithelial_matrix <- raw_matrix[
    ,colnames(raw_matrix) %in% epithelial_ids
  ]
  # fetch genes remaining from filtered infercnv output:
  last_run_genes <- gsub(
    "\"",
    "",
    system(
      paste0("awk '{print $1}' ", normal_infercnv_dir, "infercnv.12_denoised.observations.txt"),
      intern=TRUE
    )
  )[-1]
  # filter input matrix for only genes retained in last run infercnv output and write:
  epithelial_df <- epithelial_matrix[last_run_genes,]
  
#  # subset genes and cells:
#  epithelial_df <- epithelial_df[1:3000, 1:150]
  non_epithelial_ids <- as.character(
    infercnv_metadata$cell_ids[
      infercnv_metadata$cell_type != "Epithelial"
    ]
  )
  filtered_ids <- c(colnames(epithelial_df), non_epithelial_ids)
  infercnv_metadata <- infercnv_metadata[
    as.character(infercnv_metadata$cell_ids) %in% filtered_ids,
  ]
  
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
  print(paste0("Gene numbers of InferCNV gene annotation and counts matrix are now ",
    nrow(gene_annotation), " and ", nrow(epithelial_df), " respectively"))
  
  # number genes in gene annotation:
  gene_annotation$number <- seq(1, nrow(gene_annotation))

  saveRDS(
    gene_annotation, 
    paste0(common_Robject_dir, "1a.gene_annotation.Rdata")
  )
  saveRDS(
    epithelial_df, 
    paste0(common_Robject_dir, "1b.original_epithelial_df.Rdata")
  )
  saveRDS(
    infercnv_metadata, 
    paste0(common_Robject_dir, "1c.original_metadata.Rdata")
  )

} else {

  gene_annotation <- readRDS(
    paste0(common_Robject_dir, "1a.gene_annotation.Rdata")
  )
  epithelial_df <- readRDS(
    paste0(common_Robject_dir, "1b.original_epithelial_df.Rdata")
  )
  infercnv_metadata <- readRDS(
    paste0(common_Robject_dir, "1c.original_metadata.Rdata")
  )

}


################################################################################
### 2. Prepare chromosome information ###
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
### 3. Plot average fold difference from median ###
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
if (!file.exists(paste0(common_plot_dir, "1.log_original_fold_change_from_median.pdf"))) {
  p <- ggplot(log_original_fold_change_df, aes(x=number, y=count))
  p <- p + geom_point(colour = "#E8D172")
  p <- p + xlab("Genomic location")
  #p <- p + xlim(c(0, nrow(centered_original_counts_df)))
  p <- p + scale_x_continuous(
    breaks = unlist(chromosome_midpoints),
    labels = c("chr1", 2:length(chromosome_midpoints)),
    #limits = c(0,nrow(centered_original_counts_df)), 
    expand = c(0, 0)
  )
  p <- p + ylab("Log10 fold change from median")
  p <- p + scale_y_continuous(
    breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4),
    labels = c("-4", "-3", "-2", "-1", "0", "1", "2", "3", "4"),
    limits = c(-4, 4)
  )
  for (end in chromosome_ends) {
    p <- p + geom_vline(xintercept=end)
  }
  p <- p + geom_segment(
    x=0, 
    xend=max(unlist(chromosome_ends)), 
    y=log_median_original_fold_change, 
    yend=log_median_original_fold_change, 
    size=1, color="red"
  )
  p <- p + theme(
    axis.text.x=element_text(size = 24),
    axis.ticks.x=element_blank(),
    axis.text.y=element_text(size = 24),
    axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0))
  )

  pdf(paste0(common_plot_dir, "1.log10_original_fold_change_from_median.pdf"), width = 20)
    print(p)
  dev.off()

}

#save.image(paste0(Robject_dir, "temp.Rdata"))
#rm(list=ls())
#simulation_number <- 1
#Robject_dir <- paste0("/share/ScratchGeneral/jamtor/projects/thesis/Figure_2.2_accuracy_vs_coverage/results/cancer_simulation/CID4520N_cancer_sim/",
#  simulation_number, "/Rdata/")
#load(paste0(Robject_dir, "temp.Rdata"))


################################################################################
### 4. Add CNVs to df ###
################################################################################

######
#rm(initial_CNV_data)
######

# determine total gene length of genome to ensure total CNV length < 0.4 x 
# genome length to comply with InferCNV assumptions:
genome_length <- nrow(epithelial_df)

if ( !file.exists(paste0(Robject_dir, "2.initial_CNV_data.Rdata")) ) {

  print("Adding CNVs to dataset...")
  writeLines("\n")

  if (gap_or_CNV == "CNV") {

    # choose random indices for gain positions for CNVs:
    set.seed(seed_record["gain_start",])
    gene_no <- nrow(epithelial_df)
    random_gain_starts <- sample(1:gene_no, 1000)
  
    # for each length, introduce a gain of that length:
    CNV_data_list <- introduce_CNVs(feature_lengths, type = "gain", random_gain_starts, epithelial_df, 
      CNV_record = NULL, modified_df = NULL, log_modified_fold_change_df = NULL)
  
    # choose random indices for loss positions for CNVs:
    set.seed(seed_record["loss_start",])
    gene_no <- nrow(epithelial_df)
    random_loss_starts <- sample(1:gene_no, 1000)
  
    # for each length, introduce a loss of that length:
    initial_CNV_data <- introduce_CNVs(feature_lengths, type = "loss", random_starts = random_loss_starts, 
      epithelial_df, CNV_record = CNV_data_list$CNV_record, modified_df = CNV_data_list$modified_df, 
      log_modified_fold_change_df = CNV_data_list$log_modified_fold_change_df)
  
    saveRDS(initial_CNV_data, paste0(Robject_dir, "/2.initial_CNV_data.Rdata"))
    
  } else if (gap_or_CNV == "gap") {

    if (CNV_type == "gain") {

      # choose random indices for gain positions for gaps:
      set.seed(seed_record["gain_gap_start",])
      gene_no <- nrow(epithelial_df)
      random_gain_gap_starts <- sample(1:gene_no, 1000)
        
      # for each length, introduce a gap of that length:
      initial_CNV_data <- introduce_gaps(gap_lengths = feature_lengths, type = "gain", 
        random_starts = random_gain_gap_starts, epithelial_df = epithelial_df, 
        extended_gap_record = NULL, CNV_record = NULL, modified_df = NULL, 
        log_modified_fold_change_df = NULL)

    } else {

      # choose random indices for loss positions for gaps:
      set.seed(seed_record["loss_gap_start",])
      gene_no <- nrow(epithelial_df)
      random_loss_gap_starts <- sample(1:gene_no, 1000)
      # for each length, introduce a loss gap of that length:
      initial_CNV_data <- introduce_gaps(feature_lengths, type = "loss", 
        random_starts = random_loss_gap_starts, epithelial_df, 
        extended_gap_record = NULL, CNV_record = NULL, modified_df = NULL, 
        log_modified_fold_change_df = NULL)
  
    }

    saveRDS(initial_CNV_data, paste0(Robject_dir, "/2.initial_CNV_data.Rdata"))

  }

} else {
  initial_CNV_data <- readRDS(paste0(Robject_dir, "/2.initial_CNV_data.Rdata"))
}


################################################################################
### 5. Determine fold change from original median of each 
# gene and median fold change for each CNV region ###
################################################################################

if ( !file.exists(paste0(Robject_dir, "3.CNV_data.Rdata")) ) {

  print("Plotting new fold change from median for all genes...")
  writeLines("\n")


  # add non-CNV regions to CNV record:
  CNV_record_gr <- GRanges(
    seqnames = Rle("genome"),
    ranges = IRanges(start = initial_CNV_data$CNV_record$start, end = initial_CNV_data$CNV_record$end),
    strand = Rle("*"),
    multiplier = initial_CNV_data$CNV_record$multiplier,
    log_median_modified_FC = initial_CNV_data$CNV_record$log_median_modified_FC
  )
  # identify gaps between marked CNV regions as non-CNV regions and fill in values accordingly:
  non_CNV_record <- gaps(CNV_record_gr)
  non_CNV_record$multiplier = rep(1, length(non_CNV_record))
  non_CNV_record$log_median_modified_FC = rep(0, length(non_CNV_record))
  CNV_record_gr <- c(CNV_record_gr, non_CNV_record)
  # convert back to data frame and order:
  CNV_record <- data.frame(
    start = start(ranges(CNV_record_gr)),
    end = end(ranges(CNV_record_gr)),
    multiplier = CNV_record_gr$multiplier,
    log_median_modified_FC = CNV_record_gr$log_median_modified_FC
  )
  CNV_record <- CNV_record[order(CNV_record$start),]
  # fill in last non-CNV segment:
  CNV_record <- rbind(
    CNV_record,
    data.frame(
      start = CNV_record$end[nrow(CNV_record)]+1,
      end = nrow(epithelial_df),
      multiplier = 1,
      log_median_modified_FC = 0
    )
  )
  
  # add chromosome information:
  CNV_record$start_chr <- "chr1"
  CNV_record$end_chr <- "chr1"
  for (k in 1:length(chromosome_ends)) {
    if (k==1) {

      CNV_record$start_chr[
        CNV_record$start <= unlist(chromosome_ends[k])
      ] <- names(chromosome_ends)[k]

      CNV_record$end_chr[
        CNV_record$end <= unlist(chromosome_ends[k])
      ] <- names(chromosome_ends)[k]

    } else {

      CNV_record$start_chr[
        CNV_record$start <= unlist(chromosome_ends[k]) & 
        CNV_record$start > unlist(chromosome_ends[k-1])
      ] <- names(chromosome_ends)[k]

      CNV_record$end_chr[
        CNV_record$end <= unlist(chromosome_ends[k]) & 
        CNV_record$end > unlist(chromosome_ends[k-1])
      ] <- names(chromosome_ends)[k]

    }
  }
  
  if (gap_or_CNV == "CNV") {

    CNV_data <- list(
      CNV_record = CNV_record,
      modified_df = initial_CNV_data$modified_df,
      log_modified_fold_change_df = initial_CNV_data$log_modified_fold_change_df
    )

  } else if (gap_or_CNV == "gap" & CNV_type == "gain") {

    CNV_data <- list(
      extended_gap_record = initial_CNV_data$extended_gap_record,
      CNV_record = CNV_record,
      modified_df = initial_CNV_data$modified_df,
      log_modified_fold_change_df = initial_CNV_data$log_modified_fold_change_df
    )

  } else {

    CNV_data <- list(
      extended_gap_record = initial_CNV_data$extended_gap_record,
      CNV_record = CNV_record,
      modified_df = initial_CNV_data$modified_df,
      log_modified_fold_change_df = initial_CNV_data$log_modified_fold_change_df
    )
    
  }
    
  saveRDS(CNV_data, paste0(Robject_dir, "3.CNV_data.Rdata"))
  
} else {
  CNV_data <- readRDS(paste0(Robject_dir, "3.CNV_data.Rdata"))
}
  
  
################################################################################
### 6. Create modified counts and line CNV plots ###
################################################################################

# plot median fold change from original median for modified data:
if (!file.exists(paste0(plot_dir, "2a.pre_noise_log_modified_fold_change_from_median_line_only.pdf"))) {
  p <- ggplot(CNV_data$log_modified_fold_change_df, aes(x=number, y=count))
  p <- p + scale_x_continuous(
    breaks = unlist(chromosome_midpoints),
    labels = c("chr1", 2:length(chromosome_midpoints)),
    limits = c(0,length(CNV_data$log_modified_fold_change_df$count)), 
    expand = c(0, 0)
  )
  p <- p + scale_y_continuous(
    breaks = c(-3, -2, -1, 0, 1, 2, 3),
    labels = c("-3", "-2", "-1", "0", "1", "2", "3"),
    limits = c(-4, 4)
  )
  for (end in chromosome_ends) {
    p <- p + geom_vline(xintercept=end)
  }
  for (r in 1:nrow(CNV_data$CNV_record)) {
    print(r)
    # create horizontal line:
    p <- p + geom_segment(
      x=CNV_data$CNV_record$start[r], 
      xend=CNV_data$CNV_record$end[r], 
      y=CNV_data$CNV_record$log_median_modified_FC[r], 
      yend=CNV_data$CNV_record$log_median_modified_FC[r], 
      size=1, color="#37841f"
    )
    # create left vertical line:
    if (r != 1) {
      p <- p + geom_segment(
        x=CNV_data$CNV_record$start[r], 
        xend=CNV_data$CNV_record$start[r], 
        y=CNV_data$CNV_record$log_median_modified_FC[r-1], 
        yend=CNV_data$CNV_record$log_median_modified_FC[r], 
        size=1, color="#37841f"
      )
    }
  }
  p <- p + theme(
    axis.text.x=element_text(size = 24),
    axis.ticks.x=element_blank(),
    axis.text.y=element_text(size = 24),
    axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0))
  )
  pdf(paste0(plot_dir, 
    "2a.pre_noise_log_modified_fold_change_from_median_line_only.pdf"), width = 20)
    print(p)
  dev.off()
  
}

# plot counts:
if (!file.exists(paste0(plot_dir, "2b.pre_noise_log_modified_fold_change_from_median.pdf"))) {
  p <- ggplot(CNV_data$log_modified_fold_change_df, aes(x=number, y=count))
  p <- p + geom_point(colour = "#E8D172")
  p <- p + xlab("Genomic location")
  p <- p + scale_x_continuous(
    breaks = unlist(chromosome_midpoints),
    labels = c("chr1", 2:length(chromosome_midpoints)),
    limits = c(0,nrow(CNV_data$log_modified_fold_change_df)), 
    expand = c(0, 0)
  )
  p <- p + ylab("Log10 fold change")
  p <- p + scale_y_continuous(
    breaks = c(-3, -2, -1, 0, 1, 2, 3),
    labels = c("-3", "-2", "-1", "0", "1", "2", "3"),
    limits = c(-3, 3)
  )
  for (end in chromosome_ends) {
    p <- p + geom_vline(xintercept=end)
  }
  for (r in 1:nrow(CNV_data$CNV_record)) {
    print(r)
    # create horizontal line:
    p <- p + geom_segment(
      x=CNV_data$CNV_record$start[r], 
      xend=CNV_data$CNV_record$end[r], 
      y=CNV_data$CNV_record$log_median_modified_FC[r], 
      yend=CNV_data$CNV_record$log_median_modified_FC[r], 
      size=1, color="red"
    )
    # create left vertical line:
    if (r != 1) {
      p <- p + geom_segment(
        x=CNV_data$CNV_record$start[r], 
        xend=CNV_data$CNV_record$start[r], 
        y=CNV_data$CNV_record$log_median_modified_FC[r-1], 
        yend=CNV_data$CNV_record$log_median_modified_FC[r], 
        size=1, color="red"
      )
    }
  }
  p <- p + theme(
    axis.text.x=element_text(size = 24),
    axis.ticks.x=element_blank(),
    axis.text.y=element_text(size = 24),
    axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0))
  )
  pdf(paste0(plot_dir, "2b.pre_noise_log_modified_fold_change_from_median.pdf"), width = 20)
    print(p)
  dev.off()
}

#save.image(paste0(Robject_dir, "temp2.Rdata"))
#load(paste0(Robject_dir, "temp2.Rdata"))


################################################################################
### 7. Bind new epithelial and original stromal counts ###
################################################################################

# load non-epithelial_df:
non_epithelial_df <- readRDS(paste0(in_Robject_dir, "2b.non_epithelial_df.Rdata"))

# only keep genes in modified_df:
print(paste0(
  "Number of genes in non_epithelial_df before removing genes not in modified_df = ", 
  nrow(epithelial_df))
)
non_epithelial_df <- non_epithelial_df[rownames(CNV_data$modified_df),]
print(paste0(
  "Dimensions of non_epithelial_df after removing genes not in modified_df = ", 
  nrow(non_epithelial_df))
)

# bind non-epithelial df and modified_df together as new_counts:
if (gap_or_CNV == "CNV") {

  merged_CNV_data <- list(
    CNV_record = CNV_data$CNV_record,
    merged_df = cbind(
      CNV_data$modified_df,
      non_epithelial_df
    ),
    log_modified_fold_change_df = CNV_data$log_modified_fold_change_df
  )

} else if (gap_or_CNV == "gap") {

  merged_CNV_data <- list(
    extended_gap_record = CNV_data$extended_gap_record,
    CNV_record = CNV_data$CNV_record,
    merged_df = cbind(
      CNV_data$modified_df,
      non_epithelial_df
    ),
    log_modified_fold_change_df = CNV_data$log_modified_fold_change_df
  )

}

print("Dimensions of merged_df = ")
print(dim(merged_CNV_data$merged_df))


###################################################################################
### 8. Simulate noise dataset from noise counts ###
###################################################################################

if (!file.exists(paste0(noise_dir, "/noise_df.Rdata"))) {

  print(paste0("noise_df does not exist, creating now..."))

  all_cells <- read.csv(paste0(emptydrops_dir, "01_Emptydrops_out.csv"))
  filtered_cells <- read.csv(
    paste0(emptydrops_dir, "02_emptydrops_filtered_cell_ids.csv")
  )
  raw_10X <- readRDS(paste0(seurat_dir, "01_Read10X_raw_data.Rdata"))
  seurat_10X <- readRDS(paste0(seurat_dir, "03_seurat_object_processed.Rdata"))
  
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
 
  # calculate params from input to use for simulation:
  splat_params <- splatEstimate(as.matrix(noise_input))
  # reduce simulation to number of cells in simulation:
  splat_params <- setParam(splat_params, "nGenes", nrow(merged_CNV_data$merged_df))
  saveRDS(splat_params, paste0(noise_dir, "/noise_params.Rdata"))

  # simulate the data:
  sim <- splatSimulate(splat_params, batchCells = ncol(merged_CNV_data$merged_df))
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
    print(input_end_genes[g])

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

  if (!file.exists(paste0(noise_dir, "log10_noise_input_scatterplot.pdf"))) {
    pdf(paste0(noise_dir, "log10_noise_input_scatterplot.pdf"))
      p1
    dev.off()
  }


  # keep only genes in noise_counts:
  rownames(noise_counts) <- rownames(merged_CNV_data$merged_df)
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
    print(sim_end_genes[g])

    sim_df_rownames <- gsub("-|_", "", rownames(sim_df))
    sim_end_gene <- gsub("-|_", "", sim_end_genes[g])
    ind <- grep(paste0("\\b", sim_end_gene, "\\b"), sim_df_rownames)
    
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

  if (!file.exists(paste0(noise_dir, "log10_noise_sim_scatterplot.pdf"))) {
    pdf(paste0(noise_dir, "log10_noise_sim_scatterplot.pdf"))
      p2
    dev.off()
  }
  
  saveRDS(noise_counts, paste0(noise_dir, "/noise_df.Rdata"))

} else {
  print(paste0("noise_df already exists, loading..."))
  noise_counts <- readRDS(paste0(noise_dir, "/noise_df.Rdata"))
}

#save.image(paste0(Robject_dir, "temp3.Rdata"))
#load(paste0(Robject_dir, "temp3.Rdata"))


###################################################################################
### 9. Add noise to simulated counts ###
###################################################################################

if (!file.exists(paste0(Robject_dir, "4.final_CNV_data.Rdata"))) {
 
  # ensure metadata has same cells as modified_df and label stromal cells:
  infercnv_metadata <- infercnv_metadata[colnames(merged_CNV_data$merged_df),]
  epithelial_metadata <- infercnv_metadata[
    infercnv_metadata$cell_type == "Epithelial",
  ]
  
  # add noise to non-downsampled counts:
  counts_without_noise <- merged_CNV_data$merged_df
  counts_with_noise <- merged_CNV_data$merged_df
  
  # create updated CNV scatterplots post-noise addition:
  # generate mean original segment fold change vector with each value representing 
  # a gene:
  final_epithelial_df <- counts_with_noise[
    ,colnames(counts_with_noise) %in% epithelial_metadata$cell_ids[
      epithelial_metadata$cell_type == "Epithelial"
    ]
  ]
  
  final_CNV_record <- merged_CNV_data$CNV_record
  final_CNV_record$log_median_modified_FC_post_noise <- NA

  for (i in 1:nrow(final_CNV_record)) {
  
    if (final_CNV_record$start[i] != final_CNV_record$end[i]) {
      average_original_counts <- apply(
        epithelial_df[final_CNV_record$start[i]:final_CNV_record$end[i],], 1, mean
      )
    } else {
      average_original_counts <- mean(epithelial_df[final_CNV_record$start[i],])
    }
  
    # add 0.1 to all values:
    #average_original_counts[average_original_counts == 0] <- 1e-3
    average_original_counts <- average_original_counts + 1e-3
    # determine median:
    median_average_original_counts <- median(average_original_counts)
    # divide by median to get fold change from median:
    original_fold_change <- average_original_counts/median_average_original_counts
    # check median of original fold change = 1:
    median_original_fold_change <- median(original_fold_change)
    # generate mean modified fold change from original median vector with each value 
    # representing a gene:
    if (final_CNV_record$start[i] != final_CNV_record$end[i]) {
      average_modified_counts <- apply(
        final_epithelial_df[final_CNV_record$start[i]:final_CNV_record$end[i],], 1, mean
      )
    } else {
      average_modified_counts <- mean(final_epithelial_df[final_CNV_record$start[i],])
    }
    # add 0.1 to zero values:
    average_modified_counts <- average_modified_counts + 1e-3
    # divide by median to get fold change from original median and add to df for plotting:
    modified_fold_change <- average_modified_counts/median_average_original_counts
    # take log of the median of modified fold change to mark on plot:
    median_modified_fold_change <- median(modified_fold_change)    
    log_median_modified_fold_change <- log10(median_modified_fold_change)
    # add to final_CNV_record:
    final_CNV_record$log_median_modified_FC_post_noise[i] <- log_median_modified_fold_change
    # take the log10 of all fold changes:
    log_modified_fold_change <- log10(modified_fold_change)
    # add modified fold change values to log_original_fold_change_df:
    if (!exists("log_modified_fold_change_df")) {
      log_modified_fold_change_df <- log_original_fold_change_df
    }
    log_modified_fold_change_df$count[
      final_CNV_record$start[i]:final_CNV_record$end[i]
    ] <- log_modified_fold_change
  }

  # add noised data:
  if (gap_or_CNV == "CNV") {

    noised_CNV_data <- list(
      CNV_record = final_CNV_record,
      noised_df = counts_with_noise,
      log_modified_fold_change_df = merged_CNV_data$log_modified_fold_change_df,
      infercnv_metadata = infercnv_metadata
    )

  } else if (gap_or_CNV == "gap") {

    noised_CNV_data <- list(
      extended_gap_record = merged_CNV_data$extended_gap_record,
      CNV_record = final_CNV_record,
      noised_df = counts_with_noise,
      log_modified_fold_change_df = merged_CNV_data$log_modified_fold_change_df,
      infercnv_metadata = infercnv_metadata
    )

  }

  saveRDS(noised_CNV_data, paste0(Robject_dir, "4.final_CNV_data.Rdata"))

} else {
  noised_CNV_data <- readRDS(paste0(Robject_dir, "4.final_CNV_data.Rdata"))
}


###################################################################################
### 10. Plot post-noise gene expression profiles ###
###################################################################################

# plot median fold change from original median for modified data:
if (
  !file.exists(
    paste0(
      plot_dir, "3a.log_modified_fold_change_from_median_with_noise_line_only.png"
      )
    )
) {
  p <- ggplot(
    noised_CNV_data$log_modified_fold_change_df, aes(x=number, y=count)
  )
  p <- p + scale_x_continuous(
    breaks = unlist(chromosome_midpoints),
    labels = c("chr1", 2:length(chromosome_midpoints)),
    limits = c(0,length(noised_CNV_data$log_modified_fold_change_df$count)), 
    expand = c(0, 0)
  )
  p <- p + scale_y_continuous(
      breaks = c(-3, -2, -1, 0, 1, 2, 3),
      labels = c("-3", "-2", "-1", "0", "1", "2", "3"),
      limits = c(-3, 3)
    )
  for (end in chromosome_ends) {
    p <- p + geom_vline(xintercept=end)
  }
  for (r in 1:nrow(noised_CNV_data$CNV_record)) {
    print(r)
    # create horizontal line:
    p <- p + geom_segment(
      x=noised_CNV_data$CNV_record$start[r], 
      xend=noised_CNV_data$CNV_record$end[r], 
      y=noised_CNV_data$CNV_record$log_median_modified_FC_post_noise[r], 
      yend=noised_CNV_data$CNV_record$log_median_modified_FC_post_noise[r], 
      size=1, color="#37841f"
    )
    # create left vertical line:
    if (r != 1) {
      p <- p + geom_segment(
        x=noised_CNV_data$CNV_record$start[r], 
        xend=noised_CNV_data$CNV_record$start[r], 
        y=noised_CNV_data$CNV_record$log_median_modified_FC_post_noise[r-1], 
        yend=noised_CNV_data$CNV_record$log_median_modified_FC_post_noise[r], 
        size=1, color="#37841f"
      )
    }
  }
  p <- p + theme(
    axis.text.x=element_text(size = 24),
    axis.ticks.x=element_blank(),
    axis.text.y=element_text(size = 24),
    axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0))
  )
  png(
    paste0(plot_dir, 
      "3a.log_modified_fold_change_from_median_with_noise_line_only.png"
      ), 
      width = 20,
      res = 300,
      units = "in"
    )
    print(p)
  dev.off()
  
}

# plot counts:
if (
  !file.exists(
    paste0(plot_dir, "3b.log_modified_fold_change_from_median_with_noise.pdf")
    )
  ) {
  p <- ggplot(
    noised_CNV_data$log_modified_fold_change_df, aes(x=number, y=count)
  )
  p <- p + geom_point(colour = "#E8D172")
  p <- p + xlab("Genomic location")
  p <- p + scale_x_continuous(
    breaks = unlist(chromosome_midpoints),
    labels = c("chr1", 2:length(chromosome_midpoints)),
    limits = c(0,nrow(noised_CNV_data$log_modified_fold_change_df)), 
    expand = c(0, 0)
  )
  p <- p + ylab("Log10 fold change")
  p <- p + scale_y_continuous(
     breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4),
     labels = c("-4", "-3", "-2", "-1", "0", "1", "2", "3", "4"),
     limits = c(
      min(noised_CNV_data$log_modified_fold_change_df$count), 
      max(noised_CNV_data$log_modified_fold_change_df$count)
    )
   )
  for (end in chromosome_ends) {
    p <- p + geom_vline(xintercept=end)
  }
  for (r in 1:nrow(noised_CNV_data$CNV_record)) {
    print(r)
    # create horizontal line:
    p <- p + geom_segment(
      x=noised_CNV_data$CNV_record$start[r], 
      xend=noised_CNV_data$CNV_record$end[r], 
      y=noised_CNV_data$CNV_record$log_median_modified_FC_post_noise[r], 
      yend=noised_CNV_data$CNV_record$log_median_modified_FC_post_noise[r], 
      size=1, color="red"
    )
    # create left vertical line:
    if (r != 1) {
      p <- p + geom_segment(
        x=noised_CNV_data$CNV_record$start[r], 
        xend=noised_CNV_data$CNV_record$start[r], 
        y=noised_CNV_data$CNV_record$log_median_modified_FC_post_noise[r-1], 
        yend=noised_CNV_data$CNV_record$log_median_modified_FC_post_noise[r], 
        size=1, color="red"
      )
    }
  }
  p <- p + theme(
    axis.text.x=element_text(size = 24),
    axis.ticks.x=element_blank(),
    axis.text.y=element_text(size = 24),
    axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0))
  )
  pdf(
    paste0(plot_dir, "3b.log_modified_fold_change_from_median_with_noise.pdf"), 
    width = 20
  )
    print(p)
  dev.off()

}
  
if (!file.exists(paste0(out_dir, "input_matrix.txt"))) {
  system(paste0("mkdir -p ", out_dir))
  write.table(
    noised_CNV_data$noised_df, 
    paste0(out_dir, "input_matrix.txt"), 
    sep = "\t", quote = F, col.names=T, row.names=T)
  write.table(noised_CNV_data$infercnv_metadata, 
    paste0(out_dir, "metadata.txt"), 
    sep = "\t", quote = F, col.names = F, row.names = F)
}


################################################################################
### 12. Plot normal and simulated expression together ###
################################################################################

# replot normal expression profile without x-axis text/title:
p <- ggplot(log_original_fold_change_df, aes(x=number, y=count))
p <- p + geom_point(colour = "#E8D172")
p <- p + xlab("Genomic location")
#p <- p + xlim(c(0, nrow(centered_original_counts_df)))
p <- p + scale_x_continuous(
  breaks = unlist(chromosome_midpoints),
  labels = c("chr1", 2:length(chromosome_midpoints)),
  #limits = c(0,nrow(centered_original_counts_df)), 
  expand = c(0, 0)
)
p <- p + ylab("Log10 fold change from median")
if (gap_or_CNV == "CNV") {
  p <- p + scale_y_continuous(
    breaks = c(-3, 0, 3),
    labels = c("-3", "0", "3"),
    limits = c(-3, 3)
  )
} else {
  p <- p + scale_y_continuous(
    breaks = c(-1.5, 0, 1.5),
    labels = c("-1.5", "0", "1.5"),
    limits = c(-1.5, 1.5)
  )
}
for (end in chromosome_ends) {
  p <- p + geom_vline(xintercept=end)
}
p <- p + geom_segment(
  x=0, 
  xend=max(unlist(chromosome_ends)), 
  y=log_median_original_fold_change, 
  yend=log_median_original_fold_change, 
  size=1, color="red"
)
p <- p + theme(
  axis.text.x=element_blank(),
  axis.ticks.x=element_blank(),
  axis.title.x = element_blank(),
  axis.text.y=element_text(size = 30),
  axis.title.y = element_blank()
)
normal_grid <- ggplotGrob(p)
dev.off()

# replot sim expression profile:
p <- ggplot(
  noised_CNV_data$log_modified_fold_change_df, aes(x=number, y=count)
)
p <- p + geom_point(colour = "#E8D172")
p <- p + xlab("Genomic location")
p <- p + scale_x_continuous(
  breaks = unlist(chromosome_midpoints),
  labels = c("chr1", 2:20, "\n21", 22),
  limits = c(0,nrow(noised_CNV_data$log_modified_fold_change_df)), 
  expand = c(0, 0)
)
p <- p + ylab("Log10 fold change")
if (gap_or_CNV == "CNV") {
  p <- p + scale_y_continuous(
    breaks = c(-3, 0, 3),
    labels = c("-3", "0", "3"),
    limits = c(-3, 3)
  )
} else {
  p <- p + scale_y_continuous(
    breaks = c(-1.5, 0, 1.5),
    labels = c("-1.5", "0", "1.5"),
    limits = c(-1.5, 1.5)
  )
}
for (end in chromosome_ends) {
  p <- p + geom_vline(xintercept=end)
}
for (r in 1:nrow(noised_CNV_data$CNV_record)) {
  print(r)
  # create horizontal line:
  p <- p + geom_segment(
    x=noised_CNV_data$CNV_record$start[r], 
    xend=noised_CNV_data$CNV_record$end[r], 
    y=noised_CNV_data$CNV_record$log_median_modified_FC_post_noise[r], 
    yend=noised_CNV_data$CNV_record$log_median_modified_FC_post_noise[r], 
    size=1, color="red"
  )
  # create left vertical line:
  if (r != 1) {
    p <- p + geom_segment(
      x=noised_CNV_data$CNV_record$start[r], 
      xend=noised_CNV_data$CNV_record$start[r], 
      y=noised_CNV_data$CNV_record$log_median_modified_FC_post_noise[r-1], 
      yend=noised_CNV_data$CNV_record$log_median_modified_FC_post_noise[r], 
      size=1, color="red"
    )
  }
}
p <- p + theme(
  axis.text.x=element_text(size = 27),
  axis.ticks.x=element_blank(),
  axis.text.y=element_text(size = 30),
  axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 32),
  axis.title.y = element_blank()
)
sim_grid <- ggplotGrob(p)
dev.off()

# plot together:
png(
  paste0(plot_dir, "normal_vs_sim_expression_plot.png"), 
  height = 12, 
  width = 23,
  res = 300,
  units = "in"
)   
  grid.newpage()
    pushViewport(viewport(x = 0.57, y = 0.73, width = 0.84, height = 0.43))
      #grid.rect()
      grid.draw(normal_grid)
    popViewport()
    pushViewport(viewport(x = 0.57, y = 0.25, width = 0.84, height = 0.5))
      #grid.rect()
      grid.draw(sim_grid)
    popViewport()
    
    # draw label text:
    pushViewport(viewport(x = 0.07, y = 0.75, width = 0.13, height = 0.2, just = "right"))
      #grid.rect()
      grid.text("Normal\nbreast", gp=gpar(fontsize=33, fontface = "bold"), just = "left")
    popViewport()
  
    pushViewport(viewport(x = 0.07, y = 0.29, width = 0.13, height = 0.2, just = "right"))
      #grid.rect()
      grid.text("CNV\nsimulation", gp=gpar(fontsize=33, fontface = "bold"), just = "left")
    popViewport()

    pushViewport(viewport(x = 0.14, y = 0.52, width = 0.01, height = 0.7, just = "right"))
      #grid.rect()
      grid.text("Fold change from median", gp=gpar(fontsize=30), rot = 90)
    popViewport()

dev.off()


