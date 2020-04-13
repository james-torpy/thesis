#! /share/ClusterShare/software/contrib/briglo/octR/src/R-3.6.0/builddir/bin/Rscript


###############################################################################
### Simulates cancer datasets taking as inputs number of CNVs, CNV lengths, ###
### gain/loss multipliers and downsamples required                          ###
###############################################################################

args = commandArgs(trailingOnly=TRUE)

project_name <- "thesis"
subproject_name <- "Figure_2.6_random_permutation_of_normal"
sample_name <- args[1]
print(paste0("sample name = ", sample_name))
nUMI_threshold <- as.numeric(args[2])
print(paste0("nUMI_threshold = ", nUMI_threshold))
nGene_threshold <- as.numeric(args[3])
print(paste0("nGene_threshold = ", nGene_threshold))
# proportion of genes to be randomly permuated:
permutation_proportion <- as.numeric(args[4])
print(paste0("permutation proportion = ", permutation_proportion))
# range of values which could be substituted:
potential_values <- as.numeric(
  unlist(
    strsplit(
      args[5],
      split = "_"
    )
  )
)
print("possible permutation values: ")
print(potential_values)
simulation_number <- as.numeric(args[6])
print(paste0("Simulation number = ", simulation_number))
t_cells_included <- as.logical(args[7])
print(paste0("T-cells included in normal InferCNV run? ", 
  t_cells_included))
analysis_mode <- args[8]
print(paste0("Analysis mode of normal InferCNV run = ", 
  analysis_mode, "_mode"))
neutral_signal_range <- unlist(strsplit(args[9], split = "_"))
print("Neutral signal range = ")
print(neutral_signal_range)

#project_name <- "thesis"
#subproject_name <- "Figure_2.6_random_permutation_of_normal"
#sample_name <- "CID4520N"
#print(paste0("sample name = ", sample_name))
#nUMI_threshold <- 25000
#print(paste0("nUMI_threshold = ", nUMI_threshold))
#nGene_threshold <- 5000
#print(paste0("nGene_threshold = ", nGene_threshold))
# proportion of genes to be randomly permuated:
#permutation_proportion <- 0.005
#print(paste0("permutation proportion = ", permutation_proportion))
# range of values which could be substituted:
#potential_values <- as.numeric(
#  unlist(
#    strsplit(
#      "0_0.5_1.5_2_3",
#      split = "_"
#    )
#  )
#)
#print("possible permutation values: ")
#print(potential_values)
#simulation_number <- 1
#print(paste0("Simulation number = ", simulation_number))
#t_cells_included <- TRUE
#print(paste0("T-cells included in normal InferCNV run? ", 
#  t_cells_included))
#analysis_mode <- "samples"
#print(paste0("Analysis mode of normal InferCNV run = ", 
#  analysis_mode, "_mode"))
#neutral_signal_range <- unlist(strsplit("0.97_1.03", split = "_"))
#print("Neutral signal range = ")
#print(neutral_signal_range)

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
library(data.table)
library(ComplexHeatmap, lib.loc = lib_loc)
library(circlize, lib.loc = lib_loc)

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

in_dir <- paste0(project_dir, "/raw_files/seurat_objects/",
  sample_name, "/")

if (t_cells_included) {
  normal_path <- paste0(results_dir, "infercnv/t_cells_included/", 
    "/", sample_name, "/")
} else {
  normal_path <- paste0(results_dir, "infercnv/t_cells_excluded", 
    "/", sample_name, "/")
}
normal_Robject_dir <- paste0(normal_path, "Rdata/")
normal_input_dir <- paste0(normal_path, "input_files/")
normal_infercnv_dir <- paste0(normal_path, analysis_mode, "_mode/")

out_path <- paste0(results_dir, "permutated_", sample_name, "/")

common_Robject_dir <- paste0(out_path, "Rdata/")
system(paste0("mkdir -p ", common_Robject_dir))
common_plot_dir <- paste0(out_path, "plots/")
system(paste0("mkdir -p ", common_plot_dir))
common_table_dir <- paste0(out_path, "tables/")
system(paste0("mkdir -p ", common_table_dir))

permutated_out_path <- paste0(out_path, permutation_proportion, 
  "_proportion/", simulation_number, "/remove_genes_1/")
Robject_dir <- paste0(permutated_out_path, "/Rdata/")
system(paste0("mkdir -p ", Robject_dir))
plot_dir <- paste0(permutated_out_path, "/plots/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(permutated_out_path, "/tables/")
system(paste0("mkdir -p ", table_dir))

prop_table_dir <- paste0(out_path, permutation_proportion, 
  "_proportion/tables/")
system(paste0("mkdir -p ", prop_table_dir))

if (t_cells_included) {
  out_dir <- paste0(results_dir, "infercnv/t_cells_included/", 
    "permutated_", sample_name, "/", permutation_proportion, 
    "_proportion/", simulation_number, "/remove_genes_1/input_files/")
} else {
  out_dir <- paste0(results_dir, "infercnv/t_cells_excluded/",  
    "permutated_", sample_name, "/", permutation_proportion, 
    "_proportion/", simulation_number, "/remove_genes_1/input_files/")
}

system(paste0("mkdir -p ", out_dir))

print(paste0("Sample directory = ", in_dir))
print(paste0("Reference directory = ", ref_dir))

print(paste0("Robject directory = ", Robject_dir))
print(paste0("Plot dir = ", plot_dir))
print(paste0("Table directory = ", table_dir))

print(paste0("Out (InferCNV) directory = ", out_dir))

print(paste0("Generating permutated normal data set from ", sample_name))
print(paste0("Filtering out cells with less than ", nUMI_threshold, " UMIs and ",
  nGene_threshold, " genes"))


################################################################################
### 0. Define functions and choose seeds ###
################################################################################

prepare_infercnv_metadata <- dget(paste0(func_dir, "prepare_infercnv_metadata.R"))
fetch_chromosome_boundaries <- dget(paste0(func_dir, 
  "fetch_chromosome_boundaries.R"))
#introduce_permuatations <- dget(paste0(func_dir, "introduce_permutations.R"))

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
    row.names = c("permutation_indices", "potential_values"),
    seed = sample(1:999, 2)
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


###################################################################################
### 1. Determine artefact regions from normal InferCNV output to be removed from
# the dataset before random permutation ###
###################################################################################

if (
  !file.exists(
    paste0(common_Robject_dir, "0a.normal_artefact_record.Rdata")
  ) | 
  !file.exists(
    paste0(common_Robject_dir, "0b.normal_artefact_genes.Rdata")
  )
) {

  # load normal infercnv output:
  print("Loading normal InferCNV output files...")
  normal_infercnv <- as.data.frame(t(read.table(paste0(normal_infercnv_dir, 
    "infercnv.12_denoised.observations.txt"))))
  
  # determine average signal across all cells for each gene:
  normal_averages <- apply(normal_infercnv, 2, mean)
  
  for ( i in 1:(ncol(normal_infercnv)-9) ) {

    averages_per_cell <- apply(normal_infercnv[i:(i+9)], 1, mean)
    
    # if at least 10% of cells have average signal in loss range, label as loss:
    if ( length(which(averages_per_cell < neutral_signal_range[1])) > 
      floor(0.1*length(averages_per_cell)) ) {

      segment_record <- data.frame(
        start = i,
        end = i+9,
        average_signal = mean(normal_averages[i:(i+9)]),
        type = "loss"
      )

    } else if ( length(which(averages_per_cell > neutral_signal_range[2])) > 
      floor(0.1*length(averages_per_cell)) ) {

      segment_record <- data.frame(
        start = i,
        end = i+9,
        average_signal = mean(normal_averages[i:(i+9)]),
        type = "gain"
      )

    } else {

      segment_record <- data.frame(
        start = i,
        end = i+9,
        average_signal = mean(normal_averages[i:(i+9)]),
        type = "no_CNV"
      )

    }

    if (i==1) {
      artefact_record <- segment_record
    } else {
      artefact_record <- rbind(artefact_record, segment_record)
    }
    
    artefact_record$type <- as.character(artefact_record$type)

  }

  split_record <- split(artefact_record, rleid(artefact_record$type))
  final_artefact_record <- lapply(split_record, function(x) {
    return(
      data.frame(
        start = x$start[1],
        end = x$end[nrow(x)],
        type = x$type[1]
      )
    )
  })
  final_artefact_record <- do.call("rbind", final_artefact_record)

  artefact_indices <- final_artefact_record[final_artefact_record$type != "no_CNV",]

  # fetch chromosome info:
  chr_data <- fetch_chromosome_boundaries(normal_infercnv, ref_dir)

  # identify chromosomes artefacts belong to:
  for (k in 1:length(chr_data$ends)) {
    if (k==1) {

      artefact_indices$start_chr[
        artefact_indices$start <= chr_data$ends[k]
      ] <- names(chr_data$ends)[k]

      artefact_indices$end_chr[
        artefact_indices$end <= chr_data$ends[k]
      ] <- names(chr_data$ends)[k]

    } else {

      artefact_indices$start_chr[
        artefact_indices$start <= chr_data$ends[k] & 
        artefact_indices$start > chr_data$ends[k-1]
      ] <- names(chr_data$ends)[k]

      artefact_indices$end_chr[
        artefact_indices$end <= chr_data$ends[k] & 
        artefact_indices$end > chr_data$ends[k-1]
      ] <- names(chr_data$ends)[k]

    }
  }

  for (m in 1:nrow(artefact_indices)) {
    if (m==1) {
      normal_artefact <- artefact_indices$start[m]:artefact_indices$end[m]
    } else {
      normal_artefact <- c(normal_artefact,
        artefact_indices$start[m]:artefact_indices$end[m]
      )
    }   
  }

  # identify non-artefact genes:
  artefact_genes <- colnames(normal_infercnv)[normal_artefact]
  # remove duplicates:
  artefact_genes <- artefact_genes[!duplicated(artefact_genes)]

  filtered_infercnv <- normal_infercnv[
    , !(colnames(normal_infercnv) %in% artefact_genes)
  ]

  print(
    paste0(
      "Dimensions of filtered_infercnv = ", 
        paste(as.character(dim(filtered_infercnv)), collapse=",")
    )
  )

  # plot artefact-free normal heatmap:  
  # prepare df for plotting:
  plot_object <- filtered_infercnv
  colnames(plot_object) <- rep("la", ncol(plot_object))
  plot_object <- as.matrix(plot_object)
  
  # define heatmap colours:
  na_less_vector <- unlist(plot_object)
  na_less_vector <- na_less_vector[!is.na(na_less_vector)]
  heatmap_cols <- colorRamp2(c(min(na_less_vector), 1, max(na_less_vector)), 
        c("#00106B", "white", "#680700"), space = "sRGB")

  print("Generating filtered normal heatmap...")

  final_heatmap <- Heatmap(
    plot_object, name = paste0("hm"), 
    col = heatmap_cols,
    cluster_columns = F, cluster_rows = F,
    show_row_names = F, show_column_names = F,
    #column_names_gp = gpar(col = "white"),
    show_row_dend = F,
    heatmap_legend_param = list(title = "CNV\nscore", 
    at = c(round(min(na_less_vector), 1), 1, round(max(na_less_vector), 1)),
    color_bar = "continuous", grid_height = unit(1.5, "cm"), 
    grid_width = unit(1.5, "cm"), legend_direction = "horizontal",
    title_gp = gpar(fontsize = 18, fontface = "bold"), 
    labels_gp = gpar(fontsize = 12)),
    use_raster = T, raster_device = c("png")
  )

  annotated_heatmap <- grid.grabExpr(
    draw(final_heatmap, gap = unit(6, "mm"), heatmap_legend_side = "left")
  )
  dev.off()

  # plot final annotated heatmap:
  pdf(
    paste0(common_plot_dir, "filtered_normal_infercnv_plot.pdf"), 
    height = 13, width = 20
  )   
    
    grid.newpage()
    pushViewport(viewport(x = 0.01, y = 0.16, width = 0.99, height = 0.78, 
      just = c("left", "bottom")))
      grid.draw(annotated_heatmap)
      decorate_heatmap_body("hm", {
        for ( e in 1:length(chr_data$end_pos) ) {
          grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), 
            gp = gpar(lwd = 1, col = "#383838"))
          grid.text(gsub("chr", "", names(chr_data$lab_pos)[e]), chr_data$lab_pos[e], 
            unit(0, "npc") + unit(-3.5, "mm"), gp=gpar(fontsize=18))
        }
      })
    popViewport()
  
  dev.off()

  saveRDS(
    final_artefact_record, 
    paste0(common_Robject_dir, "0a.normal_artefact_record.Rdata")
  )
  saveRDS(
    artefact_genes, 
    paste0(common_Robject_dir, "0b.normal_artefact_genes.Rdata")
  )

} else {

  final_artefact_record <- readRDS(
    paste0(common_Robject_dir, "0a.normal_artefact_record.Rdata")
  )
  artefact_genes <- readRDS(
    paste0(common_Robject_dir, "0b.normal_artefact_genes.Rdata")
  )

}


###################################################################################
### 2. Fetch counts matrix from normal sample, create original infercnv input files 
# and plot nUMI and nGene ###
###################################################################################

if (
  !file.exists(paste0(Robject_dir, "1a.epithelial_df.Rdata")) |
  !file.exists(paste0(Robject_dir, "1b.non_epithelial_df.Rdata")) |
  !file.exists(paste0(Robject_dir, "1c.metadata.Rdata"))
) {
  
  input_list <- list(
    readRDS(paste0(normal_Robject_dir, "2a.epithelial_df.Rdata"))
  )
  input_list[[2]] <- readRDS(paste0(normal_Robject_dir, "2b.non_epithelial_df.Rdata"))
  input_list[[3]] <- readRDS(paste0(normal_Robject_dir, "2c.gene_annotation.Rdata"))
  names(input_list) <- c("epithelial_df", "non_epithelial_df", "gene_annotation")

  infercnv_metadata <- read.table(
    paste0(normal_input_dir, "metadata.txt"),
    sep = "\t",
    header = F,
    as.is = T
  )
  colnames(infercnv_metadata) <- c("cell_ids", "cell_types")

  # remove artefact regions for all dfs:
  filtered_list <- lapply(input_list, function(x) {
    df <- x[!(rownames(x) %in% artefact_genes),]
  })

  epithelial_df <- filtered_list$epithelial_df
  non_epithelial_df <- filtered_list$non_epithelial_df
  gene_annotation <- filtered_list$gene_annotation

  print(
    paste0(
      "Dimensions of epithelial df after filtering out initial artefacts = ", 
        paste(as.character(dim(epithelial_df)), collapse=",")
    )
  )

  print(
    paste0(
      "Dimensions of non-epithelial df = ", 
        paste(as.character(dim(non_epithelial_df)), collapse=",")
    )
  )

  print(
    paste0(
      "Dimensions of gene annotation = ", 
        paste(as.character(dim(gene_annotation)), collapse=",")
    )
  )

  print(
    paste0(
      "Dimensions of infercnv metadata = ", 
        paste(as.character(dim(infercnv_metadata)), collapse=",")
    )
  )

  saveRDS(epithelial_df, paste0(Robject_dir, 
    "/1a.epithelial_df.Rdata"))
  saveRDS(non_epithelial_df, paste0(Robject_dir, 
    "/1b.non_epithelial_df.Rdata"))
  saveRDS(gene_annotation, paste0(Robject_dir, 
    "/1c.gene_annotation.Rdata"))
  saveRDS(infercnv_metadata, paste0(Robject_dir, 
    "/1d.metadata.Rdata"))

} else {

  epithelial_df <- readRDS(paste0(Robject_dir, 
    "/1a.epithelial_df.Rdata"))
  non_epithelial_df <- readRDS(paste0(Robject_dir, 
    "/1b.non_epithelial_df.Rdata"))
  gene_annotation <- readRDS(paste0(Robject_dir, 
    "/1c.gene_annotation.Rdata"))
  infercnv_metadata <- readRDS(paste0(Robject_dir, 
    "/1d.metadata.Rdata"))

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

print(
  paste0(
    "Dimensions of log_original_fold_change_df = ", 
      paste(as.character(dim(log_original_fold_change_df)), collapse=",")
  )
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
    labels = paste0("chr", 1:length(chromosome_midpoints)),
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

  pdf(paste0(common_plot_dir, "1.log10_original_fold_change_from_median.pdf"), width = 20)
    print(p)
  dev.off()

}


################################################################################
### 4. Permutate genes randomly in df ###
################################################################################

if (
  !file.exists(paste0(Robject_dir, "2.initial_permutation_data.Rdata")) | 
  !file.exists(paste0(prop_table_dir, "gene_numbers.txt"))
) {

  if (permutation_proportion != 0) {

    # determine number of genes to permutate:
    gene_no <- nrow(epithelial_df)
    permutation_gene_no <- round(permutation_proportion*gene_no, 0)
  
    # choose random indices for positions of permutations:
    set.seed(seed_record["permutation_indices",])
    permutation_record <- data.frame(
      indices = sample(1:gene_no, permutation_gene_no)
    )
  
    # choose permutation values to multiply counts by:
    set.seed(seed_record["potential_values",])
    permutation_record$values <- sample(
      potential_values, 
      permutation_gene_no,
      replace = T
    )
  
    # order by indices:
    permutation_record <- permutation_record[order(permutation_record$indices),]
    # add gene names:
    permutation_record$genes <- rownames(epithelial_df)[permutation_record$indices]
  
    # check portion of epithelial_df corresponding to permutation indices:
    for (i in 1:3) {
      eg <- data.frame(
        head(epithelial_df[permutation_record$indices,][i,], 20)
      )
      colnames(eg) <- rownames(epithelial_df)[permutation_record$indices][i]
      print("E.g. portions of original epithelial_df to be changed: ")
      print(eg)
    }
  
    # if any genes corresponding with permutation indices are zero, 
    # substitute with 1:
    modified_df <- epithelial_df
    modified_df[permutation_record$indices,][
      modified_df[permutation_record$indices,] == 0
    ] <- 1
  
    # multiply columns of permutation indices by the corresponding 
    # permutation values:
    modified_df[permutation_record$indices,] <- 
      modified_df[permutation_record$indices,] * permutation_record$values
  
    # check portion of modified_df corresponding to permutation indices:
    for (i in 1:3) {
      print(paste0("Check the following portion was multiplied by ", 
        permutation_record$values[i], ":"))
      eg <- data.frame(
        head(modified_df[permutation_record$indices,][i,], 20)
      )
      colnames(eg) <- rownames(modified_df)[permutation_record$indices][i]
      print(eg)
    }
  
    # generate mean original fold change vector with each value representing 
    # a gene:
    average_original_counts <- apply(epithelial_df, 1, mean)
    # add 0.1 to all values:
    #average_original_counts[average_original_counts == 0] <- 1e-3
    average_original_counts <- average_original_counts + 1e-4
    # determine median:
    median_average_original_counts <- median(average_original_counts)
    # divide by median to get fold change from median:
    original_fold_change <- average_original_counts/median_average_original_counts
    # check median of original fold change = 1:
    median_original_fold_change <- median(original_fold_change)
    # generate mean modified fold change from original median vector with each value 
    # representing a gene:
    average_modified_counts <- apply(modified_df, 1, mean)
    # add 0.1 to all values to account for zeros:
    average_modified_counts <- average_modified_counts + 1e-4
    # divide by median to get fold change from original median and add to df for plotting:
    modified_fold_change <- average_modified_counts/median_average_original_counts
    # take log of the median of modified fold change to mark on plot:
    median_modified_fold_change <- median(modified_fold_change)    
    log_median_modified_fold_change <- log10(median_modified_fold_change)
    # take the log10 of all fold changes:
    log_modified_fold_change <- log10(modified_fold_change)
    # add modified fold change values to log_original_fold_change_df:
    log_modified_fold_change_df <- log_original_fold_change_df
    log_modified_fold_change_df$count <- log_modified_fold_change
  
    permutation_data <- list(
      permutation_record,
      modified_df,
      log_modified_fold_change_df,
      log_median_modified_fold_change
    )
    names(permutation_data) <- c("permutation_record", "modified_df", 
      "log_modified_fold_change_df", "log_median_modified_fold_change")
  
    print(
      paste0(
        "Dimensions of modified_df = ", 
          paste(as.character(dim(modified_df)), collapse=",")
      )
    )

    gene_numbers <- data.frame(
      row.names = c("permutated_genes", "total_genes"),
      number = c(ncol(modified_df), nrow(permutation_data$log_modified_fold_change_df))
    )
    print(gene_numbers)

  } else {

    print("Permutation proportion = 0, not permutating any genes")

    permutation_data <- list(
      epithelial_df,
      log_original_fold_change_df,
      log_median_original_fold_change
    )
    names(permutation_data) <- c("modified_df", 
      "log_modified_fold_change_df", "log_median_modified_fold_change")

    gene_numbers <- data.frame(
      row.names = c("permutated_genes", "total_genes"),
      number = c(0, ncol(epithelial_df))
    )
  
  }

  saveRDS(
    permutation_data, paste0(Robject_dir, "/2.initial_permutation_data.Rdata")
  )

  write.table(
    gene_numbers, 
    paste0(prop_table_dir, "gene_numbers.txt"),
    sep = "\t",
    quote = F,
    row.names = T,
    col.names = F
  )

} else {

  permutation_data <- readRDS(
    paste0(Robject_dir, "/2.initial_permutation_data.Rdata")
  )

  gene_numbers <- read.table(
    paste0(prop_table_dir, "gene_numbers.txt"),
    sep = "\t",
    header = F,
    as.is = T
  )

}
  
  
################################################################################
### 5. Create modified counts and line CNV plots ###
################################################################################

# plot median fold change from original median for modified data:
if (
  !file.exists(
    paste0(
      plot_dir, "2b.pre_noise_log_modified_fold_change_from_median.pdf"
    )
  )
) {

  # create colour vector to differentiate modified values:
  plot_df <- permutation_data$log_modified_fold_change_df
  plot_df$colour <- "non-modified"
  plot_df$colour[
    permutation_data$permutation_record$indices
  ] <- "modified"

  # label original values of modified genes by changing the value/
  # colour of the gene directly following it (bit hacky/cheaty),
  # or before it if it's the last gene in the genome:

  # fetch indices of modified genes to fetch original gene values:
  orig_gene_indices <- plot_df$number[plot_df$colour == "modified"]
  names(orig_gene_indices) <- rownames(plot_df)[plot_df$number %in% orig_gene_indices]
  
  # define mock genes, 1 after or before the modified genes, to change to the original value
  # of modified genes:
  mock_orig_indices <- orig_gene_indices
  mock_orig_indices[
    orig_gene_indices < nrow(epithelial_df)
  ] <- mock_orig_indices[
    mock_orig_indices < nrow(epithelial_df)
  ] + 1

  mock_orig_indices[
    mock_orig_indices == nrow(epithelial_df)
  ] <- mock_orig_indices[
    mock_orig_indices == nrow(epithelial_df)
  ] - 1

  for (m in 1:length(mock_orig_indices)) {
    # assign the gene adjacent to the modified gene the original value of the modified gene:
    plot_df$count[plot_df$number == mock_orig_indices[m]] <- 
      log_original_fold_change_df[names(mock_orig_indices)[m],]$count
    # adjust the colour of this mock original gene to 'original value'
    plot_df$colour[plot_df$number == mock_orig_indices[m]] <- "original value"
  }

  # adjust the order of the df so modified and original values are printed last:
  plot_df <- rbind(
    plot_df[plot_df$colour == "non-modified",],
    plot_df[plot_df$colour == "modified",],
    plot_df[plot_df$colour == "original value",]
  )
  plot_df$colour <- factor(
    plot_df$colour, levels = c("non-modified", "modified", "original value")
  )
  
  p <- ggplot(
    plot_df, 
    aes(x=number, y=count, color=colour)
  )
  p <- p + geom_point()
  p <- p + scale_color_manual(values=c("#F0F2F7", "#46B711", "#0737D3"))
  p <- p + xlab("Genomic location")
  p <- p + scale_x_continuous(
    breaks = unlist(chromosome_midpoints),
    labels = paste0("chr", 1:length(chromosome_midpoints)),
    limits = c(0,nrow(plot_df)), 
    expand = c(0, 0)
  )
  p <- p + ylab("Log10 fold change")
  p <- p + scale_y_continuous(
    breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4),
    labels = c("-4", "-3", "-2", "-1", "0", "1", "2", "3", "4"),
    limits = c(min(permutation_data$log_modified_fold_change_df$count), max(permutation_data$log_modified_fold_change_df$count))
  )
  for (end in chromosome_ends) {
    p <- p + geom_vline(xintercept=end)
  }
  p <- p + geom_segment(
    x=0, 
    xend=max(unlist(chromosome_ends)), 
    y=permutation_data$log_median_modified_fold_change, 
    yend=permutation_data$log_median_modified_fold_change, 
    size=1, color="red"
  )
  pdf(paste0(plot_dir, "2.log_modified_fold_change_from_median.pdf"), width = 20)
    print(p)
  dev.off()
}


################################################################################
### 6. Bind new epithelial and original stromal counts and write as infercnv 
# inputs ###
################################################################################

if (
  !file.exists(paste0(Robject_dir, "3.final_permutation_data.Rdata"))
) {

  merged_permutation_data <- list(
    permutation_record = permutation_data$permutation_record,
    merged_df = cbind(
      permutation_data$modified_df,
      non_epithelial_df
    ),
    log_modified_fold_change_df = permutation_data$log_modified_fold_change_df,
    infercnv_metadata = infercnv_metadata
  )
  
  print("Final dimensions of modified_df = ")
  print(dim(permutation_data$modified_df))
  print("Final dimensions of non_epithelial_df = ")
  print(dim(non_epithelial_df))
  print("Final dimensions of merged_df = ")
  print(dim(merged_permutation_data$merged_df))
  print("Final dimensions of infercnv_metadata = ")
  print(dim(merged_permutation_data$infercnv_metadata))
    
  if (!file.exists(paste0(out_dir, "input_matrix.txt"))) {
    print("Writing InferCNV input matrix...")
    write.table(
      merged_permutation_data$merged_df, 
      paste0(out_dir, "input_matrix.txt"), 
      sep = "\t", quote = F, col.names=T, row.names=T)
    print("Writing InferCNV input metadata...")
    write.table(as.matrix(merged_permutation_data$infercnv_metadata), 
      paste0(out_dir, "metadata.txt"), 
      sep = "\t", 
      quote = F, 
      col.names = F, 
      row.names = F
    )
  }
  
  saveRDS(
    merged_permutation_data, paste0(Robject_dir, "3.final_permutation_data.Rdata")
  )

}

