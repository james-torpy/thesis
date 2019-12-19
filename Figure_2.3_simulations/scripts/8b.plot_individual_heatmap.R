#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

RStudio <- FALSE

project_name <- "thesis"
subproject_name <- "Figure_2.3_simulations"
args = commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
include_t_cells <- as.logical(args[2])
downsample <- as.logical(args[3])
downsample_proportion <- args[4]
neutral_signal_range <- unlist(strsplit(args[5], split = "_"))
nUMI_threshold <- as.numeric(args[6])
nGene_threshold <- as.numeric(args[7])

sample_name <- "CID4520N_cancer_sim"
include_t_cells <- TRUE
downsample <- TRUE
downsample_proportion <- "1.0"
neutral_signal_range <- unlist(strsplit("0.97_1.03", split = "_"))
nUMI_threshold <- 25000
nGene_threshold <- 5000

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("T-cells included? ", as.character(include_t_cells)))
print(paste0("Downsampled?", as.character(downsample)))
print(paste0("Downsample proportion = ", as.character(downsample_proportion)))
print("Neutral signal range = ")
print(neutral_signal_range)
print(paste0("nUMI threshold = ", nUMI_threshold))
print(paste0("nGene threshold = ", as.character(nGene_threshold)))


library(Seurat)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(dplyr)
library(cowplot)
library(viridis)
library(GenomicRanges)
library(Polychrome)

if (RStudio) {
  
  library(ComplexHeatmap)
  library(circlize)
  library(scales)
  library(fpc)
  library(naturalsort)
  
  home_dir <- "/Users/jamestorpy/clusterHome/"
  
} else {
  
  lib_loc <- "/share/ScratchGeneral/jamtor/R/3.5dev/"
  library(ComplexHeatmap, lib.loc=lib_loc)
  library(circlize, lib.loc = lib_loc)
  library(scales, lib.loc = lib_loc)
  library(fpc, lib.loc = lib_loc)
  library(naturalsort, lib.loc = lib_loc)
  
  home_dir <- "/share/ScratchGeneral/jamtor/"
  
}

project_dir <- paste0(home_dir, "projects/", 
  project_name, "/", subproject_name, "/")
ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/functions/")
results_dir <- seurat_path <- paste0(project_dir, "results/")
seurat_dir <- paste0(project_dir, "raw_files/seurat_objects/", 
  sample_name, "/")

if (downsample) {

  sim_dir <- paste0(results_dir, "cancer_simulation/", sample_name,
    "/downsampling/Rdata/")

  if (downsample_proportion == "1.0") {
  
    if (include_t_cells) {
      in_dir <- paste0(results_dir, "infercnv/t_cells_included/", 
        sample_name, "/downsampling/")
    } else {
      in_dir <- paste0(results_dir, "infercnv/t_cells_excluded/", 
        sample_name, "/downsampling/")
    }
    
  } else {
  
    if (include_t_cells) {
      in_dir <- paste0(results_dir, "infercnv/t_cells_included/", 
        sample_name, "/downsampling/", downsample_proportion, 
        "_downsampling/")
      non_downsampled_dir <- paste0(results_dir, "infercnv/t_cells_included/", 
        sample_name, "/downsampling/")
    } else {
      in_dir <- paste0(results_dir, "infercnv/t_cells_excluded/", 
        sample_name, "/downsampling/", downsample_proportion, 
        "_downsampling/")
      non_downsampled_dir <- paste0(results_dir, "infercnv/t_cells_excluded/", 
        sample_name, "/downsampling/")
    }

  }

} else {
  sim_dir <- paste0(results_dir, "cancer_simulation/", sample_name,
    "/Rdata/")

  if (include_t_cells) {
    in_dir <- paste0(results_dir, "infercnv/t_cells_included/", sample_name, "/")
  } else {
    in_dir <- paste0(results_dir, "infercnv/t_cells_excluded/", sample_name, "/")
  }
}

Robject_dir <- paste0(in_dir, "Rdata/")
system(paste0("mkdir -p ", Robject_dir))
plot_dir <- paste0(in_dir, "plots/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(in_dir, "tables/")
system(paste0("mkdir -p ", table_dir))

input_dir <- paste0(in_dir, "input_files/")
CNV_dir <- paste0(results_dir, "seurat_objects/", sample_name, "/Rdata/")

print(paste0("In directory = ", in_dir))
print(paste0("Reference directory = ", ref_dir))
print(paste0("R function directory = ", func_dir))
print(paste0("R object directory = ", Robject_dir))
print(paste0("Table directory = ", table_dir))
print(paste0("Plot directory = ", plot_dir))


################################################################################
### 0. Define functions ###
################################################################################

gg_color_hue <- dget(paste0(func_dir, "gg_color_hue.R"))
prepare_infercnv_metadata <- dget(paste0(func_dir, 
  "prepare_infercnv_metadata.R"))

# create function to split list of coordinates into 21 gene slices
# and call as false/true postives or false/true negatives:
slice_and_assess_accuracy <- function(coordinate_list, neutral_signal_range) {

  segment_coords <- seq(coordinate_list$start, coordinate_list$end)
    
  # split into 21 gene segments:
  for (i in 1:ceiling(length(segment_coords)/21)) {
    if (i==1) {
      split_vector <- rep(i, 21)
    } else {
      split_vector <- c( split_vector, rep(i, 21) )
    }
  }

  if ( length(split_vector) > length(segment_coords) ) {
    split_vector <- split_vector[1:length(segment_coords)]
  }

  slice_coords <- split(segment_coords, split_vector)
  accuracy_calls <- lapply(slice_coords, function(y) {

    segment_slice <- data.frame(
      start = range(y)[1],
      end = range(y)[2]
    )
    segment_slice <- cbind(
      segment_slice, 
      subset(coordinate_list, select = c(multiplier, start_chr, end_chr, type))
    )
    # segment average_epithelial into CNV_indices segments and record the mean:
    epithelial_segment <- average_epithelial[segment_slice$start:segment_slice$end]
    segment_slice$mean_epithelial_signal <- mean(epithelial_segment)

    if (segment_slice$mean_epithelial_signal > neutral_signal_range[1] & segment_slice$mean_epithelial_signal < neutral_signal_range[2]) {
      segment_slice$CNV_call <- "neutral"
    } else if (segment_slice$mean_epithelial_signal > neutral_signal_range[2]) {
      segment_slice$CNV_call <- "gain"
    } else if (segment_slice$mean_epithelial_signal < neutral_signal_range[1]) {
      segment_slice$CNV_call <- "loss"
    }

    if (segment_slice$type == "neutral") {

      if (segment_slice$CNV_call == "gain" | segment_slice$CNV_call == "loss") {
        segment_slice$accuracy <-  "false_positive"
      } else if (segment_slice$CNV_call == "neutral") {
        segment_slice$accuracy <-  "true_negative"
      }

    } else if (segment_slice$type == "gain") {
        
      if (segment_slice$CNV_call == "gain") {
        segment_slice$accuracy <-  "true_positive"
      } else if (segment_slice$CNV_call == "neutral") {
        segment_slice$accuracy <-  "false_negative"
      } else if (segment_slice$CNV_call == "loss") {
        segment_slice$accuracy <-  "wrong_call"
      }
    
    } else if (segment_slice$type == "loss") {
    
      if (segment_slice$CNV_call == "loss") {
        segment_slice$accuracy <-  "true_positive"
      } else if (segment_slice$CNV_call == "neutral") {
        segment_slice$accuracy <-  "false_negative"
      } else if (segment_slice$CNV_call == "gain") {
        segment_slice$accuracy <-  "wrong_call"
      }
    
    }

    if (length(y) == 21) {
      segment_slice$slice_type <- "full"
    } else {
      segment_slice$slice_type <- "partial"
    }

    return(segment_slice)

  })
  res <- do.call("rbind", accuracy_calls)

  return(res)
}

# write function to expand coordinates to vector of accuracy calls:
create_extended_vector <- function(accuracy_df, column) {
  for (j in 1:nrow(accuracy_df)) {
    vector_length <- length(seq(accuracy_df$start[j], accuracy_df$end[j]))
    if (j==1) {
      result_vector <- rep(
        eval(parse(text=paste0("accuracy_df$", column, "[j]"))), vector_length
      )
    } else {
      result_vector <- c(
        result_vector,
        rep(
          eval(parse(text=paste0("accuracy_df$", column, "[j]"))), vector_length
        )
      )
    }
  }
  return(result_vector)
}

get_false_positives <- dget(paste0(func_dir, 
  "get_false_positives.R"))

fetch_chromosome_boundaries <- dget(paste0(func_dir, 
  "fetch_chromosome_boundaries.R"))


################################################################################
### 1. Load InferCNV output and create heatmap and metadata dfs ###
################################################################################

if (!file.exists(paste0(Robject_dir, "/1b.initial_epithelial_metadata.Rdata"))) {
  # load InferCNV output:
  print("Loading InferCNV output files...")
  infercnv_output <- as.data.frame(t(read.table(paste0(in_dir, 
  	"infercnv.12_denoised.observations.txt"))))

  # load metadata df:
  metadata_df <- read.table(paste0(input_dir, "metadata.txt"), header = F,
    sep = "\t", as.is = TRUE)
  colnames(metadata_df) <- c("cell_ids", "cell_type")
  row.names(metadata_df) <- metadata_df$cell_ids

  # determine the epithelial cells and only include these in heatmap:
  print(paste0("Number of heatmap rows before non-epithelial thrown: ", 
  	nrow(infercnv_output)))
  epithelial_ids <- metadata_df$cell_ids[grep("pithelial", metadata_df$cell_type)]
  epithelial_heatmap <- infercnv_output[rownames(infercnv_output) %in% epithelial_ids,]
  print(paste0("Number of heatmap rows after non-epithelial thrown: ", 
  	nrow(epithelial_heatmap)))

  # create epithelial_metadata df and only include epithelial cells in epithelial_heatmap:
  print("Creating epithelial metadata df...")
  epithelial_metadata <- metadata_df[rownames(epithelial_heatmap),]

  saveRDS(epithelial_heatmap, paste0(Robject_dir, "/1a.initial_epithelial_heatmap.Rdata"))
  saveRDS(epithelial_metadata, paste0(Robject_dir, "/1b.initial_epithelial_metadata.Rdata"))

} else {

	print("Loading heatmap and metadata dfs...")
	epithelial_heatmap <- readRDS(paste0(Robject_dir, "/1a.initial_epithelial_heatmap.Rdata"))
	epithelial_metadata <- readRDS(paste0(Robject_dir, "/1b.initial_epithelial_metadata.Rdata"))

}


################################################################################
### 2. Add annotation metadata ###
################################################################################

if (!file.exists(paste0(Robject_dir, 
  "2b.epithelial_metadata_with_cell_type_and_QC.Rdata"))) {
  # order epithelial metadata cell type cluster levels:
  epithelial_metadata$cell_type <- factor(
    epithelial_metadata$cell_type,
    levels = naturalsort(unique(epithelial_metadata$cell_type))
  )
  print(paste0(
    "Are epithelial_metadata rownames in the same order as epithelial_heatmap? ",
    identical(rownames(epithelial_heatmap), rownames(epithelial_metadata))
  ))

  # add nUMI and nGene data to epithelial_metadata:
  print("Adding QC metrics to epithelial metadata df...")
  count_df <- read.table(paste0(input_dir, "input_matrix.txt"), header = TRUE,
    sep = "\t", as.is = TRUE)
  count_df_rownames <- rownames(count_df)
  nUMI <- apply(count_df, 2, sum)
  nGene <- apply(count_df, 2, function(x) length(x[x!=0]))

  QC <- data.frame(
    row.names = colnames(count_df),
    nUMI = nUMI,
    nGene = nGene
  )
  QC <- QC[rownames(epithelial_metadata),]
  epithelial_metadata <- cbind(epithelial_metadata, QC)

  # plot distributions of nUMI and nGene to determine cutoff for high coverage cells:
  if (downsample_proportion == "1.0") {
    if (!file.exists(paste0(plot_dir, "nUMI_density_plot.png"))) {
      
      nUMI_density_plot <- density(QC$nUMI, bw="SJ")
      pdf(paste0(plot_dir, "nUMI_density_plot.pdf"))
        plot(nUMI_density_plot, main=NA, xlab = "nUMI")
      dev.off()
      png(paste0(plot_dir, "nUMI_density_plot.png"))
        plot(nUMI_density_plot, main=NA, xlab = "nUMI")
      dev.off()

      log_nUMI_density_plot <- density(log10(QC$nUMI), bw="SJ")
      pdf(paste0(plot_dir, "log_nUMI_density_plot.pdf"))
        plot(log_nUMI_density_plot, main=NA, xlab = "log10 nUMI")
      dev.off()
      png(paste0(plot_dir, "log_nUMI_density_plot.png"))
        plot(log_nUMI_density_plot, main=NA, xlab = "log10 nUMI")
      dev.off()

      nGene_density_plot <- density(QC$nGene, bw="SJ")
      pdf(paste0(plot_dir, "nGene_density_plot.pdf"))
        plot(nGene_density_plot, main=NA, xlab = "nGene")
      dev.off()
      png(paste0(plot_dir, "nGene_density_plot.png"))
        plot(nGene_density_plot, main=NA, xlab = "nGene")
      dev.off()

      log_nGene_density_plot <- density(log10(QC$nGene), bw="SJ")
      pdf(paste0(plot_dir, "log_nGene_density_plot.pdf"))
        plot(log_nGene_density_plot, main=NA, xlab = "log10 nGene")
      dev.off()
      png(paste0(plot_dir, "log_nGene_density_plot.png"))
        plot(log_nGene_density_plot, main=NA, xlab = "log10 nGene")
      dev.off()

    }
  }

  print(paste0(
    "Are epithelial_metadata rownames still in the same order as epithelial_heatmap?? ",
    identical(rownames(epithelial_heatmap), rownames(epithelial_metadata))
  ))

  saveRDS(
    epithelial_heatmap, paste0(Robject_dir, 
    "2a.epithelial_heatmap_with_cell_type_and_QC.Rdata")
  )
  saveRDS(
    epithelial_metadata, paste0(Robject_dir, 
    "2b.epithelial_metadata_with_cell_type_and_QC.Rdata")
  )
} else {
   epithelial_heatmap <- readRDS(
    paste0(Robject_dir, 
    "2a.epithelial_heatmap_with_cell_type_and_QC.Rdata")
  )
  epithelial_metadata <- readRDS(
    paste0(Robject_dir, 
    "2b.epithelial_metadata_with_cell_type_and_QC.Rdata")
  )
  print(paste0(
    "Are epithelial_metadata rownames still in the same order as epithelial_heatmap?? ",
    identical(rownames(epithelial_heatmap), rownames(epithelial_metadata))
  ))
}


################################################################################
### 3. Create QC annotations  ###
################################################################################

# create QC annotations:
nUMI_annotation <- rowAnnotation(
  correlation_annotation = anno_barplot(
    epithelial_metadata$nUMI,
    gp = gpar(
      col = "#D8B72E", 
      width = unit(4, "cm")
    ), 
    border = FALSE, 
    which = "row", 
    axis = F
  )
)
nUMI_annotation@name <- "nUMI"
nGene_annotation <- rowAnnotation(
  correlation_annotation = anno_barplot(
    epithelial_metadata$nGene, name = "nGene",
    gp = gpar(
      col = "#9ECAE1", 
      width = unit(4, "cm")
    ), 
    border = FALSE, 
    which = "row", 
    axis = F
  )
)
nGene_annotation@name <- "nGene"


################################################################################
### 4. Create simulated CNV annotation  ###
################################################################################

# create simulated CNV annotation:
simulated_CNV_plot_data <- readRDS(paste0(sim_dir, 
  "simulated_CNV_plot_data.Rdata"))
log_adjusted_fold_change_df <- 
  simulated_CNV_plot_data$log_adjusted_fold_change_df
CNV_indices <- simulated_CNV_plot_data$CNV_indices
colnames(CNV_indices) <- gsub("chromosome", "chr", colnames(CNV_indices))
zero_count_value <- simulated_CNV_plot_data$zero_count_value

# record original CNV lengths:
CNVs_only <- CNV_indices[CNV_indices$multiplier != 1,]
orig_CNV_lengths <- CNVs_only$end - CNVs_only$start

# only keep genes present in epithelial_heatmap:
log_adjusted_fold_change_df <- log_adjusted_fold_change_df[
    colnames(epithelial_heatmap),
]
log_adjusted_fold_change_df <- log_adjusted_fold_change_df[
    colnames(epithelial_heatmap),
]

# load all gene annotation:
gene_annotation <- readRDS(paste0(sim_dir, 
  "/2b.gene_annotation.Rdata"))

colnames(gene_annotation) <- c("gene", "chromosome", "start", "end")
# number genes:
gene_annotation$number <- seq_along(gene_annotation$gene)
# keep only those in epithelial_heatmap:
gene_annotation <- gene_annotation[gene_annotation$gene %in% 
  colnames(epithelial_heatmap),]
# for each row in CNV_indices, designate start as either
# 1 or gene number of last segment + 1, count genes remaining in gene 
# annotation within CNV_indices coordinates, and define end as 
# start + (no. genes-1) to account for the reduced lengths of segments
# due to filtered out genes:
for (n in 1:nrow(CNV_indices)) {
 
  print(n)

  no_genes <- length(gene_annotation$gene[
    gene_annotation$number >= CNV_indices$start[n] & 
    gene_annotation$number <= CNV_indices$end[n]
  ])

  if (no_genes > 0) {

  	if (n==1) {
  	  CNV_indices$start[n] <- 1
    } else {
  	  CNV_indices$start[n] <- CNV_indices$end[n-1] + 1
      print(CNV_indices$start[n])
    }
    CNV_indices$end[n] <- CNV_indices$start[n] + (no_genes-1)
    print(CNV_indices$end[n])

  } else {

    # if no genes in coordinate window, make end = start and
    # remove after loop:
    if (n==1) {
      CNV_indices$start[n] <- 1
    } else {
      CNV_indices$start[n] <- CNV_indices$end[n-1]
      print(CNV_indices$start[n])
    }
    CNV_indices$end[n] <- CNV_indices$start[n]
    print(CNV_indices$end[n])

  }

}
CNV_indices <- CNV_indices[CNV_indices$start != CNV_indices$end,]

# record length and midpoint of segments:
CNV_indices$length <- CNV_indices$end - CNV_indices$start
CNV_indices$midpoints <- CNV_indices$start + floor(CNV_indices$length/2)
CNV_indices$number <- seq_along(CNV_indices$start)
CNV_indices$ticks <- "exclude"

if (downsample) {
  CNV_indices$ticks[CNV_indices$multiplier != 1] <- "include"
} else {
  CNV_indices$ticks[
    c(9, 11, 13, 15, 17, 19, 21, 24, 26, 28, 30, 32, 34, 36, 38, 41, 43, 
    45, 47, 49, 51, 53, 55, 57, 59)
  ] <- "include"
}

# fetch chromosome boundary co-ordinates:
if (!file.exists(paste0(Robject_dir, "chromosome_data.Rdata"))) {
  chr_data <- fetch_chromosome_boundaries(epithelial_heatmap, ref_dir)
  saveRDS(chr_data, paste0(Robject_dir, "chromosome_data.Rdata"))
} else {
  chr_data <- readRDS(paste0(Robject_dir, "chromosome_data.Rdata"))
}

p <- ggplot(log_adjusted_fold_change_df, 
  aes(x=number, y=count))
p <- p + scale_x_continuous(
  limits = c(
    0,length(log_adjusted_fold_change_df$count)
  ), 
  expand = c(0, 0),
  breaks = CNV_indices$midpoints[CNV_indices$ticks == "include"],
  labels = CNV_indices$length[CNV_indices$ticks == "include"]
)
p <- p + scale_y_continuous(
  breaks = c(zero_count_value, -1, 0, 1),
  limits = c(zero_count_value, 1)
)
for (c_end in chr_data$ends) {
  p <- p + geom_vline(xintercept=c_end)
}
for (r in 1:nrow(CNV_indices)) {
  print(r)
  # create horizontal line:
  p <- p + geom_segment(
    x=CNV_indices$start[r], 
    xend=CNV_indices$end[r], 
    y=CNV_indices$segment_log_median_FC[r], 
    yend=CNV_indices$segment_log_median_FC[r], 
    size=1, color="#37841f"
  )

  # create left vertical line:

  if (r != 1) {
    p <- p + geom_segment(
      x=CNV_indices$start[r], 
      xend=CNV_indices$start[r], 
      y=CNV_indices$segment_log_median_FC[r-1], 
      yend=CNV_indices$segment_log_median_FC[r], 
      size=1, color="#37841f"
    )
  }
}
# create 0 line:
p <- p + geom_segment(
  x=CNV_indices$start[1],
  xend=CNV_indices$end[
    nrow(CNV_indices)
  ],
  y=0,
  yend=0
)
# remove axis labels:
p <- p + theme(
  axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank()
)
grid_sim_plot <- ggplotGrob(p)
dev.off()

if (!file.exists(paste0(plot_dir, "sim_CNV_plot.pdf"))) {
  
  pdf(paste0(plot_dir, "sim_CNV_plot.pdf"))
    grid.draw(grid_sim_plot)
  dev.off()

}


################################################################################
### 5. Determine neutral regions  ###
################################################################################

# record gain or loss in each genomic segment:
CNV_indices$type <- "neutral"
CNV_indices$type[CNV_indices$multiplier > 1] <- "gain"
CNV_indices$type[CNV_indices$multiplier < 1] <- "loss"
CNV_indices$index <- seq(1, nrow(CNV_indices), 1)
# split into list:
CNV_indices_list <- split(CNV_indices, CNV_indices$index)

# plot distribution of average epithelial CNV signal to determine thresholds
# for gains and losses:
average_epithelial <- apply(epithelial_heatmap, 2, mean)
# create density plot of average CNV values:
if (!file.exists(paste0(plot_dir, "epithelial_CNV_density_plot.png"))) {
  average_epithelial_density_plot <- density(average_epithelial, bw="SJ")
  pdf(paste0(plot_dir, "epithelial_CNV_density_plot.pdf"))
    plot(average_epithelial_density_plot, main=NA, xlab = "CNV signal")
  dev.off()
  png(paste0(plot_dir, "epithelial_CNV_density_plot.png"))
    plot(average_epithelial_density_plot, main=NA, xlab = "CNV signal")
  dev.off()
}
# determine range of mean epithelial signal across all neutral regions:
neutral_indices <- CNV_indices[CNV_indices$type == "neutral",]

for (r in 1:nrow(neutral_indices)) {
  print(r)
  if (r==1) {
    neutral_signal <- average_epithelial[
      neutral_indices$start[r]:neutral_indices$end[r]
    ]
  } else {
    neutral_signal <- c(
      neutral_signal,
      average_epithelial[
        neutral_indices$start[r]:neutral_indices$end[r]
      ]
    )
  }
}

# create density plot of neutral signal:
if (!file.exists(paste0(plot_dir, "filtered_neutral_density_plot.png"))) {
  neutral_signal_density_plot <- density(neutral_signal, bw="SJ")
  pdf(paste0(plot_dir, "neutral_signal_density_plot.pdf"))
    plot(neutral_signal_density_plot, main=NA, xlab = "Neutral signal")
  dev.off()
  png(paste0(plot_dir, "neutral_signal_density_plot.png"))
    plot(neutral_signal_density_plot, main=NA, xlab = "Neutral signal")
  dev.off()
  
  # define neutral signal as neutral_signal_range[1] < x < neutral_signal_range[2]:
  filtered_neutral_density_plot <- density(
    neutral_signal[
      neutral_signal > neutral_signal_range[1] & neutral_signal < neutral_signal_range[2]
    ],
    bw="SJ",
  )
  pdf(paste0(plot_dir, "filtered_neutral_density_plot.pdf"))
    plot(filtered_neutral_density_plot, main=NA, xlab = "Neutral signal")
  dev.off()
  png(paste0(plot_dir, "filtered_neutral_density_plot.png"))
    plot(filtered_neutral_density_plot, main=NA, xlab = "Neutral signal")
  dev.off()
}


################################################################################
### 6. Determine true positives and false negatives  ###
################################################################################

if (!file.exists(
  paste0(Robject_dir, "1.true_positives_and_false_negatives.Rdata")
)) {
  
  CNV_regions <- CNV_indices[CNV_indices$type != "neutral",]
  CNV_regions_list <- split(CNV_regions, 1:nrow(CNV_regions))
  
  # make a call whether gain, neutral or loss using InferCNV signal from the mean 
  # of each 21 gene region (determined from the appearance of 21 gene long signal
  # in InferCNV plot) and call as true/false positives/negatives:
  true_positive_and_false_negative_list <- lapply(
    CNV_regions_list, slice_and_assess_accuracy, neutral_signal_range
  )
  true_positive_and_false_negative_df <- do.call(
    "rbind", true_positive_and_false_negative_list
  )
  true_positive_and_false_negative_df <- subset(
    true_positive_and_false_negative_df, 
    select = c(
      start, end, start_chr, end_chr, type, CNV_call, accuracy, 
      slice_type
    )
  )
  colnames(true_positive_and_false_negative_df) <- c(
    "start", "end", "start_chr", "end_chr", "type", "CNV_call", "accuracy", 
    "slice_type"
  )
  
  saveRDS(true_positive_and_false_negative_df, paste0(
    Robject_dir, "1.true_positives_and_false_negatives.Rdata"
  ))
  
} else {

  true_positive_and_false_negative_df <- readRDS(paste0(
    Robject_dir, "1.true_positives_and_false_negatives.Rdata"
  ))
}


################################################################################
### 7. Determine false positives and true negatives  ###
################################################################################

# find genes with positive signal in the region and, if 20 genes around it also 
# have positive signal call as false positives:
neutral_regions <- CNV_indices[CNV_indices$type == "neutral",]
neutral_regions_list <- split(neutral_regions, 1:nrow(neutral_regions))

# determine length of neutral_regions:
neutral_regions_type_vector <- create_extended_vector(neutral_regions, "type")

if (!file.exists(
  paste0(Robject_dir, "2.all_accuracy_calls.Rdata")
)) {
 
  # fetch chromosome coordinates for genes present in epithelial_heatmap:
  gene_annotation$new_number <- seq_along(gene_annotation$gene)
  
  # determine co-ordinates of false gains and losses:
  false_positives_and_true_negatives <- lapply(
    neutral_regions_list, get_false_positives, average_epithelial, 
    min_CNV_length = 21, gene_annotation
  )

  # combine into data frame and order:
  false_positive_and_true_negative_df <- do.call(
    "rbind", false_positives_and_true_negatives
  )

  # split false_positives into 21 gene segments:
  false_positive_unsplit_df <- false_positive_and_true_negative_df[
    false_positive_and_true_negative_df$accuracy == "false_positive",
  ]
  if (nrow(false_positive_unsplit_df) > 0) {
    false_positive_unsplit_df$multiplier <- 1
    false_positive_list <- split(false_positive_unsplit_df, 1:nrow(false_positive_unsplit_df))
    false_positive_df <- do.call(
      "rbind", lapply(false_positive_list, slice_and_assess_accuracy, neutral_signal_range)
    )
    false_positive_df <- subset(false_positive_df, 
      select = c(start, end, start_chr, end_chr, type, CNV_call, accuracy, 
      slice_type))
  }
  
  # split true negatives into 21 gene segments:
  true_negative_unsplit_df <- false_positive_and_true_negative_df[
    false_positive_and_true_negative_df$accuracy == "true_negative",
  ]
  true_negative_unsplit_df$multiplier <- 1
  true_negative_list <- split(true_negative_unsplit_df, 1:nrow(true_negative_unsplit_df))
  true_negative_df <- do.call(
    "rbind", lapply(true_negative_list, slice_and_assess_accuracy, neutral_signal_range)
  )
  true_negative_df <- subset(true_negative_df, 
  select = c(start, end, start_chr, end_chr, type, CNV_call, accuracy, 
    slice_type))

  # as these are known to be true negatives, adjust those falsely called 
  # true positives:
  true_negative_df$accuracy <- "true_negative"
  true_negative_df$CNV_call <- "neutral"

  # merge with other accuracy calls and order:
  if (exists("false_positive_df")) {
    all_accuracy_calls <- rbind(
      rbind(
        true_positive_and_false_negative_df, 
        false_positive_df
      ), 
      true_negative_df
    )
  } else {
    all_accuracy_calls <- rbind(
      true_positive_and_false_negative_df, 
      true_negative_df
    )
  }
  
  all_accuracy_calls <- all_accuracy_calls[order(all_accuracy_calls$start),]

  # find gene lengths of accuracy calls and make sure they add up:
  all_length <- length(create_extended_vector(all_accuracy_calls, "accuracy"))
  true_positive_length <- length(
    create_extended_vector(
      all_accuracy_calls[all_accuracy_calls$accuracy == "true_positive",], 
      "accuracy"
    )
  )
  if (nrow(all_accuracy_calls[all_accuracy_calls$accuracy == "false_positive",]) > 0) {
    false_positive_length <- length(
      create_extended_vector(
        all_accuracy_calls[all_accuracy_calls$accuracy == "false_positive",], 
        "accuracy"
      )
    )
  } else {
    false_positive_length <- 0
  }
  true_negative_length <- length(
    create_extended_vector(
      all_accuracy_calls[all_accuracy_calls$accuracy == "true_negative",], 
      "accuracy"
    )
  )
  false_negative_length <- length(
    create_extended_vector(
      all_accuracy_calls[all_accuracy_calls$accuracy == "false_negative",], 
      "accuracy"
    )
  )

  print(paste0("Total number of genes = ", all_length))
  print(paste0("Number of true positive genes = ", true_positive_length))
  print(paste0("Number of false positive genes = ", false_positive_length))
  print(paste0("Number of true negative genes = ", true_negative_length))
  print(paste0("Number of false negative genes = ", false_negative_length))
  print(paste("Does the total no. genes add up?", 
    all_length == true_positive_length + false_positive_length + 
    true_negative_length + false_negative_length))

  # remove partial segments for counting accuracy calls:
  full_accuracy_calls <- all_accuracy_calls[
    all_accuracy_calls$slice_type == "full",
  ]
  # find numbers of accuracy calls and make sure they add up:
  all_number <- nrow(full_accuracy_calls)
  accuracy_metrics <- table(full_accuracy_calls$accuracy)
  if (!("false_positive" %in% names(accuracy_metrics))) {
    accuracy_metrics[4] <- 0
    names(accuracy_metrics)[4] <- "false_positive"
  }

  print(paste0("Total number of calls = ", all_number))
  print(paste0("Number of true positives = ", accuracy_metrics["true_positive"]))
  print(paste0("Number of false positives = ", accuracy_metrics["false_positive"]))
  print(paste0("Number of true negatives = ", accuracy_metrics["true_negative"]))
  print(paste0("Number of false negatives = ", accuracy_metrics["false_negative"]))
  print(paste("Does the total no. genes add up?", 
    all_number == accuracy_metrics["true_positive"] + 
    accuracy_metrics["false_positive"] + 
    accuracy_metrics["true_negative"] + 
    accuracy_metrics["false_negative"]))

  CNV_sensitivity <- round(
    accuracy_metrics["true_positive"]/(
      accuracy_metrics["true_positive"] + accuracy_metrics["false_negative"]
    ), 3
  )
  CNV_specificity <- round(
    accuracy_metrics["true_negative"]/(
      accuracy_metrics["true_negative"] + accuracy_metrics["false_positive"]
    ), 3
  )

  print(paste0("Sensitivity is ", CNV_sensitivity))
  print(paste0("Specificity is ", CNV_specificity))

  accuracy_metrics <- data.frame(
    type = names(accuracy_metrics),
    number = as.vector(accuracy_metrics)
  )

  accuracy_metrics <- rbind(
    accuracy_metrics,
    data.frame(
      type = c("sensitivity", "specificity"), 
      number = c(CNV_sensitivity, CNV_specificity)
    )
  )

  saveRDS(all_accuracy_calls,
    paste0(Robject_dir, "2.all_accuracy_calls.Rdata")
  )
  saveRDS(accuracy_metrics,
    paste0(Robject_dir, "accuracy_metrics.Rdata")
  )

} else {

  all_accuracy_calls <- readRDS(
    paste0(Robject_dir, "2.all_accuracy_calls.Rdata")
  )
  accuracy_metrics <- readRDS(
    paste0(Robject_dir, "accuracy_metrics.Rdata")
  )

}


################################################################################
### 8. Annotate true/false positives/negatives ###
################################################################################

# create annotation of true/false positives/negatives:
accuracy_annotation_vector <- create_extended_vector(
  all_accuracy_calls, "accuracy"
)
cols <- c("#7CBA61", "#9DC9DC", "#B488B4", 
  "#F6DC15")
names(cols) <- c("true_negative", "false_negative", "true_positive",  
  "false_positive")

accuracy_annotation <- HeatmapAnnotation(
  accuracy = accuracy_annotation_vector,
  col = list(accuracy = cols),
  annotation_name_side = "left",
  annotation_legend_param = list(title = "", 
  	labels_gp = gpar(fontsize = 12))
)
accuracy_annotation@name <- "accuracy"


################################################################################
### 9. Calculate correlation with average of non-downsampled heatmap  ###
################################################################################

if (downsample & downsample_proportion != "1.0") {

  # load InferCNV output:
  print("Loading non-downsampled InferCNV heatmap...")
  non_downsampled_output <- as.data.frame(t(read.table(paste0(non_downsampled_dir, 
    "infercnv.12_denoised.observations.txt"))))

  # load metadata df:
  non_downsampled_metadata <- read.table(paste0(non_downsampled_dir, 
    "input_files/metadata.txt"), header = F, sep = "\t", as.is = TRUE)
  colnames(non_downsampled_metadata) <- c("cell_ids", "cell_type")
  row.names(non_downsampled_metadata) <- non_downsampled_metadata$cell_ids

  # ensure only epithelial cells in non-downsampled data:
  non_downsampled_ids <- non_downsampled_metadata$cell_ids[
    grep("pithelial", non_downsampled_metadata$cell_type)
  ]
  non_downsampled_heatmap <- non_downsampled_output[
    rownames(non_downsampled_output) %in% non_downsampled_ids,
  ]

  # calculate mean of downsampled and non-downsampled data and make them
  # the same length:
  non_downsampled_mean_CNV <- apply(non_downsampled_heatmap, 2, mean)
  downsampled_mean_CNV <- apply(epithelial_heatmap, 2, mean)
  non_downsampled_mean_CNV <- non_downsampled_mean_CNV[
    names(downsampled_mean_CNV)
  ]

  # calculate correlation between non-downsampled and downsampled CNVs:
  correlation_with_original <- cor.test(
    as.numeric(downsampled_mean_CNV), 
    as.numeric(non_downsampled_mean_CNV), 
    method = "pearson"
  )
  cor_result <- data.frame(
    R_squared = correlation_with_original$estimate, 
    p_val = correlation_with_original$p.value
  )

  # add to accuracy metrics
  final_accuracy_metrics <- rbind(
    accuracy_metrics,
    data.frame(
      type = c("mean_UMI", "no_genes", "pearson_R_squared", "pearson_p_val"),
      number = c(
      	round(mean(epithelial_metadata$nUMI), 0),
        length(downsampled_mean_CNV), 
        cor_result$R_squared,
        cor_result$p_val
      )
    )
  )
  final_accuracy_metrics$number <- round(final_accuracy_metrics$number, 3)

} else {
  colnames(accuracy_metrics) <- c("type", "number")
  final_accuracy_metrics <- rbind(
    accuracy_metrics,
    data.frame(
      type = "mean_UMI",
      number = round(mean(epithelial_metadata$nUMI), 0)
    )
  )
}

write.table(
  final_accuracy_metrics,
  paste0(table_dir, "accuracy_metrics.txt"),
  sep = "\t",
  row.names = F,
  col.names = F,
  quote = F
)


################################################################################
### 10. Generate heatmap ###
################################################################################

# prepare df for plotting:
plot_object <- epithelial_heatmap
colnames(plot_object) <- rep("la", ncol(plot_object))

# define heatmap colours:
na_less_vector <- unlist(plot_object)
na_less_vector <- na_less_vector[!is.na(na_less_vector)]
heatmap_cols <- colorRamp2(c(min(na_less_vector), 1, max(na_less_vector)), 
      c("#00106B", "white", "#680700"), space = "sRGB")
print("Generating final heatmap...")

# create main CNV heatmap:
final_heatmap <- Heatmap(
  plot_object, name = paste0("hm"), 
  col = heatmap_cols,
  cluster_columns = F, cluster_rows = F,
  show_row_names = F, show_column_names = F,
  column_names_gp = gpar(col = "white"),
  show_row_dend = F,
  bottom_annotation = accuracy_annotation,
  heatmap_legend_param = list(title = "CNV\nscore", 
    at = c(round(min(na_less_vector), 1), 1, round(max(na_less_vector), 1)),
    color_bar = "continuous", grid_height = unit(1.5, "cm"), 
    grid_width = unit(1.5, "cm"), legend_direction = "horizontal",
    title_gp = gpar(fontsize = 18, fontface = "bold"), 
    labels_gp = gpar(fontsize = 12)),
  use_raster = T, raster_device = c("png")
)

ht_list <- final_heatmap + nUMI_annotation + nGene_annotation

annotated_heatmap <- grid.grabExpr(
  draw(ht_list, gap = unit(6, "mm"), heatmap_legend_side = "left",
  annotation_legend_side = "left")
)

# determine where starting co-ordinates for heatmap are based upon longest cluster name
# (0.00604 units per character):
longest_cluster_name <- max(nchar(unique(as.character(epithelial_metadata$cell_type))))
x_coord <- longest_cluster_name*0.0037


################################################################################
### 11. Plot heatmap ###
################################################################################

# plot final annotated heatmap:
pdf(paste0(plot_dir, "infercnv_plot_with_accuracy_metrics.pdf"), height = 13, width = 18)   
  
  grid.newpage()
  pushViewport(viewport(x = 0.005, y = 0.16, width = 0.99, height = 0.78, 
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

  pushViewport(viewport(x = x_coord + 0.9112, y = 0.05, 
    width = 0.823, height = 0.088, just = c("right", "bottom")))
    grid.draw(grid_sim_plot)
  popViewport()

  pushViewport(viewport(x=x_coord + 0.9155, y=0.1, width = 0.1, height = 0.1, 
    just = "bottom"))
    grid.text("nUMI", gp=gpar(fontsize=16), rot=65)
  popViewport()
  pushViewport(viewport(x=x_coord + 0.94, y=0.1, width = 0.1, height = 0.1, 
    just = "bottom"))
    grid.text("nGene", gp=gpar(fontsize=16), rot=65)
  popViewport()

  pushViewport(viewport(x=x_coord + 0.0765, y=0.86, width = 0.1, height = 0.1, 
    just = "right"))
    grid.text(paste0("No. genes = ", ncol(epithelial_heatmap)), 
      gp=gpar(fontsize=16))
  popViewport()
  pushViewport(viewport(x=x_coord + 0.071, y=0.83, width = 0.1, height = 0.1, 
    just = "right"))
    grid.text(
      paste0(
      "mean UMI = ", 
      final_accuracy_metrics$number[final_accuracy_metrics$type == "mean_UMI"]
      ), gp=gpar(fontsize=16)
    )
  popViewport()
  pushViewport(viewport(x=x_coord + 0.075, y=0.8, width = 0.1, height = 0.1, 
    just = "right"))
    grid.text(
      paste0(
      "Sensitivity = ", 
      final_accuracy_metrics$number[final_accuracy_metrics$type == "sensitivity"]
      ), gp=gpar(fontsize=16)
    )
  popViewport()
  pushViewport(viewport(x=x_coord + 0.07, y=0.75, width = 0.1, height = 0.1, 
    just = "right"))
    grid.text(
      paste0(
      "Specificity = ", 
      final_accuracy_metrics$number[final_accuracy_metrics$type == "specificity"]
      ), gp=gpar(fontsize=16)
    )
  popViewport() 

  if (exists("cor_result")) {
    pushViewport(viewport(x=x_coord + 0.072, y=0.71, width = 0.1, height = 0.1, 
      just = "right"))
      grid.text(paste0("Correlation with\noriginal = ", round(cor_result$R_squared, 3)), 
        gp=gpar(fontsize=16))
    popViewport()
    pushViewport(viewport(x=x_coord + 0.055, y=0.66, width = 0.1, height = 0.1, 
      just = "right"))
      grid.text(paste0("p.val = ", round(cor_result$p_val, 3)), 
        gp=gpar(fontsize=16))
    popViewport()
  }
  
dev.off()

print(paste0("Heatmap created, output in ", plot_dir))

# print epithelial_metadata as table:
write.table(epithelial_metadata, paste0(table_dir, "epithelial_metadata.txt"), 
  sep="\t", quote=F, row.names=F, col.name=T)

print(paste0("Heatmap created, output in ", plot_dir))

#convert pdf to png:
system(paste0("for p in ", plot_dir, "*.pdf; do echo $p; f=$(basename $p); echo $f; ",
              "new=$(echo $f | sed 's/.pdf/.png/'); echo $new; ", 
              "convert -density 150 ", plot_dir, "$f -quality 90 ", plot_dir, "$new; done"))



################################################################################
### development ###
################################################################################

#false_positive_segments <- list(
#  gain = false_positive_segments$gain,
#  loss = data.frame(
#    start = c(5400, 5560, 5800),
#    end = c(5420, 5580, 5820),
#    type = rep("neutral", 3),
#    CNV_call = c("loss", "neutral", "loss"),
#    accuracy = c("false_positive", "true_negative", "false_positive"),
#    slice_type = rep("full", 3)
#  )
#)
#
#false_positive_grs <- list(
#  gain = false_positive_grs$gain,
#  loss = GRanges(
#    seqnames = rep("all", 3),
#    ranges = IRanges(
#      start = c(5400, 5560, 5800),
#      end = c(5420, 5580, 5820)
#    ),
#    strand = Rle(rep("*", nrow(temp_coords))),
#    accuracy = c("false_positive", "true_negative", "false_positive"),
#    CNV_call = c("loss", "neutral", "loss"),
#    slice_type = rep("full", 3) 
#  )
#)
#
#o=1
#lapply(neutral_regions_list, function(neutral_region) {
#  print(o)
#  # set default accuracy and CNV_calls to be changed if necessary:
#  neutral_region$accuracy <- "true negative"
#  neutral_region$CNV_call <- "neutral"
#  neutral_region$slice_type <- "full"
#  neutral_region$slice_type[
#    neutral_region$end-neutral_region$start+1 < 21
#  ] <- "partial"
#  
#  # subset neutral_region:
#  neutral_region <- subset(neutral_region, select = c(start, end, type, 
#    CNV_call, accuracy, slice_type))
#
#  # fetch the mean signal at each gene:
#  mean_segment_signal <- average_epithelial[
#    neutral_region$start:neutral_region$end
#  ]
#  print(paste0(length(mean_segment_signal), " genes in segment"))
#
#  index_types <- c("gain", "loss")
#
#  for (i in 1:length(index_types)) {
#
#    # determine genes with a CNV signal:
#    if (index_types[i] == "gain") {
#      signal_indices <- which(mean_segment_signal > neutral_signal_range[2])
#      print(paste0("Length of gain genes = ", length(signal_indices)))
#    } else if (index_types[i] == "loss") {
#      signal_indices <- which(mean_segment_signal < neutral_signal-range[1])
#      print(paste0("Length of loss genes = ", length(signal_indices)))
#    }
#  }
#  o <<- o+1
#  return(neutral_region)
#})


