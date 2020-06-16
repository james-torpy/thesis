#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

project_name <- "thesis"
subproject_name <- "Figure_2.4_min_length_and_gap"
args = commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
t_cells_included <- as.logical(args[2])
simulation_number <- args[3]
analysis_mode <- args[4]
neutral_signal_range <- unlist(strsplit(args[5], split = "_"))
nUMI_threshold <- as.numeric(args[6])
nGene_threshold <- as.numeric(args[7])
gap_or_CNV <- args[8]
CNV_type <- args[9]
min_CNV_proportion <- as.numeric(args[10])

#sample_name <- "CID4520N_cancer_sim"
#t_cells_included <- TRUE
#simulation_number <- "1"
#analysis_mode <- "samples"
#neutral_signal_range <- unlist(strsplit("0.97_1.03", split = "_"))
#nUMI_threshold <- 25000
#nGene_threshold <- 5000
#gap_or_CNV <- "gap"
#CNV_type <- "gain"
#min_CNV_proportion <- as.numeric("0.5")

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("T-cells included? ", as.character(t_cells_included)))
print(paste0("Simulation number = ", simulation_number))
print(paste0("Analysis mode = ", analysis_mode))
print("Neutral signal range = ")
print(neutral_signal_range)
print(paste0("nUMI threshold = ", nUMI_threshold))
print(paste0("nGene threshold = ", as.character(nGene_threshold)))
print(paste0("Gap or CNV = ", gap_or_CNV))
print(paste0("CNV type = ", CNV_type))
print(paste0("Min. CNV proportion required for call = ", min_CNV_proportion))

library(Seurat)
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(cowplot)
library(viridis)
library(GenomicRanges)
library(ggplot2)
 
lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(Polychrome, lib.loc=lib_loc)
library(ComplexHeatmap, lib.loc=lib_loc)
library(circlize, lib.loc = lib_loc)
#library(scales, lib.loc = lib_loc)
library(fpc, lib.loc = lib_loc)
library(naturalsort, lib.loc = lib_loc)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", 
  project_name, "/", subproject_name, "/")
ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/functions/")
results_dir <- seurat_path <- paste0(project_dir, "results/")

orig_sample_name <- gsub("_cancer_sim", "", sample_name)
seurat_dir <- paste0(project_dir, "raw_files/seurat_objects/", 
  orig_sample_name, "/")

if (gap_or_CNV == "CNV") {
  sim_dir <- paste0(results_dir, "cancer_simulation/", sample_name,
    "/", gap_or_CNV, "/", simulation_number, "/Rdata/")
} else {
  sim_dir <- paste0(results_dir, "cancer_simulation/", sample_name,
    "/", gap_or_CNV, "/", simulation_number, "/", CNV_type, "/Rdata/")
}
general_sim_dir <- paste0(results_dir, "cancer_simulation/", sample_name,
  "/Rdata/")

if (t_cells_included) {
  in_path <- paste0(results_dir, "infercnv/t_cells_included/", 
    sample_name, "/", gap_or_CNV, "/", simulation_number, "/")
} else {
  in_path <- paste0(results_dir, "infercnv/t_cells_excluded/", 
    sample_name, "/", gap_or_CNV, "/", simulation_number, "/")
}

if (gap_or_CNV == "CNV") {
  in_dir <- paste0(in_path, analysis_mode, "_mode/")
  input_dir <- paste0(in_path, "input_files/")
} else {
  in_dir <- paste0(in_path, CNV_type, "/", analysis_mode, "_mode/")  
  input_dir <- paste0(in_path, CNV_type, "/input_files/")
}

Robject_dir <- paste0(in_dir, "Rdata/")
system(paste0("mkdir -p ", Robject_dir))
plot_dir <- paste0(in_dir, "plots/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(in_dir, "tables/")
system(paste0("mkdir -p ", table_dir))

print(paste0("In directory = ", in_dir))
print(paste0("Reference directory = ", ref_dir))
print(paste0("R function directory = ", func_dir))
print(paste0("R object directory = ", Robject_dir))
print(paste0("Table directory = ", table_dir))
print(paste0("Plot directory = ", plot_dir))


################################################################################
### 0. Define functions ###
################################################################################

fetch_chromosome_boundaries <- dget(paste0(func_dir, 
  "fetch_chromosome_boundaries.R"))

create_extended_vector <- function(df, column) {
  for (j in 1:nrow(df)) {
    vector_length <- length(seq(df$start[j], df$end[j]))
    if (j==1) {
      result_vector <- rep(
        eval(parse(text=paste0("df$", column, "[j]"))), vector_length
      )
    } else {
      result_vector <- c(
        result_vector,
        rep(
          eval(parse(text=paste0("df$", column, "[j]"))), vector_length
        )
      )
    }
  }
  return(result_vector)
}


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
### 3. Format CNV indices  ###
################################################################################

# load CNV information:
CNV_data <- readRDS(paste0(sim_dir, "4.final_CNV_data.Rdata"))

# load either CNV or extended gap indices as 'feature indices':
if (gap_or_CNV == "CNV") {

  feature_indices <- CNV_data$CNV_record
  CNV_indices <- CNV_data$CNV_record

} else {

  if (CNV_type == "gain") {
    if (CNV_data$extended_gap_record$multiplier[1] == 3) {
      feature_indices <- CNV_data$extended_gap_record
      CNV_indices <- CNV_data$CNV_record
    } else if (CNV_data$extended_gap_record$multiplier[1] == 0) {
      feature_indices <- CNV_data[[2]]$extended_gap_record
      CNV_indices <- CNV_data[[2]]$CNV_record
    }
  } else if (CNV_type == "loss") {
    if (CNV_data$extended_gap_record$multiplier[1] == 0) {
      feature_indices <- CNV_data$extended_gap_record
      CNV_indices <- CNV_data$CNV_record
    } else if (CNV_data$extended_gap_record$multiplier[1] == 3) {
      feature_indices <- CNV_data[[2]]$extended_gap_record
      CNV_indices <- CNV_data[[2]]$CNV_record
    }
  }

  # subtract 100 gene CNV regions ither side of gap:
  feature_indices$start <- feature_indices$start + 100
  feature_indices$end <- feature_indices$end - 100

}

log_modified_fold_change_df <- CNV_data$log_modified_fold_change_df

# only keep genes present in epithelial_heatmap:
log_modified_fold_change_df <- log_modified_fold_change_df[
    colnames(epithelial_heatmap),
]

# record length and midpoint of segments:
feature_indices$length <- (feature_indices$end - feature_indices$start)
feature_indices$midpoints <- feature_indices$start + floor(feature_indices$length/2)
feature_indices$number <- seq_along(feature_indices$start)

# fetch chromosome boundary co-ordinates:
if (!file.exists(paste0(Robject_dir, "chromosome_data.Rdata"))) {
  chr_data <- fetch_chromosome_boundaries(epithelial_heatmap, ref_dir)
  saveRDS(chr_data, paste0(Robject_dir, "chromosome_data.Rdata"))
} else {
  chr_data <- readRDS(paste0(Robject_dir, "chromosome_data.Rdata"))
}

# if start and end chromosome information not present for CNVs, fill in:
if ( !("start_chr" %in% colnames(feature_indices)) ) {
  feature_indices$start_chr <- "chr1"
  feature_indices$end_chr <- "chr1"
  for (k in 1:length(chr_data$ends)) {
    if (k==1) {
  
      feature_indices$start_chr[
        feature_indices$start <= chr_data$ends[k]
      ] <- names(chr_data$ends)[k]
  
      feature_indices$end_chr[
        feature_indices$end <= chr_data$ends[k]
      ] <- names(chr_data$ends)[k]
  
    } else {
  
      feature_indices$start_chr[
        feature_indices$start <= chr_data$ends[k] & 
        feature_indices$start > chr_data$ends[k-1]
      ] <- names(chr_data$ends)[k]
  
      feature_indices$end_chr[
        feature_indices$end <= unlist(chr_data$ends[k]) & 
        feature_indices$end > unlist(chr_data$ends[k-1])
      ] <- names(chr_data$ends)[k]
  
    }
  }
}


################################################################################
<<<<<<< HEAD:Figure_2.3_min_length_and_gap/scripts/4b.plot_individual_heatmap.R
### 4. Create simulated CNV annotation  ###
################################################################################

p <- ggplot(log_modified_fold_change_df, 
  aes(x=number, y=count))
p <- p + scale_x_continuous(
  limits = c(
    0,length(log_modified_fold_change_df$count)
  ), 
  expand = c(0, 0),
  breaks = feature_indices$midpoints[feature_indices$ticks == "include"],
  labels = feature_indices$length[feature_indices$ticks == "include"]
)
p <- p + scale_y_continuous(
  breaks = c(-3, -1, 0, 1),
  limits = c(-3, 1)
)
for (c_end in chr_data$ends) {
  p <- p + geom_vline(xintercept=c_end)
}
for (r in 1:nrow(CNV_indices)) {
  # create horizontal line:
  p <- p + geom_segment(
    x=CNV_indices$start[r], 
    xend=CNV_indices$end[r], 
    y=CNV_indices$log_median_modified_FC[r], 
    yend=CNV_indices$log_median_modified_FC[r], 
    size=1, color="#37841f"
  )

  # create left vertical line:
  if (r != 1) {
    p <- p + geom_segment(
      x=CNV_indices$start[r], 
      xend=CNV_indices$start[r], 
      y=CNV_indices$log_median_modified_FC[r-1], 
      yend=CNV_indices$log_median_modified_FC[r], 
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
### 5. Determine neutral value  ###
=======
### 4. Determine neutral value  ###
>>>>>>> starting_fig_2.7:Figure_2.4_min_length_and_gap/scripts/4b.plot_individual_heatmap.R
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
if (!file.exists(paste0(plot_dir, "epithelial_CNV_histogram.png"))) {
  average_epithelial_histogram <- density(average_epithelial, bw="SJ")
  pdf(paste0(plot_dir, "epithelial_CNV_histogram.pdf"))
    plot(average_epithelial_histogram, main=NA, xlab = "CNV signal")
  dev.off()
  png(paste0(plot_dir, "epithelial_CNV_histogram.png"))
    plot(average_epithelial_histogram, main=NA, xlab = "CNV signal")
  dev.off()
}

score_table <- table(round(unlist(epithelial_heatmap), 6))
neutral_value <- as.numeric(
  names(score_table)[which.max(score_table)]
)


################################################################################
### 5. Determine accurately called CNVs/gaps ###
################################################################################

if (gap_or_CNV == "CNV") {
  if (
    !file.exists(paste0(table_dir, "CNV_calls.txt")) |
    !file.exists(paste0(Robject_dir, "3.final_CNV_indices.Rdata"))
  ) {

    feature_indices <- CNV_indices
    feature_indices$length <- feature_indices$end-feature_indices$start

    # for each range, fetch mean infercnv signal:
    feature_indices$mean_signal <- NA
    for (r in 1:nrow(feature_indices)) {
      signal_segment <- average_epithelial[
        feature_indices$start[r]:feature_indices$end[r]
      ]
      feature_indices$mean_signal[r] <- mean(signal_segment)
    }

    # annotate non-CNV regions:
    feature_indices$call <- "non-CNV"

    for (r in 1:nrow(feature_indices)) {

      heatmap_segment <- epithelial_heatmap[,feature_indices[r,]$start:feature_indices[r,]$end]
      #cell_means <- apply(heatmap_segment, 1, mean)

      if (feature_indices[r,]$type == "gain") {

        proportion_genes <- apply(heatmap_segment, 2, function(x) {
          if (length(which(round(x, 6) > neutral_value)) > length(x)*min_CNV_proportion) {
            return("feature_called")
          } else {
            return("feature_not_called")
          }
        })

        if (
          length(which(proportion_genes == "feature_called")) >= 
          length(proportion_genes)*min_CNV_proportion
        ) {
            feature_indices$call[r] <- "feature_called"
            print(paste0("CNV present for gain length ", feature_indices$length[r]))
          
        } else {
          feature_indices$call[r] <- "feature_not_called"
          print(paste0("CNV NOT present for gain length ", feature_indices$length[r]))
          
        }

      } else if (feature_indices[r,]$type == "loss") {

        proportion_genes <- apply(heatmap_segment, 2, function(x) {
          if (length(which(round(x, 6) < neutral_value)) >= length(x)*min_CNV_proportion) {
            return("feature_called")
          } else {
            return("feature_not_called")
          }
        })

        if (
          length(which(proportion_genes == "feature_called")) >= 
          length(proportion_genes)*min_CNV_proportion
        ) {
            feature_indices$call[r] <- "feature_called"
            print(paste0("CNV present for loss length ", feature_indices$length[r]))
          } else {
          feature_indices$call[r] <- "feature_not_called"
          print(paste0("CNV NOT present for gain length ", feature_indices$length[r]))
          
        }
       
      }

    }

    # record calls for gain and loss of each length:
    feature_calls <- feature_indices[feature_indices$call != "non-CNV",]
    feature_calls$type <- "gain"
    feature_calls$type[feature_calls$multiplier == 0] <- "loss"
    feature_calls$length <- feature_calls$end-feature_calls$start
    feature_calls <- subset(feature_calls,
      select = c(type, length, call)
    )[order(feature_calls$length),]
    feature_calls <- feature_calls[order(feature_calls$type),]
    
    write.table(
      feature_calls, 
      paste0(table_dir, "CNV_calls.txt"),
      sep = "\t",
      quote = F,
      col.names = T,
      row.names = F
    )

    saveRDS(feature_indices, paste0(Robject_dir, "3.final_CNV_indices.Rdata"))

  } else {

    feature_calls <- read.table(
      paste0(table_dir, "CNV_calls.txt"),
      sep = "\t",
      header = T,
      as.is = T
    )
    feature_indices <- readRDS(paste0(Robject_dir, "3.final_CNV_indices.Rdata"))

  }

} else if (gap_or_CNV == "gap") {

  if (
    !file.exists(paste0(table_dir, "gap_calls.txt")) |
    !file.exists(paste0(Robject_dir, "3.final_gap_indices.Rdata"))
  ) {


    # for each range, determine whether there is a min_gap_length long gene window with 
    # mean signal within neutral_signal_ranges or the range of the opposite CNV type
    # in at least 80% of genes:
    feature_indices$call <- "feature_not_called"

    for (r in 1:nrow(feature_indices)) {

      print(r)
      for (g in feature_indices[r,]$start:feature_indices[r,]$end) {

        heatmap_segment <- epithelial_heatmap[,g:(g+(min_gap_length-1))]

        if (CNV_type == "gain") {


    for (r in 1:nrow(feature_indices)) {

      heatmap_segment <- epithelial_heatmap[,feature_indices[r,]$start:feature_indices[r,]$end]
      #cell_means <- apply(heatmap_segment, 1, mean)

      if (CNV_type == "gain") {

        proportion_genes <- apply(heatmap_segment, 2, function(x) {
          if (length(which(round(x, 6) <= neutral_value)) >= length(x)*min_CNV_proportion) {
            return("feature_called")
          } else {
            return("feature_not_called")
          }
        })

        if (
          length(which(proportion_genes == "feature_called")) >= 
          length(proportion_genes)*min_CNV_proportion
        ) {
            feature_indices$call[r] <- "feature_called"
            print(paste0("CNV present for gap length ", feature_indices$length[r]))
          
        } else {
          feature_indices$call[r] <- "feature_not_called"
          print(paste0("CNV NOT present for gap length ", feature_indices$length[r]))
          
        }

      } else if (CNV_type == "loss") {

        proportion_genes <- apply(heatmap_segment, 2, function(x) {
          if (length(which(round(x, 6) >= neutral_value)) >= length(x)*min_CNV_proportion) {
            return("feature_called")
          } else {
            return("feature_not_called")
          }
        })

        if (
          length(which(proportion_genes == "feature_called")) >= 
          length(proportion_genes)*min_CNV_proportion
        ) {
            feature_indices$call[r] <- "feature_called"
            print(paste0("CNV present for gap length ", feature_indices$length[r]))
          } else {
          feature_indices$call[r] <- "feature_not_called"
          print(paste0("CNV NOT present for gap length ", feature_indices$length[r]))
          
        }
       
      }

    }

    # order feature_indices and fill in areas of no gaps:
    feature_indices <- feature_indices[order(feature_indices$start),]
    feature_gr <- GRanges(
      seqnames = Rle("genome"),
      ranges = IRanges(start = feature_indices$start, end = feature_indices$end),
      strand = Rle("*"),
      start_chr = feature_indices$start_chr,
      end_chr = feature_indices$end_chr,
      call = feature_indices$call
    )
    feature_gaps <- gaps(feature_gr)
    feature_gaps$call <- "non-gap"
    feature_gr_full <- c(feature_gr, feature_gaps)
    feature_gr_full <- feature_gr_full[order(end(ranges(feature_gr_full)))]

    # fill in chromosome information:
    for (n in 1:length(feature_gr_full)) {
      if (feature_gr_full$call[n] == "non-gap") {

        if (n==1) {
          feature_gr_full$start_chr[n] <- "chr1"
        } else {
          feature_gr_full$start_chr[n] <- feature_gr_full$end_chr[n-1]
        }

        feature_gr_full$end_chr[n] <- feature_gr_full$start_chr[n+1]

      }
    }
    
    # add last feature neutral region:
    feature_gr_full <- c(
      feature_gr_full,
      GRanges(
        seqnames = Rle("genome"),
        ranges = IRanges(
          start = end(ranges(feature_gr_full))[length(feature_gr_full)]+1, 
          end = nrow(log_modified_fold_change_df)
        ),
        strand = Rle("*"),
        start_chr = feature_gr_full$end_chr[length(feature_gr_full)],
        end_chr = "chr22",
        call = "non-gap"
      )
    )
    
    # convert back to df and calculate lengths:
    feature_indices <- data.frame(
      start = start(ranges(feature_gr_full)),
      end = end(ranges(feature_gr_full)),
      call = feature_gr_full$call,
      length = end(ranges(feature_gr_full))-start(ranges(feature_gr_full))
    )

    # record calls for gain and loss of each length:
    feature_calls <- feature_indices[feature_indices$call != "non-gap",]
    feature_calls <- feature_calls[order(feature_calls$length),]
    
    write.table(
      feature_calls, 
      paste0(table_dir, "gap_calls.txt"),
      sep = "\t",
      quote = F,
      col.names = T,
      row.names = F
    )
  
    saveRDS(feature_indices, paste0(Robject_dir, "3.final_gap_indices.Rdata"))

  } else {

    feature_calls <- read.table(
      paste0(table_dir, "gap_calls.txt"),
      sep = "\t",
      header = T,
      as.is = T
    )
    feature_indices <- readRDS(paste0(Robject_dir, "3.final_gap_indices.Rdata"))

  }

}


################################################################################
### 4. Create simulated CNV annotation  ###
################################################################################

if (gap_or_CNV == "CNV") {
  feature_indices$ticks[feature_indices$multiplier != 1] <- "include"
} else {
  feature_indices$ticks[feature_indices$call != "non-gap"] <- "include"
}

# stagger CNV/gap labels:
stag_lab <- feature_indices$length[feature_indices$ticks == "include"]
stag_lab[c(FALSE, TRUE)] <- paste0("\n", stag_lab[c(FALSE, TRUE)])
# create CNV annotation based on fold change CNV:
p <- ggplot(log_modified_fold_change_df, 
  aes(x=number, y=count))
p <- p + scale_x_continuous(
  limits = c(
    0,length(log_modified_fold_change_df$count)
  ), 
  expand = c(0, 0),
  breaks = feature_indices$midpoints[feature_indices$ticks == "include"],
  labels = stag_lab
)
p <- p + scale_y_continuous(
  breaks = c(0, 1, 2, 3),
  limits = c(
    min(CNV_indices$multiplier), 
    max(CNV_indices$multiplier)
  ),
  labels = c("0", "1", "2", "3")
)
for (c_end in chr_data$ends) {
  p <- p + geom_vline(xintercept=c_end)
}
for (r in 1:nrow(CNV_indices)) {
  # create horizontal line:
  p <- p + geom_segment(
    x=CNV_indices$start[r], 
    xend=CNV_indices$end[r], 
    y=CNV_indices$multiplier[r], 
    yend=CNV_indices$multiplier[r], 
    size=1, color="#430F82"
  )

  # create left vertical line:
  if (r != 1) {
    p <- p + geom_segment(
      x=CNV_indices$start[r], 
      xend=CNV_indices$start[r], 
      y=CNV_indices$multiplier[r-1], 
      yend=CNV_indices$multiplier[r], 
      size=1, color="#430F82"
    )
  }
}
# create 1 line:
p <- p + geom_segment(
  x=CNV_indices$start[1],
  xend=CNV_indices$end[
    nrow(CNV_indices)
  ],
  y=1,
  yend=1
)
# create axis titles:
p <- p + ylab("Copy number\nfold change")
if (gap_or_CNV == "CNV") {
  p <- p + xlab("CNV lengths")
} else {
  p <- p + xlab("Gap lengths")
}
p <- p + theme(
  axis.title.x = element_text(size=30, margin = margin(t = 20, r = 0, b = 0, l = 0)),
  axis.text.x = element_text(size=22),
  axis.text.y = element_text(size=24),
  axis.title.y = element_text(size=30, margin = margin(t = 0, r = 30, b = 0, l = 0))
)

grid_sim_plot <- ggplotGrob(p)
dev.off()

if (!file.exists(paste0(plot_dir, "sim_CNV_annotation_plot.pdf"))) {
  pdf(paste0(plot_dir, "sim_CNV_annotation_plot.pdf"),
    width = 20)
    grid.draw(grid_sim_plot)
  dev.off()
}

################################################################################
### 6. Annotate calls ###
################################################################################

if (gap_or_CNV == "CNV") {

  # change annotation labels:
  call_annot_indices <- feature_indices
  call_annot_indices$call <- gsub("feature", "CNV", call_annot_indices$call)
  call_annot_indices$call <- gsub("_", " ", call_annot_indices$call)

  # expand call lengths < 15 genes to 15 genes:
  for (r in 1:nrow(call_annot_indices)) {
    if (
      call_annot_indices$call[r] != "non-CNV" & call_annot_indices$length[r] < 15
    ) {
      call_annot_indices$end[r] <- call_annot_indices$start[r] + 14
      call_annot_indices$start[r+1] <- call_annot_indices$end[r]+1
    }
  }

  # create annotation of all calls:
<<<<<<< HEAD:Figure_2.3_min_length_and_gap/scripts/4b.plot_individual_heatmap.R
  call_annotation_vector <- factor(create_extended_vector(feature_indices, "call"))
  levels(call_annotation_vector) <- c("CNV called", "CNV not called", "gap")

  cols <- c("#430F82", "#7CBA61", "#E7E4D3")
=======
  call_annotation_vector <- factor(create_extended_vector(call_annot_indices, "call"))
  levels(call_annotation_vector) <- c("CNV called", "CNV not called", "gap")

  cols <- c("#7CBA61", "#430F82", "#E7E4D3")
>>>>>>> starting_fig_2.7:Figure_2.4_min_length_and_gap/scripts/4b.plot_individual_heatmap.R
  names(cols) <- levels(call_annotation_vector)
  
  call_annotation <- HeatmapAnnotation(
    accuracy = call_annotation_vector,
    col = list(accuracy = cols),
    annotation_name_side = "left",
    annotation_legend_param = list(title = "", 
    	labels_gp = gpar(fontsize = 12))
  )
  call_annotation@name <- "call"

  call_cols <- structure(
    c("#7CBA61", "#430F82", "#E7E4D3"),
    names = levels(call_annotation_vector)
  )
  call_heatmap <- Heatmap(
    t(matrix(call_annotation_vector)),
    name = "call_heatmap",
    col = call_cols,
    show_heatmap_legend = FALSE
  )
  call_heatmap@name <- "call_heatmap"

  call_heatmap_obj <- grid.grabExpr(
    draw(call_heatmap, heatmap_legend_side = "left")
  )
  dev.off()

} else if (gap_or_CNV == "gap") {

  # change annotation labels:
  call_annot_indices <- feature_indices
  call_annot_indices$call <- gsub("feature", "gap", call_annot_indices$call)
  call_annot_indices$call <- gsub("_", " ", call_annot_indices$call)

  # expand call lengths < 15 genes to 15 genes:
  for (r in 1:nrow(call_annot_indices)) {
    if (
      call_annot_indices$call[r] != "non-CNV" & call_annot_indices$length[r] < 15
    ) {
      call_annot_indices$end[r] <- call_annot_indices$start[r] + 14
      call_annot_indices$start[r+1] <- call_annot_indices$end[r]+1
    }
  }

  # create annotation of all calls:
  call_annotation_vector <- factor(create_extended_vector(call_annot_indices, "call"))
  levels(call_annotation_vector) <- c("gap called", "gap not called", "non-gap")

  cols <- c("#7CBA61", "#430F82", "#E7E4D3")

  names(cols) <- levels(call_annotation_vector)
  
  call_annotation <- HeatmapAnnotation(
    accuracy = call_annotation_vector,
    col = list(accuracy = cols),
    annotation_name_side = "left",
    annotation_legend_param = list(title = "", 
      labels_gp = gpar(fontsize = 12))
  )
  call_annotation@name <- "call"

  call_cols <- structure(
    c("#7CBA61", "#430F82", "#E7E4D3"),
    names = levels(call_annotation_vector)
  )
  call_heatmap <- Heatmap(
    t(matrix(call_annotation_vector)),
    name = "call_heatmap",
    col = call_cols,
    show_heatmap_legend = FALSE
  )
  call_heatmap@name <- "call_heatmap"

  call_heatmap_obj <- grid.grabExpr(
    draw(call_heatmap, heatmap_legend_side = "left")
  )
  dev.off()

}


################################################################################
### 8. Generate annotated heatmap ###
################################################################################

# prepare df for plotting:
plot_object <- epithelial_heatmap
colnames(plot_object) <- rep("la", ncol(plot_object))
plot_object <- as.matrix(plot_object)

# define heatmap colours:
na_less_vector <- unlist(plot_object)
na_less_vector <- na_less_vector[!is.na(na_less_vector)]
heatmap_cols <- colorRamp2(c(min(na_less_vector), 1, max(na_less_vector)), 
      c("#00106B", "white", "#680700"), space = "sRGB")
print("Generating final heatmap...")

# generate heatmap legend:
signal_ranges <- round(range(unlist(plot_object)), 1)
lgd <- Legend(
  at = c(signal_ranges[1], 1, signal_ranges[2]),
  col_fun = heatmap_cols, 
  title = "CNV\nscore", 
  direction = "horizontal",
  grid_height = unit(2.5, "cm"),
  grid_width = unit(0.1, "cm"),
  labels_gp = gpar(fontsize = 22),
  title_gp = gpar(fontsize = 28, fontface = "plain")
)

final_heatmap_annotated <- Heatmap(
  plot_object, name = paste0("hm"), 
  col = heatmap_cols,
  cluster_columns = F, cluster_rows = F,
  show_row_names = F, show_column_names = F,
  #column_names_gp = gpar(col = "white"),
  show_row_dend = F,
  bottom_annotation = call_annotation,
  show_heatmap_legend = F,
  use_raster = T, raster_device = c("png")
)

annotated_heatmap <- grid.grabExpr(
  draw(final_heatmap_annotated, gap = unit(6, "mm"), heatmap_legend_side = "left",
  annotation_legend_side = "left")
)

# determine where starting co-ordinates for heatmap are based upon longest cluster name
# (0.00604 units per character):
longest_cluster_name <- max(nchar(unique(as.character(epithelial_metadata$cell_type))))
x_coord <- longest_cluster_name*0.0037


################################################################################
### 9. Plot annotated heatmap ###
################################################################################

# plot final annotated heatmap:
if (gap_or_CNV == "CNV") {
  pdf(
    paste0(plot_dir, "annotated_infercnv_plot_with_CNV_calls.pdf"), 
    height = 13, width = 20
  )   
} else {
  pdf(
    paste0(plot_dir, "annotated_infercnv_plot_with_gap_calls.pdf"), 
    height = 13, width = 20
  )   
}
  
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

  # plot legend:
  pushViewport(viewport(x = unit(2, "cm"), y = unit(14.5, "cm"), width = unit(0.1, "cm"), 
    height = unit(0.4, "cm"), just = c("right", "bottom")))
    draw(lgd, x = unit(0.1, "cm"), y = unit(0.1, "cm"), just = c("left", "bottom"))
  popViewport()

  pushViewport(viewport(x = x_coord + 0.965, y = 0.05, 
    width = 0.88, height = 0.088, just = c("right", "bottom")))
    grid.draw(grid_sim_plot)
  popViewport()

  pushViewport(viewport(x=x_coord + 0.0765, y=0.89, width = 0.1, height = 0.1, 
    just = "right"))
    grid.text(paste0("No. genes = ", ncol(epithelial_heatmap)), 
      gp=gpar(fontsize=16))
  popViewport()

  if (gap_or_CNV == "CNV") {

    pushViewport(viewport(x=x_coord + 0.071, y=0.86, width = 0.1, height = 0.1, 
    just = "right"))
      grid.text(
        paste0(
          "Correct gain calls = ", 
          length(which(feature_calls$call == "feature_called" & feature_calls$type == "gain"))
        ), gp=gpar(fontsize=16)
      )
    popViewport()
    pushViewport(viewport(x=x_coord + 0.075, y=0.83, width = 0.1, height = 0.1, 
      just = "right"))
      grid.text(
        paste0(
          "Correct loss calls = ", 
          length(which(feature_calls$call == "feature_called" & feature_calls$type == "loss"))
        ), gp=gpar(fontsize=16)
      )
    popViewport()
    pushViewport(viewport(x=x_coord + 0.07, y=0.8, width = 0.1, height = 0.1, 
      just = "right"))
      grid.text(
        paste0(
          "Incorrect gain calls = ", 
          length(which(feature_calls$call == "feature_not_called" & feature_calls$type == "gain"))
        ), gp=gpar(fontsize=16)
      )
    popViewport() 
    pushViewport(viewport(x=x_coord + 0.073, y=0.77, width = 0.1, height = 0.1, 
      just = "right"))
      grid.text(
        paste0(
          "Incorrect loss calls = ", 
          length(which(feature_calls$call == "feature_not_called" & feature_calls$type == "loss"))
        ), gp=gpar(fontsize=16)
      )
    popViewport() 

  } else if (gap_or_CNV == "gap") {

    pushViewport(viewport(x=x_coord + 0.071, y=0.86, width = 0.1, height = 0.1, 
    just = "right"))
      grid.text(
        paste0(
          "Correct calls = ", 
          length(which(feature_calls$call == "feature_called"))
        ), gp=gpar(fontsize=16)
      )
    popViewport()
    pushViewport(viewport(x=x_coord + 0.075, y=0.83, width = 0.1, height = 0.1, 
      just = "right"))
      grid.text(
        paste0(
          "Incorrect calls = ", 
          length(which(feature_calls$call == "feature_not_called"))
        ), gp=gpar(fontsize=16)
      )
    popViewport()

  }
  
dev.off()


################################################################################
### 10. Generate basic heatmap ###
################################################################################

final_heatmap_basic <- Heatmap(
  plot_object, name = paste0("hm"), 
  col = heatmap_cols,
  cluster_columns = F, cluster_rows = F,
  show_row_names = F, show_column_names = F,
  #column_names_gp = gpar(col = "white"),
  show_row_dend = F,
  show_heatmap_legend = F,
  use_raster = T, raster_device = c("png")
)

annotated_heatmap <- grid.grabExpr(
  draw(final_heatmap_basic, gap = unit(6, "mm"), heatmap_legend_side = "left",
  annotation_legend_side = "left")
)
dev.off()

# determine where starting co-ordinates for heatmap are based upon longest cluster name
# (0.00604 units per character):
longest_cluster_name <- max(nchar(unique(as.character(epithelial_metadata$cell_type))))
x_coord <- longest_cluster_name*0.0037


################################################################################
### 11. Plot basic heatmap ###
################################################################################

if (gap_or_CNV == "CNV") {
  png(
    paste0(plot_dir, "infercnv_plot_with_CNV_calls.png"), 
    height = 13, width = 23, res = 300, units = "in"
  )   
} else {
  png(
    paste0(plot_dir, "infercnv_plot_with_gap_calls.png"), 
    height = 13, width = 23, res = 300, units = "in"
  )   
}
  pushViewport(viewport(x = 0.15, y = 0.3, width = 0.8, height = 0.67, 
      just = c("left", "bottom")))
    grid.draw(annotated_heatmap)
    decorate_heatmap_body("hm", {
      for ( e in 1:length(chr_data$end_pos) ) {
        grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), 
          gp = gpar(lwd = 1, col = "#383838"))
        if (e==1) {
          grid.text(names(chr_data$lab_pos)[e], chr_data$lab_pos[e], 
          unit(0, "npc") + unit(-3.5, "mm"), gp=gpar(fontsize=24))
        } else if (e==21) {
          grid.text(paste0("\n", gsub("chr", "", names(chr_data$lab_pos)[e])), chr_data$lab_pos[e], 
          unit(0, "npc") + unit(-3.5, "mm"), gp=gpar(fontsize=24))
        } else {
          grid.text(gsub("chr", "", names(chr_data$lab_pos)[e]), chr_data$lab_pos[e], 
          unit(0, "npc") + unit(-3.5, "mm"), gp=gpar(fontsize=24))
        }
        
      }
    })
#    grid.rect()
  popViewport()

  # plot legend:
  pushViewport(viewport(x = unit(2, "cm"), y = unit(16, "cm"), width = unit(5, "cm"), 
    height = unit(7, "cm"), just = c("left", "bottom")))
    #grid.rect()
    draw(lgd)
  popViewport()

  pushViewport(viewport(x = 0.082, y = 0.02, 
      width = 0.87099, height = 0.24, just = c("left", "bottom")))
    grid.draw(grid_sim_plot)
    #grid.rect()
  popViewport()
  
dev.off()


#################################################################################
#### 12. Plot average signal vs accuracy annotation ###
#################################################################################

# create area plot presenting InferCNV signal:
area_df <- data.frame(
  index = seq_along(average_epithelial),
  average_score = average_epithelial-neutral_value,
  type = "neutral",
  stringsAsFactors = F
)

# label gains and losses:
area_df$type[area_df$average_score > 0] <- "gain"
area_df$type[area_df$average_score < 0] <- "loss"

CNV_indices$type[CNV_indices$multiplier > 1] <- "gain"
CNV_indices$type[CNV_indices$multiplier == 1] <- "neutral"
CNV_indices$type[CNV_indices$multiplier < 1] <- "loss"
scaled_CNV_indices <- data.frame(
  start = CNV_indices$start,
  end = CNV_indices$end,
  type = CNV_indices$type,
  multiplier = CNV_indices$multiplier
)
# 0 = 0 alleles = -0.04
# 0.5 = 1 allele = -0.02
# 1 = 2 alleles = 0
# 1.5 = 3 alleles = 0.02
# 2 = 4 alleles = 0.04
# 3 = 6 alleles = 0.08

scaled_CNV_indices$multiplier[
  scaled_CNV_indices$multiplier == 0
] <- -0.08
scaled_CNV_indices$multiplier[
  scaled_CNV_indices$multiplier == 0.5
] <- -0.02
scaled_CNV_indices$multiplier[
  scaled_CNV_indices$multiplier == 1
] <- 0
scaled_CNV_indices$multiplier[
  scaled_CNV_indices$multiplier == 1.5
] <- 0.02
scaled_CNV_indices$multiplier[
  scaled_CNV_indices$multiplier == 2
] <- 0.04
scaled_CNV_indices$multiplier[
  scaled_CNV_indices$multiplier == 3
] <- 0.08

# expand CNVs less than 15 genes long:
scaled_CNV_indices$length <- scaled_CNV_indices$end-scaled_CNV_indices$start
scaled_CNV_indices$midpoints = scaled_CNV_indices$start + floor(scaled_CNV_indices$length/2)
scaled_CNV_indices$keep = TRUE

for (r in 1:nrow(scaled_CNV_indices)) {
  if (scaled_CNV_indices$length[r] < 40 & scaled_CNV_indices$multiplier[r] != 0) {

  	if (scaled_CNV_indices$start[r] > 20) {
  	  scaled_CNV_indices$start[r] <- scaled_CNV_indices$midpoint[r]-20
  	  scaled_CNV_indices$end[r-1] <- scaled_CNV_indices$start[r]+1
  	} else {
  	  scaled_CNV_indices$start[r] <- 1
  	  if (r != 1) {
  	  	scaled_CNV_indices$keep[1] <- FALSE
  	  }
  	}

    scaled_CNV_indices$end[r] <- scaled_CNV_indices$midpoints[r]+20
    scaled_CNV_indices$start[r+1] <- scaled_CNV_indices$end[r]+1

  }
}

# calculate midpoints for feature_indices and label those to be included:
feature_indices$midpoints = feature_indices$start + floor(feature_indices$length/2)
feature_indices$ticks = "not_include"

if (gap_or_CNV == "CNV") {
  feature_indices$ticks[feature_indices$multiplier != 1] = "include"
} else {
  feature_indices$ticks[feature_indices$call != "non-gap"] = "include"
}


# define colours:
cols <- c("#F7B7B5", "#76C1C1", "black")
area_df$type <- factor(area_df$type, levels = c("gain", "loss", "neutral"))

# plot on barplot:
p <- ggplot(area_df, aes(x=index, y=average_score, fill = type))
p <- p + geom_bar(stat="identity")
p <- p + scale_fill_manual(values = cols)
p <- p + scale_x_continuous(
  limits = c(
    0,length(log_modified_fold_change_df$count)
  ), 
  expand = c(0, 0),
  breaks = feature_indices$midpoints[feature_indices$ticks == "include"],
  labels = stag_lab
)
p <- p + scale_y_continuous(
  limits = c(-0.09, 0.09),
  sec.axis = sec_axis(
    ~., 
    "Copy number\nfold change", 
    breaks = c(-0.08, 0, 0.08),
    labels = c("Total\nloss", "1", "3")
  )
)
p <- p + theme_cowplot(12)
p <- p + theme(
  axis.title.x = element_text(size=30, margin = margin(t = 20, r = 0, b = 0, l = 0)),
  axis.text.x = element_text(size=22, margin = margin(t = 70, r = 0, b = 0, l = 0)),
  #axis.title.x = element_blank(),
  #axis.text.x = element_blank(),
  axis.text.y = element_text(size=24),
  axis.title.y = element_text(size=30, margin = margin(t = 0, r = 30, b = 0, l = 0)),
  axis.title.y.right = element_text(size=30, margin = margin(t = 0, r = 0, b = 0, l = 0)),
  legend.position = "none"
)
p <- p + ylab("Mean CNV signal")
if (gap_or_CNV == "CNV") {
  p <- p + xlab("CNV length (genes)")
} else {
  p <- p + xlab("Gap length (genes)")
}
for (c_end in chr_data$ends) {
  p <- p + geom_vline(xintercept=c_end)
}
# create 0 line:
p <- p + geom_segment(
  x=scaled_CNV_indices$start[1],
  xend=scaled_CNV_indices$end[
    nrow(scaled_CNV_indices)
  ],
  y=0,
  yend=0
)
for (r in 1:nrow(scaled_CNV_indices)) {
  # create horizontal line:
  p <- p + geom_segment(
    x=scaled_CNV_indices$start[r], 
    xend=scaled_CNV_indices$end[r], 
    y=scaled_CNV_indices$multiplier[r], 
    yend=scaled_CNV_indices$multiplier[r], 
    size=0.75, color="#430F82"
  )

  # create left vertical line:
  if (r != 1) {
    p <- p + geom_segment(
      x=scaled_CNV_indices$start[r], 
      xend=scaled_CNV_indices$start[r], 
      y=scaled_CNV_indices$multiplier[r-1], 
      yend=scaled_CNV_indices$multiplier[r], 
      size=0.75, color="#430F82"
    )
  }
}

# convert barplot to grid object:
signal_plot <- ggplotGrob(p)
dev.off()

png(
  paste0(plot_dir, "signal_plot.png"), 
  height = 8, 
  width = 20,
  res = 300,
  units = "in"
)
  grid.newpage()
    pushViewport(viewport(x = 0.027, y = 0.001, width = 0.964, height = 0.9, 
      just = c("left", "bottom")))
      grid.draw(signal_plot)
    popViewport()
dev.off()

# add accuracy annotation to barplot:
if (gap_or_CNV == "CNV") {
  png(
    paste0(plot_dir, "signal_vs_simulated_CNV_plot.png"), 
    height = 8, 
    width = 20,
    res = 300,
    units = "in"
  )   
} else {
  png(
    paste0(plot_dir, "signal_vs_simulated_gap_plot.png"), 
    height = 8, 
    width = 20,
    res = 300,
    units = "in"
  )   
}

  grid.newpage()
  pushViewport(viewport(x = 0.027, y = 0.01, width = 0.964, height = 0.89, 
    just = c("left", "bottom")))
    grid.draw(signal_plot)
  popViewport()
  pushViewport(viewport(x = 0.106, y = 0.1753, width = 0.805, height = 0.13, 
    just = c("left", "bottom")))
    grid.draw(call_heatmap_obj)
  popViewport()
  for ( e in 1:length(chr_data$lab_pos) ) {
    pushViewport(viewport(x = 0.082 + chr_data$lab_pos[e]/1.25, y = 0.86, width = 0.05, height = 0.05, 
      just = c("left", "bottom")))
      if (e==1) {
        grid.text(names(chr_data$lab_pos)[e], gp=gpar(fontsize=13, fontface = "bold"))
      } else {
        grid.text(gsub("chr", "", names(chr_data$lab_pos)[e]), gp=gpar(fontsize=13, fontface = "bold"))
      }
    popViewport()
  }

  # draw legend text:
  pushViewport(viewport(x = unit(0.695, "cm"), y = unit(4.6, "cm"), width = 0.1, height = 0.1, 
    just = c("left", "bottom")))
    grid.text("gap called", gp=gpar(fontsize=20))
  popViewport()

  pushViewport(viewport(x = unit(3, "mm"), y = unit(3.4, "cm"), width = 0.1, height = 0.1, 
    just = c("left", "bottom")))
    grid.text("gap not", gp=gpar(fontsize=20))
  popViewport()

  pushViewport(viewport(x = 0.1, y = unit(2.6, "cm"), width = 0.1, height = 0.1, 
    just = c("right", "bottom")))
    grid.text("called", gp=gpar(fontsize=20))
  popViewport()

  # draw legend squares:
  pushViewport(viewport(x = unit(0.8, "cm"), y = unit(3.37, "cm"), width = unit(2, "cm"), height = unit(2, "cm"), 
    just = c("right", "bottom")))
    grid.rect(x = 1, y = 1, width = unit(5, "mm"), height = unit(5, "mm"),
      just = c("left", "bottom"), gp=gpar(col = "#7CBA61", fill = "#7CBA61"))
  popViewport()
  pushViewport(viewport(x = unit(0.8, "cm"), y = unit(2.13, "cm"), width = unit(2, "cm"), height = unit(2, "cm"), 
    just = c("right", "bottom")))
    grid.rect(x = 1, y = 1, width = unit(5, "mm"), height = unit(5, "mm"),
      just = c("left", "bottom"), gp=gpar(col = "#430F82", fill = "#430F82"))
  popViewport()

dev.off()



#######
## mock CNV annotation for positioning:
#p <- ggplot(log_modified_fold_change_df, 
#  aes(x=number, y=count))
#p <- p + scale_x_continuous(
#  limits = c(
#    0,length(log_modified_fold_change_df$count)
#  ), 
#  expand = c(0, 0),
#  breaks = feature_indices$midpoints[feature_indices$ticks == "include"],
#  labels = stag_lab
#)
#p <- p + scale_y_continuous(
#  breaks = c(0, 1, 2, 3),
#  limits = c(
#    min(CNV_indices$multiplier), 
#    max(CNV_indices$multiplier)
#  ),
#  labels = c("0", "1", "2", "3")
#)
#for (c_end in chr_data$ends) {
#  p <- p + geom_vline(xintercept=c_end)
#}
#
## create axis titles:
#p <- p + ylab("Copy number\nfold change")
#if (gap_or_CNV == "CNV") {
#  p <- p + xlab("CNV lengths")
#} else {
#  p <- p + xlab("Gap lengths")
#}
#p <- p + theme(
#  axis.title.x = element_text(size=30, margin = margin(t = 20, r = 0, b = 0, l = 0)),
#  axis.text.x = element_text(size=22),
#  axis.text.y = element_text(size=24),
#  axis.title.y = element_text(size=30, margin = margin(t = 0, r = 30, b = 0, l = 0))
#)
#
#grid_sim_mock <- ggplotGrob(p)
#dev.off()
#######
## mock signal plot:
#p <- ggplot(area_df, aes(x=index, y=average_score, fill = type))
#p <- p + scale_x_continuous(
#  limits = c(
#    0,length(log_modified_fold_change_df$count)
#  ), 
#  expand = c(0, 0),
#  breaks = feature_indices$midpoints[feature_indices$ticks == "include"],
#  labels = stag_lab
#)
#p <- p + scale_y_continuous(
#  limits = c(-0.09, 0.09),
#  sec.axis = sec_axis(
#    ~., 
#    "Copy number\nfold change", 
#    breaks = c(-0.08, 0, 0.08),
#    labels = c("Total\nloss", "1", "3")
#  )
#)
#p <- p + theme_cowplot(12)
#p <- p + theme(
#  axis.title.x = element_text(size=30, margin = margin(t = 20, r = 0, b = 0, l = 0)),
#  axis.text.x = element_text(size=22, margin = margin(t = 70, r = 0, b = 0, l = 0)),
#  #axis.title.x = element_blank(),
#  #axis.text.x = element_blank(),
#  axis.text.y = element_text(size=24),
#  axis.title.y = element_text(size=30, margin = margin(t = 0, r = 30, b = 0, l = 0)),
#  axis.title.y.right = element_text(size=30, margin = margin(t = 0, r = 0, b = 0, l = 0)),
#  legend.position = "none"
#)
#p <- p + ylab("Mean CNV signal")
#p <- p + xlab("CNV length (genes)")
#for (c_end in chr_data$ends) {
#  p <- p + geom_vline(xintercept=c_end)
#}
#
## convert barplot to grid object:
#mock_signal_plot <- ggplotGrob(p)
#dev.off()
#######
