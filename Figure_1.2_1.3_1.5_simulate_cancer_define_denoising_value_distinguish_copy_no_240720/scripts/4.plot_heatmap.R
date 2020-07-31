#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

project_name <- "thesis"
subproject_name <- "Figure_1.2_and_1.3_simulate_cancer_and_define_denoising_value"
args = commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
sim_name <- args[2]
denoise_value <- args[3]
analysis_mode <- args[4]
min_CNV_proportion <- as.numeric(args[5])
print("Minimum CNV proportion = ")
print(min_CNV_proportion)
min_CNV_length <- as.numeric(args[5])
print("Minimum CNV length = ")
print(min_CNV_length)

sample_name <- "CID4520N"
sim_name <- "sim4"
denoise_value <- "1.3"
analysis_mode <- "samples"
min_CNV_proportion <- as.numeric("0.5")
min_CNV_length <- 20
copy_no_signal <- "0_-0.0755_0.5_-0.0617_1.5_0.0368_2_0.0566_3_0.0673"
temp_signal <- unlist(
  strsplit(copy_no_signal, "_")
)
copy_no_signal <- as.numeric(temp_signal[c(FALSE, TRUE)])
names(copy_no_signal) <- temp_signal[c(TRUE, FALSE)]

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("Dataset name = ", sim_name))
print(paste0("Denoising value = ", denoise_value))
print(paste0("Min. no cells with CNV signal required for CNV calls = ", 
  min_CNV_proportion))
min_CNV_proportion <- as.numeric(min_CNV_proportion)

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(Repitools, lib.loc = lib_loc)
library(Seurat)
library(reshape2)
library(dplyr)
library(cowplot)
library(GenomicRanges)
library(ggplot2)
library(data.table)
library(Polychrome, lib.loc=lib_loc)
library(ComplexHeatmap, lib.loc=lib_loc)
library(circlize, lib.loc = lib_loc)
library(fpc, lib.loc = lib_loc)
library(naturalsort, lib.loc = lib_loc)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", 
  project_name, "/", subproject_name, "/")
ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/functions/")
results_dir <- paste0(project_dir, "results/")
seurat_dir <- paste0(project_dir, "raw_files/seurat_objects/", 
  sample_name, "/")
copy_dist_dir <- paste0(results_dir, "/infercnv/", sample_name, 
  "/distinguish_between_CNV_results/Rdata/")

sim_dir <- paste0(results_dir, "cancer_simulation/", sample_name,
  "/", sim_name, "/Rdata/")
general_sim_dir <- paste0(results_dir, "cancer_simulation/", sample_name,
  "/Rdata/")

in_path <- paste0(results_dir, "infercnv/", sample_name, "/", sim_name, "/")
in_dir <- paste0(in_path, denoise_value, "_denoising/", analysis_mode, 
  "_mode/")

Robject_dir <- paste0(in_dir, "Rdata/")
system(paste0("mkdir -p ", Robject_dir))
plot_dir <- paste0(in_dir, "plots/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(in_dir, "tables/")
system(paste0("mkdir -p ", table_dir))

input_dir <- paste0(in_path, "input_files/")

print(paste0("In directory = ", in_dir))
print(paste0("Reference directory = ", ref_dir))
print(paste0("R function directory = ", func_dir))
print(paste0("R object directory = ", Robject_dir))
print(paste0("Table directory = ", table_dir))
print(paste0("Plot directory = ", plot_dir))


################################################################################
### 0. Define functions ###
################################################################################

create_extended_vector <- dget(paste0(func_dir, 
  "create_extended_vector.R"))

fetch_chromosome_boundaries <- dget(paste0(func_dir, 
  "fetch_chromosome_boundaries.R"))

prepare_for_granges <- dget(paste0(func_dir, 
  "prepare_for_granges.R"))

count_calls <- dget(paste0(func_dir, 
  "count_calls.R"))

create_signal_plot <- dget(paste0(func_dir, 
  "create_signal_plot.R"))


################################################################################
### 1. Load InferCNV output and create heatmap and metadata dfs ###
################################################################################

if (!file.exists(paste0(Robject_dir, "/1b.initial_epithelial_metadata.Rdata"))) {
  
  # load InferCNV output:
  print("Loading InferCNV output files...")

  if (denoise_value == "no") {
    infercnv_output <- as.data.frame(t(read.table(paste0(in_dir, 
      "infercnv.observations.txt"))))
  } else {
    infercnv_output <- as.data.frame(t(read.table(paste0(in_dir, 
      "infercnv.observations.txt"))))
  }


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

# fetch chromosome boundary co-ordinates:
if (!file.exists(paste0(Robject_dir, "chromosome_data.Rdata"))) {
  chr_data <- fetch_chromosome_boundaries(epithelial_heatmap, ref_dir)
  saveRDS(chr_data, paste0(Robject_dir, "chromosome_data.Rdata"))
} else {
  chr_data <- readRDS(paste0(Robject_dir, "chromosome_data.Rdata"))
}


################################################################################
### 2. Calculate mean fold difference from CNV-neutral value per cell  ###
################################################################################

if (sim_name == "normal" | sim_name == "filtered_normal" | sim_name == "real_cancer") {
 
  # determine neutral value (based on sd denoising):
  score_table <- table(round(unlist(epithelial_heatmap), 6))
  if (denoise_value == "no") {
    neutral_value <- mean(unlist(epithelial_heatmap))
  } else {
    neutral_value <- as.numeric(
      names(score_table)[which.max(score_table)]
    )
  }
  CNV_fd <- apply(epithelial_heatmap, 1, function(x) mean(abs(neutral_value-x)))
  overall_CNV_score <- round(mean(CNV_fd), 5)
  write.table(
    overall_CNV_score,
    paste0(table_dir, "overall_mean_CNV_score.txt"),
    sep = "\t",
    quote = F,
    col.names = F,
    row.names = F
  )

  #CNV_means <- apply(epithelial_heatmap, 1, function(x) mean(abs(x)))
  CNV_fd <- data.frame(
    cell_id = names(CNV_fd),
    score = CNV_fd
  )
  write.table(
    CNV_fd,
    paste0(table_dir, "mean_CNV_score_per_cell.txt"),
    sep = "\t",
    quote = F,
    col.names = T,
    row.names = F
  )
  

} else {


  ################################################################################
  ### 3. Format CNV indices  ###
  ################################################################################
  
  # create simulated CNV annotation:
  simulated_CNV_plot_data <- readRDS(paste0(sim_dir, 
    "simulated_CNV_plot_data.Rdata"))
  log_modified_fold_change_df <- 
    simulated_CNV_plot_data$log_modified_fold_change_df
  CNV_indices <- simulated_CNV_plot_data$CNV_indices
  
  # if CNV-neutral regions not present, fill in:
  if ( !(1.0 %in% CNV_indices$multiplier) ) {
  
    print("Filling in CNV-neutral regions...")
  
    # order CNV_indices and fill in areas of neutral CNV:
    CNV_indices <- CNV_indices[order(CNV_indices$start),]
    CNV_gr <- GRanges(
      seqnames = Rle("genome"),
      ranges = IRanges(start = CNV_indices$start, end = CNV_indices$end),
      strand = Rle("*"),
      multiplier = CNV_indices$multiplier,
      log_median_modified_FC = CNV_indices$log_median_modified_FC
    )
    CNV_gaps <- gaps(CNV_gr)
    CNV_gaps$multiplier <- 1
    CNV_gaps$log_median_modified_FC <- 0
    CNV_gr_full <- c(CNV_gr, CNV_gaps)
    CNV_gr_full <- CNV_gr_full[order(end(ranges(CNV_gr_full)))]
    
    # add last CNV neutral region:
    CNV_gr_full <- c(
      CNV_gr_full,
      GRanges(
        seqnames = Rle("genome"),
        ranges = IRanges(
          start = end(ranges(CNV_gr_full))[length(CNV_gr_full)]+1, 
          end = nrow(log_modified_fold_change_df)
        ),
        strand = Rle("*"),
        multiplier = 1,
        log_median_modified_FC = 0
      )
    )
    
    # convert back to df:
    CNV_indices <- data.frame(
      start = start(ranges(CNV_gr_full)),
      end = end(ranges(CNV_gr_full)),
      multiplier = CNV_gr_full$multiplier,
      log_median_modified_FC = CNV_gr_full$log_median_modified_FC
    )
  
  }
  
  # record original CNV lengths:
  CNVs_only <- CNV_indices[CNV_indices$multiplier != 1,]
  orig_CNV_lengths <- CNVs_only$end - CNVs_only$start
  
  # only keep genes present in epithelial_heatmap:
  log_modified_fold_change_df <- log_modified_fold_change_df[
      colnames(epithelial_heatmap),
  ]
  
  # load all gene annotation:
  gene_annotation <- readRDS(paste0(general_sim_dir, 
    "/1e.gene_annotation.Rdata"))
  
  colnames(gene_annotation) <- c("gene", "chromosome", "start", "end")
  # number genes:
  gene_annotation$number <- seq_along(gene_annotation$gene)
  # keep only those in epithelial_heatmap:
  gene_annotation <- gene_annotation[gene_annotation$gene %in% 
    colnames(epithelial_heatmap),]
  # for each row in CNV_indices, designate start as either
  # 1 or end gene number of last segment containing genes + 1, count genes 
  # remaining in gene annotation within CNV_indices coordinates, and 
  # define end as start + (no. genes-1) to account for the reduced lengths 
  # of segments due to filtered out genes:
  for (n in 1:nrow(CNV_indices)) {
   
    print(n)
  
    # determine how many genes are present in segment:
    no_genes <- length(gene_annotation$gene[
      gene_annotation$number >= CNV_indices$start[n] & 
      gene_annotation$number <= CNV_indices$end[n]
    ])
    print(paste0("No. genes = ", no_genes))
  
    if (no_genes > 0) {
  
      if (n==1) {
  
        CNV_indices$start[n] <- 1
  
      } else {
  
        # find last segment containing genes:
        previous_ends <- CNV_indices$end[1:(n-1)]
        # if at least one end co-ordinate is not zero, use as
        # start co-ordinate for current segment, otherwise make
        # start co-ordinate = 1:
        if (any(previous_ends != 0)) {
          next_segment_up <- max(which(previous_ends != 0))
          CNV_indices$start[n] <- CNV_indices$end[next_segment_up] + 1
        } else {
          CNV_indices$start[n] <- 1
        }
       print(paste0("Start of segment is ", CNV_indices$start[n]))
  
      }
  
      CNV_indices$end[n] <- CNV_indices$start[n] + (no_genes-1)
      print(paste0("End of segment is ", CNV_indices$end[n]))
  
    } else {
  
      # mark segments without any genes for removal:
      CNV_indices$start[n] <- 0
      CNV_indices$end[n] <- 0
  
    }
  }
  CNV_indices <- CNV_indices[CNV_indices$start != 0,]
  
  # record length and midpoint of segments:
  CNV_indices$length <- (CNV_indices$end - CNV_indices$start)+1
  CNV_indices$midpoints <- CNV_indices$start + floor(CNV_indices$length/2)
  CNV_indices$number <- seq_along(CNV_indices$start)
  CNV_indices$ticks <- "exclude"
  CNV_indices$ticks[CNV_indices$multiplier != 1] <- "include"
  
  # if start and end chromosome information not present for CNVs, fill in:
  if ( !("start_chr" %in% colnames(CNV_indices)) ) {
    CNV_indices$start_chr <- "chr1"
    CNV_indices$end_chr <- "chr1"
    for (k in 1:length(chr_data$ends)) {
      if (k==1) {
    
        CNV_indices$start_chr[
          CNV_indices$start <= chr_data$ends[k]
        ] <- names(chr_data$ends)[k]
    
        CNV_indices$end_chr[
          CNV_indices$end <= chr_data$ends[k]
        ] <- names(chr_data$ends)[k]
    
      } else {
    
        CNV_indices$start_chr[
          CNV_indices$start <= chr_data$ends[k] & 
          CNV_indices$start > chr_data$ends[k-1]
        ] <- names(chr_data$ends)[k]
    
        CNV_indices$end_chr[
          CNV_indices$end <= unlist(chr_data$ends[k]) & 
          CNV_indices$end > unlist(chr_data$ends[k-1])
        ] <- names(chr_data$ends)[k]
    
      }
    }
  }
  
}


################################################################################
### 4. Create simulated CNV annotations  ###
################################################################################

if (sim_name != "normal" & sim_name != "filtered_normal" & sim_name != "real_cancer") {  
 
  if (!file.exists(paste0(Robject_dir, "grid_sim_plot.Rdata"))) {

    CNV_indices$length_lab <- CNV_indices$length
    CNV_indices$length_lab[CNV_indices$ticks == "include"][c(FALSE, TRUE)] <- paste0(
      "\n", CNV_indices$length_lab[CNV_indices$ticks == "include"][c(FALSE, TRUE)]
    )

    # create CNV annotation based on fold change CNV:
    p <- ggplot(log_modified_fold_change_df, 
      aes(x=number, y=count))
    p <- p + scale_x_continuous(
      limits = c(
        0,length(log_modified_fold_change_df$count)
      ), 
      expand = c(0, 0),
      breaks = CNV_indices$midpoints[CNV_indices$ticks == "include"],
      labels = CNV_indices$length_lab[CNV_indices$ticks == "include"]
    )
    p <- p + scale_y_continuous(
      breaks = c(0, 1, 3),
      limits = c(
        min(CNV_indices$multiplier), 
        max(CNV_indices$multiplier)
      ),
      labels = c("Total\nloss", "1", "3")
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
    p <- p + xlab("CNV lengths (no. genes)")
    # remove axis labels:
    p <- p + theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      text = element_text(size=30),
      axis.text.y = element_text(size=35),
      axis.title.y = element_text(
        size=35, 
        margin=margin(t = 0, r = 30, b = 0, l = 0)
      )
    )
    grid_sim_plot <- ggplotGrob(p)
    dev.off()

    saveRDS(grid_sim_plot, paste0(Robject_dir, "grid_sim_plot.Rdata"))
  
    if (!file.exists(paste0(plot_dir, "sim_CNV_annotation_plot.pdf"))) {
      pdf(paste0(plot_dir, "sim_CNV_annotation_plot.pdf"))
        grid.draw(grid_sim_plot)
      dev.off()
    }

  } else {
    grid_sim_plot <- readRDS(paste0(Robject_dir, "grid_sim_plot.Rdata"))
  }
  
}

  
################################################################################
### 5. Determine neutral regions  ###
################################################################################

if (sim_name != "normal" & sim_name != "filtered_normal" & sim_name != "real_cancer") {  
 
  # record gain or loss in each genomic segment:
  CNV_indices$type <- "neutral"
  CNV_indices$type[CNV_indices$multiplier > 1] <- "gain"
  CNV_indices$type[CNV_indices$multiplier < 1] <- "loss"
  CNV_indices$index <- seq(1, nrow(CNV_indices), 1)
  # split into list:
  CNV_indices_list <- split(CNV_indices, CNV_indices$index)

  # determine neutral value (based on sd denoising):
  if (!file.exists(paste0(Robject_dir, "neutral_value.Rdata"))) {

    score_table <- table(round(unlist(epithelial_heatmap), 6))
    if (denoise_value == "no") {
      neutral_value <- mean(unlist(epithelial_heatmap))
    } else {
      neutral_value <- as.numeric(
        names(score_table)[which.max(score_table)]
      )
    }
    saveRDS(neutral_value, paste0(Robject_dir, "neutral_value.Rdata"))

  } else {
    neutral_value <- readRDS(paste0(Robject_dir, "neutral_value.Rdata"))
  }
  
  # plot distribution of average epithelial CNV signal to determine thresholds
  # for gains and losses:
  average_epithelial <- apply(epithelial_heatmap, 2, mean)
  # create histogram of average CNV values:
  if (!file.exists(paste0(plot_dir, "epithelial_CNV_histogram.png"))) {
    average_epithelial_histogram <- hist(average_epithelial)
    dev.off()
    pdf(paste0(plot_dir, "epithelial_CNV_histogram.pdf"))
      plot(
        average_epithelial_histogram, 
        main=NA, 
        xlab = "CNV signal",
        ylab = "",
        yaxt='n',
        cex.axis=1.5,
        cex.lab=2
      )
    dev.off()
    png(paste0(plot_dir, "epithelial_CNV_histogram.png"))
      plot(
        average_epithelial_histogram, 
        main=NA, 
        xlab = "CNV signal",
        ylab = "",
        yaxt='n',
        cex.axis=1.5,
        cex.lab=2
      )
    dev.off()
  }

  # plot InferCNV score distribution:
  if (!file.exists(paste0(plot_dir, "epithelial_CNV_histogram.png"))) {
    pdf(paste0(plot_dir, "epithelial_CNV_histogram.pdf"))
      hist(unlist(epithelial_heatmap))
    dev.off()
    pdf(paste0(plot_dir, "epithelial_CNV_histogram.pdf"))
      hist(unlist(epithelial_heatmap))
    dev.off()
  }
  
}


################################################################################
### 6. Determine accuracy calls  ###
################################################################################

if (sim_name != "normal" & sim_name != "filtered_normal" & sim_name != "real_cancer") {  
 
  if (!file.exists(paste0(Robject_dir, "accuracy_calls.Rdata"))) {
    
    CNV_indices_list <- split(CNV_indices, 1:nrow(CNV_indices))
  
    CNV_accuracy_list <- lapply(CNV_indices_list, function(x) {
  
      # create df displaying per-gene CNV information:
      CNV_gene_df <- data.frame(
        gene_no = x$start:x$end,
        CNV_type = x$type,
        mean_signal = average_epithelial[x$start:x$end],
        infercnv_call = NA,
        accuracy_call = NA
      )
  
      for (j in 1:nrow(CNV_gene_df)) {
        scores <- round(epithelial_heatmap[,rownames(CNV_gene_df)[j]], 6)
        if (
          length(which(scores < neutral_value)) >= length(scores)*min_CNV_proportion
        ) {
          CNV_gene_df$infercnv_call[j] <- "loss"
        } else if (
          length(which(scores > neutral_value)) >= length(scores)/2*min_CNV_proportion
        ) {
          CNV_gene_df$infercnv_call[j] <- "gain"
        } else {
          CNV_gene_df$infercnv_call[j] <- "neutral"
        }
      }

      # determine accuracy call based on expected and called CNV status:
      CNV_gene_df$accuracy_call[
        CNV_gene_df$CNV_type != "neutral" & 
        CNV_gene_df$CNV_type == CNV_gene_df$infercnv_call
      ] <- "true_positive"
  
      CNV_gene_df$accuracy_call[
        CNV_gene_df$CNV_type != "neutral" & 
        CNV_gene_df$infercnv_call == "neutral"
      ] <- "false_negative"
  
      CNV_gene_df$accuracy_call[
        CNV_gene_df$CNV_type != "neutral" & 
        CNV_gene_df$infercnv_call != "neutral" & 
        CNV_gene_df$CNV_type != CNV_gene_df$infercnv_call
      ] <- "wrong_call"
  
      CNV_gene_df$accuracy_call[
        CNV_gene_df$CNV_type == "neutral" & 
        CNV_gene_df$CNV_type == CNV_gene_df$infercnv_call
      ] <- "true_negative"
  
      CNV_gene_df$accuracy_call[
        CNV_gene_df$CNV_type == "neutral" & 
        CNV_gene_df$CNV_type != CNV_gene_df$infercnv_call
      ] <- "false_positive"
  
      return(CNV_gene_df)
  
    })
  
    CNV_accuracy_df <- do.call("rbind", CNV_accuracy_list)
    
    saveRDS(
      CNV_accuracy_df, 
      paste0(Robject_dir, "accuracy_calls.Rdata")
    )
    
  } else {
  
    CNV_accuracy_df <- readRDS(
      paste0(Robject_dir, "accuracy_calls.Rdata"
    ))
  }
  
}

  
################################################################################
### 7. Determine accuracy calls  ###
################################################################################

if (sim_name != "normal" & sim_name != "filtered_normal" & sim_name != "real_cancer") {  
 
  if (!file.exists(paste0(table_dir, "accuracy_metrics.txt")) | 
    !file.exists(paste0(Robject_dir, "accuracy_annotation_vector.Rdata")) | 
    !file.exists(paste0(Robject_dir, "accuracy_heatmap_obj.Rdata"))
  ) {

    # find gene lengths of accuracy calls and make sure they add up:
    # create accuracy metrics df:
    accuracy_metrics <- data.frame(
      row.names = c("true_positive", "true_negative", 
        "false_positive", "false_negative", "wrong_call"),
      number = c(
        nrow(CNV_accuracy_df[CNV_accuracy_df$accuracy_call == "true_positive",]),
        nrow(CNV_accuracy_df[CNV_accuracy_df$accuracy_call == "true_negative",]),
        nrow(CNV_accuracy_df[CNV_accuracy_df$accuracy_call == "false_positive",]),
        nrow(CNV_accuracy_df[CNV_accuracy_df$accuracy_call == "false_negative",]),
        nrow(CNV_accuracy_df[CNV_accuracy_df$accuracy_call == "wrong_call",])
      )
    )
    
    print(paste0("Total number of genes = ", nrow(CNV_accuracy_df)))
    print(paste0(
      "Number of true positive genes = ", accuracy_metrics["true_positive",]
    ))
    print(paste0(
      "Number of true negative genes = ", accuracy_metrics["true_negative",]
    ))
    print(paste0(
      "Number of false positive genes = ", accuracy_metrics["false_positive",]
    ))
    print(paste0(
      "Number of false negative genes = ", accuracy_metrics["false_negative",]
    ))
    print(paste0(
      "Number of wrong calls = ", accuracy_metrics["wrong_call",]
    ))
  
    CNV_sensitivity <- round(
      accuracy_metrics["true_positive",]/(
        accuracy_metrics["true_positive",] + accuracy_metrics["false_negative",]
      ), 3
    )
    CNV_specificity <- round(
      accuracy_metrics["true_negative",]/(
        accuracy_metrics["true_negative",] + accuracy_metrics["false_positive",]
      ), 3
    )
    CNV_precision <- round(
      accuracy_metrics["true_positive",]/(
        accuracy_metrics["true_positive",] + accuracy_metrics["false_positive",]
      ), 3
    )
    CNV_F1 <- round(
      2*(
        (CNV_precision*CNV_sensitivity) / (CNV_precision+CNV_sensitivity)
      ), 3
    )
  
    print(paste0("Sensitivity is ", CNV_sensitivity))
    print(paste0("Specificity is ", CNV_specificity))
    print(paste0("Precision is ", CNV_precision))
    print(paste0("F1 score is ", CNV_F1))
  
    accuracy_metrics <- rbind(
      accuracy_metrics,
      data.frame(
        row.names = c("sensitivity", "specificity",
          "precision", "F1"), 
        number = c(CNV_sensitivity, CNV_specificity,
          CNV_precision, CNV_F1)
      )
    )

    write.table(
      accuracy_metrics,
      paste0(table_dir, "accuracy_metrics.txt"),
      sep = "\t",
      quote = F,
      col.names = T,
      row.names = T
    )


    ################################################################################
    ### 8. Annotate true/false positives/negatives ###
    ################################################################################
    
    # create annotation of true/false positives/negatives:
    accuracy_annotation_vector <- CNV_accuracy_df$accuracy_call
  
    cols <- c("#430F82", "#7CBA61", "#B488B4", 
      "#F6DC15", "#C02456")
    names(cols) <- c("true_negative", "false_negative", "true_positive",  
      "false_positive", "wrong_call")
    
    accuracy_annotation <- HeatmapAnnotation(
      accuracy = accuracy_annotation_vector,
      col = list(accuracy = cols),
      annotation_name_side = "left",
      annotation_legend_param = list(title = "", 
        labels_gp = gpar(fontsize = 12))
    )
    accuracy_annotation@name <- "accuracy"
  
    accuracy_cols <- structure(
      c("#430F82", "#B488B4", "#F6DC15", "#7CBA61", "#C02456"),
      names = c("true_positive", "true_negative", "false_positive", 
        "false_negative", "wrong_call")
    )
    accuracy_heatmap <- Heatmap(
      t(matrix(accuracy_annotation_vector)),
      name = "annotation_heatmap",
      col = accuracy_cols,
      show_heatmap_legend = FALSE
  
    )
    accuracy_heatmap@name <- "accuracy_heatmap"
  
    accuracy_heatmap_obj <- grid.grabExpr(
      draw(accuracy_heatmap, heatmap_legend_side = "left")
    )
    dev.off()
  

    saveRDS(
      accuracy_annotation_vector, 
      paste0(Robject_dir, "accuracy_annotation_vector.Rdata")
    )

    saveRDS(
      accuracy_heatmap_obj, 
      paste0(Robject_dir, "accuracy_heatmap_obj.Rdata")
    )

  } else {
  
    accuracy_metrics <- read.table(
      paste0(table_dir, "accuracy_metrics.txt"),
      sep = "\t",
      header = T,
      as.is = T
    )

    accuracy_annotation_vector <- readRDS( 
      paste0(Robject_dir, "accuracy_annotation_vector.Rdata")
    )

    accuracy_heatmap_obj <- readRDS(
      paste0(Robject_dir, "accuracy_heatmap_obj.Rdata")
    )

  }

}

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

# choose min and max InferCNV values for colour scheme:
for (r in 1:10) {
  if (!(round(min(na_less_vector), r) == round(max(na_less_vector), r))) {
    legend_lims <- c(round(min(na_less_vector), r), round(max(na_less_vector), r))
    break()
  }
}
  

################################################################################
### 9. Generate heatmap ###
################################################################################

if (!file.exists(paste0(plot_dir, "annotated_infercnv_plot.png"))) {

  if (sim_name == "normal" | sim_name == "filtered_normal" | sim_name == "real_cancer") {
    final_heatmap <- Heatmap(
      plot_object, 
      name = paste0("hm"),
      col = heatmap_cols,
      cluster_columns = F, 
      cluster_rows = F,
      show_row_names = F, 
      show_column_names = F,
      show_row_dend = F,
      heatmap_legend_param = list(
        title = "CNV\nscore", 
        at = legend_lims,
        color_bar = "continuous", 
        grid_height = unit(1.5, "cm"), 
        grid_width = unit(1.5, "cm"), 
        legend_direction = "horizontal",
        title_gp = gpar(fontsize = 18, fontface = "bold"), 
        labels_gp = gpar(fontsize = 12)
      ),
      use_raster = T, 
      raster_device = c("png")
    )
  } else {
      final_heatmap <- Heatmap(
        plot_object, 
        name = paste0("hm"), 
        col = heatmap_cols,
        cluster_columns = F, 
        cluster_rows = F,
        show_row_names = F, 
        show_column_names = F,
        show_row_dend = F,
        #bottom_annotation = accuracy_annotation,
        heatmap_legend_param = list(
          title = "CNV\nscore", 
          at = legend_lims,
          color_bar = "continuous", 
          grid_height = unit(1.5, "cm"), 
          grid_width = unit(1.5, "cm"), 
          legend_direction = "horizontal",
          title_gp = gpar(fontsize = 18, 
          fontface = "bold"), 
          labels_gp = gpar(fontsize = 12)
        ),
        use_raster = T, 
        raster_device = c("png")
      )
  }
  
  annotated_heatmap <- grid.grabExpr(
    draw(final_heatmap, gap = unit(6, "mm"), heatmap_legend_side = "left",
    annotation_legend_side = "left")
  )
  
  
  ################################################################################
  ### 10. Plot annotated heatmap ###
  ################################################################################
  
  # plot final annotated heatmap:
  png(
    paste0(plot_dir, "annotated_infercnv_plot.png"), 
    height = 13, 
    width = 20, 
    res = 300, 
    units = "in"
  )
  
  if (sim_name == "normal" | sim_name == "filtered_normal" | sim_name == "real_cancer") {
  
    grid.newpage()
    pushViewport(viewport(x = 0.07, y = 0.1, width = 0.92, height = 0.85, 
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
    pushViewport(viewport(x=x_coord + 0.085, y=0.89, width = 0.1, height = 0.1, 
      just = "right"))
      grid.text(paste0("Mean score =\n", overall_CNV_score), 
        gp=gpar(fontsize=16))
    popViewport()
  
  } else {
  
    grid.newpage()
    pushViewport(viewport(x = 0.01, y = 0.235, width = 0.99, height = 0.75, 
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
    pushViewport(viewport(x = x_coord + 0.964, y = 0.01, 
      width = 0.916, height = 0.2, just = c("right", "bottom")))
      grid.draw(grid_sim_plot)
    popViewport()
    pushViewport(viewport(x=x_coord + 0.071, y=0.86, width = 0.1, height = 0.1, 
    just = "right"))
      grid.text(
        paste0(
        "Sensitivity = ", 
        accuracy_metrics["sensitivity",]
        ), gp=gpar(fontsize=16)
      )
    popViewport()
    pushViewport(viewport(x=x_coord + 0.075, y=0.83, width = 0.1, height = 0.1, 
    just = "right"))
      grid.text(
        paste0(
        "Specificity = ", 
        accuracy_metrics["specificity",]
        ), gp=gpar(fontsize=16)
      )
    popViewport() 
    pushViewport(viewport(x=x_coord + 0.07, y=0.8, width = 0.1, height = 0.1, 
    just = "right"))
      grid.text(
        paste0(
        "Precision = ", 
        accuracy_metrics["precision",]
        ), gp=gpar(fontsize=16)
      )
    popViewport() 
    pushViewport(viewport(x=x_coord + 0.073, y=0.77, width = 0.1, height = 0.1, 
    just = "right"))
      grid.text(
        paste0(
        "F1 score = ", 
        accuracy_metrics["F1",]
        ), gp=gpar(fontsize=16)
      )
    popViewport() 
    pushViewport(viewport(x=x_coord + 0.072, y=0.74, width = 0.1, height = 0.1, 
    just = "right"))
      grid.text(
        paste0(
        "No. wrong calls = ", 
        accuracy_metrics["wrong_call",]
        ), gp=gpar(fontsize=16)
      )
    popViewport() 
  }
    
  dev.off()

  print(paste0("Heatmap created, output in ", plot_dir))

}


################################################################################
### 11. Create basic heatmap ###
################################################################################

if (!file.exists(paste0(plot_dir, "infercnv_plot.png"))) {

  final_heatmap <- Heatmap(
    plot_object, 
    name = paste0("hm"),
    col = heatmap_cols,
    cluster_columns = F, 
    cluster_rows = F,
    show_row_names = F, 
    show_column_names = F,
    show_row_dend = F,
    show_heatmap_legend = F,
    use_raster = T, 
    raster_device = c("png")
  )
  
  basic_heatmap <- grid.grabExpr(
    draw(final_heatmap, gap = unit(6, "mm"), heatmap_legend_side = "left",
    annotation_legend_side = "left")
  )
  dev.off()
  
  # generate heatmap legend:
  if (denoise_value == "2") {
    signal_ranges <- range(unlist(plot_object))
    lgd <- Legend(
      at = c(0.9 , 1, 1.1),
      col_fun = heatmap_cols, 
      title = "CNV\nsignal", 
      direction = "horizontal",
      grid_height = unit(2.5, "cm"),
      grid_width = unit(0.1, "cm"),
      labels_gp = gpar(fontsize = 22),
      title_gp = gpar(fontsize = 28, fontface = "plain")
    )
  } else {
    signal_ranges <- round(range(unlist(plot_object)), 2)
    lgd <- Legend(
      at = c(signal_ranges[1], 1, signal_ranges[2]),
      labels = c("loss", "", "gain"),
      col_fun = heatmap_cols, 
      title = "CNV signal", 
      direction = "horizontal",
      grid_height = unit(3, "cm"),
      grid_width = unit(4.5, "cm"),
      legend_height = unit(3, "cm"),
      legend_width = unit(4.5, "cm"),
      labels_gp = gpar(fontsize = 30),
      title_gp = gpar(fontsize = 32, fontface = "plain")
    )
  }
  
  # plot final basic heatmap:
  png(
    paste0(plot_dir, "infercnv_plot.png"), 
    height = 14, 
    width = 25, 
    res = 300, 
    units = "in"
  )
  
    if (sim_name == "normal" | sim_name == "filtered_normal" | sim_name == "real_cancer") {
  
      grid.newpage()
      pushViewport(viewport(x = 0.165, y = 0.15, width = 0.81, height = 0.85, 
        just = c("left", "bottom")))
        grid.draw(basic_heatmap)
        decorate_heatmap_body("hm", {
          for ( e in 1:length(chr_data$end_pos) ) {
            grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), 
              gp = gpar(lwd = 1, col = "#383838"))
            if (e==1) {
              grid.text(names(chr_data$lab_pos)[e], chr_data$lab_pos[e], 
              unit(0, "npc") + unit(-3.5, "mm"), gp=gpar(fontsize=26))
            } else if (e==21) {
              grid.text(paste0("\n", gsub("chr", "", names(chr_data$lab_pos)[e])), chr_data$lab_pos[e], 
              unit(0, "npc") + unit(-3.5, "mm"), gp=gpar(fontsize=26))
            } else {
              grid.text(gsub("chr", "", names(chr_data$lab_pos)[e]), chr_data$lab_pos[e], 
              unit(0, "npc") + unit(-3.5, "mm"), gp=gpar(fontsize=26))        
            }
          }
        })
      popViewport()

      # plot heatmap legend:
      pushViewport(viewport(x = unit(3.5, "cm"), y = unit(28, "cm"), width = unit(4.5, "cm"), 
        height = unit(4, "cm")))
        draw(lgd, x = unit(0.1, "cm"), y = unit(0.1, "cm"), just = c("left", "bottom"))
      popViewport()

  
    } else {
  
      # plot CNV heatmap:
      grid.newpage()
      pushViewport(viewport(x = 0.132, y = 0.28, width = 0.85, height = 0.67, 
        just = c("left", "bottom")))
        grid.draw(basic_heatmap)
        decorate_heatmap_body("hm", {
          for ( e in 1:length(chr_data$end_pos) ) {
            grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), 
              gp = gpar(lwd = 1, col = "#383838"))
            if (e==1) {
              grid.text(names(chr_data$lab_pos)[e], chr_data$lab_pos[e], 
              unit(0, "npc") + unit(-3.5, "mm"), gp=gpar(fontsize=30))
            } else if (e==21) {
              grid.text(paste0("\n", gsub("chr", "", names(chr_data$lab_pos)[e])), chr_data$lab_pos[e], 
              unit(0, "npc") + unit(-3.5, "mm"), gp=gpar(fontsize=30))
            } else {
              grid.text(gsub("chr", "", names(chr_data$lab_pos)[e]), chr_data$lab_pos[e], 
              unit(0, "npc") + unit(-3.5, "mm"), gp=gpar(fontsize=30))        
            }
          }
        })
      popViewport()
    
      # plot ground truth CNVs:
      pushViewport(viewport(x = 0.505, y = 0.13, width = 0.953, height = 0.2))
        grid.draw(grid_sim_plot)
      popViewport()

      # plot heatmap legend:
      pushViewport(viewport(x = unit(3.5, "cm"), y = unit(28, "cm"), width = unit(4.5, "cm"), 
        height = unit(4, "cm")))
        draw(lgd, x = unit(0.1, "cm"), y = unit(0.1, "cm"), just = c("left", "bottom"))
      popViewport()
  
    }
    
  dev.off()

  print(paste0("Heatmap created, output in ", plot_dir))

}

# print epithelial_metadata as table:
write.table(epithelial_metadata, paste0(table_dir, "epithelial_metadata.txt"), 
  sep="\t", quote=F, row.names=F, col.name=T)

print(paste0("Heatmap created, output in ", plot_dir))


#################################################################################
#### 12. Detect CNV peaks ###
#################################################################################

if (sim_name != "normal" & sim_name != "filtered_normal"  | sim_name == "real_cancer" & !file.exists(paste0(Robject_dir, "/CNV_data.Rdata"))) {
  
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

  ######
  # determine whether there is a gain or loss signal for each gene:
  CNV_by_gene <- apply(epithelial_heatmap, 2, function(x) {
    if (length(which(round(x, 6) > neutral_value)) >= 
      length(x)*min_CNV_proportion) {
      return("gain")
    } else if (length(which(round(x, 6) < neutral_value)) >= 
      length(x)*min_CNV_proportion) {
      return("loss")
    } else {
      return("neutral")
    }
  })

  # remove gain or loss CNVs in <10 gene stretches:
  run_length <- rle(CNV_by_gene)
  rle_df <- data.frame(
    row.names = c(
      colnames(epithelial_heatmap)[1],
      names(run_length$lengths)[1:(length(run_length$lengths)-1)]
    ),
    length = run_length$lengths,
    value = run_length$values
  )
  temp_CNV_df <- rle_df[rle_df$value != "neutral",]
  non_CNV <- temp_CNV_df[temp_CNV_df$length < min_CNV_length,]

  if (nrow(non_CNV) > 0) {
    for (r in 1:nrow(non_CNV)) {
      gene_ind <- which(names(CNV_by_gene) == rownames(non_CNV)[r])
      CNV_by_gene[gene_ind:(gene_ind + non_CNV$length[r] - 1)] <- 
        "neutral"
    }
  }

  # record CNV indices:
  CNV_df <- data.frame(
    index = 1:length(CNV_by_gene),
    type = CNV_by_gene
  )
  # split by type:
  split_CNV <- split(CNV_df, rleid(CNV_df$type))
  # convert into indices:
  split_indices <- lapply(split_CNV, function(x) {
    return(
      data.frame(
        start = x$index[1],
        end = x$index[nrow(x)],
        type = as.character(x$type[1])
      )
    )
  })
  detected_CNV <- do.call("rbind", split_indices)

  # fetch chromosome info:
  chr_data <- fetch_chromosome_boundaries(epithelial_heatmap, ref_dir)

  for (m in 1:nrow(detected_CNV)) {
    if (m==1) {
      heatmap_CNV <- detected_CNV$start[m]:detected_CNV$end[m]
    } else {
      heatmap_CNV <- c(heatmap_CNV,
        detected_CNV$start[m]:detected_CNV$end[m]
      )
    }   
  }  
  for (k in 1:length(chr_data$ends)) {
    if (k==1) {

      detected_CNV$start_chr[
        detected_CNV$start <= chr_data$ends[k]
      ] <- names(chr_data$ends)[k]

      detected_CNV$end_chr[
        detected_CNV$end <= chr_data$ends[k]
      ] <- names(chr_data$ends)[k]

    } else {

      detected_CNV$start_chr[
        detected_CNV$start <= chr_data$ends[k] & 
        detected_CNV$start > chr_data$ends[k-1]
      ] <- names(chr_data$ends)[k]

      detected_CNV$end_chr[
        detected_CNV$end <= chr_data$ends[k] & 
        detected_CNV$end > chr_data$ends[k-1]
      ] <- names(chr_data$ends)[k]

    }
  }

  # isolate CNVs only and join adjacent ones <20 genes apart:
  CNV_only <- detected_CNV[detected_CNV$type != "neutral",]
  CNV_only$keep <- TRUE
  for (i in 2:nrow(CNV_only)) {
    if ( (CNV_only$start[i] - CNV_only$end[i-1]) <= 20 ) {
      CNV_only$keep[i] <- FALSE
      CNV_only$keep[i-1] <- FALSE
      new_row <- data.frame(
        start = CNV_only$start[i-1],
        end = CNV_only$end[i],
        type = CNV_only$type[i],
        start_chr = CNV_only$start_chr[i-1],
        end_chr = CNV_only$end_chr[i],
        keep = TRUE
      )
      CNV_only <- rbind(CNV_only, new_row)
    }
  }
  # remove old indices of CNVs that have been concatenated, and keep column:
  CNV_only <- CNV_only[CNV_only$keep,]
  CNV_only <- subset(
    CNV_only[order(CNV_only$start),], select = -keep
  )

}


#################################################################################
#### 12. Predict copy number of CNV peaks ###
#################################################################################

if (sim_name != "normal" & sim_name != "filtered_normal"  | sim_name == "real_cancer" & !file.exists(paste0(Robject_dir, "/CNV_data.Rdata"))) {
 
  # fetch signal values for each CNV and estimate copy number for all CNVs
  # by comparing this distribution to copy number signal across all simulations:
  copy_dist <- readRDS(
  	paste0(
      copy_dist_dir, 
      "all_true_positive_CNV_signal_per_multiplier.Rdata"
  	)
  )

  # fetch signal values and estimate copy number for all CNVs:
  for (r in 1:nrow(CNV_only)) {
    
    # fetch signal values:
    CNV_signal <- area_df$average_score[(CNV_only$start[r]):(CNV_only$end[r])]
    
    # determine which copy number distribution CNV_signal belongs to 
    # (highest p-value of wilcox test):
    wilcox_p <- lapply(copy_dist, function(x) {
      return(wilcox.test(CNV_signal, x$signal)$p.value)
    })

    if (r==1) {
      wilcox_ps <- list(wilcox_p)
    } else {
      wilcox_ps[[r]] <- wilcox_p
    }

    CNV_only$copy_no[r] <- names(which.max(unlist(wilcox_p)))
    
  }

  # define text colour and midpoints for plotting:
  CNV_only$midpoint <- CNV_only$start + (floor((CNV_only$end - CNV_only$start)/2))
  CNV_only$estimate_colour[CNV_only$copy_no > 1] <- "#BF3667"
  CNV_only$estimate_colour[CNV_only$copy_no < 1] <- "#58B9DB"

  # split into gain and loss estimate dfs:
  all_peak_estimates <- CNV_only
  peak_estimates <- list(
    gain = CNV_only[CNV_only$type == "gain",],
    loss = CNV_only[CNV_only$type == "loss",]
  )

  # insert newlines where needed for labelling:
  peak_estimates <- lapply(peak_estimates, function(x) {

    newline_record <- data.frame(
      row.names = 1:3,
      in_label = rep(FALSE, 3)
    )
    x$estimate_lab <- x$copy_no
    x$estimate_lab[c(FALSE, TRUE)] <- paste0(
      "\n", x$estimate_lab[c(FALSE, TRUE)]
    )
    
    return(x)

  })

}  


#################################################################################
#### 13. Determine which calls were accurate and save ###
#################################################################################

if (sim_name != "normal" & sim_name != "filtered_normal"  | sim_name == "real_cancer" & !file.exists(paste0(Robject_dir, "/CNV_data.Rdata"))) {
 
  # annotate all genes with chromosome position:
  for (l in 1:length(chr_data$lengths)) {
    if (l==1) {
      chrs <- rep(names(chr_data$lengths)[l], chr_data$lengths[l])
    } else {
      chrs <- c(chrs, rep(names(chr_data$lengths)[l], chr_data$lengths[l]))
    }
  }
  chr_positions <- data.frame(
    gene = rownames(area_df),
    chr = chrs
  )
  chr_positions <- do.call(
    "rbind",
    lapply(
      split(chr_positions, chr_positions$chr),
      function(x) {
        x$index <- 1:nrow(x)
        return(x)
      }
    )
  )
  chr_positions <- chr_positions[naturalorder(chr_positions$chr),]

  # isolate CNV indices only and annotate chromosomal start and end positions:
  CNV_only_indices <- CNV_indices[CNV_indices$multiplier != 1,]
  CNV_only_indices <- subset(
    CNV_only_indices, 
    select = c(start, end, multiplier, start_chr, end_chr)
  )
  colnames(CNV_only_indices) <- gsub(
    "multiplier", "copy_no", colnames(CNV_only_indices)
  )

  # split CNVs over 2 chromosomes as they cannot be added to GRanges objects:
  CNV_only_indices <- prepare_for_granges(CNV_only_indices)
  peak_estimates$gain <- prepare_for_granges(peak_estimates$gain)
  peak_estimates$loss <- prepare_for_granges(peak_estimates$loss)
  
  # add known CNVs to granges object:
  known_gr <- GRanges(
	  seqnames = Rle(CNV_only_indices$start_chr),
    ranges = IRanges(
      start = CNV_only_indices$chr_start, 
      end = CNV_only_indices$chr_end
    ),
    strand = Rle("*"),
    end_chr = CNV_only_indices$end_chr,
    copy_no = CNV_only_indices$copy_no
  )

  estimated_vs_known_gains <- count_calls(peak_estimates$gain, known_gr)
  estimated_vs_known_losses <- count_calls(peak_estimates$loss, known_gr)

  saveRDS(
    list(
      gains = estimated_vs_known_gains,
      losses = estimated_vs_known_losses
    ),
    paste0(Robject_dir, "estimated_vs_known_counts.Rdata")
  )

}

#save.image(paste0(Robject_dir, "pre_signal_plots.Rdata"))

#project_name <- "thesis"
#subproject_name <- "Figure_1.2_and_1.3_simulate_cancer_and_define_denoising_value"
#sample_name <- "CID4520N"
#sim_name <- "sim2"
#denoise_value <- "1.3"
#analysis_mode <- "samples"
#
#lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
#library(Repitools, lib.loc = lib_loc)
#library(Seurat)
#library(reshape2)
#library(dplyr)
#library(cowplot)
#library(GenomicRanges)
#library(ggplot2)
#library(data.table)
#library(Polychrome, lib.loc=lib_loc)
#library(ComplexHeatmap, lib.loc=lib_loc)
#library(circlize, lib.loc = lib_loc)
#library(fpc, lib.loc = lib_loc)
#library(naturalsort, lib.loc = lib_loc)
#
#home_dir <- "/share/ScratchGeneral/jamtor/"
#project_dir <- paste0(home_dir, "projects/", 
#  project_name, "/", subproject_name, "/")
#results_dir <- paste0(project_dir, "results/")
#in_path <- paste0(results_dir, "infercnv/", sample_name, "/", sim_name, "/")
#in_dir <- paste0(in_path, denoise_value, "_denoising/", analysis_mode, 
#  "_mode/")
#
#Robject_dir <- paste0(in_dir, "Rdata/")
#load(paste0(Robject_dir, "pre_signal_plots.Rdata"))
#
#create_extended_vector <- dget(paste0(func_dir, 
#  "create_extended_vector.R"))
#
#fetch_chromosome_boundaries <- dget(paste0(func_dir, 
#  "fetch_chromosome_boundaries.R"))
#
#prepare_for_granges <- dget(paste0(func_dir, 
#  "prepare_for_granges.R"))
#
#count_calls <- dget(paste0(func_dir, 
#  "count_calls.R"))
#
#create_signal_plot <- dget(paste0(func_dir, 
#  "create_signal_plot.R"))
#

#################################################################################
#### 14. Create average signal plots ###
#################################################################################

if (sim_name != "normal" & sim_name != "filtered_normal"  | sim_name == "real_cancer" & !file.exists(paste0(Robject_dir, "/CNV_data.Rdata"))) {
  if (!file.exists(paste0(plot_dir, "signal_vs_simulated_CNV_plot_with_copy_no_estimates.png"))) {

    # generate signal plots with and without accuracy annotations:
    create_signal_plot(
      CNV_indices,
      area_df,
      log_modified_fold_change_df,
      peak_estimates,
      chr_data,
      CNV_accuracy_df,
      accuracy_annotation = TRUE,
      copy_no_estimates = FALSE,
      func_dir,
      plot_dir
    )

    create_signal_plot(
      CNV_indices,
      area_df,
      log_modified_fold_change_df,
      peak_estimates,
      chr_data,
      CNV_accuracy_df,
      accuracy_annotation = FALSE,
      copy_no_estimates = FALSE,
      func_dir,
      plot_dir
    )

    # generate signal plot with copy number estimates:
    create_signal_plot(
      CNV_indices,
      area_df,
      log_modified_fold_change_df,
      peak_estimates,
      chr_data,
      CNV_accuracy_df,
      accuracy_annotation = FALSE,
      copy_no_estimates = TRUE,
      func_dir,
      plot_dir
    )

  }


  #################################################################################
  #### 17. Save data for distinguishing between different gain/loss copy 
  # number values ###
  #################################################################################
  
  CNV_data <- list(
    CNV_indices = CNV_indices,
    average_signal = area_df$average_score,
    accuracy_annotation_vector = accuracy_annotation_vector,
    genes = names(average_epithelial)
  )
  
  saveRDS(CNV_data, paste0(Robject_dir, "/CNV_data.Rdata"))

}

# plot decoy CNV data file for snakemake:
if (sim_name == "filtered_normal") {
  CNV_data <- NA
  saveRDS(CNV_data, paste0(Robject_dir, "/CNV_data.Rdata"))
}



###################################################################################
## regenerate plot with gap for accuracy annotation:
#    p <- ggplot(area_df, aes(x=index, y=average_score, fill = type))
#    p <- p + geom_bar(stat="identity")
#    p <- p + scale_fill_manual(values = cols)
#    p <- p + scale_x_continuous(
#      limits = c(
#        0,length(log_modified_fold_change_df$count)
#      ), 
#      expand = c(0, 0),
#      breaks = CNV_indices$midpoints[CNV_indices$ticks == "include"],
#      labels = stag_lab
#    )
#    p <- p + scale_y_continuous(
#      limits = c(-0.09, 0.09),
#      sec.axis = sec_axis(
#        ~., 
#        "Copy number\nfold change", 
#        breaks = c(-0.08, 0, 0.08),
#        labels = c("Total\nloss", "1", "3")
#      )
#    )
#    p <- p + theme_cowplot(12)
#    p <- p + theme(
#      axis.title.x = element_text(size=30, margin = margin(t = 20, r = 0, b = 0, l = 0)),
#      axis.text.x = element_text(size=23, margin = margin(t = 90, r = 0, b = 0, l = 0)),
#      axis.ticks.x = element_blank(),
#      axis.text.y = element_text(size=30),
#      axis.title.y = element_text(size=30, margin = margin(t = 0, r = 30, b = 0, l = 0)),
#      axis.title.y.right = element_text(size=30, margin = margin(t = 0, r = 0, b = 0, l = 0)),
#      legend.position = "none"
#    )
#    p <- p + ylab("Mean CNV signal")
#    p <- p + xlab("CNV length (genes)")
#    for (c_end in chr_data$ends) {
#      p <- p + geom_vline(xintercept=c_end)
#    }
#    # create 0 line:
#    p <- p + geom_segment(
#      x=scaled_CNV_indices$start[1],
#      xend=scaled_CNV_indices$end[
#        nrow(scaled_CNV_indices)
#      ],
#      y=0,
#      yend=0
#    )
#    for (r in 1:nrow(scaled_CNV_indices)) {
#      # create horizontal line:
#      p <- p + geom_segment(
#        x=scaled_CNV_indices$start[r], 
#        xend=scaled_CNV_indices$end[r], 
#        y=scaled_CNV_indices$multiplier[r], 
#        yend=scaled_CNV_indices$multiplier[r], 
#        size=0.75, color="#430F82"
#      )
#    
#      # create left vertical line:
#      if (r != 1) {
#        p <- p + geom_segment(
#          x=scaled_CNV_indices$start[r], 
#          xend=scaled_CNV_indices$start[r], 
#          y=scaled_CNV_indices$multiplier[r-1], 
#          yend=scaled_CNV_indices$multiplier[r], 
#          size=0.75, color="#430F82"
#        )
#      }
#    }
#    # convert barplot to grid object:
#    signal_plot <- ggplotGrob(p)
#    dev.off()
#
#    png(
#      paste0(plot_dir, "signal_vs_simulated_CNV_plot.png"), 
#      height = 8, 
#      width = 22,
#      res = 300,
#      units = "in"
#    )   
#      grid.newpage()
#  
#        # draw signal plot:
#        pushViewport(viewport(x = 0.06, y = 0.001, width = 0.93, height = 0.9, 
#          just = c("left", "bottom")))
#          grid.draw(signal_plot_no_annot)
#        popViewport()
#
#        # draw accuracy annotation:
#        pushViewport(viewport(x = 0.14, y = 0.082, width = 0.77, height = 0.13, 
#          just = c("left", "bottom")))
#          #grid.draw(accuracy_heatmap_obj)
#          grid.rect()
#        popViewport()
#
#        # draw accuracy legend:
#        pushViewport(viewport(x = unit(1, "cm"), y = unit(2, "cm"), 
#                          width = unit(5.5, "cm"), height = unit(4.5, "cm"), 
#                          just = c("left", "bottom")))
#    
#          #grid.rect()
#          if (all_accuracy == "false_positive_true_positive_wrong_call") {
#            true_pos_false_pos_wrong_legend()
#          } else if (all_accuracy == "false_negative_false_positive_true_negative_true_positive") {
#            true_pos_neg_false_pos_neg_legend()
#          }
#        
#        popViewport()
#
#        # label chromosomes:
#        for ( e in 1:length(chr_data$lab_pos) ) {
#          pushViewport(
#            viewport(
#              x = 0.12 + chr_data$lab_pos[e]/1.315, 
#              y = 0.88, width = 0.05, height = 0.05, 
#              just = c("left", "bottom")
#            )
#          )
#            if (e==1) {
#              grid.text(
#                paste0(
#                  "\n", names(chr_data$lab_pos)[e]
#                ), gp=gpar(fontsize=22))
#            } else if (e==21) {
#              grid.text(
#                gsub("chr", "", names(chr_data$lab_pos)[e]), 
#                gp=gpar(fontsize=22)
#              )
#            } else {
#              grid.text(
#                paste0(
#                  "\n", gsub("chr", "", names(chr_data$lab_pos)[e])
#                ), gp=gpar(fontsize=22)
#              )
#            }
#          popViewport()
#        }
#      
#    dev.off()
#