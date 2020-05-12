#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

project_name <- "thesis"
subproject_name <- "Figure_2.2_and_2.3_define_denoising_value"
args = commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
sim_name <- args[2]
denoise_value <- args[3]
analysis_mode <- args[4]
min_CNV_proportion <- as.numeric(args[5])

#sample_name <- "CID4520N"
#sim_name <- "sim3"
#denoise_value <- "no"
#analysis_mode <- "samples"
#min_CNV_proportion <- as.numeric("0.5")

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("Dataset name = ", sim_name))
print(paste0("Denoising value = ", denoise_value))
print(paste0("Min. no cells with CNV signal required for CNV calls = ", 
  min_CNV_proportion))
min_CNV_proportion <- as.numeric(min_CNV_proportion)

library(Seurat)
library(reshape2)
library(dplyr)
library(cowplot)
library(GenomicRanges)
library(ggplot2)
library(data.table)
  
lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
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
results_dir <- seurat_path <- paste0(project_dir, "results/")
seurat_dir <- paste0(project_dir, "raw_files/seurat_objects/", 
  sample_name, "/")

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

 # write legend functions:
true_pos_neg_false_pos_neg_legend <- function() {
  # plot legend text:
  pushViewport(viewport(x = 0.48, y = 0.87, 
                        width = unit(1, "cm"), height = unit(0.5, "cm"), 
                        just = c("left")))
      grid.text("true positive", gp=gpar(fontsize=18))
  popViewport()
  pushViewport(viewport(x = 0.5, y = 0.655, 
                        width = unit(1, "cm"), height = unit(0.5, "cm"), 
                        just = c("left")))
    grid.text("true negative", gp=gpar(fontsize=18))
  popViewport()
  pushViewport(viewport(x = 0.5, y = 0.435, 
                        width = unit(1, "cm"), height = unit(0.5, "cm"), 
                        just = c("left")))
    grid.text("false positive", gp=gpar(fontsize=18))
  popViewport()
  pushViewport(viewport(x = 0.52, y = 0.215, 
                        width = unit(1, "cm"), height = unit(0.5, "cm"), 
                        just = c("left")))
    grid.text("false negative", gp=gpar(fontsize=18))
  popViewport()
  
  # plot legend squares:
  pushViewport(viewport(x = 0, y = 0.81, 
                        width = unit(1, "cm"), height = unit(0.5, "cm"), 
                        just = c("left")))
    grid.rect(width = unit(5, "mm"), height = unit(5, "mm"),
            just = c("left", "bottom"), gp=gpar(col = "#430F82", fill = "#430F82"))
  popViewport()
  pushViewport(viewport(x = 0, y = 0.595, 
                        width = unit(1, "cm"), height = unit(0.5, "cm"), 
                        just = c("left")))
    grid.rect(width = unit(5, "mm"), height = unit(5, "mm"),
              just = c("left", "bottom"), gp=gpar(col = "#B488B4", fill = "#B488B4"))
  popViewport()
  pushViewport(viewport(x = 0, y = 0.38, 
                        width = unit(1, "cm"), height = unit(0.5, "cm"), 
                        just = c("left")))
    grid.rect(width = unit(5, "mm"), height = unit(5, "mm"),
            just = c("left", "bottom"), gp=gpar(col = "#F6DC15", fill = "#F6DC15"))
  popViewport()
  pushViewport(viewport(x = 0, y = 0.165, 
                        width = unit(1, "cm"), height = unit(0.5, "cm"), 
                        just = c("left")))
    grid.rect(width = unit(5, "mm"), height = unit(5, "mm"),
            just = c("left", "bottom"), gp=gpar(col = "#7CBA61", fill = "#7CBA61"))
  popViewport()
}
true_pos_false_pos_wrong_legend <- function() {
  # plot legend text:
  pushViewport(viewport(x = 0.48, y = 0.87, 
                        width = unit(1, "cm"), height = unit(0.5, "cm"), 
                        just = c("left")))
      grid.text("true positive", gp=gpar(fontsize=18))
  popViewport()
  pushViewport(viewport(x = 0.5, y = 0.655, 
                        width = unit(1, "cm"), height = unit(0.5, "cm"), 
                        just = c("left")))
    grid.text("false positive", gp=gpar(fontsize=18))
  popViewport()
  pushViewport(viewport(x = 0.43, y = 0.435, 
                        width = unit(1, "cm"), height = unit(0.5, "cm"), 
                        just = c("left")))
    grid.text("wrong call", gp=gpar(fontsize=18))
  popViewport()
  
  # plot legend squares:
  pushViewport(viewport(x = 0, y = 0.81, 
                        width = unit(1, "cm"), height = unit(0.5, "cm"), 
                        just = c("left")))
    grid.rect(width = unit(5, "mm"), height = unit(5, "mm"),
            just = c("left", "bottom"), gp=gpar(col = "#430F82", fill = "#430F82"))
  popViewport()
  pushViewport(viewport(x = 0, y = 0.595, 
                        width = unit(1, "cm"), height = unit(0.5, "cm"), 
                        just = c("left")))
    grid.rect(width = unit(5, "mm"), height = unit(5, "mm"),
              just = c("left", "bottom"), gp=gpar(col = "#F6DC15", fill = "#F6DC15"))
  popViewport()
  pushViewport(viewport(x = 0, y = 0.38, 
                        width = unit(1, "cm"), height = unit(0.5, "cm"), 
                        just = c("left")))
    grid.rect(width = unit(5, "mm"), height = unit(5, "mm"),
            just = c("left", "bottom"), gp=gpar(col = "#C02456", fill = "#C02456"))
  popViewport()
}


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
 
if (sim_name == "normal" | sim_name == "filtered_normal") {
 
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
  
  
  ################################################################################
  ### 4. Create simulated CNV annotations  ###
  ################################################################################
  
  if (!file.exists(paste0(Robject_dir, "grid_sim_plot.Rdata"))) {

    # create CNV annotation based on fold change CNV:
    p <- ggplot(log_modified_fold_change_df, 
      aes(x=number, y=count))
    p <- p + scale_x_continuous(
      limits = c(
        0,length(log_modified_fold_change_df$count)
      ), 
      expand = c(0, 0),
      breaks = CNV_indices$midpoints[CNV_indices$ticks == "include"],
      labels = CNV_indices$length[CNV_indices$ticks == "include"]
    )
    p <- p + scale_y_continuous(
      breaks = c(0, 0.5, 1, 2, 3),
      limits = c(
        min(CNV_indices$multiplier), 
        max(CNV_indices$multiplier)
      ),
      labels = c("0", "0.5", "1", "2", "3")
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
    # create y-axis title:
    p <- p + ylab("Simulated copy number fold change")
    # remove axis labels:
    p <- p + theme(
      axis.title.x=element_blank(),
      axis.text.x = element_text(size=24),
      text = element_text(size=24),
      axis.text.y = element_text(size=20)
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
    pdf(paste0(plot_dir, "epithelial_CNV_histogram.pdf"))
      plot(average_epithelial_histogram, main=NA, xlab = "CNV signal")
    dev.off()
    png(paste0(plot_dir, "epithelial_CNV_histogram.png"))
      plot(average_epithelial_histogram, main=NA, xlab = "CNV signal")
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
  
  
  ################################################################################
  ### 6. Determine accuracy calls  ###
  ################################################################################

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
  
  
  ################################################################################
  ### 7. Determine accuracy calls  ###
  ################################################################################
  
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
  
    #  pdf(paste0(plot_dir, "accuracy_heatmap.pdf"), height = 2, width = 20)
    #  grid.newpage()
    #    pushViewport(viewport(x = 0.07, y = 0.1, width = 0.92, height = 0.85, 
    #      just = c("left", "bottom")))
    #      grid.draw(accuracy_heatmap_obj)
    #    popViewport()
    #  dev.off()

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
# determine where starting co-ordinates for heatmap are based upon longest cluster name
# (0.00604 units per character):
longest_cluster_name <- max(nchar(unique(as.character(epithelial_metadata$cell_type))))
x_coord <- longest_cluster_name*0.0037
  

################################################################################
### 9. Generate heatmap ###
################################################################################

if (!file.exists(paste0(plot_dir, "annotated_infercnv_plot.png"))) {

  if (sim_name == "normal" | sim_name == "filtered_normal") {
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
        bottom_annotation = accuracy_annotation,
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
  
  if (sim_name == "normal" | sim_name == "filtered_normal") {
  
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
      title = "CNV\nscore", 
      direction = "horizontal",
      grid_height = unit(2.5, "cm"),
      grid_width = unit(0.1, "cm"),
      labels_gp = gpar(fontsize = 22),
      title_gp = gpar(fontsize = 28, fontface = "plain")
    )
  } else {
    signal_ranges <- round(range(unlist(plot_object)), 2)
    lgd <- Legend(
      at = c(
        round(signal_ranges[1], 1), 
        1, 
        round(signal_ranges[2], 1)
      ),
      col_fun = heatmap_cols, 
      title = "CNV\nscore", 
      direction = "horizontal",
      grid_height = unit(2.5, "cm"),
      grid_width = unit(0.1, "cm"),
      labels_gp = gpar(fontsize = 22),
      title_gp = gpar(fontsize = 28, fontface = "plain")
    )
  }
  
  # plot final basic heatmap:
  png(
    paste0(plot_dir, "infercnv_plot.png"), 
    height = 12, 
    width = 23, 
    res = 300, 
    units = "in"
  )
  
    if (sim_name == "normal" | sim_name == "filtered_normal") {
  
      grid.newpage()
      pushViewport(viewport(x = 0.15, y = 0.15, width = 0.8, height = 0.8, 
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
      # plot legend:
      pushViewport(viewport(x = unit(2, "cm"), y = unit(14.5, "cm"), width = unit(0.1, "cm"), 
        height = unit(0.4, "cm"), just = c("right", "bottom")))
        draw(lgd, x = unit(0.1, "cm"), y = unit(0.1, "cm"), just = c("left", "bottom"))
      popViewport()
  
    } else {
  
      grid.newpage()
      pushViewport(viewport(x = 0.035, y = 0.22, width = 0.95, height = 0.7, 
        just = c("left", "bottom")))
        grid.draw(basic_heatmap)
        decorate_heatmap_body("hm", {
          for ( e in 1:length(chr_data$end_pos) ) {
            grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), 
              gp = gpar(lwd = 1, col = "#383838"))
            grid.text(gsub("chr", "", names(chr_data$lab_pos)[e]), chr_data$lab_pos[e], 
              unit(0, "npc") + unit(-3.5, "mm"), gp=gpar(fontsize=24))
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
      popViewport()
    
      pushViewport(viewport(x = x_coord + 0.945, y = 0.001, 
        width = 0.973, height = 0.2, just = c("right", "bottom")))
        grid.draw(grid_sim_plot)
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
#### 12. Plot average signal vs accuracy annotation ###
#################################################################################

if (sim_name != "normal" & 
  sim_name != "filtered_normal" & 
	!file.exists(paste0(Robject_dir, "/CNV_data.Rdata"))) {
  
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

  if (!file.exists(paste0(plot_dir, "signal_vs_simulated_CNV_plot.png"))) {

    # define colours:
    cols <- c("#F7B7B5", "#76C1C1", "black")
    
    scaled_CNV_indices <- data.frame(
      start = CNV_indices$start,
      end = CNV_indices$end,
      type = CNV_indices$type,
      multiplier = CNV_indices$multiplier
    )
  
    # 0 = 0 alleles = -0.08
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
  
    # stagger CNV labels for plotting:
    stag_lab <- CNV_indices$length[CNV_indices$ticks == "include"]
    stag_lab[c(FALSE, TRUE)] <- paste0("\n", stag_lab[c(FALSE, TRUE)])
    # plot on barplot:
    area_df$type <- factor(area_df$type, levels = c("gain", "loss", "neutral"))
    p <- ggplot(area_df, aes(x=index, y=average_score, fill = type))
    p <- p + geom_bar(stat="identity")
    p <- p + scale_fill_manual(values = cols)
    p <- p + scale_x_continuous(
      limits = c(
        0,length(log_modified_fold_change_df$count)
      ), 
      expand = c(0, 0),
      breaks = CNV_indices$midpoints[CNV_indices$ticks == "include"],
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
      axis.title.x = element_text(size=25, margin = margin(t = 20, r = 0, b = 0, l = 0)),
      axis.text.x = element_text(size=18, margin = margin(t = 70, r = 0, b = 0, l = 0)),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size=24),
      axis.title.y = element_text(size=25, margin = margin(t = 0, r = 30, b = 0, l = 0)),
      axis.title.y.right = element_text(size=25, margin = margin(t = 0, r = 0, b = 0, l = 0)),
      legend.position = "none"
    )
    p <- p + ylab("Mean CNV signal")
    p <- p + xlab("CNV length (genes)")
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
  
    # determine accuracy calls present:
    all_accuracy <- unique(CNV_accuracy_df$accuracy_call)
    all_accuracy <- paste0(all_accuracy[order(all_accuracy)], collapse = "_")
  
    # add accuracy annotation to barplot:
    png(
      paste0(plot_dir, "signal_vs_simulated_CNV_plot.png"), 
      height = 8, 
      width = 22,
      res = 300,
      units = "in"
    )   
      grid.newpage()
  
        pushViewport(viewport(x = 0.06, y = 0.001, width = 0.93, height = 0.9, 
          just = c("left", "bottom")))
          grid.draw(signal_plot)
        popViewport()
        pushViewport(viewport(x = 0.13, y = 0.15, width = 0.793, height = 0.13, 
          just = c("left", "bottom")))
          grid.draw(accuracy_heatmap_obj)
        popViewport()
        for ( e in 1:length(chr_data$lab_pos) ) {
          pushViewport(viewport(x = 0.105 + chr_data$lab_pos[e]/1.27, y = 0.86, width = 0.05, height = 0.05, 
            just = c("left", "bottom")))
            if (e==1) {
              grid.text(names(chr_data$lab_pos)[e], gp=gpar(fontsize=13, fontface = "bold"))
            } else {
              grid.text(gsub("chr", "", names(chr_data$lab_pos)[e]), gp=gpar(fontsize=13, fontface = "bold"))
            }
          popViewport()
        }
      
        pushViewport(viewport(x = unit(1, "cm"), y = unit(2, "cm"), 
                          width = unit(5.5, "cm"), height = unit(4.5, "cm"), 
                          just = c("left", "bottom")))
    
          #grid.rect()
          if (all_accuracy == "false_positive_true_positive_wrong_call") {
            true_pos_false_pos_wrong_legend()
          } else if (all_accuracy == "false_negative_false_positive_true_negative_true_positive_wrong_call") {
            true_pos_neg_false_pos_neg_legend()
          }
        
        popViewport()
    
    dev.off()
  }


  #################################################################################
  #### 13. Save data for distinguishing between different gain/loss copy 
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




