#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

RStudio <- FALSE

project_name <- "thesis"
subproject_name <- "Figure_2.6_random_permutation_of_normal"
args = commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
permutation_proportion <- args[2]
simulation_number <- args[3]
t_cells_included <- as.logical(args[4])
analysis_mode <- args[5]
neutral_signal_range <- unlist(strsplit(args[6], split = "_"))

#sample_name <- "permutated_CID4520N"
#permutation_proportion <- 0.01
#simulation_number <- "1"
#t_cells_included <- TRUE
#analysis_mode <- "samples"
#neutral_signal_range <- unlist(strsplit("0.97_1.03", split = "_"))

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("Proportion of genes permutated = ", permutation_proportion))
print(paste0("Simulation number = ", simulation_number))
print(paste0("T-cells included? ", as.character(t_cells_included)))
print(paste0("Analysis mode = ", analysis_mode))
print("Neutral signal range = ")
print(neutral_signal_range)

library(Seurat)
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(cowplot)
library(viridis)
library(GenomicRanges)
library(ggplot2)

if (RStudio) {
  
  library(ComplexHeatmap)
  library(circlize)
  library(scales)
  library(fpc)
  library(naturalsort)
  
  home_dir <- "/Users/jamestorpy/clusterHome/"
  
} else {
  
  lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
  library(Polychrome, lib.loc=lib_loc)
  library(ComplexHeatmap, lib.loc=lib_loc)
  library(circlize, lib.loc = lib_loc)
  library(fpc, lib.loc = lib_loc)
  library(naturalsort, lib.loc = lib_loc)
  
  home_dir <- "/share/ScratchGeneral/jamtor/"
  
}

project_dir <- paste0(home_dir, "projects/", 
  project_name, "/", subproject_name, "/")
ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/functions/")
results_dir <- seurat_path <- paste0(project_dir, "results/")

#orig_sample_name <- gsub("permutated_", "", sample_name)
#seurat_dir <- paste0(project_dir, "raw_files/seurat_objects/", 
#  orig_sample_name, "/")

perm_dir <- paste0(results_dir, sample_name, "/", permutation_proportion, 
  "_proportion/", simulation_number, "/Rdata/")

#general_sim_dir <- paste0(results_dir, "cancer_simulation/", sample_name,
#  "/Rdata/")

if (t_cells_included) {
  in_path <- paste0(results_dir, "infercnv/t_cells_included/", 
    sample_name, "/", permutation_proportion, "_proportion/", 
    simulation_number, "/")
  
} else {
  in_path <- paste0(results_dir, "infercnv/t_cells_excluded/", 
    sample_name, "/", permutation_proportion, "_proportion/", 
    simulation_number, "/")
}
in_dir <- paste0(in_path, analysis_mode, "_mode/")
input_dir <- paste0(in_path, "input_files/")

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

#create_extended_vector <- function(df, column) {
#  for (j in 1:nrow(df)) {
#    vector_length <- length(seq(df$start[j], df$end[j]))
#    if (j==1) {
#      result_vector <- rep(
#        eval(parse(text=paste0("df$", column, "[j]"))), vector_length
#      )
#    } else {
#      result_vector <- c(
#        result_vector,
#        rep(
#          eval(parse(text=paste0("df$", column, "[j]"))), vector_length
#        )
#      )
#    }
#  }
#  return(result_vector)
#}


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
### 3. Format permutation indices  ###
################################################################################

# load permutation information:
permutation_data <- readRDS(paste0(perm_dir, "3.final_permutation_data.Rdata"))

# fetch chromosome boundary co-ordinates:
if (!file.exists(paste0(Robject_dir, "chromosome_data.Rdata"))) {
  chr_data <- fetch_chromosome_boundaries(epithelial_heatmap, ref_dir)
  saveRDS(chr_data, paste0(Robject_dir, "chromosome_data.Rdata"))
} else {
  chr_data <- readRDS(paste0(Robject_dir, "chromosome_data.Rdata"))
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
### 6. Determine accurately called CNVs/gaps ###
################################################################################

if (gap_or_CNV == "CNV") {
  if (
    !file.exists(paste0(table_dir, "CNV_calls.txt")) |
    !file.exists(paste0(Robject_dir, "3.final_CNV_indices.Rdata"))
  ) {

    feature_indices <- CNV_indices
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

    # annotate all loss ranges with mean signal >= lower range of 
    # neutral_ranges as "feature_not_called":
    feature_indices$call[
      feature_indices$multiplier == 0 & 
      feature_indices$mean_signal >= neutral_signal_range[1]
    ] <- "feature_not_called"
    # annotate all loss ranges with mean signal < lower range of 
    # neutral_ranges as "feature_called":
    feature_indices$call[
      feature_indices$multiplier == 0 & 
      feature_indices$mean_signal < neutral_signal_range[1]
    ] <- "feature_called"
   
    # annotate all gain ranges with mean signal <= upper range of 
    # neutral_ranges as "feature_not_called":
    feature_indices$call[
      feature_indices$multiplier == 3 & 
      feature_indices$mean_signal <= neutral_signal_range[2]
    ] <- "feature_not_called"
    # annotate all gain ranges with mean signal > upper range of 
    # neutral_ranges as "feature_called":
    feature_indices$call[
      feature_indices$multiplier == 3 & 
      feature_indices$mean_signal > neutral_signal_range[2]
    ] <- "feature_called"

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

    # for each range, determine whether there is a 40 gene window with mean
    # signal within neutral_signal_ranges or the range of the opposite CNV type:
    feature_indices$call <- "feature_not_called"
    for (r in 1:nrow(feature_indices)) {

      for (g in feature_indices[r,]$start:feature_indices[r,]$end) {

        segment_mean <- mean(average_epithelial[g:g+39])

        if (CNV_type == "gain") {
          if (segment_mean <= neutral_signal_range[2]) {
            feature_indices$call[r] <- "feature_called"
            break()
          }
        } else if (CNV_type == "loss") {
          if (segment_mean >= neutral_signal_range[1]) {
            feature_indices$call[r] <- "feature_called"
            break()
          }
        }

      }

    }

#    # for each range, fetch mean infercnv signal:
#    feature_indices$mean_signal <- NA
#    for (r in 1:nrow(feature_indices)) {
#      signal_segment <- average_epithelial[
#        feature_indices$start[r]:feature_indices$end[r]
#      ]
#      feature_indices$mean_signal[r] <- mean(signal_segment)
#    }
#
#    # annotate non-called gap regions:
#    feature_indices$call <- "feature_not_called"
#
#    # annotate all gap ranges which have a gap of at least 40 genes within 
#    # neutral_signal_ranges or the range of the opposite CNV type as "feature_called":
#    if (CNV_type == "gain") {
#      feature_indices$call[
#        feature_indices$mean_signal <= neutral_signal_range[2]
#      ] <- "feature_called"
#    } else if (CNV_type == "loss") {
#      feature_indices$call[
#        feature_indices$mean_signal >= neutral_signal_range[1]
#      ] <- "feature_called"
#    }

    # order feature_indices and fill in areas of no gaps:
    feature_indices <- feature_indices[order(feature_indices$start),]
    feature_gr <- GRanges(
      seqnames = Rle("genome"),
      ranges = IRanges(start = feature_indices$start, end = feature_indices$end),
      strand = Rle("*"),
      call = feature_indices$call
    )
    feature_gaps <- gaps(feature_gr)
    feature_gaps$call <- "non-gap"
    feature_gr_full <- c(feature_gr, feature_gaps)
    feature_gr_full <- feature_gr_full[order(end(ranges(feature_gr_full)))]
    
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
### 7. Annotate calls ###
################################################################################

if (gap_or_CNV == "CNV") {

  # change annotation labels:
  call_annot_indices <- feature_indices
  call_annot_indices$call <- gsub("feature", "CNV", call_annot_indices$call)
  call_annot_indices$call <- gsub("_", " ", call_annot_indices$call)
  # create annotation of all calls:
  call_annotation_vector <- factor(create_extended_vector(feature_indices, "call"))
  levels(call_annotation_vector) <- c("CNV called", "VNC not called", "non-gap")

  cols <- c("#F6DC15", "#7CBA61", "#E7E4D3")
  names(cols) <- levels(call_annotation_vector)
  
  call_annotation <- HeatmapAnnotation(
    accuracy = call_annotation_vector,
    col = list(accuracy = cols),
    annotation_name_side = "left",
    annotation_legend_param = list(title = "", 
    	labels_gp = gpar(fontsize = 12))
  )
  call_annotation@name <- "call"

} else if (gap_or_CNV == "gap") {

  # change annotation labels:
  call_annot_indices <- feature_indices
  call_annot_indices$call <- gsub("feature", "gap", call_annot_indices$call)
  call_annot_indices$call <- gsub("_", " ", call_annot_indices$call)
  # create annotation of all calls:
  call_annotation_vector <- factor(create_extended_vector(call_annot_indices, "call"))
  levels(call_annotation_vector) <- c("gap called", "gap not called", "non-gap")

  cols <- c("#F6DC15", "#7CBA61", "#E7E4D3")
  names(cols) <- levels(call_annotation_vector)
  
  call_annotation <- HeatmapAnnotation(
    accuracy = call_annotation_vector,
    col = list(accuracy = cols),
    annotation_name_side = "left",
    annotation_legend_param = list(title = "", 
      labels_gp = gpar(fontsize = 12))
  )
  call_annotation@name <- "call"

}


################################################################################
### 8. Generate heatmap ###
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

final_heatmap <- Heatmap(
  plot_object, name = paste0("hm"), 
  col = heatmap_cols,
  cluster_columns = F, cluster_rows = F,
  show_row_names = F, show_column_names = F,
  #column_names_gp = gpar(col = "white"),
  show_row_dend = F,
  bottom_annotation = call_annotation,
  heatmap_legend_param = list(title = "CNV\nscore", 
  at = c(round(min(na_less_vector), 1), 1, round(max(na_less_vector), 1)),
  color_bar = "continuous", grid_height = unit(1.5, "cm"), 
  grid_width = unit(1.5, "cm"), legend_direction = "horizontal",
  title_gp = gpar(fontsize = 18, fontface = "bold"), 
  labels_gp = gpar(fontsize = 12)),
  use_raster = T, raster_device = c("png")
)

annotated_heatmap <- grid.grabExpr(
  draw(final_heatmap, gap = unit(6, "mm"), heatmap_legend_side = "left",
  annotation_legend_side = "left")
)
dev.off()

# determine where starting co-ordinates for heatmap are based upon longest cluster name
# (0.00604 units per character):
longest_cluster_name <- max(nchar(unique(as.character(epithelial_metadata$cell_type))))
x_coord <- longest_cluster_name*0.0037


################################################################################
### 11. Plot heatmap ###
################################################################################

# plot final annotated heatmap:
if (gap_or_CNV == "CNV") {
  pdf(
    paste0(plot_dir, "infercnv_plot_with_CNV_calls.pdf"), 
    height = 13, width = 20
  )   
} else {
  pdf(
    paste0(plot_dir, "infercnv_plot_with_gap_calls.pdf"), 
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

##convert pdf to png:
#system(paste0("for p in ", plot_dir, "*.pdf; do echo $p; f=$(basename $p); echo $f; ",
#              "new=$(echo $f | sed 's/.pdf/.png/'); echo $new; ", 
#              "convert -density 150 ", plot_dir, "$f -quality 90 ", plot_dir, "$new; done"))
#
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          