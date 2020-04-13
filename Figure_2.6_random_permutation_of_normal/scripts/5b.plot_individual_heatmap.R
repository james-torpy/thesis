#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

RStudio <- FALSE

project_name <- "thesis"
subproject_name <- "Figure_2.6_random_permutation_of_normal"
args = commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
permutation_proportion <- args[2]
simulation_number <- args[3]
analysis_mode <- args[4]
neutral_signal_range <- unlist(strsplit(args[5], split = "_"))
min_artefact_proportion <- as.numeric(args[6])
min_artefact_length <- as.numeric(args[7])

#sample_name <- "permutated_CID4520N"
#permutation_proportion <- 0.3
#simulation_number <- "1"
#analysis_mode <- "samples"
#neutral_signal_range <- unlist(strsplit("0.97_1.03", split = "_"))
#min_artefact_proportion <- 0.4
#min_artefact_length <- 20

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("Proportion of genes permutated = ", permutation_proportion))
print(paste0("Simulation number = ", simulation_number))
print(paste0("Analysis mode = ", analysis_mode))
print("Neutral signal range = ")
print(neutral_signal_range)
print(paste0("Min proportion of cells for artefact detection = ", min_artefact_proportion))
print(paste0("Min length for artefact detection = ", min_artefact_length))

library(Seurat)
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(cowplot)
library(viridis)
library(GenomicRanges)
library(ggplot2)
library(data.table)

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

perm_dir <- paste0(results_dir, sample_name, "/", permutation_proportion, 
  "_proportion/", simulation_number, "/Rdata/")

in_path <- paste0(results_dir, "infercnv/", sample_name, "/", 
  permutation_proportion, "_proportion/", simulation_number, "/")

in_dir <- paste0(in_path, analysis_mode, "_mode/")
input_dir <- paste0(in_path, "input_files/")

out_dir <- paste0(in_dir, min_artefact_proportion, 
  "_min_artefact_proportion/", min_artefact_length, 
  "_min_artefact_length/")

Robject_dir <- paste0(out_dir, "Rdata/")
system(paste0("mkdir -p ", Robject_dir))
plot_dir <- paste0(out_dir, "plots/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(out_dir, "tables/")
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
      result_vector <- as.character(
      	rep(
          eval(parse(text=paste0("df$", column, "[j]"))), vector_length
        )
      )
    } else {
      result_vector <- c(
        result_vector,
        as.character(
          rep(
            eval(parse(text=paste0("df$", column, "[j]"))), vector_length
          )
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
### 3. Load permutation indices and create annotation ###
################################################################################

# fetch chromosome boundary co-ordinates:
chr_data <- fetch_chromosome_boundaries(epithelial_heatmap, ref_dir)

# load permutation information:
permutation_data <- readRDS(paste0(perm_dir, "3.final_permutation_data.Rdata"))

# subset annotation_data and reset indices column according to what genes are 
# present in epithelial_heatmap:
annotation_data <- permutation_data$permutation_record
annotation_data <- annotation_data[
  which(
    annotation_data$genes %in% colnames(epithelial_heatmap)
  ), 
]
for (i in 1:nrow(annotation_data)) {
  annotation_data$indices[i] <- grep(
  	paste0("\\b", annotation_data$genes[i], "\\b"), 
  	colnames(epithelial_heatmap)
  )
}

# subset to indices and values only:
annotation_data <- subset(annotation_data, select = -genes)

# plot permutated gene annotation:
p <- ggplot(annotation_data, aes(x=indices, y=values))
p <- p + geom_point(size = 0.8, colour = "#46B711")
p <- p + scale_x_continuous(
  limits = c(1, ncol(epithelial_heatmap)), 
  expand = c(0, 0)
)
p <- p + scale_y_continuous(
  breaks = c(0, 0.5, 1, 1.5, 2, 3),
  limits = c(0, 3)
)
for (c_end in chr_data$ends) {
  p <- p + geom_vline(xintercept=c_end)
}

# create line at 1:
p <- p + geom_segment(
  x=1,
  xend=annotation_data$indices[
    nrow(annotation_data)
  ],
  y=1,
  yend=1
)
# remove axis labels:
p <- p + theme(
  axis.title.x=element_blank(),
  axis.text.x=element_blank(),
  axis.ticks.x=element_blank(),
  axis.title.y=element_blank()
)
grid_perm_plot <- ggplotGrob(p)
dev.off()

if (!file.exists(paste0(plot_dir, "permutated_gene_plot.pdf"))) {
  pdf(
  	paste0(plot_dir, "permutated_gene_plot.pdf"), 
  	width = 15, height = 10
  )
    grid.draw(grid_perm_plot)
  dev.off()
  png(
    paste0(plot_dir, "permutated_gene_plot.png"), 
    width = 15, height = 10,
    units = "in",
    res = 300
  )
    grid.draw(grid_perm_plot)
  dev.off()
}


################################################################################
### 4. Determine presence of artefacts in infercnv output ###
################################################################################

if (!file.exists(paste0(table_dir, "final_artefact_record.txt")) | 
	!file.exists(paste0(table_dir, "artefact_counts.txt")) |
	!file.exists(paste0(table_dir, "artefact_lengths.txt")) |
	!file.exists(paste0(Robject_dir, "artefact_annotation.Rdata"))) {

  # determine average signal across all cells for each gene:
  heatmap_averages <- apply(epithelial_heatmap, 2, mean)
  
  for ( i in 1:(ncol(epithelial_heatmap)-(min_artefact_length-1)) ) {

   averages_per_cell <- apply(epithelial_heatmap[i:(i+(min_artefact_length-1))], 1, mean)
    
    # if at least 10% of cells have average signal in loss range, label as loss
    # and vice versa for gains:
    if ( length(which(averages_per_cell < neutral_signal_range[1])) > 
      floor(min_artefact_proportion*length(averages_per_cell)) ) {

      segment_record <- data.frame(
        start = i,
        end = i+(min_artefact_length-1),
        average_signal = mean(heatmap_averages[i:(i+(min_artefact_length-1))]),
        type = "loss"
      )

    } else if ( length(which(averages_per_cell > neutral_signal_range[2])) > 
      floor(min_artefact_proportion*length(averages_per_cell)) ) {

      segment_record <- data.frame(
        start = i,
        end = i+(min_artefact_length-1),
        average_signal = mean(heatmap_averages[i:(i+(min_artefact_length-1))]),
        type = "gain"
      )

    } else {

      segment_record <- data.frame(
        start = i,
        end = i+(min_artefact_length-1),
        average_signal = mean(heatmap_averages[i:(i+(min_artefact_length-1))]),
        type = "no_artefact"
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
  final_artefact_record$type <- as.character(final_artefact_record$type)

  # fix overlapping indices by prioritising artefact annotation over 
  # no artefact:
  no_artefact_indices <- which(final_artefact_record$type == "no_artefact")
  for (ind in no_artefact_indices) {
  	# if the first non-artefact region is not at the start of the genome,
  	#  make start of non-artefact region one after the end of the previous 
  	# artefact region
  	if (final_artefact_record$start[ind] != 1) {
  	  final_artefact_record$start[ind] <- final_artefact_record$end[(ind-1)] + 1
  	}
  	# if the last non-artefact region is not at the end of the genome,
  	# make end of non-artefact region one before the start of the next 
  	# artefact region: 
  	if (final_artefact_record$end[ind] != ncol(epithelial_heatmap)) {
  	  if (final_artefact_record$start[ind] < final_artefact_record$start[(ind+1)]) {
        final_artefact_record$end[ind] <- final_artefact_record$start[(ind+1)] - 1  	  	
  	  } else {
  	  	final_artefact_record$type[ind] <- "remove"
  	  	final_artefact_record$start[(ind+1)] <- final_artefact_record$end[(ind-1)] + 1
  	  }
    }
  }
  final_artefact_record <- final_artefact_record[
    final_artefact_record$type != "remove",
  ]

  artefact_indices <- final_artefact_record[
    final_artefact_record$type != "no_artefact",
  ]

  if (nrow(artefact_indices) > 0) {
    for (m in 1:nrow(artefact_indices)) {
      if (m==1) {
        heatmap_artefact <- artefact_indices$start[m]:artefact_indices$end[m]
      } else {
        heatmap_artefact <- c(heatmap_artefact,
          artefact_indices$start[m]:artefact_indices$end[m]
        )
      }   
    }

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
  
    # count gains and losses:
    artefact_count <- data.frame(
    	row.names = c("gain", "loss"),
    	count = c(
    	  length(which(final_artefact_record$type == "gain")),
    	  length(which(final_artefact_record$type == "loss"))
    	)
    )
    
    # create artefact annotation:
    artefact_vector <- factor(
      create_extended_vector(final_artefact_record, "type"),
      levels = c("no_artefact", "gain", "loss")
    )
    cols <- c("#E5E4DF", "#BF3667", "#58B9DB")
    names(cols) <- levels(artefact_vector)
    
    artefact_annotation <- HeatmapAnnotation(
      artefact_call = artefact_vector,
      col = list(artefact_call = cols),
      annotation_name_side = "left",
      annotation_legend_param = list(title = "", 
      	labels_gp = gpar(fontsize = 12))
    )
    artefact_annotation@name <- "artefact"

    # annotate artefact lengths:
    artefact_lengths <- data.frame(
      length = artefact_indices$end - artefact_indices$start
    )
    artefact_lengths$midpoint <- floor(
    	artefact_indices$start + artefact_lengths$length/2
    )
    artefact_lengths$position <- 
      artefact_lengths$midpoint/ncol(epithelial_heatmap)


    ################################################################################
    ### 5. Save artefact metadata ###
    ################################################################################
  
    # write artefact results as tables:
    write.table(
   	  final_artefact_record, 
   	  paste0(table_dir, "final_artefact_record.txt"),
   	  sep = "\t",
   	  quote = F,
   	  row.names = F,
   	  col.names = T
    )
  
    write.table(
   	  artefact_count, 
   	  paste0(table_dir, "artefact_counts.txt"),
   	  sep = "\t",
   	  quote = F,
   	  row.names = T,
   	  col.names = T
    )

    write.table(
   	  artefact_lengths, 
   	  paste0(table_dir, "artefact_lengths.txt"),
   	  sep = "\t",
   	  quote = F,
   	  row.names = T,
   	  col.names = T
    )
  
    saveRDS(
    	artefact_annotation, 
    	paste0(Robject_dir, "artefact_annotation.Rdata")
    )

  } else if (nrow(artefact_indices) == 0) {

    artefact_count <- data.frame(
      row.names = c("gain", "loss"),
      count = c(0, 0)
    )

    write.table(
      artefact_count, 
      paste0(table_dir, "artefact_counts.txt"),
      sep = "\t",
      quote = F,
      row.names = T,
      col.names = T
    )

    artefact_lengths <- data.frame(
      length = 0,
      midpoint = 0
    )

    write.table(
   	  artefact_lengths, 
   	  paste0(table_dir, "artefact_lengths.txt"),
   	  sep = "\t",
   	  quote = F,
   	  row.names = T,
   	  col.names = T
    )

  }

} else {

  if (file.exists(paste0(table_dir, "final_artefact_record.txt"))) {
    final_artefact_record <- read.table(
      paste0(table_dir, "final_artefact_record.txt"),
      sep = "\t",
      header = T,
      as.is = T
    )
  }

  artefact_count <- read.table(
 	  paste0(table_dir, "artefact_counts.txt"),
 	  sep = "\t",
 	  header = T,
 	  as.is = T
  )
  artefact_lengths <- read.table(
  	paste0(table_dir, "artefact_lengths.txt"),
  	  sep = "\t",
 	  header = T,
 	  as.is = T
  )

  if (file.exists(paste0(table_dir, "final_artefact_record.txt"))) {
    artefact_annotation <- readRDS(
    	paste0(Robject_dir, "artefact_annotation.Rdata")
    )
  }
  
}


################################################################################
### 6. Generate heatmap ###
################################################################################

if (!file.exists(paste0(plot_dir, "infercnv_plot_with_artefact_calls.pdf"))) {
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
  
  if (sum(artefact_count) > 0) {
    final_heatmap <- Heatmap(
      plot_object, name = paste0("hm"), 
      col = heatmap_cols,
      cluster_columns = F, cluster_rows = F,
      show_row_names = F, show_column_names = F,
      #column_names_gp = gpar(col = "white"),
      show_row_dend = F,
      bottom_annotation = artefact_annotation,
      heatmap_legend_param = list(title = "CNV\nscore", 
      at = c(round(min(na_less_vector), 1), 1, round(max(na_less_vector), 1)),
      color_bar = "continuous", grid_height = unit(1.5, "cm"), 
      grid_width = unit(1.5, "cm"), legend_direction = "horizontal",
      title_gp = gpar(fontsize = 18, fontface = "bold"), 
      labels_gp = gpar(fontsize = 12)),
      use_raster = T, raster_device = c("png")
    )
  } else {
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
  }
  
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
  ### 7. Plot heatmap ###
  ################################################################################
  
  # plot final annotated heatmap:
  pdf(
    paste0(plot_dir, "infercnv_plot_with_artefact_calls.pdf"), 
    height = 13, width = 20
  )   
    
    grid.newpage()
    if (sum(artefact_count) > 0) {
      pushViewport(viewport(x = 0.01, y = 0.16, width = 0.99, height = 0.78, 
        just = c("left", "bottom")))
    } else {
      pushViewport(viewport(x = 0.065, y = 0.16, width = 0.934, height = 0.78, 
        just = c("left", "bottom")))
    }
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
  
    for (i in 1:nrow(artefact_lengths)) {
      pushViewport(viewport(x = x_coord + (1/artefact_lengths$position[i]*0.027) + 
        artefact_lengths$position[i], y = 0.085, 
        width = 0.1, height = 0.1, just = c("right", "bottom")))
        grid.text(artefact_lengths$length[i], gp=gpar(fontsize=14))
      popViewport()
    }
    
    pushViewport(viewport(x=x_coord + 0.07, y=0.135, width = 0.05, height = 0.1, 
      just = "right"))
      grid.text("artefact lengths", gp=gpar(fontsize=14))
    popViewport()
  
  
    pushViewport(viewport(x = x_coord + 0.964, y = 0.001, 
      width = 0.909, height = 0.12, just = c("right", "bottom")))
      grid.draw(grid_perm_plot)
    popViewport()
  
    pushViewport(viewport(x=x_coord + 0.0765, y=0.89, width = 0.1, height = 0.1, 
      just = "right"))
      grid.text(paste0("No. gain calls = ", artefact_count["gain",]), 
        gp=gpar(fontsize=16))
    popViewport()
    pushViewport(viewport(x=x_coord + 0.071, y=0.86, width = 0.1, height = 0.1, 
    just = "right"))
      grid.text(paste0("No. loss calls = ", artefact_count["loss",]), 
        gp=gpar(fontsize=16))
    popViewport()
    
  dev.off()

}


################################################################################
### 7. Plot average signal vs genes permutated ###
################################################################################

if (!file.exists(paste0(plot_dir, "annotated_artefact_level_plot.png"))) {

  # calculate mean signal per gene:
  mean_signal <- apply(epithelial_heatmap, 2, mean)
  
  # split mean signal by 10 gene intervals:
  number_divisions <- ceiling(length(mean_signal)/10)
  for (i in 1:number_divisions) {
    if (i==1) {
      split_vector <- rep(i, 10)
    } else {
      split_vector <- c(split_vector, rep(i, 10))
    }
  }
  
  # calculate how many elements less in the last split group and remove:
  missing_elements <- (number_divisions*10)-length(mean_signal)
  split_vector <- split_vector[1:(length(split_vector)-missing_elements)]
  
  # split mean signal:
  signal_per_10 <- unlist(
    lapply(
      split(mean_signal, split_vector), mean
    )
  )
  # add start and end co-ordinates and scale signal values:
  signal_df <- data.frame(
    start = ((as.numeric(names(signal_per_10))-1)*10)+1,
    end = as.numeric(names(signal_per_10))*10,
    mean_signal = round(signal_per_10, 3)
  )
  # adjust last end co-ordinate for missing values:
  signal_df$end[nrow(signal_df)] <- length(mean_signal)
  
  # collapse entries with the same mean signal:
  split_df <- split(signal_df, rleid(signal_df$mean_signal))
  collapsed_signal_df <- do.call(
    "rbind",
    lapply(split_df, function(x) {
      return(
        data.frame(
          start = x$start[1],
          end = x$end[nrow(x)],
          mean_signal = x$mean_signal[1]
        )
      )
    }) 
  )
  
  # label loss and gain signal:
  collapsed_signal_df$type[
    collapsed_signal_df$mean_signal == 1
  ] <- "neutral"
  collapsed_signal_df$type[
    collapsed_signal_df$mean_signal > 1
  ] <- "gain"
  collapsed_signal_df$type[
    collapsed_signal_df$mean_signal < 1
  ] <- "loss"
  
  # determine midpoints of chromosomes:
  chr_data$midpoints <- chr_data$ends-floor(chr_data$lengths/2)
  
  # define ylim as 0.01 less and more than max and min signal values rounded to 
  # 2 decimal places:
  #y_lim <- round(range(collapsed_signal_df$mean_signal), 2)
  #y_lim[1] <- y_lim[1]-0.02
  #y_lim[2] <- y_lim[2]+0.01
  y_lim <- c(0.95, 1.05)
  
  # change values of annotation_data to match those of collapsed_signal_df:
  perm_data <- annotation_data
  perm_data$plot_values <- perm_data$values
  perm_data$plot_values[perm_data$plot_values == 3] <- y_lim[2]
  perm_data$plot_values[perm_data$plot_values == 2] <- 1 + 2*((y_lim[2]-1)/4)
  perm_data$plot_values[perm_data$plot_values == 1.5] <- 1 + ((y_lim[2]-1)/4)
  perm_data$plot_values[perm_data$plot_values == 0.5] <- 1 - 2*(abs(y_lim[1]-1)/4)
  perm_data$plot_values[perm_data$plot_values == 0] <- y_lim[1]
  
  # plot CNV signal across genome:
  p <- ggplot(perm_data, aes(x=indices, y=plot_values))
  ## highlight neutral region:
  #p <- p + geom_rect(
  #  data=NULL,
  #  aes(xmin=-Inf, xmax=Inf, ymin=0.97, ymax=1.03),
  #  fill="lightgreen"
  #)
  p <- p + geom_point(size = 2, colour = "#46B711")
  p <- p + scale_y_continuous(
    breaks = c(y_lim[1], 1, y_lim[2]),
    limits = c(y_lim[1], y_lim[2]),
    labels = c(y_lim[1], 1, y_lim[2])
  )
  p <- p + scale_x_continuous(
    limits = c(
      0, collapsed_signal_df$end[nrow(collapsed_signal_df)]
    ), 
    expand = c(0, 0),
    breaks = chr_data$midpoints,
    labels = seq(1, 22)
  )
  
  for (c_end in chr_data$ends) {
    p <- p + geom_vline(xintercept=c_end, color="#A0A0A0", size = 0.2)
  }
  
  for (r in 1:nrow(collapsed_signal_df)) {
  
    # make gains red, losses blue and neutral regions black:
    if (collapsed_signal_df$type[r] == "neutral") {
      temp_col <- "#777777"
    } else if (collapsed_signal_df$type[r] == "gain") {
      temp_col <- "#BF3667"
    } else if (collapsed_signal_df$type[r] == "loss") {
      temp_col <- "#58B9DB"
    }
  
    # create horizontal line:
    p <- p + geom_segment(
      x=collapsed_signal_df$start[r], 
      xend=collapsed_signal_df$end[r], 
      y=collapsed_signal_df$mean_signal[r], 
      yend=collapsed_signal_df$mean_signal[r], 
      size=0.6, color=temp_col
    )
  
    if (collapsed_signal_df$type[r] != "neutral") {
      # create left vertical line:
      if (r != 1) {
        p <- p + geom_segment(
          x=collapsed_signal_df$start[r], 
          xend=collapsed_signal_df$start[r], 
          y=1, 
          yend=collapsed_signal_df$mean_signal[r], 
          size=0.6, color=temp_col
        )
      }
      # create left vertical line:
      if (r != nrow(collapsed_signal_df)) {
        p <- p + geom_segment(
          x=collapsed_signal_df$end[r], 
          xend=collapsed_signal_df$end[r], 
          y=1, 
          yend=collapsed_signal_df$mean_signal[r], 
          size=0.6, color=temp_col
        )
      }
    }
  }
  
  # remove axis labels:
  p <- p + theme(
    axis.title.x=element_blank(),
  #  axis.text.x=element_blank(),
    axis.ticks.x=element_blank()
  )
  
  # label y-axis:
  p <- p + ylab("Mean CNV signal")

  print(paste0("Creating annotated artefact level plot..."))

  artefact_level_plot <- ggplotGrob(p)
  dev.off()
  
#  pdf(paste0(plot_dir, "artefact_level_plot.pdf"), width = 20)
#    grid.draw(artefact_level_plot)
#  dev.off()
#  png(paste0(plot_dir, "artefact_level_plot.png"), 
#    width = 20,
#    height = 10,
#    units = "in",
#    res = 300
#  )
#    grid.draw(artefact_level_plot)
#  dev.off()
  
  # plot artefact plot with detected artefact annotation:
  annotation_grid <- grid.grabExpr(
    draw(artefact_annotation)
  )
  dev.off()
  
#  pdf(paste0(plot_dir, "annotated_artefact_level_plot.pdf"), width = 20)
#  
#    grid.newpage()
#      pushViewport(viewport(x = 0.065, y = 0.16, width = 0.934, height = 0.78, 
#        just = c("left", "bottom")))
#        grid.draw(artefact_level_plot)
#      popViewport()
#      pushViewport(viewport(x = 0.106, y = 0.1, 
#        width = 0.89, height = 0.12, just = c("left", "bottom")))
#        grid.draw(annotation_grid)
#      popViewport()
#  
#  dev.off()

  print(paste0("Creating annotated artefact level plot file..."))
  png(
    paste0(plot_dir, "annotated_artefact_level_plot.png"), 
    width = 20,
    height = 10,
    units = "in",
    res = 300
  )
  
    grid.newpage()
      pushViewport(viewport(x = 0.065, y = 0.16, width = 0.934, height = 0.78, 
        just = c("left", "bottom")))
        grid.draw(artefact_level_plot)
      popViewport()
      pushViewport(viewport(x = 0.106, y = 0.1, 
        width = 0.89, height = 0.12, just = c("left", "bottom")))
        grid.draw(annotation_grid)
      popViewport()
  
  dev.off()

}


