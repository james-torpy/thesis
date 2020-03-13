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
}


################################################################################
### 4. Determine presence of artefacts in infercnv output ###
################################################################################

if (!file.exists(paste0(table_dir, "final_artefact_record.txt")) | 
	!file.exists(paste0(table_dir, "artefact_counts.txt"))) {

  # determine average signal across all cells for each gene:
  heatmap_averages <- apply(epithelial_heatmap, 2, mean)
  
  for ( i in 1:(ncol(epithelial_heatmap)-9) ) {
    print(i)
    # determine average signal for every combination of 10 adjacent genes:
    if (i == 1) {
      artefact_record <- data.frame(
        start = i,
        end = i+9,
        average_signal = mean(heatmap_averages[i:(i+9)]),
        type = "no artefact"
      )
    } else {
      artefact_record <- rbind(
        artefact_record,
        data.frame(
          start = i,
          end = i+9,
          average_signal = mean(heatmap_averages[i:(i+9)]),
          type = "no artefact"
        )
      )
    }
    artefact_record$type <- as.character(artefact_record$type)
  
    # report segment type:
    if (artefact_record$average_signal[i] < neutral_signal_range[1]) {
      artefact_record$type[i] <- "loss"
    } else if (artefact_record$average_signal[i] > neutral_signal_range[2]) {
      artefact_record$type[i] <- "gain"
    }
  
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

  # fix overlapping indices by prioritising artefact annotation over 
  # no artefact:
  no_artefact_indices <- which(final_artefact_record$type == "no artefact")
  for (ind in no_artefact_indices) {
  	# if the first non-artefact region is not at the start of the genome,
  	#  make start of non-artefact region one after the end of the previous 
  	# artefact region
  	if (final_artefact_record$start[ind] != 1) {
  	  final_artefact_record$start[ind] <- final_artefact_record$end[ind-1] + 1
  	}
  	# if the last non-artefact region is not at the end of the genome,
  	# make end of non-artefact region one before the start of the next 
  	# artefact region: 
  	if (final_artefact_record$end[ind] != ncol(epithelial_heatmap)) {
      final_artefact_record$end[ind] <- final_artefact_record$start[ind+1] - 1
    }
  }

  artefact_indices <- final_artefact_record[
    final_artefact_record$type != "no artefact",
  ]
  for (m in 1:nrow(artefact_indices)) {
    if (m==1) {
      heatmap_artefact <- artefact_indices$start[m]:artefact_indices$end[m]
    } else {
      heatmap_artefact <- c(heatmap_artefact,
        artefact_indices$start[m]:artefact_indices$end[m]
      )
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
    levels = c("no artefact", "gain", "loss")
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

} else {

  final_artefact_record <- read.table(
  	paste0(table_dir, "final_artefact_record.txt"),
  	sep = "\t",
  	header = T,
  	as.is = T
  )

  artefact_count <- read.table(
  	paste0(table_dir, "artefact_counts.txt"),
  	sep = "\t",
  	header = T,
  	as.is = T
  )

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
  bottom_annotation = artefact_annotation,
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
pdf(
  paste0(plot_dir, "infercnv_plot_with_artefact_calls.pdf"), 
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

  pushViewport(viewport(x = x_coord + 0.964, y = 0.035, 
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

##convert pdf to png:
#system(paste0("for p in ", plot_dir, "*.pdf; do echo $p; f=$(basename $p); echo $f; ",
#              "new=$(echo $f | sed 's/.pdf/.png/'); echo $new; ", 
#              "convert -density 150 ", plot_dir, "$f -quality 90 ", plot_dir, "$new; done"))
#
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          