#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

RStudio <- FALSE

project_name <- "thesis"
subproject_name <- "Figure_2.3_simulations"
args = commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
include_t_cells <- as.logical(args[2])
downsample_proportion <- args[3]

sample_name <- "CID4520N_cancer_sim"
include_t_cells <- TRUE
downsample_proportion <- "1"

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("T-cells included?", as.character(include_t_cells)))

library(Seurat)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(dplyr)
library(cowplot)
library(viridis)

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
sim_dir <- paste0(results_dir, "cancer_simulation/", sample_name,
  "/Rdata/")

if (downsample_proportion == "1") {

  if (include_t_cells) {
    in_dir <- paste0(results_dir, "infercnv/t_cells_included/", sample_name, "/")
  } else {
    in_dir <- paste0(results_dir, "infercnv/t_cells_excluded/", sample_name, "/")
  }
  
} else {

  if (include_t_cells) {
    in_dir <- paste0(results_dir, "infercnv/t_cells_included/", sample_name, "/",
      downsample_proportion, "_downsampling/")
  } else {
    in_dir <- paste0(results_dir, "infercnv/t_cells_excluded/", sample_name, "/",
      downsample_proportion, "_downsampling/")
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
  nUMI <- apply(count_df, 2, sum)
  count_df_rownames <- rownames(count_df)
  nGene <- apply(count_df, 2, function(x) length(x[x!=0]))

  QC <- data.frame(
    row.names = colnames(count_df),
    nUMI = nUMI,
    nGene = nGene
  )
  QC <- QC[rownames(epithelial_metadata),]
  epithelial_metadata <- cbind(epithelial_metadata, QC)
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

######


################################################################################
### 3. Calculate sensitivity and specificity  ###
################################################################################

# load known CNVs coordinates from simulated data:  
CNV_indices <- readRDS(paste0(sim_dir, "/CNV_indices.Rdata"))
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
average_epithelial_density_plot <- density(average_epithelial, bw="SJ")
pdf(paste0(plot_dir, "epithelial_CNV_density_plot.pdf"))
  plot(average_epithelial_density_plot, main=NA, xlab = "CNV signal")
dev.off()
png(paste0(plot_dir, "epithelial_CNV_density_plot.png"))
  plot(average_epithelial_density_plot, main=NA, xlab = "CNV signal")
dev.off()

# record gain and losses in InferCNV data:
CNV_region_calls <- lapply(CNV_indices_list, function(x) {
  # segment average_epithelial into CNV_indices segments and record the mean:
  epithelial_segment <- average_epithelial[x$start:x$end]
  x$mean_epithelial_signal <- mean(epithelial_segment)

  # make a call whether gain, neutral or loss using InferCNV signal from this 
  # region based on whether at least 75% is above, equal to or below 1 
  # respectively:
  x$epithelial_CNV_call <- "neutral"
  if (
    length(which(epithelial_segment > 0.95 & epithelial_segment < 1.05)) > (0.75*length(epithelial_segment))
  ) {
    x$epithelial_CNV_call <- "neutral"
  } else if (
    length(which(epithelial_segment > 1.05)) > (0.75*length(epithelial_segment))
  ) {
    x$epithelial_CNV_call <- "gain"
  } else if (
    length(which(epithelial_segment < 0.95)) > (0.75*length(epithelial_segment))
  ) {
    x$epithelial_CNV_call <- "loss"
  }
  return(x)

}

# determine number of true and false positives, true and false negatives:
CNV_accuracy_metrics <- lapply(CNV_region_calls, function(x) {

  x$accuracy <- "true_negative"

  if (CNV_indices$type == "neutral") {

    if (CNV_indices$epithelial_CNV_call == "neutral") {
      x$accuracy <- "true_negative"
    } else if (CNV_indices$epithelial_CNV_call != "neutral") {
      x$accuracy <- "false_positive"
    }

  } else if (CNV_indices$type == "gain") {

    if (CNV_indices$epithelial_CNV_call == "gain") {
      x$accuracy <-  "true_positive"
    } else if (CNV_indices$epithelial_CNV_call == "neutral") {
      x$accuracy <-  "false_negative"
    } else if (CNV_indices$epithelial_CNV_call == "loss") {
      x$accuracy <-  "wrong_call"
    }

  } else if (CNV_indices$type == "loss") {

    if (CNV_indices$epithelial_CNV_call == "loss") {
      x$accuracy <-  "true_positive"
    } else if (CNV_indices$epithelial_CNV_call == "neutral") {
      x$accuracy <-  "false_negative"
    } else if (CNV_indices$epithelial_CNV_call == "gain") {
      x$accuracy <-  "wrong_call"
    }
  }

})

######


################################################################################
### 4. Create QC annotations  ###
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
### 5. Create simulated CNV annotation  ###
################################################################################

# create simulated CNV annotation:
simulated_CNV_plot_data <- readRDS(paste0(sim_dir, 
  "simulated_CNV_plot_data.Rdata"))
log_adjusted_fold_change_df <- 
  simulated_CNV_plot_data$log_adjusted_fold_change_df
CNV_indices <- simulated_CNV_plot_data$CNV_indices
zero_count_value <- simulated_CNV_plot_data$zero_count_value

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
    }
    CNV_indices$end[n] <- CNV_indices$start[n] + (no_genes-1)
  }

}

# record length and midpoint of segments:
CNV_indices$length <- CNV_indices$end - CNV_indices$start
CNV_indices$midpoints <- CNV_indices$start + floor(CNV_indices$length/2)
CNV_indices$number <- seq_along(CNV_indices$start)
CNV_indices$ticks <- "exclude"
CNV_indices$ticks[
  c(9, 11, 13, 15, 17, 19, 21, 24, 26, 28, 30, 32, 34, 36, 38, 41, 43, 
  	45, 47, 49, 51, 53, 55, 57, 59)
] <- "include"

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
for (end in chr_data$ends) {
  p <- p + geom_vline(xintercept=end)
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

pdf(paste0(plot_dir, "sim_CNV_plot.pdf"))
  grid.draw(grid_sim_plot)
dev.off()


################################################################################
### 6. Generate heatmap ###
################################################################################

# fetch chromosome boundary co-ordinates:
if (!file.exists(paste0(Robject_dir, "chromosome_data.Rdata"))) {
  chr_data <- fetch_chromosome_boundaries(epithelial_heatmap, ref_dir)
  saveRDS(chr_data, paste0(Robject_dir, "chromosome_data.Rdata"))
} else {
  chr_data <- readRDS(paste0(Robject_dir, "chromosome_data.Rdata"))
}
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
  show_row_names = F, show_column_names = T,
  column_names_gp = gpar(col = "white"),
  show_row_dend = F,
  heatmap_legend_param = list(title = "CNV\nscore", color_bar = "continuous", 
    grid_height = unit(1.5, "cm"), grid_width = unit(1.5, "cm"), legend_direction = "horizontal",
    title_gp = gpar(fontsize = 18, fontface = "bold"), labels_gp = gpar(fontsize = 12)),
  use_raster = T, raster_device = c("png")
)


ht_list <- final_heatmap + nUMI_annotation + nGene_annotation

annotated_heatmap <- grid.grabExpr(
  draw(ht_list, gap = unit(6, "mm"), heatmap_legend_side = "left")
)

# determine where starting co-ordinates for heatmap are based upon longest cluster name
# (0.00604 units per character):
longest_cluster_name <- max(nchar(unique(as.character(epithelial_metadata$cell_type))))
x_coord <- longest_cluster_name*0.0037


################################################################################
### 6. Plot heatmap ###
################################################################################

# plot final annotated heatmap:
pdf(paste0(plot_dir, "infercnv_plot2.pdf"), height = 13, width = 18)   
  
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
    width = 0.853, height = 0.088, just = c("right", "bottom")))
  grid.draw(grid_sim_plot)
  popViewport()

  pushViewport(viewport(x=x_coord + 0.9155, y=0.1, width = 0.1, height = 0.1, 
    just = "bottom"))
    grid.text("nUMI", rot=65)
  popViewport()
  pushViewport(viewport(x=x_coord + 0.94, y=0.1, width = 0.1, height = 0.1, 
    just = "bottom"))
    grid.text("nGene", rot=65)
  popViewport()
    
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
