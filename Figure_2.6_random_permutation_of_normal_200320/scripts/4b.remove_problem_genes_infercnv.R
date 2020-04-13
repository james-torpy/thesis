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
# proportion of genes to be randomly permuated:
permutation_proportion <- as.numeric(args[2])
print(paste0("permutation proportion = ", permutation_proportion))
simulation_number <- as.numeric(args[3])
print(paste0("Simulation number = ", simulation_number))
t_cells_included <- as.logical(args[4])
print(paste0("T-cells included in normal InferCNV run? ", 
  t_cells_included))
analysis_mode <- args[5]
print(paste0("Analysis mode of normal InferCNV run = ", 
  analysis_mode, "_mode"))
neutral_signal_range <- unlist(strsplit(args[6], split = "_"))
print("Neutral signal range = ")
print(neutral_signal_range)
numcores <- as.numeric(args[7])
print(paste0("Number cores = ", numcores))
remove_genes_number <- as.numeric(args[8])
print(paste0("Number of times repeated after this run = ", remove_genes_number))


#project_name <- "thesis"
#subproject_name <- "Figure_2.6_random_permutation_of_normal"
#sample_name <- "permutated_CID4520N"
#print(paste0("sample name = ", sample_name))
# proportion of genes to be randomly permuated:
#permutation_proportion <- 0.005
#print(paste0("permutation proportion = ", permutation_proportion))
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
#numcores <- 6
#print(paste0("Number cores = ", numcores))
#remove_genes_number <- 3
#print(paste0("Number of times repeated after this run = ", repeat_number))


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
library(infercnv, lib.loc=lib_loc)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", 
  project_name, "/", subproject_name, "/")
ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/functions/")
results_dir <- paste0(project_dir, "results/")

normal_artefact_dir <- paste0(results_dir, sample_name, "/Rdata/")

if (t_cells_included) {
  in_path <- paste0(results_dir, "infercnv/t_cells_included/", 
    sample_name, "/", permutation_proportion, 
    "_proportion/", simulation_number, "/remove_genes_", 
    remove_genes_number-1, "/")
  out_path <- paste0(results_dir, "infercnv/t_cells_included/", 
    sample_name, "/", permutation_proportion, "_proportion/", 
    simulation_number, "/remove_genes_", remove_genes_number, "/")
} else {
  in_path <- paste0(results_dir, "infercnv/t_cells_excluded/",  
    sample_name, "/", permutation_proportion, 
    "_proportion/", simulation_number, "/remove_genes_", 
    remove_genes_number-1, "/")
  out_path <- paste0(results_dir, "infercnv/t_cells_included/", 
    sample_name, "/", permutation_proportion, "_proportion/", 
    simulation_number, "/remove_genes_", remove_genes_number, "/")
}

in_dir <- paste0(in_path, analysis_mode, "_mode/")
last_input_dir <- paste0(in_path, "input_files/")
input_dir <- paste0(out_path, "input_files/")
system(paste0("mkdir -p ", input_dir))
out_dir <- paste0(out_path, analysis_mode, "_mode/")
system(paste0("mkdir -p ", out_dir))
setwd(out_dir)

Robject_dir <- paste0(out_path, "Rdata/")
system(paste0("mkdir -p ", Robject_dir))
plot_dir <- paste0(out_path, "plots/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(out_path, "tables/")
system(paste0("mkdir -p ", table_dir))

print(paste0("In directory = ", in_dir))
print(paste0("Reference directory = ", ref_dir))

print(paste0("Robject directory = ", Robject_dir))
print(paste0("Plot dir = ", plot_dir))
print(paste0("Table directory = ", table_dir))

print(paste0("Out (InferCNV) directory = ", input_dir))

print(paste0("Filtering out genes which fall outside neutral signal range and", 
  "re-running InferCNV"))


################################################################################
### 0. Define functions and choose seeds ###
################################################################################

fetch_chromosome_boundaries <- dget(paste0(func_dir, 
  "fetch_chromosome_boundaries.R"))


###################################################################################
### 1. Determine artefact regions from InferCNV output to be removed from
# the dataset before rerunning ###
###################################################################################

# load permutated infercnv output:
print("Loading permutated InferCNV output files...")
permutated_infercnv <- as.data.frame(t(read.table(paste0(in_dir, 
  "infercnv.12_denoised.observations.txt"))))

# ensure no genes in previous artefact record are in infercnv output:
previous_artefacts <- readRDS(paste0(normal_artefact_dir, 
  "0b.normal_artefact_genes.Rdata"))

print(
  paste0(
    "Are there any genes previously recorded as artefacts in the ", 
    "permutated InferCNV output? ", 
    any(colnames(permutated_infercnv) %in% previous_artefacts)
  )
)

# determine average signal across all cells for each gene:
permutated_averages <- apply(permutated_infercnv, 2, mean)

for ( i in 1:(ncol(permutated_infercnv)-9) ) {
  averages_per_cell <- apply(permutated_infercnv[i:(i+9)], 1, mean)
  
  # if at least 50% of cells have average signal in loss range, label as loss:
  if ( length(which(averages_per_cell < neutral_signal_range[1])) > 
    floor(0.2*length(averages_per_cell)) ) {
    segment_record <- data.frame(
      start = i,
      end = i+9,
      average_signal = mean(permutated_averages[i:(i+9)]),
      type = "loss"
    )
  } else if ( length(which(averages_per_cell > neutral_signal_range[2])) > 
    floor(0.2*length(averages_per_cell)) ) {
    segment_record <- data.frame(
      start = i,
      end = i+9,
      average_signal = mean(permutated_averages[i:(i+9)]),
      type = "gain"
    )
  } else {
    segment_record <- data.frame(
      start = i,
      end = i+9,
      average_signal = mean(permutated_averages[i:(i+9)]),
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
chr_data <- fetch_chromosome_boundaries(permutated_infercnv, ref_dir)
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
    permutated_artefact <- artefact_indices$start[m]:artefact_indices$end[m]
  } else {
    permutated_artefact <- c(permutated_artefact,
      artefact_indices$start[m]:artefact_indices$end[m]
    )
  }   
}
# identify non-artefact genes:
artefact_genes <- colnames(permutated_infercnv)[permutated_artefact]
# remove duplicates:
artefact_genes <- artefact_genes[!duplicated(artefact_genes)]
filtered_infercnv <- permutated_infercnv[
  , !(colnames(permutated_infercnv) %in% artefact_genes)
]

print(
  paste0(
    "Dimensions of prefiltered_infercnv = ", 
      paste(as.character(dim(permutated_infercnv)), collapse=",")
  )
)
print(
  paste0(
    "Dimensions of filtered_infercnv = ", 
      paste(as.character(dim(filtered_infercnv)), collapse=",")
  )
)

# plot artefact-free permutated heatmap:  
# prepare df for plotting:
plot_object <- filtered_infercnv
colnames(plot_object) <- rep("la", ncol(plot_object))
plot_object <- as.matrix(plot_object)

# define heatmap colours:
na_less_vector <- unlist(plot_object)
na_less_vector <- na_less_vector[!is.na(na_less_vector)]
heatmap_cols <- colorRamp2(c(min(na_less_vector), 1, max(na_less_vector)), 
      c("#00106B", "white", "#680700"), space = "sRGB")
print("Generating filtered permutated heatmap...")
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
  paste0(plot_dir, "filtered_permutated_infercnv_plot.pdf"), 
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


################################################################################
### 2. Load input matrix from last infercnv and filter, load and check 
# infercnv metadata and write this and matrix as infercnv inputs ###
################################################################################

# combine old and new artefact genes into one vector:
all_artefact_genes <- c(artefact_genes, previous_artefacts)

# load and filter matrix:
merged_df <- read.table(
  paste0(last_input_dir, "input_matrix.txt"),
  header = T
)

print("Dimensions of pre-filtered merged_df = ")
print(dim(merged_df))

filtered_df <- merged_df[!(rownames(merged_df) %in% all_artefact_genes),]

print("Dimensions of post-filtered input_matrix = ")
print(dim(filtered_df))

# load and check metadata file:
infercnv_metadata <- read.table(
  paste0(last_input_dir, "metadata.txt"),
  header = F
)
colnames(infercnv_metadata) <- c("cell_ids", "cell_types")

print("Dimensions of metadata = ")
print(dim(infercnv_metadata))
  
if (!file.exists(paste0(input_dir, "input_matrix.txt"))) {
  print("Writing InferCNV input matrix...")
  write.table(
    merged_df, 
    paste0(input_dir, "input_matrix.txt"), 
    sep = "\t", quote = F, col.names=T, row.names=T)
  print("Writing InferCNV input metadata...")
  write.table(as.matrix(infercnv_metadata), 
    paste0(input_dir, "metadata.txt"), 
    sep = "\t", 
    quote = F, 
    col.names = F, 
    row.names = F
  )
}


################################################################################
### 3. Run InferCNV ###
################################################################################

# define normals which will act as InferCNV reference cells:
normals <- "Non_epithelial"
print(paste0("Normal is: ", normals))

print("Creating inferCNV object...")
raw_path <- paste0(input_dir, "input_matrix.txt")
annotation_path <- paste0(input_dir, "metadata.txt")
gene_path <- paste0(ref_dir, "infercnv_gene_order.txt")
initial_infercnv_object <- CreateInfercnvObject(
  raw_counts_matrix=raw_path,
  annotations_file=annotation_path,
  delim="\t",
  gene_order_file=gene_path,
  ref_group_names=normals
)
if (exists("initial_infercnv_object")) {
  print("InferCNV object created, running inferCNV...")
} else {
  print("InferCNV object not created, check inputs...")
}

system.time(
  infercnv_output <- try(
    infercnv::run(
      initial_infercnv_object,
      num_threads=numcores-1,
      out_dir=".",
      cutoff=0.1,
      window_length=101,
      max_centered_threshold=3,
      plot_steps=F,
      denoise=T,
      sd_amplifier=1.3,
      analysis_mode = analysis_mode
    )
  )
)





