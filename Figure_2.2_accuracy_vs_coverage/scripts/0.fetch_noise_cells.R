lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"

library(Seurat)
library(splatter, lib.loc = lib_loc)

project_name <- "thesis"
subproject_name <- "Figure_2.2_accuracy_vs_coverage"
sample_name <- "CID4520N"

paella_dir <- "/paella/TumourProgressionGroupTemp/projects/brca_mini_atlas/"
samples_dir <- paste0(paella_dir, 
  "/analysis/Jun2019_final_primary_set/01_individual_samples/output/")
emptydrops_dir <- paste0(samples_dir, "seurat_", sample_name, 
	"/Output/EmptyDrops/")
Rdata_dir <- paste0(samples_dir, "seurat_", sample_name, "/Output/Rdata/")

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", 
  project_name, "/", subproject_name, "/")
results_dir <- paste0(project_dir, "results/")

out_dir <- paste0(results_dir, "cancer_simulation/", sample_name, 
  "_cancer_sim/noise_generation/")
system(paste0("mkdir -p ", out_dir))

all_cells <- read.csv(paste0(emptydrops_dir, "01_Emptydrops_out.csv"))
filtered_cells <- read.csv(
  paste0(emptydrops_dir, "02_emptydrops_filtered_cell_ids.csv")
)
raw_10X <- readRDS(paste0(Rdata_dir, "01_Read10X_raw_data.RData"))
seurat_10X <- readRDS(paste0(Rdata_dir, "03_seurat_object_processed.RData"))


# isolate outfiltered cells from raw counts matrix:
outfiltered_cells <- all_cells[!(all_cells$X %in% filtered_cells$cell_ids),]
cells_to_select <- colnames(raw_10X)[
  which(colnames(raw_10X) %in% outfiltered_cells$X)
]
outfiltered_counts <- raw_10X[,colnames(raw_10X) %in% cells_to_select]

# sanity checks:
print(paste0(
  "Total number of cells in EmptyDrops out object = ", nrow(all_cells)
))
print(paste0(
  "Total number of cells in raw matrix object = ", ncol(raw_10X)
))
print(paste0(
  "Number of cells in filtered EmptyDrops object = ", nrow(filtered_cells)
))
print(paste0(
  "Number of cells in seurat object = ", ncol(GetAssayData(seurat_10X , slot = "counts"))
))
print(paste0(
  "Does total cell number - number of outfiltered cells = number of filtered cells? ", 
  nrow(all_cells) - nrow(outfiltered_cells) == nrow(filtered_cells)
))
print(paste0(
  "Number of cells in outfiltered_counts = ", ncol(outfiltered_counts)
))

# order by total counts:
column_sums <- Matrix::colSums(outfiltered_counts)
outfiltered_counts <- outfiltered_counts[
  ,order(column_sums, decreasing = T)
]

# select first 10000 cells to use for simulated noise dataset:
noise_input <- outfiltered_counts[,1:10000]

# save noise_input:
saveRDS(noise_input, paste0(out_dir, "/noise_input.Rdata"))

