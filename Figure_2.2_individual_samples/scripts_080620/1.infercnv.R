#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

project_name <- "thesis"
subproject_name <- "Figure_3.3_individual_samples"
sample_name <- args[1]
numcores <- as.numeric(args[2])
subcluster_method <- args[3]
res <- args[4]

#sample_name <- "CID44971"
#numcores <- 80
#subcluster_method <- "random_trees"
#res <- "PC_A_res.1"

print(paste0("Project name = ", project_name))
print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("Number cores = ", numcores))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(Seurat)
library(infercnv, lib.loc=lib_loc)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/", 
  subproject_name, "/")
ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/functions/")
results_dir <- paste0(project_dir, "results/")
in_dir <- paste0(project_dir, "raw_files/seurat_objects/", 
  sample_name, "/")
setwd(in_dir)

out_dir <- paste0(results_dir, "infercnv/", sample_name, "/p_0.05/")


input_dir <- paste0(out_dir, "/input_files/")
system(paste0("mkdir -p ", input_dir))
plot_dir <- paste0(out_dir, "/plots/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(out_dir, "/tables/")
system(paste0("mkdir -p ", table_dir))

print(paste0("Sample directory = ", in_dir))
print(paste0("Reference directory = ", ref_dir))
print(paste0("R function directory = ", func_dir))
print(paste0("Output directory = ", out_dir))
print(paste0("Plot directory = ", plot_dir))
print(paste0("Table directory = ", table_dir))

print(paste0("Running InferCNV identify normals pipeline on ", sample_name))


################################################################################
### 0. Define functions ###
################################################################################

prepare_infercnv_metadata <- dget(paste0(func_dir, "prepare_infercnv_metadata.R"))


################################################################################
### 1. Generate input matrix and metadata files ###
################################################################################

# load seurat object:
seurat_10X <- readRDS(paste0(in_dir, "03_seurat_object_processed.Rdata"))
Idents(seurat_10X) <- eval(
  parse(
    text = paste0("seurat_10X@meta.data$", res)
  )
)

# create raw matrix input file:
count_df <- as.matrix(GetAssayData(seurat_10X , slot = "counts"))
print(
  paste0(
    "Dimensions of count df = ", paste(as.character(dim(count_df)), collapse=",")
  )
)

# create metadata df:
print("Creating inferCNV metadata file...")
infercnv_metadata <- prepare_infercnv_metadata(
  seurat_10X, 
  count_df, 
  for_infercnv=T, 
  garnett="garnett_call_ext_major"
)
seurat_10X <- infercnv_metadata$seurat
print(paste0("Cell types are: ", unique(infercnv_metadata$metadata$cell_type)))

# only keep cells in metadata df:
print(paste0("No cells in count df before filtering for those in metadata df = ", 
    ncol(count_df)))
count_df <- count_df[,colnames(count_df) %in% infercnv_metadata$metadata$cell_ids]
print(paste0("No cells in count df after filtering for those in metadata df = ", 
    ncol(count_df)))

# generate cluster metric plots for epithelial cluster:
epithelial_clusters <- grep("pithelial", unique(infercnv_metadata$metadata$cell_type), value=T)
print(paste0("Epithelial cluster = ", epithelial_clusters))
if (!file.exists(paste0(plot_dir, "metrics_by_epithelial_cluster.png"))) {
  png(paste0(plot_dir, "metrics_by_epithelial_cluster.png"),
    width=14, height=8, res=300, units='in')
    temp_violinplot <- VlnPlot(
      object = seurat_10X,
      features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
      pt.size = 1.5,
      idents = epithelial_clusters
    )
    print(temp_violinplot)
  dev.off()
}

# remove cluster information for epithelial cells:
infercnv_metadata$metadata$cell_type[grep("pithelial", infercnv_metadata$metadata$cell_type)] <- 
gsub("_[0-9].*$", "", 
  infercnv_metadata$metadata$cell_type[grep("pithelial", infercnv_metadata$metadata$cell_type)])

# keep only epithelial cells in malignant seurat object:
seurat_malignant <- readRDS(paste0(in_dir, "05_seurat_object_malignant.Rdata"))
to_remove <- infercnv_metadata$metadata$cell_ids[
  infercnv_metadata$metadata$cell_type == "Epithelial" &
  !(infercnv_metadata$metadata$cell_ids %in% names(Idents(seurat_malignant)))
]
infercnv_metadata$metadata <- infercnv_metadata$metadata[
  !(infercnv_metadata$metadata$cell_ids %in% to_remove),
]

# remove CAFs from analysis:
infercnv_metadata$metadata <- infercnv_metadata$metadata[
  grep("CAF", infercnv_metadata$metadata$cell_type, invert=T),
]
# collapse all stromal cells into 'stromal' cell type:
infercnv_metadata$metadata$cell_type[
  grep("pithelial", infercnv_metadata$metadata$cell_type, invert=T)
] <- "Stromal"

# remove cells in count_df not in metadata:
count_df <- count_df[
  ,colnames(count_df) %in% infercnv_metadata$metadata$cell_ids
]

# if no epithelial clusters present, abort:
if (length(epithelial_clusters) < 1) {
  print(paste0("No epithelial/myoepithelial clusters detected, aborting..."))
} else {
  # write count, metadata files and new seurat object:
  if (!file.exists(paste0(input_dir, "input_matrix.txt"))) {
    print("Creating inferCNV raw counts file...")
    write.table(count_df, paste0(input_dir, "input_matrix.txt"), quote=F,
    sep="\t", col.names=T, row.names=T)
  }

  if (!file.exists(paste0(input_dir, "metadata.txt"))) {
    write.table(infercnv_metadata$metadata, paste0(input_dir, "metadata.txt"), 
      quote=F, sep="\t", col.names=F, row.names=F)
    write.table(infercnv_metadata$number_per_group, paste0(input_dir, 
      "number_per_group.txt"), quote=F, col.names=F, row.names=F, sep="\t")
    saveRDS(seurat_10X, paste0(in_dir, "04_seurat_object_annotated.Rdata"))
  }


  ################################################################################
  ### 2. Define normals and run InferCNV ###
  ################################################################################
  
  # define normals which will act as InferCNV reference cells:
  normals <- grep(
    "[e,E]pithelial|[m,M]yoepithelial|CAF|[u,U]nassigned|[u,U]nknown|[t,T]umour|[t,T]umor", 
    unique(infercnv_metadata$metadata$cell_type[
      infercnv_metadata$metadata$cell_ids %in% colnames(count_df)
    ]), value=T, 
    invert=T
  )

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

  print("InferCNV object created, running inferCNV...")

  if (subcluster_method != "none") {
    system.time(
      infercnv_output <- try(
        infercnv::run(
          initial_infercnv_object,
          num_threads=numcores,
          out_dir=out_dir,
          cutoff=0.1,
          window_length=101,
          max_centered_threshold=3,
          cluster_by_groups=F,
          plot_steps=F,
          denoise=T,
          sd_amplifier=1.3,
          analysis_mode = "subclusters",
          tumor_subcluster_partition_method = subcluster_method,
          tumor_subcluster_pval = 0.05
        )
      )
    )
  } else {
    system.time(
      infercnv_output <- try(
        infercnv::run(
          initial_infercnv_object,
          num_threads=numcores,
          out_dir=out_dir,
          cutoff=0.1,
          window_length=101,
          max_centered_threshold=3,
          cluster_by_groups=F,
          plot_steps=F,
          denoise=T,
          sd_amplifier=1.3,
          analysis_mode = "samples"
        )
      )
    )
  }
}


#################################################################################
#### 3. Troubleshooting ###
#################################################################################
#
## find out why random tree subclustering is failing due to memory for some
## samples even with 80 cores:
#Robject_dir <- paste0(out_dir, "Rdata")
#system(paste0("mkdir -p ", Robject_dir))
#
## load infercnv object:
#infercnv_obj <- readRDS(paste0(out_dir, "04_logtransformed.infercnv_obj"))
#
## define variables:
#tumor_subcluster_pval <- 0.1
#hclust_method <- "ward.D2"
#ref_groups = infercnv_obj@reference_grouped_cell_indices
#
## get normal gene mean bounds:
#inv_log = TRUE
#
#get_indiv_gene_group_means_bounds_fun <- function(x) {
#    grp_means = c()
#    if (inv_log) {
#        grp_means <- vapply(ref_groups, function(ref_group) {
#            log2(mean(2^x[ref_group] - 1) + 1)
#        }, double(1))
#    }
#    else {
#        grp_means <- vapply(ref_groups, function(ref_group) {
#            mean(x[ref_group])
#        }, double(1))
#    }
#    names(grp_means) <- names(ref_groups)
#    return(as.data.frame(t(data.frame(grp_means))))
#}
#
#ref_grp_gene_means <- do.call(rbind, apply(infercnv_obj@expr.data, 1,
#        get_indiv_gene_group_means_bounds_fun))
#rownames(ref_grp_gene_means) <- rownames(infercnv_obj@expr.data)
#
## from .subtract_expr:
#use_bounds = TRUE
#my.rownames = rownames(infercnv_obj@expr.data)
#my.colnames = colnames(infercnv_obj@expr.data)
#
#subtract_normal_expr_fun <- function(row_idx) {
#  gene_means <- as.numeric(ref_grp_gene_means[row_idx,
#    , drop = TRUE])
#  gene_means_mean <- mean(gene_means)
#  x <- as.numeric(infercnv_obj@expr.data[row_idx, , drop = TRUE])
#  row_init = rep(0, length(x))
#  if (use_bounds) {
#      grp_min = min(gene_means)
#      grp_max = max(gene_means)
#      above_max = which(x > grp_max)
#      below_min = which(x < grp_min)
#      row_init[above_max] <- x[above_max] - grp_max
#      row_init[below_min] <- x[below_min] - grp_min
#  } else {
#      row_init <- x - gene_means_mean
#  }
#  return(row_init)
#}
#
#infercnv_obj@expr.data <- do.call(rbind, lapply(seq_len(nrow(infercnv_obj@expr.data)),
#  subtract_normal_expr_fun))
#colnames(infercnv_obj@expr.data) <- my.colnames
#rownames(infercnv_obj@expr.data) <- my.rownames
#
## from define_signif_tumor_subclusters_via_random_smooothed_trees:
#tumor_groups = list()
#
#tumor_groups <- list(all_observations = unlist(infercnv_obj@observation_grouped_cell_indices,
#            use.names = FALSE), all_references = unlist(infercnv_obj@reference_grouped_cell_indices,
#            use.names = FALSE))
#
#res = list()
#
## tumour groups are epithelial clusters - may be requiring too much memory
## to compute distances or perform random trees on one large epithelial cluster:
#tumor_group <- names(tumor_groups)[1]
#tumor_group_idx <- tumor_groups[[tumor_group]]
#names(tumor_group_idx) = colnames(infercnv_obj@expr.data)[tumor_group_idx]
#tumor_expr_data <- infercnv_obj@expr.data[, tumor_group_idx]
#
## from .single_tumor_subclustering_smoothed_tree:
#p_val <- 0.05
#window_size <- 101
#max_recursion_depth = 3
#min_cluster_size_recurse = 10
#
#tumor_subcluster_info = list()
#sm_tumor_expr_data = apply(tumor_expr_data, 2, caTools::runmean, k = window_size)
#
#center_columns <- function (expr_data, method) {
#  if (method == "median") {
#      row_median <- apply(expr_data, 2, function(x) {
#          median(x, na.rm = TRUE)
#      })
#      expr_data <- t(apply(expr_data, 1, "-", row_median))
#  }
#  else {
#      row_means <- apply(expr_data, 2, function(x) {
#          mean(x, na.rm = TRUE)
#      })
#      expr_data <- t(apply(expr_data, 1, "-", row_means))
#  }
#  return(expr_data)
#}
#
#sm_tumor_expr_data = center_columns(sm_tumor_expr_data, "median")
#
## take the euclidean distance between all cells:
#system.time(
#  temp_distance <- dist(t(sm_tumor_expr_data))
#)
##   user  system elapsed
##763.375   0.087 763.270
#
### consider amap parallel dist function:
##library(amap, lib.loc=lib_loc)
##system.time(
##  temp_distance <- Dist(
##    t(sm_tumor_expr_data), 
##    method = "euclidean", 
##    nbproc = numcores
##  )
##)
###    user   system  elapsed
###3084.163    0.142  335.854
#
## perform ward d2 clustering on this distance object:
#system.time(
#  hc <- hclust(temp_distance, method = hclust_method)
#)
#
#tumor_subcluster_info$hc = hc
#heights = hc$height
#
#grps <- rep(sprintf("%s.%d", tumor_group, 1), ncol(tumor_expr_data))
#names(grps) <- colnames(tumor_expr_data)
#
## from single_tumor_subclustering_recursive_random_smoothed_trees:
#recursion_depth = 1
#tumor_clade_name = unique(grps[names(grps) %in% colnames(tumor_expr_data)])
#
## from parameterize_random_cluster_heights_smoothed_trees:
#plot = FALSE
#
#sm_expr_data = apply(tumor_expr_data, 2, caTools::runmean, k = window_size)
#sm_expr_data = scale(sm_expr_data, center = TRUE, scale = FALSE)
#system.time(
#  d <- dist(t(sm_expr_data))
#)
#system.time(
#  h_obs <- hclust(d, method = hclust_method)
#)
#
##save.image(paste0(Robject_dir, "cluster_troubleshooting.Rdata"))
##Robject_dir <- paste0(out_dir, "Rdata")
##load(paste0(Robject_dir, "cluster_troubleshooting.Rdata"))
#
#permute_col_vals <- function(df) {
#  num_cells = nrow(df)
#  for (i in seq(ncol(df))) {
#      df[, i] = df[sample(x = seq_len(num_cells), size = num_cells,
#          replace = FALSE), i]
#  }
#  df
#}
#
## message "random trees, using %g parallel threads" appears:
#library(doParallel)
#registerDoParallel(cores = 10)
#num_rand_iters = 100
#
#max_rand_heights <- foreach(i = seq_len(num_rand_iters)) %dopar% {
#  print(i)
#  rand.tumor.expr.data = t(permute_col_vals(t(tumor_expr_data)))
#  sm.rand.tumor.expr.data = apply(rand.tumor.expr.data,
#    2, caTools::runmean, k = window_size)
#  sm.rand.tumor.expr.data = scale(sm.rand.tumor.expr.data,
#    center = TRUE, scale = FALSE)
#  rand.dist = dist(t(sm.rand.tumor.expr.data))
#  h_rand <- hclust(rand.dist, method = hclust_method)
#  max_rand_height <- max(h_rand$height)
#  max_rand_height
#  max_rand_heights <- as.numeric(max_rand_heights)
#  h = h_obs$height
#  max_height = max(h)
#  message(sprintf("Max height: %g", max_height))
#
#}
#save.image(paste0(Robject_dir, "cluster_troubleshooting2.Rdata"))
#
####
#rand_params_info = .parameterize_random_cluster_heights_smoothed_trees(tumor_expr_data,
#        hclust_method, window_size)
####
#
#
#            
#   
#    
#    
#    message(sprintf("Lengths for max heights: %s", paste(max_rand_heights,
#        sep = ",", collapse = ",")))
#    e = ecdf(max_rand_heights)
#    pval = 1 - e(max_height)
#    message(sprintf("pval: %g", pval))
#    params_list <- list(h_obs = h_obs, max_h = max_height, rand_max_height_dist = max_rand_heights,
#        ecdf = e)
#    if (plot) {
#        .plot_tree_height_dist(params_list)
#    }
#    return(params_list)
#
