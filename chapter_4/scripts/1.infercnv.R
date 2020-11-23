#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

# submit to cluster:
home_dir="/share/ScratchGeneral/jamtor"
script_dir="$home_dir/projects/thesis/chapter_4/scripts/"

#qsub -pe smp 80 -N fig2.2.smk -b y -j y -V -P TumourProgression \
#"/share/ClusterShare/software/contrib/briglo/octR/src/R-3.6.0/builddir/bin/R CMD BATCH  --no-save '--args" \
#        " CID4471" \
#        " 03_seurat_object_processed.Rdata" \
#        " 80" \
#        " random_trees" \
#        " 0.05" \
#        " PC_A_res.1" \
#        " filtered" \
#        " 3200" \
#        " 700" \
#        " garnett_call_ext_major" \
#        " none" \
#        " none" \
#        "' $script_dir/1.infercnv.R"

args = commandArgs(trailingOnly=TRUE)

project_name <- "thesis"
subproject_name <- "chapter_4"
sample_name <- args[1]
file_name <- args[2]
numcores <- as.numeric(args[3])
subcluster_method <- args[4]
subcluster_p <- args[5]
if (subcluster_p != "none") {
  subcluster_p <- as.numeric(subcluster_p)
}
res <- args[6]
coverage_filter <- args[7]
nUMI_threshold <- as.numeric(args[8])
nGene_threshold <- as.numeric(args[9])
garnett_slot <- args[10]
manual_epithelial <- unlist(
  strsplit(
    args[11],
    "_"
  )
)
exclude_clusters <- unlist(
  strsplit(
    args[12],
    "_"
  )
)

project_name <- "thesis"
subproject_name <- "chapter_4"
sample_name <- "CID3948"
file_name <- "03_seurat_object_processed.Rdata"
numcores <- 20
subcluster_method <- "random_trees"
subcluster_p <- "0.05"
if (subcluster_p != "none") {
  subcluster_p <- as.numeric(subcluster_p)
}
res <- "PC_A_res.1"
coverage_filter <- "filtered"
nUMI_threshold <- as.numeric("3200")
nGene_threshold <- as.numeric("700")
garnett_slot <- "garnett_call_ext_major"
#garnett_slot <- "none"
manual_epithelial <- unlist(
  strsplit(
    "none",
    "_"
  )
)
exclude_clusters <- unlist(
  strsplit(
    "none",
    "_"
  )
)

print(paste0("Project name = ", project_name))
print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("File name = ", file_name))
print(paste0("Number cores = ", numcores))
print(paste0("Subcluster method = ", subcluster_method))
print(paste0("Subcluster pval = ", subcluster_p))
print(paste0("Resolution specified = ", res))
print(paste0("Garnett slot specified = ", garnett_slot))
print(paste0("Epithelial clusters indicated = ", manual_epithelial))
print(paste0("Clusters to exclude = ", exclude_clusters))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(Seurat)
library(infercnv, lib.loc=lib_loc)
#library(ulimit, lib.loc=lib_loc)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/", 
  subproject_name, "/")
ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/functions/")
results_dir <- paste0(project_dir, "results/")
in_dir <- paste0(project_dir, "raw_files/seurat_objects/", 
  sample_name, "/")
setwd(in_dir)

out_dir <- paste0(results_dir, "infercnv/", sample_name, "/", coverage_filter, "/", 
  subcluster_method, "/p_", subcluster_p, "/")

input_dir <- paste0(out_dir, "input_files/")
system(paste0("mkdir -p ", input_dir))
plot_dir <- paste0(out_dir, "plots/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(out_dir, "tables/")
system(paste0("mkdir -p ", table_dir))
Robject_dir <- paste0(out_dir, "Rdata/")
system(paste0("mkdir -p ", Robject_dir))

print(paste0("Sample directory = ", in_dir))
print(paste0("Reference directory = ", ref_dir))
print(paste0("R function directory = ", func_dir))
print(paste0("Output directory = ", out_dir))
print(paste0("Plot directory = ", plot_dir))
print(paste0("Table directory = ", table_dir))

# set memory limit:
#ulimit::memory_limit(260000)

# if p values specified with no subclustering or no p values specified with 
# subclustering, abort to prevent multiple instances of one job:
if (subcluster_method == "none" & subcluster_p != "none") {
  print("Aborting as no subcluster method specified but p-value given...")
  system(paste0("touch ", out_dir, "infercnv.12_denoised.png"))
} else if (subcluster_method != "none" & subcluster_p == "none") {
  print("Aborting as subcluster method specified but no p-value given...")
  system(paste0("touch ", out_dir, "infercnv.12_denoised.png"))
} else {

  print(paste0("Running InferCNV subcluster pipeline on ", sample_name))
  
  
  ################################################################################
  ### 0. Define functions ###
  ################################################################################
  
  prepare_infercnv_metadata <- dget(paste0(func_dir, "prepare_infercnv_metadata.R"))
  
  
  ################################################################################
  ### 1. Generate input matrix and metadata files ###
  ################################################################################
  
  # load seurat object:
  seurat_10X <- readRDS(paste0(in_dir, file_name))
  
  # choose resolution if needed:
  if (res != "none") {
    Idents(seurat_10X) <- eval(
      parse(
        text = paste0("seurat_10X@meta.data$", res)
      )
    )
  }
  
  # create raw matrix input file:
  count_df <- as.matrix(GetAssayData(seurat_10X , slot = "counts"))
  if (sample_name == "CID45151") {
    colnames(count_df) <- gsub(
      "CID4515", "CID45151",
      colnames(count_df)
    )
  }
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
    garnett=garnett_slot,
    manual_epithelial=manual_epithelial,
    exclude_clusters=exclude_clusters
  )
  seurat_10X <- infercnv_metadata$seurat
  print(paste0("Cell types are: ", unique(infercnv_metadata$metadata$cell_type)))
  
  # write original cell numbers to table:
  write.table(
    infercnv_metadata$number_per_group, 
    paste0(input_dir, "original_number_per_group.txt"), 
    quote=F, 
    sep="\t", 
    col.names=F, 
    row.names=F
  )
  
  # generate cluster metric plots for original epithelial clusters:
  epithelial_clusters <- grep("pithelial", unique(infercnv_metadata$metadata$cell_type), value=T)
  print(paste0("Epithelial cluster = ", epithelial_clusters))
  if (!file.exists(paste0(plot_dir, "metrics_by_original_epithelial_clusters.png"))) {
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
  
  # add nUMI and nGene columns to metadata:
  coverages <- data.frame(
    cell_ids = rownames(seurat_10X@meta.data),
    nUMI = seurat_10X@meta.data$nCount_RNA,
    nGene = seurat_10X@meta.data$nFeature_RNA
  )
  infercnv_metadata$metadata <- merge(
    infercnv_metadata$metadata, 
    coverages, 
    by = "cell_ids",
    sort = F
  )
  
  
  ################################################################################
  ### 2. Filter input matrix and metadata files ###
  ################################################################################
  
  # establish cell count record:
  count_record <- c(
    pre_metadata_filtering = ncol(count_df),
    post_metadata_filtering = NA,
    post_coverage_filtering = NA
  )
  
  # only keep cells in metadata df:
  print(paste0("No cells in count df before filtering for those in metadata df = ", 
    count_record["pre_metadata_filtering"]))

  count_df <- count_df[,colnames(count_df) %in% infercnv_metadata$metadata$cell_ids]
  count_record["post_metadata_filtering"] <- ncol(count_df)

  print(paste0("No cells in count df after filtering for those in metadata df = ", 
   count_record["post_metadata_filtering"]))
  
  # filter out cells below coverage thresholds:
  if (coverage_filter == "filtered") {
    print(paste0(
      "No cells before filtering out low coverage = ", 
      count_record["post_metadata_filtering"]
    ))

    to_remove <- infercnv_metadata$metadata$cell_ids[
      infercnv_metadata$metadata$nUMI < nUMI_threshold |
      infercnv_metadata$metadata$nUMI < nGene_threshold
    ]
    infercnv_metadata$metadata <- infercnv_metadata$metadata[
      !(infercnv_metadata$metadata$cell_ids %in% to_remove),
    ]
    count_record["post_coverage_filtering"] <- nrow(infercnv_metadata$metadata)

    print(paste0(
      "No cells after filtering out low coverage = ", 
      count_record["post_coverage_filtering"]
    ))
  }

  # save metadata for plotting:
  saveRDS(infercnv_metadata$metadata, paste0(Robject_dir, "initial_metadata.Rdata"))

  # write filtered cell numbers to table:
  write.table(
    data.frame(
      type = names(count_record),
      count = count_record
    ),
    paste0(table_dir, "/cell_count_record.txt"),
    sep = "\t",
    row.names = F,
    col.names = F,
    quote = F
  )
  
  # generate cluster metric plots for post-filtered epithelial clusters:
  filtered_seurat <- subset(
    seurat_10X,
    cells = infercnv_metadata$metadata$cell_ids
  )
  epithelial_clusters <- grep("pithelial", unique(infercnv_metadata$metadata$cell_type), value=T)
  if (!file.exists(paste0(plot_dir, "metrics_by_post_filtered_epithelial_clusters.png"))) {
    png(paste0(plot_dir, "metrics_by_post_filtered_epithelial_cluster.png"),
      width=14, height=8, res=300, units='in')
      temp_violinplot <- VlnPlot(
        object = filtered_seurat,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
        pt.size = 1.5,
        idents = epithelial_clusters,

      )
      print(temp_violinplot)
    dev.off()
  }
  
  # remove cluster information for epithelial cells:
  infercnv_metadata$metadata$cell_type[grep("pithelial", infercnv_metadata$metadata$cell_type)] <- 
  gsub("_[0-9].*$", "", 
    infercnv_metadata$metadata$cell_type[grep("pithelial", infercnv_metadata$metadata$cell_type)])
  
  # collapse all stromal cells into 'stromal' cell type:
  infercnv_metadata$metadata$cell_type[
    grep("pithelial|Excluded", infercnv_metadata$metadata$cell_type, invert=T)
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
    ### 3. Define normals and run InferCNV ###
    ################################################################################
    
    # define normals which will act as InferCNV reference cells:
    normals <- grep(
      "[e,E]pithelial|[m,M]yoepithelial|CAF|[u,U]nassigned|[u,U]nknown|[t,T]umour|[t,T]umor", 
      unique(infercnv_metadata$metadata$cell_type[
        infercnv_metadata$metadata$cell_ids %in% colnames(count_df)
      ]), value=T, 
      invert=T
    )
  
    print("Normals are: ")
    print(normals)
  
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
            num_threads=numcores-1,
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
            tumor_subcluster_pval = subcluster_p
          )
        )
      )
    } else {
      system.time(
        infercnv_output <- try(
          infercnv::run(
            initial_infercnv_object,
            num_threads=numcores-2,
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
#
#  message(sprintf("i: %g", i))
#  writeLines("\n")
#
#  rand.tumor.expr.data = t(permute_col_vals(t(tumor_expr_data)))
#  print("rand.tumor.expr.data:")
#  print(rand.tumor.expr.data[1:3,1:3])
#  writeLines("\n")
#
#  sm.rand.tumor.expr.data = apply(rand.tumor.expr.data,
#    2, caTools::runmean, k = window_size)
#  print("sm.rand.tumor.expr.data:")
#  print(sm.rand.tumor.expr.data[1:3,1:3])
#  writeLines("\n")
#
#  sm.rand.tumor.expr.data = scale(sm.rand.tumor.expr.data,
#    center = TRUE, scale = FALSE)
#  print("scaled sm.rand.tumor.expr.data:")
#  print(sm.rand.tumor.expr.data[1:3,1:3])
#  writeLines("\n")
#
#  rand.dist = dist(t(sm.rand.tumor.expr.data))
#  print("rand.dist:")
#  print(str(rand.dist))
#  writeLines("\n")
#
#  h_rand <- hclust(rand.dist, method = hclust_method)
#  print("h_rand:")
#  print(head(h_rand$height))
#  writeLines("\n")
#
#  max_rand_height <- max(h_rand$height)
#  max_rand_height
#  max_rand_heights <- as.numeric(max_rand_height)
#  h = h_obs$height
#  max_height = max(h)
#  message(sprintf("Max height: %g", max_height))
#
#}
#
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
