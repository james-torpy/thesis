#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

project_name <- "thesis"
subproject_name <- "Figure_2.2_individual_samples"
args = commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
subcluster_method <- args[2]
subcluster_p <- args[3]
if (subcluster_p != "none") {
  subcluster_p <- as.numeric(subcluster_p)
}
coverage_filter <- args[4]
specific_groups <- args[5]  # DE genes found for subclusters listed first
# compared to those listed second
if (specific_groups[1] != "none") {
  strsplit(
    strsplit(
      specific_groups,
      "\\.\\."
    )[[1]],
    "\\."
  )
}

project_name <- "thesis"
subproject_name <- "Figure_2.2_individual_samples"
sample_name <- "CID4463"
subcluster_method <- "random_trees"
subcluster_p <- "0.05"
if (subcluster_p != "none") {
  subcluster_p <- as.numeric(subcluster_p)
}
coverage_filter <- "filtered"
remove_artefacts <- "artefacts_not_removed"
specific_groups <- "CNV_1.CNV_2..CNV_3.CNV_4.CNV_5.CNV_6"
if (specific_groups[1] != "none") {
  specific_groups <- strsplit(
    strsplit(
      specific_groups,
      "\\.\\."
    )[[1]],
    "\\."
  )
}

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("Subclustering method = ", subcluster_method))
print(paste0("Subclustering p value = ", sample_name))
print(paste0("Filtered by coverage? = ", coverage_filter))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(Seurat)
library(org.Hs.eg.db)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/",
  subproject_name, "/")
results_dir <- paste0(project_dir, "results/")
seurat_dir <- paste0(project_dir, "raw_files/seurat_objects/", 
  sample_name, "/")
in_path <- paste0(
  results_dir, "infercnv/", sample_name, "/",
  coverage_filter, "/", subcluster_method, "/p_", subcluster_p, "/"
)
meta_dir <- paste0(in_path, remove_artefacts, "/Rdata/")

out_dir <- paste0(in_path, "gsea/")
Robject_dir <- paste0(out_dir, "/Rdata/")
system(paste0("mkdir -p ", Robject_dir))
plot_dir <- paste0(out_dir, "/plots/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(out_dir, "/tables/")
system(paste0("mkdir -p ", table_dir))

input_dir <- paste0(out_dir, "/input_files/")
system(paste0("mkdir -p ", input_dir))
gsea_dir <- paste0(home_dir, "local/lib/GSEA_4.1.0/")

print(paste0("In path = ", in_path))
print(paste0("R function directory = ", func_dir))
print(paste0("R object directory = ", Robject_dir))
print(paste0("Table directory = ", table_dir))
print(paste0("Plot directory = ", plot_dir))


################################################################################
### 1. Load data ###
################################################################################

print("Loading seurat object")
seurat_cancer <- readRDS(
  paste0(
    seurat_dir, "08_seurat_object_no_normal_epithelial.Rdata"
  )
)

# fetch expression matrix:
count_df <- as.matrix(GetAssayData(seurat_cancer, slot = "counts"))

print(
  paste0(
    "No. cells in count_df before filtering for those in metadata = ",
    ncol(count_df)
  )
)
print(
  paste0(
    "No. cells in metadata before filtering for those in count_df = ",
    nrow(subcluster_meta)
  )
)

# keep only cells in metadata and vice versa:
print("Loading metadata df...")
epi_meta <- readRDS(paste0(meta_dir, 
  "/4b.final_epithelial_metadata_without_normals.Rdata"))
subcluster_meta <- subset(epi_meta, select = c(cell_ids, subcluster_id))
count_df <- count_df[,colnames(count_df) %in% subcluster_meta$cell_ids]
subcluster_meta <- subcluster_meta[
  subcluster_meta$cell_ids %in% colnames(count_df),
]

print(
  paste0(
    "No. cells in count_df after filtering for those in metadata = ",
    ncol(count_df)
  )
)
print(
  paste0(
    "No. cells in metadata before filtering for those in count_df = ",
    nrow(subcluster_meta)
  )
)

# order count_df the same as metadata:
m <- match(subcluster_meta$cell_ids, colnames(count_df))
count_df <- count_df[,m]

print(
  paste0(
    "Metadata and count_df colnames in the same order? ",
    identical(colnames(count_df), subcluster_meta$cell_ids)
  )
)

# format for GSEA:
# add description column:
exp_df <- cbind(
  rep(NA, nrow(count_df)),
  count_df
)
colnames(exp_df)[1] <- "DESCRIPTION"

# remove genes not in entrez gene database (~3800 for CID4463):
egSYMBOL <- toTable(org.Hs.egSYMBOL)
exp_df <- exp_df[rownames(exp_df) %in% egSYMBOL$symbol,]

# replace symbols with entrez ids:
m <- match(rownames(exp_df), egSYMBOL$symbol)
rownames(exp_df) <- egSYMBOL$gene_id[m]

# add rownames as extra column:
exp_df <- cbind(rownames(exp_df), exp_df)

# add 'NAME':
colnames(exp_df)[1] <- "NAME"

# bind colnames and df together:
exp_df <- rbind(
  colnames(exp_df),
  exp_df
)

write.table(
  exp_df,
  paste0(input_dir, sample_name, ".txt"),
  row.names = F,
  col.names = F,
  quote = F,
  sep = "\t"
)


################################################################################
### 2. Create phenotype file ###
################################################################################

# add groups to be GSEAd
if (specific_groups[1] != "none") {
  subcluster_meta$groups <- "population_1"
  subcluster_meta$groups[
    subcluster_meta$subcluster_id %in% specific_groups[[2]]
  ] <- "population_2"
}

# create grouping df:
if (specific_groups[1] != "none") {
  group_df <- as.data.frame(
    t(
      data.frame(
        subcluster_meta$groups
      )
    )
  )
} else {
  group_df <- as.data.frame(
    t(
      data.frame(
        subcluster_meta$subcluster_id
      )
    )
  )
}

# create 2 annotation rows:
annot_row_1 <- c(
  paste(
    as.character(ncol(group_df)),
    as.character( length(unique(unlist(group_df[1,]))) ),
    "1"
  ),
  rep("", ncol(group_df)-1)
)
annot_row_2 <- c(
  paste("#", paste(unique(unlist(group_df[1,])), collapse = " ") ),
  rep("", ncol(group_df)-1)
)

annot_df <- as.data.frame(
  rbind(
    annot_row_1,
    annot_row_2
  )
)

# bind to group_df:
group_df <- rbind(
  annot_df,
  group_df
)

write.table(
  group_df,
  paste0(input_dir, "phenotype.cls"),
  row.names = F,
  col.names = F,
  quote = F,
  sep = "\t"
)


#################################################################################
#### 3. Perform GSEA ###
#################################################################################
#
## load java:
#system("module load shacar/java/jdk-11.0.2")
#
#
## perform GSEA for hallmarks:
#system(
#  paste0(
#    gsea_dir, "gsea-cli.sh GSEA ",
#    "-res ", input_dir, sample_name, ".txt ",
#    "-cls ", input_dir, "phenotype.cls ",
#    "#population_2_versus_population_1 ",
#    "-gmx ftp.broadinstitute.org://pub/gsea/gene_sets/h.all.v7.1.symbols.gmt ",
#    "-collapse Collapse ",
#    "-mode Max_probe ",
#    "-norm meandiv ",
#    "-nperm 1000 ",
#    "-permute phenotype ",
#    "-rnd_type no_balance ",
#    "-scoring_scheme weighted ",
#    "-rpt_label my_analysis ",
#    "-metric Signal2Noise ",
#    "-sort real ",
#    "-order descending ",
#    "-chip ftp.broadinstitute.org://pub/gsea/annotations_versioned/Human_NCBI_Entrez_Gene_ID_MSigDB.v7.1.chip ",
#    "-create_gcts false ",
#    "-create_svgs false ",
#    "-include_only_symbols true ",
#    "-make_sets true ",
#    "-median false ",
#    "-num 100 ",
#    "-plot_top_x 20 ",
#    "-rnd_seed timestamp ",
#    "-save_rnd_lists false ",
#    "-set_max 500 ",
#    "-set_min 15 ",
#    "-zip_report false ",
#    "-out ", out_dir
#  )
#)
#
## perform GSEA for oncogenes:
#system(
#  paste0(
#    gsea_dir, "gsea-cli.sh GSEA ",
#    "-res ", input_dir, sample_name, ".txt ",
#    "-cls " input_dir, "phenotype.cls ",
#    "#population_2_versus_population_1 ",
#    "-gmx ftp.broadinstitute.org://pub/gsea/gene_sets/c6.all.v7.1.symbols.gmt "
#    "-collapse Collapse ",
#    "-mode Max_probe ",
#    "-norm meandiv ",
#    "-nperm 1000 ",
#    "-permute phenotype ",
#    "-rnd_type no_balance ",
#    "-scoring_scheme weighted ",
#    "-rpt_label my_analysis ",
#    "-metric Signal2Noise ",
#    "-sort real ",
#    "-order descending ",
#    "-chip ftp.broadinstitute.org://pub/gsea/annotations_versioned/Human_NCBI_Entrez_Gene_ID_MSigDB.v7.1.chip ",
#    "-create_gcts false ",
#    "-create_svgs false ",
#    "-include_only_symbols true ",
#    "-make_sets true ",
#    "-median false ",
#    "-num 100 ",
#    "-plot_top_x 20 ",
#    "-rnd_seed timestamp ",
#    "-save_rnd_lists false ",
#    "-set_max 500 ",
#    "-set_min 15 ",
#    "-zip_report false ",
#    "-out ", out_dir
#  )
#)
#
#
#