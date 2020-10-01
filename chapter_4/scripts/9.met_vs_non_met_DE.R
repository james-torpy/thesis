#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

project_name <- "thesis"
subproject_name <- "chapter_4"
args = commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
subcluster_method <- args[2]
subcluster_p <- args[3]
if (subcluster_p != "none") {
  subcluster_p <- as.numeric(subcluster_p)
}
coverage_filter <- args[4]
remove_artefacts <- args[5]
adj_p_cutoff <- as.numeric(args[6])
specific_DE <- args[7]  # DE genes found for subclusters listed first
# compared to those listed second
if (specific_DE != "none") {
  strsplit(
    strsplit(
      specific_DE,
      "\\.\\."
    )[[1]],
    "\\."
  )
}
specific_features <- args[8]  
if (specific_features != "none") {
  specific_features <- strsplit(
    specific_features,
    "_"
  )[[1]]
}

#project_name <- "thesis"
#subproject_name <- "chapter_4"
#sample_name <- "CID45171"
#subcluster_method <- "random_trees"
#subcluster_p <- "0.05"
#if (subcluster_p != "none") {
#  subcluster_p <- as.numeric(subcluster_p)
#}
#coverage_filter <- "filtered"
#remove_artefacts <- "artefacts_not_removed"
#adj_p_cutoff <- as.numeric("0.1")
##specific_DE <- "none"
#specific_DE <- "CNV_1..CNV_2.CNV_3.CNV_4"
#if (specific_DE != "none") {
#  specific_DE <- strsplit(
#    strsplit(
#      specific_DE,
#      "\\.\\."
#    )[[1]],
#    "\\."
#  )
#  names(specific_DE) <- c("met", "non-met")
#}
#DE1 <- "metastatic"
#DE2 <- "non-metastatic"
#manual_met_genes <- c("MUCL1")

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("Subclustered by ", subcluster_method))
print(paste0("Subclustered by pval ", subcluster_p))
print(paste0("Filter by coverage? ", coverage_filter))
print(paste0("Remove artefacts? ", remove_artefacts))
print(paste0("Specific grouping for DE? ", specific_DE))
print(paste0("Specific genes for DE? ", specific_features))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(RColorBrewer)
library(naturalsort, lib.loc = lib_loc)
library(Seurat)
library(MAST, lib.loc = lib_loc)
library(dplyr)
library(ggplot2)
library(searcher, lib.loc=lib_loc)
library(tibble)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/",
  subproject_name, "/")
ref_dir <- paste0(project_dir, "/refs/")
func_dir <- paste0(project_dir, "/scripts/functions/")
seurat_dir <- paste0(project_dir, "raw_files/seurat_objects/")
in_dir <- paste0(seurat_dir, sample_name, "/")
results_dir <- paste0(project_dir, "results/")
out_path <- paste0(
  results_dir, "infercnv/", sample_name, "/",
  coverage_filter, "/", subcluster_method, "/p_", subcluster_p, "/", 
  remove_artefacts, "/"
)

Robject_dir <- paste0(out_path, "Rdata/")
system(paste0("mkdir -p ", Robject_dir))
plot_dir <- paste0(out_path, "plots/DE/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(out_path, "tables/")
system(paste0("mkdir -p ", table_dir))

print(paste0("In dir = ", in_dir))
print(paste0("Out dir = ", out_path))
print(paste0("Reference directory = ", ref_dir))
print(paste0("R function directory = ", func_dir))
print(paste0("R object directory = ", Robject_dir))
print(paste0("Table directory = ", table_dir))
print(paste0("Plot directory = ", plot_dir))


################################################################################
### 0. Define colours and functions ###
################################################################################

extra_colours <- c(brewer.pal(12, "Set3"), brewer.pal(12, "Paired"))
col_palette <- c(brewer.pal(8, "Dark2")[c(3:6,8)], brewer.pal(12, "Set3"), brewer.pal(8, "Accent"),
  "#660200", "#918940", "black", "#9ECAE1", "#3F007D", "#67000D", "#FDD0A2", "#08519C",
  "#DBB335", "#5AA050", "#807DBA", "#1A6000", "#F16913", "#FD8D3C", "#DEEBF7", 
  "#7F2704", "#DADAEB", "#FC9272", "#BCBDDC", extra_colours, "#b2182b", "#85929E", 
  "#9B59B6", "#74add1","#1b7837", "#b8e186", "#fed976","#e7298a", "#18ffff", "#ef6c00",
  "#A93226", "#E611ED","orange", "#b8bc53", "#5628ce", "#fa909c", "#8ff331","#270e26")
col_palette <- col_palette[-7]

subcluster_cols <- read.table(
  paste0(ref_dir, "CNV_colour_palette.txt"),
  header = F,
  stringsAsFactors = F,
  comment.char = ""
)[,1]

plot_DE <- dget(paste0(func_dir, "plot_DE.R"))


################################################################################
### 1. Load and format data ###
################################################################################

# load seurat object:
seurat_10X <- readRDS(paste0(in_dir, "04_seurat_object_annotated.Rdata"))

# load metadata and subset:
epi_meta <- readRDS(
  paste0(
    Robject_dir, 
    "5b.final_epithelial_metadata_without_normals.Rdata"
  )
)

subcluster_meta <- subset(epi_meta, select = c("cell_ids", "subcluster_id"))

# change seurat object Idents to subcluster ids:
temp_idents <- as.character(Idents(seurat_10X))
names(temp_idents) <- as.character(names(Idents(seurat_10X)))
m <- match(subcluster_meta$cell_ids, names(temp_idents))
temp_idents[m] <- as.character(subcluster_meta$subcluster_id)
Idents(seurat_10X) <- factor(
  as.character(temp_idents),
  levels = naturalsort(unique(temp_idents))
)

# subset object keeping only subcluster cells:
seurat_sub <- subset(
  seurat_10X,
  idents = as.character(unique(subcluster_meta$subcluster_id))
)


################################################################################
### 2. DE between met and non-met populations ###
################################################################################

# Find upregulated genes in met population (to be used for met gene list):
met_DE_genes <- plot_DE(
  seurat_object = seurat_sub,
  min_pct = 0.1,
  logfc_thresh = 0.4,
  return_thresh = 0.05,
  only.pos = TRUE,
  plot_dir,
  filter_sig = TRUE,
  CNA_assoc_only = TRUE,
  table_dir,
  Robject_dir,
  ref_dir,
  file_prefix = "met_vs_non-met_CNA_assoc_subpop_upreg",
  specific_DE = specific_DE,
  group1 = DE1,
  group2 = DE2,
  specific_genes = "none",
  return_table = TRUE
)

# filter metastasis genes for those in GSEA met gene sets:
if (!file.exists(paste0(ref_dir, "msigdbr_met_genes.Rdata"))) {
	library(msigdbr, lib.loc = lib_loc)
  # fetch all cancer associated GSEA gene set categories:
  m_df = msigdbr(species = "Homo sapiens")
  cancer_cats_temp <- m_df[m_df$gs_cat == "C2" & m_df$gs_subcat == "CGP",]
  cancer_cats <- rbind(
    cancer_cats_temp, 
    m_df[m_df$gs_cat == "C6",]
  )
  
  met_indices <- grep("metasta", cancer_cats$gs_name, ignore.case = T)
  remove_indices <- grep(
    "strom|lymphatic_vessel", 
    cancer_cats$gs_name, 
    ignore.case = TRUE
  )
  met_indices <- met_indices[!(met_indices %in% remove_indices)]
  msig_met_genes <- data.frame(
    gene = cancer_cats$human_gene_symbol[met_indices],
    set = cancer_cats$gs_name[met_indices],
    cat = cancer_cats$gs_cat[met_indices],
    subcat = cancer_cats$gs_subcat[met_indices]
  )
  saveRDS(msig_met_genes, paste0(ref_dir, "msigdbr_met_genes.Rdata"))

  # add non-cancer negative controls:
  ctl_indices <- which(
    m_df$gs_name %in% 
      c(
        "GO_REGULATION_OF_HAIR_FOLLICLE_DEVELOPMENT", 
        "KEGG_TASTE_TRANSDUCTION", 
        "REACTOME_HIV_TRANSCRIPTION_INITIATION"
      )
    )
  msig_ctl_genes <- data.frame(
    gene = m_df$human_gene_symbol[ctl_indices],
    set = m_df$gs_name[ctl_indices],
    cat = m_df$gs_cat[ctl_indices],
    subcat = m_df$gs_subcat[ctl_indices]
  )

  msig_ctl_and_met <- list(
    ctls = msig_ctl_genes,
    met = msig_met_genes
  )
  saveRDS(msig_ctl_and_met, paste0(ref_dir, "msigdbr_met_and_ctl_genes.Rdata"))

} else {
  msig_met_genes <- readRDS(paste0(ref_dir, "msigdbr_met_genes.Rdata"))
}

# check whether each upregulated met gene falls under a msigdb term:
gene_list <- as.list(as.character(met_DE_genes$gene))
met_assoc <- lapply(gene_list, function(y) {
  return(msig_met_genes[msig_met_genes$gene == y,])
})
res <- do.call("rbind", met_assoc[lapply(met_assoc, nrow) > 0])
msig_DE_genes <- res[order(res$gene),]

# take only genes reported as upregulated
up_indices <- grep("up", msig_DE_genes$set, ignore.case = TRUE)
down_indices <- grep("down|dn", msig_DE_genes$set, ignore.case = TRUE)

msig_met_DE_genes <- msig_DE_genes[up_indices,]

write.table(
  msig_met_DE_genes,
  paste0(ref_dir, "msigdbr_met_genes.txt"),
  row.names = F,
  col.names = T,
  quote = F
)

# only take manual genes included in DE list:
manual_met_genes <- manual_met_genes[
	manual_met_genes %in% met_DE_genes$gene
]
# select and order met genes to plot:
plot_met_genes <- met_DE_genes[
	met_DE_genes$gene %in% unique(
    c(manual_met_genes, as.character(msig_met_DE_genes$gene))
  ),
]
plot_met_genes <- arrange(plot_met_genes, desc(avg_logFC))

write.table(
  plot_met_genes,
  paste0(ref_dir, "CNA_met_assoc_DE_genes.txt"),
  row.names = F,
  col.names = T,
  quote = F
)

# plot selected DE genes for matched met sample:
plot_DE(
  seurat_object = seurat_sub,
  min_pct = 0.1,
  logfc_thresh = 0.4,
  return_thresh = 0.05,
  only.pos = TRUE,
  plot_dir,
  filter_sig = TRUE,
  CNA_assoc_only = TRUE,
  table_dir,
  Robject_dir,
  ref_dir,
  file_prefix = "met_and_CNA_assoc_subpop_DE",
  specific_DE = specific_DE,
  group1 = DE1,
  group2 = DE2,
  specific_genes = plot_met_genes$gene
)

# loose DE for use with iDEA:
plot_DE(
  seurat_object = seurat_sub,
  min_pct = 0,
  logfc_thresh = 0,
  return_thresh = 1,
  only.pos = FALSE,
  plot_dir,
  filter_sig = FALSE,
  CNA_assoc_only = FALSE,
  table_dir,
  Robject_dir,
  ref_dir,
  file_prefix = "all_gene_met_vs_non-met_subpop_DE",
  specific_DE = specific_DE,
  group1 = DE1,
  group2 = DE2
)

