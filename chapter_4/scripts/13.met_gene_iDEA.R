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
up_annot_only <- as.logical(args[7])
subset_annots <- as.logical(args[8])

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
#up_annot_only <- TRUE
#subset_annots <- TRUE

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("Subclustered by ", subcluster_method))
print(paste0("Subclustered by pval ", subcluster_p))
print(paste0("Filter by coverage? ", coverage_filter))
print(paste0("Remove artefacts? ", remove_artefacts))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(RColorBrewer)
library(naturalsort, lib.loc = lib_loc)
library(dplyr)
library(tibble)
library(iDEA, lib.loc=lib_loc)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/",
  subproject_name, "/")
ref_dir <- paste0(project_dir, "/refs/")
func_dir <- paste0(project_dir, "/scripts/functions/")
results_dir <- paste0(project_dir, "results/")
in_path <- paste0(results_dir, "/infercnv/", sample_name, "/", 
  coverage_filter, "/", subcluster_method, "/p_", subcluster_p, "/", 
  remove_artefacts, "/")

in_table_dir <- paste0(in_path, "tables/")

Robject_dir <- paste0(in_path, "Rdata/iDEA/")
system(paste0("mkdir -p ", Robject_dir))
plot_dir <- paste0(in_path, "plots/iDEA/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(in_path, "tables/iDEA/")
system(paste0("mkdir -p ", table_dir))

print(paste0("In dir = ", in_table_dir))
print(paste0("Out path = ", in_path))
print(paste0("Reference directory = ", ref_dir))
print(paste0("R function directory = ", func_dir))
print(paste0("R object directory = ", Robject_dir))
print(paste0("Table directory = ", table_dir))
print(paste0("Plot directory = ", plot_dir))


################################################################################
### 0. Define colours and functions ###
################################################################################



################################################################################
### 1. Load and format data ###
################################################################################

# load msigdbr metastasis gene sets:
if (!file.exists(paste0(ref_dir, "msigdbr_met_and_ctl_genes.Rdata"))) {
  library(msigdbr, lib.loc = lib_loc)
  # fetch all cancer associated msigdbr gene set categories:
  m_df = msigdbr(species = "Homo sapiens")
  cancer_cats_temp <- m_df[m_df$gs_cat == "C2" & m_df$gs_subcat == "CGP",]
  cancer_cats <- rbind(
    cancer_cats_temp, 
    m_df[m_df$gs_cat == "C6",]
  )
  
  met_indices <- grep("metasta", cancer_cats$gs_name, ignore.case = T)
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

  msig_ctl_and_met <- rbind(
    msig_ctl_genes,
    msig_met_genes
  )
  saveRDS(msig_met_genes, paste0(ref_dir, "msigdbr_met_and_ctl_genes.Rdata"))

} else {
  msig_ctl_and_met  <- readRDS(paste0(ref_dir, "msigdbr_met_and_ctl_genes.Rdata"))
}

# keep only upregulated gene sets if necessary:
if (up_annot_only) {
  keep_indices <- grep(
    "down|dn|repress", 
    msig_ctl_and_met$set, 
    ignore.case = T,
    invert = T
  )
  msig_ctl_and_met <- msig_ctl_and_met[keep_indices,]
}

# split by gene set:
msig_ctl_and_met$set <- as.character(msig_ctl_and_met$set)
msig_list <- split(msig_ctl_and_met, msig_ctl_and_met$set)

# load custom metastasis gene set and format:
subpop_met_genes <- read.table(
  paste0(ref_dir, "msigdbr_met_DE_CNA_genes.txt"),
  header = TRUE,
  stringsAsFactors = FALSE
)

# uppend custom met gene set:
met_list <- append(list(subpop_met_genes), msig_list)
names(met_list)[1] <- "SUBPOP_MET"

# load previous DE and variance results and format:
DE <- read.table(
  paste0(in_table_dir, "all_gene_met_vs_non-met_subpop_DE.txt"),
  header = TRUE,
  stringsAsFactors = FALSE
)

DE <- DE %>% 
  column_to_rownames("gene") %>%
  subset(select = c(avg_logFC, variance))

if (!file.exists(paste0(Robject_dir, "met_annotation_df.Rdata"))) {
  # create annotation df:
  annot_df <- data.frame(
    gene = unique(rownames(DE))
  )
  
  for (l in 1:length(met_list)) {
  
    met_list[[l]]$member = 1
    to_add <- subset(met_list[[l]], select = c(gene, member))
    to_add <- to_add[!(duplicated(to_add$gene)),]
  
    # add extra column indicating genes belonging to set:
    annot_df <- merge(annot_df, to_add, by = "gene", all.x = T)
    annot_df$member[is.na(annot_df$member)] <- 0
    colnames(annot_df)[l+1] <- names(met_list)[l]
  
  }
  
  annot_df <- annot_df %>%
    column_to_rownames("gene")

  saveRDS(annot_df, paste0(Robject_dir, "met_annotation_df.Rdata"))

} else {
  annot_df <- readRDS(paste0(Robject_dir, "met_annotation_df.Rdata"))
}

if (subset_annots) {
  annot_df <- cbind(
  	annot_df[ ,
  	  colnames(annot_df) %in% c(
        "GO_REGULATION_OF_HAIR_FOLLICLE_DEVELOPMENT", 
        "KEGG_TASTE_TRANSDUCTION", 
        "REACTOME_HIV_TRANSCRIPTION_INITIATION"
      )
  	],
  	annot_df[,1:7]
  )
}


################################################################################
### 2.  ###
################################################################################

### to troubleshoot ###

#Following every 'CreateiDEAObject':
#Warning messages:
#1: In selectChildren(pids[!fin], -1) :
#  cannot wait for child 56433 as it does not exist
#2: In parallel::mccollect(...) : 1 parallel job did not deliver a result

###

idea <- CreateiDEAObject(
  DE, 
  annot_df, 
  max_var_beta = 100, 
  min_precent_annot = 0.0025, 
  num_core=1
)

# fit model:
idea <- iDEA.fit(idea,
  fit_noGS=FALSE,
  init_beta=NULL, 
  init_tau=c(-2,0.5),
  min_degene=5,
  em_iter=15,
  mcmc_iter=1000, 
  fit.tol=1e-5,
  modelVariant = F,
  verbose=TRUE
)

# compute GSEA and DE posterior probabilities:
idea <- iDEA.louis(idea)



