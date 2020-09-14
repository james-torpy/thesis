#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

project_name <- "thesis"
subproject_name <- "Figure_2.2_individual_samples"
sample_name <- "CID4463"
subcluster_method <- "random_trees"
subcluster_p <- "0.05"
coverage_filter <- "filtered"
remove_artefacts <- "artefacts_not_removed"
CNA_assoc <- TRUE
key_terms <- c("cancer", "tumour", "carcinoma", "malignant", 
  "metastasis", "resistance", "oncogene")

sample_name <- args[1]
subcluster_method <- args[2]
subcluster_p <- args[3]
if (subcluster_p != "none") {
  subcluster_p <- as.numeric(subcluster_p)
}
coverage_filter <- args[4]
remove_artefacts <- args[5]
CNA_assoc <- as.logical(args[6])
key_terms <- args[7]  
key_terms <- strsplit(
  key_terms,
  "_"
)[[1]]

lib_loc = "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(vctrs, lib.loc = lib_loc)
library(rlang, lib.loc = lib_loc)
library(msigdbr, lib.loc = lib_loc)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/", 
	subproject_name, "/")
results_dir <- paste0(project_dir, "results/")
func_dir <- paste0(project_dir, "scripts/functions/")
in_path <- paste0(results_dir, "infercnv/")
in_dir <- paste0(in_path, sample_name, "/", coverage_filter, "/",
  subcluster_method, "/p_", subcluster_p, "/", remove_artefacts, "/")
out_dir <- paste0(in_dir, "/tables/")

ref_dir <- paste0(project_dir, "refs/")

if (sample_name == "common_DE") {
  if (file.exists(paste0(in_dir, "sig_subpop_DE_CNA_assoc_only.txt"))) {
    
    print(paste0("Loading DE data for ", sample_name))
    
    DE_genes <- read.table(
      paste0(
        in_dir, 
        "tables/subpop_DE_min_pct_0.1_logfc_0_p_val_1_CNA_assoc_only.txt"
      ),
      header = TRUE
    )
  } else {
    print("DE data does not exist for ", sample_name, "...")
  }
        
} else {
  if (
    file.exists(
      paste0(
        in_dir, 
        "tables/subpop_DE_min_pct_0.5_logfc_0.7_p_val_0.01_CNA_assoc_only.txt"
      )
    )
  ) {
    print(paste0("Loading DE data for ", sample_name))
    
    DE_genes <- read.table(
      paste0(
        in_dir, 
        "tables/subpop_DE_min_pct_0.5_logfc_0.7_p_val_0.01_CNA_assoc_only.txt"
      ),
      header = TRUE
    )
  } else {
    print(paste0(
      "DE data does not exist for ", sample_names[i], "..."
    ))
  }
}

if (!file.exists(paste0(ref_dir, "gsea_brca_genes.Rdata"))) {
 
  # fetch all cancer associated GSEA gene set categories:
  m_df = msigdbr(species = "Homo sapiens")
  cancer_cats_temp <- m_df[m_df$gs_cat == "C2" & m_df$gs_subcat == "CGP",]
  cancer_cats <- rbind(
    cancer_cats_temp, 
    m_df[m_df$gs_cat == "C6",]
  )
  
  # fetch all GSEA sets with key terms in category name:
  if (exists("cancer_genes")) {
    rm(cancer_genes)
  }
  for (t in 1:length(key_terms)) {
  
    if (t==1) {
  
      cancer_genes <- data.frame(
        gene = cancer_cats$human_gene_symbol[
          grep(key_terms[t], cancer_cats$gs_name, ignore.case = T)
        ],
        set = cancer_cats$gs_name[
          grep(key_terms[t], cancer_cats$gs_name, ignore.case = T)
        ],
        cat = cancer_cats$gs_cat[
          grep(key_terms[t], cancer_cats$gs_name, ignore.case = T)
        ],
        subcat = cancer_cats$gs_subcat[
          grep(key_terms[t], cancer_cats$gs_name, ignore.case = T)
        ]
      )
  
    } else {
  
      cancer_genes <- rbind(
        cancer_genes,
        data.frame(
          gene = cancer_cats$human_gene_symbol[
            grep(key_terms[t], cancer_cats$gs_name, ignore.case = T)
          ],
          set = cancer_cats$gs_name[
            grep(key_terms[t], cancer_cats$gs_name, ignore.case = T)
          ],
          cat = cancer_cats$gs_cat[
            grep(key_terms[t], cancer_cats$gs_name, ignore.case = T)
          ],
          subcat = cancer_cats$gs_subcat[
            grep(key_terms[t], cancer_cats$gs_name, ignore.case = T)
          ]
        )
      )
  
    }
    
  }
  
  # keep only brca- or metastasis- associated genes:
  brca_genes <- cancer_genes[
    grep(
      "breast|metastasis|resistance", 
      cancer_genes$set, 
      ignore.case = T
    ),
  ]
  brca_genes <- brca_genes[
    grep(
      "PREDNISOLONE|INSULIN|CLUSTER|VINCRISTINE|SMID|FARMER|CHARAFE|BERTUCCI", 
      brca_genes$set, 
      ignore.case = T, 
      invert = T
    ),
  ]

  saveRDS(brca_genes, paste0(ref_dir, "gsea_brca_genes.Rdata"))

} else {
  brca_genes <- readRDS(paste0(ref_dir, "gsea_brca_genes.Rdata"))
}

# check whether each DE gene falls under a GSEA term:
gene_list <- as.list(as.character(DE_genes$gene))
brca_assoc <- lapply(gene_list, function(y) {
  return(brca_genes[brca_genes$gene == y,])
})
res <- do.call("rbind", brca_assoc[lapply(brca_assoc, nrow) > 0])
cancer_DE_genes <- res[order(res$gene),]

write.table(
  cancer_DE_genes,
  paste0(out_dir, "msigdb_cancer_subpop_DE_genes.txt"),
  row.names = F,
  col.names = T,
  quote = F
)

