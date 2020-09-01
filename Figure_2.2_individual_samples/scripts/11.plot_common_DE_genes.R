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
remove_artefacts <- args[5]
subset_samples <- as.logical(args[6])

project_name <- "thesis"
subproject_name <- "Figure_2.2_individual_samples"
args = commandArgs(trailingOnly=TRUE)
subcluster_method <- "random_trees"
subcluster_p <- "0.05"
if (subcluster_p != "none") {
  subcluster_p <- as.numeric(subcluster_p)
}
coverage_filter <- "filtered"
remove_artefacts <- "artefacts_not_removed"
subset_samples <- TRUE
min_common_prop <- 0.2

if (subset_samples) {

  # determine order of samples by subtype:
  ER <- c("CID4463")
  HER2 <- c("CID45171")
  TNBC <- c("CID4515")
  sample_names <- c(ER, HER2, TNBC)

  subtype_df <- data.frame(
    subtype = c(
      rep("ER", length(ER)),
      rep("HER2", length(HER2)),
      rep("TNBC", length(TNBC))
    ),
    sample = sample_names
  )

} else {

  # determine order of samples by subtype:
#  ER <- c("CID3941", "CID3948", "CID4067", "CID4290A", "CID4398", "CID4461", 
#  	"CID4463", "CID4530N", "CID4535")
  ER <- c("CID3941", "CID4067", "CID4290A", "CID4461", 
    "CID4463", "CID4530N", "CID4535")
  HER2 <- c("CID3921", "CID3586", "CID3963", "CID4066", "CID45171")
  TNBC <- c("CID44041", "CID4465", "CID44971", "CID44991", "CID4513",
  	"CID4515", "CID4523")
  sample_names <- c(ER, HER2, TNBC)

  subtype_df <- data.frame(
    subtype = c(
      rep("ER", length(ER)),
      rep("HER2", length(HER2)),
      rep("TNBC", length(TNBC))
    ),
    sample = sample_names
  )

}

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(ComplexHeatmap, lib.loc = lib_loc)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/",
  subproject_name, "/")
ref_dir <- paste0(project_dir, "/refs/")
func_dir <- paste0(project_dir, "/scripts/functions/")
results_dir <- paste0(project_dir, "results/")
in_path <- paste0(results_dir, "infercnv/")
out_dir <- paste0(results_dir, "infercnv/common_DE/", 
  coverage_filter, "/", subcluster_method, "/p_", subcluster_p, "/", 
  remove_artefacts, "/")

if (subset_samples) {
  Robject_dir <- paste0(out_dir, "Rdata_sub/")
  plot_dir <- paste0(out_dir, "plots_sub/")
  table_dir <- paste0(out_dir, "tables_sub/")
} else {
  Robject_dir <- paste0(out_dir, "Rdata/")
  plot_dir <- paste0(out_dir, "plots/")
  table_dir <- paste0(out_dir, "tables/")
}

system(paste0("mkdir -p ", Robject_dir))
system(paste0("mkdir -p ", plot_dir))
system(paste0("mkdir -p ", table_dir))

print(paste0("In path = ", in_path))
print(paste0("Out dir = ", out_dir))
print(paste0("Reference directory = ", ref_dir))
print(paste0("R function directory = ", func_dir))
print(paste0("R object directory = ", Robject_dir))
print(paste0("Plot directory = ", plot_dir))


################################################################################
### 0. Load functions and colour palette ###
################################################################################

subtype_cols <- read.table(
  paste0(ref_dir, "subtype_colour_palette.txt"),
  header = F,
  comment.char = "",
  stringsAsFactors = F
)
m <- match(c("ER", "HER2", "TNBC"), subtype_cols$V1)
subtype_cols <- subtype_cols$V2[m]
names(subtype_cols) <- c("ER", "HER2", "TNBC")


################################################################################
### 1. Load DE data ###
################################################################################

if (!file.exists(paste0(Robject_dir, "group_DE_data.Rdata"))) {

  for (s in 1:length(sample_names)) {

    print(paste0("Loading ", sample_names[s], " CNV data..."))

    sample_table_dir <- paste0(
  	  in_path, sample_names[s], "/", coverage_filter, "/", subcluster_method, 
      "/p_", subcluster_p, "/", remove_artefacts, 
      "/tables/"
    )

    DE_data <- read.table(
  		paste0(sample_table_dir, "subpop_DE.txt"),
      header = T
	  )

    all_gene_DE_data <- read.table(
      paste0(sample_table_dir, "all_gene_subpop_DE.txt"),
      header = T
    )

    # only keep genes with adjusted p value < 0.1:
    filtered_DE <- DE_data[DE_data$p_val_adj < 0.1,]
  
    if (s==1) {

      all_DE_data <- list(filtered_DE)
      names(all_DE_data) <- sample_names[s]

      all_gene_data <- list(all_gene_DE_data)
      names(all_gene_data) <- sample_names[s]

    } else {

      all_DE_data[[s]] <- filtered_DE
      names(all_DE_data)[s] <- sample_names[s]

      all_gene_data[[s]] <- all_gene_DE_data
      names(all_gene_data)[s] <- sample_names[s]

    }

  }

  saveRDS(all_DE_data, paste0(Robject_dir, "group_DE_data.Rdata"))
  saveRDS(all_gene_data, paste0(Robject_dir, "group_all_gene_DE_data.Rdata"))

} else {

  all_DE_data <- readRDS(paste0(Robject_dir, "group_DE_data.Rdata"))
  all_gene_data <- readRDS(paste0(Robject_dir, "group_all_gene_DE_data.Rdata"))
  
}


################################################################################
### 2. Identify DE genes present in a proportion of samples ###
################################################################################

# find genes present in all elements of all_DE_data list:
common_DE <- unique(
  unlist(
    lapply(all_DE_data, function(x) {
      gene_vec <- as.character(
        unlist(
          lapply(all_DE_data, function(y) {
            return(
              unique(
                x$gene[
                  which(x$gene %in% y$gene)
                ]
              )
            )
          })
        )
      )
      gene_tab <- table(gene_vec)
      return(
        names(gene_tab)[
          which(gene_tab >= ceiling(length(all_DE_data)*min_common_prop))
        ]
      )
    })
  )
)


################################################################################
### 3. Prepare DE genes present in a proportion of samples ###
################################################################################

if (length(common_DE) > 0) {

  # filter out non-common genes in all_gene_data:
  common_DE_data <- lapply(all_gene_data, function(x) {

    common_df <- x[x$gene %in% common_DE,]
    common_df$sig <- FALSE
    common_df$sig[common_df$p_val_adj < 0.1] <- TRUE

    # split by gene:
    split_df <- split(common_df, as.character(common_df$gene))

    # deal with duplicate genes:
    res_df <- do.call(
      "rbind",lapply(split_df, function(y) {

        if (nrow(y) == 1) {
  
          return(
            subset(
              y,
              select = c(gene, p_val_adj, avg_logFC, sig)
            )
          )
  
        } else {
  
          if (all(y$sig == TRUE) | all(y$sig == FALSE)) {
            return(
              subset(
                aggregate(
                  .~gene,
                  common_df, 
                  mean
                ),
                select = c(gene, p_val_adj, avg_logFC, sig)
              )
            )
  
          } else {

            return(
              subset(
                y[y$sig,],
                select = c(gene, p_val_adj, avg_logFC, sig)
              )
            )
  
          }
  
        }
      })
    )

    return(res_df)

  })

}


################################################################################
### 4. Plot common DE genes ###
################################################################################

final_heatmap <- Heatmap(
  as.matrix(plot_object), name = paste0("hm"), 
  col = heatmap_cols,
  cluster_columns = F, cluster_rows = F,
  show_row_names = F, show_column_names = T,
  column_names_gp = gpar(col = "white"),
  show_row_dend = F,
  show_heatmap_legend = F,
  use_raster = T, raster_device = c("png")
)


