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
subset_samples <- FALSE
min_common_prop <- 0.2
min_pct_all_gene <- 0.5

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
library(circlize, lib.loc = lib_loc)
library(dplyr)
library(textshape, lib.loc = lib_loc)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/",
  subproject_name, "/")
ref_dir <- paste0(project_dir, "/refs/")
func_dir <- paste0(project_dir, "/scripts/functions/")
results_dir <- paste0(project_dir, "results/")
in_path <- paste0(results_dir, "infercnv/")
out_dir <- paste0(results_dir, "infercnv/common_DE/", 
  coverage_filter, "/", subcluster_method, "/p_", subcluster_p, "/", 
  remove_artefacts, "/tables/")

if (subset_samples) {
  Robject_dir <- paste0(out_dir, "Rdata_sub/")
  plot_dir <- paste0(out_dir, "plots_sub/")
  table_dir <- paste0(out_dir, "tables_sub/")
} else {
  Robject_dir <- paste0(out_dir, "Rdata/")
  plot_dir <- paste0(out_dir, "plots/")
  table_dir <- out_dir
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

if (!file.exists(paste0(Robject_dir, "common_subpop_DE_data.Rdata"))) {

  for (s in 1:length(sample_names)) {

    print(paste0("Loading ", sample_names[s], " CNV data..."))

    sample_table_dir <- paste0(
  	  in_path, sample_names[s], "/", coverage_filter, "/", subcluster_method, 
      "/p_", subcluster_p, "/", remove_artefacts, 
      "/tables/"
    )

    DE_data <- read.table(
  		paste0(
        sample_table_dir, 
        "subpop_DE_min_pct_0.5_logfc_0.7_p_val_0.01_CNA_assoc_only.txt"
      ),
      header = T
	  )

    all_gene_DE_data <- read.table(
      paste0(
      	sample_table_dir, 
      	"subpop_DE_min_pct_0.1_logfc_0_p_val_1_CNA_assoc_only.txt"
      ),
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

  saveRDS(
    all_DE_data, 
    paste0(Robject_dir, "common_subpop_DE_data.Rdata")
  )
  saveRDS(
    all_gene_data, 
    paste0(Robject_dir, "common_all_gene_subpop DE_data.Rdata")
  )

} else {

  all_DE_data <- readRDS(
    paste0(Robject_dir, "common_subpop_DE_data.Rdata")
  )
  all_gene_data <- readRDS(
    paste0(Robject_dir, "common_all_gene_subpop DE_data.Rdata")
  )
  
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

# save common DE genes to use for looking up articles:
write.table(
  data.frame(gene = common_DE),
  paste0(table_dir, "sig_subpop_DE_CNA_assoc_only.txt"),
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)


################################################################################
### 3. Prepare DE genes present in a proportion of samples ###
################################################################################

if (length(common_DE) > 0) {

  if (subset_samples) {
    common_DE <- common_DE[1:50]
  }

  # filter out non-common genes in all_gene_data:
  common_DE_data <- lapply(all_gene_data, function(x) {

    common_df <- x[x$gene %in% common_DE,]
    common_df$gene <- as.character(common_df$gene)
    common_df$sig <- FALSE
    common_df$sig[common_df$p_val_adj < 0.1] <- TRUE

    # split by gene:
    split_by_gene <- split(common_df, as.character(common_df$gene))

    # keep only values with largest distance from 0:
    res_df <- do.call(
      "rbind",lapply(split_by_gene, function(y) {

        if (nrow(y) == 1) {
  
          return(
            subset(
              y,
              select = c(gene, p_val_adj, avg_logFC, sig)
            )
          )
  
        } else {
  
          if (all(y$sig == TRUE) | all(y$sig == FALSE)) {

            # determine difference of each logFC from 0:
            y$diff <- abs(0-y$avg_logFC)

            return(
              subset(
                merge(aggregate(diff ~ gene, max, data = y), y),
                select = c(gene, p_val_adj, avg_logFC, sig)
              )
            )
  
          } else {

            y <- y[y$sig == TRUE,]

            return(
              subset(
                merge(aggregate(diff ~ gene, max, data = y), y),
                select = c(gene, p_val_adj, avg_logFC, sig)
              )
            )
  
          }
  
        }
      })
    )
    
    return(combined_df)

  })

  # merge together common DE data:
  logFC <- data.frame(
    gene = common_DE,
    stringsAsFactors = FALSE
  )

  for (s in 1:length(common_DE_data)) {

    logFC <- merge(
      logFC,
      subset(common_DE_data[[s]], select = c(gene, avg_logFC)),
      by = "gene",
      all.x = TRUE
    )
    colnames(logFC)[s+1] <- names(common_DE_data)[s]

  }

  logFC <- logFC %>%
    column_to_rownames(loc = 1)

  p_vals <- data.frame(
    gene = common_DE
  )

  for (s in 1:length(common_DE_data)) {

    p_vals <- merge(
      p_vals,
      subset(common_DE_data[[s]], select = c(gene, p_val_adj)),
      by = "gene",
      all.x = TRUE
    )
    colnames(p_vals)[s+1] <- names(common_DE_data)[s]

  }

  p_vals <- p_vals %>%
    column_to_rownames(loc = 1)

  # change significant values to the ComplexHeatmap code for dots, and all
  # others to NA:
  p_val_dots <- p_vals
  p_val_dots[p_val_dots < 0.1] <- "*"
  p_val_dots[p_val_dots != "*"] <- " "
  p_val_dots[is.na(p_val_dots)] <- " "


  ################################################################################
  ### 4. Plot common DE genes ###
  ################################################################################
  
  # define heatmap colours:
  na_less_vector <- unlist(logFC)
  na_less_vector <- na_less_vector[!is.na(na_less_vector)]
  heatmap_cols <- colorRamp2(c(min(na_less_vector), 0, max(na_less_vector)), 
    c("#0D3084", "white", "#870C0C"), space = "sRGB")

  final_heatmap <- Heatmap(
    as.matrix(logFC), 
    na_col = "grey",
    name = paste0("hm"), 
    col = heatmap_cols,
    cluster_columns = F, cluster_rows = F,
    show_row_names = T, show_column_names = T,
    show_row_dend = T,
    show_heatmap_legend = T,
    use_raster = T, raster_device = c("png"),
    cell_fun = function(j, i, x, y, w, h, heatmap_cols) { # add text to each grid
        grid.text(p_val_dots[i, j], x, y, gp=gpar(fontsize=30, col = "#1F0B99"))
    }
  )

  ht_list <- final_heatmap

  annotated_heatmap <- grid.grabExpr(
    draw(ht_list
      #, gap = unit(6, "mm"), heatmap_legend_side = "left"
    )
  )
  dev.off()

  pdf(paste0(plot_dir, "common_gene_DE_CNA_assoc_only.pdf"), 
    height = 13, width = 20) 
  
    grid.newpage()
  
      # plot heatmap:
      pushViewport(viewport(x = 0.155, y = 0.04, width = 0.82, height = 0.95, 
        just = c("left", "bottom")))
        #grid.rect()
        grid.draw(annotated_heatmap)
      popViewport()

  dev.off()

}




