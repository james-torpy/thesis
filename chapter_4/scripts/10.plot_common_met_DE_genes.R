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
subset_samples <- as.logical(args[6])
pos_only <- as.logical(args[7])
external_list <- as.logical(args[8])
external_list_name <- as.logical(args[9])

project_name <- "thesis"
subproject_name <- "chapter_4"
args = commandArgs(trailingOnly=TRUE)
subcluster_method <- "random_trees"
subcluster_p <- "0.05"
if (subcluster_p != "none") {
  subcluster_p <- as.numeric(subcluster_p)
}
coverage_filter <- "filtered"
remove_artefacts <- "artefacts_not_removed"
subset_samples <- FALSE
pos_only <- TRUE
external_list <- "none"
external_list_name <- "none"
#external_list <- c("ANGPT1", "CASP1", "FGFR1", "GM19589", "Il1R1",
#  "MGP", "MFSS1", "NEURL1B", "PHEX", "PLSCR4", "SERPINE2", "SLPI")
#external_list_name <- "Wagenblast"

if (subset_samples) {

  # determine order of samples by subtype:
  ER <- c("CID4463")
  HER2 <- c("CID45171")
  TNBC <- c("CID44991")
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
out_dir <- paste0(results_dir, "infercnv/met_gene_expression/", 
  coverage_filter, "/", subcluster_method, "/p_", subcluster_p, "/", 
  remove_artefacts, "/")

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

filter_common_DE_genes <- dget(paste0(func_dir, "filter_common_DE_genes.R"))

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
### 1. Load DE data and met gene list ###
################################################################################

if (!file.exists(paste0(Robject_dir, "common_subpop_met_DE_data.Rdata"))) {

  for (s in 1:length(sample_names)) {

    print(paste0("Loading ", sample_names[s], " CNV data..."))

    sample_table_dir <- paste0(
  	  in_path, sample_names[s], "/", coverage_filter, "/", subcluster_method, 
      "/p_", subcluster_p, "/", remove_artefacts, 
      "/tables/"
    )


    all_gene_DE_data <- read.table(
      paste0(
      	sample_table_dir, 
      	"all_CNA_assoc_subpop_DE.txt"
      ),
      header = T
    )

    if (s==1) {

      all_gene_data <- list(all_gene_DE_data)
      names(all_gene_data) <- sample_names[s]

    } else {

      all_gene_data[[s]] <- all_gene_DE_data
      names(all_gene_data)[s] <- sample_names[s]

    }

  }

  saveRDS(
    all_gene_data, 
    paste0(Robject_dir, "common_subpop_met_DE_data.Rdata")
  )

} else {

  all_gene_data <- readRDS(
    paste0(Robject_dir, "common_subpop_met_DE_data.Rdata")
  )
  
}

# load metastasis- and CNA- associated genes:
met_genes <- read.table(
  paste0(ref_dir, "CNA_met_assoc_DE_genes.txt"),
  header = T,
  stringsAsFactors = F
)$gene


################################################################################
### 2. Prepare DE genes present in a proportion of samples ###
################################################################################

# merge together met DE data and filter for metastasis genes:
if (pos_only) {

  met_DE_data <- lapply(
    all_gene_data, 
    filter_common_DE_genes, 
    met_genes, 
    "up",
    func_dir
  )

  empty_df <- data.frame(
    gene = sort(
      met_genes
    ),
    stringsAsFactors = FALSE
  )

  plot_data <- list(
    avg_logFC = empty_df,
    p_val_adj = empty_df
  )

} else {

  met_DE_data <- lapply(
    all_gene_data, 
    filter_common_DE_genes, 
    met_genes, 
    "bi",
    func_dir
  )

  empty_df <- data.frame(
    gene = sort(
      c(
        paste0(met_genes, "_down"),
        paste0(met_genes$gene, "_up")
      )
    ),
    stringsAsFactors = FALSE
  )

  plot_data <- list(
    avg_logFC = empty_df,
    p_val_adj = empty_df
  )

}

for (l in 1:length(plot_data)) {

  for (s in 1:length(met_DE_data)) {
    temp_df <- cbind(
      met_DE_data[[s]]$gene, 
      as.data.frame(
        eval(parse(text = paste0("met_DE_data[[s]]$", names(plot_data)[l])))
      )
    )
    colnames(temp_df) <- c("gene", names(plot_data)[l])
    plot_data[[l]] <- merge(
      plot_data[[l]],
      temp_df,
      by = "gene",
      all.x = TRUE
    )
    colnames(plot_data[[l]])[s+1] <- names(met_DE_data)[s]
  }

  # remove rows with all NAs:
  plot_data[[l]] <- plot_data[[l]][
    apply(plot_data[[l]], 1, function(x) {
      !(all(is.na(x)))
    }),
  ]

  plot_data[[l]] <- plot_data[[l]] %>%
    column_to_rownames(loc = 1)

}

# ensure p_vals rows match to logFC:
plot_data$p_val_adj <- plot_data$p_val_adj[rownames(plot_data$avg_logFC),]
print(paste0(
  "Are logFC rownames identical to p_vals rownames? ", 
  identical(rownames(plot_data$avg_logFC), rownames(plot_data$p_val_adj)))
)
# change significant values to the ComplexHeatmap code for dots, and all
# others to NA:
p_val_dots <- plot_data$p_val_adj
p_val_dots[p_val_dots < 0.1] <- "*"
p_val_dots[p_val_dots != "*"] <- " "
p_val_dots[is.na(p_val_dots)] <- " "


################################################################################
### 3. Plot common DE genes ###
################################################################################

# define heatmap colours:
na_less_vector <- unlist(plot_data$avg_logFC)
na_less_vector <- na_less_vector[!is.na(na_less_vector)]
if (pos_only) {
  heatmap_cols <- colorRamp2(c(0, max(na_less_vector)), 
      c("white", "#870C0C"), space = "sRGB")
} else {
  heatmap_cols <- colorRamp2(c(min(na_less_vector), 0, max(na_less_vector)), 
    c("#0D3084", "white", "#870C0C"), space = "sRGB")
}
final_heatmap <- Heatmap(
  as.matrix(plot_data$avg_logFC), 
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

if (pos_only) {
  pdf(paste0(plot_dir, "met_and_CNA_assoc_upreg_gene_heatmap.pdf"), 
    height = 13, width = 20) 
} else {
  pdf(paste0(plot_dir, "met_and_CNA_assoc_DE_gene_heatmap.pdf"), 
    height = 13, width = 20) 
}

  grid.newpage()

    # plot heatmap:
    pushViewport(viewport(x = 0.155, y = 0.04, width = 0.82, height = 0.95, 
      just = c("left", "bottom")))
      grid.draw(annotated_heatmap)
    popViewport()
dev.off()



