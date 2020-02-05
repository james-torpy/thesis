#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

#args = commandArgs(trailingOnly=TRUE)

project_name <- "thesis"
subproject_name <- "Figure_3.5_matched_samples"
sample_name <- "CID45171"
include_t_cells <- TRUE
subpopulations <- list(
  subpop_1 <- 5,
  subpop_2 <- c(seq(0,4), 6, 7)
)

print(paste0("Project name = ", project_name))
print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name is = ", sample_name))
print(paste0("Include T cells? ", as.character(include_t_cells)))
for (i in 1:length(subpopulations)) {
  subpop_clusters <- paste(subpopulations[[i]], collapse=", ")
  print(paste0("Subpopulation ", i, 
    " is composed of the following clusters: ", subpop_clusters))
}


lib_loc <- "/share/ScratchGeneral/jamtor/R/3.5dev/"
library(Seurat)
library(ggplot2)
library(ggrepel)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", 
  project_name, "/", subproject_name, "/")
ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/functions/")
results_dir <- seurat_path <- paste0(project_dir, "results/")

paella_dir <- "/paella/TumourProgressionGroupTemp/projects/brca_mini_atlas/"
reclustered_dir <- paste0(
  paella_dir, 
  "analysis/celltype_reclustering_per_sample.DLR_working/garnett_call_ext_major/Epithelial/seurat/individual/",
  sample_name, "/RObj/"
)

if (include_t_cells) {
  out_path <- paste0(results_dir, "DE/t_cells_included/")
  out_dir <- paste0(out_path, sample_name, "/")
} else {
  out_path <- paste0(results_dir, "DE/t_cells_excluded/")
  out_dir <- paste0(out_path, sample_name, "/")
}

Robject_dir <- paste0(out_dir, "Robjects/")
system(paste0("mkdir -p ", Robject_dir))
plot_dir <- paste0(out_dir, "plots/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(out_dir, "tables/")
system(paste0("mkdir -p ", table_dir))

print(paste0("Sample directory = ", in_dir))
print(paste0("Reference directory = ", ref_dir))
print(paste0("R function directory = ", func_dir))
print(paste0("Output directory = ", out_dir))
print(paste0("Plot directory = ", plot_dir))
print(paste0("Table directory = ", table_dir))

print(paste0("Running DE between ", length(subpopulations), " subpopulations"))


################################################################################
### 1. Load seurat object and fetch cell IDs of each subpopulation ###
################################################################################

# load epithelial reclustered object:
reclustered_10X <- readRDS(paste0(reclustered_dir, "seurat.rds"))
# load correct annotation:
Idents(reclustered_10X) <- reclustered_10X@meta.data$SUBSET_D_res.0.8

# fetch all cell ids of each subpopulation:
subpopulation_cell_ids <- lapply(subpopulations, function(x) {
  return(
    names(Idents(reclustered_10X))[Idents(reclustered_10X) %in% x]
  )
})


################################################################################
### 2. Perform DE between subpopulations ###
################################################################################

# reassign Idents to cells:
temp_idents <- as.character(Idents(reclustered_10X))
names(temp_idents) <- as.character(names(Idents(reclustered_10X)))

temp_idents[
  names(temp_idents) %in% subpopulation_cell_ids[[1]]
] <- "subpopulation_1"
temp_idents[
  names(temp_idents) %in% subpopulation_cell_ids[[2]]
] <- "subpopulation_2"

Idents(reclustered_10X) <- factor(temp_idents)
names(Idents(reclustered_10X)) <- names(temp_idents)

all_DE <- FindAllMarkers(
  only.pos = T,
  object = reclustered_10X,
  min.pct = 0.5, 
  logfc.threshold = 0.7, 
  test.use = 'MAST'
)

all_DE_sorted <- arrange(all_DE, cluster, desc(avg_logFC))

write.table(all_DE_sorted, paste0(table_dir, "DE_subpop_1_vs_2.txt"))

heatmap_genes <- all_DE_sorted %>% 
  group_by(cluster) %>% 
  top_n(10, avg_logFC)

pdf(paste0(plot_dir, "top_DE_subpop_1_vs_2_heatmap.pdf"))
  print(DoHeatmap(
    reclustered_10X,
    features = heatmap_genes$gene,
    group.by = "ident"
  ))
dev.off()



if ( !file.exists(paste0(Robject_dir, "DE_subpop_1_vs_2.RData")) | 
  !file.exists(paste0(plot_dir, "DE_subpop_1_vs_2_volcano_plot.pdf")) ) {
  DE1 <- FindMarkers(
    reclustered_10X, 
    ident.1 = "subpopulation_1",
    ident.2 = "subpopulation_2",
    logfc.threshold = 0,
  )

  DE1$genes <- rownames(DE1)

  # change logFC to log2 FC:
  DE1$avg_logFC <- log2(
    exp(DE1$avg_logFC)
  )
  colnames(DE1) <- gsub("logFC", "logFC2", colnames(DE1))
  # annotate significantly DE genes:
  DE1$significant <- "non_significant"
  DE1$significant[DE1$p_val_adj < 0.05] <- "significant"

  # check DE genes:
  # log(0.5) = -0.7, log(1) = 0, log(2) = 0.7
  print(paste0(
    "No. genes reported as significantly DE in DE1: ", 
    length(which(DE1$significant == "significant"))
  ))
  upreg_DE1 <- DE1[DE1$avg_logFC2 >=1 & DE1$significant == "significant",]
  print(paste0("No. genes upregulated at least 2x in DE1: ", nrow(upreg_DE1)))
  downreg_DE1 <- DE1[DE1$avg_logFC2 <= -1 & DE1$significant == "significant",]
  print(paste0("No. genes downregulated at least 0.5x in DE1: ", nrow(downreg_DE1)))

  DE1_top_up <- top_n(DE1[DE1$significant != "non_significant",], 10, avg_logFC2)
  DE1_top_up <- DE1_top_up[order(DE1_top_up$avg_logFC2, decreasing = T),]
  DE1$significant[DE1$genes %in% DE1_top_up$genes] <- "top_upregulated"
  DE1_top_down <- top_n(DE1[DE1$significant != "non_significant",], -10, avg_logFC2)
  DE1_top_down <- DE1_top_down[order(DE1_top_down$avg_logFC2, decreasing = T),]
  DE1$significant[DE1$genes %in% DE1_top_down$genes] <- "top_downregulated"
  DE1$significant <- factor(
    DE1$significant, 
    levels = c("non-significant", "significant", "top_upregulated", "top_downregulated")
  )
  
  # volcano plot:
  p <- ggplot(data=DE1, aes( x=avg_logFC2, y=-log10(p_val_adj), color=significant))
  p <- p + geom_point(data=DE1)
  p <- p + geom_text_repel(
    data=DE1[DE1$significant == "top_upregulated" | DE1$significant == "top_downregulated",], 
    aes(label=genes)
  )
  p <- p + theme(legend.position = "none")
  p <- p + labs(x="log fold change", y="-log10 adjusted p.val")
  p <- p + scale_colour_manual(values = c("grey", "forestgreen", "darkorchid1", "darkorchid1"))
  pdf(paste0(plot_dir, "DE_subpop_1_vs_2_volcano_plot.pdf"))
    p
  dev.off()

  # heatmap of top 10 up and downregulated genes:
  top_DE <- rbind(
    DE1_top_up,
    DE1_top_down
  )
  pdf(paste0(plot_dir, "top_DE_subpop_1_vs_2_heatmap.pdf"))
    print(DoHeatmap(
      reclustered_10X,
      features = top_DE$genes,
      group.by = "ident"
    ))
  dev.off()
  
  saveRDS(DE1, paste0(Robject_dir, "DE_subpop_1_vs_2.RData"))

} else {
  DE1 <- readRDS(paste0(Robject_dir, "DE_subpop_1_vs_2.RData"))
}


if (!file.exists(paste0(Robject_dir, "DE_MAST_subpop_1_vs_2.RData"))) {
 
  DE1_MAST <- FindMarkers(
    reclustered_10X, 
    ident.1 = "subpopulation_1",
    ident.2 = "subpopulation_2",
    test.use = "MAST"
  )

  DE1_MAST$genes <- rownames(DE1_MAST)

  sig_DE1_MAST <- DE1_MAST[DE1_MAST$p_val_adj < 0.05,]
  print(paste0("No. genes reported as significantly DE in DE1_MAST: ", nrow(sig_DE1_MAST)))
  upreg_DE1_MAST <- sig_DE1_MAST[sig_DE1_MAST$avg_logFC >= 1,]
  print(paste0("No. genes reported as upregulated in DE1_MAST: ", nrow(upreg_DE1_MAST)))
  downreg_DE1_MAST <- sig_DE1_MAST[sig_DE1_MAST$avg_logFC <= 1,]
  print(paste0("No. genes reported as downregulated in DE1_MAST: ", nrow(downreg_DE1_MAST)))
  
  DE1_MAST_top_up <- top_n(sig_DE1_MAST, 10, avg_logFC)
  DE1_MAST_top_down <- top_n(sig_DE1_MAST, -10, avg_logFC)
  DE1_MAST_lab <- rbind(DE1_MAST_top_up, DE1_MAST_top_down)
  DE1_MAST$top <- DE1_MAST$genes %in% DE1_MAST_lab$genes
  DE1_MAST_lab$top <- TRUE
  
  p <- ggplot(data=DE1_MAST, aes( x=avg_logFC, y=-log10(p_val_adj), color=top))
  p <- p + geom_point(data=DE1_MAST)
  p <- p + geom_text_repel(data=DE1_MAST_lab, aes(label=genes))
  p <- p + theme(legend.position = "none")
  p <- p + labs(x="log fold change", y="-log10 adjusted p.val")
  p <- p + scale_colour_manual(values = c("forestgreen", "darkorchid1"))
  pdf(paste0(plot_dir, "DE_MAST_subpop_1_vs_2_volcano_plot.pdf"))
    p
  dev.off()
} else {
  DE1_MAST <- readRDS(paste0(Robject_dir, "DE_MAST_subpop_1_vs_2.RData"))
}
