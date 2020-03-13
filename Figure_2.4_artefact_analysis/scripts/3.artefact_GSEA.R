#! /share/ClusterShare/software/contrib/briglo/octR/src/R-3.6.0/builddir/bin/Rscript

###############################################################################
### Simulates cancer datasets taking as inputs number of CNVs, CNV lengths, ###
### gain/loss multipliers and downsamples required                          ###
###############################################################################

args = commandArgs(trailingOnly=TRUE)

project_name <- "thesis"
subproject_name <- "Figure_2.4_artefact_analysis"
sample_name <- args[1]
print(paste0("sample name = ", sample_name))
t_cells_included <- as.logical(args[2])
print(paste0("T-cells included in normal InferCNV run? ", 
  t_cells_included))
analysis_mode <- args[3]
print(paste0("Analysis mode of normal InferCNV run = ", 
  analysis_mode, "_mode"))

#project_name <- "thesis"
#subproject_name <- "Figure_2.4_artefact_analysis"
#sample_name <- "CID4520N"
#print(paste0("sample name = ", sample_name))
#t_cells_included <- TRUE
#print(paste0("T-cells included in normal InferCNV run? ", 
#  t_cells_included))
#analysis_mode <- "samples"
#print(paste0("Analysis mode of normal InferCNV run = ", 
#  analysis_mode, "_mode"))

print(paste0("Project name = ", project_name))
print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
#library(Seurat)
library(cowplot)
library(rlist, lib.loc = lib_loc)
#library(GenomicRanges)
#library(naturalsort, lib.loc = lib_loc)
#library(splatter, lib.loc = lib_loc)
library(ggplot2)
library(org.Hs.eg.db)
#library(data.table)
#library(ComplexHeatmap, lib.loc = lib_loc)
#library(circlize, lib.loc = lib_loc)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", 
  project_name, "/", subproject_name, "/")
#ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/functions/")
results_dir <- paste0(project_dir, "results/")

if (t_cells_included) {
  in_dir <- paste0(results_dir, "infercnv/t_cells_included/", 
    "/", sample_name, "/", analysis_mode, "_mode/tables/")
  out_path <- paste0(results_dir, GSEA, "/t_cells_included/", 
    "/", sample_name, "/", analysis_mode, "_mode/")
} else {
  in_dir <- paste0(results_dir, "infercnv/t_cells_excluded", 
    "/", sample_name, "/", analysis_mode, "_mode/tables/")
  out_path <- paste0(results_dir, GSEA, "/t_cells_excluded/", 
    "/", sample_name, "/", analysis_mode, "_mode/")
}

Robject_dir <- paste0(out_path, "Rdata/")
system(paste0("mkdir -p ", Robject_dir))
plot_dir <- paste0(out_path, "plots/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(out_path, "tables/")
system(paste0("mkdir -p ", table_dir))

print(paste0("Sample directory = ", in_dir))
print(paste0("Robject directory = ", Robject_dir))
print(paste0("Plot dir = ", plot_dir))
print(paste0("Table directory = ", table_dir))

print(paste0("Plotting InferCNV heatmap of ", sample_name))


##############################################################
### 0. Load data ###
##############################################################

infercnv_output <- as.data.frame(t(read.table(paste0(in_dir, 
  "/infercnv.12_denoised.observations.txt"), header=T, 
  as.is=T)))

na_less_vector <- unlist(infercnv_output)
na_less_vector <- na_less_vector[!is.na(na_less_vector)]

# generate heatmap from all output data:
all_chr_heatmap <- Heatmap(
  infercnv_output, name = paste0("hm"), 
  col = colorRamp2(c(min(na_less_vector), 1, 
  max(na_less_vector)), 
  c("#00106B", "white", "#680700"), space = "sRGB"),
  cluster_columns = F, cluster_rows = F,
  show_row_names = F, show_column_names = T,
  column_names_gp = gpar(col = "white"),
  show_row_dend = FALSE,
  use_raster = T, raster_device = c("png")
)
pdf(paste0(in_dir, "all_chr_heatmap.pdf"))
  print(all_chr_heatmap)
dev.off()

# isolate cluster with highest expression of artifacts and 
# generate heatmap:
cell_types <- read.table(paste0(input_files, 
	"metadata.txt"), header=F, as.is=T)
chosen_cells <- cell_types$V1[cell_types$V2 == "Epithelial_12"]
output_subset <- infercnv_output[chosen_cells,]

# generate heatmap from chosen cluster data:
all_chr_subset_heatmap <- Heatmap(
  output_subset, name = paste0("hm"), 
  col = colorRamp2(c(min(na_less_vector), 1, 
  max(na_less_vector)), 
  c("#00106B", "white", "#680700"), space = "sRGB"),
  cluster_columns = F, cluster_rows = F,
  show_row_names = F, show_column_names = T,
  column_names_gp = gpar(col = "white"),
  show_row_dend = FALSE,
  use_raster = T, raster_device = c("png")
)
pdf(paste0(in_dir, "all_chr_subset_heatmap.pdf"))
  print(all_chr_subset_heatmap)
dev.off()


##############################################################
### 1. GSEA on chromosome 12 artifact ###
##############################################################

# isolate chr12 from chosen cluster data:
gene_order <- read.table(paste0(ref_dir, 
	"/infercnv_gene_order.txt"))
chr12_df <- output_subset[colnames(output_subset) %in% 
	gene_order$V1[gene_order$V2 == "chr12"]]

chr12_heatmap <- Heatmap(
  chr12_df, name = paste0("hm"), 
  col = colorRamp2(c(min(na_less_vector), 1, 
  max(na_less_vector)), 
     c("#00106B", "white", "#680700"), space = "sRGB"),
  cluster_columns = F, cluster_rows = F,
  show_row_names = F, show_column_names = T,
  column_names_gp = gpar(col = "white"),
  show_row_dend = FALSE,
  use_raster = T, raster_device = c("png")
)
pdf(paste0(in_dir, "chr12_heatmap.pdf"))
  print(chr12_heatmap)
dev.off()

# isolate genes in gain artifact and generate heatmap:
chr12_gene_sums <- apply(chr12_df, 2, sum)
pdf(paste0(in_dir, "chr12_hist.pdf"))
  hist(chr12_gene_sums)
 dev.off()

chr12_gain_genes <- names(chr12_gene_sums)[chr12_gene_sums > 
  21.1]
chr12_gain_df <- chr12_df[,colnames(chr12_df) %in% 
  chr12_gain_genes]

chr12_gain_heatmap <- Heatmap(
  chr12_gain_df, name = paste0("hm"), 
  col = colorRamp2(c(min(na_less_vector), 1, 
  max(na_less_vector)), 
  c("#00106B", "white", "#680700"), space = "sRGB"),
  cluster_columns = F, cluster_rows = F,
  show_row_names = F, show_column_names = T,
  column_names_gp = gpar(col = "white"),
  show_row_dend = FALSE,
  use_raster = T, raster_device = c("png")
)
pdf(paste0(in_dir, "chr12_gain_heatmap.pdf"))
  print(chr12_gain_heatmap)
dev.off()

# find entrez IDs for gene list:
egSYMBOL <- toTable(org.Hs.egSYMBOL)
m <- match(chr12_gain_genes, egSYMBOL$symbol)
# determine alternate symbols for genes not found in egSYMBOL:
alt_symbols_needed <- chr12_gain_genes[which(is.na(m))]
chr12_gain_genes <- gsub("ATP5G2", "ATP5MC2", chr12_gain_genes)
chr12_gain_genes <- chr12_gain_genes[!(chr12_gain_genes %in% 
  c("RP11-983P16.4", "RP11-834C11.4"))]
m <- match(chr12_gain_genes, egSYMBOL$symbol)
chr12_gain_ids <- egSYMBOL$gene_id[m]

# perform GSEA:
chr12_GO <- goana(chr12_gain_ids)
top_chr12_GO <- topGO(chr12_GO)
write.table(top_chr12_GO, paste0(out_dir, 
  "/chr12_subset_GO_results.txt"))

# perform KEGG:
chr12_KEGG <- kegga(chr12_gain_ids)
top_chr12_KEGG <- topKEGG(chr12_KEGG)
write.table(top_chr12_KEGG, paste0(out_dir, 
  "/chr12_subset_KEGG_results.txt"))


##############################################################
### 2. GSEA on chromosome 17 artifact ###
##############################################################

# isolate chr17 from chosen cluster data:
chr17_df <- output_subset[colnames(output_subset) %in% 
	gene_order$V1[gene_order$V2 == "chr17"]]

chr17_heatmap <- Heatmap(
  chr17_df, name = paste0("hm"), 
  col = colorRamp2(c(min(na_less_vector), 1, 
  max(na_less_vector)), 
     c("#00106B", "white", "#680700"), space = "sRGB"),
  cluster_columns = F, cluster_rows = F,
  show_row_names = F, show_column_names = T,
  column_names_gp = gpar(col = "white"),
  show_row_dend = FALSE,
  use_raster = T, raster_device = c("png")
)
pdf(paste0(in_dir, "chr17_heatmap.pdf"))
  print(chr17_heatmap)
dev.off()

# isolate genes in gain artifact and generate heatmap:
chr17_gene_sums <- apply(chr17_df, 2, sum)
pdf(paste0(in_dir, "chr17_hist.pdf"))
  hist(chr17_gene_sums)
 dev.off()

chr17_gain_genes <- names(chr17_gene_sums)[chr17_gene_sums > 
  21.1]
chr17_gain_df <- chr17_df[,colnames(chr17_df) %in% 
  chr17_gain_genes]

chr17_gain_heatmap <- Heatmap(
  chr17_gain_df, name = paste0("hm"), 
  col = colorRamp2(c(min(na_less_vector), 1, 
  max(na_less_vector)), 
  c("#00106B", "white", "#680700"), space = "sRGB"),
  cluster_columns = F, cluster_rows = F,
  show_row_names = F, show_column_names = T,
  column_names_gp = gpar(col = "white"),
  show_row_dend = FALSE,
  use_raster = T, raster_device = c("png")
)
pdf(paste0(in_dir, "chr17_gain_heatmap.pdf"))
  print(chr17_gain_heatmap)
dev.off()

# find entrez IDs for gene list:
m <- match(chr17_gain_genes, egSYMBOL$symbol)
# determine alternate symbols for genes not found in egSYMBOL:
alt_symbols_needed <- chr17_gain_genes[which(is.na(m))]
alt_symbols <- c("CAVIN1", "RETREG3", "CENPX", "HEXD", "CYBC1")
for (i in 1:length(alt_symbols_needed)) {
  chr17_gain_genes <- gsub(alt_symbols_needed[i], alt_symbols[i], 
  chr17_gain_genes)
}
m <- match(chr17_gain_genes, egSYMBOL$symbol)
chr17_gain_ids <- unique(egSYMBOL$gene_id[m])

# perform GSEA:
chr17_GO <- goana(chr17_gain_ids)
top_chr17_GO <- topGO(chr17_GO)
write.table(top_chr17_GO, paste0(out_dir, 
  "/chr17_subset_GO_results.txt"))

# perform KEGG:
chr17_KEGG <- kegga(chr17_gain_ids)
top_chr17_KEGG <- topKEGG(chr17_KEGG)
write.table(top_chr17_KEGG, paste0(out_dir, 
  "/chr17_subset_KEGG_results.txt"))

# perform GSEA:
both_GO <- goana(c(chr12_gain_ids, chr17_gain_ids))
top_both_GO <- topGO(both_GO)
write.table(top_both_GO, paste0(out_dir, 
  "/chr12_and_17_subset_GO_results.txt"))

# perform KEGG:
both_KEGG <- kegga(c(chr12_gain_ids, chr17_gain_ids))
top_both_KEGG <- topKEGG(both_KEGG)
write.table(top_both_KEGG, paste0(out_dir, 
  "/chr12_and_17_subset_KEGG_results.txt"))


##############################################################
### 3. Visualise how many of the amplified genes are part of 
# the indicated pathways ###
##############################################################

library(biomaRt, 
  lib.loc="/share/ScratchGeneral/jamtor/R/3.5dev")
mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")

# Find how many chr12 amplified genes are part of 
# 'keratin filament' GO term:
keratin_filament_genes <- getBM(attributes=c('hgnc_symbol'),
  filters = 'go_parent_term', values = 'GO:0045095', 
  mart = mart)$hgnc_symbol
chr12_keratin_filament_genes <- 
  chr12_gain_genes[chr12_gain_genes %in% 
  keratin_filament_genes]

keratin_filament_results <- 
  data.frame(genes=keratin_filament_genes, 
  gene_in_artifact="no", stringsAsFactors=F)
keratin_filament_results$gene_in_artifact[
  keratin_filament_results$genes %in% 
  chr12_keratin_filament_genes] <- "yes"
keratin_filament_results <- 
  keratin_filament_results[
    order(keratin_filament_results$gene_in_artifact, 
    decreasing=T),
  ]
write.table(keratin_filament_results, paste0(out_dir,
  "chr12_gain_artifact_genes_in_keratin_filament_GO_pathway"),
  quote=F, col.names=T, row.names=F)


# Find how many chr17 amplified genes are part of 
# 'epithelial_cell_differentiation' GO term:
epithelial_differentiation_genes <- 
  getBM(attributes=c('hgnc_symbol'), 
  filters = 'go_parent_term', values = 'GO:0030855', 
  mart = mart)$hgnc_symbol
chr17_epithelial_differentiation_genes <- 
  chr17_gain_genes[chr17_gain_genes %in% 
  epithelial_differentiation_genes]

epithelial_differentiation_results <- 
  data.frame(genes=epithelial_differentiation_genes, 
  gene_in_artifact="no", stringsAsFactors=F)
epithelial_differentiation_results$gene_in_artifact[
  epithelial_differentiation_results$genes %in% 
  chr17_epithelial_differentiation_genes] <- "yes"
epithelial_differentiation_results <- 
  epithelial_differentiation_results[
    order(epithelial_differentiation_results$gene_in_artifact, 
    decreasing=T),
  ]
write.table(epithelial_differentiation_results, paste0(out_dir,
  "chr17_gain_artifact_genes_in_epithelial_differentiation_GO_pathway"),
  quote=F, col.names=T, row.names=F)

# Find how many chr17 amplified genes are part of 
# 'estrogen signalling pathway' KEGG term:
estrogen_signalling_kegg <- unlist(
  keggGet("hsa04915")[[1]]$GENE
)
estrogen_signalling_ids <- 
  estrogen_signalling_kegg[
    seq(1, length(estrogen_signalling_kegg), 2)
  ]
m <- match(estrogen_signalling_ids, egSYMBOL$gene_id)
estrogen_signalling_genes <- egSYMBOL$symbol[m]

chr17_estrogen_signalling_genes <- 
  chr17_gain_genes[chr17_gain_genes %in% estrogen_signalling_genes]

estrogen_signalling_results <- 
  data.frame(genes=estrogen_signalling_genes, 
  gene_in_artifact="no", stringsAsFactors=F)
estrogen_signalling_results$gene_in_artifact[
  estrogen_signalling_results$genes %in% 
  chr17_estrogen_signalling_genes] <- "yes"
estrogen_signalling_results <- 
  estrogen_signalling_results[
    order(estrogen_signalling_results$gene_in_artifact, 
    decreasing=T),
  ]
write.table(estrogen_signalling_results, paste0(out_dir,
  "chr17_gain_artifact_genes_in_estrogen_signalling_GO_pathway"),
  quote=F, col.names=T, row.names=F)


# Find how many chr17 amplified genes are part of 
# 'prolactin signalling pathway' KEGG term:
prolactin_signalling_kegg <- unlist(
  keggGet("hsa04917")[[1]]$GENE
)
prolactin_signalling_ids <- 
  prolactin_signalling_kegg[
    seq(1, length(prolactin_signalling_kegg), 2)
  ]
m <- match(prolactin_signalling_ids, egSYMBOL$gene_id)
prolactin_signalling_genes <- egSYMBOL$symbol[m]

chr17_prolactin_signalling_genes <- 
  chr17_gain_genes[chr17_gain_genes %in% prolactin_signalling_genes]

prolactin_signalling_results <- 
  data.frame(genes=prolactin_signalling_genes, 
  gene_in_artifact="no", stringsAsFactors=F)
prolactin_signalling_results$gene_in_artifact[
  prolactin_signalling_results$genes %in% 
  chr17_prolactin_signalling_genes] <- "yes"
prolactin_signalling_results <- 
  prolactin_signalling_results[
    order(prolactin_signalling_results$gene_in_artifact, 
    decreasing=T),
  ]
write.table(prolactin_signalling_results, paste0(out_dir,
  "chr17_gain_artifact_genes_in_prolactin_signalling_KEGG_pathway"),
  quote=F, col.names=T, row.names=F)


##############################################################
### 4. Randomly sample 100 bp regions of genome  100 times 
# and plot distribution of 'hits' of any genes in keratin 
# filament, epithelial cell differentiation, estrogen 
# or prolactin signalling pathways ###
##############################################################

keratin_filament_genes <- read.table(paste0(out_dir, 
  "chr12_gain_artefact_genes_in_keratin_filament_GO_pathway.txt"), 
  header=T, as.is=T)$genes
epithelial_pathway_genes <- read.table(paste0(out_dir, 
  "chr17_gain_artefact_genes_in_epithelial_differentiation_GO_pathway.txt"), 
  header=T, as.is=T)$genes
estrogen_signalling_genes <- read.table(paste0(out_dir, 
  "chr17_gain_artefact_genes_in_estrogen_signalling_KEGG_pathway.txt"), 
  header=T, as.is=T)$genes
prolactin_signalling_genes <- read.table(paste0(out_dir, 
  "chr17_gain_artefact_genes_in_prolactin_signalling_KEGG_pathway.txt"), 
  header=T, as.is=T)$genes
epithelial_estrogen_prolactin_genes <- unique(c(epithelial_pathway_genes, 
  estrogen_signalling_genes, prolactin_signalling_genes))

epithelial_genes <- list(keratin_filament_genes, 
  epithelial_pathway_genes, estrogen_signalling_genes, 
  prolactin_signalling_genes, epithelial_estrogen_prolactin_genes)

epithelial_pathway_names <- c("keratin_filament_genes", "epithelial_pathway_genes", 
                              "estrogen_signalling_genes", "prolactin_signalling_genes", 
                              "epithelial_estrogen_prolactin_genes")

infercnv_output_width <- ncol(infercnv_output)
set.seed(666)
random_ints <- sample.int(infercnv_output_width-99, 
  size = 8000, replace = F)

epithelial_path_hits <- lapply(epithelial_genes, function(x) {
  for ( i in 1:length(random_ints) ) {
    genes_in_window <- colnames(infercnv_output)[
        random_ints[i]:(random_ints[i]+99)
    ]
    if (i==1) {
      gene_hits <- c(length(
        genes_in_window[genes_in_window %in% x]
      ))
    } else {
      gene_hits[i] <- length(
        genes_in_window[genes_in_window %in% x]
      )
    }
  }


  return(gene_hits)
})
names(epithelial_pathway_names) <- epithelial_pathway_names

saveRDS(epithelial_path_hits, paste0(out_dir, "epithelial_pathway_random_sample_hits.rds"))
#epithelial_path_hits <- readRDS(paste0(out_dir, "epithelial_pathway_random_sample_hits.rds"))

plot_colours <- c("#b2182b", "#1F618D", "#053061",
                 "#bebada", "#f4a582", "#1c9099",
                 "#85929E", "#9B59B6", "#74add1",
                 "#1b7837", "#b8e186", "#fed976",
                 "#e7298a", "#18ffff", "#ef6c00",
                 "#A93226", "black","orange",
                 "#b8bc53", "#5628ce", "#fa909c",
                 "#8ff331","#270e26") 

i=1
barplots <- lapply(epithelial_path_hits, function(x) {
  df <- as.data.frame(table(x))
  df$percentage <- round(df$Freq/8000*100, 1)
  p <- ggplot(df, aes(x=x, y=Freq))
  p <- p + geom_bar(stat="identity", fill=plot_colours[2])
  p <- p + geom_text(aes(label=paste0(percentage, "%"), vjust=-0.5))
  p <- p + ggtitle(gsub("_", " ", epithelial_pathway_names[i]))
  p <- p + theme(plot.title = element_text(hjust = 0.5))
  p <- p + xlab("No. pathway genes in 100 bp window")
  p <- p + ylab("Proportion of trials")
  p <- p + theme_bw() + theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")
  )
  
  png(paste0(out_dir, "hits_in_", epithelial_pathway_names[i], 
    "_from_8000_randomly_sampled_100bp_genomic_windows.png"))
    print(p)
  dev.off()
  
  i <<- i+1
  return(p)
})






