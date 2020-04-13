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
analysis_mode <- args[2]
print(paste0("Analysis mode of normal InferCNV run = ", 
  analysis_mode, "_mode"))

project_name <- "thesis"
subproject_name <- "Figure_2.4_artefact_analysis"
sample_name <- "CID4520N"
print(paste0("sample name = ", sample_name))
analysis_mode <- "samples"
print(paste0("Analysis mode of normal InferCNV run = ", 
  analysis_mode, "_mode"))

print(paste0("Project name = ", project_name))
print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
#library(Seurat)
library(limma)
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

in_dir <- paste0(results_dir, "infercnv/", sample_name, "/", 
  analysis_mode, "_mode/Rdata/")

out_path <- paste0(results_dir, "GSEA/", sample_name, "/", 
  analysis_mode, "_mode/")

Robject_dir <- paste0(out_path, "Rdata/")
system(paste0("mkdir -p ", Robject_dir))
plot_dir <- paste0(out_path, "plots/")
system(paste0("mkdir -p ", plot_dir))
table_path <- paste0(out_path, "tables/")
system(paste0("mkdir -p ", table_path))

print(paste0("In directory = ", in_dir))
print(paste0("Robject directory = ", Robject_dir))
print(paste0("Plot dir = ", plot_dir))
print(paste0("Table directory = ", table_path))

print(paste0("Performing GSEA of artefacts from ", sample_name, "..."))


########################################################################
### 1. Load artefact data ###
########################################################################

artefact_record <- readRDS(paste0(in_dir, "1b.artefact_genes.Rdata"))


##############################################################
### 2. Artefact GSEA ###
##############################################################

artefact_GSEA <- lapply(artefact_record, function(x) {

  # find entrez IDs for gene list:
  egSYMBOL <- toTable(org.Hs.egSYMBOL)
  m <- match(x$genes, egSYMBOL$symbol)

  # determine genes not found in egSYMBOL and remove:
  genes_not_identified <- as.character(x$genes[which(is.na(m))])
  GSEA_genes <- as.character(x$genes[!(x$genes %in% genes_not_identified)])
  m <- match(GSEA_genes, egSYMBOL$symbol)
  GSEA_genes <- data.frame(
    gene_id = GSEA_genes,
    entrez_id = egSYMBOL$gene_id[m],
    stringsAsFactors = F
  )
  
  # perform GO:
  artefact_GO <- goana(GSEA_genes$entrez_id)
  top_GO <- topGO(artefact_GO)
  
  # perform KEGG:
  artefact_KEGG <- kegga(GSEA_genes$entrez_id)
  top_KEGG <- topKEGG(artefact_KEGG)
  
  return(
    list(
      genes = x,
      GO = top_GO,
      KEGG =  top_KEGG,
      removed_genes = genes_not_identified
    )
  )

})


##############################################################
### 2. Plot results ###
##############################################################

for (i in 1:length(artefact_GSEA)) {

  table_dir <- paste0(table_path, "artefact_", i, "/")
  system(paste0("mkdir -p ", table_dir))

  write.table(
    artefact_GSEA[[i]]$genes,
    paste0(table_dir, "indices_and_genes.txt"),
    sep = "\t",
    quote = F,
    col.names = T,
    row.names = F
  )

  write.table(
    artefact_GSEA[[i]]$GO,
    paste0(table_dir, "GO_results.txt"),
    sep = "\t",
    quote = F,
    col.names = T,
    row.names = F
  )

  write.table(
    artefact_GSEA[[i]]$KEGG,
    paste0(table_dir, "KEGG_results.txt"),
    sep = "\t",
    quote = F,
    col.names = T,
    row.names = F
  )

  write.table(
    artefact_GSEA[[i]]$removed_genes,
    paste0(table_dir, "removed_genes.txt"),
    sep = "\t",
    quote = F,
    col.names = T,
    row.names = F
  )

}



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






