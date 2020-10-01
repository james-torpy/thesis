#! /share/ClusterShare/software/contrib/briglo/octR/src/R-3.6.0/builddir/bin/Rscript

###############################################################################
### Simulates cancer datasets taking as inputs number of CNVs, CNV lengths, ###
### gain/loss multipliers and downsamples required                          ###
###############################################################################

args = commandArgs(trailingOnly=TRUE)

project_name <- "thesis"
subproject_name <- "Figure_2.6_artefact_analysis"
sample_name <- args[1]
print(paste0("sample name = ", sample_name))
analysis_mode <- args[2]
print(paste0("Analysis mode of normal InferCNV run = ", 
  analysis_mode, "_mode"))
remove_cell_types <- args[6]
print("Cell types removed:")
print(remove_cell_types)

#project_name <- "thesis"
#subproject_name <- "Figure_2.6_artefact_analysis"
#sample_name <- "CID4520N"
#print(paste0("sample name = ", sample_name))
#analysis_mode <- "samples"
#print(paste0("Analysis mode of normal InferCNV run = ", 
#  analysis_mode, "_mode"))
#remove_cell_types <- "No"
#print("Cell types removed:")
#print(remove_cell_types)

print(paste0("Project name = ", project_name))
print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(Seurat)
library(limma)
library(cowplot)
library(rlist, lib.loc = lib_loc)
library(ggplot2)
library(org.Hs.eg.db)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", 
  project_name, "/", subproject_name, "/")
#ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/functions/")
results_dir <- paste0(project_dir, "results/")

in_dir <- paste0(results_dir, "infercnv/", sample_name, "/", 
  remove_cell_types, "_removed/", analysis_mode, "_mode/Rdata/")
raw_dir <- paste0(project_dir, "raw_files/seurat_objects/", 
  sample_name, "/")
out_path <- paste0(results_dir, "GSEA/", sample_name, "/", 
  analysis_mode, "_mode/")

Robject_dir <- paste0(out_path, "Rdata/")
system(paste0("mkdir -p ", Robject_dir))
plot_dir <- paste0(out_path, "plots/")
system(paste0("mkdir -p ", plot_dir))
table_path <- paste0(out_path, "tables/")
system(paste0("mkdir -p ", table_path))

col_dir <- paste0(home_dir, "R/colour_palettes/")

print(paste0("In directory = ", in_dir))
print(paste0("Robject directory = ", Robject_dir))
print(paste0("Plot dir = ", plot_dir))
print(paste0("Table directory = ", table_path))

print(paste0("Finding cell types enriched for artefact gene expression in ", 
  sample_name, "..."))


###################################################################################
### 1. Load data ###
###################################################################################

artefact_indices <- readRDS(paste0(in_dir, "1a.artefact_indices.Rdata"))
artefact_only <- artefact_indices[artefact_indices$type != "neutral",]
artefact_by_gene <- readRDS(paste0(in_dir, "1b.artefact_by_gene.Rdata"))
seurat_object <- readRDS(paste0(raw_dir, "/04_seurat_object_filtered.Rdata"))


###################################################################################
### 2. Artefact Module Score ###
###################################################################################

# load filtered seurat object, take all genes of each artefact and use
# AddModuleScore to assess which clusters express this gene set the most,
# then remove any outlying clusters with high expression and run InferCNV

# group cell ids by cell type and combine B and plasma cells:
Idents(seurat_object) <- gsub("_[0-9].*", "", Idents(seurat_object))
Idents(seurat_object) <- gsub(
  "B_cells|Plasma_cells", "B_and_plasma_cells", 
  Idents(seurat_object)
)

cell_types <- split(
  names(Idents(seurat_object)), 
  as.character(Idents(seurat_object))
)

# set seed:
if (!file.exists(paste0(Robject_dir, "modscore_seeds.Rdata"))) {
  modscore_seeds <- sample(1:999, nrow(artefact_only))
  saveRDS(modscore_seeds, paste0(Robject_dir, "modscore_seeds.Rdata"))
} else {
  modscore_seeds <- readRDS(paste0(Robject_dir, "modscore_seeds.Rdata"))
}

# for each artefact, calculate mean modscore per cell type:
for (i in 1:nrow(artefact_only)) {

  artefact_genes <- names(artefact_by_gene)[
    artefact_only[i,]$start:artefact_only[i,]$end
  ]

  seurat_modscore <- AddModuleScore(
    object = seurat_object,
    features = list(as.character(artefact_genes)),
    seed = modscore_seeds[i]
  )

  if (i==1) {
    modscores <- list(seurat_modscore@meta.data$Cluster1)
    names(modscores[[i]]) <- rownames(seurat_modscore@meta.data)
  } else {
    modscores[[i]] <- seurat_modscore@meta.data$Cluster1
    names(modscores[[i]]) <- rownames(seurat_modscore@meta.data)
  }

  # plot distribution of modscores:
  if (!file.exists(
    paste0(plot_dir, "artefact_", i, "_modscore_density_plot.pdf")
  )) {
    modscore_density <- density(modscores[[i]])
    pdf(paste0(plot_dir, "artefact_", i, "_modscore_density_plot.pdf"))
      plot(modscore_density, main=NA, xlab = "modscore")
    dev.off()
  }

  # combine artefact mean modscores per cell:
  cell_type_score <- lapply(cell_types, function(x) {
    return(
      data.frame(
        artefact_chr = artefact_only$start_chr[i],
        mean_module_score = mean(modscores[[i]][x]),
        n = length(x),
        SE = ( sd(modscores[[i]][x]) )/( sqrt(length(x)) )
      )
    )
  })
  score_df <- do.call("rbind", cell_type_score)
  
  if (i==1) {
    modscore_df <- do.call("rbind", cell_type_score)
  } else {
    modscore_df <- rbind(
      modscore_df,
      do.call("rbind", cell_type_score)
    )
  }

}
modscore_df$cell_type <- gsub(
  "_", " ", gsub("[0-9]", "", rownames(modscore_df))
)

# plot modscores on barplot:
all_cols <- read.table(
  paste0(col_dir, "colour_palette_1.txt"), 
  as.is = T,
  comment.char = ""
)[,1]
cols <- all_cols[c(2, 1, 3:nrow(artefact_only))]

p <- ggplot(modscore_df, 
  aes(x = artefact_chr, y = mean_module_score, 
    colour = cell_type, fill = cell_type)
) 
p <- p + geom_bar(
  stat = "identity", position = position_dodge(0.8)
)
p <- p + geom_errorbar(
  aes(
    ymin=mean_module_score-SE, 
    ymax=mean_module_score+SE
  ),
  position = position_dodge(0.8),
  width = 0.6
)
p <- p + scale_color_manual(values = cols)
p <- p + scale_fill_manual(values = cols)
p <- p + ylab("Module score")
p <- p + xlab("Artefact location")
p <- p + theme(
  axis.title.x = element_text(size=24, margin = margin(t = 20, r = 0, b = 0, l = 0)),
  axis.text.x = element_text(size=22),
  axis.title.y = element_text(size=24, margin = margin(t = 0, r = 20, b = 0, l = 0)),
  axis.text.y = element_text(size=22),
  legend.title = element_blank(),
  legend.text = element_text(size=18),
  legend.key.size = unit(1, "cm")
)

# function to increase vertical spacing between legend keys
# @clauswilke
draw_key_polygon3 <- function(data, params, size) {
  lwd <- min(data$size, min(size) / 4)

  grid::rectGrob(
    width = grid::unit(0.6, "npc"),
    height = grid::unit(0.6, "npc"),
    gp = grid::gpar(
      col = data$colour,
      fill = alpha(data$fill, data$alpha),
      lty = data$linetype,
      lwd = lwd * .pt,
      linejoin = "mitre"
    ))
}

# register new key drawing function, 
# the effect is global & persistent throughout the R session
GeomBar$draw_key = draw_key_polygon3

pdf(
  paste0(plot_dir, "artefact_module_scores.pdf"),
  width = 15
)
  p
dev.off()
png(
  paste0(plot_dir, "artefact_module_scores.png"),
  width = 15,
  height = 10,
  res = 300,
  units = "in"
)
  p
dev.off()


###################################################################################
### 2. Identify cells with outlying positive expression of modules ###
###################################################################################

# fetch all module scores for each artefact/cell type:
tukey_results <- lapply(modscores, function(x) {

  # create df of module scores labelled by cell type:
  for (j in 1:length(cell_types)) {
    if (j==1) {
      mscore_df <- data.frame(
        cell_type = names(cell_types)[j],
        module_score = x[cell_types[[j]]]
      )
    } else {
      mscore_df <- rbind(
        mscore_df,
        data.frame(
          cell_type = names(cell_types)[j],
          module_score = x[cell_types[[j]]]
        )
      )
    }
  }

  # create linear model across all cell types:
  mscore_lm <- lm(module_score ~ cell_type, data = mscore_df)
  # perform ANOVA:
  mscore_av <- aov(mscore_lm)
  print(summary(mscore_av))

  # perform and return Tukey test:
  mscore_tukey <- TukeyHSD(mscore_av)

  return(
    list(
      ANOVA = summary(mscore_av),
      Tukey = mscore_tukey$cell_type
    )
  )

})

for (k in 1:length(tukey_results)) {

  saveRDS(
    tukey_results[[k]]$ANOVA, 
    paste0(Robject_dir, "artefact_", k, "ANOVA_results.Rdata")
  )

  table_dir <- paste0(table_path, "artefact_", k, "/")
  write.table(
    tukey_results[[k]]$Tukey,
    paste0(table_dir, "Tukey_results.txt"),
    sep = "\t",
    quote = F,
    col.names = T,
    row.names = T
  )

}

