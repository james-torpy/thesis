#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

project_name <- "thesis"
subproject_name <- "chapter_4"
args = commandArgs(trailingOnly=TRUE)
subcluster_method <- args[1]
subcluster_p <- args[2]
if (subcluster_p != "none") {
  subcluster_p <- as.numeric(subcluster_p)
}
coverage_filter <- args[3]
remove_artefacts <- args[4]
subset_data <- as.logical(args[5])
subset_samples <- as.logical(args[6])
na_colour <- args[7]
gene_proportion_threshold <- as.numeric(args[8])
QC_annot <- as.logical(args[9])
include_metabric <- as.logical(args[10])

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
subset_data <- FALSE
subset_samples <- FALSE
na_colour <- "white"
gene_proportion_threshold <- 0.5
QC_annot <- TRUE
include_metabric <- TRUE
out_suffix <- "18_tumours"

print(paste0("Subproject name = ", subproject_name))
print(paste0("Subclustered by ", subcluster_method))
print(paste0("Subclustered by pval ", subcluster_p))
print(paste0("Filter by coverage? ", coverage_filter))
print(paste0("Remove artefacts? ", remove_artefacts))

heatmap_prefix <- "combined_infercnv_heatmap"

if (subset_samples) {

  # determine order of samples by subtype:
  ER <- c("CID3941")
  HER2 <- c("CID3586")
  TNBC <- c("CID44041")
  sample_names <- c(ER, HER2, TNBC)

} else {

  if (out_suffix == "14_tumours") {
    # determine order of samples by subtype:
    ER <- c("CID4067", "CID4290A", "CID4463", "CID4530N", "CID4535")
    HER2 <- c("CID3921", "CID3586", "CID4066", "CID45171")
    TNBC <- c("CID44971", "CID44991", "CID4513", "CID4515", "CID4523")
    sample_names <- c(ER, HER2, TNBC)
  } else if (out_suffix == "18_tumours") {
    # determine order of samples by subtype:
    ER <- c("CID3941", "CID3948", "CID4067", "CID4290A",
      "CID4463", "CID4530N", "CID4535")
    HER2 <- c("CID3921", "CID3586", "CID3963", "CID4066", "CID45171")
    TNBC <- c("CID44041", "CID44971", "CID44991", "CID4513",
      "CID4515", "CID4523")
    sample_names <- c(ER, HER2, TNBC)
  }

}

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample names = ", sample_names))
print(paste0("Subset data? ", as.character(subset_data)))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(naturalsort, lib.loc = lib_loc)
library(cluster, lib.loc = lib_loc)
library(ComplexHeatmap, lib.loc=lib_loc)
library(circlize, lib.loc = lib_loc)
library(scales, lib.loc = lib_loc)
library(RColorBrewer)
library(reshape2)
library(fpc, lib.loc = lib_loc)
library(dplyr)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/",
  subproject_name, "/")
ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/functions/")
results_dir <- paste0(project_dir, "results/")
in_path <- paste0(results_dir, "infercnv/")
out_dir <- paste0(results_dir, "infercnv/combined_infercnv/", 
  coverage_filter, "/", subcluster_method, "/p_", subcluster_p, "/", 
  remove_artefacts, "/", out_suffix, "/")

if (subset_data & subset_samples) {
  Robject_dir <- paste0(out_dir, "Rdata_double_sub/")
  plot_dir <- paste0(out_dir, "plots_double_sub/")
} else if (subset_data) {
  Robject_dir <- paste0(out_dir, "Rdata_sub/")
  plot_dir <- paste0(out_dir, "plots_sub/")
} else if (subset_samples) {
  Robject_dir <- paste0(out_dir, "Rdata_sample_sub/")
  plot_dir <- paste0(out_dir, "plots_sample_sub/")
} else {
  Robject_dir <- paste0(out_dir, "Rdata/")
  plot_dir <- paste0(out_dir, "plots/")
}

system(paste0("mkdir -p ", Robject_dir))
system(paste0("mkdir -p ", plot_dir))

print(paste0("Reference directory = ", ref_dir))
print(paste0("R function directory = ", func_dir))
print(paste0("Output directory = ", out_dir))
print(paste0("R object directory = ", Robject_dir))
print(paste0("Plot directory = ", plot_dir))

print("Plotting combined CNV heatmap of: ")
print(sample_names)


################################################################################
### 0. Define functions and colours ###
################################################################################

fetch_chromosome_boundaries <- dget(paste0(func_dir, 
  "fetch_chromosome_boundaries.R"))
annotate_METABRIC_CNV_freq <- dget(paste0(func_dir, 
  "annotate_METABRIC_CNV_freq.R"))

subtype_cols <- read.table(
  paste0(ref_dir, "subtype_colour_palette.txt"),
  header = F,
  comment.char = "",
  stringsAsFactors = F
)
colnames(subtype_cols) <- c("subtype", "col")

all_cols <- read.table(
  paste0(ref_dir, "all_colour_palette.txt"),
  header = F,
  comment.char = "",
  stringsAsFactors = F
)[,1]

sample_cols <- all_cols[c(28:35, 41:43, 47:51, 53:78, 19:27)]


################################################################################
### 1. Load heatmap data and combine ###
################################################################################

print("Loading heatmap and metadata dfs for each sample...")

if (!file.exists(paste0(Robject_dir, "group_epi_heatmap.Rdata"))) {

  for (s in 1:length(sample_names)) {

    print(paste0("Loading ", sample_names[s], " heatmap df..."))

    sample_Robject_dir <- paste0(
  	  in_path, sample_names[s], "/", coverage_filter, "/", subcluster_method, 
      "/p_", subcluster_p, "/", remove_artefacts, 
      "/Rdata/"
    )

    epi_heatmap <- readRDS(
      paste0(
      	sample_Robject_dir, "5a.final_epithelial_heatmap_without_normals.Rdata"
      )
    )

    print(paste0("Loading ", sample_names[s], " metadata..."))

    epi_meta <- readRDS(
      paste0(
      	sample_Robject_dir, "5b.final_epithelial_metadata_without_normals.Rdata"
      )
    )
    epi_meta <- subset(
      epi_meta, 
      select = c(cell_ids, cell_type, nUMI, nGene, subcluster_id)
    )
    rownames(epi_meta) <- epi_meta$cell_ids

    # order by CNV subcluster:
    epi_meta <- epi_meta[
      naturalorder(epi_meta$subcluster_id),
    ]
    epi_meta$subcluster_id <- factor(
      epi_meta$subcluster_id,
      levels = naturalsort(unique(epi_meta$subcluster_id))
    )
    epi_heatmap <- epi_heatmap[epi_meta$cell_ids,]

    print(paste0(
      "Are epi_meta rownames in the same order as epi_heatmap? ",
      identical(rownames(epi_heatmap), rownames(epi_meta))
    ))
  
    if (s==1) {
      epi_heatmaps <- list(epi_heatmap)
      names(epi_heatmaps) <- sample_names[s]
      gene_list <- colnames(epi_heatmap)
      epi_metas <- list(epi_meta)
      names(epi_metas) <- sample_names[s]

    } else {
      epi_heatmaps[[s]] <- epi_heatmap
      names(epi_heatmaps)[s] <- sample_names[s]
      gene_list <- c(gene_list, colnames(epi_heatmap))
      epi_metas[[s]] <- epi_meta
      names(epi_metas)[s] <- sample_names[s]
    }

  }
  
  # get complete gene list from samples and order:
  print("Adding missing genes and collating heatmap and metadata dfs...")
  gene_order <- read.table(paste0(ref_dir, "infercnv_gene_order.txt"))
  gene_list <- unique(gene_list)
  gene_list <- as.character(gene_order$V1[gene_order$V1 %in% gene_list])

  # keep only genes in at least n% of samples:
  for (k in 1:length(epi_heatmaps)) {
  	if (k==1) {
  	  gene_present <- data.frame(
  	  	row.names = gene_list,
  	  	present = gene_list %in% colnames(epi_heatmaps[[k]])
  	  )
  	  colnames(gene_present)[k] <- names(epi_heatmaps)[k]
  	} else {
  	  gene_present$present <- gene_list %in% colnames(epi_heatmaps[[k]])
  	  colnames(gene_present)[k] <- names(epi_heatmaps)[k]
  	}
  }
  gene_freqs <- apply(gene_present, 1, function(x) length(which(x)))
  keep_genes <- names(gene_freqs[
  	gene_freqs >= floor(gene_proportion_threshold*length(epi_heatmaps))
  ])
  
  # add all missing genes as columns to all heatmap dfs:
  complete_epi_heatmaps <- lapply(epi_heatmaps, function(x) {
    
    missing_genes <- keep_genes[!(keep_genes %in% colnames(x))]
    missing_genes_df <- data.frame(matrix(NA, nrow = nrow(x), 
      ncol = length(missing_genes)))
    colnames(missing_genes_df) <- missing_genes
  
    complete_df <- cbind(x, missing_genes_df)
    m <- match(keep_genes, colnames(complete_df))
    complete_df <- complete_df[,m]
  
    return(complete_df)
  
  })
  
  # collate all dfs and check rows line up:
  group_epi_heatmap <- do.call(rbind, complete_epi_heatmaps)
  rownames(group_epi_heatmap) <- gsub("^.*\\.C", "C", rownames(group_epi_heatmap))
  print("Are all rows present in group heatmap df?")
  identical(nrow(group_epi_heatmap), sum(unlist(lapply(complete_epi_heatmaps, nrow))))
  print("Are all genes present in group heatmap df?")
  identical(colnames(group_epi_heatmap), keep_genes)
  
  # subset group_epi_heatmap if needed:
  if (subset_data) {
    group_epi_heatmap <- group_epi_heatmap[,1:300]
  }
  
  # collate heatmap metadata:
  group_epi_meta <- do.call("rbind", epi_metas)
  rownames(group_epi_meta) <- gsub(
    "^.*\\.", "", rownames(group_epi_meta)
  )

  print(paste0(
    "Are epi_meta rownames still in the same order as epi_heatmap?? ",
    identical(rownames(group_epi_heatmap), rownames(group_epi_meta))
  ))

  # select only group_epi_meta rows in group_epi_heatmap and order:
  group_epi_meta <- group_epi_meta[rownames(group_epi_heatmap),]
  m <- match(rownames(group_epi_heatmap), rownames(group_epi_meta))
  group_epi_meta <- group_epi_meta[m,]
  
  print(paste0(
    "Are epi_meta rownames still in the same order as epi_heatmap?? ",
    identical(rownames(group_epi_heatmap), rownames(group_epi_meta))
  ))

  # split heatmap by sample and remove any samples with <2 cells:
  sample_rownames <- split(
    rownames(group_epi_heatmap), gsub("_.*$", "", rownames(group_epi_heatmap))
  )
  split_heatmap <- lapply(sample_rownames, function(x) group_epi_heatmap[x,])
  for (m in 1:length(split_heatmap)) {
    if (nrow(split_heatmap[[m]]) < 2) {
      if (exists("remove_heatmaps")) {
        remove_heatmaps <- c(remove_heatmaps, names(split_heatmap)[m])
      } else {
        remove_heatmaps <- c(names(split_heatmap)[m])
      }
    }
  }
  if (exists("remove_heatmaps")) {
    split_heatmap <- split_heatmap[!(names(split_heatmap) %in% remove_heatmaps)]
    sample_names <- sample_names[!(sample_names %in% remove_heatmaps)]
    for (r in remove_heatmaps) {
      group_epi_meta <- group_epi_meta[
        grep(r, rownames(group_epi_meta), invert=T),
      ]
      group_epi_heatmap <- group_epi_heatmap[
        grep(r, rownames(group_epi_heatmap), invert=T),
      ]
    }
    print(paste0(
      "Are epi_meta rownames still in the same order as epi_heatmap?? ",
      identical(rownames(group_epi_heatmap), rownames(group_epi_meta))
    ))
  }


  ################################################################################
  ### 2. Add subtype and sample columns to group_epi_meta ###
  ################################################################################
  
  # add sample and subtype columns to epi_meta:
  group_epi_meta$sample <- gsub(
    "_.*$", "", group_epi_meta$cell_id
  )
  group_epi_meta$subtype <- NA
  group_epi_meta$subtype[
    group_epi_meta$sample %in% ER
  ] <- "ER"
  group_epi_meta$subtype[
    group_epi_meta$sample %in% HER2
  ] <- "HER2"
  group_epi_meta$subtype[
    group_epi_meta$sample %in% TNBC
  ] <- "TNBC"
  
  # reorder group_annotations_df by subtype and apply to group_epi_heatmap:
  group_epi_meta <- rbind(
    group_epi_meta[group_epi_meta$subtype=="ER",],
    group_epi_meta[group_epi_meta$subtype=="HER2",],
    group_epi_meta[group_epi_meta$subtype=="TNBC",]
  )
  
  group_epi_heatmap <- group_epi_heatmap[rownames(group_epi_meta),]
  
  print(paste0(
      "Are epi_meta rownames still in the same order as epi_heatmap?? ",
      identical(rownames(epi_heatmap), rownames(epi_meta))
  ))
  
  saveRDS(group_epi_heatmap, paste0(Robject_dir, "group_epi_heatmap.Rdata"))
  saveRDS(
    group_epi_meta, 
    paste0(Robject_dir, "group_epi_meta.Rdata")
  )

} else {

  group_epi_heatmap <- readRDS(paste0(Robject_dir, "group_epi_heatmap.Rdata"))
  group_epi_meta <- readRDS(
    paste0(Robject_dir, "group_epi_meta.Rdata")
  )
}
print(paste0("Are group_epi_heatmap and group_metadata_df rownames identical? ",
  identical(rownames(group_epi_heatmap), rownames(group_epi_meta))))


################################################################################
### 3. Create heatmap annotations ###
################################################################################

# order group_epi_meta first by subtype then sample:
group_epi_meta_split <- split(
  group_epi_meta, group_epi_meta$subtype
)
group_epi_meta_split <- list(
  ER=group_epi_meta_split$ER,
  HER2=group_epi_meta_split$HER2,
  TNBC=group_epi_meta_split$TNBC
)
ordered_group_epi_meta_split <- lapply(group_epi_meta_split, function(x) {
  return(x[order(x$sample),])
})
group_epi_meta <- do.call("rbind", ordered_group_epi_meta_split)
rownames(group_epi_meta) <- gsub(
  "^.*\\.", "", rownames(group_epi_meta)
)

# order group_epi_heatmap the same way as group_epi_meta:
group_epi_heatmap <- group_epi_heatmap[rownames(group_epi_meta),]
print(paste0("Are group_epi_heatmap and group_metadata_df rownames identical? ",
  identical(rownames(group_epi_heatmap), rownames(group_epi_meta))))

# create subtype annotation:
subtype_annotation_df <- subset(group_epi_meta, select = subtype)
subtype_annotation_df$subtype
subtype_annotation_df$subtype <- factor(subtype_annotation_df$subtype, 
  levels = c("ER", "HER2", "TNBC"))
print(paste0("Are subtype_annotation_df and group_epi_heatmap rownames identical? ",
  identical(rownames(group_epi_heatmap), rownames(subtype_annotation_df))))

# determine colours and create annotation:
m <- match(
  as.character(unique(subtype_annotation_df$subtype)),
  subtype_cols$subtype
)
subtype_annot_cols <- subtype_cols$col[m]
names(subtype_annot_cols) <- as.character(unique(subtype_annotation_df$subtype))

subtype_annotation <- Heatmap(
  as.matrix(subtype_annotation_df), 
  col = subtype_annot_cols, 
  name = "subtype_annotation", 
  width = unit(4, "mm"), 
  show_row_names = F, show_column_names = F,
  show_heatmap_legend = FALSE,
  row_split = group_epi_meta$subtype,
  row_gap = unit(1, "cm"),
  row_title = NULL
)

# create sample annotation:
sample_annotation_df <- subset(group_epi_meta, select = sample)
sample_annotation_df$sample <- factor(sample_annotation_df$sample,
  levels=unique(sample_annotation_df$sample))
print(paste0("Are sample_annotation_df and group_epi_heatmap rownames identical? ",
  identical(rownames(group_epi_heatmap), rownames(sample_annotation_df))))
# determine colours and create annotation:
sample_annot_cols <- sample_cols[1:length(unique(sample_annotation_df$sample))]
names(sample_annot_cols) <- sample_names
sample_annotation <- Heatmap(
  as.matrix(sample_annotation_df), 
  col = sample_annot_cols, 
  name = "sample_annotation", 
  width = unit(4, "mm"), 
  show_row_names = F, show_column_names = F,
  show_heatmap_legend = FALSE,
  row_split = group_epi_meta$subtype,
  row_gap = unit(1, "cm"),
  row_title = NULL
)

# create QC annotations:
nUMI_annotation <- rowAnnotation(
  correlation_annotation = anno_barplot(
    group_epi_meta$nUMI,
    gp = gpar(
      col = "#D8B72E", 
      width = unit(4, "cm")
    ), 
    border = FALSE, 
    which = "row", 
    axis = F
  ), show_annotation_name = FALSE
)
nUMI_annotation@name <- "nUMI"
nGene_annotation <- rowAnnotation(
  correlation_annotation = anno_barplot(
    group_epi_meta$nGene, name = "nGene",
    gp = gpar(
      col = "#9ECAE1", 
      width = unit(4, "cm")
    ), 
    border = FALSE, 
    which = "row", 
    axis = F
  ), show_annotation_name = FALSE
)
nGene_annotation@name <- "nGene"

# determine co-ordinates of vertical lines at chromosome borders:
chr_data <- fetch_chromosome_boundaries(group_epi_heatmap, ref_dir)
saveRDS(chr_data, paste0(Robject_dir, "chromosome_data.Rdata"))

# create heatmap annotation for METABRIC genome-wide CNV frequency:
METABRIC_CNV_freq <- read.table(
  paste0(ref_dir, "metabric_cnv_freq.txt"), 
  header=T, 
  as.is=T, 
  fill=T
)

print(paste0("Generating METABRIC CNV freq plot..."))
metabric_plot <- annotate_METABRIC_CNV_freq(
  group_epi_heatmap, 
  cnv_frequencies = METABRIC_CNV_freq,
  temp_subtype = "all", 
  chr_ends = chr_data$ends, 
  chr_lengths = chr_data$lengths,
  check_correlation = FALSE,
  func_dir = func_dir
)
dev.off()


#############################################################################
### 5. Prepare heatmap and annotations ###
################################################################################

if (!exists("split_heatmap")) {
  sample_rownames <- split(
    rownames(group_epi_heatmap), gsub("_.*$", "", rownames(group_epi_heatmap))
  )
  split_heatmap <- lapply(sample_rownames, function(x) group_epi_heatmap[x,])
}

# scale all samples the same:
print(paste0("Summaries of split heatmap df before rescaling to min/max values: "))
print(lapply(split_heatmap, function(x) summary(unlist(x))))

rescaled_split_heatmap <- lapply(split_heatmap, function(x) {
  return(as.data.frame(rescale(as.matrix(x), c(-1, 1))))
})

print(paste0("Summaries of split heatmap df after rescaling: "))
print(lapply(rescaled_split_heatmap, function(x) summary(unlist(x))))

rescaled_epi_heatmap <- do.call("rbind", rescaled_split_heatmap)
rownames(rescaled_epi_heatmap) <- gsub("^.*\\.", "", rownames(rescaled_epi_heatmap))
rescaled_epi_heatmap <- rescaled_epi_heatmap[rownames(group_epi_meta),]
print(paste0("Are rescaled_epi_heatmap and group_epi_meta rownames identical? ",
  identical(rownames(rescaled_epi_heatmap), rownames(group_epi_meta))))

saveRDS(rescaled_epi_heatmap, paste0(Robject_dir, "final_heatmap.Rdata"))

# prepare df for plotting:
plot_object <- rescaled_epi_heatmap
colnames(plot_object) <- rep("la", ncol(plot_object))

# define heatmap colours:
na_less_vector <- unlist(plot_object)
na_less_vector <- na_less_vector[!is.na(na_less_vector)]
heatmap_cols <- colorRamp2(c(min(na_less_vector), 0, max(na_less_vector)), 
      c("#00106B", "white", "#680700"), space = "sRGB")

print("Generating final heatmap...")

final_heatmap <- Heatmap( 
  as.matrix(plot_object),
  name = paste0("hm"),
  na_col = na_colour,
  col = heatmap_cols,
  cluster_columns = F, cluster_rows = F,
  show_row_names = F, show_column_names = T,
  column_names_gp = gpar(col = "white"),
  show_row_dend = FALSE,
  show_heatmap_legend = FALSE,
  use_raster = T, raster_device = c("png"),
  row_split = group_epi_meta$subtype,
  row_gap = unit(1, "cm"),
  row_title = NULL
)

# determine co-ordinates of horizontal lines at group borders:
spl_subtypes <- split(group_epi_meta, group_epi_meta$subtype)

hlines <- lapply(spl_subtypes, function(x) {

  spl_groups <- split(x$sample, x$sample)
  spl_groups <- spl_groups[unique(x$sample)]
  if (length(spl_groups) > 1) {
    for ( n in 1:(length(spl_groups)-1) ) {
      if (n==1) {
        hl <- c(length(spl_groups[[n]])/length(x$sample))
      } else {
        hl[n] <- hl[n-1] + length(spl_groups[[n]])/length(x$sample)
      }
    }
    hl <- 1-hl
  } else {
    hl <- c(length(spl_groups[[1]])/length(x$sample))
  }
  hl <- c(1, hl, 0)

  return(hl)

})

ht_list <- subtype_annotation + sample_annotation
ht_list <- ht_list + final_heatmap
ht_list <- ht_list + nUMI_annotation + nGene_annotation

annotated_heatmap <- grid.grabExpr(
  draw(ht_list, gap = unit(3, "mm"), heatmap_legend_side = "left")
)
dev.off()

# determine where starting co-ordinates for heatmap are based upon longest cluster name
# (0.00604 units per character):
longest_cluster_name <- max(nchar(unique(as.character(group_epi_meta$subtype))))
x_coord <- longest_cluster_name*0.0037

# generate heatmap legend:
signal_ranges <- round(range(unlist(plot_object), na.rm = TRUE), 2)
lgd <- Legend(
  at = c(signal_ranges[1], 0, signal_ranges[2]),
  col_fun = heatmap_cols, 
  title = "Scaled CNV\nvalue", 
  direction = "horizontal",
  grid_height = unit(3, "cm"),
  grid_width = unit(4.5, "cm"),
  legend_height = unit(3, "cm"),
  legend_width = unit(4.5, "cm"),
  labels_gp = gpar(fontsize = 20),
  title_gp = gpar(fontsize = 22, fontface = "plain")
)


###############################################################################
### 7. Plot heatmap and annotations ###
################################################################################

pdf(paste0(plot_dir, heatmap_prefix, "_rescaled.pdf"), height = 13, width = 24)

  grid.newpage()
    pushViewport(viewport(x = 0.11, y = 0.1, 
                        width = 0.89, height = 0.89, just = c("left", "bottom")))
      
      grid.draw(annotated_heatmap)

      # add lines to slice 1:
      decorate_heatmap_body("hm", {
          for ( e in 1:length(chr_data$end_pos) ) {
            grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), gp = gpar(lwd = 3, 
              col = "#383838"))
          }
          for (m in 1:length(hlines[[1]])) {
            grid.lines(
            	c(0, 1), 
            	c(hlines[[1]][m], hlines[[1]][m]), 
            	gp = gpar(lwd = 2, col = "#383838")
            )
          }
        }, slice = 1 )

      # add lines to slice 2:
      decorate_heatmap_body("hm", {
        for ( e in 1:length(chr_data$end_pos) ) {
          grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), gp = gpar(lwd = 3, 
            col = "#383838"))
        }
        for (m in 1:length(hlines[[2]])) {
          grid.lines(
          	c(0, 1), 
          	c(hlines[[2]][m], hlines[[2]][m]), 
          	gp = gpar(lwd = 2, col = "#383838")
          )
        }
      }, slice = 2 )

      decorate_heatmap_body("hm", {
          for ( e in 1:length(chr_data$end_pos) ) {
            grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), gp = gpar(lwd = 3, 
              col = "#383838"))
            grid.text(gsub("chr", "", names(chr_data$lab_pos)[e]), chr_data$lab_pos[e], 
              unit(0, "npc") + unit(-3.1, "mm"), gp=gpar(fontsize=18))
          }
          for ( m in 1:length(hlines[[3]]) ) {
            grid.lines(
            	c(0, 1), 
            	c(hlines[[3]][m], hlines[[3]][m]), 
            	gp = gpar(lwd = 2, col = "#383838")
            )
          }
        }, slice = 3
      )

    popViewport()

    # plot heatmap legend:
    pushViewport(viewport(x = unit(3.5, "cm"), y = unit(28.3, "cm"), width = unit(4.5, "cm"), 
      height = unit(4, "cm")))
      draw(lgd, x = unit(0.1, "cm"), y = unit(0.1, "cm"), just = c("left", "bottom"))
    popViewport()

    # plot subtype legend:
    pushViewport(viewport(x = unit(4, "cm"), y = unit(20, "cm"), width = unit(6, "cm"), 
      height = unit(10, "cm")))
      #grid.rect()
      # add title:
      pushViewport(viewport(x = 0.06, y = 1, width = unit(2, "cm"), height = unit(0.5, "cm")))
        grid.text("Subtypes", gp=gpar(fontsize=20), just = "left")
        #grid.rect()
      popViewport()
      for (l in 1:length(subtype_annot_cols)) {
        # add labels:
        pushViewport(viewport(
          x = 0.21, 
          y = 0.99-(0.08*l), 
          width = unit(2, "cm"), 
          height = unit(0.5, "cm")
        ))
          #grid.rect()
          grid.text(names(subtype_annot_cols)[l], gp=gpar(fontsize=18), just = "left")
        popViewport()
        # add dots:
        pushViewport(viewport(
          x = 0.1, 
          y = 0.99-(0.08*l), 
          width = unit(0.3, "cm"), 
          height = unit(0.3, "cm")
        ))
          grid.circle(gp=gpar(col = subtype_annot_cols[l], fill = subtype_annot_cols[l]))
        popViewport()
      }
    popViewport()

    # plot sample legend:
    pushViewport(viewport(x = unit(4, "cm"), y = unit(15.7, "cm"), width = unit(6, "cm"), 
      height = unit(10, "cm")))
      #grid.rect()
      # add title:
      pushViewport(viewport(x = 0.06, y = 1, width = unit(2, "cm"), height = unit(0.5, "cm")))
        grid.text("Sample IDs", gp=gpar(fontsize=20), just = "left")
        #grid.rect()
      popViewport()
      for (l in 1:length(sample_annot_cols)) {
        # add labels:
        pushViewport(viewport(
          x = 0.21, 
          y = 0.99-(0.08*l), 
          width = unit(2, "cm"), 
          height = unit(0.5, "cm")
        ))
          #grid.rect()
          grid.text(names(sample_annot_cols)[l], gp=gpar(fontsize=18), just = "left")
        popViewport()
        # add dots:
        pushViewport(viewport(
          x = 0.1, 
          y = 0.99-(0.08*l), 
          width = unit(0.3, "cm"), 
          height = unit(0.3, "cm")
        ))
          grid.circle(gp=gpar(col = sample_annot_cols[l], fill = sample_annot_cols[l]))
        popViewport()
      }
    popViewport()

    # plot METABRIC CNV frequencies:
    pushViewport(viewport(
      x = x_coord+0.103, 
      y = 0.091,
      width = 0.84, 
      height = 0.08, 
      just = c("left", "top")
    ))
      grid.draw(metabric_plot)
    popViewport()

    # add frequency label:
    pushViewport(viewport(
      x = 0.075,
      y = 0.05
    ))
      grid.text("METABRIC freq.", gp=gpar(fontsize=18))
    popViewport()

    # add QC annotation values:
    pushViewport(viewport(x=x_coord + 0.945, y=0.14, width = 0.1, height = 0.1, just = "top"))
      grid.text("nUMI", rot=65, gp=gpar(fontsize=20))
    popViewport()
    pushViewport(viewport(x=x_coord + 0.965, y=0.14, width = 0.1, height = 0.1, just = "top"))
      grid.text("nGene", rot=65, gp=gpar(fontsize=20))
    popViewport()

dev.off()

print(paste0("Heatmap created, output in ", plot_dir))

#convert pdf to png:
system(paste0("for p in ", plot_dir, "*.pdf; do echo $p; f=$(basename $p); echo $f; ",
"new=$(echo $f | sed 's/.pdf/.png/'); echo $new; ", 
"convert -density 150 ", plot_dir, "$f -quality 90 ", plot_dir, "$new; done"))


