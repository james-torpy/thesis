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
min_CNV_length <- as.numeric(args[6])
min_CNV_proportion <- as.numeric(args[7])
normals_removed <- as.logical(args[8])

#project_name <- "thesis"
#subproject_name <- "chapter_4"
#sample_name <- "CID44041"
#subcluster_method <- "random_trees"
#subcluster_p <- "0.05"
#if (subcluster_p != "none") {
#  subcluster_p <- as.numeric(subcluster_p)
#}
#coverage_filter <- "filtered"
#remove_artefacts <- "artefacts_not_removed"
#min_CNV_length <- 20
#min_CNV_proportion <- as.numeric("0.5")
#normals_removed <- TRUE

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg38, lib.loc = lib_loc)
library(ggplot2)
library(RColorBrewer)
library(ComplexHeatmap, lib.loc=lib_loc)
library(circlize, lib.loc = lib_loc)
library(naturalsort, lib.loc = lib_loc)
library(cowplot)
library(Seurat)
library(dplyr)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/",
  subproject_name, "/")
ref_dir <- paste0(project_dir, "/refs/")
func_dir <- paste0(project_dir, "/scripts/functions/")
results_dir <- paste0(project_dir, "results/")
in_path <- paste0(
  results_dir, "infercnv/", sample_name, "/",
  coverage_filter, "/", subcluster_method, "/p_", subcluster_p, "/", 
  remove_artefacts, "/"
)

Robject_dir <- paste0(in_path, "Rdata/")
system(paste0("mkdir -p ", Robject_dir))
plot_dir <- paste0(in_path, "plots/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(in_path, "tables/")
system(paste0("mkdir -p ", table_dir))

print(paste0("In path = ", in_path))
print(paste0("Reference directory = ", ref_dir))
print(paste0("R function directory = ", func_dir))
print(paste0("R object directory = ", Robject_dir))
print(paste0("Table directory = ", table_dir))
print(paste0("Plot directory = ", plot_dir))


################################################################################
### 0. Define functions ###
################################################################################

fetch_chromosome_boundaries <- dget(paste0(func_dir, 
  "fetch_chromosome_boundaries.R"))
create_extended_vector <- dget(paste0(func_dir, 
  "create_extended_vector.R"))
detect_CNVs <- dget(paste0(func_dir, "detect_CNVs.R"))
subtype_signal_plot <- dget(paste0(func_dir, "subtype_signal_plot.R"))
print_subcluster_signal <- dget(paste0(func_dir, "print_subcluster_signal.R"))

add_chromosome_indices <- function(indices_df, chr_positions) {
  indices_df$chr_start <- chr_positions$index[indices_df$start]
  indices_df$chr_end <- chr_positions$index[indices_df$end]
  return(indices_df)
}

add_overall_indices <- function(indices_df, chr_positions) {
  for (r in 1:nrow(indices_df)) {
    indices_df$start[r] <- chr_positions$overall_index[
      chr_positions$chr == indices_df$chr[r] &
      chr_positions$index == indices_df$chr_start[r]
    ]
  }
  for (r in 1:nrow(indices_df)) {
    indices_df$end[r] <- chr_positions$overall_index[
      chr_positions$chr == indices_df$chr[r] &
      chr_positions$index == indices_df$chr_end[r]
    ]
  }
  return(indices_df)
}

# create function to check whether vector is sequential:
is.sequential <- function(x){
  all(abs(diff(x)) == 1)
}

extra_colours <- c(brewer.pal(12, "Set3"), brewer.pal(12, "Paired"))
col_palette <- c(brewer.pal(8, "Dark2")[c(3:6,8)], brewer.pal(12, "Set3"), brewer.pal(8, "Accent"),
  "#660200", "#918940", "black", "#9ECAE1", "#3F007D", "#67000D", "#FDD0A2", "#08519C",
  "#DBB335", "#5AA050", "#807DBA", "#1A6000", "#F16913", "#FD8D3C", "#DEEBF7", 
  "#7F2704", "#DADAEB", "#FC9272", "#BCBDDC", extra_colours, "#b2182b", "#85929E", 
  "#9B59B6", "#74add1","#1b7837", "#b8e186", "#fed976","#e7298a", "#18ffff", "#ef6c00",
  "#A93226", "#E611ED","orange", "#b8bc53", "#5628ce", "#fa909c", "#8ff331","#270e26")
col_palette <- col_palette[-7]


################################################################################
### 1. Load and format data ###
################################################################################

print("Loading heatmap and metadata dfs...")

epithelial_heatmap <- readRDS(paste0(Robject_dir, 
  "/5a.final_epithelial_heatmap_without_normals.Rdata"))

if (normals_removed) {
  epithelial_metadata <- readRDS(paste0(Robject_dir, 
    "/5b.final_epithelial_metadata_without_normals.Rdata"))
} else {
  epithelial_metadata <- readRDS(paste0(Robject_dir, 
    "/5b.final_epithelial_metadata_with_normals.Rdata"))
}

print(paste0(
  "Are epithelial_metadata rownames in the same order as epithelial_heatmap?? ",
  identical(rownames(epithelial_heatmap), rownames(epithelial_metadata))
))

# reduce levels of subcluster ids:
epithelial_metadata$subcluster_id <- factor(
  epithelial_metadata$subcluster_id,
  levels = unique(
    as.character(epithelial_metadata$subcluster_id)
  )
)

# determine neutral CNV value as most frequent infercnv score and check:
if (!file.exists(paste0(Robject_dir, "score_table.Rdata"))) {
  score_table <- table(round(unlist(epithelial_heatmap), 6))
  saveRDS(score_table, paste0(Robject_dir, "score_table.Rdata"))
} else {
  score_table <- readRDS(paste0(Robject_dir, "score_table.Rdata"))
}

neutral_value <- as.numeric(
  names(score_table)[which.max(score_table)]
)
scores_around_neutral <- score_table[
  (grep(neutral_value, names(score_table))-5):
  (grep(neutral_value, names(score_table))+5)
]

# isolate each subpopulation matrix:
subpop_ids <- split(
  rownames(epithelial_heatmap), 
  epithelial_metadata$subcluster_id
)
subpop_matrices <- lapply(subpop_ids, function(x) {
  return(epithelial_heatmap[x,])
})

# save subpop signal data for future plots:
subpop_signal_data <- list(
  subpop_matrices = subpop_matrices,
  neutral_value = neutral_value
)
saveRDS(
  subpop_signal_data, 
  paste0(Robject_dir, "subpop_signal_data.Rdata")
)

# load gene co-ordinates and only keep chromosomes 1:22:
gene_coords <- read.table(
  paste0(ref_dir, "infercnv_gene_order.txt"),
  sep = "\t",
  header = F,
  stringsAsFactors = F,
  col.names = c("gene_id", "chr", "start", "end")
)
gene_coords <- gene_coords[!(gene_coords$chr %in% c("chrX", "chrY", "chrM")),]

# fetch genomic lengths for each chromosome:
hg38 <- BSgenome.Hsapiens.UCSC.hg38
chr_genomic <- list(
  lengths = seqlengths(hg38)[
    grep("_|M|X|Y", names(seqlengths(hg38)), invert = T)
  ]
)

# calculate starts and ends:
chr_names <- names(chr_genomic$lengths)
chr_genomic$ends <- cumsum(as.numeric(chr_genomic$lengths))
names(chr_genomic$ends) <- chr_names

chr_genomic$starts[1] <- 1
for (j in 2:length(chr_names)) {
  chr_genomic$starts[j] <- chr_genomic$ends[j-1]+1
}
names(chr_genomic$starts) <- chr_names

# add chromosome starts -1 onto gene coordinates:
coord_list <- split(gene_coords, gene_coords$chr)
coord_list <- lapply(coord_list, function(x) {
  x$start <- x$start + (chr_genomic$starts[unique(x$chr)])-1
  x$end <- x$end + (chr_genomic$starts[unique(x$chr)])-1
  return(x)
})
gene_coords <- do.call("rbind", coord_list)


################################################################################
### 2. Identify CNVs in each subpopulation and calculate genomic lengths ###
################################################################################

# fetch chromosome boundary co-ordinates:
if (!file.exists(paste0(Robject_dir, "chromosome_data.Rdata"))) {
  chr_data <- fetch_chromosome_boundaries(epithelial_heatmap, ref_dir)
  saveRDS(chr_data, paste0(Robject_dir, "chromosome_data.Rdata"))
} else {
  chr_data <- readRDS(paste0(Robject_dir, "chromosome_data.Rdata"))
}

# determine CNV co-ordinates and start and end genomic positions for each subpop:
if (!file.exists(paste0(Robject_dir, "CNV_indices_and_lengths.Rdata"))) {
  
  subpop_CNV_data <- lapply(subpop_matrices, detect_CNVs, min_CNV_proportion, 
    neutral_value)

  # label CNVs < min_CNV_length as artefacts:
  subpop_CNV_data <- lapply(subpop_CNV_data, function(x) {
    x$call <- as.character(x$call)
    x$call[x$call != "neutral" & x$length < min_CNV_length] <- 
      paste0(
        x$call[x$call != "neutral" & x$length < min_CNV_length],
        "_artefact"
      )
    return(x)
  })

  saveRDS(subpop_CNV_data, paste0(Robject_dir, "all_indices_and_lengths.Rdata"))
  
  # create CNV only list:
  subpop_CNV_only <- lapply(subpop_CNV_data, function(df) {
    return(df[df$call != "neutral",])
  })

  saveRDS(subpop_CNV_only, paste0(Robject_dir, "CNV_indices_and_lengths.Rdata"))
  
  # save mean, min and max CNV lengths:
  all_CNV_only <- do.call("rbind", subpop_CNV_only)
  CNV_lengths <- data.frame(
    metric = c(
      "gene_mean", "gene_min", "gene_max", 
      "genomic_mean", "genomic_min", "genomic_max"
    ),
    value = c(
      round(mean(all_CNV_only$length), 0), 
      min(all_CNV_only$length),
      max(all_CNV_only$length),
      round(mean(all_CNV_only$genomic_length), 1),
      min(all_CNV_only$genomic_length),
      max(all_CNV_only$genomic_length)
    )
  )
  write.table(
    CNV_lengths,
    paste0(table_dir, "/CNV_length_stats.txt"),
    quote = F,
    sep = "\t",
    row.names = F,
    col.names = T
  )

} else {
  subpop_CNV_data <- readRDS(
    paste0(Robject_dir, "all_indices_and_lengths.Rdata")
  )
  subpop_CNV_only <- readRDS(
    paste0(Robject_dir, "CNV_indices_and_lengths.Rdata")
  )
}


################################################################################
### 3. Convert to GRanges objects for comparison of subclusters ###
################################################################################

# annotate all genes with chromosome position:
for (l in 1:length(chr_data$lengths)) {
  if (l==1) {
    chrs <- rep(names(chr_data$lengths)[l], chr_data$lengths[l])
  } else {
    chrs <- c(chrs, rep(names(chr_data$lengths)[l], chr_data$lengths[l]))
  }
}
chr_positions <- data.frame(
  gene = colnames(epithelial_heatmap),
  chr = chrs
)
chr_positions <- do.call(
  "rbind",
  lapply(
    split(chr_positions, chr_positions$chr),
    function(x) {
      x$index <- 1:nrow(x)
      return(x)
    }
  )
)
chr_positions <- chr_positions[naturalorder(chr_positions$chr),]
chr_positions$overall_index <- 1:nrow(chr_positions)

# add chromosomal indices:
subpop_CNV_only <- lapply(subpop_CNV_only, add_chromosome_indices, chr_positions)


################################################################################
### 4. Plot CNV signal without CNV genomic length annotations ###
################################################################################

for (i in 1:length(subpop_matrices)) {
  print(i)

  signal_plot <- subtype_signal_plot(
    subpop_matrices[[i]],
    neutral_value,
    subpop_CNV_data[[i]]
  )

  if (i==1) {
    signal_plots <- list(signal_plot)
  } else {
    signal_plots[[i]] <- signal_plot
  }

}
names(signal_plots) <- names(subpop_matrices)

print_subcluster_signal(
  signal_plots,
  chr_data,
  plot_dir
)


################################################################################
### 5. Plot CNV signal with CNV genomic length annotations ###
################################################################################

for (i in 1:length(subpop_matrices)) {
  print(i)

  signal_plot <- subtype_signal_plot(
    subpop_matrices[[i]],
    neutral_value,
    subpop_CNV_data[[i]],
    include_lengths = TRUE
  )

  if (i==1) {
    signal_plots_with_length <- list(signal_plot)
  } else {
    signal_plots_with_length[[i]] <- signal_plot
  }

}
names(signal_plots_with_length) <- names(subpop_matrices)
dev.off()

print_subcluster_signal(
  signal_plots_with_length,
  chr_data,
  plot_dir,
  include_lengths = TRUE
)


#################################################################################
#### 6. Fetch CNV_associated DE and surrounding genes ###
#################################################################################
#
#if (file.exists(paste0(
#  table_dir, 
#  "met_and_CNA_assoc_subpop_DE.txt"
#))) {
#
#  # load DE genes:
#  DE_genes <- read.table(
#    paste0(
#      table_dir, 
#      "subpop_DE_min_pct_0.5_logfc_0.7_p_val_0.01_CNA_assoc_only.txt"
#    ),
#    header = TRUE
#  )
#  
#  # determine top DE genes:
#  top_genes <- DE_genes %>% 
#      group_by(cluster) %>% 
#      top_n(10, avg_logFC)
#  
#  # load cancer-associated DE genes:
#  cancer_DE_genes <- read.table(
#    paste0(
#      table_dir, 
#      "msigdb_cancer_subpop_DE_genes.txt"
#    ),
#    header = TRUE
#  )
#  cancer_DE_genes <- unique(as.character(cancer_DE_genes$gene))
#  
#  # keep top cancer DE genes only:
#  top_cancer <- top_genes[top_genes$gene %in% cancer_DE_genes,]
#  top_cancer <- as.data.frame(top_cancer)
#  top_cancer <- subset(top_cancer, select = c(gene, avg_logFC))
#  top_cancer$direction[top_cancer$avg_logFC > 0] <- "up"
#  top_cancer$direction[top_cancer$avg_logFC < 0] <- "down"
#  top_cancer <- top_cancer[order(top_cancer$gene),]
#  
#  # remove duplicates:
#  top_cancer$keep <- TRUE
#  for (j in 1:nrow(top_cancer)) {
#    if (j!=1) {
#    	if ( top_cancer$gene[j] == top_cancer$gene[j-1] & 
#    	  top_cancer$direction[j] == top_cancer$direction[j-1] ) {
#    	  top_cancer$keep[j] <- FALSE
#    	}
#    }
#  }
#  top_cancer <- top_cancer[top_cancer$keep,]
#  
#  # find closest matching gene in CNV heatmap:
#  DE_gene_list <- as.list(top_cancer$gene)
#  
#  closest_indices <- unlist(
#    lapply(DE_gene_list, function(x) {
#  
#      # if gene is in CNV data, record this co-ordinate:
#      if (x %in% colnames(epithelial_heatmap)) {
#        return(which(colnames(epithelial_heatmap) == x))
#      } else {
#    
#        # find where gene lies in gene coords annotation:
#        coord_ind <- which(gene_coords$gene_id == x)
#    
#        # iterate forward to find closest present gene:
#        for (i in (coord_ind+1):nrow(gene_coords)) {
#          if (gene_coords$gene[i] %in% colnames(epithelial_heatmap)) {
#            fwd_ind <- i
#            break()
#          }
#        }
#    
#        # iterate backwards to find closest present gene:
#        for (i in (coord_ind-1):1) {
#          if (gene_coords$gene[i] %in% colnames(epithelial_heatmap)) {
#            rev_ind <- i
#            break()
#          }
#        }
#    
#        # return closest index to the gene:
#        fwd_ind_diff <- fwd_ind - coord_ind
#        rev_ind_diff <- coord_ind - rev_ind
#        if (fwd_ind_diff <= rev_ind_diff) {
#          return(fwd_ind)
#        } else {
#          return(rev_ind)
#        }
#       
#      }
#    
#    })
#  
#  )
#  
#  DE_gene_ind <- data.frame(
#    gene = as.character(top_cancer$gene),
#    index = closest_indices,
#    direction = top_cancer$direction
#  )
#  
#  # stagger labels if too close to others:
#  DE_gene_ind <- DE_gene_ind[order(DE_gene_ind$index),]
#  DE_gene_ind$stagger <- 0
#  
#  stagger_labels <- function(df, stag_no) {
#    for (i in 1:length(df$index)) {
#      # for each gene, check all preceeding genes:
#      if (i!=1) {
#        if ( (df$index[i] - df$index[i-1]) < 400 ) {
#          df$stagger[i] <- stag_no
#        }
#      }
#    }
#    return(df)
#  }
#  DE_gene_ind <- stagger_labels(DE_gene_ind, 1)
#  
#  for (i in 2:10) {
#    split_df <- split(DE_gene_ind, DE_gene_ind$stagger)
#    split_df[[length(split_df)]] <- stagger_labels(split_df[[length(split_df)]], i)
#    DE_gene_ind <- do.call("rbind", split_df)
#    DE_gene_ind <- DE_gene_ind[order(DE_gene_ind$index),]
#  }
#  
#  
#  ################################################################################
#  ### 7. Plot CNV signal with CNA-assoc DE gene annotations ###
#  ################################################################################
#  
#  for (i in 1:length(subpop_matrices)) {
#    print(i)
#  
#    if (i==length(subpop_matrices)) {
#      signal_plot <- subtype_signal_plot(
#        subpop_matrices[[i]],
#        neutral_value,
#        subpop_CNV_data[[i]],
#        include_lengths = FALSE,
#        annotate_genes = DE_gene_ind
#      )
#    } else {
#      signal_plot <- subtype_signal_plot(
#        subpop_matrices[[i]],
#        neutral_value,
#        subpop_CNV_data[[i]],
#        include_lengths = FALSE
#      )
#    }
#  
#    if (i==1) {
#      signal_plots_with_DE <- list(signal_plot)
#    } else {
#      signal_plots_with_DE[[i]] <- signal_plot
#    }
#  
#  }
#  names(signal_plots_with_DE) <- names(subpop_matrices)
#  
#  print_subcluster_signal(
#    signal_plots_with_DE,
#    chr_data,
#    plot_dir,
#    include_lengths = FALSE,
#    annotate_genes = DE_gene_ind,
#    no_genes = ncol(epithelial_heatmap)
#  )
#  
#}
#  #