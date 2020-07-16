#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

#### Generate CNV heatmaps with following annotations: ###
# nUMI
# nGene
# Expression clusters
# CNV subclusters

#### Generate following t-SNEs/UMAPS: ###
# Best matching resolution clusters annotated by cell type
# Broad cell marker feature plot
# Epithelial vs myoepithelial feature plot
# Luminal vs basal feature plot

project_name <- "thesis"
subproject_name <- "Figure_3.3_individual_samples"
args = commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
min_CNV_length <- args[2]
min_CNV_proportion <- as.numeric(args[3])
loose_min_proportion <- as.numeric(args[4])
malignant_res <- args[5]
random_tree_p <- as.numeric(args[6])
subcluster_by <- args[7]
merge_similar_subclusters <- as.logical(args[8])

sample_name <- "CID4515"
min_CNV_length <- 20
min_CNV_proportion <- as.numeric("0.5")
loose_min_proportion <- as.numeric("0.4")
malignant_res <- "SUBSET_D_res.0.8"
random_tree_p <- 0.05
subcluster_by <- "random_trees"
merge_similar_subclusters <- FALSE

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

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/",
  subproject_name, "/")
ref_dir <- paste0(project_dir, "/refs/")
func_dir <- paste0(project_dir, "/scripts/functions/")
results_dir <- paste0(project_dir, "results/")
in_path <- paste0(results_dir, "infercnv/", sample_name, "/p_", 
  random_tree_p, "/clustered_by_", subcluster_by, "/")
raw_dir <- paste0(project_dir, "raw_files/")
seurat_dir <- paste0(raw_dir, "seurat_objects/", sample_name, "/")

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

add_chromosome_indices <- function(indices_df, chr_positions) {
  indices_df$chr_start <- chr_positions$index[indices_df$start]
  indices_df$chr_end <- chr_positions$index[indices_df$end]
  return(indices_df)
}

# create function to check whether vector is sequential:
is.sequential <- function(x){
  all(abs(diff(x)) == 1)
}

# create function to detect CNVs from subpopulations
detect_CNV <- function(df, min_proportion, neutral_value) {
  
  # if at least 50% of cells have signal above or below neutral_value, 
  # label gene as loss/gain respectively:
  for (c in 1:ncol(df)) {
  
    # determine if loss:
    if ( length(which(
        round(df[,c], 6) < neutral_value
      )) >= min_proportion*length(df[,c]) ) {
  
      if (c==1) {
        call_vector <- c("loss")
      } else {
        call_vector[c] <- "loss"
      }
  
    } else if ( length(which(
        round(df[,c], 6) > neutral_value
      )) >= min_proportion*length(df[,c]) ) {
  
      if (c==1) {
        call_vector <- c("gain")
      } else {
        call_vector[c] <- "gain"
      }
  
    } else {
  
      if (c==1) {
        call_vector <- c("neutral")
      } else {
        call_vector[c] <- "neutral"
      }
  
    }
  
  }
  
  # label stretches of losses/gains as CNVs:
  # fetch indices of CNVs
  CNV_indices <- data.table(
    call = rle(call_vector)$values,
    length = rle(call_vector)$lengths
  )
  for (r in 1:nrow(CNV_indices)) {
  
    if (r==1) {
      CNV_indices$start[r] <- 1
      CNV_indices$end[r] <- CNV_indices$length[r]
    } else {
      CNV_indices$start[r] <- CNV_indices$end[r-1]+1
      CNV_indices$end[r] <- CNV_indices$start[r] + (CNV_indices$length[r]-1)
    }
  
  }
  
  # find chromosomes CNVs belong to:
  for (k in 1:length(chr_data$ends)) {
    if (k==1) {

      CNV_indices$start_chr[
        CNV_indices$start <= chr_data$ends[k]
      ] <- names(chr_data$ends)[k]

      CNV_indices$end_chr[
        CNV_indices$end <= chr_data$ends[k]
      ] <- names(chr_data$ends)[k]

    } else {

      CNV_indices$start_chr[
        CNV_indices$start <= chr_data$ends[k] & 
        CNV_indices$start > chr_data$ends[k-1]
      ] <- names(chr_data$ends)[k]

      CNV_indices$end_chr[
        CNV_indices$end <= chr_data$ends[k] & 
        CNV_indices$end > chr_data$ends[k-1]
      ] <- names(chr_data$ends)[k]

    }
  }

  # make CNVs overlapping 2 chromosomes two distinct ranges:
  CNV_indices$keep <- TRUE
  for (n in 1:nrow(CNV_indices)) {
    if (CNV_indices$call[n] != "neutral" & 
        CNV_indices$start_chr[n] != CNV_indices$end_chr[n]) {
      # mark row for removal:
      CNV_indices$keep[n] <- FALSE
      # identify end of first chromosome:
      chr_end <- chr_data$ends[CNV_indices$start_chr[n]]
      new_rows <- data.frame(
        call = rep(CNV_indices$call[n], 2),
        length = rep(NA, 2),
        start = c(CNV_indices$start[n], chr_end+1),
        end = c(chr_end, CNV_indices$end[n]),
        start_chr = c(
          as.character(CNV_indices$start_chr[n]),
          as.character(CNV_indices$end_chr[n])
        ),
        end_chr = c(
          as.character(CNV_indices$start_chr[n]), 
          as.character(CNV_indices$end_chr[n])
        ),
        keep = rep(TRUE, 2)
      )
      new_rows$length <- new_rows$end-new_rows$start+1
      CNV_indices <- rbind(
        CNV_indices,
        new_rows
      )
    }
  }
  # remove CNV indices overlapping 2 chromosomes:
  CNV_indices <- CNV_indices[CNV_indices$keep,]
  CNV_indices <- subset(CNV_indices, select=-c(keep, end_chr))
  colnames(CNV_indices) <- gsub("start_chr", "chr", colnames(CNV_indices))

  # identify start and end genes of each CNV:
  CNV_indices$start_gene <- colnames(df)[CNV_indices$start]
  CNV_indices$end_gene <- colnames(df)[CNV_indices$end]
  
  # calculate genomic start of CNVs:
  merge_coords <- subset(gene_coords, select = c(gene_id, start))
  colnames(merge_coords) <- c("start_gene", "genomic_start")
  CNV_indices <- merge(
    CNV_indices,
    merge_coords,
    by = "start_gene"
  )
  # calculate genomic start of CNVs:
  merge_coords <- subset(gene_coords, select = c(gene_id, end))
  colnames(merge_coords) <- c("end_gene", "genomic_end")
  CNV_indices <- merge(
    CNV_indices,
    merge_coords,
    by = "end_gene"
  )
  # calculate genomic length of CNVs in mb:
  CNV_indices$genomic_length <- 
    round((CNV_indices$genomic_end-CNV_indices$genomic_start)/1000000, 1)

  # calculate midpoints for labelling:
  CNV_indices$midpoints <- 
    CNV_indices$start + floor((CNV_indices$end-CNV_indices$start)/2)
  
  return(CNV_indices[order(CNV_indices$start),])

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
  "/3a.filtered_heatmap_good_coverage_only.Rdata"))
epithelial_metadata <- readRDS(paste0(Robject_dir, 
  "/3b.filtered_metadata_good_coverage_only.Rdata"))
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
  
  subpop_CNV_data <- lapply(subpop_matrices, detect_CNV, min_CNV_proportion, 
    neutral_value)
  
  # create CNV only list:
  subpop_CNV_only <- lapply(subpop_CNV_data, function(df) {
    CNV_only <- df[df$call != "neutral",]
    # remove CNVs less than min_CNV_length:
    return(CNV_only[CNV_only$length >= min_CNV_length,])
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

# add chromosomal indices:
subpop_CNV_only <- lapply(subpop_CNV_only, add_chromosome_indices, chr_positions)

# convert all CNV data to genomic ranges objects:
subpop_CNV_gr <- lapply(subpop_CNV_only, function(x) {
  return(
    GRanges(
      seqnames = Rle(x$chr),
      ranges = IRanges(
        start = as.integer(x$chr_start),
        end = as.numeric(x$chr_end)
      ),
      strand = Rle("*"),
      call = x$call,
      length = x$length,
      genomic_start = x$genomic_start,
      genomic_end = x$genomic_end,
      genomic_length = x$genomic_length,
      midpoints = x$midpoints
    )
  )
})


################################################################################
### 4. Compare CNV profiles across subclusters and merge if no difference ###
################################################################################

#save.image(paste0(Robject_dir, "pre_CNV_verification.Rdata"))

#home_dir <- "/share/ScratchGeneral/jamtor/"
#project_dir <- paste0(home_dir, "projects/", project_name, "/",
#  subproject_name, "/")
#results_dir <- paste0(project_dir, "results/")
#in_path <- paste0(results_dir, "infercnv/", sample_name, "/")
#Robject_dir <- paste0(in_path, "Rdata/")
#load(paste0(Robject_dir, "pre_CNV_verification.Rdata"))
if (merge_similar_subclusters) {
  if (!file.exists(
  paste0(Robject_dir, "/4.subcluster_revised_epithelial_metadata.Rdata")
  )) {
  
    # for each subpopulation, findOverlaps looser CNV profiles of all others
    # and merge if subpop has CNVs not contained in looser profiles of others
    epithelial_metadata$revised_subcluster_id <- epithelial_metadata$subcluster_id
  
    # take each subpopulation CNV indices as query and overlap it with all other subpopulation
    # indices:
    for (q in 1:length(subpop_CNV_gr)) {
      for (s in 1:length(subpop_CNV_gr)) {
  
        olaps <- suppressWarnings(findOverlaps(subpop_CNV_gr[[q]], subpop_CNV_gr[[s]]))
  
        # if queryHits contains all query entries, record as merge:
        if ( identical(unique(queryHits(olaps)), 1:length(subpop_CNV_gr[[q]])) ) {
  
          if (!exists("merge_df")) {
            merge_df <- data.frame(
              query = names(subpop_CNV_gr)[q],
              subject = names(subpop_CNV_gr)[s],
              stringsAsFactors = F
            )
          } else {
            merge_df <- rbind(
              merge_df,
              data.frame(
                query = names(subpop_CNV_gr)[q],
                subject = names(subpop_CNV_gr)[s],
                stringsAsFactors = F
              )
            )
          }
           
         
        } else {
         
          whole_range <- 1:length(subpop_CNV_gr[[q]])
          mismatch_gr <- subpop_CNV_gr[[q]][setdiff(whole_range, queryHits(olaps))]
  
          if (!exists("mismatch_df")) {
            mismatch_df <- data.frame(
              in_subpop = names(subpop_CNV_gr)[q],
              not_in_subpop = names(subpop_CNV_gr)[s],
              chr = seqnames(mismatch_gr),
              start = start(ranges(mismatch_gr)),
              end = end(ranges(mismatch_gr)),
              call = mismatch_gr$call,
              genomic_length = mismatch_gr$genomic_length
            )
          } else {
            mismatch_df <- rbind(
              mismatch_df,
              data.frame(
                in_subpop = names(subpop_CNV_gr)[q],
                not_in_subpop = names(subpop_CNV_gr)[s],
                chr = seqnames(mismatch_gr),
                start = start(ranges(mismatch_gr)),
                end = end(ranges(mismatch_gr)),
                call = mismatch_gr$call,
                genomic_length = mismatch_gr$genomic_length
              )
            )
          }
          
        }
  
      }
    }
    # remove rows of merge_df in which subpopulations have been compared to themselves:
    merge_df <- merge_df[!(merge_df$query == merge_df$subject),]
  
    # save merge_df and merge subpopulations with all CNVs in common:
    if (nrow(merge_df) > 0) {
  
      write.table(
        merge_df,
        paste0(table_dir, "merged_subpopulations.txt"),
        sep = "\t",
        quote = F,
        row.names = F,
        col.names = T
      )
  
      for (r in 1:nrow(merge_df)) {
  
        # if both row and reverse of row is present in df, merge:
        if (nrow(merge(rev(merge_df[r,]), merge_df)) > 0) {
          if (merge_df[r,]$query %in% epithelial_metadata$revised_subcluster_id) {
  
            print(paste0("Merging ", merge_df[r,]$query, " and ", merge_df[r,]$subject, 
              " due to no detected difference between CNV profiles..."))
  
            epithelial_metadata$revised_subcluster_id[
              epithelial_metadata$revised_subcluster_id == merge_df[r,]$subject
            ] <- merge_df[r,]$query
          }
        }
        
      }
    }
    # remove any merged subpopulations secondary cluster ID:
    subpop_CNV_gr <- subpop_CNV_gr[
      names(subpop_CNV_gr) %in% epithelial_metadata$revised_subcluster_id
    ]
    
    saveRDS(epithelial_metadata, paste0(Robject_dir, 
      "/4.subcluster_revised_epithelial_metadata.Rdata"))
  
    if (exists("mismatch_df")) {
      write.table(
        mismatch_df,
        paste0(table_dir, "subpopulation_differential_CNVs.txt"),
        sep = "\t",
        quote = F,
        row.names = F,
        col.names = T
      )
    }
  
  } else {
    epithelial_metadata <- readRDS(paste0(Robject_dir, 
      "/4.subcluster_revised_epithelial_metadata.Rdata"))
  }
}



################################################################################
### 5. Add annotations and prepare heatmap ###
################################################################################

if (merge_similar_subclusters) {
  if (!file.exists(paste0(Robject_dir, "cluster_colours.Rdata"))) {

    seurat_malignant <- readRDS(paste0(seurat_dir, "05_seurat_object_malignant.Rdata"))
    Idents(seurat_malignant) <- paste0(
      seurat_malignant@meta.data$garnett_call_ext_major,
      " ",
      eval(parse(text = paste0("seurat_malignant@meta.data$", malignant_res)))
    )
    levels(Idents(seurat_malignant)) <- naturalsort(levels(Idents(seurat_malignant)))
    
    # define cluster annotation colours:
    cluster_number <- length(unique(Idents(seurat_malignant)))
    cluster_cols <- col_palette[1:cluster_number]
    names(cluster_cols) <- unique(Idents(seurat_malignant))
  
    saveRDS(cluster_cols, paste0(Robject_dir, "cluster_colours.Rdata"))
  
  } else {
    cluster_cols <- readRDS(paste0(Robject_dir, "cluster_colours.Rdata"))
  }
  
  # create group annotation:
  group_annotation_df <- subset(epithelial_metadata, select = cell_type)
  group_annotation <- Heatmap(
    as.matrix(group_annotation_df), 
    col = cluster_cols, 
    name = "Expression\nclusters", 
    width = unit(4, "mm"), 
    show_row_names = F, show_column_names = F,
    show_heatmap_legend = T
  )
  # create QC annotations:
  nUMI_annotation <- rowAnnotation(
    nUMI = anno_barplot(
      epithelial_metadata$nUMI,
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
    nGene = anno_barplot(
      epithelial_metadata$nGene, name = "nGene",
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
  
  # create subcluster annotation:
  if (merge_similar_subclusters) {
    subcluster_annot_df <- subset(epithelial_metadata, select = revised_subcluster_id)
    colnames(subcluster_annot_df) <- gsub(
    	"revised_subcluster_id", "subcluster_df", colnames(subcluster_annot_df)
    )
  } else {
    subcluster_annot_df <- subset(epithelial_metadata, select = subcluster_id)
  }
  
  # change subcluster names to be sequential:
  subcluster_annot_df$subcluster_id <- factor(subcluster_annot_df$subcluster_id)
  levels(subcluster_annot_df$subcluster_id) <- paste0(
    "CNV ", 
    1:length(unique(subcluster_annot_df$subcluster_id))
  )
  subcluster_cols <- rev(col_palette)[
    1:length(unique(subcluster_annot_df$subcluster_id))
  ]
  names(subcluster_cols) <- unique(subcluster_annot_df$subcluster_id)
  subcluster_annotation <- Heatmap(as.matrix(subcluster_annot_df), 
    col = subcluster_cols, 
    name = "subcluster_annotation", width = unit(6, "mm"), 
    show_row_names = F, show_column_names = F, 
    show_heatmap_legend = T,
    heatmap_legend_param = list(labels_gp = gpar(fontsize = 12))
  )
  
  # fetch chromosome boundary co-ordinates:
  if (!file.exists(paste0(Robject_dir, "chromosome_data.Rdata"))) {
    chr_data <- fetch_chromosome_boundaries(epithelial_heatmap, ref_dir)
    saveRDS(chr_data, paste0(Robject_dir, "chromosome_data.Rdata"))
  } else {
    chr_data <- readRDS(paste0(Robject_dir, "chromosome_data.Rdata"))
  }
  
  # prepare df for plotting:
  plot_object <- epithelial_heatmap
  colnames(plot_object) <- rep("la", ncol(plot_object))
  # define heatmap colours:
  na_less_vector <- unlist(plot_object)
  na_less_vector <- na_less_vector[!is.na(na_less_vector)]
  heatmap_cols <- colorRamp2(c(min(na_less_vector), 1, max(na_less_vector)), 
        c("#00106B", "white", "#680700"), space = "sRGB")
  
  print("Generating final heatmap...")
  # create main CNV heatmap:
  final_heatmap <- Heatmap(
    as.matrix(plot_object), name = paste0("hm"), 
    col = heatmap_cols,
    cluster_columns = F, cluster_rows = F,
    show_row_names = F, show_column_names = T,
    column_names_gp = gpar(col = "white"),
    show_row_dend = F,
    show_heatmap_legend = F,
    heatmap_legend_param = list(labels_gp = gpar(col = "red", fontsize = 12)),
    use_raster = T, raster_device = c("png")
  )
  
  ht_list <- subcluster_annotation + group_annotation + final_heatmap + 
    nUMI_annotation + nGene_annotation
  
  annotated_heatmap <- grid.grabExpr(
    draw(ht_list, gap = unit(6, "mm"), heatmap_legend_side = "left")
  )
  dev.off()
  
  
  ################################################################################
  ### 6. Plot heatmap ###
  ################################################################################
  
  # determine where starting co-ordinates for heatmap are based upon longest cluster name
  # (0.00604 units per character):
  longest_cluster_name <- max(nchar(unique(as.character(epithelial_metadata$cell_type))))
  x_coord <- longest_cluster_name*0.0037
  
  # generate heatmap legend:
  signal_ranges <- round(range(unlist(plot_object)), 2)
  lgd <- Legend(
    at = c(signal_ranges[1], 1, signal_ranges[2]),
    col_fun = heatmap_cols, 
    title = "CNV signal", 
    direction = "horizontal",
    grid_height = unit(2.5, "cm"),
    grid_width = unit(0.1, "cm"),
    labels_gp = gpar(fontsize = 16),
    title_gp = gpar(fontsize = 22, fontface = "plain")
  )
  
  # plot final annotated heatmap:
  png(paste0(plot_dir, "infercnv_plot_verified_subclusters.png"), 
    height = 13, width = 20, res = 300, units = "in") 
  
    grid.newpage()
      pushViewport(viewport(x = 0.155, y = 0.065, width = 0.752, height = 0.78, 
        just = c("left", "bottom")))
        grid.draw(annotated_heatmap)
        decorate_heatmap_body("hm", {
          for ( e in 1:length(chr_data$end_pos) ) {
            grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), 
              gp = gpar(lwd = 3, col = "#383838"))
            if (e==1) {
              grid.text(gsub("chr", "", names(chr_data$lab_pos)[e]), chr_data$lab_pos[e], 
              unit(0, "npc") + unit(-3.5, "mm"), gp=gpar(fontsize=20))
            } else if (e==21) {
              grid.text(paste0("\n", gsub("chr", "", names(chr_data$lab_pos)[e])), chr_data$lab_pos[e], 
              unit(0, "npc") + unit(-3.5, "mm"), gp=gpar(fontsize=20))
            } else {
              grid.text(gsub("chr", "", names(chr_data$lab_pos)[e]), chr_data$lab_pos[e], 
              unit(0, "npc") + unit(-3.5, "mm"), gp=gpar(fontsize=20))        
            }
          }
        })
      popViewport()
      # plot legend:
      pushViewport(viewport(x = unit(2, "cm"), y = unit(15, "cm"), width = unit(0.1, "cm"), 
        height = unit(0.4, "cm"), just = c("right", "bottom")))
        draw(lgd, x = unit(0.1, "cm"), y = unit(0.1, "cm"), just = c("left", "bottom"))
      popViewport()
      # label annotations:
      pushViewport(viewport(x=x_coord + 0.807, y=0.0001, width = 0.1, height = 0.1, 
        just = "bottom"))
        grid.text("nUMI", rot=65)
      popViewport()
      pushViewport(viewport(x=x_coord + 0.825, y=0.0001, width = 0.1, height = 0.1, 
        just = "bottom"))
        grid.text("nGene", rot=65)
      popViewport()
      
  dev.off()
  
  print(paste0("Heatmap created, output in ", plot_dir))

}



#################################################################################
#### 5. Identify CNVs common across all subpopulations to fade them out ###
#################################################################################
#common_CNV_gr <- lapply(subpop_CNV_gr, function(x) {
#  for (k in 1:length(subpop_CNV_gr)) {
#    # keep only ranges overlapping with all subpopulations:
#    olap <- findOverlaps(x, subpop_CNV_gr[[k]])
#    x <- x[queryHits(olap)]
#  }
#  return(x)
#})

################################################################################
### 6. Plot CNV signal with CNV genomic length annotations ###
################################################################################

#for (i in 1:length(subpop_matrices)) {
#
#  average_signal <- round(apply(subpop_matrices[[i]], 2, mean), 6)
#
#  # create area plot presenting InferCNV signal:
#  area_df <- data.frame(
#    index = seq_along(average_signal),
#    average_score = average_signal-neutral_value,
#    type = "neutral",
#    stringsAsFactors = F
#  )
#  
#  # label gains and losses:
#  for (r in 1:nrow(subpop_CNV_data[[i]])) {
#    area_df$type[
#      subpop_CNV_data[[i]]$start[r]:subpop_CNV_data[[i]]$end[r]
#    ] <- as.character(subpop_CNV_data[[i]]$call[r])
#  }
#  area_df$type[
#    area_df$type == "neutral" & area_df$average_score > 0
#  ] <- "gain_artefact"
#  area_df$type[
#    area_df$type == "neutral" & area_df$average_score < 0
#  ] <- "loss_artefact"
#  area_df$type <- factor(
#    area_df$type, levels = c(
#      "loss_artefact", "gain_artefact", "loss", "gain", "neutral"
#    )
#  )
#  # define colours:
#  cols <- c("#D6EAE8", "#F4E9E9", "#76C1C1", "#F7B7B5", "black")
#  
#  # prepare label data:
#  CNV_only <- subpop_CNV_data[[i]][subpop_CNV_data[[i]]$call != "neutral",]
#  CNV_only$length_labels <- CNV_only$genomic_length
#  CNV_only$midpoints <- CNV_only$midpoints + 10
#  for (r in 2:nrow(CNV_only)) {
#    if (CNV_only$start[r] - CNV_only$end[r-1] < 40) {
#      CNV_only$midpoints[r] <- CNV_only$midpoints[r]+30
#    }
#  }
##  CNV_only$length_labels[c(FALSE, TRUE, TRUE)] <- 
##    paste0("\n", CNV_only$length_labels[c(FALSE, TRUE, TRUE)])
##  CNV_only$length_labels[c(FALSE, FALSE, TRUE)] <- 
##    paste0("\n", CNV_only$length_labels[c(FALSE, FALSE, TRUE)])
#  
#  # plot on barplot:
#  p <- ggplot(area_df, aes(x=index, y=average_score, fill = type))
#  p <- p + geom_bar(stat="identity")
#  p <- p + scale_fill_manual(values = cols)
#  p <- p + scale_x_continuous(
#    label = CNV_only$length_labels,
#    breaks = CNV_only$midpoints,
#    expand = c(0,0)
#  )
#  p <- p + scale_y_continuous(
#    limits = c(-0.09, 0.09)
#  )
#  p <- p + theme_cowplot(12)
#  p <- p + theme(
#    panel.grid.major = element_blank(),
#    panel.grid.minor = element_blank(),
#    legend.position = "none",
#    text = element_text(size=20),
#    #axis.text.y = element_text(size=20),
#    axis.title.y = element_blank(),
#    axis.text.y = element_blank(),
#    axis.ticks.y = element_blank(),
#    axis.text.x = element_text(size=12, angle=55, hjust=1),
#    axis.title.x = element_blank(),
#    #axis.text.x = element_blank()
#    axis.ticks.x = element_blank()
#  )
#  #p <- p + ylab("Mean CNV signal")
#  #p <- p + xlab("Length (kb)")
#  for (c_end in chr_data$ends) {
#    p <- p + geom_vline(xintercept=c_end)
#  }
#  # create 0 line:
#  p <- p + geom_segment(
#    x=subpop_CNV_data[[i]]$start[1],
#    xend=subpop_CNV_data[[i]]$end[
#      nrow(subpop_CNV_data[[i]])
#    ],
#    y=0,
#    yend=0
#  )
#
#  # convert barplot to grid object and plot:
#  if (i==1) {
#    signal_plots <- list(ggplotGrob(p))
#    dev.off()
#  } else {
#    signal_plots[[i]] <- ggplotGrob(p)
#    dev.off()
#  }
#
#}
#names(signal_plots) <- names(subpop_matrices)
#
#pdf(
#  paste0(plot_dir, "CNV_signal_plots.pdf"), 
#  height = 15, 
#  width = 20
#)   
#  # draw signal plots and titles:
#  grid.newpage()
#
#  for (p in 1:length())
#
#  # plot 1:
#  pushViewport(viewport(x = 0.07, y = 0.87, width = 0.88, height = 0.12, 
#    just = c("left", "bottom")))
#    grid.draw(signal_plots[[1]])
#  popViewport()
#  # title 1:
#  pushViewport(viewport(x = unit(0.1, "cm"), y = unit(0.8, "npc"), width = 0.1, height = 0.1, 
#    just = c("left", "bottom")))
#    grid.text(names(signal_plots)[1], gp=gpar(fontsize=18))
#  popViewport()
# 
#  # plot 2:
#  pushViewport(viewport(x = 0.07, y = 0.75, width = 0.88, height = 0.12, 
#    just = c("left", "bottom")))
#    grid.draw(signal_plots[[2]])
#  popViewport()
#  # title 2:
#  pushViewport(viewport(x = unit(0.1, "cm"), y = unit(0.7, "npc"), width = 0.1, height = 0.1, 
#    just = c("left", "bottom")))
#    grid.text(names(signal_plots)[2], gp=gpar(fontsize=18))
#  popViewport()
#
#  # plot 3:
#  pushViewport(viewport(x = 0.07, y = 0.63, width = 0.88, height = 0.12, 
#    just = c("left", "bottom")))
#    grid.draw(signal_plots[[3]])
#  popViewport()
#  # title 3:
#  pushViewport(viewport(x = unit(0.1, "cm"), y = unit(0.6, "npc"), width = 0.1, height = 0.1, 
#    just = c("left", "bottom")))
#    grid.text(names(signal_plots)[3], gp=gpar(fontsize=18))
#  popViewport()
# 
#  # plot 4:
#  pushViewport(viewport(x = 0.07, y = 0.51, width = 0.88, height = 0.12, 
#    just = c("left", "bottom")))
#    grid.draw(signal_plots[[4]])
#  popViewport()
#  # title 4:
#  pushViewport(viewport(x = unit(0.1, "cm"), y = unit(0.5, "npc"), width = 0.1, height = 0.1, 
#    just = c("left", "bottom")))
#    grid.text(names(signal_plots)[4], gp=gpar(fontsize=18))
#  popViewport()
#  
#  # plot 5:
#  pushViewport(viewport(x = 0.07, y = 0.39, width = 0.88, height = 0.12, 
#    just = c("left", "bottom")))
#    grid.draw(signal_plots[[5]])
#  popViewport()
#  # title 5:
#  pushViewport(viewport(x = unit(0.1, "cm"), y = unit(0.4, "npc"), width = 0.1, height = 0.1, 
#    just = c("left", "bottom")))
#    grid.text(names(signal_plots)[5], gp=gpar(fontsize=18))
#  popViewport()
# 
#  # plot 6:
#  pushViewport(viewport(x = 0.07, y = 0.27, width = 0.88, height = 0.12, 
#    just = c("left", "bottom")))
#    grid.draw(signal_plots[[6]])
#  popViewport()
#  # title 6:
#  pushViewport(viewport(x = unit(0.1, "cm"), y = unit(0.3, "npc"), width = 0.1, height = 0.1, 
#    just = c("left", "bottom")))
#    grid.text(names(signal_plots)[6], gp=gpar(fontsize=18))
#  popViewport()
# 
#  # plot 7:
#  pushViewport(viewport(x = 0.07, y = 0.15, width = 0.88, height = 0.12, 
#    just = c("left", "bottom")))
#    grid.draw(signal_plots[[7]])
#  popViewport()
#  # title 7:
#  pushViewport(viewport(x = unit(0.1, "cm"), y = unit(0.2, "npc"), width = 0.1, height = 0.1, 
#    just = c("left", "bottom")))
#    grid.text(names(signal_plots)[7], gp=gpar(fontsize=18))
#  popViewport()
#  
#  # plot 8:
#  pushViewport(viewport(x = 0.07, y = 0.03, width = 0.88, height = 0.12, 
#    just = c("left", "bottom")))
#    grid.draw(signal_plots[[8]])
#  popViewport()
#  # title 8:
#  pushViewport(viewport(x = unit(0.1, "cm"), y = unit(0.1, "npc"), width = 0.1, height = 0.1, 
#    just = c("left", "bottom")))
#    grid.text(names(signal_plots)[8], gp=gpar(fontsize=18))
#  popViewport()
#  
#  # draw chromosome labels:
#  for ( e in 1:length(chr_data$lab_pos) ) {
#    pushViewport(viewport(x = 0.05 + chr_data$lab_pos[e]/1.155, y = 0.97, width = 0.05, height = 0.05, 
#      just = c("left", "bottom")))
#      if (e==1) {
#        grid.text(names(chr_data$lab_pos)[e], gp=gpar(fontsize=15))
#      } else {
#        grid.text(gsub("chr", "", names(chr_data$lab_pos)[e]), gp=gpar(fontsize=15))
#      }
#    popViewport()
#  }
#  # draw x-axis label:
#  pushViewport(viewport(x = 0.45, y = 0.001, width = 0.05, height = 0.05, 
#    just = c("left", "bottom")))
#    grid.text("CNV length (kb)", gp=gpar(fontsize=20))
#  popViewport()
#
#dev.off()
#
#