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

sample_name <- "CID4515"
min_CNV_length <- 20
min_CNV_proportion <- as.numeric("0.5")

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

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/",
  subproject_name, "/")
ref_dir <- paste0(project_dir, "/refs/")
func_dir <- paste0(project_dir, "/scripts/functions/")
results_dir <- paste0(project_dir, "results/")
in_path <- paste0(results_dir, "infercnv/", sample_name, "/")

Robject_dir <- paste0(in_path, "Rdata/")
system(paste0("mkdir -p ", Robject_dir))
plot_dir <- paste0(in_path, "plots/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(in_path, "tables/")
system(paste0("mkdir -p ", table_dir))

print(paste0("In directory = ", in_dir))
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
  "/2a.epithelial_heatmap_with_CNV_subclusters.Rdata"))
epithelial_metadata <- readRDS(paste0(Robject_dir, 
  "/2b.epithelial_metadata_with_CNV_subclusters.Rdata"))
print(paste0(
  "Are epithelial_metadata rownames in the same order as epithelial_heatmap?? ",
  identical(rownames(epithelial_heatmap), rownames(epithelial_metadata))
))

# determine neutral CNV value as most frequent infercnv score and check:
score_table <- table(round(unlist(epithelial_heatmap), 6))
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

# make start and end genomic co-ordinates positions for entire genome,
# not just chromosomal indices:

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

# determine CNV co-ordinates and start and end genomic positions for each subpop:
subpop_CNV_data <- lapply(subpop_matrices, function(x) {

  # if at least 50% of cells have signal above or below neutral_value, 
  # label gene as loss/gain respectively:
  for (c in 1:ncol(x)) {
  
    # determine if loss:
    if ( length(which(
        round(x[,c], 6) < neutral_value
      )) >= min_CNV_proportion*length(x[,c]) ) {
  
      if (c==1) {
        call_vector <- c("loss")
      } else {
        call_vector[c] <- "loss"
      }
  
    } else if ( length(which(
        round(x[,c], 6) > neutral_value
      )) >= min_CNV_proportion*length(x[,c]) ) {
  
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
  
  # fetch chromosome boundary co-ordinates:
  if (!file.exists(paste0(Robject_dir, "chromosome_data.Rdata"))) {
    chr_data <- fetch_chromosome_boundaries(epithelial_heatmap, ref_dir)
    saveRDS(chr_data, paste0(Robject_dir, "chromosome_data.Rdata"))
  } else {
    chr_data <- readRDS(paste0(Robject_dir, "chromosome_data.Rdata"))
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
        start_chr = c(CNV_indices$start_chr[n], CNV_indices$end_chr[n]),
        end_chr = c(CNV_indices$start_chr[n], CNV_indices$end_chr[n]),
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
  CNV_indices$start_gene <- colnames(x)[CNV_indices$start]
  CNV_indices$end_gene <- colnames(x)[CNV_indices$end]
  
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

})

# create CNV only list:
subpop_CNV_only <- lapply(subpop_CNV_data, function(x) {
  CNV_only <- x[x$call != "neutral",]
  # remove CNVs less than min_CNV_length:
  return(CNV_only[CNV_only$length >= 20,])
})

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


################################################################################
### 3. Identify CNVs common across all subpopulations to fade them out ###
################################################################################

# convert all CNV data to genomic ranges objects:
subpop_CNV_gr <- lapply(subpop_CNV_only, function(x) {
  return(
    GRanges(
      seqnames = Rle(x$chr),
      ranges = IRanges(
        start = as.integer(x$start),
        end = as.numeric(x$end)
      ),
      strand = Rle("*"),
      call = x$call,
      genomic_start = x$genomic_start,
      genomic_end = x$genomic_end,
      genomic_length = x$genomic_length,
      midpoints = x$midpoints
    )
  )
})

common_CNV_gr <- lapply(subpop_CNV_gr, function(x) {
  for (k in 1:length(subpop_CNV_gr)) {
    # keep only ranges overlapping with all subpopulations:
    olap <- findOverlaps(x, subpop_CNV_gr[[k]])
    x <- x[queryHits(olap)]
  }
  return(x)
})


################################################################################
### 4. Plot CNV signal with CNV genomic length annotations ###
################################################################################

for (i in 1:length(subpop_matrices)) {

  average_signal <- round(apply(subpop_matrices[[i]], 2, mean), 6)

  # create area plot presenting InferCNV signal:
  area_df <- data.frame(
    index = seq_along(average_signal),
    average_score = average_signal-neutral_value,
    type = "neutral",
    stringsAsFactors = F
  )
  
  # label gains and losses:
  for (r in 1:nrow(subpop_CNV_data[[i]])) {
    area_df$type[
      subpop_CNV_data[[i]]$start[r]:subpop_CNV_data[[i]]$end[r]
    ] <- as.character(subpop_CNV_data[[i]]$call[r])
  }
  area_df$type[
    area_df$type == "neutral" & area_df$average_score > 0
  ] <- "gain_artefact"
  area_df$type[
    area_df$type == "neutral" & area_df$average_score < 0
  ] <- "loss_artefact"
  area_df$type <- factor(
    area_df$type, levels = c(
      "loss_artefact", "gain_artefact", "loss", "gain", "neutral"
    )
  )
  # define colours:
  cols <- c("#D6EAE8", "#F4E9E9", "#76C1C1", "#F7B7B5", "black")
  
  # prepare label data:
  CNV_only <- subpop_CNV_data[[i]][subpop_CNV_data[[i]]$call != "neutral",]
  CNV_only$length_labels <- CNV_only$genomic_length
  CNV_only$midpoints <- CNV_only$midpoints + 10
  for (r in 2:nrow(CNV_only)) {
    if (CNV_only$start[r] - CNV_only$end[r-1] < 40) {
      CNV_only$midpoints[r] <- CNV_only$midpoints[r]+30
    }
  }
#  CNV_only$length_labels[c(FALSE, TRUE, TRUE)] <- 
#    paste0("\n", CNV_only$length_labels[c(FALSE, TRUE, TRUE)])
#  CNV_only$length_labels[c(FALSE, FALSE, TRUE)] <- 
#    paste0("\n", CNV_only$length_labels[c(FALSE, FALSE, TRUE)])
  
  # plot on barplot:
  p <- ggplot(area_df, aes(x=index, y=average_score, fill = type))
  p <- p + geom_bar(stat="identity")
  p <- p + scale_fill_manual(values = cols)
  p <- p + scale_x_continuous(
    label = CNV_only$length_labels,
    breaks = CNV_only$midpoints,
    expand = c(0,0)
  )
  p <- p + scale_y_continuous(
    limits = c(-0.09, 0.09)
  )
  p <- p + theme_cowplot(12)
  p <- p + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    text = element_text(size=20),
    #axis.text.y = element_text(size=20),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size=12, angle=55, hjust=1),
    axis.title.x = element_blank(),
    #axis.text.x = element_blank()
    axis.ticks.x = element_blank()
  )
  #p <- p + ylab("Mean CNV signal")
  #p <- p + xlab("Length (kb)")
  for (c_end in chr_data$ends) {
    p <- p + geom_vline(xintercept=c_end)
  }
  # create 0 line:
  p <- p + geom_segment(
    x=subpop_CNV_data[[i]]$start[1],
    xend=subpop_CNV_data[[i]]$end[
      nrow(subpop_CNV_data[[i]])
    ],
    y=0,
    yend=0
  )

  ######

  pdf(
    paste0(plot_dir, "CNV_signal_plot_1.pdf"), 
      height = 15, 
      width = 20
    )   
    p
  dev.off()
  
  ######
  
  # convert barplot to grid object and plot:
  if (i==1) {
    signal_plots <- list(ggplotGrob(p))
    dev.off()
  } else {
    signal_plots[[i]] <- ggplotGrob(p)
    dev.off()
  }

}
names(signal_plots) <- names(subpop_matrices)

pdf(
  paste0(plot_dir, "CNV_signal_plots.pdf"), 
  height = 15, 
  width = 20
)   
  # draw signal plots and titles:
  grid.newpage()

  # plot 1:
  pushViewport(viewport(x = 0.07, y = 0.87, width = 0.88, height = 0.12, 
    just = c("left", "bottom")))
    grid.draw(signal_plots[[1]])
  popViewport()
  # title 1:
  pushViewport(viewport(x = unit(0.1, "cm"), y = unit(0.8, "npc"), width = 0.1, height = 0.1, 
    just = c("left", "bottom")))
    grid.text(names(signal_plots)[1], gp=gpar(fontsize=18))
  popViewport()
 
  # plot 2:
  pushViewport(viewport(x = 0.07, y = 0.75, width = 0.88, height = 0.12, 
    just = c("left", "bottom")))
    grid.draw(signal_plots[[2]])
  popViewport()
  # title 2:
  pushViewport(viewport(x = unit(0.1, "cm"), y = unit(0.7, "npc"), width = 0.1, height = 0.1, 
    just = c("left", "bottom")))
    grid.text(names(signal_plots)[2], gp=gpar(fontsize=18))
  popViewport()

  # plot 3:
  pushViewport(viewport(x = 0.07, y = 0.63, width = 0.88, height = 0.12, 
    just = c("left", "bottom")))
    grid.draw(signal_plots[[3]])
  popViewport()
  # title 3:
  pushViewport(viewport(x = unit(0.1, "cm"), y = unit(0.6, "npc"), width = 0.1, height = 0.1, 
    just = c("left", "bottom")))
    grid.text(names(signal_plots)[3], gp=gpar(fontsize=18))
  popViewport()
 
  # plot 4:
  pushViewport(viewport(x = 0.07, y = 0.51, width = 0.88, height = 0.12, 
    just = c("left", "bottom")))
    grid.draw(signal_plots[[4]])
  popViewport()
  # title 4:
  pushViewport(viewport(x = unit(0.1, "cm"), y = unit(0.5, "npc"), width = 0.1, height = 0.1, 
    just = c("left", "bottom")))
    grid.text(names(signal_plots)[4], gp=gpar(fontsize=18))
  popViewport()
  
  # plot 5:
  pushViewport(viewport(x = 0.07, y = 0.39, width = 0.88, height = 0.12, 
    just = c("left", "bottom")))
    grid.draw(signal_plots[[5]])
  popViewport()
  # title 5:
  pushViewport(viewport(x = unit(0.1, "cm"), y = unit(0.4, "npc"), width = 0.1, height = 0.1, 
    just = c("left", "bottom")))
    grid.text(names(signal_plots)[5], gp=gpar(fontsize=18))
  popViewport()
 
  # plot 6:
  pushViewport(viewport(x = 0.07, y = 0.27, width = 0.88, height = 0.12, 
    just = c("left", "bottom")))
    grid.draw(signal_plots[[6]])
  popViewport()
  # title 6:
  pushViewport(viewport(x = unit(0.1, "cm"), y = unit(0.3, "npc"), width = 0.1, height = 0.1, 
    just = c("left", "bottom")))
    grid.text(names(signal_plots)[6], gp=gpar(fontsize=18))
  popViewport()
 
  # plot 7:
  pushViewport(viewport(x = 0.07, y = 0.15, width = 0.88, height = 0.12, 
    just = c("left", "bottom")))
    grid.draw(signal_plots[[7]])
  popViewport()
  # title 7:
  pushViewport(viewport(x = unit(0.1, "cm"), y = unit(0.2, "npc"), width = 0.1, height = 0.1, 
    just = c("left", "bottom")))
    grid.text(names(signal_plots)[7], gp=gpar(fontsize=18))
  popViewport()
  
  # plot 8:
  pushViewport(viewport(x = 0.07, y = 0.03, width = 0.88, height = 0.12, 
    just = c("left", "bottom")))
    grid.draw(signal_plots[[8]])
  popViewport()
  # title 8:
  pushViewport(viewport(x = unit(0.1, "cm"), y = unit(0.1, "npc"), width = 0.1, height = 0.1, 
    just = c("left", "bottom")))
    grid.text(names(signal_plots)[8], gp=gpar(fontsize=18))
  popViewport()
  
  # draw chromosome labels:
  for ( e in 1:length(chr_data$lab_pos) ) {
    pushViewport(viewport(x = 0.05 + chr_data$lab_pos[e]/1.155, y = 0.97, width = 0.05, height = 0.05, 
      just = c("left", "bottom")))
      if (e==1) {
        grid.text(names(chr_data$lab_pos)[e], gp=gpar(fontsize=15))
      } else {
        grid.text(gsub("chr", "", names(chr_data$lab_pos)[e]), gp=gpar(fontsize=15))
      }
    popViewport()
  }
  # draw x-axis label:
  pushViewport(viewport(x = 0.45, y = 0.001, width = 0.05, height = 0.05, 
    just = c("left", "bottom")))
    grid.text("CNV length (kb)", gp=gpar(fontsize=20))
  popViewport()

dev.off()







################################################################################
### 4. Remove CNVs below min CNV length ###
################################################################################

# annotate artefact CNVs:

# create separate subpop matrix list without CNV artefacts:






######
#  pdf(
#    paste0(plot_dir, "signal_plot_grid_test.pdf"), 
#    height = 8, 
#    width = 20
#  )   
#  
#    grid.newpage()
#    pushViewport(viewport(x = 0.02, y = 0.28, width = 0.925, height = 0.65, 
#      just = c("left", "bottom")))
#      grid.draw(signal_plot)
#    popViewport()
#  dev.off()
  
#  # create CNV call annotation:
#  call_vector <- create_extended_vector(subpop_CNV_data[[i]], "call")
#  call_cols <- structure(
#    c("#EF0F0F", "#069ABC", "#E7E4D3"),
#    names = c("gain", "loss", "neutral")
#  )
#  call_heatmap <- Heatmap(
#    t(matrix(call_vector)),
#    name = "annotation_heatmap",
#    col = call_cols,
#    show_heatmap_legend = FALSE
#  )
#  call_heatmap@name <- "call_heatmap"
#  
#  call_heatmap_obj <- grid.grabExpr(
#    draw(call_heatmap, heatmap_legend_side = "left")
#  )
#  dev.off()
#  pdf(paste0(plot_dir, "call_heatmap.pdf"), height = 2, width = 20)
#  grid.newpage()
#    pushViewport(viewport(x = 0.07, y = 0.1, width = 0.92, height = 0.85, 
#      just = c("left", "bottom")))
#      grid.draw(call_heatmap_obj)
#    popViewport()
#  dev.off()

#  # add call annotation to barplot:
#  pdf(
#    paste0(plot_dir, "CNV_signal_plot.pdf"), 
#    height = 8, 
#    width = 20
#  )   
#  
#    grid.newpage()
#    pushViewport(viewport(x = 0.02, y = 0.28, width = 0.925, height = 0.65, 
#      just = c("left", "bottom")))
#      grid.draw(signal_plot)
#    popViewport()
#    pushViewport(viewport(x = 0.076, y = 0.26, width = 0.87, height = 0.1, 
#      just = c("left", "bottom")))
#      grid.draw(call_heatmap_obj)
#    popViewport()
#    for ( e in 1:length(chr_data$lab_pos) ) {
#      pushViewport(viewport(x = 0.05 + chr_data$lab_pos[e]/1.155, y = 0.9, width = 0.05, height = 0.05, 
#        just = c("left", "bottom")))
#        if (e==1) {
#          grid.text(names(chr_data$lab_pos)[e], gp=gpar(fontsize=15))
#        } else {
#          grid.text(gsub("chr", "", names(chr_data$lab_pos)[e]), gp=gpar(fontsize=15))
#        }
#      popViewport()
#    }
#    pushViewport(viewport(x = 0.45, y = 0.13, width = 0.05, height = 0.05, 
#      just = c("left", "bottom")))
#      grid.text("CNV length (kb)", gp=gpar(fontsize=20))
#    popViewport()
#
#    # draw legend text:
#    pushViewport(viewport(x = unit(0.1, "cm"), y = unit(5.7, "cm"), width = 0.1, height = 0.1, 
#      just = c("left", "bottom")))
#      grid.text("gain", gp=gpar(fontsize=18))
#    popViewport()
#  
#    pushViewport(viewport(x = unit(1.5, "mm"), y = unit(4.58, "cm"), width = 0.1, height = 0.1, 
#      just = c("left", "bottom")))
#      grid.text("loss", gp=gpar(fontsize=18))
#    popViewport()
# 
#    # draw legend squares:
#    pushViewport(viewport(x = unit(1.1, "cm"), y = unit(4.37, "cm"), width = unit(2, "cm"), height = unit(2, "cm"), 
#      just = c("right", "bottom")))
#      grid.rect(x = 1, y = 1, width = unit(5, "mm"), height = unit(5, "mm"),
#        just = c("left", "bottom"), gp=gpar(col = "#EF0F0F", fill = "#EF0F0F"))
#    popViewport()
#    pushViewport(viewport(x = unit(1.1, "cm"), y = unit(3.33, "cm"), width = unit(2, "cm"), height = unit(2, "cm"), 
#      just = c("right", "bottom")))
#      grid.rect(x = 1, y = 1, width = unit(5, "mm"), height = unit(5, "mm"),
#        just = c("left", "bottom"), gp=gpar(col = "#069ABC", fill = "#069ABC"))
#    popViewport()
#
#  dev.off()





