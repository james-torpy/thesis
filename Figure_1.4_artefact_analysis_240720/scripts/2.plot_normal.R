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
remove_cell_types <- args[3]
print("Cell types removed:")
print(remove_cell_types)
min_artefact_proportion <- as.numeric(args[4])
print("Minimum artefact proportion = ")
print(min_artefact_proportion)
min_artefact_length <- as.numeric(args[5])
print("Minimum artefact length = ")
print(min_artefact_length)

project_name <- "thesis"
subproject_name <- "Figure_1.4_artefact_analysis"
sample_name <- "CID4520N"
print(paste0("sample name = ", sample_name))
analysis_mode <- "samples"
print(paste0("Analysis mode of normal InferCNV run = ", 
  analysis_mode, "_mode"))
remove_cell_types <- "No"
print("Cell types removed:")
min_artefact_proportion <- 0.5
print("Minimum artefact proportion = ")
print(min_artefact_proportion)
min_artefact_length <- 20
print("Minimum artefact length = ")
print(min_artefact_length)

print(paste0("Project name = ", project_name))
print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38, lib.loc = lib_loc)
library(cowplot)
library(rlist, lib.loc = lib_loc)
library(GenomicRanges)
library(naturalsort, lib.loc = lib_loc)
library(splatter, lib.loc = lib_loc)
library(ggplot2)
library(org.Hs.eg.db)
library(data.table)
library(ComplexHeatmap, lib.loc = lib_loc)
library(circlize, lib.loc = lib_loc)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", 
  project_name, "/", subproject_name, "/")
ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/functions/")
results_dir <- paste0(project_dir, "results/")

in_dir <- paste0(results_dir, "infercnv/", sample_name, "/", 
  remove_cell_types, "_removed/", analysis_mode, "_mode/")
input_dir <- paste0(results_dir, "infercnv/", sample_name, 
  "/input_files/")

Robject_dir <- paste0(in_dir, "Rdata/")
system(paste0("mkdir -p ", Robject_dir))
plot_dir <- paste0(in_dir, "plots/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(in_dir, "tables/")
system(paste0("mkdir -p ", table_dir))

print(paste0("Sample directory = ", in_dir))
print(paste0("Robject directory = ", Robject_dir))
print(paste0("Plot dir = ", plot_dir))
print(paste0("Table directory = ", table_dir))

print(paste0("Plotting InferCNV heatmap of ", sample_name))


################################################################################
### 0. Define functions and choose seeds ###
################################################################################

fetch_chromosome_boundaries <- dget(paste0(func_dir, 
  "fetch_chromosome_boundaries.R"))

create_extended_vector <- function(df, column) {
  for (j in 1:nrow(df)) {
    vector_length <- length(seq(df$start[j], df$end[j]))
    if (j==1) {
      result_vector <- as.character(
        rep(
          eval(parse(text=paste0("df$", column, "[j]"))), vector_length
        )
      )
    } else {
      result_vector <- c(
        result_vector,
        as.character(
          rep(
            eval(parse(text=paste0("df$", column, "[j]"))), vector_length
          )
        )
      )
    }
  }
  return(result_vector)
}

# accuracy annotation legend function:
artefact_legend <- function() {
  # plot legend text:
  pushViewport(viewport(x = 0.5, y = 0.655, 
                        width = unit(1, "cm"), height = unit(0.5, "cm"), 
                        just = c("left")))
    grid.text("gain artefact", gp=gpar(fontsize=18))
  popViewport()
  pushViewport(viewport(x = 0.5, y = 0.435, 
                        width = unit(1, "cm"), height = unit(0.5, "cm"), 
                        just = c("left")))
    grid.text("loss artefact", gp=gpar(fontsize=18))
  popViewport()
  
  # plot legend squares:
  pushViewport(viewport(x = 0, y = 0.595, 
                        width = unit(1, "cm"), height = unit(0.5, "cm"), 
                        just = c("left")))
    grid.rect(width = unit(5, "mm"), height = unit(5, "mm"),
              just = c("left", "bottom"), gp=gpar(col = "#BF3667", fill = "#BF3667"))
  popViewport()
  pushViewport(viewport(x = 0, y = 0.38, 
                        width = unit(1, "cm"), height = unit(0.5, "cm"), 
                        just = c("left")))
    grid.rect(width = unit(5, "mm"), height = unit(5, "mm"),
            just = c("left", "bottom"), gp=gpar(col = "#58B9DB", fill = "#58B9DB"))
  popViewport()

}


###################################################################################
### 1. Load InferCNV output and screen for artefacts ###
###################################################################################

if (
  !file.exists(
    paste0(Robject_dir, "1a.artefact_indices.Rdata")
  ) | 
  !file.exists(
    paste0(Robject_dir, "1b.artefact_by_gene.Rdata")
  ) |
  !file.exists(
    paste0(Robject_dir, "1c.artefact_annotation.Rdata")
  )
) {

  # load normal infercnv output:
  print("Loading normal InferCNV output files...")
  epithelial_heatmap <- as.data.frame(t(read.table(paste0(in_dir, 
    "infercnv.12_denoised.observations.txt"))))

  # determine neutral value (based on sd denoising):
  score_table <- table(round(unlist(epithelial_heatmap), 6))
  neutral_value <- as.numeric(
    names(score_table)[which.max(score_table)]
  )

  # determine whether there is a gain or loss signal for each gene:
  artefact_by_gene <- apply(epithelial_heatmap, 2, function(x) {
    if (length(which(round(x, 6) > neutral_value)) >= 
      length(x)*min_artefact_proportion) {
      return("gain")
    } else if (length(which(round(x, 6) < neutral_value)) >= 
      length(x)*min_artefact_proportion) {
      return("loss")
    } else {
      return("neutral")
    }
  })

  # remove gain or loss artefacts in <10 gene stretches:
  run_length <- rle(artefact_by_gene)
  rle_df <- data.frame(
    row.names = c(
      colnames(epithelial_heatmap)[1],
      names(run_length$lengths)[1:(length(run_length$lengths)-1)]
    ),
    length = run_length$lengths,
    value = run_length$values
  )
  temp_artefact_df <- rle_df[rle_df$value != "neutral",]
  non_artefact <- temp_artefact_df[temp_artefact_df$length < min_artefact_length,]

  if (nrow(non_artefact) > 0) {
    for (r in 1:nrow(non_artefact)) {
      gene_ind <- which(names(artefact_by_gene) == rownames(non_artefact)[r])
      artefact_by_gene[gene_ind:(gene_ind + non_artefact$length[r] - 1)] <- 
        "neutral"
    }
  }

  # record artefact indices:
  artefact_df <- data.frame(
    index = 1:length(artefact_by_gene),
    type = artefact_by_gene
  )
  # split by type:
  split_artefact <- split(artefact_df, rleid(artefact_df$type))
  # convert into indices:
  split_indices <- lapply(split_artefact, function(x) {
    return(
      data.frame(
        start = x$index[1],
        end = x$index[nrow(x)],
        type = as.character(x$type[1])
      )
    )
  })
  artefact_indices <- do.call("rbind", split_indices)

  # fetch chromosome info:
  chr_data <- fetch_chromosome_boundaries(epithelial_heatmap, ref_dir)

  if (nrow(artefact_indices) > 0) {
    for (m in 1:nrow(artefact_indices)) {
      if (m==1) {
        heatmap_artefact <- artefact_indices$start[m]:artefact_indices$end[m]
      } else {
        heatmap_artefact <- c(heatmap_artefact,
          artefact_indices$start[m]:artefact_indices$end[m]
        )
      }   
    }  

    for (k in 1:length(chr_data$ends)) {
      if (k==1) {
  
        artefact_indices$start_chr[
          artefact_indices$start <= chr_data$ends[k]
        ] <- names(chr_data$ends)[k]
  
        artefact_indices$end_chr[
          artefact_indices$end <= chr_data$ends[k]
        ] <- names(chr_data$ends)[k]
  
      } else {
  
        artefact_indices$start_chr[
          artefact_indices$start <= chr_data$ends[k] & 
          artefact_indices$start > chr_data$ends[k-1]
        ] <- names(chr_data$ends)[k]
  
        artefact_indices$end_chr[
          artefact_indices$end <= chr_data$ends[k] & 
          artefact_indices$end > chr_data$ends[k-1]
        ] <- names(chr_data$ends)[k]
  
      }
    }
  
    # create artefact annotation:
    artefact_annotation_vector <- factor(
      artefact_by_gene,
      levels = c("neutral", "gain", "loss")
    )
    artefact_cols <- structure(
      c("#E5E4DF", "#BF3667", "#58B9DB"),
      names = levels(artefact_annotation_vector)
    )
    artefact_heatmap <- Heatmap(
      t(matrix(artefact_annotation_vector)),
      name = "artefact_heatmap",
      col = artefact_cols,
      show_heatmap_legend = FALSE
    )
    artefact_heatmap@name <- "artefact_heatmap"
    artefact_heatmap_obj <- grid.grabExpr(
      draw(artefact_heatmap, heatmap_legend_side = "left")
    )
    dev.off()

    saveRDS(
      artefact_indices, 
      paste0(Robject_dir, "1a.artefact_indices.Rdata")
    )
    saveRDS(
      artefact_by_gene, 
      paste0(Robject_dir, "1b.artefact_by_gene.Rdata")
    )
    saveRDS(
      artefact_heatmap_obj, 
      paste0(Robject_dir, "1c.artefact_annotation.Rdata")
    )

  }

} else {

  artefact_indices <- readRDS(
    paste0(Robject_dir, "1a.artefact_indices.Rdata")
  )
  artefact_by_gene <- readRDS(
    paste0(Robject_dir, "1b.artefact_by_gene.Rdata")
  )
  artefact_heatmap_obj <- readRDS(
    paste0(Robject_dir, "1c.artefact_annotation.Rdata")
  )

  # load normal infercnv output:
  print("Loading normal InferCNV output files...")
  epithelial_heatmap <- as.data.frame(t(read.table(paste0(in_dir, 
    "infercnv.12_denoised.observations.txt"))))

  # fetch chromosome info:
  chr_data <- fetch_chromosome_boundaries(epithelial_heatmap, ref_dir)

}


###################################################################################
### 2. Convert gene lengths to genomic lengths ###
###################################################################################

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

# identify start and end genes of each artefact:
artefact_indices$start_gene <- colnames(epithelial_heatmap)[
  artefact_indices$start
]
artefact_indices$end_gene <- colnames(epithelial_heatmap)[
  artefact_indices$end
]
# calculate genomic start of CNVs:
merge_coords <- subset(gene_coords, select = c(gene_id, start))
colnames(merge_coords) <- c("start_gene", "genomic_start")
artefact_indices <- merge(
  artefact_indices,
  merge_coords,
  by = "start_gene"
)
# calculate genomic start of CNVs:
merge_coords <- subset(gene_coords, select = c(gene_id, end))
colnames(merge_coords) <- c("end_gene", "genomic_end")
artefact_indices <- merge(
  artefact_indices,
  merge_coords,
  by = "end_gene"
)
# calculate genomic length of CNVs in mb:
artefact_indices$genomic_length <- 
  round((artefact_indices$genomic_end-artefact_indices$genomic_start)/1000000, 1)
# calculate midpoints for labelling:
artefact_indices$midpoints <- 
  artefact_indices$start + floor((artefact_indices$end-artefact_indices$start)/2)

# isolate artefacts only:
artefact_only <- artefact_indices[artefact_indices$type != "neutral",]
artefact_only <- artefact_only[order(artefact_only$genomic_start),]

# insert newlines where needed for labelling:
newline_record <- data.frame(
  row.names = 1:3,
  in_label = rep(FALSE, 3)
)
artefact_only$length_lab <- artefact_only$genomic_length
for (i in 2:nrow(artefact_only)) {
  if ( (artefact_only$start[i] - artefact_only$end[i-1]) < 200 ) {
  	if ( length(grep("\n", artefact_only$length_lab[i-1])) != 0 ) {

  	  if ( length(strsplit(artefact_only$length_lab[i-1], "\n")[[1]]) == 3 ) {

  	    artefact_only$length_lab[i] <- paste0("\n\n\n", artefact_only$length_lab[i])
  	    newline_record["3",1] <- TRUE

   	  } else if ( length(strsplit(artefact_only$length_lab[i-1], "\n")[[1]]) == 2 ) {
 
   	    artefact_only$length_lab[i] <- paste0("\n\n", artefact_only$length_lab[i])
   	    newline_record["2",1] <- TRUE
 
   	  } 

  	} else {
  	  artefact_only$length_lab[i] <- paste0("\n", artefact_only$length_lab[i])
	  newline_record["1",1] <- TRUE
  	}
  }
}


###################################################################################
### 3. Plot CNV heatmap ###
###################################################################################

# prepare df for plotting:
plot_object <- epithelial_heatmap
colnames(plot_object) <- rep("la", ncol(plot_object))
plot_object <- as.matrix(plot_object)

# define heatmap colours:
na_less_vector <- unlist(plot_object)
na_less_vector <- na_less_vector[!is.na(na_less_vector)]
heatmap_cols <- colorRamp2(c(min(na_less_vector), 1, max(na_less_vector)), 
      c("#00106B", "white", "#680700"), space = "sRGB")
print("Generating filtered normal heatmap...")

# generate heatmap legend:
rounding_no <- 1
signal_ranges <- round(range(unlist(plot_object)), rounding_no)
while (signal_ranges[1] == signal_ranges[2]) {
  rounding_no <- rounding_no + 1
  signal_ranges <- round(range(unlist(plot_object)), rounding_no)
}

lgd <- Legend(
  at = c(signal_ranges[1], 1, signal_ranges[2]),
  labels = c("loss", "", "gain"),
  col_fun = heatmap_cols, 
  title = "CNV signal", 
  direction = "horizontal",
  grid_height = unit(3, "cm"),
  grid_width = unit(4.5, "cm"),
  legend_height = unit(3, "cm"),
  legend_width = unit(4.5, "cm"),
  labels_gp = gpar(fontsize = 30),
  title_gp = gpar(fontsize = 32, fontface = "plain")
)

final_heatmap_basic <- Heatmap(
  plot_object, name = paste0("hm"), 
  col = heatmap_cols,
  cluster_columns = F, cluster_rows = F,
  show_row_names = F, show_column_names = F,
  #column_names_gp = gpar(col = "white"),
  show_row_dend = F,
  show_heatmap_legend = F,
  use_raster = T, raster_device = c("png")
)

annotated_heatmap <- grid.grabExpr(
  draw(final_heatmap_basic, gap = unit(6, "mm"), heatmap_legend_side = "left")
)
dev.off()

png(
  paste0(plot_dir, "infercnv_plot.png"), 
  height = 13, width = 23, res = 300, units = "in"
)   
  pushViewport(viewport(x = 0.147, y = 0.15, width = 0.805, height = 0.8, 
      just = c("left", "bottom")))
    grid.draw(annotated_heatmap)
    decorate_heatmap_body("hm", {
      for ( e in 1:length(chr_data$end_pos) ) {
        grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), 
          gp = gpar(lwd = 3, col = "#383838"))
        if (e==1) {
          grid.text(names(chr_data$lab_pos)[e], chr_data$lab_pos[e], 
          unit(0, "npc") + unit(-3.5, "mm"), gp=gpar(fontsize=26))
        } else if (e==21) {
          grid.text(paste0("\n", gsub("chr", "", names(chr_data$lab_pos)[e])), chr_data$lab_pos[e], 
          unit(0, "npc") + unit(-3.5, "mm"), gp=gpar(fontsize=26))
        } else {
          grid.text(gsub("chr", "", names(chr_data$lab_pos)[e]), chr_data$lab_pos[e], 
          unit(0, "npc") + unit(-3.5, "mm"), gp=gpar(fontsize=26))
        }
        
      }
    })
    #grid.rect()
  popViewport()

  # plot legend:
  pushViewport(viewport(x = unit(2, "cm"), y = unit(16, "cm"), width = unit(5, "cm"), 
    height = unit(7, "cm"), just = c("left", "bottom")))
    #grid.rect()
    draw(lgd)
  popViewport()

dev.off()


#################################################################################
#### 4. Plot average signal vs artefact annotation ###
#################################################################################

# determine neutral value (based on sd denoising):
if (!exists("neutral_value")) {
  print("Loading normal InferCNV output files...")
  epithelial_heatmap <- as.data.frame(t(read.table(paste0(in_dir, 
    "infercnv.12_denoised.observations.txt"))))
  score_table <- table(round(unlist(epithelial_heatmap), 6))
  neutral_value <- as.numeric(
    names(score_table)[which.max(score_table)]
  )
}

# create area plot presenting InferCNV signal:
average_epithelial <- apply(epithelial_heatmap, 2, mean)

area_df <- data.frame(
  index = seq_along(average_epithelial),
  average_score = average_epithelial-neutral_value,
  type = "neutral",
  stringsAsFactors = F
)

# label gains and losses:
area_df$type[area_df$average_score > 0] <- "gain"
area_df$type[area_df$average_score < 0] <- "loss"

# define colours:
cols <- c("#F7B7B5", "#76C1C1", "black")
area_df$type <- factor(area_df$type, levels = c("gain", "loss", "neutral"))

# plot on barplot:
p <- ggplot(area_df, aes(x=index, y=average_score, fill = type))
p <- p + geom_bar(stat="identity")
p <- p + scale_fill_manual(values = cols)
p <- p + scale_x_continuous(
  limits = c(
    0,ncol(epithelial_heatmap)
  ), 
  expand = c(0, 0),
  breaks = artefact_only$midpoints,
  labels = c(artefact_only$length_lab)
)
p <- p + scale_y_continuous(
  limits = c(-0.09, 0.09)
)
p <- p + theme_cowplot(12)
p <- p + theme(
  axis.title.x = element_text(size=30, margin = margin(t = 30, r = 0, b = 0, l = 0)),
  axis.text.x = element_text(size=23, margin = margin(t = 70, r = 0, b = 0, l = 0)),
  axis.ticks.x = element_blank(),
  axis.title.y = element_text(size=30, margin = margin(t = 0, r = 30, b = 0, l = 0)),
  axis.text.y = element_text(size=24),
  legend.position = "none"
)
p <- p + xlab("Artefact genomic length (mb)")
p <- p + ylab("Mean CNV signal")
for (c_end in chr_data$ends) {
  p <- p + geom_vline(xintercept=c_end)
}

# convert barplot to grid object:
signal_plot <- ggplotGrob(p)
dev.off()

# plot signal and accuracy annotation together:
png(
  paste0(plot_dir, "signal_vs_simulated_CNV_plot.png"), 
  height = 8, 
  width = 22,
  res = 300,
  units = "in"
)   
  grid.newpage()

  	# plot signal:
    pushViewport(viewport(x = 0.52, y = 0.5, width = 0.93, height = 0.9))
      grid.draw(signal_plot)
    popViewport()

    # plot artefact annotation:
    if (newline_record["3",]) {
      pushViewport(viewport(x = 0.556, y = 0.365, width = 0.857, height = 0.13))
    } else if (newline_record["2",]) {
      pushViewport(viewport(x = 0.556, y = 0.3, width = 0.857, height = 0.13))
    } else if (newline_record["1",]) {
      pushViewport(viewport(x = 0.556, y = 0.3, width = 0.857, height = 0.13))
    } else {
      pushViewport(viewport(x = 0.556, y = 0.3, width = 0.857, height = 0.13))
    }
      grid.draw(artefact_heatmap_obj)
    popViewport()
    
    pushViewport(viewport(x = 0.555, y = 0.925, width = 0.85, height = 0.07))
      #grid.rect()
      for ( e in 1:length(chr_data$lab_pos) ) {
        if (e==1) {
          pushViewport(viewport(x = chr_data$lab_pos[e], y = 0.5, width = 0.01, height = 0.5))
            grid.text(names(chr_data$lab_pos)[e], gp=gpar(fontsize=22))
          popViewport()
        } else if (e==21) {
          pushViewport(viewport(x = chr_data$lab_pos[e], y = 0.9, width = 0.01, height = 0.5))
            grid.text(gsub("chr", "", names(chr_data$lab_pos)[e]), gp=gpar(fontsize=22))
          popViewport()
        } else {
          pushViewport(viewport(x = chr_data$lab_pos[e], y = 0.5, width = 0.01, height = 0.5))
            grid.text(gsub("chr", "", names(chr_data$lab_pos)[e]), gp=gpar(fontsize=22))
          popViewport()
        }
      }
    popViewport()
    pushViewport(viewport(x = unit(1, "cm"), y = unit(3.5, "cm"), 
                      width = unit(5.5, "cm"), height = unit(4.5, "cm"), 
                      just = c("left", "bottom")))
      #grid.rect()
      artefact_legend()
    popViewport()
dev.off()


# create chromosome border plot:
p <- ggplot(area_df, aes(x=index, y=average_score, fill = type))
p <- p + scale_x_continuous(
  limits = c(
    0,ncol(epithelial_heatmap)
  ), 
  expand = c(0, 0),
)
p <- p + scale_y_continuous(
  limits = c(-0.09, 0.09)
)
p <- p + theme_cowplot(12)
p <- p + theme(
  axis.line=element_blank(),
  axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  legend.position = "none"
)
for (c_end in chr_data$ends) {
  p <- p + geom_vline(xintercept=c_end)
}

# convert barplot to grid object:
chr_border_plot <- ggplotGrob(p)

plot_data <- list(
  signal_plot = signal_plot,
  artefact_heatmap_obj = artefact_heatmap_obj,
  chr_data = chr_data,
  chr_border_plot = chr_border_plot
)
saveRDS(plot_data, paste0(Robject_dir, "signal_plot_data.Rdata"))


#################################################################################
#### 4. Record first and last 10 genes of each artefact ###
#################################################################################

for (i in 1:nrow(artefact_only)) {

  df <- data.frame(
    first_ten = names(
      artefact_genes[
        (artefact_only$start[i]):(artefact_only$start[i]+9)
      ]
    ),
    last_ten = names(
      artefact_genes[
        (artefact_only$end[i]):(artefact_only$end[i]-9)
      ]
    )
  )

  if (i==1) {
    artefact_gene_record <- list(df)
  } else {
    artefact_gene_record[[i]] <- df
  }
  names(artefact_gene_record)[i] <- artefact_only$start_chr[i]

}

######

## fetch all genes in artefacts:
#all_artefact_genes <- lapply(artefact_gene_record, function(x) {
#  gene_coords$gene_id[
#    which(gene_coords$gene == x$first_ten[1]):which(gene_coords$gene == x$last_ten[10])
#  ]
#})
#
#saveRDS(
#  all_artefact_genes,
#  paste0(Robject_dir, "artefact_gene_record.Rdata")
#)


#######
#
## mock plot:
## plot on barplot:
#p <- ggplot(area_df, aes(x=index, y=average_score, fill = type))
#p <- p + scale_x_continuous(
#  limits = c(
#    0,ncol(epithelial_heatmap)
#  ), 
#  expand = c(0, 0),
#  breaks = artefact_only$midpoints,
#  labels = artefact_only$length_lab
#)
#p <- p + scale_y_continuous(
#  limits = c(-0.09, 0.09)
#)
#p <- p + theme_cowplot(12)
#p <- p + theme(
#  axis.title.x = element_text(size=30, margin = margin(t = 30, r = 0, b = 0, l = 0)),
#  axis.text.x = element_text(size=22, margin = margin(t = 70, r = 0, b = 0, l = 0)),
#  axis.ticks.x = element_blank(),
#  axis.title.y = element_text(size=30, margin = margin(t = 0, r = 30, b = 0, l = 0)),
#  axis.text.y = element_text(size=24),
#  legend.position = "none"
#)
#p <- p + xlab("Artefact genomic length (mb)")
#p <- p + ylab("Mean CNV signal")
#for (c_end in chr_data$ends) {
#  p <- p + geom_vline(xintercept=c_end)
#}
#
## convert barplot to grid object:
#mock_signal_plot <- ggplotGrob(p)
#dev.off()
#
#png(
#  paste0(plot_dir, "mock_signal_vs_simulated_CNV_plot.png"), 
#  height = 8, 
#  width = 22,
#  res = 300,
#  units = "in"
#)   
#
#  grid.newpage()
#    pushViewport(viewport(x = 0.52, y = 0.5, width = 0.93, height = 0.9))
#      grid.draw(mock_signal_plot)
#    popViewport()
#    pushViewport(viewport(x = 0.556, y = 0.3, width = 0.857, height = 0.13))
#      grid.draw(artefact_heatmap_obj)
#    popViewport()
#    pushViewport(viewport(x = 0.555, y = 0.925, width = 0.85, height = 0.07))
#      #grid.rect()
#      for ( e in 1:length(chr_data$lab_pos) ) {
#        if (e==1) {
#          pushViewport(viewport(x = chr_data$lab_pos[e], y = 0.5, width = 0.01, height = 0.5))
#            grid.text(names(chr_data$lab_pos)[e], gp=gpar(fontsize=16, fontface = "bold"))
#          popViewport()
#        } else if (e==21) {
#          pushViewport(viewport(x = chr_data$lab_pos[e], y = 0.9, width = 0.01, height = 0.5))
#            grid.text(gsub("chr", "", names(chr_data$lab_pos)[e]), gp=gpar(fontsize=16, fontface = "bold"))
#          popViewport()
#        } else {
#          pushViewport(viewport(x = chr_data$lab_pos[e], y = 0.5, width = 0.01, height = 0.5))
#            grid.text(gsub("chr", "", names(chr_data$lab_pos)[e]), gp=gpar(fontsize=16, fontface = "bold"))
#          popViewport()
#        }
#      }
#    popViewport()
#    pushViewport(viewport(x = unit(1, "cm"), y = unit(0.8, "cm"), 
#                      width = unit(5.5, "cm"), height = unit(4.5, "cm"), 
#                      just = c("left", "bottom")))
#      #grid.rect()
#      artefact_legend()
#    popViewport()
#dev.off()


#png(
#  paste0(plot_dir, "signal_plot.png"), 
#  height = 8, 
#  width = 20,
#  res = 300,
#  units = "in"
#)
#  grid.newpage()
#    pushViewport(viewport(x = 0.027, y = 0.001, width = 0.964, height = 0.9, 
#      just = c("left", "bottom")))
#      grid.draw(signal_plot)
#    popViewport()
#dev.off()



