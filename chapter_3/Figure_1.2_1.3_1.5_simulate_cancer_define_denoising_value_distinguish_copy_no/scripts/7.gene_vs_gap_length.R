#! /share/ClusterShare/software/contrib/briglo/octR/src/R-3.6.0/builddir/bin/RRscript

project_name <- "thesis"
subproject_name <- "Figure_2.2_and_2.3_define_denoising_value"
args = commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
analysis_mode <- args[2]
sim_names <- args[3]

project_name <- "thesis"
subproject_name <- "Figure_2.2_and_2.3_define_denoising_value"
args = commandArgs(trailingOnly=TRUE)
sample_name <- "CID4520N"
analysis_mode <- "samples"
sim_names <- paste0("sim", 1:30, collapse = ".")

# split multi-element variable vectors:
sim_names <- unlist(
  strsplit(
    sim_names,
    split = "\\."
  )
)
sim_names <- sim_names[grep("normal", sim_names, invert = T)]

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("Analysis mode = ", analysis_mode))
print("Simulation names = ")
print(sim_names)

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(rlang, lib.loc = lib_loc)
library(dplyr, lib.loc = lib_loc)
library(ggpubr, lib.loc = lib_loc)
#library(ggplot2)
library(dplyr)
library(cowplot)
theme_set(theme_minimal())
library(scales)
library(BSgenome.Hsapiens.UCSC.hg38, lib.loc = lib_loc)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", 
  project_name, "/", subproject_name, "/")
results_dir <- seurat_path <- paste0(project_dir, "results/")
ref_dir <- paste0(project_dir, "refs/")
in_path <- paste0(results_dir, "infercnv/", sample_name, "/")

out_path <- paste0(in_path, "gene_vs_genomic_length/")
Robject_dir <- paste0(out_path, "Rdata/")
system(paste0("mkdir -p ", Robject_dir))
plot_dir <- paste0(out_path, "plots/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(out_path, "tables/")
system(paste0("mkdir -p ", table_dir))


#################################################################################
#### 1. Load genomic data ###
#################################################################################

if (!file.exists(paste0(ref_dir, "infercnv_genomic_coordinates.Rdata"))) {

  #load gene co-ordinates and only keep chromosomes 1:22:
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
  
  saveRDS(gene_coords, paste0(ref_dir, "infercnv_genomic_coordinates.Rdata"))

} else {
  gene_coords <- readRDS(paste0(ref_dir, "infercnv_genomic_coordinates.Rdata"))
}


#################################################################################
#### 2. Load CNV indices for each simulation ###
#################################################################################

for (i in 1:length(sim_names)) {

  if (i==1) {

    CNV_data <- list(
      readRDS(
        paste0(in_path, sim_names[i], "/no_denoising/", analysis_mode, 
          "_mode/Rdata/CNV_data.Rdata")
      )
    )

  } else {

    CNV_data[[i]] <- readRDS(
      paste0(in_path, sim_names[i], "/no_denoising/", analysis_mode, 
        "_mode/Rdata/CNV_data.Rdata")
    )

  }

}


###################################################################################
### 3. Convert gene lengths to genomic lengths ###
###################################################################################

genomic_indices <- lapply(CNV_data, function(x) {

  CNV_indices <- x$CNV_indices

  # identify start and end genes of each CNV:
  CNV_indices$start_gene <- x$genes[CNV_indices$start]
  CNV_indices$end_gene <- x$genes[CNV_indices$end]

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
  
  # isolate CNVs only:
  CNV_only <- CNV_indices[CNV_indices$type != "neutral",]
  CNV_only <- CNV_only[order(CNV_only$genomic_start),]
  
  # insert newlines where needed for labelling:
  newline_record <- data.frame(
    row.names = 1:3,
    in_label = rep(FALSE, 3)
  )
  CNV_only$length_lab <- CNV_only$genomic_length
  for (i in 2:nrow(CNV_only)) {
    if ( (CNV_only$start[i] - CNV_only$end[i-1]) < 200 ) {
    	if ( length(grep("\n", CNV_only$length_lab[i-1])) != 0 ) {
  
    	  if ( length(strsplit(CNV_only$length_lab[i-1], "\n")[[1]]) == 3 ) {
  
    	    CNV_only$length_lab[i] <- paste0("\n\n\n", CNV_only$length_lab[i])
    	    newline_record["3",1] <- TRUE
  
     	  } else if ( length(strsplit(CNV_only$length_lab[i-1], "\n")[[1]]) == 2 ) {
   
     	    CNV_only$length_lab[i] <- paste0("\n\n", CNV_only$length_lab[i])
     	    newline_record["2",1] <- TRUE
   
     	  } 
  
    	} else {
    	  CNV_only$length_lab[i] <- paste0("\n", CNV_only$length_lab[i])
  	  newline_record["1",1] <- TRUE
    	}
    }
  }

  return(CNV_only)
  
})
all_CNV_data <- do.call("rbind", genomic_indices)
length_data <- subset(all_CNV_data, select = c(length, genomic_length))
  
  
###################################################################################
### 4. Plot length in genes vs in mb ###
###################################################################################

p <- ggplot(data = length_data, aes(x=length, y=genomic_length))
p <- p + geom_point(color="#E5E4DF")
p <- p + geom_smooth(method = "lm", colour = "black")
p <- p + scale_x_continuous(
  limits = c(0, 410),
  expand = c(0, 0),
  breaks = seq(0, 400, 50)
)
p <- p + scale_y_continuous(
  limits = c(0, 200),
  expand = c(0, 0),
  breaks = seq(0, 200, 50)
)
p <- p + theme_cowplot(12)
p <- p + xlab("Length (genes)")
p <- p + ylab("Genomic length (mb)")
p <- p + theme(
  axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
  axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
)

png(
  paste0(plot_dir, "gene_vs_genomic_length_scatter.png"), 
  width = 7, height = 5, unit = "in", res = 300
)
  p
dev.off()



  
  
  
  