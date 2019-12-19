#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

project_name <- "thesis"
subproject_name <- "Figure_2.3_simulations"
numcores <- 10
sample_name <- args[1]
subset_data <- as.logical(args[2])
downsample <- as.logical(args[3])
downsample_proportions <- as.numeric(
  unlist(
    strsplit(
      args[4],
      split = "_"
    )
  )
)
nUMI_threshold <- as.numeric(args[5])
nGene_threshold <- as.numeric(args[6])

#sample_name <- "CID4520N"
#subset_data <- FALSE
#downsample <- TRUE
#nUMI_threshold <- 25000
#downsample_proportions <- as.numeric(
#  unlist(
#    strsplit(
#      "0.05_0.1_0.15_0.2_0.3_0.4_0.5_0.6_0.7_0.8_0.9",
#      split = "_"
#    )
#  )
#)
#nGene_threshold <- 5000

RStudio <- FALSE

print(paste0("Project name = ", project_name))
print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("Subset data? ", as.character(subset_data)))

library(Seurat)
library(cowplot)
library(rlist)
library(naturalsort, lib.loc = "/share/ScratchGeneral/jamtor/R/3.5dev/")

if (RStudio) {
  home_dir <- "/Users/jamestorpy/clusterHome/"
} else {
  home_dir <- "/share/ScratchGeneral/jamtor/"
}
project_dir <- paste0(home_dir, "projects/", 
  project_name, "/", subproject_name, "/")
ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/functions/")
results_dir <- seurat_path <- paste0(project_dir, "results/")
in_dir <- paste0(project_dir, "raw_files/seurat_objects/", 
  sample_name, "/")
annot_dir <- paste0(results_dir, "infercnv/t_cells_included/", 
  sample_name, "/input_files/")

if (downsample) {
  out_path <- paste0(results_dir, "cancer_simulation/", sample_name, 
  	"_cancer_sim/downsampling/")
  out_dir <- paste0(results_dir, "infercnv/t_cells_included/", sample_name, 
    "_cancer_sim/downsampling/")
} else {
  out_path <- paste0(results_dir, "cancer_simulation/", sample_name, 
  	"_cancer_sim/")
  out_dir <- paste0(results_dir, "infercnv/t_cells_included/", sample_name, 
    "_cancer_sim/")
}
input_dir <- paste0(out_dir, "input_files/")
system(paste0("mkdir -p ", input_dir))

metadata_dir <- paste0(results_dir, "infercnv/t_cells_included/", sample_name, 
  "/input_files/")

Robject_dir <- paste0(out_path, "Rdata/")
system(paste0("mkdir -p ", Robject_dir))
plot_dir <- paste0(out_path, "plots/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(out_path, "tables/")
system(paste0("mkdir -p ", table_dir))


print(paste0("Sample directory = ", in_dir))
print(paste0("Reference directory = ", ref_dir))
print(paste0("Output directory = ", out_dir))
print(paste0("Plot directory = ", plot_dir))
print(paste0("Table directory = ", table_dir))

print(paste0("Generating simulated cancer data set from ", sample_name))

if (downsample) {
  # define CNVs for each chromosome:
  chromosome_CNV_splits_1 <- data.frame(
    row.names = paste0("chr", c(1, 2, 3, 4, 5,   6, 7)),
             no_parts =       c(1, 1, 1, 2, 2,   1, 2),
           multiplier =       c(1, 3, 0, 2, 0.5, 1, 1)
  )

  chromosome_CNV_intervals <- data.frame(
    lengths =    c(400, 400, 300, 300, 200,  200,
      1300, 150, 150, 125, 125, 100, 100, 75, 75),
    multiplier = c(3,   0,    2,  0.5,  1,    1,
      1,    3,    0,   3,   0,   3,   0,  3,  1),
    gaps = c(rep(300, 3), 600, rep(300, 2), 0, rep(300, 7), 0)
  )

  chromosome_CNV_splits_2 <- data.frame(
    row.names = paste0("chr", c(17, 18, 19, 20, 21, 22)),
             no_parts =       rep(1, 6),
           multiplier =       rep(1, 6)
  )

} else {
  # define CNVs for each chromosome:
  chromosome_CNV_splits_1 <- data.frame(
    row.names = paste0("chr", c(1, 2, 3, 4, 5,   6, 7)),
             no_parts =       c(1, 1, 1, 2, 2,   1, 2),
           multiplier =       c(1, 3, 0, 2, 0.5, 1, 1.5)
  )

  chromosome_CNV_intervals <- data.frame(
    lengths =    c(400, 400, 300, 300, 200,  200,    #1
      1300, 150, 150, 125, 125, 100, 100, 75, 75,    #2
      1600, 200, 200, 200, 200, 200,                 #3
      200, 200, 200, 200, 200),                      #4
    multiplier = c(3,   0,    2,  0.5,  1.5,  0.5,   #1
      1,   3,     0,   3,   0,   3,   0,  3,  0,     #2
      1,    3,   3,   3,    3,  3,                   #3
      0,   0,   0,    0,  0),                        #4
    gaps =       c(rep(300, 3), 600, rep(300, 2), 0, rep(300, 7), 0,  #1-4                     #1-5
      0,    200, 150, 100,  75,  300,                                 #5
      200, 150, 100,  75, 0)                                          #6
  )
}


################################################################################
### 0. Define functions ###
################################################################################

# add function to take a vector, split it into even parts define a multiplier 
# for every alternate coordinate:
split_into_parts <- function(vector, n, multiplier_value) {

  if (n==1) {

    split_coords <- data.frame(
      start = vector[1],
      end = vector[length(vector)],
      multiplier = multiplier_value
    )

  } else {

    split_vectors <- split(vector, cut(vector, n, labels = FALSE))
    split_coords <- as.data.frame(
      do.call(
        "rbind",
        lapply(split_vectors, function(x) {
          return(c(min(x), max(x)))
        })
      )
    )
    colnames(split_coords) <- c("start", "end")
    split_coords$multiplier <- 1

  }

  if (n>1 & n<9) {
    split_coords$multiplier[
      rep(c(TRUE, FALSE), length(split_coords$multiplier)/2)
    ] <- multiplier_value
  } else {
    split_coords$multiplier[
      c(
        1,
        (length(split_coords$multiplier)/4)+1,
        (2*(length(split_coords$multiplier)/4))+1,
        (3*(length(split_coords$multiplier)/4))+1
      )
    ] <- multiplier_value
  }


  return(split_coords)

}


###################################################################################
### 1. Fetch counts matrix from normal sample and plot nUMI and nGene           ###
###################################################################################

if (!file.exists(paste0(Robject_dir, "/1.original_epithelial_df.Rdata"))) {
  # load seurat object:
  seurat_10X <- readRDS(paste0(in_dir, "04_seurat_object_annotated.Rdata"))
  
  # create raw matrix input file and subset if necessary:
  count_df <- as.matrix(GetAssayData(seurat_10X , slot = "counts"))
  print(
    paste0(
      "Dimensions of count df = ", paste(as.character(dim(count_df)), collapse=",")
    )
  )
  
  # isolate epithelial cells from count_df:
  cell_annotation <- read.table(
    paste0(annot_dir, "metadata.txt"),
    header = F,
    sep = "\t"
  )
  epithelial_df <- count_df[
    ,as.character(cell_annotation$V1[cell_annotation$V2 == "Epithelial"])
  ]

  print(
    paste0(
      "Dimensions of epithelial df = ", paste(as.character(dim(epithelial_df)), collapse=",")
    )
  )
  
  if (downsample) {
  	# create density plots of nUMI and nGene:
    QC <- data.frame(
      row.names = colnames(epithelial_df),
      nUMI = apply(epithelial_df, 2, sum),
      nGene = apply(epithelial_df, 2, function(x) length(x[x!=0]))
    )
    QC <- QC[colnames(epithelial_df),]
  
    nUMI_density_plot <- density(QC$nUMI)
    pdf(paste0(plot_dir, "nUMI_density_plot.pdf"))
      plot(nUMI_density_plot, main=NA, xlab = "nUMI")
    dev.off()
    png(paste0(plot_dir, "nUMI_density_plot.png"))
      plot(nUMI_density_plot, main=NA, xlab = "nUMI")
    dev.off()
  
    nGene_density_plot <- density(QC$nGene)
    pdf(paste0(plot_dir, "nGene_density_plot.pdf"))
      plot(nGene_density_plot, main=NA, xlab = "nGene")
    dev.off()
    png(paste0(plot_dir, "nGene_density_plot.png"))
      plot(nGene_density_plot, main=NA, xlab = "nGene")
    dev.off()
  
    log_nUMI_density_plot <- density(log10(QC$nUMI))
    pdf(paste0(plot_dir, "log10_nUMI_density_plot.pdf"))
      plot(log_nUMI_density_plot, main=NA, xlab = "log10 nUMI")
    dev.off()
    png(paste0(plot_dir, "log10_nUMI_density_plot.png"))
      plot(log_nUMI_density_plot, main=NA, xlab = "log10 nUMI")
    dev.off()
  
    log_nGene_density_plot <- density(log10(QC$nGene))
    pdf(paste0(plot_dir, "log10_nGene_density_plot.pdf"))
      plot(log_nGene_density_plot, main=NA, xlab = "log10 nGene")
    dev.off()
    png(paste0(plot_dir, "log10_nGene_density_plot.png"))
      plot(log_nGene_density_plot, main=NA, xlab = "log10 nGene")
    dev.off()
  
    # filter out cells with nUMI < nUMI_threshold and nGene < nGene_threshold
    print(paste0("Number of cells before filtering out low coverage: ",
      nrow(QC)))
    cells_to_keep <- rownames(QC)[QC$nUMI > nUMI_threshold & QC$nGene > nGene_threshold]
    print(paste0("Number of cells after filtering out low coverage: ",
      length(cells_to_keep)))
    epithelial_df <- epithelial_df[
      ,colnames(epithelial_df) %in% cells_to_keep
    ]

    saveRDS(epithelial_df, paste0(Robject_dir, "/1.original_epithelial_df.Rdata"))
    
    p <- ggplot(QC, aes(x=nUMI, y=nGene))
    p <- p + geom_point()
    p <- p + xlab("nUMI")
    p <- p + ylab("nGene")
    p <- p + theme(legend.title = element_blank())
    pdf(paste0(plot_dir, "QC_quad_plot.pdf"), 
      width = 10, height = 6)
      print(p)
    dev.off()
    png(paste0(plot_dir, "QC_quad_plot.png"), 
      width = 450, height = 270)
      print(p)
    dev.off()
  
    # create total count density quad plot:
    print("Determining total count...")
    total_counts <- apply(epithelial_df, 2, sum)
    total_count_density_plot <- density(total_counts, bw="SJ")
    pdf(paste0(plot_dir, "total_count_density_plot.pdf"))
      plot(total_count_density_plot, main=NA, xlab = "Total counts")
    dev.off()
    png(paste0(plot_dir, "total_count_density_plot.png"))
      plot(total_count_density_plot, main=NA, xlab = "Total counts")
    dev.off()
  
    # create log10 total count density quad plot:
    log_total_count_density_plot <- density(log10(total_counts), bw="SJ")
    pdf(paste0(plot_dir, "log10_total_count_density_plot.pdf"))
      plot(log_total_count_density_plot, main=NA, xlab = "Total counts")
    dev.off()
    png(paste0(plot_dir, "log10_total_count_density_plot.png"))
      plot(log_total_count_density_plot, main=NA, xlab = "Total counts")
    dev.off()
  }

} else {
  epithelial_df <- readRDS(paste0(Robject_dir, "/1.original_epithelial_df.Rdata"))
}


################################################################################
### 2. Format annotation and counts matrices ###
################################################################################

if (!file.exists(paste0(Robject_dir, 
  "/2a.epithelial_df.Rdata")) | !file.exists(paste0(Robject_dir, 
  "/2b.gene_annotation.Rdata"))) {

  # load in gene annotation and deterime chromosome lengths:
  gene_annotation <- read.table(
    paste0(ref_dir, "infercnv_gene_order.txt"),
    header = F,
    sep = "\t",
    as.is = T
  )
  colnames(gene_annotation) <- c("symbol", "chromosome", "start", "end")
  
  # subset gene annotation and epithelial_df so they contain the same genes:
  rownames(gene_annotation) <- gene_annotation$symbol

  genes_not_in_gene_annotation <- rownames(epithelial_df)[
    !(rownames(epithelial_df) %in% rownames(gene_annotation))
  ]
  print(paste0(length(genes_not_in_gene_annotation), 
    " genes not present in InferCNV gene annotation but present in ", sample_name, 
    " counts matrix"))
  
  genes_not_in_epithelial_df <- rownames(gene_annotation)[
    !(rownames(gene_annotation) %in% rownames(epithelial_df))
  ]
  print(paste0(length(genes_not_in_epithelial_df), " genes not present in ", sample_name, 
    " counts matrix but present in InferCNV gene annotation"))
  
  print("Removing genes not present in both InferCNV gene annotation and counts matrix...")
  gene_annotation <- gene_annotation[
    rownames(gene_annotation) %in% rownames(epithelial_df),
  ]
  epithelial_df <- epithelial_df[rownames(epithelial_df) %in% rownames(gene_annotation),]

  print(paste0("Gene numbers of InferCNV gene annotation and counts matrix are now ",
    nrow(gene_annotation), " and ", nrow(epithelial_df), " respectively"))

  # remove genes in chromosomes M, X, Y:
  print("Removing genes from chromosomes M, X, Y...")
  gene_annotation <- gene_annotation[
    !(gene_annotation$chromosome %in% c("chrM", "chrX", "chrY")),
  ]
  epithelial_df <- epithelial_df[
    gene_annotation$symbol[!(gene_annotation$chromosome %in% c("chrM", "chrX", "chrY"))],
  ]

  print(paste0("Gene numbers of InferCNV gene annotation and counts matrix are now ",
    nrow(gene_annotation), " and ", nrow(epithelial_df), " respectively"))
  
  # number genes in gene annotation:
  gene_annotation$number <- seq(1, nrow(gene_annotation))

  saveRDS(epithelial_df, paste0(Robject_dir, "/2a.epithelial_df.Rdata"))
  saveRDS(gene_annotation, paste0(Robject_dir, 
    "/2b.gene_annotation.Rdata"))

} else {

  epithelial_df <- readRDS(paste0(Robject_dir, "/2a.epithelial_df.Rdata"))
  gene_annotation <- readRDS(paste0(Robject_dir, 
    "/2b.gene_annotation.Rdata"))

}


################################################################################
### 3. Prepare chromosome information ###
################################################################################

# determine chromosome information:
chromosome_lengths <- lapply(split(gene_annotation, gene_annotation$chromosome), 
  nrow)
chromosome_lengths <- chromosome_lengths[
  naturalsort(names(chromosome_lengths))
]

# subset to first 4 chromosomes:
if (subset_data) {
  chromosome_lengths <- chromosome_lengths[1:4]
  gene_annotation <- gene_annotation[
    gene_annotation$chromosome %in% c("chr1", "chr2", "chr3", "chr4"),
  ]
  epithelial_df <- epithelial_df[gene_annotation$symbol,]
}

for (n in 1:length(chromosome_lengths)) {
  print(n)
  if (n==1) {
    chromosome_coords <- list(1:chromosome_lengths[[n]])
  } else {
    # define starting coordinate as last coordinate of previous chromosome + 1:
    start_coord <- chromosome_coords[[n-1]][length(chromosome_coords[[n-1]])]+1
    chromosome_coords[[n]] <- start_coord:(start_coord+chromosome_lengths[[n]]-1)
  }
}
names(chromosome_coords) <- names(chromosome_lengths)

# determine chromosome ends and midpoints:
chromosome_ends <- lapply(chromosome_coords, max)
chromosome_midpoints <- lapply(chromosome_coords, function(x) x[floor(length(x)/2)])


################################################################################
### 4. Plot average fold difference from mean ###
################################################################################

# generate average original counts vector with each value representing a gene:
average_original_counts <- apply(epithelial_df, 1, mean)
# add 0.1 to all values:
average_original_counts <- average_original_counts + 3e-4
#average_original_counts[average_original_counts == 0] <- 3e-4
# determine mean:
mean_average_original_counts <- mean(average_original_counts)

# divide by mean to get fold change from mean and add to df for plotting:
original_fold_change <- average_original_counts/mean_average_original_counts
original_fold_change_df <- data.frame(
  number = seq(1, length(original_fold_change)),
  count = original_fold_change
)
# take mean of original fold change to mark on plot:
mean_original_fold_change <- mean(original_fold_change)

# plot fold change values:
if (!file.exists(paste0(plot_dir, "0a.original_fold_change_from_mean.pdf"))) {
  p <- ggplot(original_fold_change_df, aes(x=number, y=count))
  p <- p + geom_point(colour = "#E8D172")
  p <- p + xlab("Genomic coordinates")
  #p <- p + xlim(c(0, nrow(centered_original_counts_df)))
  p <- p + scale_x_continuous(
    breaks = unlist(chromosome_midpoints),
    labels = paste0("chr", 1:length(chromosome_midpoints)), 
    expand = c(0, 0)
  )
  p <- p + ylab("Fold change")
  p <- p + scale_y_continuous(
    breaks = c(min(original_fold_change), 1, 20),
    labels = c("Zero counts", "1", "20"),
    limits = c(min(original_fold_change), 20)
  )
  for (end in chromosome_ends) {
    p <- p + geom_vline(xintercept=end)
  }
  p <- p + geom_segment(
    x=0, 
    xend=max(unlist(chromosome_ends)), 
    y=mean_original_fold_change,
    yend=mean_original_fold_change, 
    size=1, color="red"
  )

  pdf(paste0(plot_dir, "0a.original_fold_change_from_mean.pdf"), 
  	width = 20)
    print(p)
  dev.off()

}

# take the log10:
log_original_fold_change <- log10(original_fold_change)
# tabulate for plotting:
log_original_fold_change_df <- data.frame(
  number = seq(1, length(log_original_fold_change)),
  count = log_original_fold_change
)
# calculate the mean of the counts vector:
log_mean_original_fold_change <- log10(mean_original_fold_change)

# plot counts:
if (!file.exists(paste0(plot_dir, "0b.log_original_fold_change_from_mean.pdf"))) {
  p <- ggplot(log_original_fold_change_df, aes(x=number, y=count))
  p <- p + geom_point(colour = "#E8D172")
  p <- p + xlab("Genomic coordinates")
  #p <- p + xlim(c(0, nrow(centered_original_counts_df)))
  p <- p + scale_x_continuous(
    breaks = unlist(chromosome_midpoints),
    labels = paste0("chr", 1:length(chromosome_midpoints)),
    #limits = c(0,nrow(centered_original_counts_df)), 
    expand = c(0, 0)
  )
  p <- p + ylab("Log10 fold change from mean")
  p <- p + scale_y_continuous(
    breaks = c(min(log_original_fold_change), -1, 0, 1, 2, 3, 4),
    labels = c("Zero counts", "-1", "0", "1", "2", "3", "4"),
  )
  for (end in chromosome_ends) {
    p <- p + geom_vline(xintercept=end)
  }
  p <- p + geom_segment(
    x=0, 
    xend=max(unlist(chromosome_ends)), 
    y=log_mean_original_fold_change, 
    yend=log_mean_original_fold_change, 
    size=1, color="red"
  )

  pdf(paste0(plot_dir, "0b.log10_original_fold_change_from_mean.pdf"), width = 20)
    print(p)
  dev.off()
}


################################################################################
### 5. Plot average fold difference from median ###
################################################################################

# generate average original counts vector with each value representing a gene:
average_original_counts <- apply(epithelial_df, 1, mean)
# add 0.1 to zero values:
average_original_counts <- average_original_counts + 3e-4
# determine median:
median_average_original_counts <- median(average_original_counts)

# divide by median to get fold change from median and add to df for plotting:
original_fold_change <- average_original_counts/median_average_original_counts
original_fold_change_df <- data.frame(
  number = seq(1, length(original_fold_change)),
  count = original_fold_change
)
# take mean of original fold change to mark on plot:
median_original_fold_change <- median(original_fold_change)

# plot fold change values:
if (!file.exists(paste0(plot_dir, "0c.original_fold_change_from_median.pdf"))) {
  p <- ggplot(original_fold_change_df, aes(x=number, y=count))
  p <- p + geom_point(colour = "#E8D172")
  p <- p + xlab("Genomic coordinates")
  #p <- p + xlim(c(0, nrow(centered_original_counts_df)))
  p <- p + scale_x_continuous(
    breaks = unlist(chromosome_midpoints),
    labels = paste0("chr", 1:length(chromosome_midpoints)),
    #limits = c(0,nrow(centered_original_counts_df)), 
    expand = c(0, 0)
  )
  p <- p + ylab("Fold change")
  p <- p + scale_y_continuous(
    breaks = c(0.0005, 1, 20),
    labels = c("Zero counts", "1", "20"),
    limits = c(0.0005, 25)
  )
  for (end in chromosome_ends) {
    p <- p + geom_vline(xintercept=end)
  }
  p <- p + geom_segment(
    x=0, 
    xend=max(unlist(chromosome_ends)), 
    y=median_original_fold_change,
    yend=median_original_fold_change, 
    size=1, color="red"
  )

  pdf(paste0(plot_dir, "0c.original_fold_change_from_median.pdf"), width = 20)
    print(p)
  dev.off()
}

# take the log10:
log_original_fold_change <- log10(original_fold_change)
# tabulate for plotting:
log_original_fold_change_df <- data.frame(
  number = seq(1, length(log_original_fold_change)),
  count = log_original_fold_change
)
# calculate the mean of the counts vector:
log_median_original_fold_change <- log10(median_original_fold_change)

# plot counts:
if (!file.exists(paste0(plot_dir, "1.log_original_fold_change_from_median.pdf"))) {
  p <- ggplot(log_original_fold_change_df, aes(x=number, y=count))
  p <- p + geom_point(colour = "#E8D172")
  p <- p + xlab("Genomic coordinates")
  #p <- p + xlim(c(0, nrow(centered_original_counts_df)))
  p <- p + scale_x_continuous(
    breaks = unlist(chromosome_midpoints),
    labels = paste0("chr", 1:length(chromosome_midpoints)),
    #limits = c(0,nrow(centered_original_counts_df)), 
    expand = c(0, 0)
  )
  p <- p + ylab("Log10 fold change from median")
  p <- p + scale_y_continuous(
    breaks = c(min(log_original_fold_change), -1, 0, 1, 2, 3, 4),
    labels = c("Zero counts", "-1", "0", "1", "2", "3", "4"),
    #limits = c(-5, 8)
  )
  for (end in chromosome_ends) {
    p <- p + geom_vline(xintercept=end)
  }
  p <- p + geom_segment(
    x=0, 
    xend=max(unlist(chromosome_ends)), 
    y=log_median_original_fold_change, 
    yend=log_median_original_fold_change, 
    size=1, color="red"
  )

  pdf(paste0(plot_dir, "1.log10_original_fold_change_from_median.pdf"), width = 20)
    print(p)
  dev.off()
}


################################################################################
### 6. Generate map of CNVs per chromosome ###
################################################################################

# add first intervals to split chromosomes into CNV regions:
for (i in 1:nrow(chromosome_CNV_splits_1)) {

  print(paste0("Creating CNV profile for ", 
    rownames(chromosome_CNV_splits_1)[i]))

  temp_df <- split_into_parts(
    eval(parse(text = paste0("chromosome_coords$", 
      rownames(chromosome_CNV_splits_1)[i]))),
    n=chromosome_CNV_splits_1$no_parts[i],
    multiplier_value=chromosome_CNV_splits_1$multiplier[i]
  )
  temp_df$start_chromosome <- rownames(chromosome_CNV_splits_1)[i]
  temp_df$end_chromosome <- rownames(chromosome_CNV_splits_1)[i]

  if (i==1) {
    CNV_indices <- temp_df
  } else {
    CNV_indices <- rbind(CNV_indices, temp_df)
  }
}
CNV_indices_nrow <- nrow(CNV_indices)

# add intervals of small CNV regions interspersed with CNV 
# neutral regions:
for (j in 1:nrow(chromosome_CNV_intervals)) {
  print(j)

  last_position <- CNV_indices$end[nrow(CNV_indices)]

  # define next CNA interval:
  CNA_segment <- data.frame(
    start = last_position+1,
    end = last_position+1+chromosome_CNV_intervals$lengths[j],
    multiplier = chromosome_CNV_intervals$multiplier[j]
  )
  # identify which chromosome contains the start interval:
  l=1
  start_containing_chromosome <- lapply(chromosome_coords, function(x) {
    if (CNA_segment$start %in% x) {
      l <<- l+1
      return(names(chromosome_coords)[l-1])
    } else {
      l <<- l+1
      return(NULL)
    }
  })
  start_containing_chromosome <- unlist(start_containing_chromosome)
  CNA_segment$start_chromosome <- start_containing_chromosome

  # identify which chromosome contains the end interval:
  l=1
  end_containing_chromosome <- lapply(chromosome_coords, function(x) {
    if (CNA_segment$end %in% x) {
      l <<- l+1
      return(names(chromosome_coords)[l-1])
    } else {
      l <<- l+1
      return(NULL)
    }
  })
  end_containing_chromosome <- unlist(end_containing_chromosome)
  CNA_segment$end_chromosome <- end_containing_chromosome

  CNV_indices <- rbind(CNV_indices, CNA_segment)

  # define CNV neutral interval of 300 genes:
  CNV_segment <- data.frame(
    start = CNA_segment$end+1,
    end = CNA_segment$end+1+chromosome_CNV_intervals$gaps[j],
    multiplier = 1
  )
  # identify which chromosome contains the start interval:
  l=1
  start_containing_chromosome <- lapply(chromosome_coords, function(x) {
    if (CNV_segment$start %in% x) {
      l <<- l+1
      return(names(chromosome_coords)[l-1])
    } else {
      l <<- l+1
      return(NULL)
    }
  })
  start_containing_chromosome <- unlist(start_containing_chromosome)
  CNV_segment$start_chromosome <- start_containing_chromosome
  # identify which chromosome contains the end interval:
  l=1
  end_containing_chromosome <- lapply(chromosome_coords, function(x) {
    if (CNV_segment$end %in% x) {
      l <<- l+1
      return(names(chromosome_coords)[l-1])
    } else {
      l <<- l+1
      return(NULL)
    }
  })
  end_containing_chromosome <- unlist(end_containing_chromosome)
  CNV_segment$end_chromosome <- end_containing_chromosome

  # add next CNV neutral interval:
  CNV_indices <- rbind(CNV_indices, CNV_segment)

}

# remove CNV_indices entries with length = 0:
CNV_indices <- CNV_indices[CNV_indices$start - CNV_indices$end != 0,]

# fill in rest of last chromosome with copy number neutral region:
end_position <- CNV_indices$end[nrow(CNV_indices)]
end_chromosome <- CNV_indices$end_chromosome[nrow(CNV_indices)]
CNV_segment <- data.frame(
  start = end_position+1,
  end = eval(parse(text=paste0("chromosome_ends$", end_chromosome))),
  multiplier = 1,
  start_chromosome = end_chromosome,
  end_chromosome = end_chromosome
)
CNV_indices <- rbind(CNV_indices, CNV_segment)

# add second intervals to split chromosomes into CNV regions:
if (exists("chromosome_CNV_splits_2")) {
  for (i in 1:nrow(chromosome_CNV_splits_2)) {

    print(paste0("Creating CNV profile for ", 
      rownames(chromosome_CNV_splits_2)[i]))
  
    temp_df <- split_into_parts(
      eval(parse(text = paste0("chromosome_coords$", 
        rownames(chromosome_CNV_splits_2)[i]))),
      n=chromosome_CNV_splits_2$no_parts[i],
      multiplier_value=chromosome_CNV_splits_2$multiplier[i]
    )
    temp_df$start_chromosome <- rownames(chromosome_CNV_splits_2)[i]
    temp_df$end_chromosome <- rownames(chromosome_CNV_splits_2)[i]
  
    CNV_indices <- rbind(CNV_indices, temp_df)
  }
  CNV_indices_nrow <- nrow(CNV_indices)
}


################################################################################
### 7. Apply CNV changes to epithelial_df ###
################################################################################
  
if (!file.exists(paste0(Robject_dir, "/2.simulated_epithelial_df.Rdata")) | 
  !file.exists(paste0(Robject_dir, "/CNV_indices.Rdata")) |
  !file.exists(paste0(Robject_dir, "final_fold_change.Rdata"))) {
  
  print("Adding CNVs to counts matrix...")

  for (r in 1:nrow(CNV_indices)) {
    
    print(r)

    # generate mean original segment counts vector with each value representing 
    # a gene:
    average_original_counts <- apply(
      epithelial_df[CNV_indices$start[r]:CNV_indices$end[r],], 1, mean
    )
    # add 0.1 to all values:
  	#average_original_counts[average_original_counts == 0] <- 1e-3
    average_original_counts <- average_original_counts + 1e-3
  	# determine median:
  	median_average_original_counts <- median(average_original_counts)
  	# divide by median to get fold change from median and add to df for plotting:
  	original_fold_change <- average_original_counts/median_average_original_counts
  	# take median of original fold change to mark on plot:
  	median_original_fold_change <- median(original_fold_change)

    # multiply adjusted of epithelial_df according to CNV_indices$multiplier:
    epithelial_df[CNV_indices$start[r]:CNV_indices$end[r],] <- 
      epithelial_df[CNV_indices$start[r]:CNV_indices$end[r],]*CNV_indices$multiplier[r]

    # generate mean adjusted counts vector with each value representing 
    # a gene:
    average_adjusted_counts <- apply(
      epithelial_df[CNV_indices$start[r]:CNV_indices$end[r],], 1, mean
    )
    # add 0.1 to zero values:
	  average_adjusted_counts <- average_adjusted_counts + 1e-3
	  # divide by median to get fold change from original median and add to df for plotting:
	  adjusted_fold_change <- average_adjusted_counts/median_average_original_counts
	  # take median of adjusted fold change to mark on plot:
	  median_adjusted_fold_change <- median(adjusted_fold_change)    

	  # take the log10:
    log_adjusted_fold_change <- log10(adjusted_fold_change)
    # if any of log_adjusted_fold_change < zero_count_value, 
    # adjust to zero_count_value:
    if (r > 3) {
      log_adjusted_fold_change[
        log_adjusted_fold_change < zero_count_value
      ] <- zero_count_value
    }

    # add to table for plotting:
    if (r==1) {
      log_adjusted_fold_change_df <- data.frame(
        number = seq(1, length(log_adjusted_fold_change)),
        count = log_adjusted_fold_change
      )
    } else {
      log_adjusted_fold_change_df <- rbind(
        log_adjusted_fold_change_df,
        data.frame(
          number = seq(
            nrow(log_adjusted_fold_change_df)+1, 
            nrow(log_adjusted_fold_change_df)+length(log_adjusted_fold_change)
          ),
          count = log_adjusted_fold_change
        )
      )
    }

    # calculate the mean of the counts vector and add to CNV_indices:
    log_median_adjusted_fold_change <- log10(median_adjusted_fold_change)
    CNV_indices$segment_log_median_FC[r] <- log_median_adjusted_fold_change

    # define loss value:
    if (r==3) {
      zero_count_value <- log_median_adjusted_fold_change
    }
    # if log_median_adjusted_fold_change < zero_count_value, adjust to zero_count_value:
    if (r > 3) {
      if (CNV_indices$multiplier[r] == 0) {
        CNV_indices$segment_log_median_FC[r] <- zero_count_value
      }
    }
  }


  ################################################################################
  ### 8. Create adjusted line only CNV plot ###
  ################################################################################

  if (!file.exists(paste0(plot_dir, 
    "2a.log_adjusted_fold_change_from_median_line_only.pdf"))) {
    p <- ggplot(log_adjusted_fold_change_df, aes(x=number, y=count))
    p <- p + scale_x_continuous(
      breaks = unlist(chromosome_midpoints),
      labels = paste0("chr", 1:length(chromosome_midpoints)),
      limits = c(0,length(log_adjusted_fold_change_df$count)), 
      expand = c(0, 0)
    )
    p <- p + scale_y_continuous(
      breaks = c(zero_count_value, -1, 0, 1),
      labels = c("Zero counts", "-1", "0", "1"),
      limits = c(zero_count_value, 1)
    )
    for (end in chromosome_ends) {
      p <- p + geom_vline(xintercept=end)
    }
    for (r in 1:nrow(CNV_indices)) {
      print(r)
      # create horizontal line:
      p <- p + geom_segment(
        x=CNV_indices$start[r], 
        xend=CNV_indices$end[r], 
        y=CNV_indices$segment_log_median_FC[r], 
        yend=CNV_indices$segment_log_median_FC[r], 
        size=1, color="#37841f"
      )
      # create left vertical line:
      if (r != 1) {
        p <- p + geom_segment(
          x=CNV_indices$start[r], 
          xend=CNV_indices$start[r], 
          y=CNV_indices$segment_log_median_FC[r-1], 
          yend=CNV_indices$segment_log_median_FC[r], 
          size=1, color="#37841f"
        )
      }
    }
  
    pdf(paste0(plot_dir, 
      "2a.log_adjusted_fold_change_from_median_line_only.pdf"), width = 20)
      print(p)
    dev.off()
  }


  ################################################################################
  ### 9. Create adjusted counts CNV plot ###
  ################################################################################

  # plot counts:
  if (!file.exists(paste0(plot_dir, "2b.log_adjusted_fold_change_from_median.pdf"))) {
    p <- ggplot(log_adjusted_fold_change_df, aes(x=number, y=count))
    p <- p + geom_point(colour = "#E8D172")
    p <- p + xlab("Genomic coordinates")
    p <- p + scale_x_continuous(
      breaks = unlist(chromosome_midpoints),
      labels = paste0("chr", 1:length(chromosome_midpoints)),
      limits = c(0,nrow(log_adjusted_fold_change_df)), 
      expand = c(0, 0)
    )
    p <- p + ylab("Log10 fold change")
    p <- p + scale_y_continuous(
      breaks = c(zero_count_value, -1, 0, 1, 2, 3, 4),
      labels = c("Zero counts", "-1", "0", "1", "2", "3", "4"),
      limits = c(zero_count_value, max(log_adjusted_fold_change_df$count))
    )
    for (end in chromosome_ends) {
      p <- p + geom_vline(xintercept=end)
    }
    for (r in 1:nrow(CNV_indices)) {
      print(r)
      # create horizontal line:
      p <- p + geom_segment(
        x=CNV_indices$start[r], 
        xend=CNV_indices$end[r], 
        y=CNV_indices$segment_log_median_FC[r], 
        yend=CNV_indices$segment_log_median_FC[r], 
        size=1, color="red"
      )
      # create left vertical line:
      if (r != 1) {
        p <- p + geom_segment(
          x=CNV_indices$start[r], 
          xend=CNV_indices$start[r], 
          y=CNV_indices$segment_log_median_FC[r-1], 
          yend=CNV_indices$segment_log_median_FC[r], 
          size=1, color="red"
        )
      }
    }
  
    pdf(paste0(plot_dir, "2b.log_adjusted_fold_change_from_median.pdf"), width = 20)
      print(p)
    dev.off()
  }

  simulated_CNV_plot_data <- list(
    log_adjusted_fold_change_df,
    zero_count_value,
    CNV_indices
  )
  names(simulated_CNV_plot_data) <- c("log_adjusted_fold_change_df", 
    "zero_count_value", "CNV_indices")
  saveRDS(simulated_CNV_plot_data, 
    paste0(Robject_dir, "simulated_CNV_plot_data.Rdata")
  )
  saveRDS(epithelial_df, paste0(Robject_dir, "/2.simulated_epithelial_df.Rdata"))
  saveRDS(CNV_indices, paste0(Robject_dir, "/CNV_indices.Rdata"))
  
} else {
  
  epithelial_df <- readRDS(paste0(Robject_dir, "/2.simulated_epithelial_df.Rdata"))
  CNV_indices <- readRDS(paste0(Robject_dir, "/CNV_indices.Rdata"))

}


################################################################################
### 6. Bind new epithelial and original stromal counts, check against metadata
# and save as infercnv input files ###
################################################################################

if (!exists("count_df")) {
  # load seurat object:
  seurat_10X <- readRDS(paste0(in_dir, "04_seurat_object_annotated.Rdata"))
  # create raw matrix input file and subset if necessary:
  count_df <- as.matrix(GetAssayData(seurat_10X , slot = "counts"))
}

# only keep genes in both count and epithelial_df:
print(paste0(
  "Number of genes in count_df before removing genes not in epithelial_df = ", 
  nrow(count_df))
)
count_df <- count_df[rownames(epithelial_df),]
print(paste0(
  "Dimensions of count_df after removing genes not in epithelial_df = ", 
  nrow(count_df))
)

if (!exists("cell_annotation")) {
  cell_annotation <- read.table(
    paste0(annot_dir, "metadata.txt"),
    header = F,
    sep = "\t"
  )
}

print(
  paste0(
    "Columns in stromal part of count_df = ", 
    ncol(
      count_df[
        ,as.character(cell_annotation$V1[cell_annotation$V2 != "Epithelial"])
      ]
    )
  )
)
print(paste0("Columns in epithelial_df = ", ncol(epithelial_df)))

# bind stromal count_df and epithelial_df together as new_counts:
new_counts <- cbind(
  epithelial_df,
  count_df[
    ,as.character(cell_annotation$V1[cell_annotation$V2 != "Epithelial"])
  ]
)
print(paste0("Columns in new_counts = ", ncol(new_counts)))

# only keep cells in cell_annotation:
new_counts <- new_counts[
  ,colnames(new_counts) %in% as.character(cell_annotation$V1)
]
print(paste0(
  "Columns in new_counts after removing cells not in cell_annotation = ", 
  ncol(new_counts))
)

write.table(new_counts, paste0(out_dir, "input_files/input_matrix.txt"), quote=F,
  sep="\t", col.names=T, row.names=T)

# copy original metadata and ensure it has same cells as epithelial_df:
metadata_df <- read.table(paste0(metadata_dir, "metadata.txt"),
  header = F, sep = "\t", as.is = T)
rownames(metadata_df) <- metadata_df$V1
epithelial_metadata <- metadata_df[colnames(epithelial_df),]
metadata_df <- rbind(metadata_df[metadata_df$V2 == "Stromal",],
  epithelial_metadata)

write.table(metadata_df, paste0(input_dir, "metadata.txt"), sep = "\t",
  quote = F, col.names = F, row.names = F)


################################################################################
### 7. Downsample new counts and save ###
################################################################################

if (downsample) {
  library(DropletUtils)

  print("Total counts before downsampling:")
  print(sum(as.vector(new_counts)))
  
  for (proportion in downsample_proportions) {
  
    print(paste0("Downsampling counts by ", proportion, "..."))
    downsampled_counts <- downsampleMatrix(new_counts, proportion, bycol=F)
    print("Total counts after downsampling:")
    print(sum(as.vector(downsampled_counts)))
  
    downsample_dir <- paste0(out_dir, proportion, "_downsampling/input_files/")
    system(paste0("mkdir -p ", downsample_dir))
    write.table(downsampled_counts, paste0(downsample_dir, "input_matrix.txt"), quote=F,
      sep="\t", col.names=T, row.names=T)
  
    write.table(metadata_df, paste0(downsample_dir, "metadata.txt"), sep = "\t",
      quote = F, col.names = F, row.names = F)
  
  }
}


################################################################################
### 7. Convert PDF to PNG ###
################################################################################

system(paste0("for p in ", plot_dir, "*.pdf; do echo $p; f=$(basename $p); echo $f; ",
              "new=$(echo $f | sed 's/.pdf/.png/'); echo $new; ", 
              "convert -density 150 ", plot_dir, "$f -quality 90 ", plot_dir, "$new; done"))


