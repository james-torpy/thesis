#! /share/ClusterShare/software/contrib/briglo/octR/src/R-3.6.0/builddir/bin/Rscript


###############################################################################
### Simulates cancer datasets taking as inputs number of CNVs, CNV lengths, ###
### gain/loss multipliers and downsamples required                          ###
###############################################################################

args = commandArgs(trailingOnly=TRUE)

project_name <- "thesis"
subproject_name <- 
  "Figure_1.2_1.3_1.5_simulate_cancer_define_denoising_value_distinguish_copy_no"
sample_name <- args[1]
print(paste0("sample name = ", sample_name))
nUMI_threshold <- as.numeric(args[2])
print(paste0("nUMI_threshold = ", nUMI_threshold))
nGene_threshold <- as.numeric(args[3])
print(paste0("nGene_threshold = ", nGene_threshold))
# range of number of CNVs per simulation, chosen randomly
# defines number of CNV loop iterations with one CNV added per iteration:
CNV_no_range <- as.numeric(
  unlist(
    strsplit(
      args[4],
      split = "_"
    )
  )
)
print(paste0("Range of CNVs possible = ", paste(CNV_no_range, collapse = "_")))

# range of lengths of CNVs, chosen randomly:
CNV_lengths <- as.numeric(
  unlist(
    strsplit(
      args[5],
      split = "_"
    )
  )
)
print(paste0("Possible CNV lengths = ", paste(CNV_lengths, collapse = "_")))

# vector of possible gain/loss multipliers, chosen randomly:
CNV_multipliers <- as.numeric(
  unlist(
    strsplit(
      args[6],
      split = "_"
    )
  )
)
print(paste0("Possible CNV multipliers = ", paste(CNV_multipliers, collapse = "_")))

sim_name <- args[7]
print(paste0("Simulation name = ", sim_name))

noise_cell_no <- as.numeric(args[8])
print(paste0("Noise input cell number = ", noise_cell_no))

#project_name <- "thesis"
#subproject_name <- 
#  "Figure_1.2_1.3_1.5_simulate_cancer_define_denoising_value_distinguish_copy_no"
#sample_name <- "CID4520N"
#nUMI_threshold <- 25000
#nGene_threshold <- 5000
## range of number of CNVs per simulation, chosen randomly
## defines number of CNV loop iterations with one CNV added per iteration:
#CNV_no_range <- as.numeric(
#  unlist(
#    strsplit(
#      "10_50",
#      split = "_"
#    )
#  )
#)
## range of lengths of CNVs, chosen randomly:
#CNV_lengths <- as.numeric(
#  unlist(
#    strsplit(
#      "50_75_100_150_150_200_200_250_250_300_300_350_400_450_500_550_600_650_700_750_800",
#      split = "_"
#    )
#  )
#)
## vector of possible gain/loss multipliers, chosen randomly:
#CNV_multipliers <- as.numeric(
#  unlist(
#    strsplit(
#      "3_2_1.5_0.5_0",
#      split = "_"
#    )
#  )
#)
#sim_name <- "sim3"
#noise_cell_no <- 5000

print(paste0("Project name = ", project_name))
print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(Seurat)
library(cowplot)
library(rlist, lib.loc = lib_loc)
library(GenomicRanges)
library(naturalsort, lib.loc = lib_loc)
library(splatter, lib.loc = lib_loc)
library(ggplot2)
library(org.Hs.eg.db)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", 
  project_name, "/", subproject_name, "/")
ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/functions/")
results_dir <- paste0(project_dir, "results/")
annot_dir <- paste0(results_dir, "infercnv/", sample_name, 
  "/normal/input_files/")

raw_dir <- paste0(project_dir, "raw_files/")
seurat_dir <- paste0(raw_dir, "seurat_objects/", sample_name, "/")
emptydrops_dir <- paste0(raw_dir, "/emptydrops/", sample_name, "/")
normal_dir <- paste0(results_dir, "infercnv/", sample_name, 
  "/normal/input_files/")
annot_dir <- paste0(results_dir, "infercnv/", sample_name, 
  "/normal/Rdata/")

sim_path <- paste0(results_dir, "cancer_simulation/", sample_name, "/")
noise_dir <- paste0(sim_path, "noise_generation/")
system(paste0("mkdir -p ", noise_dir))

common_Robject_dir <- paste0(sim_path, "Rdata/")
system(paste0("mkdir -p ", common_Robject_dir))
common_plot_dir <- paste0(sim_path, "plots/")
system(paste0("mkdir -p ", common_plot_dir))
common_table_dir <- paste0(sim_path, "tables/")
system(paste0("mkdir -p ", common_table_dir))

Robject_dir <- paste0(sim_path, "/", sim_name, "/Rdata/")
system(paste0("mkdir -p ", Robject_dir))
plot_dir <- paste0(sim_path, "/", sim_name, "/plots/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(sim_path, "/", sim_name, "/tables/")
system(paste0("mkdir -p ", table_dir))

out_dir <- paste0(results_dir, "infercnv/", sample_name, "/", 
  sim_name, "/input_files/")
system(paste0("mkdir -p ", out_dir))

filtered_dir <- paste0(results_dir, "infercnv/", sample_name, 
  "/filtered_normal/input_files/")
system(paste0("mkdir -p ", filtered_dir))

print(paste0("Sample directory = ", normal_dir))
print(paste0("Reference directory = ", ref_dir))
print(paste0("Output directory = ", out_dir))
print(paste0("Plot directory = ", plot_dir))
print(paste0("Table directory = ", table_dir))

print(paste0("Generating simulated cancer data set from ", sample_name))
print(paste0("Filtering out cells with less than ", nUMI_threshold, " UMIs and ",
  nGene_threshold, " genes"))


################################################################################
### 0. Define functions and choose seeds ###
################################################################################

# choose and record seeds:
if (file.exists(paste0(table_dir, "random_seed_record.txt"))) {
  seed_record <- read.table(
    paste0(table_dir, "random_seed_record.txt"),
    sep = "\t",
    header = T,
    as.is = T
  )
} else {
  seed_record <- data.frame(
    row.names = c("CNV_no", "start_position", "length", "multiplier"),
    seed = sample(1:999, 4)
  )
  write.table(
    seed_record, 
    paste0(table_dir, "random_seed_record.txt"),
    sep = "\t",
    row.names = T,
    col.names = T,
    quote = F
  )
}

set.seed(seed_record["CNV_no",])
CNV_no <- sample(seq(CNV_no_range[1], CNV_no_range[2]), 1)
print(paste0(
  "Adding ", CNV_no, " CNVs ",
  "with lengths between ", CNV_lengths[1], " and ", CNV_lengths[length(CNV_lengths)]
))


###################################################################################
### 1. Prepare input data ###
###################################################################################

if (
  !file.exists(paste0(common_Robject_dir, "/1a.epithelial_df.Rdata")) | 
  !file.exists(paste0(common_Robject_dir, "/1b.filtered_epithelial_df.Rdata")) | 
  !file.exists(paste0(common_Robject_dir, "/1c.non_epithelial_df.Rdata")) | 
  !file.exists(paste0(common_Robject_dir, "/1d.metadata.Rdata")) |
  !file.exists(paste0(common_Robject_dir, "/1e.gene_annotation.Rdata")) |
  !file.exists(paste0(filtered_dir, "input_matrix.txt"))
) {

  # load normal counts matrix and split into epithelial and non-epithelial dfs:
  count_df <- read.table(
    paste0(normal_dir, "input_matrix.txt"),
    sep = "\t",
    header = T,
    as.is = T
  )
  
  metadata <- read.table(
    paste0(normal_dir, "metadata.txt"),
    sep = "\t",
    header = F,
    as.is = T
  )
  colnames(metadata) <- c("cell_ids", "cell_types")
  
  epithelial_df <- count_df[
    , colnames(count_df) %in% metadata$cell_ids[metadata$cell_types == "Epithelial"]
  ]
  non_epithelial_df <- count_df[
    , colnames(count_df) %in% metadata$cell_ids[metadata$cell_types == "Non_epithelial"]
  ]
  
  print("count_df dimensions:")
  print(dim(count_df))
  print("epithelial_df dimensions:")
  print(dim(epithelial_df))
  print("non_epithelial_df dimensions:")
  print(dim(non_epithelial_df))

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
  filtered_epithelial_df <- epithelial_df[
    ,colnames(epithelial_df) %in% cells_to_keep
  ]
  
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
  
  # load gene annotation data:
  gene_annotation <- readRDS(paste0(annot_dir, "2c.gene_annotation.Rdata"))

  # bind non_epithelial_df and modified_df together as new_counts:
  print("filtered_epithelial_df dimensions:")
  print(dim(filtered_epithelial_df))
  
  filtered_counts <- cbind(
    filtered_epithelial_df,
    non_epithelial_df
  )
  print("filtered_counts dimensions: ")
  print(dim(filtered_counts))

  filtered_metadata <- metadata[
    metadata$cell_ids %in% colnames(filtered_counts),
  ]

  # write as infercnv input files:
  if (!file.exists(paste0(filtered_dir, "input_matrix.txt"))) {
    write.table(filtered_counts, paste0(filtered_dir, "input_matrix.txt"), quote=F,
      sep="\t", col.names=T, row.names=T)
    write.table(filtered_metadata, paste0(filtered_dir, "metadata.txt"), sep = "\t",
    quote = F, col.names = F, row.names = F)
  }

  # save data:
  saveRDS(epithelial_df, paste0(common_Robject_dir, "/1a.epithelial_df.Rdata"))
  saveRDS(filtered_epithelial_df, paste0(common_Robject_dir, "/1b.filtered_epithelial_df.Rdata"))
  saveRDS(non_epithelial_df, paste0(common_Robject_dir, "/1c.non_epithelial_df.Rdata"))
  saveRDS(metadata, paste0(common_Robject_dir, "/1d.metadata.Rdata"))
  saveRDS(gene_annotation, paste0(common_Robject_dir, "/1e.gene_annotation.Rdata"))
  
} else {

  # load data:
  epithelial_df <- readRDS(paste0(common_Robject_dir, "/1a.epithelial_df.Rdata"))
  filtered_epithelial_df <- readRDS(paste0(common_Robject_dir, "/1b.filtered_epithelial_df.Rdata"))
  non_epithelial_df <- readRDS(paste0(common_Robject_dir, "/1c.non_epithelial_df.Rdata"))
  metadata <- readRDS(paste0(common_Robject_dir, "/1d.metadata.Rdata"))
  gene_annotation <- readRDS(paste0(common_Robject_dir, "/1e.gene_annotation.Rdata"))
  
}

if (sim_name != "filtered_normal") {

  ################################################################################
  ### 2. Prepare chromosome information ###
  ################################################################################
  
  # determine chromosome information:
  chromosome_lengths <- lapply(split(gene_annotation, gene_annotation$chromosome), 
    nrow)
  chromosome_lengths <- chromosome_lengths[
    naturalsort(names(chromosome_lengths))
  ]
  
  for (n in 1:length(chromosome_lengths)) {
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
  ### 3. Plot average fold difference from median ###
  ################################################################################
  
  # generate average original counts vector with each value representing a gene:
  average_original_counts <- apply(filtered_epithelial_df, 1, mean)
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
  
  # take mean of original fold change:
  median_original_fold_change <- median(original_fold_change)
  
  # take the log10:
  log_original_fold_change <- log10(original_fold_change)
  # tabulate for plotting:
  log_original_fold_change_df <- data.frame(
    number = seq(1, length(log_original_fold_change)),
    count = log_original_fold_change
  )
  
  # calculate the median of the counts vector:
  log_median_original_fold_change <- log10(median_original_fold_change)
  
  # plot counts:
  if (!file.exists(paste0(common_plot_dir, "1.log_original_fold_change_from_median.pdf"))) {
    p <- ggplot(log_original_fold_change_df, aes(x=number, y=count))
    p <- p + geom_point(colour = "#E8D172")
    p <- p + xlab("Genomic location")
    #p <- p + xlim(c(0, nrow(centered_original_counts_df)))
    p <- p + scale_x_continuous(
      breaks = unlist(chromosome_midpoints),
      labels = c("chr1", 2:length(chromosome_midpoints)),
      #limits = c(0,nrow(centered_original_counts_df)), 
      expand = c(0, 0)
    )
    p <- p + ylab("Log10 fold change from median")
    p <- p + scale_y_continuous(
      breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4),
      labels = c("-4", "-3", "-2", "-1", "0", "1", "2", "3", "4"),
      limits = c(-4, 4)
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
    p <- p + theme(
      text=element_text(size = 24),
      axis.text.x=element_text(size = 20),
      axis.ticks.x=element_blank(),
      axis.text.y=element_text(size = 20),
      axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0))
    )
  
  	# create full pdf and png:
    pdf(paste0(common_plot_dir, "1.log10_original_fold_change_from_median.pdf"), width = 20)
      print(p)
    dev.off()

    png(paste0(common_plot_dir, "1.log10_original_fold_change_from_median.png"), 
      width = 20,
      height = 7,
      res = 300,
      units = "in")
      print(p)
    dev.off()

    # create title-less plot:
    p <- p + theme(
      text=element_text(size = 30),
      axis.text.x=element_blank(),
      axis.text.y=element_text(size = 30),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )

    png(paste0(common_plot_dir, "1.log10_original_fold_change_from_median_minimal.png"), 
      width = 22,
      height = 7,
      res = 300,
      units = "in")
      print(p)
    dev.off()

#    # create label-less png to reduce load time:
#    p <- p + theme(
#      axis.text.x=element_blank(),
#      axis.title.x=element_blank(),
#      axis.ticks.x=element_blank(),
#      axis.text.y=element_blank(),
#      axis.title.y=element_blank(),
#      axis.ticks.y=element_blank()
#    )
#
#    png(paste0(common_plot_dir, "1.log10_original_fold_change_from_median.png"), 
#      width = 20,
#      height = 10,
#      res = 300,
#      units = "in")
#      print(p)
#    dev.off()
#
#    # create label only pdf to copy to png:
#    p <- ggplot(log_original_fold_change_df, aes(x=number, y=count))
#    p <- p + xlab("Genomic location")
#    #p <- p + xlim(c(0, nrow(centered_original_counts_df)))
#    p <- p + scale_x_continuous(
#      breaks = unlist(chromosome_midpoints),
#      labels = paste0("chr", 1:length(chromosome_midpoints)),
#      #limits = c(0,nrow(centered_original_counts_df)), 
#      expand = c(0, 0)
#    )
#    p <- p + ylab("Log10 fold change from median")
#    p <- p + scale_y_continuous(
#      breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4),
#      labels = c("-4", "-3", "-2", "-1", "0", "1", "2", "3", "4"),
#      limits = c(-4, 4)
#    )
#    for (end in chromosome_ends) {
#      p <- p + geom_vline(xintercept=end)
#    }
#    p <- p + geom_segment(
#      x=0, 
#      xend=max(unlist(chromosome_ends)), 
#      y=log_median_original_fold_change, 
#      yend=log_median_original_fold_change, 
#      size=1, color="red"
#    )
#  
#    pdf(paste0(common_plot_dir, "1.log10_original_fold_change_from_median.pdf"), width = 20)
#      print(p)
#    dev.off()
  
  }
  
  
  ################################################################################
  ### 4. Add CNVs to df, determine fold change from original median of each 
  # gene and median fold change for each CNV region ###
  ################################################################################
  
  # determine total gene length of genome to ensure total CNV length < 0.4 x 
  # genome length to comply with InferCNV assumptions:
  genome_length <- nrow(filtered_epithelial_df)
  
  if (
    !file.exists(paste0(Robject_dir, "/2a.pre_noise_simulated_epithelial_df.Rdata")) | 
    !file.exists(paste0(Robject_dir, "/2b.pre_noise_log_modified_fold_change_df.Rdata")) |
    !file.exists(paste0(Robject_dir, "/2c.CNV_record.Rdata")) | 
    !file.exists(paste0(Robject_dir, "simulated_CNV_plot_data.Rdata"))
  ) {
  
    print("Adding CNVs to dataset...")
    writeLines("\n")
  
    # choose random indices for start positions:
    set.seed(seed_record["start_position",])
    gene_no <- nrow(filtered_epithelial_df)
    random_starts <- sample(1:gene_no, 1000)
  
    # choose random indices for lengths:
    set.seed(seed_record["length",])
    random_lengths <- CNV_lengths[
      sample(1:length(CNV_lengths), 1000, replace=T)
    ]
  
     # choose random indices for multipliers:
    set.seed(seed_record["multiplier",])
    random_multipliers <- CNV_multipliers[
      sample(1:length(CNV_multipliers), CNV_no, replace = T)
    ]
  
    r=1
    for (i in 1:CNV_no) {
  
      if (exists("CNV_record")) {
        total_CNV_length <- sum(CNV_record$end-CNV_record$start)
      } else {
        total_CNV_length <- 0
      }
  
  
      if ( total_CNV_length < (genome_length*0.4) ) {
  
        print(paste0(
          "Total CNV length still less than genome length x 0.4, adding another CNV..."
        ))
        writeLines("\n")
  
        repeat {
          # choose start position at random:
          start_position <- random_starts[r]
        
          # choose CNV length at random:
          CNV_length <- random_lengths[r]
          end_position <- start_position+CNV_length
          CNV_region <- start_position:end_position
        
          # check this region does not overlap end of genome or any prior CNVs and if not, break loop:
          if (end_position < nrow(filtered_epithelial_df)) {
            if (exists("CNV_record")) {
      
              print("CNV_record exists, checking region does not overlap with any prior CNV regions...")
              
              no_break_indicator <- FALSE
      
              for (n in 1:nrow(CNV_record) ) {
                prior_region <- CNV_record$start[n]:CNV_record$end[n]
                if ( length(Reduce(intersect, list(prior_region, CNV_region))) > 0 ) {
                  no_break_indicator <- TRUE
                }
              }
              if (no_break_indicator) {
                print("Overlaps found between prior CNV region/s, trying different region...")
              } else {
                print("No overlaps found, keeping region...")
                writeLines("\n")
                break
              }
        
            } else {
              print("CNV_record does not exist, keeping region...")
              break
              writeLines("\n")
            }
    
          } else {
            print("Region overlaps with end of genome, trying different region...")
          }
  
          r <<- r+1
        }
      
        # choose CNV multiplier at random:
        CNV_multiplier <- random_multipliers[i]
    
        # record region so no future CNVs overlap:
        if (!(exists("CNV_record"))) {
          CNV_record <- data.frame(
            start = start_position,
            end = end_position,
            multiplier = CNV_multiplier,
            log_median_modified_FC = NA
          )
        } else {
          CNV_record <- rbind(
            CNV_record, 
            data.frame(
              start = start_position,
              end = end_position,
              multiplier = CNV_multiplier,
              log_median_modified_FC = NA
            )
          )
        }
      
        # add CNV to modified_df if exists, or original filtered_epithelial_df if not:
        if (!exists("modified_df")) {
          modified_df <- filtered_epithelial_df
        }
        modified_df[CNV_record$start[i]:CNV_record$end[i],] <- 
          modified_df[CNV_record$start[i]:CNV_record$end[i],]*CNV_multiplier
  
        # generate mean original segment fold change vector with each value representing 
        # a gene:
        average_original_counts <- apply(
          filtered_epithelial_df[CNV_record$start[i]:CNV_record$end[i],], 1, mean
        )
        # add 0.1 to all values:
        #average_original_counts[average_original_counts == 0] <- 1e-3
        average_original_counts <- average_original_counts + 1e-4
        # determine median:
        median_average_original_counts <- median(average_original_counts)
        # divide by median to get fold change from median:
        original_fold_change <- average_original_counts/median_average_original_counts
        # check median of original fold change = 1:
        median_original_fold_change <- median(original_fold_change)
    
        # generate mean modified fold change from original median vector with each value 
        # representing a gene:
        average_modified_counts <- apply(
          modified_df[CNV_record$start[i]:CNV_record$end[i],], 1, mean
        )
        # add 0.1 to all values to account for zeros:
        average_modified_counts <- average_modified_counts + 1e-4
        # divide by median to get fold change from original median and add to df for plotting:
        modified_fold_change <- average_modified_counts/median_average_original_counts
        # take log of the median of modified fold change to mark on plot:
        median_modified_fold_change <- median(modified_fold_change)    
        log_median_modified_fold_change <- log10(median_modified_fold_change)
    
        # add to CNV_record:
        CNV_record$log_median_modified_FC[i] <- log_median_modified_fold_change

        ###
#        if (is.na(CNV_record$log_median_modified_FC[i])) {
#          print(i)
#          print("NA found")
#          break()
#        }
        ###

        # take the log10 of all fold changes:
        log_modified_fold_change <- log10(modified_fold_change)
    
        # add modified fold change values to log_original_fold_change_df:
        if (!exists("log_modified_fold_change_df")) {
          log_modified_fold_change_df <- log_original_fold_change_df
        }
        log_modified_fold_change_df$count[
          CNV_record$start[i]:CNV_record$end[i]
        ] <- log_modified_fold_change
  
      } else {
        print(paste0(
          "Total CNV length is too close to genome length x 0.4, no CNVs added..."
        ))
        writeLines("\n")
      }
    }
  
    # add non-CNV regions to CNV record:
    CNV_record_gr <- GRanges(
      seqnames = Rle("genome"),
      ranges = IRanges(start = CNV_record$start, end = CNV_record$end),
      strand = Rle("*"),
      multiplier = CNV_record$multiplier,
      log_median_modified_FC = CNV_record$log_median_modified_FC
    )
    # identify gaps between marked CNV regions as non-CNV regions and fill in values accordingly:
    non_CNV_record <- gaps(CNV_record_gr)
    non_CNV_record$multiplier = rep(1, length(non_CNV_record))
    non_CNV_record$log_median_modified_FC = rep(0, length(non_CNV_record))
    CNV_record_gr <- c(CNV_record_gr, non_CNV_record)
    # convert back to data frame and order:
    final_CNV_record <- data.frame(
      start = start(ranges(CNV_record_gr)),
      end = end(ranges(CNV_record_gr)),
      multiplier = CNV_record_gr$multiplier,
      log_median_modified_FC = CNV_record_gr$log_median_modified_FC
    )
    final_CNV_record <- final_CNV_record[order(final_CNV_record$start),]
    # fill in last non-CNV segment:
    final_CNV_record <- rbind(
      final_CNV_record,
      data.frame(
        start = final_CNV_record$end[nrow(final_CNV_record)]+1,
        end = nrow(filtered_epithelial_df),
        multiplier = 1,
        log_median_modified_FC = 0
      )
    )
    
    # add chromosome information:
    final_CNV_record$start_chr <- "chr1"
    final_CNV_record$end_chr <- "chr1"
    for (k in 1:length(chromosome_ends)) {
      if (k==1) {
  
        final_CNV_record$start_chr[
          final_CNV_record$start <= unlist(chromosome_ends[k])
        ] <- names(chromosome_ends)[k]
  
        final_CNV_record$end_chr[
          final_CNV_record$end <= unlist(chromosome_ends[k])
        ] <- names(chromosome_ends)[k]
  
      } else {
  
        final_CNV_record$start_chr[
          final_CNV_record$start <= unlist(chromosome_ends[k]) & 
          final_CNV_record$start > unlist(chromosome_ends[k-1])
        ] <- names(chromosome_ends)[k]
  
        final_CNV_record$end_chr[
          final_CNV_record$end <= unlist(chromosome_ends[k]) & 
          final_CNV_record$end > unlist(chromosome_ends[k-1])
        ] <- names(chromosome_ends)[k]
  
      }
    }
    
    saveRDS(modified_df, paste0(Robject_dir, "/2a.pre_noise_simulated_epithelial_df.Rdata"))
    saveRDS(
      log_modified_fold_change_df, 
      paste0(Robject_dir, "/2b.pre_noise_log_modified_fold_change_df.Rdata")
    )
    saveRDS(final_CNV_record, paste0(Robject_dir, "/2c.CNV_record.Rdata"))
  
    simulated_CNV_plot_data <- list(
      log_modified_fold_change_df,
      final_CNV_record
    )
    names(simulated_CNV_plot_data) <- c("log_modified_fold_change_df", 
      "CNV_indices")
    saveRDS(simulated_CNV_plot_data, 
      paste0(Robject_dir, "simulated_CNV_plot_data.Rdata")
    )
  
  } else {
  
    modified_df <- readRDS(paste0(Robject_dir, 
      "/2a.pre_noise_simulated_epithelial_df.Rdata"))
    log_modified_fold_change_df <- readRDS(paste0(Robject_dir, 
      "/2b.pre_noise_log_modified_fold_change_df.Rdata"))
    final_CNV_record <- readRDS(paste0(Robject_dir, "/2c.CNV_record.Rdata"))
    simulated_CNV_plot_data <- readRDS(paste0(Robject_dir, "simulated_CNV_plot_data.Rdata"))
  
  }
  
  
  ################################################################################
  ### 5. Create modified counts and line CNV plots ###
  ################################################################################
  
  # plot median fold change from original median for modified data:
  if (!file.exists(paste0(plot_dir, "2a.pre_noise_log_modified_fold_change_from_median_line_only.pdf"))) {
    p <- ggplot(log_modified_fold_change_df, aes(x=number, y=count))
    p <- p + scale_x_continuous(
      breaks = unlist(chromosome_midpoints),
      labels = c("chr1", 2:length(chromosome_midpoints)),
      limits = c(0,length(log_modified_fold_change_df$count)), 
      expand = c(0, 0)
    )
    p <- p + scale_y_continuous(
      breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4),
      labels = c("-4", "-3", "-2", "-1", "0", "1", "2", "3", "4"),
      limits = c(
        min(log_modified_fold_change_df$count), 
        max(log_modified_fold_change_df$count)
      )
    )
    for (end in chromosome_ends) {
      p <- p + geom_vline(xintercept=end)
    }
    for (r in 1:nrow(final_CNV_record)) {
      # create horizontal line:
      p <- p + geom_segment(
        x=final_CNV_record$start[r], 
        xend=final_CNV_record$end[r], 
        y=final_CNV_record$log_median_modified_FC[r], 
        yend=final_CNV_record$log_median_modified_FC[r], 
        size=1, color="#37841f"
      )
      # create left vertical line:
      if (r != 1) {
        p <- p + geom_segment(
          x=final_CNV_record$start[r], 
          xend=final_CNV_record$start[r], 
          y=final_CNV_record$log_median_modified_FC[r-1], 
          yend=final_CNV_record$log_median_modified_FC[r], 
          size=1, color="#37841f"
        )
      }
    }
    p <- p + theme(
      text=element_text(size = 24),
      axis.text.x=element_text(size = 20),
      axis.text.y=element_text(size = 20),
      axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0))
    )

    pdf(paste0(plot_dir, 
      "2a.pre_noise_log_modified_fold_change_from_median_line_only.pdf"), width = 20)
      print(p)
    dev.off()
    
  }
  
  # plot counts:
  if (!file.exists(paste0(plot_dir, "2b.pre_noise_log_modified_fold_change_from_median.pdf"))) {
    p <- ggplot(log_modified_fold_change_df, aes(x=number, y=count))
    p <- p + geom_point(colour = "#E8D172")
    p <- p + xlab("Genomic location")
    p <- p + scale_x_continuous(
      breaks = unlist(chromosome_midpoints),
      labels = c("chr1", 2:length(chromosome_midpoints)),
      limits = c(0,nrow(log_modified_fold_change_df)), 
      expand = c(0, 0)
    )
    p <- p + ylab("Log10 fold change")
    p <- p + scale_y_continuous(
      breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4),
      labels = c("-4", "-3", "-2", "-1", "0", "1", "2", "3", "4"),
      limits = c(
          min(log_modified_fold_change_df$count), 
          max(log_modified_fold_change_df$count)
        )
      )
    for (end in chromosome_ends) {
      p <- p + geom_vline(xintercept=end)
    }
    for (r in 1:nrow(final_CNV_record)) {
      # create horizontal line:
      p <- p + geom_segment(
        x=final_CNV_record$start[r], 
        xend=final_CNV_record$end[r], 
        y=final_CNV_record$log_median_modified_FC[r], 
        yend=final_CNV_record$log_median_modified_FC[r], 
        size=1, color="red"
      )
      # create left vertical line:
      if (r != 1) {
        p <- p + geom_segment(
          x=final_CNV_record$start[r], 
          xend=final_CNV_record$start[r], 
          y=final_CNV_record$log_median_modified_FC[r-1], 
          yend=final_CNV_record$log_median_modified_FC[r], 
          size=1, color="red"
        )
      }
    }
    p <- p + theme(
      text=element_text(size = 24),
      axis.text.x=element_text(size = 20),
      axis.text.y=element_text(size = 20),
      axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0))
    )

    pdf(paste0(plot_dir, "2b.pre_noise_log_modified_fold_change_from_median.pdf"), width = 20)
      print(p)
    dev.off()
  }
  
  
  ################################################################################
  ### 6. Bind new epithelial and original stromal counts and check against 
  # metadata ###
  ################################################################################
  
  # bind non_epithelial_df and modified_df together as new_counts:
  print("modified_df dimensions:")
  print(dim(modified_df))
  
  new_counts <- cbind(
    modified_df,
    non_epithelial_df
  )
  print("New_counts dimensions: ")
  print(dim(new_counts))
  
  
  ###################################################################################
  ### 7. Simulate noise dataset from noise counts ###
  ###################################################################################
  
  if (!file.exists(paste0(noise_dir, "/noise_df.Rdata"))) {
  
    all_cells <- read.csv(paste0(emptydrops_dir, "01_Emptydrops_out.csv"))
    filtered_cells <- read.csv(
      paste0(emptydrops_dir, "02_emptydrops_filtered_cell_ids.csv")
    )
    raw_10X <- readRDS(paste0(seurat_dir, "01_Read10X_raw_data.Rdata"))
    seurat_10X <- readRDS(paste0(seurat_dir, "03_seurat_object_processed.Rdata"))
    
    # isolate outfiltered cells from raw counts matrix:
    outfiltered_cells <- all_cells[!(all_cells$X %in% filtered_cells$cell_ids),]
    cells_to_select <- colnames(raw_10X)[
      which(colnames(raw_10X) %in% outfiltered_cells$X)
    ]
    outfiltered_counts <- raw_10X[,colnames(raw_10X) %in% cells_to_select]
    
    # sanity checks:
    print(paste0(
      "Total number of cells in EmptyDrops out object = ", nrow(all_cells)
    ))
    print(paste0(
      "Total number of cells in raw matrix object = ", ncol(raw_10X)
    ))
    print(paste0(
      "Number of cells in filtered EmptyDrops object = ", nrow(filtered_cells)
    ))
    print(paste0(
      "Number of cells in seurat object = ", ncol(GetAssayData(seurat_10X , slot = "counts"))
    ))
    print(paste0(
      "Does total cell number - number of outfiltered cells = number of filtered cells? ", 
      nrow(all_cells) - nrow(outfiltered_cells) == nrow(filtered_cells)
    ))
    print(paste0(
      "Number of cells in outfiltered_counts = ", ncol(outfiltered_counts)
    ))
  
    # order by total counts:
    column_sums <- Matrix::colSums(outfiltered_counts)
    outfiltered_counts <- outfiltered_counts[
      ,order(column_sums, decreasing = T)
    ]
  
    # plot distribution of outfiltered counts:
    count_sums <- colSums(outfiltered_counts)
    outfiltered_counts <- outfiltered_counts[,count_sums != 0]
    outfiltered_counts <- outfiltered_counts[,1:noise_cell_no]
    outfiltered_means <- apply(outfiltered_counts, 2, mean)
    if (!file.exists(paste0(noise_dir, "outfiltered_count_density_plot.pdf"))) {
      outfiltered_count_density_plot <- density(as.vector(outfiltered_means))
      pdf(paste0(noise_dir, "outfiltered_count_density_plot.pdf"))
        plot(outfiltered_count_density_plot, main=NA, xlab = "nUMI")
      dev.off()
    }
  
    # select for cells to use for simulated noise dataset:
    noise_input <- outfiltered_counts[,1:noise_cell_no]
   
    if (!file.exists(paste0(noise_dir, "/noise_params.Rdata"))) {
  
      # calculate params from input to use for simulation:
      splat_params <- splatEstimate(as.matrix(noise_input))
      # reduce simulation to number of cells in simulation:
      splat_params <- setParam(splat_params, "nGenes", nrow(new_counts))
  
      saveRDS(splat_params, paste0(noise_dir, "/noise_params.Rdata"))
  
    } else {
  
      splat_params <- readRDS(paste0(noise_dir, "/noise_params.Rdata"))
  
    }
  
    # simulate the data:
    sim <- splatSimulate(splat_params, batchCells = ncol(new_counts))
    noise_counts <- counts(sim)
    rownames(noise_counts) <- rownames(new_counts)
   
    print("Table of simulated noise counts:")
    table(unlist(noise_counts))

    plot_noise <- function(
      noise_df,
      ref_dir,
      plot_dir
    ) {

      # find chromosome ends to mark on plot:
      # load infercnv gene annotation:
      gene_annot <- read.table(
        paste0(ref_dir, "infercnv_gene_order.txt"),
        header = F
      )
      colnames(gene_annot) <- c("gene", "chr", "start", "end")
  
      # calculate mean of each input gene:
      input_means <- apply(noise_df, 1, mean)
      
      # only keep input_means genes in gene_annot and vice versa:
      gene_annot <- gene_annot[gene_annot$gene %in% names(input_means),]
      input_means <- input_means[names(input_means) %in% gene_annot$gene]
  
      # reorder input_means:
      input_means <- input_means[as.character(gene_annot$gene)]
  
      # remove chrX, Y:
      to_remove <- gene_annot$gene[
        gene_annot$chr == "chrX" | gene_annot$chr == "chrY" | 
        gene_annot$chr == "chrM"
      ]
      gene_annot <- gene_annot[!(gene_annot$gene %in% to_remove), ]
      input_means <- input_means[!(names(input_means) %in% to_remove)]
  
      # split annotation and find end indices:
      split_annot <- split(gene_annot, gene_annot$chr)
      end_genes <- unlist(
        lapply(split_annot, function(x) {
          return(tail(as.character(x$gene), 1))
        })
      )
      # determine chr end indices:
      end_indices <- sort(match(end_genes, names(input_means)))
  
      # determine middle gene for each chromosome:
      mid_genes <- unlist(
         lapply(split_annot, function(x) {
           return(as.character(x$gene)[floor(nrow(x)/2)])
         })
       )
      # determine chr mid indices:
      mid_indices <- sort(match(mid_genes, names(input_means)))
  
      input_df <- data.frame(
        number = 1:length(input_means),
        input = input_means
      )
  
      # log 10:
      input_df$log_input <- log10(input_df$input+1e-4)
  
      # split input_means into n gene bins:
      bin_no <- 500
      div_no <- floor(length(input_df$log_input)/bin_no)
      split_vec <- c(
        rep(1:div_no, each=bin_no),
        rep(div_no+1, length(input_df$log_input) - (div_no*bin_no))
      )
      split_means <- split(input_df$log_input, split_vec)
  
      # calcuate mean of each bin and log10:
      bin_means <-  unlist(lapply(split_means, mean))
      bin_medians <-  unlist(lapply(split_means, median))
      
      p1 <- ggplot(input_df, aes(x=number, y=log_input))
      p1 <- p1 + geom_point(colour="#E7E4D3")
      p1 <- p1 + xlab("Genomic location")
      p1 <- p1 + ylab("log10 mean count")
  
      for (ind in end_indices) {
        
        # create chromosome end vertical line:
        p1 <- p1 + geom_segment(
          x=ind,
          xend=ind,
          y=min(input_df$log_input), 
          yend=max(input_df$log_input), 
          size=0.5, color="black"
        )
    
      }
  
      # mark bin means:
      for (b in 1:length(bin_means)) {
        if (b==length(bin_means)) {
          p1 <- p1 + geom_segment(
            x=((b-1)*bin_no)+1,
            xend=length(input_df$log_input),
            y=bin_means[b], 
            yend=bin_means[b], 
            size=0.5, color="red"
          )
        } else {
          p1 <- p1 + geom_segment(
            x=((b-1)*bin_no)+1,
            xend=b*bin_no,
            y=bin_means[b], 
            yend=bin_means[b], 
            size=0.5, color="red"
          )
        }
      }
  
      # mark bin vertical lines:
      for (b in 1:length(bin_means)) {
        if (b!=length(bin_means)) {
          p1 <- p1 + geom_segment(
            x=b*bin_no,
            xend=b*bin_no,
            y=bin_means[b], 
            yend=bin_means[b+1], 
            size=0.5, color="red"
          )
        }
      }
    
      # create chromosome labels:
      p1 <- p1 + scale_x_continuous(
        breaks = mid_indices,
        labels = c("chr1", as.character(2:length(input_mid_ind))),
        expand = c(0, 0)
      )

    }

    # plot mean noise input across the genome:  
    noise_input_plot <- plot_noise(
      noise_input,
      ref_dir,
      plot_dir
    )
  
    png(
      paste0(noise_dir, "log10_noise_input_scatterplot.png"),
      width = 20,
      height = 7,
      res = 300,
      unit = "in"
    )
      noise_input_plot
    dev.off()
  
    # plot mean simulated noise across the genome:  
    sim_noise_plot <- plot_noise(
      noise_counts,
      ref_dir,
      plot_dir
    )
  
    png(
      paste0(noise_dir, "log10_simulated_noise_scatterplot.png"),
      width = 20,
      height = 7,
      res = 300,
      unit = "in"
    )
      sim_noise_plot
    dev.off()

    saveRDS(noise_counts, paste0(noise_dir, "/noise_df.Rdata"))
  
  } else {
    noise_counts <- readRDS(paste0(noise_dir, "/noise_df.Rdata"))
  }
  
  
  ###################################################################################
  ### 8. Add noise to simulated counts ###
  ###################################################################################
  
  # ensure metadata has same cells as modified_df and label stromal cells:
  epithelial_metadata <- metadata[
    grep("Epithelial", metadata$cell_types),
  ]
  epithelial_metadata <- epithelial_metadata[
    epithelial_metadata$cell_ids %in% colnames(modified_df),
  ]
  
  non_epithelial_metadata <- metadata[
    grep("Non_epithelial", metadata$cell_types),
  ]
  
  metadata_df <- rbind(
    non_epithelial_metadata,
    epithelial_metadata
  )
  
  # add noise to counts:
  without_noise <- new_counts
  print("E.g. summary of new counts before noise added:") 
  print(summary(without_noise[,1]))
  with_noise <- without_noise + noise_counts
  print("E.g. summary of new counts after noise added:") 
  print(summary(with_noise[,1]))
  
  # create updated CNV scatterplots post-noise addition:
  # generate mean original segment fold change vector with each value representing 
  # a gene:
  final_epithelial_df <- with_noise[
    ,colnames(with_noise) %in% epithelial_metadata$cell_ids
  ]
  
  final_CNV_record$log_median_modified_FC_post_noise <- NA
  for (i in 1:nrow(final_CNV_record)) {
 
    if (final_CNV_record$start[i] != final_CNV_record$end[i]) {
      average_original_counts <- apply(
        filtered_epithelial_df[final_CNV_record$start[i]:final_CNV_record$end[i],], 1, mean)
    } else {
      average_original_counts <- mean(as.numeric(filtered_epithelial_df[final_CNV_record$start[i],]))
    }
  
    # add 0.1 to all values:
    #average_original_counts[average_original_counts == 0] <- 1e-3
    average_original_counts <- average_original_counts + 1e-3
    # determine median:
    median_average_original_counts <- median(average_original_counts)
    # divide by median to get fold change from median:
    original_fold_change <- average_original_counts/median_average_original_counts
    # check median of original fold change = 1:
    median_original_fold_change <- median(original_fold_change)
    # generate mean modified fold change from original median vector with each value 
    # representing a gene:
    if (final_CNV_record$start[i] != final_CNV_record$end[i]) {
      average_modified_counts <- apply(
        final_epithelial_df[final_CNV_record$start[i]:final_CNV_record$end[i],], 1, mean
      )
    } else {
      average_modified_counts <- mean(as.numeric(final_epithelial_df[final_CNV_record$start[i],]))
    }
    # add 0.1 to zero values:
    average_modified_counts <- average_modified_counts + 1e-3
    # divide by median to get fold change from original median and add to df for plotting:
    modified_fold_change <- average_modified_counts/median_average_original_counts
    # take log of the median of modified fold change to mark on plot:
    median_modified_fold_change <- median(modified_fold_change)    
    log_median_modified_fold_change <- log10(median_modified_fold_change)
    # add to final_CNV_record:
    final_CNV_record$log_median_modified_FC_post_noise[i] <- log_median_modified_fold_change
    # take the log10 of all fold changes:
    log_modified_fold_change <- log10(modified_fold_change)
    # add modified fold change values to log_original_fold_change_df:
    if (!exists("log_modified_fold_change_df")) {
      log_modified_fold_change_df <- log_original_fold_change_df
    }
    log_modified_fold_change_df$count[
      final_CNV_record$start[i]:final_CNV_record$end[i]
    ] <- log_modified_fold_change

  }

  
  ###################################################################################
  ### 9. Plot post-noise gene expression profiles ###
  ###################################################################################
  
  # plot median fold change from original median for modified data:
  if (!file.exists(paste0(plot_dir, "3a.log_modified_fold_change_from_median_with_noise_line_only.pdf"))) {
    p <- ggplot(log_modified_fold_change_df, aes(x=number, y=count))
    p <- p + scale_x_continuous(
      breaks = unlist(chromosome_midpoints),
      labels = c("chr1", 2:length(chromosome_midpoints)),
      limits = c(0,length(log_modified_fold_change_df$count)), 
      expand = c(0, 0)
    )
    p <- p + scale_y_continuous(
        breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4),
        labels = c("-4", "-3", "-2", "-1", "0", "1", "2", "3", "4"),
        limits = c(
          min(log_modified_fold_change_df$count), 
          max(log_modified_fold_change_df$count)
        )
      )
    for (end in chromosome_ends) {
      p <- p + geom_vline(xintercept=end)
    }
    for (r in 1:nrow(final_CNV_record)) {
      # create horizontal line:
      p <- p + geom_segment(
        x=final_CNV_record$start[r], 
        xend=final_CNV_record$end[r], 
        y=final_CNV_record$log_median_modified_FC_post_noise[r], 
        yend=final_CNV_record$log_median_modified_FC_post_noise[r], 
        size=1, color="#37841f"
      )
      # create left vertical line:
      if (r != 1) {
        p <- p + geom_segment(
          x=final_CNV_record$start[r], 
          xend=final_CNV_record$start[r], 
          y=final_CNV_record$log_median_modified_FC_post_noise[r-1], 
          yend=final_CNV_record$log_median_modified_FC_post_noise[r], 
          size=1, color="#37841f"
        )
      }
    }
    p <- p + theme(
      text=element_text(size = 24),
      axis.text.x=element_text(size = 20),
      axis.text.y=element_text(size = 20),
      axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0))
    )

    pdf(paste0(plot_dir, 
      "3a.log_modified_fold_change_from_median_with_noise_line_only.pdf"), width = 20)
      print(p)
    dev.off()
    
  }
  
  # plot counts:
  if (!file.exists(paste0(plot_dir, "3b.log_modified_fold_change_from_median_with_noise.pdf"))) {
    p <- ggplot(log_modified_fold_change_df, aes(x=number, y=count))
    p <- p + geom_point(colour = "#E8D172")
    p <- p + xlab("Genomic location")
    p <- p + scale_x_continuous(
      breaks = unlist(chromosome_midpoints),
      labels = c("chr1", 2:length(chromosome_midpoints)),
      limits = c(0,nrow(log_modified_fold_change_df)), 
      expand = c(0, 0)
    )
    p <- p + ylab("Log10 fold change")
    p <- p + scale_y_continuous(
       breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4),
       labels = c("-4", "-3", "-2", "-1", "0", "1", "2", "3", "4"),
       limits = c(min(log_modified_fold_change_df$count), max(log_modified_fold_change_df$count))
     )
    for (end in chromosome_ends) {
      p <- p + geom_vline(xintercept=end)
    }
    for (r in 1:nrow(final_CNV_record)) {
      # create horizontal line:
      p <- p + geom_segment(
        x=final_CNV_record$start[r], 
        xend=final_CNV_record$end[r], 
        y=final_CNV_record$log_median_modified_FC_post_noise[r], 
        yend=final_CNV_record$log_median_modified_FC_post_noise[r], 
        size=1, color="red"
      )
      # create left vertical line:
      if (r != 1) {
        p <- p + geom_segment(
          x=final_CNV_record$start[r], 
          xend=final_CNV_record$start[r], 
          y=final_CNV_record$log_median_modified_FC_post_noise[r-1], 
          yend=final_CNV_record$log_median_modified_FC_post_noise[r], 
          size=1, color="red"
        )
      }
    }
    p <- p + theme(
      text=element_text(size = 30),
      axis.text.x=element_text(size = 30),
      axis.text.y=element_text(size = 30),
      axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0))
    )

    pdf(paste0(plot_dir, "3b.log_modified_fold_change_from_median_with_noise.pdf"), width = 20)
      print(p)
    dev.off()

    p <- p + theme(
      text=element_text(size = 30),
      axis.text.x=element_text(size = 30),
      axis.text.y=element_text(size = 30),
      axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
      axis.title.y = element_blank()
    )

    png(paste0(plot_dir, "3b.log_modified_fold_change_from_median_with_noise_minimal.png"), 
      width = 22,
      height = 7,
      res = 300,
      units = "in")
      print(p)
    dev.off()

  }
  
  if (!file.exists(paste0(out_dir, "input_matrix.txt"))) {
    write.table(with_noise, paste0(out_dir, "input_matrix.txt"), quote=F,
      sep="\t", col.names=T, row.names=T)
    write.table(metadata_df, paste0(out_dir, "metadata.txt"), sep = "\t",
    quote = F, col.names = F, row.names = F)
  }

  print("Simulation infercnv output created.")

} else {
  print("Filtered normal infercnv output created.")
}




