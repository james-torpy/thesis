#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

project_name <- "thesis"
subproject_name <- "Figure_2.2_individual_samples"
args = commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
subcluster_method <- args[2]
subcluster_p <- args[3]
if (subcluster_p != "none") {
  subcluster_p <- as.numeric(subcluster_p)
}
coverage_filter <- args[4]
remove_artefacts <- args[5]
subset_data <- as.logical(args[6])
subset_samples <- as.logical(args[7])
na_colour <- args[8]
gene_proportion_threshold <- as.numeric(args[9])
QC_annot <- as.logical(args[10])
subset_samples <- as.logical(args[11])

project_name <- "thesis"
subproject_name <- "Figure_2.2_individual_samples"
args = commandArgs(trailingOnly=TRUE)
subcluster_method <- "random_trees"
subcluster_p <- "0.05"
if (subcluster_p != "none") {
  subcluster_p <- as.numeric(subcluster_p)
}
coverage_filter <- "filtered"
remove_artefacts <- "artefacts_not_removed"
subset_samples <- FALSE

if (subset_samples) {

  # determine order of samples by subtype:
  ER <- c("CID3941")
  HER2 <- c("CID3586")
  TNBC <- c("CID44041")
  sample_names <- c(ER, HER2, TNBC)

  subtype_df <- data.frame(
    subtype = c(
      rep("ER", length(ER)),
      rep("HER2", length(HER2)),
      rep("TNBC", length(TNBC))
    ),
    sample = sample_names
  )

} else {

  # determine order of samples by subtype:
  ER <- c("CID3941", "CID3948", "CID4067", "CID4290A", "CID4398", "CID4461", 
  	"CID4463", "CID4530N", "CID4535")
  HER2 <- c("CID3921", "CID3586", "CID3963", "CID4066", "CID45171")
  TNBC <- c("CID44041", "CID4465", "CID44971", "CID44991", "CID4513",
  	"CID4515", "CID4523")
  sample_names <- c(ER, HER2, TNBC)

  subtype_df <- data.frame(
    subtype = c(
      rep("ER", length(ER)),
      rep("HER2", length(HER2)),
      rep("TNBC", length(TNBC))
    ),
    sample = sample_names
  )

}

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(data.table)
library(ComplexHeatmap, lib.loc = lib_loc)
library(grid)
library(ggplot2)
library(cowplot)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/",
  subproject_name, "/")
ref_dir <- paste0(project_dir, "/refs/")
func_dir <- paste0(project_dir, "/scripts/functions/")
results_dir <- paste0(project_dir, "results/")
in_path <- paste0(results_dir, "infercnv/")
out_dir <- paste0(results_dir, "infercnv/combined_CNV_stats/", 
  coverage_filter, "/", subcluster_method, "/p_", subcluster_p, "/", 
  remove_artefacts, "/")

if (subset_samples) {
  Robject_dir <- paste0(out_dir, "Rdata_sub/")
  plot_dir <- paste0(out_dir, "plots_sub/")
  table_dir <- paste0(out_dir, "tables_sub/")
} else {
  Robject_dir <- paste0(out_dir, "Rdata/")
  plot_dir <- paste0(out_dir, "plots/")
  table_dir <- paste0(out_dir, "tables/")
}

system(paste0("mkdir -p ", Robject_dir))
system(paste0("mkdir -p ", plot_dir))
system(paste0("mkdir -p ", table_dir))

print(paste0("In path = ", in_path))
print(paste0("Out dir = ", out_dir))
print(paste0("Reference directory = ", ref_dir))
print(paste0("R function directory = ", func_dir))
print(paste0("R object directory = ", Robject_dir))
print(paste0("Plot directory = ", plot_dir))


################################################################################
### 0. Load functions and colour palette ###
################################################################################

plot_densities <- dget(
  paste0(func_dir, "plot_densities.R")
)

subtype_cols <- read.table(
  paste0(ref_dir, "subtype_colour_palette.txt"),
  header = F,
  comment.char = "",
  stringsAsFactors = F
)
m <- match(c("ER", "HER2", "TNBC"), subtype_cols$V1)
subtype_cols <- subtype_cols$V2[m]
names(subtype_cols) <- c("ER", "HER2", "TNBC")


################################################################################
### 1. Load CNV data ###
################################################################################

if (!file.exists(paste0(Robject_dir, "group_CNV_data.Rdata"))) {

  for (s in 1:length(sample_names)) {

    print(paste0("Loading ", sample_names[s], " CNV data..."))

    sample_Robject_dir <- paste0(
  	  in_path, sample_names[s], "/", coverage_filter, "/", subcluster_method, 
      "/p_", subcluster_p, "/", remove_artefacts, 
      "/Rdata/"
    )

    CNV_data <- readRDS(
  		paste0(sample_Robject_dir, "CNV_indices_and_lengths.Rdata")
	)
  
    if (s==1) {
      all_CNV_data <- list(CNV_data)
      names(all_CNV_data) <- sample_names[s]
    } else {
      all_CNV_data[[s]] <- CNV_data
      names(all_CNV_data)[s] <- sample_names[s]
    }

  }

  saveRDS(all_CNV_data, paste0(Robject_dir, "group_CNV_data.Rdata"))

} else {

  all_CNV_data <- readRDS(paste0(Robject_dir, "group_CNV_data.Rdata"))
  
}


################################################################################
### 2. Fetch subpop numbers, CNV numbers and lengths ###
################################################################################

subpop_no <- unlist(
  lapply(all_CNV_data, function(x) {
    return(length(unique(names(x))))
  })
)

# count CNVs in each subcluster:
CNV_counts <- unlist(
  lapply(all_CNV_data, function(x) {
    lapply(x, function(y) {
      nrow(y)
    })
  })
)
CNV_mean_counts <- unlist(
  lapply(all_CNV_data, function(x) {
    round(mean(unlist(lapply(x, function(y) {
      nrow(y)
    }))), 0)
  })
)

# record CNV lengths in each dataset:
CNV_lengths <- lapply(all_CNV_data, function(x) {
  all_subs <- do.call("rbind", x)
  return(all_subs$length)
})

CNV_genomic_lengths <- lapply(all_CNV_data, function(x) {
  all_subs <- do.call("rbind", x)
  return(all_subs$genomic_length)
})

# collate data and save:
CNV_dist_data <- data.frame(
  sample_name = factor(sample_names, levels = sample_names),
  subtype = subtype_df$subtype,
  no_subpops = subpop_no,
  mean_CNV_count = CNV_mean_counts,
  mean_CNV_length = round(unlist(lapply(CNV_lengths, mean)), 0),
  mean_CNV_genomic_length = round(unlist(lapply(CNV_genomic_lengths, mean)), 0)
)

write.table(
  CNV_dist_data,
  paste0(table_dir, "CNV_metadata.txt"),
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)


################################################################################
### 3. Plot data ###
################################################################################

plot_densities(
  CNV_counts,
  "counts",
  subtype_df,
  plot_dir
)

mean(CNV_counts)

# collate all lengths and rename each as sample of origin:
length_vec <- unlist(CNV_lengths)
for (l in 1:length(CNV_lengths)) {
  if (l==1) {
    vec_names <- rep(names(CNV_lengths)[l], length(CNV_lengths[[l]]))
  } else {
    vec_names <- c(
      vec_names, 
      rep(names(CNV_lengths)[l], length(CNV_lengths[[l]]))
    )
  }
}
names(length_vec) <- vec_names

plot_densities(
  length_vec,
  "lengths",
  subtype_df,
  plot_dir
)

mean(length_vec)

# collate all lengths and rename each as sample of origin:
genomic_length_vec <- unlist(CNV_genomic_lengths)
for (l in 1:length(CNV_genomic_lengths)) {
  if (l==1) {
    vec_names <- rep(names(CNV_genomic_lengths)[l], length(CNV_genomic_lengths[[l]]))
  } else {
    vec_names <- c(
      vec_names, 
      rep(names(CNV_genomic_lengths)[l], length(CNV_genomic_lengths[[l]]))
    )
  }
}
names(genomic_length_vec) <- vec_names

plot_densities(
  genomic_length_vec,
  "genomic_lengths",
  subtype_df,
  plot_dir
)

mean(genomic_length_vec)

p <- ggplot(CNV_dist_data, aes(x=sample_name, y=no_subpops, fill=subtype))
p <- p + geom_bar(stat = "identity")
p <- p + xlab("Subtype")
p <- p + ylab("No. subpopulations")
p <- p + theme_cowplot(12)
p <- p + scale_fill_manual(values = subtype_cols)
#p <- p + scale_fill_manual(values = cols, labels = c("LumA", "LumB", "HER2", "TNBC", "Metaplastic"))
p <- p + labs(fill = "Subtype")
p <- p + theme(
  axis.text.x = element_text(angle = 45, hjust = 1)
)

pdf(
  paste0(plot_dir, "subcluster_numbers.pdf"),
  height = 5,
  width = 7
)
  p
dev.off()

#convert pdf to png:
system(paste0("for p in ", plot_dir, "*.pdf; do echo $p; f=$(basename $p); echo $f; ",
"new=$(echo $f | sed 's/.pdf/.png/'); echo $new; ", 
"convert -density 150 ", plot_dir, "$f -quality 90 ", plot_dir, "$new; done"))


################################################################################
### 4. Determine significant differences ###
################################################################################

# split subtype no by subtype
split_subpop_no <- split(subpop_no, as.character(subtype_df$subtype))
names(split_subpop_no) <- c("ER", "HER2", "TNBC")

# use one-tailed t-test to determine whether sig diff between TNBC subpop no 
# and rest:
subpop_t <- t.test(
  split_subpop_no$TNBC, 
  c(split_subpop_no$ER, split_subpop_no$HER2),
  alternative = "greater"
)

# use one-tailed t-test to determine whether sig diff between ER CNV count no 
# and rest:
length_t <- t.test(
  spl_counts$ER, 
  c(spl_counts$HER, spl_counts$TNBC),
  alternative = "less"
)



