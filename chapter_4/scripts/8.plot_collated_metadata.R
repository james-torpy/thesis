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
subset_samples <- as.logical(args[5])

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
subset_samples <- FALSE
out_suffix <- "14_tumours"

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
library(dplyr, lib.loc = lib_loc)
library(ggpubr, lib.loc = lib_loc)
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
  remove_artefacts, "/", out_suffix, "/")

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
plot_boxplot <- dget(
  paste0(func_dir, "plot_boxplot.R")
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

    # remove artefacts:
    CNV_data <- lapply(CNV_data, function(x) {
      return(x[grep("artefact", x$call, invert = T),])
    })
  
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

# calculate summary stats:
stats_data <- split(CNV_dist_data, CNV_dist_data$subtype)
stats_data <- append(stats_data, list(CNV_dist_data))
names(stats_data)[length(stats_data)] <- "all"

all_sum_stats <- lapply(stats_data, function(x) {

  sum_stats <- round(
    as.data.frame(
      t(
        apply(
          subset(x, select = -c(sample_name, subtype)),
          2,
          summary
        )
      )
    )
  )

  sum_stats$Mean <- round(
    apply(
      subset(x, select = -c(sample_name, subtype)),
      2,
      mean
    ),
    1
  )
  
  sum_stats$std.dev <- round(
    apply(
      subset(x, select = -c(sample_name, subtype)),
      2,
      sd
    ),
    2
  )

  if (length(unique(x$subtype)) == 1) {
    sum_stats$subtype <- rep(unique(x$subtype), nrow(sum_stats))
  } else {
    sum_stats$subtype <- rep("all", nrow(sum_stats))
  }
  
  return(sum_stats)

})

all_stats_df <- do.call("rbind", all_sum_stats)

write.table(
  all_stats_df,
  paste0(table_dir, "CNA_metadata_summaries.txt"),
  col.names = TRUE,
  row.names = TRUE,
  quote = FALSE
)


################################################################################
### 3. Plot data ###
################################################################################

plot_densities(
  CNV_counts,
  "counts",
  subtype_cols,
  subtype_df,
  plot_dir
)
mean(CNV_counts)

plot_boxplot(
  CNV_counts,
  "counts",
  subtype_cols,
  subtype_df,
  plot_dir,
  lib_loc
)

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
  subtype_cols,
  subtype_df,
  plot_dir
)
mean(length_vec)

plot_boxplot(
  length_vec,
  "lengths",
  subtype_cols,
  subtype_df,
  plot_dir,
  lib_loc
)

plot_violin(
  length_vec,
  "lengths",
  subtype_cols,
  subtype_df,
  plot_dir,
  lib_loc
)

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
  subtype_cols,
  subtype_df,
  plot_dir
)

plot_boxplot(
  genomic_length_vec,
  "genomic_lengths",
  subtype_cols,
  subtype_df,
  plot_dir,
  lib_loc
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

names(CNV_counts) <- gsub(".CNV.*$", "", names(CNV_counts))

names(CNV_lengths) <- paste0(names(CNV_lengths), "_")
CNV_lengths <- unlist(CNV_lengths)
names(CNV_lengths) <- gsub("_.*$", "", names(CNV_lengths))

names(CNV_genomic_lengths) <- paste0(names(CNV_genomic_lengths), "_")
CNV_genomic_lengths <- unlist(CNV_genomic_lengths)
names(CNV_genomic_lengths) <- gsub("_.*$", "", names(CNV_genomic_lengths))

all_data_types <- list(
  nos = subpop_no,
  counts = CNV_counts,
  lengths = CNV_lengths,
  genomic_lengths = CNV_genomic_lengths
)
all_data_nos <- lapply(all_data_types, function(x) round(mean(x), 2))

t_tests <- lapply(all_data_types, function(x) {

  test_df <- data.frame(
    value = x,
    sample = names(x)
  )
  test_df <- merge(test_df, subtype_df, by="sample")

  # split subtype no by subtype
  split_data <- split(test_df, as.character(test_df$subtype))
  
  # use one-tailed t-test to determine whether sig diff between TNBC/ER subpop no 
  # and rest:

  t_res <- list(
    ER_vs_rest = t.test(
      split_data$ER$value, 
      c(split_data$HER2$value, split_data$TNBC$value)
    ),
    HER2_vs_rest = t.test(
      split_data$HER2$value, 
      c(split_data$ER$value, split_data$TNBC$value)
    ),
    TNBC_vs_rest_greater = t.test(
      split_data$TNBC$value, 
      c(split_data$ER$value, split_data$HER2$value),
      alternative = "greater"
    ),
    ER_vs_HER2 = t.test(
      split_data$ER$value, 
      split_data$HER2$value
    ),
    TNBC_vs_ER_greater = t.test(
      split_data$TNBC$value, split_data$ER$value,
      alternative = "greater"
    ),
    TNBC_vs_HER2_greater = t.test(
      split_data$TNBC$value, split_data$HER2$value,
      alternative = "greater"
    )
  )

  t_df <- data.frame(
    comparison = c(
      "ER_vs_rest",
      "HER2_vs_rest",
      "TNBC_vs_rest_greater",
      "ER_vs_HER2",
      "TNBC_vs_ER_greater",
      "TNBC_vs_HER2_greater"
    ),
    p_val = c(
      round(t_res$ER_vs_rest$p.val, 5),
      round(t_res$HER2_vs_rest$p.val, 5),
      round(t_res$TNBC_vs_rest_greater$p.val, 5),
      round(t_res$ER_vs_HER2$p.val, 5),
      round(t_res$TNBC_vs_ER_greater$p.val, 5),
      round(t_res$TNBC_vs_HER2_greater$p.val, 5)
    ),
    mean1 = c(
      round(t_res$ER_vs_rest$estimate[1], 1),
      round(t_res$HER2_vs_rest$estimate[1], 1),
      round(t_res$TNBC_vs_rest_greater$estimate[1], 1),
      round(t_res$ER_vs_HER2$estimate[1], 1),
      round(t_res$TNBC_vs_ER_greater$estimate[1], 1),
      round(t_res$TNBC_vs_HER2_greater$estimate[1], 1)
    ),
    mean2 = c(
      round(t_res$ER_vs_rest$estimate[2], 1),
      round(t_res$HER2_vs_rest$estimate[2], 1),
      round(t_res$TNBC_vs_rest_greater$estimate[2], 1),
      round(t_res$ER_vs_HER2$estimate[2], 1),
      round(t_res$TNBC_vs_ER_greater$estimate[2], 1),
      round(t_res$TNBC_vs_HER2_greater$estimate[2], 1)
    )
  )

  sd_df <- data.frame(
    subtype = c("ER", "HER2", "TNBC"),
    std_dev = c(
      round(sd(split_data$ER$value), 1),
      round(sd(split_data$HER2$value), 1),
      round(sd(split_data$TNBC$value), 1)
    )
  )

  return(
    list(
      t_tests = t_df,
      std_devs = sd_df
    )
  )

})

t_tests$n <- c(
  ER = length(split_data$ER$value),
  HER2 = length(split_data$HER2$value),
  TNBC = length(split_data$TNBC$value)
)

for (i in 1:length(t_tests)) {
  for (j in 1:length(t_tests[[i]])) {
    write.table(
      t_tests[[i]][[j]],
      paste0(table_dir, names(t_tests)[i], "_", names(t_tests[[i]])[j], ".txt"),
      quote = F,
      row.names = F,
      col.names = T
    )
  }
}



















