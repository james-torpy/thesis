#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

RStudio <- FALSE

project_name <- "thesis"
subproject_name <- "Figure_2.3_min_length_and_gap"
print(paste0("Subproject name = ", subproject_name))
args = commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
print(paste0("Sample name = ", sample_name))
t_cells_included <- as.logical(args[2])
print(paste0("T-cells included? ", as.character(t_cells_included)))
simulation_numbers <- args[3]
if ( length(grep("_", simulation_numbers)) > 0 ) {
  simulation_numbers <- as.character(
    as.numeric(
      strsplit(simulation_numbers, "_")[[1]][1]
    ):as.numeric(
      strsplit(simulation_numbers, "_")[[1]][2]
    )
  )
}
print("Simulation numbers are: ")
print(simulation_numbers)
analysis_mode <- args[4]
print(paste0("Analysis mode = ", analysis_mode))

#sample_name <- "CID4520N_cancer_sim"
#t_cells_included <- TRUE
#simulation_numbers <- as.character(
#  as.numeric(
#    strsplit("1_30", "_")[[1]][1]
#  ):as.numeric(
#    strsplit("1_30", "_")[[1]][2]
#  )
#)
#analysis_mode <- "samples"

library(ggplot2)
theme_set(theme_minimal())
library(dplyr)
library(cowplot)

if (RStudio) {
  home_dir <- "/Users/jamestorpy/clusterHome/"  
} else { 
  lib_loc <- "/share/ScratchGeneral/jamtor/R/3.5dev/" 
  home_dir <- "/share/ScratchGeneral/jamtor/" 
}
library(naturalsort, lib.loc = lib_loc)

project_dir <- paste0(home_dir, "projects/", 
  project_name, "/", subproject_name, "/")
results_dir <- seurat_path <- paste0(project_dir, "results/")
func_dir <- paste0(project_dir, "scripts/functions/")
col_dir <- paste0(home_dir, "/R/colour_palettes/")

if (t_cells_included) {
  in_path <- paste0(results_dir, "infercnv/t_cells_included/", 
    sample_name, "/")
} else {
  in_path <- paste0(results_dir, "infercnv/t_cells_excluded/", 
    sample_name, "/")
}

plot_dir <- paste0(in_path, "final_plots/")
system(paste0("mkdir -p ", plot_dir))

# create list to process CNV, gain_gap and loss_gap data in parallel :
type_list <- list("CNV", "gain_gap", "loss_gap")

for (i in 1:length(simulation_numbers)) {

  calls_list <- lapply(type_list, function(x) {

    print(paste0("Fetching calls for ", x, " simulation ", 
      simulation_numbers[i], "..."))

    if (x == "CNV") {

      in_dir <- paste0(in_path, "/", x, "/", simulation_numbers[i],
        "/", analysis_mode, "_mode/tables/")
      calls <- read.table(
        paste0(in_dir, "CNV_calls.txt"),
        sep = "\t",
        header = T,
        as.is = T
      )

      ######
      # replace call entries with underscored 'feature' calls if needed:
      calls$call <- gsub("CNV|gap", "feature", calls$call)
      calls$call <- gsub(" ", "_", calls$call)
      ######

    } else {

      CNV_type <- strsplit(x, "_")[[1]][1]
      in_dir <- paste0(in_path, "/gap/", simulation_numbers[i],
        "/", CNV_type, "/", analysis_mode, "_mode/tables/")
      calls <- read.table(
        paste0(in_dir, "gap_calls.txt"),
        sep = "\t",
        header = T,
        as.is = T
      )
      ######
      # reorder columns if needed:
      calls <- data.frame(
        length = calls$length,
        call = calls$call
      )
      ######
      # replace call entries with underscored 'feature' calls if needed:
      calls$call <- gsub("CNV|gap", "feature", calls$call)
      calls$call <- gsub(" ", "_", calls$call)
      ######

    }
  
    calls$simulation_no <- simulation_numbers[i]

    return(calls)

  })

  # split CNV calls into gains and losses and add as separate elements:
  calls_list <- list(
    gain_CNV = split(calls_list[[1]], calls_list[[1]]$type)[[1]],
    loss_CNV = split(calls_list[[1]], calls_list[[1]]$type)[[2]],
    gain_gap = calls_list[[2]],
    loss_gap = calls_list[[3]] 
  )
  calls_list$gain_CNV <- subset(calls_list$gain_CNV, select = -type)
  calls_list$loss_CNV <- subset(calls_list$loss_CNV, select = -type)

  # add each element of calls_list to each element of plot_calls:
  writeLines("\n")
  print(paste0("Adding calls for simulation ", 
      simulation_numbers[i], " to all simulation calls list..."))
  writeLines("\n")

  if (i==1) {
    plot_calls <- calls_list
  } else {
    for (j in 1:length(plot_calls)) {
      plot_calls[[j]] <- rbind(plot_calls[[j]], calls_list[[j]])
    }
  }

}

# split df based upon feature length and count each call:
plot_call_count <- lapply(plot_calls, function(x) {

  split_calls <- split(x, x$length)
  call_counts <- lapply(split_calls, function(y) {
    counts_df <- data.frame(
      number = c(
        called = length(which(y$call == "feature_called")),
        not_called = length(which(y$call == "feature_not_called"))
      )
    )
    counts_df["total", ] <- counts_df["called",] + 
      counts_df["not_called",]
    counts_df[
      "percentage_correct_calls", 
    ] <- round(counts_df["called", ]/counts_df["total", ]*100, 1)
    return(counts_df)
  })
  names(call_counts) <- paste0("length_", names(call_counts))

  return(call_counts)

})

# write tables:
for (l in 1:length(plot_call_count)) {
  for (m in 1:length(plot_call_count[[l]])) {
    table_dir <- paste0(in_path, "final_tables/")
    system(paste0("mkdir -p ", table_dir))
    write.table(
      plot_call_count[[l]][[m]],
      paste0(table_dir, 
        names(plot_call_count)[l], "_", 
        names(plot_call_count[[l]])[m], "_calls"
      ),
      sep = "\t",
      quote = F,
      row.names = T,
      col.names = F
    )
  }
}

# keep only proportion of correct calls:
correct_call_proportions <- lapply(plot_call_count, function(x) {
  correct_proportions <- lapply(x, function(y) y["percentage_correct_calls",])
  return(
    data.frame(
      percent_correct_calls = do.call("c", correct_proportions),
      length = factor(
      	gsub("length_", "", names(correct_proportions)),
      	levels = naturalsort(gsub("length_", "", names(correct_proportions)))
      )
    )
  )
})

cols <- rep(
  read.table(
    paste0(col_dir, "colour_palette_1.txt"),
    sep = "\t",
    header = F,
    as.is = T,
    comment.char = ""
  )[,1][1:2],
  2
)

for (n in 1:length(correct_call_proportions)) {
  col <- cols[n]
  p <- ggplot(
    correct_call_proportions[[n]], 
    aes(x = length, y = percent_correct_calls)
  ) 
  p <- p + geom_bar(stat = "identity", width = 0.7, fill = col)
  p <- p + scale_color_manual(values = cols)
  if (length(grep("CNV", names(correct_call_proportions)[n])) > 0) {
    p <- p + xlab("CNV length")
    p <- p + ylab(
      paste0("Correct detections (% of total CNV number)")
    )
  } else {
    p <- p + xlab("Gap length")
    p <- p + ylab(
      paste0("Correct detections (% of total gap number)")
    )
  }
    
  pdf(
    paste0(
      plot_dir, "proportion_correct_calls_for_", 
      names(correct_call_proportions)[n], "s.pdf"), 
    width = 7, height = 5
  )
    print(p)
  dev.off()
  
  png(
    paste0(
      plot_dir, "proportion_correct_calls_for_", 
      names(correct_call_proportions)[n], "s.png"), 
    width = 7, height = 5, units = "in", res = 300
  )
    print(p)
  dev.off()
  
}

