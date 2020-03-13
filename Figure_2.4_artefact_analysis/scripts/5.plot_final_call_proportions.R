#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

project_name <- "thesis"
subproject_name <- "Figure_2.6_random_permutation_of_normal"
sample_name <- args[1]
print(paste0("Sample name = ", sample_name))
t_cells_included <- as.logical(args[2])
print(paste0("T-cells included? ", as.character(t_cells_included)))
permutation_proportions <- strsplit(args[3], "_")[[1]]
print(paste0("Permutation proportions ", permutation_proportions))
simulation_numbers <- as.character(
  as.numeric(
    strsplit(args[4], "_")[[1]][1]
  ):as.numeric(
    strsplit(args[4], "_")[[1]][2]
  )
)
print("Simulation numbers: ")
print(simulation_numbers)
analysis_mode <- args[5]
print(paste0("Analysis mode = ", analysis_mode))

#project_name <- "thesis"
#subproject_name <- "Figure_2.6_random_permutation_of_normal"
#sample_name <- "permutated_CID4520N"
#t_cells_included <- TRUE
##permutation_proportions <- strsplit("0.01_0.05_0.1_0.2_0.3_0.4", "_")[[1]]
#permutation_proportions <- strsplit("0.01_0.05", "_")[[1]]
#simulation_numbers <- as.character(
#  as.numeric(
#    strsplit("1_2", "_")[[1]][1]
#  ):as.numeric(
#    strsplit("1_2", "_")[[1]][2]
#  )
#)
#analysis_mode <- "samples"

library(ggplot2)
theme_set(theme_minimal())
library(dplyr)
library(cowplot)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", 
  project_name, "/", subproject_name, "/")
results_dir <- paste0(project_dir, "results/")
func_dir <- paste0(project_dir, "scripts/functions/")
col_dir <- paste0(home_dir, "/R/colour_palettes/")

if (t_cells_included) {
  in_path <- paste0(results_dir, "infercnv/t_cells_included/", 
    sample_name, "/")
} else {
  in_path <- paste0(results_dir, "infercnv/t_cells_excluded/", 
    sample_name, "/")
}

out_dir <- paste0(in_path, "final_plots/")
system(paste0("mkdir -p ", out_dir))

# convert permutation_proportions to list:
permutation_list <- as.list(permutation_proportions)

for (i in 1:length(simulation_numbers)) {

  artefact_list <- lapply(permutation_list, function(x) {

    print(paste0("Fetching counts for simulation ", 
      simulation_numbers[i], " with ", x, " of genes permutated..."))
    
    in_dir <- paste0(in_path, x, "_proportion/", simulation_numbers[i], "/", 
      analysis_mode, "_mode/tables/")
    
    artefact_counts <- read.table(
      paste0(in_dir, "artefact_counts.txt"),
      sep = "\t",
      header = T,
      as.is = T
    )

    artefact_df <- data.frame(
      Proportion = rep(as.numeric(x), 2),
      Type = c("gain", "loss"),
      Count = c(
        artefact_counts["gain",],
        artefact_counts["loss",]
      ),
      Sim_no = rep(simulation_numbers[i], 2)
    )

    return(artefact_df)

  })

  if (i==1) {
    plot_counts <- do.call("rbind", artefact_list)
  } else {
    plot_counts <- rbind(
      plot_counts,
      do.call("rbind", artefact_list)
    )
  }

}

# split df based upon upon permutation proportion:
split_counts1 <- split(plot_counts, plot_counts$Proportion)
split_counts2 <- lapply(split_counts1, function(x) {
  return(split(x, x$Type))
})
split_counts3 <- do.call("c", split_counts2)

final_counts <- lapply(split_counts3, function(x) {
  
  res_df <- data.frame(
    Proportion = x$Proportion[1],
    Type = x$Type[1],
    Count = mean(x$Count),
    SE = sd(x$Count)/sqrt(nrow(x))
  )

})
final_count_df <- do.call("rbind", final_counts)
levels(final_count_df$Type) <- c("gain", "loss")

cols <- c("#BF3667", "#58B9DB")

p <- ggplot(final_count_df, aes(x = Proportion, y = Count)) 
p <- p + geom_line(aes(color = Type))
p <- p + geom_errorbar(aes(ymin=Count-SE, ymax=Count+SE))
p <- p + scale_color_manual(values = cols)
p <- p + scale_x_continuous(
  breaks = c(0, as.numeric(permutation_proportions))
)
p <- p + scale_y_continuous(
  breaks = seq(
    min(plot_counts$Count),
    max(plot_counts$Count)
  )
)
p <- p + theme(
  panel.grid.minor = element_blank(),
  axis.title.x = element_text(size = 9),
  axis.title.y = element_blank(),
  legend.title = element_blank()
)
p <- p + xlab("Number of artefacts")

pdf(
  paste0(out_dir, "artefact_results.pdf"), 
  width = 7, height = 5
)
  print(p)
dev.off()

png(
  paste0(out_dir, "artefact_results.png"), 
  width = 7, height = 5, units = "in", res = 300
)
  print(p)
dev.off()

print(paste0("All output plots in", out_dir))