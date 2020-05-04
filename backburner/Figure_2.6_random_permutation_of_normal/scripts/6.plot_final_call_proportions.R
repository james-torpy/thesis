#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

project_name <- "thesis"
subproject_name <- "Figure_2.6_random_permutation_of_normal"
sample_name <- args[1]
print(paste0("Sample name = ", sample_name))
permutation_proportions <- strsplit(args[2], "_")[[1]]
print(paste0("Permutation proportions ", permutation_proportions))
simulation_numbers <- as.character(
  as.numeric(
    strsplit(args[3], "_")[[1]][1]
  ):as.numeric(
    strsplit(args[3], "_")[[1]][2]
  )
)
print("Simulation numbers: ")
print(simulation_numbers)
analysis_mode <- args[4]
print(paste0("Analysis mode = ", analysis_mode))
min_artefact_proportion <- args[5]
print(
  paste0(
    "Min proportion of cells for artefact detection = ", 
    min_artefact_proportion
  )
)
min_artefact_length <- args[6]
print(
  paste0(
    "Min length for artefact detection = ", min_artefact_length
  )
)

#project_name <- "thesis"
#subproject_name <- "Figure_2.6_random_permutation_of_normal"
#sample_name <- "permutated_CID4520N"
##permutation_proportions <- strsplit("0.01_0.05_0.1_0.2_0.3_0.4", "_")[[1]]
#permutation_proportions <- strsplit("0.3_0.4", "_")[[1]]
#simulation_numbers <- as.character(
#  as.numeric(
#    strsplit("1_30", "_")[[1]][1]
#  ):as.numeric(
#    strsplit("1_30", "_")[[1]][2]
#  )
#)
#analysis_mode <- "samples"
#min_artefact_proportion <- 0.4
#min_artefact_length <- 20

library(ggplot2)
theme_set(theme_minimal())
library(dplyr)
library(cowplot)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", 
  project_name, "/", subproject_name, "/")
results_dir <- paste0(project_dir, "results/")
func_dir <- paste0(project_dir, "scripts/functions/")
col_dir <- paste0(home_dir, "R/colour_palettes/")
perm_path <- paste0(results_dir, "permutated_CID4520N/")

in_path <- paste0(results_dir, "infercnv/", sample_name, "/") 

plot_dir <- paste0(
  in_path, 
  "final_plots/",
  "/", min_artefact_proportion, "_min_artefact_proportion", 
  "/", min_artefact_length, "_min_artefact_length/"
)
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(
  in_path, 
  "final_tables/",
  "/", min_artefact_proportion, "_min_artefact_proportion", 
  "/", min_artefact_length, "_min_artefact_length/"
)
system(paste0("mkdir -p ", table_dir))


####################################################################################
### 1. Load data ###
####################################################################################

# convert permutation_proportions to list:
permutation_list <- as.list(permutation_proportions)

for (i in 1:length(simulation_numbers)) {

  artefact_list <- lapply(permutation_list, function(x) {

    print(paste0("Fetching counts for simulation ", 
      simulation_numbers[i], " with ", x, " of genes permutated..."))
    
    in_dir <- paste0(in_path, x, "_proportion/", simulation_numbers[i], "/", 
      analysis_mode, "_mode/", min_artefact_proportion, "_min_artefact_proportion/",
      min_artefact_length, "_min_artefact_length/tables/")
    
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
    Count = round(mean(x$Count), 0),
    SE = sd(x$Count)/sqrt(nrow(x))
  )

})
final_count_df <- do.call("rbind", final_counts)
levels(final_count_df$Type) <- c("gain", "loss")
final_count_df$Proportion <- as.character(final_count_df$Proportion)

if (!file.exists(paste0(table_dir, "final_counts.txt"))) {
  write.table(
    final_count_df, 
    paste0(table_dir, "final_counts.txt"),
    sep = "\t",
    row.names = F,
    col.names = T,
    quote = F
  )
}


####################################################################################
### 2. Record metadata ###
####################################################################################

# fetch artefact lengths and permutation numbers:
for (i in 1:length(simulation_numbers)) {

  permutation_metadata <- lapply(permutation_list, function(x) {
    
    # load artefact lengths:
    print(paste0("Fetching artefact lengths for simulation ", 
      simulation_numbers[i], " with ", x, " of genes permutated..."))

    in_dir <- paste0(in_path, x, "_proportion/", simulation_numbers[i], "/", 
      analysis_mode, "_mode/", min_artefact_proportion, "_min_artefact_proportion/",
      min_artefact_length, "_min_artefact_length/tables/")

    artefact_lengths <- read.table(
      paste0(in_dir, "artefact_lengths.txt"),
      sep = "\t",
      header = T,
      as.is = T
    )

    # load permuatation records:
    print(paste0("Fetching number genes permutated for simulation ", 
      simulation_numbers[i], " with ", x, " of genes permutated..."))
    perm_dir <- paste0(perm_path, x, "_proportion/", simulation_numbers[i], 
      "/Rdata/")
    permutation_record <- readRDS(paste0(perm_dir, 
    	"3.final_permutation_data.Rdata"))

    return(
      data.frame(
        proportion = x,
      	mean_artefact_length = mean(artefact_lengths$length),
      	permutation_number = nrow(permutation_record$permutation_record)
      )
    )

  })
  names(permutation_metadata) <- unlist(permutation_list)

  if (i==1) {
  	metadata_df <- do.call("rbind", permutation_metadata)
  } else {
  	metadata_df <- rbind(
  	  metadata_df,
  	  do.call("rbind", permutation_metadata)
  	)
  }
  
}

split_metadata <- split(metadata_df, metadata_df$proportion)
mean_metadata <- lapply(split_metadata, function(x) {
  return(
    data.frame(
      proportion = x$proportion[1],
      mean_artefact_length = mean(x$mean_artefact_length),
      permutation_number = x$permutation_number[1]
    )
  )
})
final_metadata <- do.call("rbind", mean_metadata)

write.table(
  final_metadata,
  paste0(table_dir, "mean_artefect_length_and_no_genes_permuatated.txt"),
  sep = "\t",
  quote = F,
  row.names = T,
  col.names = T
)


####################################################################################
### 3. Plot artefact numbers ###
####################################################################################

cols <- c("#BF3667", "#58B9DB")

p <- ggplot(final_count_df, aes(x = Proportion, y = Count, group = Type)) 
p <- p + geom_line(aes(color = Type))
p <- p + geom_errorbar(aes(ymin = Count-SE, ymax = Count+SE, color = Type))
p <- p + scale_color_manual(values = cols)
p <- p + scale_x_discrete(
  labels = paste0(
    final_count_df$Proportion[c(TRUE, FALSE)], 
    "\n(", final_metadata$permutation_number, ")")
)
p <- p + scale_y_continuous(
  breaks = seq(
    min(plot_counts$Count),
    max(plot_counts$Count),
    5
  )
)
p <- p + theme(
  panel.grid.minor = element_blank(),
  axis.title.x = element_text(size = 11),
  axis.title.y = element_text(size = 11),
  axis.text.x = element_text(size = 10),
  axis.text.y = element_text(size = 10),
  legend.title = element_blank()
)
p <- p + ylab("Number of artefacts detected")
p <- p + xlab("Proportion of genes permutated")

pdf(
  paste0(plot_dir, "artefact_results.pdf"), 
  width = 7, height = 5
)
  print(p)
dev.off()

png(
  paste0(plot_dir, "artefact_results.png"), 
  width = 7, height = 5, units = "in", res = 300
)
  print(p)
dev.off()

print(paste0("All output plots in", plot_dir))