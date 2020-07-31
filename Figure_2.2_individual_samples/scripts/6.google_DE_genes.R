
project_name <- "thesis"
subproject_name <- "Figure_2.2_individual_samples"
sample_name <- "CID4463"
subcluster_method <- "random_trees"
subcluster_p <- "0.05"
if (subcluster_p != "none") {
  subcluster_p <- as.numeric(subcluster_p)
}
coverage_filter <- "filtered"
remove_artefacts <- "artefacts_not_removed"
remove_normals <- "normals_removed"

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))

library(searcher)
library(dplyr)

home_dir <- "/Users/jamestorpy/clusterHome/"
project_dir <- paste0(home_dir, "projects/", project_name, "/",
                      subproject_name, "/")
results_dir <- paste0(project_dir, "results/")
out_path <- paste0(
  results_dir, "infercnv/", sample_name, "/",
  coverage_filter, "/", subcluster_method, "/p_", subcluster_p, "/", 
  remove_artefacts, "/"
)
table_dir <- paste0(out_path, "tables/")

# load DE genes:
DE <- read.table(
  paste0(table_dir, "subpop_DE.txt"),
  header=T
)
top_DE <- DE %>% 
  group_by(cluster) %>% 
  top_n(10, avg_logFC)
gene_list <- top_DE$gene

# google DE genes:
for (g in gene_list) {
  
  print(paste0("Googling ", g))
  
  search_site(
    paste0(g, " + genecards + neuro"),
    site = "google",
    rlang = TRUE
  )
  
}

top_subset <- top_DE[!(top_DE$cluster %in% c("CNV_1", "CNV_2")),]
gene_subset <- top_subset$gene

for (g in gene_subset) {
  
  print(paste0("Googling ", g))
  
  search_site(
    paste0(g, " + neuroendocrine"),
    site = "google",
    rlang = TRUE
  )
  
}
