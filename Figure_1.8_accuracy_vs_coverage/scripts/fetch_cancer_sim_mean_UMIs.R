in_path <- "/share/ScratchGeneral/jamtor/projects/thesis/Figure_2.2_accuracy_vs_coverage/results/infercnv/t_cells_included/CID4520N_cancer_sim/both/"

for ( i in 1:30 ) {
  if (i==1) {
    dfs <- list(read.table(paste0(in_path, i, "/no_downsampling/input_files/input_matrix.txt")))
  } else {
  	dfs[i] <- read.table(paste0(in_path, i, "/no_downsampling/input_files/input_matrix.txt"))
  }
}

mean_umis <- lapply(dfs, function(x) {
  nUMI <- rowSums(x)
  return(mean(nUMI))
})