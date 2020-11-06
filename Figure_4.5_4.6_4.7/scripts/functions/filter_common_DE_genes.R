filter_common_DE_genes <- function(
  gene_data,
  selected_genes,
  direction,
  func_dir,
  ref_dir,
  Robject_dir
) {

  keep_top_common_DE <- dget(paste0(func_dir, "keep_top_common_DE.R"))

  sample_name <- unique(gene_data$sample_name)
  print(paste0("Fetching common DE genes from ", sample_name))

  if (selected_genes[1] != "none") {
    DE_df <- gene_data[gene_data$gene %in% selected_genes,]
  } else {
    DE_df <- gene_data
  }
  
  DE_df$gene <- as.character(DE_df$gene)
  DE_df$sig <- FALSE
  DE_df$sig[DE_df$p_val_adj < 0.1] <- TRUE

  # split by gene:
  split_by_gene <- split(DE_df, as.character(DE_df$gene))
  res_df <- do.call(
    "rbind",
    lapply(split_by_gene, keep_top_common_DE, direction)
  )
  res_df$symbol <- res_df$gene
  
  return(res_df)

}