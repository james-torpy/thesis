filter_common_DE_genes <- function(
  gene_data,
  met_genes,
  direction = "up",
  func_dir
) {

  keep_top_common_DE <- dget(paste0(func_dir, "keep_top_common_DE.R"))

  met_df <- gene_data[gene_data$gene %in% met_genes,]
  met_df$gene <- as.character(met_df$gene)
  met_df$sig <- FALSE
  met_df$sig[met_df$p_val_adj < 0.1] <- TRUE

  if (direction == "bi") {

    # split by DE direction:
    split_df <- split(met_df, met_df$avg_logFC > 0)
    names(split_df) <- c("down", "up")
  
    # deal with duplicate genes:
    for (k in 1:length(names(split_df))) {
      # split by gene:
      split_by_gene <- split(split_df[[k]], as.character(split_df[[k]]$gene))
      temp_res_df <- do.call(
        "rbind",
        lapply(split_by_gene, keep_top_common_DE, direction == "bi")
      )
      if (k==1) {
        res_df <- temp_res_df
      } else {
        res_df <- rbind(res_df, temp_res_df)
      }
    }

    res_df$gene[res_df$avg_logFC < 0] <- paste0(
      res_df$gene[res_df$avg_logFC < 0], "_down"
    )
    res_df$gene[res_df$avg_logFC > 0] <- paste0(
      res_df$gene[res_df$avg_logFC > 0], "_up"
    )

  } else {

    # split by gene:
    split_by_gene <- split(met_df, as.character(met_df$gene))
    res_df <- do.call(
      "rbind",
      lapply(split_by_gene, keep_top_common_DE, "up")
    )

  }
  
  return(res_df)

}