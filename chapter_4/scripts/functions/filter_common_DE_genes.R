filter_common_DE_genes <- function(
  gene_data,
  selected_genes,
  direction = "bi",
  CNA_assoc_only = TRUE,
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

  if (direction == "bi") {

    # split by DE direction:
    split_df <- split(DE_df, DE_df$avg_logFC > 0)
    names(split_df) <- c("down", "up")
  
    # deal with duplicate genes:
    for (k in 1:length(names(split_df))) {
      # split by gene:
      split_by_gene <- split(split_df[[k]], as.character(split_df[[k]]$gene))
      temp_res_df <- do.call(
        "rbind",
        lapply(
          split_by_gene, 
          keep_top_common_DE, 
          direction = names(split_df)[k]
        )
      )
      if (k==1) {
        res_df <- temp_res_df
      } else {
        res_df <- rbind(res_df, temp_res_df)
      }
    }
    res_df$symbol <- res_df$gene

    res_df$gene[res_df$avg_logFC < 0] <- paste0(
      res_df$gene[res_df$avg_logFC < 0], "_down"
    )
    res_df$gene[res_df$avg_logFC > 0] <- paste0(
      res_df$gene[res_df$avg_logFC > 0], "_up"
    )

  } else {

    # split by gene:
    split_by_gene <- split(DE_df, as.character(DE_df$gene))
    res_df <- do.call(
      "rbind",
      lapply(split_by_gene, keep_top_common_DE, "up")
    )
    res_df$symbol <- res_df$gene

  }

  if (CNA_assoc_only) {
  
    # define sample Robject dir:
    sample_Robject_dir <- gsub("common_gene_expression", sample_name, Robject_dir)

    # load CNA and gene coordinate data:
    gene_coords <- read.table(paste0(ref_dir, "infercnv_gene_order.txt"))
    colnames(gene_coords) <- c("gene", "chr", "start", "end")
    CNA_data <- readRDS(paste0(sample_Robject_dir, "CNV_indices_and_lengths.Rdata"))
    
    # identify genes within CNA areas in at least one subpopulation:
    CNA_genes <- lapply(CNA_data, function(x) {
      for (i in 1:nrow(x)) {
        gene_seg <- gene_coords$gene[
          (which(gene_coords$gene == x$start_gene[i])):(which(gene_coords$gene == x$end_gene[i]))
        ]
        if (i==1) {
          subset_CNA_genes <- as.character(gene_seg)
        } else {
          subset_CNA_genes <- c(subset_CNA_genes, as.character(gene_seg))
        }
      }
      return(subset_CNA_genes)      
    })
    CNA_genes <- unique(unlist(CNA_genes))

    # filter out genes not in CNAs:
    res_df <- res_df[res_df$symbol %in% CNA_genes,]
    
  }
  
  return(res_df)

}