plot_DE <- function(
  seurat_object,
  min_pct,
  logfc_thresh,
  return_thresh,
  plot_dir,
  filter_sig = TRUE,
  CNA_assoc_only = TRUE,
  table_dir,
  Robject_dir,
  ref_dir
) {
  
  if (CNA_assoc_only) {
    filename <- paste0(
        table_dir, 
        "subpop_DE_min_pct_", min_pct, 
        "_logfc_", logfc_thresh, 
        "_p_val_", return_thresh,
        "_CNA_assoc_only.txt"
      )
  } else {
    filename <- paste0(
      table_dir, 
      "subpop_DE_min_pct_", min_pct, 
      "_logfc_", logfc_thresh, 
      "_p_val_", return_thresh,
      ".txt"
    )
  }


  # fetch scaled data slot of Seurat:
  scaled_genes <- GetAssayData(seurat_object, slot = "scale.data", assay = "RNA")

  if (!file.exists(filename)) {

    # find DE genes using MAST:
    DE_res <- FindAllMarkers(
      only.pos = F,
      object = seurat_object,
      min.pct = min_pct, 
      #min.pct = 0.8,
      logfc.threshold = logfc_thresh, 
      #logfc.threshold = 0.5,
      test.use = 'MAST',
      return.thresh = return_thresh
      #return.thresh = 0.01
    )
  
    # write table:
    if (nrow(DE_res) > 0) {

      # sort with genes with most FC at top:
      DE_sorted <- arrange(DE_res, cluster, desc(avg_logFC))
  
      # filter for significant genes if needed:
      if (filter_sig) {
        DE_sorted <- DE_sorted[DE_sorted$p_val_adj < 0.1,]
        # calculate and add variances of genes for iDEA:
        DE_sorted_counts <- seurat_sub@assays$RNA[DE_sorted$gene]
        DE_sorted$variance <- apply(DE_sorted_counts, 1, function(x) (sd(x))^2)
      }
  
      if (CNA_assoc_only) {
  
        # load CNA and gene coordinate data:
        gene_coords <- read.table(paste0(ref_dir, "infercnv_gene_order.txt"))
        colnames(gene_coords) <- c("gene", "chr", "start", "end")
        CNA_data <- readRDS(paste0(Robject_dir, "CNV_indices_and_lengths.Rdata"))
        
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
        DE_sorted <- DE_sorted[DE_sorted$gene %in% CNA_genes,]
  
        write.table(
          DE_sorted, 
          filename,
          col.names = TRUE,
          row.names = FALSE,
          quote = FALSE
        )
  
      } else {
  
        write.table(
          DE_sorted, 
          filename,
          col.names = TRUE,
          row.names = FALSE,
          quote = FALSE
        )
  
      }
  
    } else {
  
      if (CNA_assoc_only) {
  
        DE_sorted <- read.table(
          paste0(
            table_dir, 
            "subpop_DE_min_pct_", min_pct, 
            "_logfc_", logfc_thresh, 
            "_p_val_", return_thresh,
            "_CNA_assoc_only.txt"
          ),
          header = TRUE
        )
  
      } else {
  
        DE_sorted <- read.table(
          paste0(
            table_dir, 
            "subpop_DE_min_pct_", min_pct, 
            "_logfc_", logfc_thresh, 
            "_p_val_", return_thresh,
            ".txt"
          ),
          header = TRUE
        )
  
      }

    }

  } else {
  
    if (CNA_assoc_only) {

      DE_sorted <- read.table(
        paste0(
          table_dir, 
          "subpop_DE_min_pct_", min_pct, 
          "_logfc_", logfc_thresh, 
          "_p_val_", return_thresh,
          "_CNA_assoc_only.txt"
        ),
        header = TRUE
      )

    } else {

      DE_sorted <- read.table(
        paste0(
          table_dir, 
          "subpop_DE_min_pct_", min_pct, 
          "_logfc_", logfc_thresh, 
          "_p_val_", return_thresh,
          ".txt"
        ),
        header = TRUE
      )

    }

  }

  if (nrow(DE_sorted) > 0) {

    heatmap_genes <- DE_sorted %>% 
    group_by(cluster) %>% 
    top_n(10, avg_logFC)
  
    if (any(heatmap_genes$gene %in% rownames(scaled_genes))) {

      hmap <- DoHeatmap(
        seurat_object,
        features = as.character(heatmap_genes$gene),
        group.by = "ident"
      ) + theme(
        text = element_text(size = 20)
      )

      if (CNA_assoc_only) {

        png(
          paste0(
            plot_dir, 
            "subpop_DE_min_pct_", min_pct, 
            "_logfc_", logfc_thresh, 
            "_p_val_", return_thresh,
            "_DE_heatmap_CNA_assoc_only.png"
          ),
          height = 9,
          width = 15,
          res = 300,
          units = "in"
        )
          print(hmap)
        dev.off()

      } else {

        png(
          paste0(
            plot_dir, 
            "subpop_DE_min_pct_", min_pct, 
            "_logfc_", logfc_thresh, 
            "_p_val_", return_thresh,
            "_DE_heatmap.png"
          ),
          height = 9,
          width = 15,
          res = 300,
          units = "in"
        )
          print(hmap)
        dev.off()

      }

    } else {
      print("No DE genes in scaled gene list")
    }

  } else {
    print("No DE genes identified")
  }
 
}
