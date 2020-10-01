plot_DE <- function(
  seurat_object,
  min_pct,
  logfc_thresh,
  return_thresh, 
  only.pos = FALSE,
  plot_dir,
  filter_sig = TRUE,
  CNA_assoc_only = TRUE,
  table_dir,
  Robject_dir,
  ref_dir,
  file_prefix,
  specific_DE = "none",
  group1 = "none",
  group2 = "none",
  specific_genes = "none",
  return_table = FALSE
) {

  if (specific_DE[1] != "none") {
    # change cluster labels to 'group1' or 'group2':
    idents <- as.character(Idents(seurat_object))
    for (id in specific_DE[[1]]) {
      idents[idents == id] <- group1
    }
    idents[idents != group1] <- group2
    Idents(seurat_object) <- factor(idents)
  }

  if (!file.exists(paste0(table_dir, file_prefix, ".txt"))) {

    # find DE genes using MAST:
    if (specific_DE[1] == "none") {
      DE_res <- FindAllMarkers(
        only.pos = only.pos,
        object = seurat_object,
        min.pct = min_pct, 
        logfc.threshold = logfc_thresh,
        test.use = 'MAST',
        return.thresh = return_thresh
      )
    } else {
      DE_res <- FindMarkers(
        only.pos = only.pos,
        object = seurat_object,
        ident.1 = group1,
        ident.2 = group2,
        min.pct = min_pct, 
        logfc.threshold = logfc_thresh, 
        test.use = 'MAST'
      )
      DE_res$gene <- rownames(DE_res)
    }
  
    # write table:
    if (nrow(DE_res) > 0) {

      # sort with genes with most FC at top:
      if ("cluster" %in% colnames(DE_res)) {
        DE_sorted <- arrange(DE_res, cluster, desc(avg_logFC))
      } else {
        DE_sorted <- arrange(DE_res, desc(avg_logFC))
      }
      
      # filter for significant genes if needed:
      if (filter_sig) {
        DE_sorted <- DE_sorted[DE_sorted$p_val_adj < 0.1,]
      }

      # calculate and add variances of genes for iDEA:
      DE_sorted_counts <- seurat_object@assays$RNA[DE_sorted$gene]
      DE_sorted$variance <- apply(DE_sorted_counts, 1, function(x) (sd(x))^2)
  
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
          paste0(table_dir, file_prefix, ".txt"),
          col.names = TRUE,
          row.names = FALSE,
          quote = FALSE
        )
  
      } else {
  
        write.table(
          DE_sorted, 
          paste0(table_dir, file_prefix, ".txt"),
          col.names = TRUE,
          row.names = FALSE,
          quote = FALSE
        )
  
      }
  
    } else {

      # write dummy file for Snakemake:
      write.table(
        "No DE genes detected", 
        paste0(table_dir, file_prefix, ".txt"),
        col.names = TRUE,
        row.names = FALSE,
        quote = FALSE
      )

    }

  } else {
  
    DE_sorted <- read.table(
      paste0(table_dir, file_prefix, ".txt"),
      header = TRUE
    )

  }

  if (nrow(DE_sorted) > 0) {

    if ("cluster" %in% colnames(DE_sorted)) {
      heatmap_genes <- DE_sorted %>% 
      group_by(cluster) %>% 
      top_n(10, avg_logFC)
    } else {
      heatmap_genes <- DE_sorted
    }
    
    if (specific_genes[1] == "none") {

      if (!file.exists(paste0(
        Robject_dir, file_prefix, "_scaled_gene_seurat_object.Rdata"
      ))) {

        seurat_object <- ScaleData(
          object = seurat_object, 
          features = heatmap_genes$gene,
          vars.to.regress = c("nFeature_RNA", "nCount_RNA")
        )
        saveRDS(seurat_object, 
          paste0(
          Robject_dir, file_prefix, "_scaled_gene_seurat_object.Rdata"
          )
        )

      } else {

        seurat_object <- readRDS(paste0(
          Robject_dir, file_prefix, "_scaled_gene_seurat_object.Rdata"
        ))

      }

      if (length(heatmap_genes$gene) < 50) {
        hmap <- DoHeatmap(
          seurat_object,
          features = as.character(heatmap_genes$gene),
          group.by = "ident"
        ) + theme(
          text = element_text(size = 20)
        )
      } else {
        hmap <- DoHeatmap(
          seurat_object,
          features = as.character(heatmap_genes$gene),
          group.by = "ident"
        ) + theme(
          text = element_text(size = 10)
        )  
      }

      png(
        paste0(plot_dir, file_prefix, ".png"),
        height = 9,
        width = 15,
        res = 300,
        units = "in"
      )
        print(hmap)
      dev.off()

    } else {

      seurat_object <- ScaleData(
        object = seurat_object, 
        features = as.character(specific_genes),
        vars.to.regress = c("nFeature_RNA", "nCount_RNA")
      )

      if (length(specific_genes) < 50) {
        hmap <- DoHeatmap(
          seurat_object,
          features = as.character(specific_genes),
          group.by = "ident"
        ) + theme(
          text = element_text(size = 20)
        )
      } else {
        hmap <- DoHeatmap(
          seurat_object,
          features = as.character(specific_genes),
          group.by = "ident"
        ) + theme(
          text = element_text(size = 10)
        )  
      }

      png(
        paste0(plot_dir, file_prefix, ".png"),
        height = 9,
        width = 15,
        res = 300,
        units = "in"
      )
        print(hmap)
      dev.off()

    } 

    if (return_table) {
      return(DE_sorted)
    }

  } else {
    print("No DE genes identified")

    # create dummy file for Snakemake:
    system(
      paste0(
        "touch ", plot_dir, file_prefix, ".png"
      )
    )
  }

}
