DE_heatmap <- function(
  DE,
  seurat_object,
  file_prefix,
  Robject_dir,
  no_genes_returned,
  data_type = "non_CNA",
  specific_genes = "none"
) {

  if (data_type == "CNA_assoc") {
    file_prefix <- paste0(file_prefix, "_CNA_assoc")
  }
  if ("cluster" %in% colnames(DE)) {
    heatmap_genes <- DE %>% 
    group_by(cluster) %>% 
    top_n(10, avg_logFC)
  } else {
    up_genes <- DE %>% 
      top_n(no_genes_returned/2, avg_logFC)
    down_genes <- DE %>% 
      top_n(-(no_genes_returned/2), avg_logFC)
    heatmap_genes <- rbind(up_genes, down_genes)
  }
  
  if (specific_genes[1] == "none") {

    # scale DE genes:
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

    # adjust group labels:
    levels(Idents(seurat_object)) <- gsub(
      "_", " ", levels(Idents(seurat_object))
    )

    # generate heatmap:
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

    pdf(
      paste0(plot_dir, file_prefix, ".pdf"),
      height = 9,
      width = 15
    )
      print(hmap)
    dev.off()

    png(
      paste0(plot_dir, file_prefix, ".png"),
      height = 9,
      width = 15,
      res = 300,
      units = "in"
    )
      print(hmap)
    dev.off()

    return(heatmap_genes)

  } else {

    DE <- DE[order(DE$avg_logFC, decreasing = TRUE),]
    # filter specific genes for only DE genes:
    specific_DE_genes <- DE$gene[DE$gene %in% specific_genes]
    # reduce to required no of genes:
    if (length(specific_DE_genes) > no_genes_returned) {
      specific_DE_genes <- specific_DE_genes[1:no_genes_returned]
    }

    # adjust group labels:
    levels(Idents(seurat_object)) <- gsub(
      "_", " ", levels(Idents(seurat_object))
    )

    # generate heatmap:
    if (length(specific_DE_genes) > 0) {
      seurat_object <- ScaleData(
        object = seurat_object, 
        features = as.character(specific_DE_genes),
        vars.to.regress = c("nFeature_RNA", "nCount_RNA")
      )

      if (length(specific_DE_genes) < 50) {
        hmap <- DoHeatmap(
          seurat_object,
          features = as.character(specific_DE_genes),
          group.by = "ident"
        ) + theme(
          text = element_text(size = 20)
        )
      } else {
        hmap <- DoHeatmap(
          seurat_object,
          features = as.character(specific_DE_genes),
          group.by = "ident"
        ) + theme(
          text = element_text(size = 10)
        )  
      }

      pdf(
        paste0(plot_dir, file_prefix, ".pdf"),
        height = 9,
        width = 15
      )
        print(hmap)
      dev.off()

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
      print("No CNA-associated genes are in the top DE genes")

      # create dummy file for snakemake:
      system(paste0("touch ", plot_dir, file_prefix, ".png"))
    }

    return(specific_DE_genes)

  }

}
