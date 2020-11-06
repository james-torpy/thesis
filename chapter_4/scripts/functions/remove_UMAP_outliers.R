remove_UMAP_outliers <- function(seurat_object) {

  embeddings <- eval(parse(
    text = paste0("seurat_object@reductions$UMAP", epi_PC, "@cell.embeddings")
  ))
  mean_embeddings <- apply(embeddings, 2, mean)
  std_dev_embeddings <- apply(embeddings, 2,sd)
  
  x_outliers <- rownames(embeddings)[
    embeddings[,1] > (mean_embeddings[1] + 
      (outlier_sd_multiplier*std_dev_embeddings[1])) | 
      embeddings[,1] < (mean_embeddings[1] - 
      (outlier_sd_multiplier*std_dev_embeddings[1]))
  ]
  
  no_outlier <- rownames(embeddings)[
    !(rownames(embeddings) %in% x_outliers)
  ]
  
  y_outliers <- rownames(embeddings)[
    embeddings[,2] > (mean_embeddings[2] + (3*std_dev_embeddings[2])) | 
    embeddings[,2] < (mean_embeddings[2] - (3*std_dev_embeddings[2]))
  ]
  no_outlier <- no_outlier[
    !(no_outlier %in% y_outliers)
  ]
  
  if (!exists("no_outlier")) {
    no_outlier <- rownames(embeddings)
  }
  
  # remove epithelial clusters with <5 cells:
  split_idents <- split(
    Idents(seurat_object),
    Idents(seurat_object)
  )
  small_clusters <- levels(Idents(seurat_object))[
    unlist(
      lapply(split_idents, function(x) length(x) < 5)
    )
  ]
  no_outlier <- no_outlier[
    no_outlier %in% names(Idents(seurat_object))[
      !(Idents(seurat_object) %in% small_clusters)
    ]
  ]
  
  # only include cells in epithelial_metadata:
  include_cells <- rownames(seurat_object@meta.data)
  include_cells <- include_cells[
    include_cells %in% epithelial_metadata$cell_ids
  ]
  
  no_outlier <- no_outlier[
    no_outlier %in% epithelial_metadata$cell_ids
  ]
  
  # subset seurat object to remove outliers, rename clusters and save:
  return(
    list(
      subset(
        seurat_object,
        cells = include_cells
      ),
      subset(
        seurat_object,
        cells = no_outlier
      ),
      include_cells
    )
  )
}