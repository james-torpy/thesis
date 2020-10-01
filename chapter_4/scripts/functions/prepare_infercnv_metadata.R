prepare_infercnv_metadata <- function(
  seurat_object, 
  count_df, 
  for_infercnv, 
  garnett, 
  manual_epithelial,
  exclude_clusters) {

  numbers_only <- function(x) !grepl("\\D", x)

  # annotate cell identities using garnett_call_ext_major metadata column:
  if (
    numbers_only(Idents(seurat_object)[1])
  ) {

    annotated_idents <- as.character(Idents(seurat_10X))

    if (exclude_clusters[1] != "none") {
      annotated_idents[annotated_idents %in% exclude_clusters] <- paste0(
        "Exclude_", annotated_idents[annotated_idents %in% exclude_clusters]
      )  
    }

    if (manual_epithelial[1] != "none") {
    
      
      ident_ids <- names(Idents(seurat_10X))
      annotated_idents[annotated_idents %in% manual_epithelial] <- paste0(
        "Epithelial_", annotated_idents[annotated_idents %in% manual_epithelial]
      )
      annotated_idents[grep("Epithelial", annotated_idents, invert=T)] <- paste0(
        "Stromal_", annotated_idents[grep("Epithelial", annotated_idents, invert=T)]
      )
    
    } else {
    
      annotated_idents <- gsub(
        " ", 
        "_", 
        paste0(
          eval(parse(text=paste0("seurat_object@meta.data$", garnett))), " ", 
          Idents(seurat_object)
        )
      )
    
    }
    
  
    # remove sample id if present:
    annotated_idents <- gsub("CID.*_", "", annotated_idents)
    Idents(seurat_object) <- factor(annotated_idents)
  
    # create infercnv metadata:
    temp_metadata <- data.frame(cell_ids = names(Idents(seurat_object)), 
      cell_type = Idents(seurat_object), stringsAsFactors = F)
    temp_metadata <- temp_metadata[order(temp_metadata$cell_type),]
  
    # throw unassigned cells and CAFs (as they may be malignant):
    temp_metadata <- temp_metadata[grep("[u,U]nknown|[u,U]nassigned|CAF", 
      temp_metadata$cell_type, invert=T),]
  
    # only include cells present in count_df:
    print(paste0("No cells in metadata df before filtering for those in count df = ", 
      nrow(temp_metadata)))
    temp_metadata <- temp_metadata[temp_metadata$cell_ids %in% colnames(count_df),]
    print(paste0("No cells in metadata df after filtering for those in count df = ", 
      nrow(temp_metadata)))

    if (for_infercnv) {
      # label cells in clusters of < 2 cells as 'outliers' as these will break InferCNV:
      temp_metadata$cell_type <- as.character(temp_metadata$cell_type)
      cluster_list <- split(temp_metadata, temp_metadata$cell_type)
      metadata_outliers_labelled <- do.call(
        "rbind", lapply(cluster_list, function(x) {
          if (nrow(x) < 2) {
            x$cell_type <- gsub("_[0-9].*$", "_outlier", x$cell_type)
          }
          return(x)
        })
      )
  
      # remove outlier cell types with <2 cells:
      cluster_list2 <- split(metadata_outliers_labelled, 
        metadata_outliers_labelled$cell_type)
      i=1
      metadata_final <- do.call(
        "rbind", lapply(cluster_list2, function(x) {
          if (nrow(x) < 2) {
            print(paste0(cluster_list[[i]]$cell_type[1], 
              " removed as contained <2 cells"))
            i <<- i+1
            return(NULL)
          } else {
            i <<- i+1
            return(x)
          }
        })
      )
    } else {
      metadata_final <- temp_metadata
    }
    rownames(metadata_final) <- metadata_final$cell_ids 
    
    # record number per cell type:
    number_per_cell_type <- as.data.frame(table(metadata_final$cell_type))
  
    # create and label results list:
    result_list <- list(metadata_final, number_per_cell_type, seurat_object)
    names(result_list) <- c("metadata", "number_per_group", "seurat")
  
    return(result_list)
  } else {
  # create infercnv metadata:
  temp_metadata <- data.frame(cell_ids = names(Idents(seurat_object)), 
    cell_type = Idents(seurat_object), stringsAsFactors = F)
  temp_metadata <- temp_metadata[order(temp_metadata$cell_type),]

  # throw unassigned cells and CAFs (as they may be malignant):
  temp_metadata <- temp_metadata[grep("[u,U]nknown|[u,U]nassigned|CAF", 
    temp_metadata$cell_type, invert=T),]

  # only include cells present in count_df:
  print(paste0("No cells in metadata df before filtering for those in count df = ", 
    nrow(temp_metadata)))
  temp_metadata <- temp_metadata[temp_metadata$cell_ids %in% colnames(count_df),]
  print(paste0("No cells in metadata df after filtering for those in count df = ", 
    nrow(temp_metadata)))
  if (for_infercnv) {
    # label cells in clusters of < 2 cells as 'outliers' as these will break InferCNV:
    temp_metadata$cell_type <- as.character(temp_metadata$cell_type)
    cluster_list <- split(temp_metadata, temp_metadata$cell_type)
    metadata_outliers_labelled <- do.call(
      "rbind", lapply(cluster_list, function(x) {
        if (nrow(x) < 2) {
          x$cell_type <- gsub("_[0-9].*$", "_outlier", x$cell_type)
        }
        return(x)
      })
    )

    # remove outlier cell types with <2 cells:
    cluster_list2 <- split(metadata_outliers_labelled, 
      metadata_outliers_labelled$cell_type)
    i=1
    metadata_final <- do.call(
      "rbind", lapply(cluster_list2, function(x) {
        if (nrow(x) < 2) {
          print(paste0(cluster_list[[i]]$cell_type[1], 
            " removed as contained <2 cells"))
          i <<- i+1
          return(NULL)
        } else {
          i <<- i+1
          return(x)
        }
      })
    )
  } else {
    metadata_final <- temp_metadata
  }
  rownames(metadata_final) <- metadata_final$cell_ids 
  
  # record number per cell type:
  number_per_cell_type <- as.data.frame(table(metadata_final$cell_type))

  # create and label results list:
  result_list <- list(metadata_final, number_per_cell_type, seurat_object)
  names(result_list) <- c("metadata", "number_per_group", "seurat")

  return(result_list)
  }
}