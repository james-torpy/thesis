create_matched_heatmap_annotations <- function(
  combined_data,
  type_cols,
  subcluster_cols,
  func_dir,
  ref_dir
) {

  # load function to annotate chromosome boundaries:
  fetch_chromosome_boundaries <- dget(paste0(func_dir, 
  "fetch_chromosome_boundaries.R"))

  # order metadata by sample:
  combined_data$metadata <- combined_data$metadata[
    order(combined_data$metadata$sample),
  ]
  rownames(combined_data$metadata) <- combined_data$metadata$cell_ids
  
  # order heatmap the same way as metadata:
  combined_data$heatmap <- combined_data$heatmap[rownames(combined_data$metadata),]
  print(paste0("Are heatmap and metadata rownames identical? ",
    identical(rownames(combined_data$heatmap), rownames(combined_data$metadata))))
  
  # create type annotation:
  type_annotation_df <- subset(combined_data$metadata, select = type)
  type_annotation_df$type <- factor(type_annotation_df$type,
    levels=unique(type_annotation_df$type))
  print(paste0("Are type_annotation_df and heatmap rownames identical? ",
    identical(rownames(combined_data$heatmap), rownames(type_annotation_df))))
  # determine colours and create annotation:
  names(type_cols) <- unique(combined_data$metadata$type)
  type_annotation <- Heatmap(
    as.matrix(type_annotation_df), 
    col = type_cols, 
    name = "type_annotation", 
    width = unit(4, "mm"), 
    show_row_names = F, show_column_names = F,
    show_heatmap_legend = FALSE,
    row_gap = unit(1, "cm"),
    row_title = NULL
  )

  # create subcluster annotation:
  subcluster_annotation_df <- subset(combined_data$metadata, select = subcluster_id)
  subcluster_annotation_df$subcluster_id <- factor(subcluster_annotation_df$subcluster_id,
    levels=unique(subcluster_annotation_df$subcluster_id))
  print(paste0("Are subcluster_annotation_df and heatmap rownames identical? ",
    identical(rownames(combined_data$heatmap), rownames(subcluster_annotation_df))))
  # determine colours and create annotation:
  subcluster_annot_cols <- subcluster_cols[1:length(unique(subcluster_annotation_df$subcluster_id))]
  names(subcluster_annot_cols) <- unique(combined_data$metadata$subcluster_id)
  subcluster_annotation <- Heatmap(
    as.matrix(subcluster_annotation_df), 
    col = subcluster_annot_cols, 
    name = "subcluster_annotation", 
    width = unit(4, "mm"), 
    show_row_names = F, show_column_names = F,
    show_heatmap_legend = FALSE,
    row_gap = unit(1, "cm"),
    row_title = NULL
  )
  
  # create QC and GIN annotations:
  GIN_values <- rescale(combined_data$metadata$CNA_value, c(0,1))
  GIN_annotation <- rowAnnotation(
    correlation_annotation = anno_barplot(
      GIN_values, name = "GIN",
      gp = gpar(
        col = "#D95F02", 
        width = unit(4, "cm")
      ), 
      border = FALSE, 
      which = "row", 
      axis = F
    ), show_annotation_name = FALSE
  )
  GIN_annotation@name <- "GIN"
  nUMI_annotation <- rowAnnotation(
    correlation_annotation = anno_barplot(
      combined_data$metadata$nUMI,
      gp = gpar(
        col = "#D8B72E", 
        width = unit(4, "cm")
      ), 
      border = FALSE, 
      which = "row", 
      axis = F
    ), show_annotation_name = FALSE
  )
  nUMI_annotation@name <- "nUMI"
  nGene_annotation <- rowAnnotation(
    correlation_annotation = anno_barplot(
      combined_data$metadata$nGene, name = "nGene",
      gp = gpar(
        col = "#9ECAE1", 
        width = unit(4, "cm")
      ), 
      border = FALSE, 
      which = "row", 
      axis = F
    ), show_annotation_name = FALSE
  )
  nGene_annotation@name <- "nGene"
  
  # determine co-ordinates of vertical lines at chromosome borders:
  chr_data <- fetch_chromosome_boundaries(combined_data$heatmap, ref_dir)

  return(
    list(
      type = type_annotation,
      subcluster = subcluster_annotation,
      GIN = GIN_annotation,
      nUMI = nUMI_annotation,
      nGene = nGene_annotation,
      chr_data = chr_data,
      subcluster_cols = subcluster_annot_cols
    )
  )

}
