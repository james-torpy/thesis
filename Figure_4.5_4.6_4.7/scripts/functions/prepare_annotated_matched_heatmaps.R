
prepare_annotated_matched_heatmaps <- function(
  combined_data,
  combined_annotations
) {

  # scale all samples the same:
  print(paste0("Summary of heatmap df gene before rescaling to min/max values: "))
  print(
    lapply(combined_data$heatmap, function(x) {
      summary(unlist(x))
    })[[1]]
  )
  
  combined_data$rescaled_heatmap <- as.data.frame(
    rescale(
      as.matrix(combined_data$heatmap), 
      c(-1, 1)
    )
  )
  
  print(paste0("Summary of heatmap df gene after rescaling: "))
  print(
    lapply(combined_data$rescaled_heatmap, function(x) {
      summary(unlist(x))
    })[[1]]
  )

  # prepare df for plotting:
  plot_object <- combined_data$rescaled_heatmap
  colnames(plot_object) <- rep("la", ncol(plot_object))
  
  # define heatmap colours:
  na_less_vector <- unlist(plot_object)
  na_less_vector <- na_less_vector[!is.na(na_less_vector)]
  heatmap_cols <- colorRamp2(c(min(na_less_vector), 0, max(na_less_vector)), 
        c("#00106B", "white", "#680700"), space = "sRGB")
  
  print("Generating final heatmap...")
  
  final_heatmap <- Heatmap( 
    as.matrix(plot_object),
    name = paste0("hm"),
    na_col = na_colour,
    col = heatmap_cols,
    cluster_columns = F, cluster_rows = F,
    show_row_names = F, show_column_names = T,
    column_names_gp = gpar(col = "white"),
    show_row_dend = FALSE,
    show_heatmap_legend = FALSE,
    use_raster = T, raster_device = c("png"),
    row_gap = unit(1, "cm"),
    row_title = NULL
  )
  
  # determine co-ordinates of horizontal lines at subpop borders:
  combined_data$metadata$subcluster_id <- as.character(
    combined_data$metadata$subcluster_id
  )
  hlines <- 1 - (
    cumsum(
      rle(combined_data$metadata$subcluster_id)$lengths/
        nrow(combined_data$metadata)
    )
  )

  ht_list <- combined_annotations$type
  ht_list <- ht_list + combined_annotations$subcluster
  ht_list <- ht_list + final_heatmap
  ht_list <- ht_list + combined_annotations$GIN
  ht_list <- ht_list + combined_annotations$nUMI + combined_annotations$nGene
  
  annotated_heatmap <- grid.grabExpr(
    draw(ht_list, gap = unit(3, "mm"), heatmap_legend_side = "left")
  )
  dev.off()

  # determine where starting co-ordinates for heatmap are based upon longest cluster name
  # (0.00604 units per character):
  longest_cluster_name <- max(nchar(unique(as.character(combined_data$metadata$sample))))
  x_coord <- longest_cluster_name*0.0037
  
  # generate heatmap legend:
  signal_ranges <- round(range(unlist(plot_object), na.rm = TRUE), 2)
  lgd <- Legend(
    at = c(signal_ranges[1], 0, signal_ranges[2]),
    col_fun = heatmap_cols, 
    title = "Scaled CNA\nvalue", 
    direction = "horizontal",
    grid_height = unit(3, "cm"),
    grid_width = unit(4.5, "cm"),
    legend_height = unit(3, "cm"),
    legend_width = unit(4.5, "cm"),
    labels_gp = gpar(fontsize = 20),
    title_gp = gpar(fontsize = 22, fontface = "plain")
  )

  return(
    list(
      heatmap = annotated_heatmap,
      hlines = hlines,
      x_coord = x_coord,
      legend = lgd
    )
  )

}
