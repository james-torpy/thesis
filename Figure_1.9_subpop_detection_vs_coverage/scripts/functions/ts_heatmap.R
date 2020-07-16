ts_heatmap <- function(temp_heatmap, temp_metadata, plot_dir) {

  # define heatmap colours:
  na_less_vector <- unlist(temp_heatmap)
  na_less_vector <- na_less_vector[!is.na(na_less_vector)]
  temp_cols <- colorRamp2(c(min(na_less_vector), 1, max(na_less_vector)), 
        c("#00106B", "white", "#680700"), space = "sRGB")
  
  # prepare df for plotting:
  temp_object <- temp_heatmap
  colnames(temp_object) <- rep("la", ncol(temp_object))
  
  plot_heatmap <- Heatmap(
    as.matrix(temp_object), name = paste0("hm"), 
    col = temp_cols,
    cluster_columns = F, cluster_rows = F,
    show_row_names = F, show_column_names = T,
    column_names_gp = gpar(col = "white"),
    show_row_dend = F,
    show_heatmap_legend = F,
    use_raster = T, raster_device = c("png")
  )
  
  final_heatmap <- grid.grabExpr(
    draw(plot_heatmap, gap = unit(6, "mm"), heatmap_legend_side = "left")
  )
  dev.off()
  
  # plot final annotated heatmap:
  png(paste0(plot_dir, "ts_heatmap.png"), 
    height = 13, width = 20, res = 300, units = "in") 
  
    grid.newpage()
  
      # plot heatmap:
      pushViewport(viewport(x = 0.155, y = 0.065, width = 0.8, height = 0.85, 
        just = c("left", "bottom")))
        #grid.rect()
        grid.draw(final_heatmap)
        
      popViewport()
      
  dev.off()
  
}