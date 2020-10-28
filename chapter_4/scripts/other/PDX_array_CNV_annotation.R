  pdx_name <- gsub("CID", "PDX", sample_name)
  PDX_CNVs <- all_array_CNVs[,colnames(all_array_CNVs) %in% pdx_name]
  names(PDX_CNVs) <- rownames(all_array_CNVs)

  PDX_CNVs[PDX_CNVs > 6] <- 6
  PDX_CNVs <- PDX_CNVs/2
  PDX_CNVs <- PDX_CNVs[colnames(heatmap_df)]
  PDX_CNVs[is.na(PDX_CNVs)] <- 2
  names(PDX_CNVs) <- colnames(heatmap_df)
  PDX_CNVs <- as.data.frame(PDX_CNVs)
  PDX_CNVs <- as.data.frame(t(PDX_CNVs))

  saveRDS(PDX_CNVs, paste0(Robject_dir, "PDX_CNVs.RData"))

  PDX_heatmap <- Heatmap(
  	PDX_CNVs, name=paste0("array_heatmap"),
  	col = colorRamp2(c(min(PDX_CNVs[1,]), 1, max(PDX_CNVs[1,])), 
      c("#00106B", "white", "#680700"), space = "sRGB"),
  	  cluster_columns = F, cluster_rows = F, show_row_dend = FALSE,
  	  show_row_names = F, show_column_names = F,
  	  heatmap_legend_param = list(title = "SNP array\nPDX CNV", legend_direction = "horizontal"),
  	  use_raster = T, raster_device = c("png")
  )
  grid_PDX_heatmap <- grid.grabExpr(draw(PDX_heatmap, heatmap_legend_side = "left"))