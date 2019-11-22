group_annotation <- Heatmap(
    group_annotation_df, 
    col = cluster_cols, 
    name = "group_annotation", 
    width = unit(4, "mm"), 
    show_row_names = F, show_column_names = F,
    heatmap_legend_param = list(
      title = "Cluster", title_gp = gpar(fontsize = 18, fontface = "bold", lineheight = 2), 
      labels_gp = gpar(fontsize = 18), 
      at = as.character(levels(group_annotation_df$reclustered_cell_type))
    )
  )

test1_heatmap <- grid.grabExpr(
  draw(
  	group_annotation, gap = unit(6, "mm"), 
  	heatmap_legend_side = "left"
  )
)
test1_draw <- draw(
  	group_annotation, gap = unit(6, "mm"), 
  	heatmap_legend_side = "left"
  )

pdf(paste0(plot_dir, "/test1.pdf"))
  grid.newpage()
    pushViewport(viewport(x = 0.005, y = 0.075, width = 0.99, 
    height = 0.77, just = c("left", "bottom")))
    grid.draw(test1_heatmap)
  popViewport()
dev.off()


group_annotation@matrix_param$gap <- unit(5, "mm")

test2_heatmap <- grid.grabExpr(
  draw(
  	group_annotation, gap = unit(6, "mm"), 
  	heatmap_legend_side = "left"
  )
)
 
pdf(paste0(plot_dir, "/test2.pdf"))
  grid.newpage()
    pushViewport(viewport(x = 0.005, y = 0.075, width = 0.99, 
    height = 0.77, just = c("left", "bottom")))
    grid.draw(test2_heatmap)
  popViewport()
dev.off()


group_annotation <- Heatmap(
    group_annotation_df, 
    col = cluster_cols, 
    name = "group_annotation", 
    width = unit(4, "mm"), 
    show_row_names = F, show_column_names = F,
    heatmap_legend_param = list(
      title = "Cluster", title_gp = gpar(fontsize = 18, fontface = "bold", lineheight = 2), 
      labels_gp = gpar(fontsize = 18), 
      at = as.character(levels(group_annotation_df$reclustered_cell_type))
    )
  )

test3_draw <- draw(
  group_annotation, gap = unit(6, "mm"), 
  heatmap_legend_side = "left"
)
test3_draw@heatmap_legend_param$size <- unit(4, "cm")
test3_heatmap <- grid.grabExpr(test3_draw)

pdf(paste0(plot_dir, "/test3.pdf"))
  grid.newpage()
    pushViewport(viewport(x = 0.005, y = 0.075, width = 0.99, 
    height = 0.77, just = c("left", "bottom")))
    grid.draw(test3_heatmap)
  popViewport()
dev.off()
