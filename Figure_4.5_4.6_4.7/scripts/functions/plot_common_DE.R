plot_common_DE <- function(
  DE_plot_data,
  data_type
) {

  # define heatmap and significance star colours:
  na_less_vector <- unlist(DE_plot_data$avg_logFC)
  na_less_vector <- na_less_vector[!is.na(na_less_vector)]

  if (length(grep("down", data_type)) > 0) {
    heatmap_cols <- colorRamp2(c(min(na_less_vector), 0), 
      c("#0D3084", "white"), space = "sRGB")
    star_col <- "#870C0C"
  } else if (length(grep("up", data_type)) > 0) {
    heatmap_cols <- colorRamp2(c(0, max(na_less_vector)), 
      c("white", "#870C0C"), space = "sRGB")
    star_col <- "#1F0B99"
  }

  mat <- as.matrix(DE_plot_data$avg_logFC)
  rownames(mat) <- gsub("_.*$", "", rownames(mat))
  pdots <- DE_plot_data$p_val_dots
  # order subtype_df the same as matrix column:
  m <- match(colnames(mat), subtype_df$sample)
  subtype_df <- subtype_df[m,]
  # create subtype annotation:
  subtype_annot <- HeatmapAnnotation(
    subtype = as.character(subtype_df$subtype),
    col = list(
      subtype = c(
        subtype_cols["ER"], 
        subtype_cols["HER2"], 
        subtype_cols["TNBC"]
      )
    ),
    show_legend = F,
    show_annotation_name = F
  )
  # create heatmap legend:
  if (length(grep("up", data_type)) > 0) {
    lgd <- Legend(
      at = seq(
        round(min(na_less_vector)), 
        round(max(na_less_vector))
      ),
      col_fun = heatmap_cols, 
      title = "logFC", 
      direction = "vertical",
      grid_height = unit(1, "cm"),
      grid_width = unit(1, "cm"),
      #legend_height = unit(4.5, "cm"),
      #legend_width = unit(2, "cm"),
      labels_gp = gpar(fontsize = 20),
      title_gp = gpar(fontsize = 22, fontface = "plain")
    )
  } else if (length(grep("down", data_type)) > 0) {
   lgd <- Legend(
      at = seq(
        round(min(na_less_vector)), 
        round(max(na_less_vector))
      ),
      col_fun = heatmap_cols, 
      title = "logFC", 
      direction = "vertical",
      grid_height = unit(1, "cm"),
      grid_width = unit(1, "cm"),
      #legend_height = unit(4.5, "cm"),
      #legend_width = unit(2, "cm"),
      labels_gp = gpar(fontsize = 20),
      title_gp = gpar(fontsize = 22, fontface = "plain")
    )
  }
   
  # create main heatmap:
  final_heatmap <- Heatmap(
    mat, 
    na_col = "white",
    name = paste0("hm"), 
    col = heatmap_cols,
    cluster_columns = F, cluster_rows = F,
    show_row_names = T, show_column_names = T,
    row_names_gp = gpar(fontsize = 17),
    column_names_gp = gpar(fontsize = 17),
    show_row_dend = T,
    show_heatmap_legend = F,
    use_raster = T, raster_device = c("png"),
    cell_fun = function(j, i, x, y, w, h, heatmap_cols) { # add text to each grid
        grid.text(pdots[i, j], x, y, gp=gpar(fontsize=30, col = star_col))
    },
    bottom_annotation = subtype_annot
  )
  # make grid object:
  ht_list <- final_heatmap
  annotated_heatmap <- grid.grabExpr(
    draw(ht_list
      #, gap = unit(6, "mm"), heatmap_legend_side = "left"
    )
  )
  dev.off()
  
  if (data_type == "down_all") {
    pdf(paste0(plot_dir, "downreg_gene_heatmap.pdf"), 
      height = 25, width = 20) 
  } else if (data_type == "up_all") {
    pdf(paste0(plot_dir, "upreg_gene_heatmap.pdf"), 
      height = 25, width = 20) 
  } else if (data_type == "down_CNA") {
    pdf(paste0(plot_dir, "downreg_CNA_assoc_gene_heatmap.pdf"), 
      height = 8, width = 25) 
  } else if (data_type == "up_CNA") {
      pdf(paste0(plot_dir, "upreg_CNA_assoc_gene_heatmap.pdf"), 
        height = 8, width = 25) 
  }
    grid.newpage()
  
      # plot heatmap:
      pushViewport(viewport(x = 0.08, y = 0.04, width = 0.9, height = 0.95, 
        just = c("left", "bottom")))
        grid.draw(annotated_heatmap)
      popViewport()

      # plot heatmap legend:
      pushViewport(viewport(x = unit(2.5, "cm"), y = unit(15, "cm"), width = unit(1, "cm"), 
        height = unit(2, "cm")))
        draw(lgd, x = unit(0.1, "cm"), y = unit(0.1, "cm"))
      popViewport()

#      # plot NA legend square:
#      pushViewport(viewport(x = unit(2, "cm"), y = unit(33, "cm"), 
#        width = unit(3.3, "cm"), height = unit(1, "cm")))
#        #grid.rect()
#        grid.rect(width = unit(1, "cm"), height = unit(1, "cm"),
#          just = "right", gp=gpar(col = "#C5E0BA", fill = "#C5E0BA"))
#        if (length(grep("down", data_type)) > 0) {
#          grid.text("No\ndownreg.", x=0.75, gp=gpar(fontsize=20))
#        } else if (length(grep("up", data_type)) > 0) {
#          grid.text("No\nupreg.", x=0.75, gp=gpar(fontsize=20))
#        }
#      popViewport()

      # print subtype labels:
      pushViewport(viewport(x = unit(14.5, "cm"), y = unit(0.5, "cm"), width = unit(3, "cm"), 
        height = unit(1, "cm")))
        grid.text("ER+", gp=gpar(fontsize=25, fontface = "bold", col = subtype_cols[1]))
      popViewport()
      pushViewport(viewport(x = unit(33, "cm"), y = unit(0.5, "cm"), width = unit(3, "cm"), 
        height = unit(1, "cm")))
        grid.text("HER2+", gp=gpar(fontsize=25, fontface = "bold", col = subtype_cols[2]))
      popViewport()
      pushViewport(viewport(x = unit(51, "cm"), y = unit(0.5, "cm"), width = unit(3, "cm"), 
        height = unit(1, "cm")))
        grid.text("TNBC", gp=gpar(fontsize=25, fontface = "bold", col = subtype_cols[3]))
      popViewport()

  dev.off()
  
}