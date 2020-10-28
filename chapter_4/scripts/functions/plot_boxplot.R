plot_boxplot <- function(
  CNV_data,
  data_type,
  subtype_cols,
  subtype_df,
  plot_dir,
  lib_loc
) {

  library(ggpubr, lib.loc = lib_loc)

  # change names of lists to sample ids and arrange into df:
  CNV_df <- data.frame(
    sample = gsub("\\..*$", "", names(CNV_data)),
    value = CNV_data
  )
  plot_df <- merge(CNV_df, subtype_df, by = "sample")
  plot_df$log2_value <- log2(plot_df$value)

  my_comparisons = list( c("ER", "HER2"), c("ER", "TNBC"), c("HER2", "TNBC") )

  p <- ggboxplot(
    plot_df, 
    x = "subtype", 
    y = "value",
    fill = "subtype",
    palette = subtype_cols,
    bxp.errorbar = TRUE
  )
  p <- p + xlab("Subtype")
  p <- p + ylab(data_type)
  p <- p + theme(
    legend.position = "none"
  )
  p <- p + stat_compare_means(
    comparisons = my_comparisons,
    method = "t.test",
    label = "p.signif"
  )
  p <- p + theme(
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
  )

  pdf(paste0(plot_dir, data_type, "_boxplot_subtypes.pdf"))
    print(p)
  dev.off()

  p <- ggboxplot(
    plot_df, 
    x = "subtype", 
    y = "log2_value",
    fill = "subtype",
    palette = subtype_cols,
    bxp.errorbar = TRUE
  )
  p <- p + xlab("Subtype")
  p <- p + ylab(data_type)
  p <- p + theme(
    legend.position = "none"
  )
  p <- p + stat_compare_means(
    comparisons = my_comparisons,
    method = "t.test",
    label = "p.signif"
  )
  p <- p + theme(
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
  )

  png(paste0(plot_dir, "log2_", data_type, "_boxplot_subtypes.pdf"))
    print(p)
  dev.off()

}



