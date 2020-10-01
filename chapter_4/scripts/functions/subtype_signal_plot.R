subtype_signal_plot <- function(
  plot_df,
  neutral_score,
  CNV_indices,
  include_lengths = FALSE,
  annotate_genes = NULL
) {

  average_signal <- round(apply(plot_df, 2, mean), 6)

  # create area plot presenting InferCNV signal:
  area_df <- data.frame(
    index = seq_along(average_signal),
    average_score = average_signal-neutral_score,
    type = "neutral",
    stringsAsFactors = F
  )
  
  # label gains and losses:
  for (r in 1:nrow(CNV_indices)) {
    area_df$type[
      CNV_indices$start[r]:CNV_indices$end[r]
    ] <- as.character(CNV_indices$call[r])
  }
  area_df$type[
    area_df$type == "neutral" & area_df$average_score > 0
  ] <- "gain_artefact"
  area_df$type[
    area_df$type == "neutral" & area_df$average_score < 0
  ] <- "loss_artefact"
  area_df$type <- factor(
    area_df$type, levels = c(
      "loss_artefact", "gain_artefact", "loss", "gain", "neutral"
    )
  )
  # define colours:
  cols <- c("#D6EAE8", "#F4E9E9", "#76C1C1", "#F7B7B5", "black")
  
  # prepare label data:
  CNV_only <- CNV_indices[CNV_indices$call != "neutral",]
  CNV_only <- CNV_only[grep("artefact", CNV_only$call, invert=TRUE),]
  CNV_only$length_labels <- CNV_only$genomic_length
  # stagger labels:
  CNV_only$length_labels[c(FALSE, TRUE)] <- 
    paste0("\n", CNV_only$length_labels[c(FALSE, TRUE)])
  CNV_only$midpoints <- CNV_only$midpoints + 10
  for (r in 2:nrow(CNV_only)) {
    if (CNV_only$start[r] - CNV_only$end[r-1] < 40) {
      CNV_only$midpoints[r] <- CNV_only$midpoints[r]+30
    }
  }

  # plot on barplot:
  p <- ggplot(area_df, aes(x=index, y=average_score, fill = type))
  p <- p + geom_bar(stat="identity")
  p <- p + scale_fill_manual(values = cols)
  if (include_lengths) {
    p <- p + scale_x_continuous(
      label = CNV_only$length_labels,
      breaks = CNV_only$midpoints,
      expand = c(0,0)
    )
  } else if (!is.null(annotate_genes)) {
    p <- p + scale_x_continuous(
      label = annotate_genes$gene,
      breaks = annotate_genes$index,
      expand = c(0,0)
    )
  } else {
    p <- p + scale_x_continuous(
      expand = c(0,0)
    )
  }
  p <- p + scale_y_continuous(
    limits = c(-0.2, 0.2)
  )
  p <- p + theme_cowplot(12)
  if (include_lengths) {
    p <- p + theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size=15),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  } else if (!is.null(annotate_genes)) {
     p <- p + theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    )
  } else {
    p <- p + theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  }
  
  #p <- p + ylab("Mean CNV signal")
  #p <- p + xlab("Length (kb)")
  for (c_end in chr_data$ends) {
    p <- p + geom_vline(xintercept=c_end)
  }
  # create 0 line:
  p <- p + geom_segment(
    x=CNV_indices$start[1],
    xend=CNV_indices$end[
      nrow(CNV_indices)
    ],
    y=0,
    yend=0
  )

  # convert barplot to grid object:
  signal_plot <- ggplotGrob(p)

  return(signal_plot)

}