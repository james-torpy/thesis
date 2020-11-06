subtype_signal_plot <- function(
  plot_df,
  neutral_score,
  CNA_indices,
  annotate_genes = NULL,
  include_artefacts = FALSE,
  include_lengths = FALSE
) {

  # load packages and functions:
  library(cowplot)

  average_signal <- round(apply(plot_df, 2, mean), 6)

  # create area plot presenting InferCNV signal:
  area_df <- data.frame(
    index = seq_along(average_signal),
    average_score = average_signal-neutral_score,
    type = "neutral",
    stringsAsFactors = F
  )
  
  # label gains and losses:
  for (r in 1:nrow(CNA_indices)) {
    area_df$type[
      CNA_indices$start[r]:CNA_indices$end[r]
    ] <- as.character(CNA_indices$call[r])
  }
  area_df$type[
    area_df$type == "neutral" & area_df$average_score > 0
  ] <- "gain_artefact"
  area_df$type[
    area_df$type == "neutral" & area_df$average_score < 0
  ] <- "loss_artefact"

  if (include_artefacts) {
    if ("neutral" %in% area_df$type) {
      # make factor:
      area_df$type <- factor(
        area_df$type,
        levels = c("gain", "neutral", "loss", "gain_artefact", "loss_artefact")
      )
      # define colours:
      cols <- c("#F7B7B5", "black", "#76C1C1", "#F4E9E9", "#D6EAE8")
    } else {
      # make factor:
      area_df$type <- factor(
        area_df$type,
        levels = c("gain", "loss", "gain_artefact", "loss_artefact")
      )
      # define colours:
      cols <- c("#F7B7B5", "#76C1C1", "#F4E9E9", "#D6EAE8")
    }
    
  } else {

    # change all artefact signal to neutral signal:
    area_df$average_score[grep("artefact", area_df$type)] <- 0
    area_df$type[grep("artefact", area_df$type)] <- "neutral"

    if ("neutral" %in% area_df$type) {
      
      # make factor:
      area_df$type <- factor(
        area_df$type,
        levels = c("gain", "neutral", "loss")
      )
      # define colours:
      cols <- c("#F7B7B5", "black", "#76C1C1")

    } else {

      # make factor:
      area_df$type <- factor(
        area_df$type,
        levels = c("gain", "loss")
      )
      # define colours:
      cols <- c("#F7B7B5", "#76C1C1")

    }

  }
  
  # prepare label data:
  CNA_only <- CNA_indices[CNA_indices$call != "neutral",]
  CNA_only <- CNA_only[grep("artefact", CNA_only$call, invert=TRUE),]
  CNA_only$length_labels <- CNA_only$genomic_length
  # stagger labels:
  CNA_only$length_labels[c(FALSE, TRUE)] <- 
    paste0("\n", CNA_only$length_labels[c(FALSE, TRUE)])
  CNA_only$midpoints <- CNA_only$midpoints + 10
  for (r in 2:nrow(CNA_only)) {
    if (CNA_only$start[r] - CNA_only$end[r-1] < 40) {
      CNA_only$midpoints[r] <- CNA_only$midpoints[r]+30
    }
  }

  # plot on barplot:
  p <- ggplot(area_df, aes(x=index, y=average_score, fill = type))
  p <- p + geom_bar(stat="identity")
  p <- p + scale_fill_manual(values = cols)
  if (include_lengths) {
    p <- p + scale_x_continuous(
      label = CNA_only$length_labels,
      breaks = CNA_only$midpoints,
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
  
  #p <- p + ylab("Mean CNA signal")
  #p <- p + xlab("Length (kb)")
  for (c_end in chr_data$ends) {
    p <- p + geom_vline(xintercept=c_end)
  }
  # create 0 line:
  p <- p + geom_segment(
    x=CNA_indices$start[1],
    xend=CNA_indices$end[
      nrow(CNA_indices)
    ],
    y=0,
    yend=0
  )

  # convert barplot to grid object:
  signal_plot <- ggplotGrob(p)

#  ######
#  pdf(
#    paste0(plot_dir, "test.pdf"), 
#    height = 14, 
#    width = 20
#  )   
#
#  # draw signal plots and titles:
#  grid.newpage()
# 
#      
#        # plot:
#        pushViewport(viewport(
#          x = 0.52, y = 0.91, 
#          width = 0.94, height = 0.14
#        ))
#          grid.draw(signal_plot)
#        popViewport()
#  dev.off()
#  ######

  return(signal_plot)

}