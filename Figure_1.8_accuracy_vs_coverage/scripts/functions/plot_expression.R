plot_expression <- function(
  log_fc_df, 
  chromosome_midpoints, 
  chromosome_ends,
  y_lims = c(-2.5, 2),
  CNV_record = "none",
  line_only = FALSE
) { 

  p <- ggplot(log_fc_df, aes(x=number, y=count))
  if (!line_only) {
    p <- p + geom_point(colour = "#E8D172")
  }
  p <- p + xlab("Genomic location")
  p <- p + scale_x_continuous(
    breaks = unlist(chromosome_midpoints),
    labels = paste0("chr", 1:length(chromosome_midpoints)),
    expand = c(0, 0),
    limits = c(min(log_fc_df$number), max(log_fc_df$number))
  )
  p <- p + ylab("Log10 fold change from median")
  p <- p + scale_y_continuous(
    breaks = seq(y_lims[1], y_lims[2]),
    labels = as.character(seq(y_lims[1], y_lims[2])),
    limits = y_lims
  )
  for (end in chromosome_ends) {
    p <- p + geom_vline(xintercept=end)
  }
  if (CNV_record == "none") {
    p <- p + geom_segment(
      x=0, 
      xend=max(unlist(chromosome_ends)), 
      y=log_median_original_fold_change, 
      yend=log_median_original_fold_change, 
      size=1, color="red"
    )
  } else {
    for (r in 1:nrow(CNV_record)) {
      # if CNV is differentiating between subpops, highlight in yellow:
      if ("differentiating" %in% colnames(CNV_record)) {
        if (CNV_record$differentiating[r]) {
          vert_line_col <- "#37841f"
          hor_line_col <- "#37841f"
        } else {
          vert_line_col <- "red"
          hor_line_col <- "red"
        }
        if (length(CNV_record$differentiating[r-1]) > 0){
          if (CNV_record$differentiating[r-1]) {
            vert_line_col <- "#37841f"
            hor_line_col <- "red"
          }
        }
      } else {
        vert_line_col <- "red"
        hor_line_col <- "red"
      }
      # create horizontal line:
      p <- p + geom_segment(
        x=CNV_record$start[r], 
        xend=CNV_record$end[r], 
        y=CNV_record$log_median_modified_FC[r], 
        yend=CNV_record$log_median_modified_FC[r], 
        size=1, color=hor_line_col
      )
      # create left vertical line:
      if (r != 1) {
        p <- p + geom_segment(
          x=CNV_record$start[r], 
          xend=CNV_record$start[r], 
          y=CNV_record$log_median_modified_FC[r-1], 
          yend=CNV_record$log_median_modified_FC[r], 
          size=1, color=vert_line_col
        )
      }
    }
  }
  
  return(p)
}