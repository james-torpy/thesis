
create_signal_plot <- function(
  all_indices, 
  plot_df, 
  log_mod_fc_df,
  peak_data,
  chr_coords, 
  accuracy_df,
  accuracy_annotation = TRUE,
  copy_no_estimates = FALSE,
  func_dir,
  out_dir
) {

#  all_indices <- CNV_indices
#  plot_df <- area_df
#  log_mod_fc_df <- log_modified_fold_change_df
#  peak_data <- peak_estimates
#  chr_coords <- chr_data
#  accuracy_df <- CNV_accuracy_df
#  out_dir <- plot_dir
#  accuracy_annotation = FALSE
#  copy_no_estimates = TRUE
  
  # load legend functions:
  true_pos_neg_false_pos_neg_legend <- dget(paste0(func_dir, 
    "true_pos_neg_false_pos_neg_legend.R"))
  
  true_pos_false_pos_wrong_legend <- dget(paste0(func_dir,
    "true_pos_false_pos_wrong_legend.R"))

  if (accuracy_annotation) {
    print("Plotting CNV signal with accuracy annotation...")
  } else {
    print("Plotting CNV signal without accuracy annotation...")
  }

  # define colours:
  cols <- c("#F7B7B5", "#76C1C1", "black")
  all_indices$type[all_indices$multiplier > 1] <- "gain"
  all_indices$type[all_indices$multiplier == 1] <- "neutral"
  all_indices$type[all_indices$multiplier < 1] <- "loss"
  
  scaled_CNV_indices <- data.frame(
    start = all_indices$start,
    end = all_indices$end,
    type = all_indices$type,
    multiplier = all_indices$multiplier
  )

  # 0 = 0 alleles = -0.08
  # 0.5 = 1 allele = -0.02
  # 1 = 2 alleles = 0
  # 1.5 = 3 alleles = 0.02
  # 2 = 4 alleles = 0.04
  # 3 = 6 alleles = 0.08
  scaled_CNV_indices$multiplier[
    scaled_CNV_indices$multiplier == 0
  ] <- -0.08
  scaled_CNV_indices$multiplier[
    scaled_CNV_indices$multiplier == 0.5
  ] <- -0.02
  scaled_CNV_indices$multiplier[
    scaled_CNV_indices$multiplier == 1
  ] <- 0
  scaled_CNV_indices$multiplier[
    scaled_CNV_indices$multiplier == 1.5
  ] <- 0.02
  scaled_CNV_indices$multiplier[
    scaled_CNV_indices$multiplier == 2
  ] <- 0.04
  scaled_CNV_indices$multiplier[
    scaled_CNV_indices$multiplier == 3
  ] <- 0.08

  # expand CNVs less than 15 genes long:
  scaled_CNV_indices$length <- scaled_CNV_indices$end-scaled_CNV_indices$start
  scaled_CNV_indices$midpoints = scaled_CNV_indices$start + floor(scaled_CNV_indices$length/2)
  scaled_CNV_indices$keep = TRUE
  
  for (r in 1:nrow(scaled_CNV_indices)) {
    if (scaled_CNV_indices$length[r] < 40 & scaled_CNV_indices$multiplier[r] != 0) {
  
      if (scaled_CNV_indices$start[r] > 20) {
        scaled_CNV_indices$start[r] <- scaled_CNV_indices$midpoint[r]-20
        scaled_CNV_indices$end[r-1] <- scaled_CNV_indices$start[r]+1
      } else {
        scaled_CNV_indices$start[r] <- 1
        if (r != 1) {
          scaled_CNV_indices$keep[1] <- FALSE
        }
      }
  
      scaled_CNV_indices$end[r] <- scaled_CNV_indices$midpoints[r]+20
      scaled_CNV_indices$start[r+1] <- scaled_CNV_indices$end[r]+1
  
    }
  }

  # stagger CNV labels for plotting:
  stag_lab <- CNV_indices$length[CNV_indices$ticks == "include"]
  stag_lab[c(FALSE, TRUE)] <- paste0("\n", stag_lab[c(FALSE, TRUE)])
  # plot on barplot:
  plot_df$type <- factor(plot_df$type, levels = c("gain", "loss", "neutral"))
  p <- ggplot(plot_df, aes(x=index, y=average_score, fill = type))
  p <- p + geom_bar(stat="identity")
  p <- p + scale_fill_manual(values = cols)
  p <- p + scale_x_continuous(
    limits = c(
      0,length(log_mod_fc_df$count)
    ), 
    expand = c(0, 0),
    breaks = CNV_indices$midpoints[CNV_indices$ticks == "include"],
    labels = stag_lab
  )
  p <- p + scale_y_continuous(
    limits = c(-0.097, 0.094),
    breaks = c(-0.05, 0, 0.05),
    labels = c("-0.05", "0", "0.05"),
    sec.axis = sec_axis(
      ~., 
      "Copy number\nfold change", 
      breaks = c(-0.08, 0, 0.08),
      labels = c("Total\nloss", "1", "3")
    )
  )
  p <- p + theme_cowplot(12)
  if (accuracy_annotation) {

    p <- p + theme(
      axis.title.x = element_text(size=30, margin = margin(t = 20, r = 0, b = 0, l = 0)),
      axis.text.x = element_text(size=23, margin = margin(t = 70, r = 0, b = 0, l = 0)),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size=30),
      axis.title.y = element_text(size=30, margin = margin(t = 0, r = 30, b = 0, l = 0)),
      axis.title.y.right = element_text(size=30, margin = margin(t = 0, r = 0, b = 0, l = 0)),
      legend.position = "none"
    )

  } else {

    if (copy_no_estimates) {

      p <- p + theme(
        axis.title.x = element_text(size=30, margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=30),
        axis.title.y = element_text(size=30, margin = margin(t = 0, r = 30, b = 0, l = 0)),
        axis.title.y.right = element_text(size=30, margin = margin(t = 0, r = 0, b = 0, l = 0)),
        legend.position = "none"
      )

      # add gain copy numbr estimates:
      for (m in 1:nrow(peak_data$gain)) {
        p <- p + annotate(
          geom = "text", 
          x = peak_data$gain$midpoint[m], 
          y = 0.094, 
          label = peak_data$gain$estimate_lab[m],
          color = "#BF3667", 
          size = 9, 
          fontface = "bold"
        )
      }

      # add loss copy numbr estimates:
      for (m in 1:nrow(peak_data$loss)) {
        p <- p + annotate(
          geom = "text", 
          x = peak_data$loss$midpoint[m], 
          y = -0.086, 
          label = peak_data$loss$estimate_lab[m],
          color = "#58B9DB", 
          size = 9,
          fontface = "bold"
        )
      }

    } else {

      p <- p + theme(
        axis.title.x = element_text(size=30, margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=23, margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=30),
        axis.title.y = element_text(size=30, margin = margin(t = 0, r = 30, b = 0, l = 0)),
        axis.title.y.right = element_text(size=30, margin = margin(t = 0, r = 0, b = 0, l = 0)),
        legend.position = "none"
      )

    }

  }
  p <- p + ylab("Mean CNV signal")
  p <- p + xlab("CNV length (genes)")
  for (c_end in chr_coords$ends) {
    p <- p + geom_vline(xintercept=c_end)
  }
  # create 0 line:
  p <- p + geom_segment(
    x=scaled_CNV_indices$start[1],
    xend=scaled_CNV_indices$end[
      nrow(scaled_CNV_indices)
    ],
    y=0,
    yend=0
  )
  for (r in 1:nrow(scaled_CNV_indices)) {
    # create horizontal line:
    p <- p + geom_segment(
      x=scaled_CNV_indices$start[r], 
      xend=scaled_CNV_indices$end[r], 
      y=scaled_CNV_indices$multiplier[r], 
      yend=scaled_CNV_indices$multiplier[r], 
      size=0.75, color="#430F82"
    )
  
    # create left vertical line:
    if (r != 1) {
      p <- p + geom_segment(
        x=scaled_CNV_indices$start[r], 
        xend=scaled_CNV_indices$start[r], 
        y=scaled_CNV_indices$multiplier[r-1], 
        yend=scaled_CNV_indices$multiplier[r], 
        size=0.75, color="#430F82"
      )
    }
  }
  # convert barplot to grid object:
  signal_plot <- ggplotGrob(p)
  dev.off()

#######
#
#  png(
#      paste0(out_dir, "signal_plot.png"), 
#      height = 8, 
#      width = 22,
#      res = 300,
#      units = "in"
#    )
#
#    grid.newpage()
#
#      # draw signal plot:
#      pushViewport(viewport(x = 0.06, y = 0.001, width = 0.93, height = 0.95, 
#        just = c("left", "bottom")))
#        grid.draw(signal_plot)
#      popViewport()
#
#  dev.off()
#
#######

  # determine accuracy calls present:
  all_accuracy <- unique(accuracy_df$accuracy_call)
  all_accuracy <- paste0(all_accuracy[order(all_accuracy)], collapse = "_")

  if (accuracy_annotation) {

    print(paste0(
      "Writing plot to file: ", 
      out_dir, 
      "signal_vs_simulated_CNV_plot.png"
    ))
  
    # plot average signal:
    png(
      paste0(out_dir, "signal_vs_simulated_CNV_plot.png"), 
      height = 8, 
      width = 22,
      res = 300,
      units = "in"
    )   

  } else {

    if (copy_no_estimates) {

      print(paste0(
        "Writing plot to file: ", 
        out_dir, 
        "signal_vs_simulated_CNV_plot_with_copy_no_estimates.png"
      ))
    
      # plot average signal:
      png(
        paste0(out_dir, "signal_vs_simulated_CNV_plot_with_copy_no_estimates2.png"), 
        height = 8, 
        width = 22,
        res = 300,
        units = "in"
      )  

    } else {

      print(paste0(
        "Writing plot to file: ", 
        out_dir, 
        "signal_vs_simulated_CNV_plot_no_accuracy_annotation.png"
      ))
    
      # plot average signal:
      png(
        paste0(out_dir, "signal_vs_simulated_CNV_plot_no_accuracy_annotation.png"), 
        height = 8, 
        width = 22,
        res = 300,
        units = "in"
      )  

    }
    
  }

    grid.newpage()

      if (copy_no_estimates) {

        # draw signal plot:
        pushViewport(viewport(x = 0.02, y = 0.001, width = 0.92, height = 0.95, 
          just = c("left", "bottom")))
          grid.draw(signal_plot)
        popViewport()

        # label chromosomes:
        for ( e in 1:length(chr_coords$lab_pos) ) {
          pushViewport(
            viewport(
              x = 0.12 + chr_coords$lab_pos[e]/1.315, 
              y = 0.94, width = 0.05, height = 0.05, 
              just = c("left", "bottom")
            )
          )
            if (e==1) {
              grid.text(
                paste0(
                  "\n", names(chr_coords$lab_pos)[e]
                ), gp=gpar(fontsize=22))
            } else if (e==21) {
              grid.text(
                gsub("chr", "", names(chr_coords$lab_pos)[e]), 
                gp=gpar(fontsize=22)
              )
            } else {
              grid.text(
                paste0(
                  "\n", gsub("chr", "", names(chr_coords$lab_pos)[e])
                ), gp=gpar(fontsize=22)
              )
            }
          popViewport()
        }

      } else {

        # draw signal plot:
        pushViewport(viewport(x = 0.06, y = 0.001, width = 0.93, height = 0.9, 
          just = c("left", "bottom")))
          grid.draw(signal_plot)
        popViewport()

        # label chromosomes:
        for ( e in 1:length(chr_coords$lab_pos) ) {
          pushViewport(
            viewport(
              x = 0.12 + chr_coords$lab_pos[e]/1.315, 
              y = 0.91, width = 0.05, height = 0.05, 
              just = c("left", "bottom")
            )
          )
            if (e==1) {
              grid.text(
                paste0(
                  "\n", names(chr_coords$lab_pos)[e]
                ), gp=gpar(fontsize=22))
            } else if (e==21) {
              grid.text(
                gsub("chr", "", names(chr_coords$lab_pos)[e]), 
                gp=gpar(fontsize=22)
              )
            } else {
              grid.text(
                paste0(
                  "\n", gsub("chr", "", names(chr_coords$lab_pos)[e])
                ), gp=gpar(fontsize=22)
              )
            }
          popViewport()
        }

      }

      if (accuracy_annotation) {

        # draw accuracy annotation:
        pushViewport(viewport(x = 0.14, y = 0.17, width = 0.77, height = 0.13, 
          just = c("left", "bottom")))
          grid.draw(accuracy_heatmap_obj)
        popViewport()

        # draw accuracy legend:
        pushViewport(viewport(x = unit(1, "cm"), y = unit(2, "cm"), 
                          width = unit(5.5, "cm"), height = unit(4.5, "cm"), 
                          just = c("left", "bottom")))
          if (all_accuracy == "false_positive_true_positive_wrong_call") {
            true_pos_false_pos_wrong_legend()
          } else if (all_accuracy == "false_negative_false_positive_true_negative_true_positive") {
            true_pos_neg_false_pos_neg_legend()
          }
        popViewport()

      }

      if (copy_no_estimates) {

        # label gain copy number estimates:
        pushViewport(viewport(x = 0.082, y = 0.89, width = 0.05, height = 0.05))
          #grid.rect()
          grid.text("Gain copy\nestimate:", gp=gpar(fontsize=24, fontface="bold", col = "#BF3667"))
        popViewport()
    
#        pushViewport(viewport(x = 0.526, y = 0.867, width = 0.785, height = 0.05))
#          #grid.rect()
#          for (k in 1:nrow(peak_data$gain)) {
#            pushViewport(viewport(x = 0.013 + (peak_data$gain$midpoint[k]/(nrow(area_df)+10)), width = 0.01, height = 0.8))
#              #grid.rect()
#              grid.text(
#                peak_data$gain$estimate_lab[k], 
#                gp=gpar(fontsize=16, fontface="bold", col = peak_data$gain$estimate_colour[k])
#              )
#            popViewport()
#          }
#        popViewport()
    
        # label loss copy number estimates:
        pushViewport(viewport(x = 0.082, y = 0.18, width = 0.05, height = 0.05))
          #grid.rect()
          grid.text("Loss copy\nestimate:", gp=gpar(fontsize=24, fontface="bold", col = "#58B9DB"))
        popViewport()
    
#        pushViewport(viewport(x = 0.526, y = 0.25, width = 0.785, height = 0.05))
#          #grid.rect()
#          for (k in 1:nrow(peak_data$loss)) {
#            pushViewport(viewport(x = 0.013 + (peak_data$loss$midpoint[k]/(nrow(area_df)+10)), width = 0.01, height = 0.8))
#              #grid.rect()
#              grid.text(
#                peak_data$loss$estimate_lab[k], 
#                gp=gpar(fontsize=16, fontface="bold", col = peak_data$loss$estimate_colour[k])
#              )
#            popViewport()
#          }
#        popViewport()

      }
    
  dev.off()

}