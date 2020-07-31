print_subcluster_signal <- function(
  s_plots,
  chr_coords,
  plot_dir,
  include_lengths = FALSE
) {

  if (include_lengths) {
    pdf(
      paste0(plot_dir, "CNV_signal_plots_with_lengths.pdf"), 
      height = 14, 
      width = 20
    )   

  } else {
    pdf(
      paste0(plot_dir, "CNV_signal_plots.pdf"), 
      height = 10, 
      width = 20
    )   
  }
  
  # draw signal plots and titles:
  grid.newpage()

    for (p in 1:length(s_plots)) {
  
      if (include_lengths) {
        # plot:
        pushViewport(viewport(
          x = 0.52, y = 0.91-(0.15*(p-1)), 
          width = 0.94, height = 0.14
        ))
          grid.draw(s_plots[[p]])
        popViewport()
    
        # title:
        pushViewport(viewport(
          x = 0.03, y = 0.927-(0.15*(p-1)), 
          width = 0.1, height = 0.1
        ))
          grid.text(
            gsub("_", " ", names(s_plots)[p])
            , gp=gpar(fontsize=18)
          )
        popViewport()
      } else {
        # plot:
        pushViewport(viewport(
          x = 0.52, y = 0.9-(0.16*(p-1)), 
          width = 0.94, height = 0.14
        ))
          grid.draw(s_plots[[p]])
        popViewport()
    
        # title:
        pushViewport(viewport(
          x = 0.03, y = 0.905-(0.16*(p-1)), 
          width = 0.1, height = 0.1
        ))
          grid.text(
            gsub("_", " ", names(s_plots)[p])
            , gp=gpar(fontsize=18)
          )
        popViewport()
      }  
  
    }

    # draw chromosome labels:
    for ( e in 1:length(chr_coords$lab_pos) ) {
      pushViewport(viewport(x = 0.037 + chr_coords$lab_pos[e]/1.085, y = 0.96, width = 0.05, height = 0.05, 
        just = c("left", "bottom")))
        if (e==1) {
          grid.text(
            paste0(
              "\n",
              names(chr_coords$lab_pos)[e]
            ), gp=gpar(fontsize=18)
          )
        } else if (e==21) {
          grid.text(
            gsub("chr", "", names(chr_coords$lab_pos)[e]), 
              gp=gpar(fontsize=18)
          )
        } else {
          grid.text(
            paste0(
              "\n",
              gsub("chr", "", names(chr_coords$lab_pos)[e])
            ), 
            gp=gpar(fontsize=18)
          )
        }
      popViewport()
    }

    if (include_lengths) {
      # draw x-axis label:
      pushViewport(viewport(x = 0.5, y = 0.03, width = 0.05, height = 0.05, 
        just = c("left", "bottom")))
        grid.text("CNV length (kb)", gp=gpar(fontsize=20))
      popViewport()

    }
   
  dev.off()
}