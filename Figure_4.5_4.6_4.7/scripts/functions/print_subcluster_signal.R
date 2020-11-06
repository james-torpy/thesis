print_subcluster_signal <- function(
  s_plots,
  chr_coords,
  plot_dir,
  annotate_genes = NULL,
  no_genes = NULL,
  include_lengths = FALSE,
  with_artefacts = FALSE
) {

  library(grid)

  if (include_lengths) {
    filename_pre <- paste0("CNV_signal_plots_with_lengths")
  } else if (!is.null(annotate_genes)) {
    filename_pre <- paste0("CNV_signal_plots_with_DE_genes")
  } else {
    filename_pre <- paste0("CNV_signal_plots") 
  }

  if (with_artefacts) {
    filename_full <- paste0(filename_pre, "_with_artefacts.pdf")
  } else {
    filename_full <- paste0(filename_pre, ".pdf")
  }
 
  
  pdf(
    paste0(plot_dir, filename_full), 
    height = 14, 
    width = 20
  )   

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

      } else if (!is.null(annotate_genes)) {

        if (p == length(s_plots)) {

          # plot:
          pushViewport(viewport(
            x = 0.52, y = 0.91-(0.15*(p-1)), 
            width = 0.94, height = 0.13
          ))
            grid.draw(s_plots[[p]])
          popViewport()
      
          # title:
          pushViewport(viewport(
            x = 0.03, y = 0.937-(0.15*(p-1)), 
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
            x = 0.52, y = 0.91-(0.15*(p-1)), 
            width = 0.94, height = 0.13
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

        }
        
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

    # annotate genes if needed:
    if (!is.null(annotate_genes)) {

      pushViewport(viewport(
        x = 0.521, y = 0.06, 
        width = 0.931, height = 0.07
      ))
        grid.rect()

      for (j in 1:nrow(annotate_genes)) {
        # adjust x and y for stagger value:
        pushViewport(viewport(
          x = (annotate_genes$index[j]/no_genes), 
          y = 1 - (0.2*annotate_genes$stagger[j]), 
          width = 0.01, 
          height = 1
        ))          
          if (annotate_genes$direction[j] == "up") {
            grid.text(annotate_genes$gene[j], gp=gpar(fontsize=16, col="red"))
          } else {
            grid.text(annotate_genes$gene[j], gp=gpar(fontsize=16, col="blue"))
          }
        popViewport()
      }
        
      popViewport()

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
      pushViewport(viewport(x = 0.5, y = 0.05, width = 0.05, height = 0.05, 
        just = c("left", "bottom")))
        grid.text("CNV length (kb)", gp=gpar(fontsize=20))
      popViewport()
    } else if (!is.null(annotate_genes)) {
      pushViewport(viewport(x = 0.52, y = 0.035, width = 0.05, height = 0.05, 
        just = c("top")))
        grid.text("DE CNA assoc. genes", gp=gpar(fontsize=20))
      popViewport()
    }
   
  dev.off()
}