true_pos_neg_false_pos_neg_legend <- function() {
  # plot legend text:
  pushViewport(viewport(x = 0.48, y = 0.87, 
                        width = unit(1, "cm"), height = unit(0.5, "cm"), 
                        just = c("left")))
      grid.text("true positive", gp=gpar(fontsize=18))
  popViewport()
  pushViewport(viewport(x = 0.5, y = 0.655, 
                        width = unit(1, "cm"), height = unit(0.5, "cm"), 
                        just = c("left")))
    grid.text("true negative", gp=gpar(fontsize=18))
  popViewport()
  pushViewport(viewport(x = 0.5, y = 0.435, 
                        width = unit(1, "cm"), height = unit(0.5, "cm"), 
                        just = c("left")))
    grid.text("false positive", gp=gpar(fontsize=18))
  popViewport()
  pushViewport(viewport(x = 0.52, y = 0.215, 
                        width = unit(1, "cm"), height = unit(0.5, "cm"), 
                        just = c("left")))
    grid.text("false negative", gp=gpar(fontsize=18))
  popViewport()
  
  # plot legend squares:
  pushViewport(viewport(x = 0, y = 0.81, 
                        width = unit(1, "cm"), height = unit(0.5, "cm"), 
                        just = c("left")))
    grid.rect(width = unit(5, "mm"), height = unit(5, "mm"),
            just = c("left", "bottom"), gp=gpar(col = "#430F82", fill = "#430F82"))
  popViewport()
  pushViewport(viewport(x = 0, y = 0.595, 
                        width = unit(1, "cm"), height = unit(0.5, "cm"), 
                        just = c("left")))
    grid.rect(width = unit(5, "mm"), height = unit(5, "mm"),
              just = c("left", "bottom"), gp=gpar(col = "#B488B4", fill = "#B488B4"))
  popViewport()
  pushViewport(viewport(x = 0, y = 0.38, 
                        width = unit(1, "cm"), height = unit(0.5, "cm"), 
                        just = c("left")))
    grid.rect(width = unit(5, "mm"), height = unit(5, "mm"),
            just = c("left", "bottom"), gp=gpar(col = "#F6DC15", fill = "#F6DC15"))
  popViewport()
  pushViewport(viewport(x = 0, y = 0.165, 
                        width = unit(1, "cm"), height = unit(0.5, "cm"), 
                        just = c("left")))
    grid.rect(width = unit(5, "mm"), height = unit(5, "mm"),
            just = c("left", "bottom"), gp=gpar(col = "#7CBA61", fill = "#7CBA61"))
  popViewport()
}