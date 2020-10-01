true_pos_false_pos_wrong_legend <- function() {
  # plot legend text:
  pushViewport(viewport(x = 0.48, y = 0.87, 
                        width = unit(1, "cm"), height = unit(0.5, "cm"), 
                        just = c("left")))
      grid.text("true positive", gp=gpar(fontsize=18))
  popViewport()
  pushViewport(viewport(x = 0.5, y = 0.655, 
                        width = unit(1, "cm"), height = unit(0.5, "cm"), 
                        just = c("left")))
    grid.text("false positive", gp=gpar(fontsize=18))
  popViewport()
  pushViewport(viewport(x = 0.43, y = 0.435, 
                        width = unit(1, "cm"), height = unit(0.5, "cm"), 
                        just = c("left")))
    grid.text("wrong call", gp=gpar(fontsize=18))
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
              just = c("left", "bottom"), gp=gpar(col = "#F6DC15", fill = "#F6DC15"))
  popViewport()
  pushViewport(viewport(x = 0, y = 0.38, 
                        width = unit(1, "cm"), height = unit(0.5, "cm"), 
                        just = c("left")))
    grid.rect(width = unit(5, "mm"), height = unit(5, "mm"),
            just = c("left", "bottom"), gp=gpar(col = "#C02456", fill = "#C02456"))
  popViewport()
}