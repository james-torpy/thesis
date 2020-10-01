create_legend <- function(annot_type, labs, sort_labs = TRUE, cols, lib_loc) {

  require(naturalsort, lib.loc = lib_loc)

  if (sort_labs) {
    l_labels <- gsub("_", " ", naturalsort(unique(labs)))
  } else {
    l_labels <- gsub("_", " ", unique(labs)) 
  }

  # add title:
  pushViewport(viewport(x = 0.06, y = 1, width = unit(2, "cm"), height = unit(0.5, "cm")))
    if (annot_type != "Normal vs cancer") {
      grid.text(
        paste0(annot_type, "\nclusters"), 
        gp=gpar(fontsize=20), just = "left"
      )
    } else {
      grid.text(
        paste0(annot_type, "\nclusters"),
        gp=gpar(fontsize=20), just = "left"
      )
    }
    #grid.rect()
  popViewport()
  for (l in 1:length(l_labels)) {
    # add labels:
    pushViewport(viewport(
      x = 0.25, 
      y = 0.95-(0.12*l), 
      width = unit(2, "cm"), 
      height = unit(0.5, "cm")
    ))
      #grid.rect()
      grid.text(l_labels[l], gp=gpar(fontsize=18), just = "left")
    popViewport()
    # add boxes:
    pushViewport(viewport(
      x = 0.14, 
      y = 0.95-(0.12*l), 
      width = unit(0.7, "cm"), 
      height = unit(0.7, "cm")
    ))
      grid.rect(gp=gpar(col = cols[l], fill = cols[l]))
    popViewport()
  }

}