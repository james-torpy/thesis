#! /share/ClusterShare/software/contrib/briglo/octR/src/R-3.6.0/builddir/bin/Rscript

###############################################################################
### Simulates cancer datasets taking as inputs number of CNVs, CNV lengths, ###
### gain/loss multipliers and downsamples required                          ###
###############################################################################

args = commandArgs(trailingOnly=TRUE)

project_name <- "thesis"
subproject_name <- "Figure_1.4_artefact_analysis"
sample_name <- args[1]
print(paste0("sample name = ", sample_name))
analysis_mode <- args[2]
print(paste0("Analysis mode of normal InferCNV run = ", 
  analysis_mode, "_mode"))

project_name <- "thesis"
subproject_name <- "Figure_1.4_artefact_analysis"
sample_name <- "CID4520N"
print(paste0("sample name = ", sample_name))
analysis_mode <- "samples"
print(paste0("Analysis mode of normal InferCNV run = ", 
  analysis_mode, "_mode"))
remove_cell_types <- c("No", "T_cells", "Endothelial", "Myeloid_cells", 
  "Endothelial.Myeloid_cells")
print("Plotting signal from InferCNV output without: ")
print(remove_cell_types)
remove_cell_type_labels <- c("No cells", "T cells", "Endothelial cells", "Myeloid cells", 
  "Endothelial and\nmyeloid cells")

print(paste0("Project name = ", project_name))
print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(cowplot)
library(ggplot2)
library(ComplexHeatmap, lib.loc = lib_loc)
library(circlize, lib.loc = lib_loc)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", 
  project_name, "/", subproject_name, "/")
func_dir <- paste0(project_dir, "scripts/functions/")
results_dir <- paste0(project_dir, "results/")

in_dir <- paste0(results_dir, "infercnv/", sample_name, "/")
out_dir <- paste0(in_dir, "final_plot/")
system(paste0("mkdir -p ", out_dir))

print(paste0("In directory = ", in_dir))
print(paste0("Out directory = ", out_dir))

print(paste0("Plotting selected artefact signal plots from ", 
  sample_name, "..."))


###############################################################################
### 1. Load plot data ###
###############################################################################

remove_cell_list <- as.list(remove_cell_types)
names(remove_cell_list) <- remove_cell_list

plot_data <- lapply(remove_cell_list, function(x) {
  return(
    readRDS(
      paste0(in_dir, x, "_removed/", analysis_mode, 
        "_mode/Rdata/signal_plot_data.Rdata")
    )
  )
})

# plot signal plots together:
pdf(
  paste0(out_dir, "selected_artefact_annotations.pdf"), 
  height = 10, 
  width = 22
)   
  grid.newpage()

    # plot one viewport per plot:
    # plot 1:
    pushViewport(viewport(x = 0.525, y = 0.85, width = 0.97, height = 0.18))
      #grid.rect()
      # plot label:
      pushViewport(viewport(x = 0.025, y = 0.5, width = 0.08, height = 0.4))
        #grid.rect()
        grid.text(paste0(remove_cell_type_labels[1], "\nremoved"), 
          gp=gpar(fontsize=26, fontface = "bold"))
      popViewport()
      # plot artefact annotation:
      pushViewport(viewport(x = 0.55, y = 0.5, width = 0.9, height = 0.9))
        grid.draw(plot_data[[1]]$artefact_heatmap_obj)
        #grid.rect()
      popViewport()
      # plot chromosome borders:
      pushViewport(viewport(x = 0.55, y = 0.5, width = 0.9, height = 0.9))
        grid.draw(plot_data[[1]]$chr_border_plot)
        #grid.rect()
      popViewport()
    popViewport()

    # plot 2:
    pushViewport(viewport(x = 0.525, y = 0.66, width = 0.97, height = 0.18))
      #grid.rect()
      # plot label:
      pushViewport(viewport(x = 0.025, y = 0.5, width = 0.08, height = 0.4))
        #grid.rect()
        grid.text(paste0(remove_cell_type_labels[2], "\nremoved"), 
          gp=gpar(fontsize=26, fontface = "bold"))
      popViewport()
      # plot artefact annotation:
      pushViewport(viewport(x = 0.55, y = 0.5, width = 0.9, height = 0.9))
        grid.draw(plot_data[[2]]$artefact_heatmap_obj)
        #grid.rect()
      popViewport()
      # plot chromsome borders:
      pushViewport(viewport(x = 0.55, y = 0.5, width = 0.9, height = 0.9))
        grid.draw(plot_data[[2]]$chr_border_plot)
        #grid.rect()
      popViewport()
    popViewport()

    # plot 3:
    pushViewport(viewport(x = 0.525, y = 0.47, width = 0.97, height = 0.18))
      #grid.rect()
      # plot label:
      pushViewport(viewport(x = 0.025, y = 0.5, width = 0.08, height = 0.4))
        #grid.rect()
        grid.text(paste0(remove_cell_type_labels[3], "\nremoved"), 
          gp=gpar(fontsize=26, fontface = "bold"))
      popViewport()
      # plot artefact annotation:
      pushViewport(viewport(x = 0.55, y = 0.5, width = 0.9, height = 0.9))
        grid.draw(plot_data[[3]]$artefact_heatmap_obj)
        #grid.rect()
      popViewport()
      # plot chromsome borders:
      pushViewport(viewport(x = 0.55, y = 0.5, width = 0.9, height = 0.9))
        grid.draw(plot_data[[3]]$chr_border_plot)
        #grid.rect()
      popViewport()
    popViewport()

    # plot 4:
    pushViewport(viewport(x = 0.525, y = 0.28, width = 0.97, height = 0.18))
      #grid.rect()
      # plot label:
      pushViewport(viewport(x = 0.025, y = 0.5, width = 0.08, height = 0.4))
        #grid.rect()
        grid.text(paste0(remove_cell_type_labels[4], "\nremoved"), 
          gp=gpar(fontsize=26, fontface = "bold"))
      popViewport()
      # plot artefact annotation:
      pushViewport(viewport(x = 0.55, y = 0.5, width = 0.9, height = 0.9))
        grid.draw(plot_data[[4]]$artefact_heatmap_obj)
        #grid.rect()
      popViewport()
      # plot chromsome borders:
      pushViewport(viewport(x = 0.55, y = 0.5, width = 0.9, height = 0.9))
        grid.draw(plot_data[[4]]$chr_border_plot)
        #grid.rect()
      popViewport()
    popViewport()

    # plot 5:
    pushViewport(viewport(x = 0.525, y = 0.09, width = 0.97, height = 0.18))
      #grid.rect()
      # plot label:
      pushViewport(viewport(x = 0.025, y = 0.5, width = 0.08, height = 0.4))
        #grid.rect()
        grid.text(paste0(remove_cell_type_labels[5], "\nremoved"), 
          gp=gpar(fontsize=26, fontface = "bold"))
      popViewport()
      # plot artefact annotation:
      pushViewport(viewport(x = 0.55, y = 0.5, width = 0.9, height = 0.9))
        grid.draw(plot_data[[5]]$artefact_heatmap_obj)
        #grid.rect()
      popViewport()
      # plot chromsome borders:
      pushViewport(viewport(x = 0.55, y = 0.5, width = 0.9, height = 0.9))
        grid.draw(plot_data[[5]]$chr_border_plot)
        #grid.rect()
      popViewport()
    popViewport()

    # plot chromosome annotation:
    pushViewport(viewport(x = 0.555, y = 0.95, width = 0.85, height = 0.07))
      #grid.rect()
      for ( e in 1:length(plot_data[[1]]$chr_data$lab_pos) ) {
        if (e==1) {
          pushViewport(viewport(x = plot_data[[1]]$chr_data$lab_pos[e], 
            y = 0.5, width = 0.01, height = 0.5))
            grid.text(names(plot_data[[1]]$chr_data$lab_pos)[e], 
              gp=gpar(fontsize=26))
          popViewport()
        } else if (e==21) {
          pushViewport(viewport(x = plot_data[[1]]$chr_data$lab_pos[e], 
            y = 0.9, width = 0.01, height = 0.5))
            grid.text(gsub("chr", "", names(plot_data[[1]]$chr_data$lab_pos)[e]), 
              gp=gpar(fontsize=26))
          popViewport()
        } else {
          pushViewport(viewport(x = plot_data[[1]]$chr_data$lab_pos[e], 
            y = 0.5, width = 0.01, height = 0.5))
            grid.text(gsub("chr", "", names(plot_data[[1]]$chr_data$lab_pos)[e]), 
              gp=gpar(fontsize=26))
          popViewport()
        }
      }
    popViewport()

    # plot boxes:
    # chr 2 artefact:
    pushViewport(viewport(x = 0.255, y = 0.47, width = 0.015, height = 0.93))
      grid.rect(gp = gpar(col = "#D95F02", fill = "transparent", lwd = 4, lty = "longdash"))
    popViewport()
    # chr 5 artefact:
    pushViewport(viewport(x = 0.385, y = 0.47, width = 0.015, height = 0.93))
      grid.rect(gp = gpar(col = "#D95F02", fill = "transparent", lwd = 4, lty = "longdash"))
    popViewport()
    # chr 6 artefact:
    pushViewport(viewport(x = 0.41, y = 0.47, width = 0.015, height = 0.93))
      grid.rect(gp = gpar(col = "#D95F02", fill = "transparent", lwd = 4, lty = "longdash"))
    popViewport()
    # chr 11 artefact:
    pushViewport(viewport(x = 0.591, y = 0.47, width = 0.015, height = 0.93))
      grid.rect(gp = gpar(col = "#D95F02", fill = "transparent", lwd = 4, lty = "longdash"))
    popViewport()
    # chr 12 artefact:
    pushViewport(viewport(x = 0.664, y = 0.47, width = 0.015, height = 0.93))
      grid.rect(gp = gpar(col = "#D95F02", fill = "transparent", lwd = 4, lty = "longdash"))
    popViewport()
    # chr 19 artefact:
    pushViewport(viewport(x = 0.875, y = 0.47, width = 0.015, height = 0.93))
      grid.rect(gp = gpar(col = "#D95F02", fill = "transparent", lwd = 4, lty = "longdash"))
    popViewport()

dev.off()

#convert pdf to png:
system(paste0("for p in ", out_dir, "*.pdf; do echo $p; f=$(basename $p); echo $f; ",
"new=$(echo $f | sed 's/.pdf/.png/'); echo $new; ", 
"convert -density 150 ", out_dir, "$f -quality 90 ", out_dir, "$new; done"))

