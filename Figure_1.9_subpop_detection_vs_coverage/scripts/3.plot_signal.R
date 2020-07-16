#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

#### Generate CNV heatmaps with following annotations: ###
# CNV subclusters

project_name <- "thesis"
subproject_name <- "Figure_1.9_subpop_detection_vs_coverage"
args = commandArgs(trailingOnly=TRUE)

sample_name <- args[1]
simulation_number <- args[2]
downsample_proportion <- args[3]
remove_artefacts <- args[4]
QC_annot <- as.logical(args[5])
plot_references <- as.logical(args[6])

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("Simulation number ", simulation_number))
print(paste0("Downsample proportion ", downsample_proportion))
print(paste0("Remove artefacts? ", remove_artefacts))
print(paste0("Print QC annotations? ", QC_annot))
print(paste0("Plot reference cells? ", plot_references))

sample_name <- "CID4520N_cancer_sim"
simulation_number <- "3"
downsample_proportion <- "no"
remove_artefacts <- "artefacts_removed"
plot_references <- FALSE

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(scales, lib.loc = lib_loc)
library(cluster)
library(ggplot2)
library(Seurat)
library(infercnv, lib.loc = lib_loc)
library(RColorBrewer)
library(ComplexHeatmap, lib.loc=lib_loc)
library(circlize, lib.loc = lib_loc)
library(reshape2)
library(fpc, lib.loc = lib_loc)
library(dplyr)
library(naturalsort, lib.loc = lib_loc)
library(cowplot)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/",
  subproject_name, "/")
ref_dir <- paste0(project_dir, "/refs/")
func_dir <- paste0(project_dir, "/scripts/functions/")
results_dir <- paste0(project_dir, "results/")
raw_dir <- paste0(project_dir, "raw_files/")
seurat_dir <- paste0(raw_dir, "seurat_objects/", sample_name, "/")
in_dir <- paste0(results_dir, "infercnv/", sample_name, "/", 
  simulation_number, "/", downsample_proportion, "_downsampling/")
input_dir <- paste0(in_dir, "/input_files/")
sim_Robject_dir <- paste0(results_dir, "cancer_simulation/",
  sample_name, "/", simulation_number, "/Rdata/")

out_dir <- paste0(in_dir, remove_artefacts, "/")

Robject_dir <- paste0(out_dir, "Rdata/")
system(paste0("mkdir -p ", Robject_dir))
plot_dir <- paste0(out_dir, "plots/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(out_dir, "tables/")
system(paste0("mkdir -p ", table_dir))

print(paste0("In directory = ", in_dir))
print(paste0("Reference directory = ", ref_dir))
print(paste0("R function directory = ", func_dir))
print(paste0("R object directory = ", Robject_dir))
print(paste0("Table directory = ", table_dir))
print(paste0("Plot directory = ", plot_dir))


################################################################################
### 0. Define functions ###
################################################################################

fetch_chromosome_boundaries <- dget(paste0(func_dir, 
  "fetch_chromosome_boundaries.R"))
get_subpops <- dget(paste0(func_dir, "get_subpops.R"))
ts_heatmap <- dget(paste0(func_dir, "ts_heatmap.R"))

extra_colours <- c(brewer.pal(12, "Set3"), brewer.pal(12, "Paired"))
col_palette <- c(brewer.pal(8, "Dark2")[c(3:6,8)], brewer.pal(12, "Set3"), brewer.pal(8, "Accent"),
  "#660200", "#918940", "#9ECAE1", "#3F007D", "#67000D", "#FDD0A2", "#08519C",
  "#DBB335", "#5AA050", "#807DBA", "#1A6000", "#F16913", "#FD8D3C", "#DEEBF7", 
  "#7F2704", "#DADAEB", "#FC9272", "#BCBDDC", extra_colours, "#b2182b", "#85929E", 
  "#9B59B6", "#74add1","#1b7837", "#b8e186", "#fed976","#e7298a", "#18ffff", "#ef6c00",
  "#A93226", "#E611ED","orange", "#b8bc53", "#5628ce", "#fa909c", "#8ff331")
col_palette <- col_palette[-7]


################################################################################
### 1. Load InferCNV output and create heatmap and metadata dfs ###
################################################################################

if (!file.exists(paste0(
  Robject_dir, "/1b.initial_epithelial_metadata.Rdata")
)) {

  # load InferCNV output:
  print("Loading InferCNV output files...")
  epithelial_heatmap <- as.data.frame(t(read.table(paste0(in_dir, 
  	"infercnv.12_denoised.observations.txt"))))

  ######
  # subset:
  # epithelial_heatmap <- epithelial_heatmap[1:200,1:500]
  ######

  # load and isolate epithelial metadata:
  initial_metadata <- read.table(
    paste0(input_dir, "metadata.txt"),
    sep = "\t",
    header = F,
  )
  colnames(initial_metadata) <- c("cell_ids", "cell_types")
  epithelial_metadata <- initial_metadata[
    grep("Epithelial", initial_metadata$cell_type),
  ]
  # replace dashes with full stops for cell ids if needed:
  epithelial_metadata$cell_ids <- gsub(
    "-", ".", epithelial_metadata$cell_ids
  )

  # remove cells not present in heatmap:
  print(paste0(
    "No cells in metadata before filtering for those in heatmap only = ", 
    nrow(epithelial_metadata)
  ))
  epithelial_metadata <- epithelial_metadata %>%
    filter(cell_ids %in% rownames(epithelial_heatmap))
  print(paste0(
    "No cells in metadata after filtering for those in heatmap only = ", 
    nrow(epithelial_metadata)
  ))

  saveRDS(epithelial_heatmap, paste0(Robject_dir, 
    "/1a.initial_epithelial_heatmap.Rdata"))
  saveRDS(epithelial_metadata, paste0(Robject_dir, 
    "/1b.initial_epithelial_metadata.Rdata"))

} else {

  print("Loading heatmap and metadata dfs...")
  epithelial_heatmap <- readRDS(paste0(Robject_dir, 
    "/1a.initial_epithelial_heatmap.Rdata"))
  epithelial_metadata <- readRDS(paste0(Robject_dir, 
    "/1b.initial_epithelial_metadata.Rdata"))
 
}


################################################################################
### 2. Remove artefacts ###
################################################################################

if (remove_artefacts == "artefacts_removed") {

  if (!file.exists(paste0(Robject_dir, 
    "/2.epithelial_heatmap_artefacts_removed.Rdata"))) {

    # load known artefact genes:
    artefact_genes <- readRDS(paste0(ref_dir, "artefact_gene_record.Rdata"))
    
    # determine diploid signal
    signal_table <- table(round(unlist(epithelial_heatmap), 6))
    diploid_signal <- as.numeric(names(signal_table)[which.max(signal_table)])
    
    # locate artefact within heatmap and make all signal within diploid:
    for (a in 1:length(artefact_genes)) {
    	print(a)
      # identify first and last gene in artefact present in epithelial heatmap:
      first_genes <- as.character(
        artefact_genes[[a]]$first_ten[
          artefact_genes[[a]]$first_ten %in% colnames(epithelial_heatmap)
        ]
      )
      first_gene <- first_genes[1]
    
      last_genes <- as.character(
        artefact_genes[[a]]$last_ten[
          artefact_genes[[a]]$last_ten %in% colnames(epithelial_heatmap)
        ]
      )
      last_gene <- last_genes[length(last_genes)]
    
      heatmap_coords <- c(
        grep(paste0("^", first_gene, "$"), colnames(epithelial_heatmap)),
        grep(paste0("^", last_gene, "$"), colnames(epithelial_heatmap))
      )
  
      # add gene tails equal to 1/3 of the artefact length to start and end of 
      # artefact coordinates to account for residue signal:
      tail_length <- round((heatmap_coords[2]-heatmap_coords[1])/3, 0)
      heatmap_coords[1] <- heatmap_coords[1] - tail_length
      heatmap_coords[2] <- heatmap_coords[2] + tail_length
      coord_vector <- heatmap_coords[1]:heatmap_coords[2]
  
      epithelial_heatmap[,coord_vector] <- diploid_signal
    
    }

    saveRDS(epithelial_heatmap, paste0(Robject_dir, 
      "/2.epithelial_heatmap_artefacts_removed.Rdata"))
    
  } else {

  	epithelial_heatmap <- readRDS(paste0(Robject_dir, 
      "/2.epithelial_heatmap_artefacts_removed.Rdata"))

  }
}


################################################################################
### 3. Add subcluster annotation and order by metadata ###
################################################################################

# load final output object and fetch subcluster data:  
final_infercnv_obj <- readRDS(paste0(in_dir, "run.final.infercnv_obj"))
epithelial_metadata <- get_subpops(final_infercnv_obj, epithelial_metadata)

# add and order by known subpops:
known_subpops <- readRDS(paste0(sim_Robject_dir, "subpop_ids.Rdata"))

#####

test_cells <- known_subpops[[1]]
test_metadata <- epithelial_metadata[
  epithelial_metadata$cell_ids %in% test_cells,
]

#####


known_df <- data.frame(
  cell_ids = do.call("c", known_subpops),
  subpop_ids = c(
    rep(names(known_subpops)[1], length(known_subpops[[1]])),
    rep(names(known_subpops)[2], length(known_subpops[[2]]))
  )
)
epithelial_metadata <- merge(epithelial_metadata, known_df, by="cell_ids")
epithelial_metadata <- epithelial_metadata[
  naturalorder(epithelial_metadata$subcluster_ids),
]
epithelial_metadata <- epithelial_metadata[
  naturalorder(epithelial_metadata$subpop_ids),
]

#######
#epithelial_metadata <- epithelial_metadata[
#  epithelial_metadata$subcluster_ids == "CNV_5",
#]
#######

# adjust order of heatmap:
epithelial_heatmap <- epithelial_heatmap[epithelial_metadata$cell_ids,]

print(paste0(
  "Is heatmap order the same as metadata? ", 
  identical(rownames(epithelial_heatmap), epithelial_metadata$cell_ids)
))

#######
#ts_heatmap(epithelial_heatmap, epithelial_metadata, plot_dir)
#######


################################################################################
### 4. Create heatmap and annotations ###
################################################################################

# create known subpop annotation:
# define cluster annotation colours:
subpop_number <- length(unique(epithelial_metadata$subpop_ids))
subpop_cols <- col_palette[1:subpop_number]
names(subpop_cols) <- unique(epithelial_metadata$subpop_ids)
subpop_annot_df <- subset(epithelial_metadata, select = subpop_ids)
subpop_annot <- Heatmap(
  as.matrix(subpop_annot_df), 
  col = subpop_cols,
  name = "Simulated\nsubpopulation", 
  width = unit(4, "mm"), 
  show_row_names = F, show_column_names = F,
  show_heatmap_legend = FALSE
)

# create CNV subcluster annotation:
subcluster_annot_df <- subset(epithelial_metadata, select = subcluster_ids)
subcluster_cols <- rev(col_palette)[
  1:length(unique(subcluster_annot_df$subcluster_ids))
]
names(subcluster_cols) <- unique(subcluster_annot_df$subcluster_ids)
CNV_subcluster_annot <- Heatmap(as.matrix(subcluster_annot_df), 
  col = subcluster_cols, 
  name = "CNV\nsubclusters", width = unit(6, "mm"), 
  show_row_names = F, show_column_names = F, 
  show_heatmap_legend = F,
  heatmap_legend_param = list(
    title = "",
    labels_gp = gpar(fontsize = 16)
  )
)


################################################################################
### 5. Create epithelial and reference heatmap ###
################################################################################

# determine heatmap colours
if (plot_references) {

  # load InferCNV output:
  print("Loading InferCNV reference output files...")
  ref_heatmap <- as.data.frame(t(read.table(paste0(in_dir, 
    "infercnv.12_denoised.references.txt"))))

  both_heatmap <- rbind(epithelial_heatmap, ref_heatmap)

  # prepare df for plotting:
  col_heatmap <- both_heatmap
  # define heatmap colours:
  na_less_vector <- unlist(col_heatmap)
  na_less_vector <- na_less_vector[!is.na(na_less_vector)]
  heatmap_cols <- colorRamp2(c(min(na_less_vector), 1, max(na_less_vector)), 
        c("#00106B", "white", "#680700"), space = "sRGB")

  # prepare df for plotting:
  ref_object <- ref_heatmap
  colnames(ref_object) <- rep("la", ncol(ref_object))
  
  print("Generating reference heatmap...")
  # create main CNV heatmap:
  final_ref_heatmap <- Heatmap(
    as.matrix(ref_object), name = paste0("ref_hm"), 
    col = heatmap_cols,
    cluster_columns = F, cluster_rows = F,
    show_row_names = F, show_column_names = T,
    column_names_gp = gpar(col = "white"),
    show_row_dend = F,
    show_heatmap_legend = F,
    use_raster = T, raster_device = c("png")
  )

  # converrt to grid object:
  grid_ref_heatmap <- grid.grabExpr(
    draw(final_ref_heatmap, gap = unit(6, "mm"))
  )
  dev.off()

} else {
 
  # define heatmap colours:
  na_less_vector <- unlist(epithelial_heatmap)
  na_less_vector <- na_less_vector[!is.na(na_less_vector)]
  heatmap_cols <- colorRamp2(c(min(na_less_vector), 1, max(na_less_vector)), 
        c("#00106B", "white", "#680700"), space = "sRGB")

}

# prepare df for plotting:
plot_object <- epithelial_heatmap
colnames(plot_object) <- rep("la", ncol(plot_object))

print("Generating final heatmap...")
# create main CNV heatmap:
final_heatmap <- Heatmap(
  as.matrix(plot_object), name = paste0("hm"), 
  col = heatmap_cols,
  cluster_columns = F, cluster_rows = F,
  show_row_names = F, show_column_names = T,
  column_names_gp = gpar(col = "white"),
  show_row_dend = F,
  show_heatmap_legend = F,
  use_raster = T, raster_device = c("png")
)

chr_data <- fetch_chromosome_boundaries(epithelial_heatmap, ref_dir)
 
# generate heatmap legend:
signal_ranges <- round(range(unlist(plot_object)), 2)
lgd <- Legend(
  at = c(signal_ranges[1], 1, signal_ranges[2]),
  col_fun = heatmap_cols, 
  title = "CNV signal", 
  direction = "horizontal",
  grid_height = unit(2.5, "cm"),
  grid_width = unit(0.1, "cm"),
  labels_gp = gpar(fontsize = 16),
  title_gp = gpar(fontsize = 22, fontface = "plain")
)

# organise chromosome labels:
chr_labels <- chr_data$lab_pos
names(chr_labels) <- gsub("chr", "", names(chr_labels))
names(chr_labels) <- gsub("^1$", "chr1", names(chr_labels))
names(chr_labels) <- gsub("21", "\n21", names(chr_labels))


################################################################################
### 6. Create and plot annotated heatmap ###
################################################################################

ht_list <- subpop_annot + CNV_subcluster_annot + final_heatmap

annotated_heatmap <- grid.grabExpr(
  draw(ht_list, gap = unit(6, "mm"), heatmap_legend_side = "left")
)
dev.off()

#save.image(paste0(Robject_dir, "temp.Rdata"))
#load(paste0(Robject_dir, "temp.Rdata"))

if (plot_references) {

  # plot final annotated heatmap:
  png(paste0(plot_dir, "infercnv_plot.png"), 
    height = 13, width = 20, res = 300, units = "in") 
  
    grid.newpage()

      # create reference heatmap viewport:
      pushViewport(viewport(x = 0.56, y = 0.68, width = 0.765, height = 0.3,
        just = "bottom"))
        #grid.rect()
        # plot heatmap:
        grid.draw(grid_ref_heatmap)
        decorate_heatmap_body("ref_hm", {
          for ( e in 1:length(chr_data$end_pos) ) {
            grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), 
              gp = gpar(lwd = 1, col = "#383838"))
          }
        })
      popViewport()

      # create reference heatmap label:
      pushViewport(viewport(x = 0.1, y = 0.95, 
        width = unit(2, "cm"), height = unit(0.5, "cm")))
        #grid.rect()
        grid.text("Reference cells:", gp=gpar(fontsize=20))
      popViewport()

      # create epithelial heatmap viewport:
      pushViewport(viewport(x = 0.5, y = 0, width = 1, height = 0.75,
        just = "bottom"))

         # plot heatmap:
        pushViewport(viewport(x = 0.155, y = 0.065, width = 0.85, height = 0.85, 
          just = c("left", "bottom")))
          #grid.rect()
          grid.draw(annotated_heatmap)
          decorate_heatmap_body("hm", {
            for ( e in 1:length(chr_data$end_pos) ) {
            grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), 
              gp = gpar(lwd = 1, col = "#383838"))
            grid.text(names(chr_labels)[e], chr_labels[e], 
                unit(0, "npc") + unit(-3.5, "mm"), gp=gpar(fontsize=18))
            }
          })
        popViewport()
   
        # plot heatmap legend:
        pushViewport(viewport(x = unit(2, "cm"), y = unit(19, "cm"), width = unit(0.1, "cm"), 
          height = unit(0.4, "cm"), just = c("right", "bottom")))
          draw(lgd, x = unit(0.1, "cm"), y = unit(0.1, "cm"), just = c("left", "bottom"))
        popViewport()
     
        # label annotations:
        pushViewport(viewport(x=0.95, y=0.1, width = 0.1, height = 0.1, 
          just = "top"))
          grid.text("nUMI", rot=65, gp=gpar(fontsize=18))
        popViewport()
        pushViewport(viewport(x=0.98, y=0.1, width = 0.1, height = 0.1, 
          just = "top"))
          grid.text("nGene", rot=65, gp=gpar(fontsize=18))
        popViewport()

      popViewport()
        
  dev.off()

  pdf(paste0(plot_dir, "infercnv_plot.pdf"), 
    height = 13, width = 20) 
  
    grid.newpage()

    # create reference heatmap viewport:
    pushViewport(viewport(x = 0.56, y = 0.68, width = 0.765, height = 0.3,
      just = "bottom"))
      #grid.rect()
      # plot heatmap:
      grid.draw(grid_ref_heatmap)
      decorate_heatmap_body("ref_hm", {
        for ( e in 1:length(chr_data$end_pos) ) {
          grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), 
            gp = gpar(lwd = 1, col = "#383838"))
        }
      })
    popViewport()
    # create reference heatmap label:
    pushViewport(viewport(x = 0.1, y = 0.95, 
      width = unit(2, "cm"), height = unit(0.5, "cm")))
      #grid.rect()
      grid.text("Reference cells:", gp=gpar(fontsize=20))
    popViewport()
    # create epithelial heatmap viewport:
    pushViewport(viewport(x = 0.5, y = 0, width = 1, height = 0.75,
      just = "bottom"))
       # plot heatmap:
      pushViewport(viewport(x = 0.155, y = 0.065, width = 0.85, height = 0.85, 
        just = c("left", "bottom")))
        #grid.rect()
        grid.draw(annotated_heatmap)
        decorate_heatmap_body("hm", {
          for ( e in 1:length(chr_data$end_pos) ) {
          grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), 
            gp = gpar(lwd = 1, col = "#383838"))
          grid.text(names(chr_labels)[e], chr_labels[e], 
              unit(0, "npc") + unit(-3.5, "mm"), gp=gpar(fontsize=18))
          }
        })
      popViewport()
 
      # plot heatmap legend:
      pushViewport(viewport(x = unit(2, "cm"), y = unit(19, "cm"), width = unit(0.1, "cm"), 
        height = unit(0.4, "cm"), just = c("right", "bottom")))
        draw(lgd, x = unit(0.1, "cm"), y = unit(0.1, "cm"), just = c("left", "bottom"))
      popViewport()
  
      # label annotations:
      pushViewport(viewport(x=0.95, y=0.1, width = 0.1, height = 0.1, 
        just = "top"))
        grid.text("nUMI", rot=65, gp=gpar(fontsize=18))
      popViewport()
      pushViewport(viewport(x=0.98, y=0.1, width = 0.1, height = 0.1, 
        just = "top"))
        grid.text("nGene", rot=65, gp=gpar(fontsize=18))
      popViewport()
    popViewport()
        
  dev.off()

} else {

  # plot final annotated heatmap:
  png(paste0(plot_dir, "infercnv_plot.png"), 
    height = 13, width = 20, res = 300, units = "in") 
  
    grid.newpage()
  
      # plot heatmap:
      pushViewport(viewport(x = 0.155, y = 0.065, width = 0.8, height = 0.85, 
        just = c("left", "bottom")))
        #grid.rect()
        grid.draw(annotated_heatmap)
        decorate_heatmap_body("hm", {
          for ( e in 1:length(chr_data$end_pos) ) {
          grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), 
            gp = gpar(lwd = 1, col = "#383838"))
          grid.text(names(chr_labels)[e], chr_labels[e], 
              unit(0, "npc") + unit(-3.5, "mm"), gp=gpar(fontsize=18))
          }
        })
      popViewport()

      # plot subpop legend:
      l_labels <- gsub("_", " ", naturalsort(unique(epithelial_metadata$subpop_ids)))
      pushViewport(viewport(x = unit(4.5, "cm"), y = unit(8, "cm"), width = unit(6, "cm"), 
        height = unit(10, "cm")))
        #grid.rect()
        # add title:
        pushViewport(viewport(x = 0.06, y = 1, width = unit(2, "cm"), height = unit(0.5, "cm")))
          grid.text("Simulated\nsubpopulation", gp=gpar(fontsize=20), just = "left")
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
            grid.rect(gp=gpar(col = subpop_cols[l], fill = subpop_cols[l]))
          popViewport()
        }
      popViewport()

      # plot CNV legend:
      l_labels <- gsub("_", " ", unique(epithelial_metadata$subcluster_id))

      pushViewport(viewport(x = unit(4.4, "cm"), y = unit(26, "cm"), width = unit(6, "cm"), 
        height = unit(10, "cm")))
        #grid.rect()

        # add title:
        pushViewport(viewport(x = 0.06, y = 1, width = unit(2, "cm"), height = unit(0.5, "cm")))
          grid.text("CNV\nclusters", gp=gpar(fontsize=20), just = "left")
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
            grid.rect(gp=gpar(col = subcluster_cols[l], fill = subcluster_cols[l]))
          popViewport()
        }
      popViewport()

      # plot heatmap legend:
      pushViewport(viewport(x = unit(2, "cm"), y = unit(15.5, "cm"), width = unit(0.1, "cm"), 
        height = unit(0.4, "cm"), just = c("right", "bottom")))
        draw(lgd, x = unit(0.1, "cm"), y = unit(0.1, "cm"), just = c("left", "bottom"))
      popViewport()

  dev.off()
  
  # plot final annotated heatmap:
  pdf(paste0(plot_dir, "infercnv_plot.pdf"), 
    height = 13, width = 20) 
  
    grid.newpage()
  
    # plot heatmap:
    pushViewport(viewport(x = 0.155, y = 0.065, width = 0.8, height = 0.85, 
      just = c("left", "bottom")))
      #grid.rect()
      grid.draw(annotated_heatmap)
      decorate_heatmap_body("hm", {
        for ( e in 1:length(chr_data$end_pos) ) {
        grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), 
          gp = gpar(lwd = 1, col = "#383838"))
        grid.text(names(chr_labels)[e], chr_labels[e], 
            unit(0, "npc") + unit(-3.5, "mm"), gp=gpar(fontsize=18))
        }
      })
    popViewport()

    # plot subpop legend:
    l_labels <- gsub("_", " ", naturalsort(unique(epithelial_metadata$subpop_ids)))
    pushViewport(viewport(x = unit(4.5, "cm"), y = unit(8, "cm"), width = unit(6, "cm"), 
      height = unit(10, "cm")))
      #grid.rect()
      # add title:
      pushViewport(viewport(x = 0.06, y = 1, width = unit(2, "cm"), height = unit(0.5, "cm")))
        grid.text("Simulated\nsubpopulation", gp=gpar(fontsize=20), just = "left")
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
          grid.rect(gp=gpar(col = subpop_cols[l], fill = subpop_cols[l]))
        popViewport()
      }
    popViewport()

    # plot CNV legend:
    l_labels <- gsub("_", " ", unique(epithelial_metadata$subcluster_id))

    pushViewport(viewport(x = unit(4.4, "cm"), y = unit(26, "cm"), width = unit(6, "cm"), 
      height = unit(10, "cm")))
      #grid.rect()
      # add title:
      pushViewport(viewport(x = 0.06, y = 1, width = unit(2, "cm"), height = unit(0.5, "cm")))
        grid.text("CNV\nclusters", gp=gpar(fontsize=20), just = "left")
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
          grid.rect(gp=gpar(col = subcluster_cols[l], fill = subcluster_cols[l]))
        popViewport()
      }
    popViewport()

    # plot heatmap legend:
    pushViewport(viewport(x = unit(2, "cm"), y = unit(15.5, "cm"), width = unit(0.1, "cm"), 
      height = unit(0.4, "cm"), just = c("right", "bottom")))
      draw(lgd, x = unit(0.1, "cm"), y = unit(0.1, "cm"), just = c("left", "bottom"))
    popViewport()

  dev.off()

}

print(paste0("Heatmap created, output in ", plot_dir))
  
  