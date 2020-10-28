#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

project_name <- "thesis"
subproject_name <- "chapter_4"
args = commandArgs(trailingOnly=TRUE)
patient_name <- args[1]
subcluster_method <- args[2]
subcluster_p <- args[3]
if (subcluster_p != "none") {
  subcluster_p <- as.numeric(subcluster_p)
}
coverage_filter <- args[4]
remove_artefacts <- args[5]
subset_data <- as.logical(args[6])
QC_annot <- as.logical(args[7])

project_name <- "thesis"
subproject_name <- "chapter_4"
args = commandArgs(trailingOnly=TRUE)
patient_names <- c("CID4515", "CID4517")
subcluster_method <- "random_trees"
subcluster_p <- "0.05"
if (subcluster_p != "none") {
  subcluster_p <- as.numeric(subcluster_p)
}
coverage_filter <- "filtered"
remove_artefacts <- "artefacts_not_removed"
subset_data <- FALSE
na_colour <- "white"
QC_annot <- TRUE

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("Subclustered by ", subcluster_method))
print(paste0("Subclustered by pval ", subcluster_p))
print(paste0("Filter by coverage? ", coverage_filter))
print(paste0("Remove artefacts? ", remove_artefacts))

# determine order of samples:

sample_names <- paste0(patient_name[1], 1:2)

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample names = ", sample_names))
print(paste0("Subset data? ", as.character(subset_data)))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(naturalsort, lib.loc = lib_loc)
library(cluster, lib.loc = lib_loc)
library(ComplexHeatmap, lib.loc=lib_loc)
library(circlize, lib.loc = lib_loc)
library(scales, lib.loc = lib_loc)
library(RColorBrewer)
library(reshape2)
library(fpc, lib.loc = lib_loc)
library(dplyr)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/",
  subproject_name, "/")
ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/functions/")
results_dir <- paste0(project_dir, "results/")
in_path <- paste0(results_dir, "infercnv/")
out_path <- paste0(results_dir, "infercnv/matched_mets/", coverage_filter, 
  "/", subcluster_method, "/p_", subcluster_p, "/", remove_artefacts, "/")

if (subset_data) {
  Robject_dir <- paste0(out_path, "Rdata_sub/")
  table_dir <- paste0(out_path, "tables_sub/")
} else {
  Robject_dir <- paste0(out_path, "Rdata/")
  table_dir <- paste0(out_path, "tables/")
}
system(paste0("mkdir -p ", Robject_dir))
system(paste0("mkdir -p ", table_dir))

print(paste0("Reference directory = ", ref_dir))
print(paste0("R function directory = ", func_dir))
print(paste0("Common R object directory = ", Robject_dir))

print("Plotting primary and met CNV heatmaps of: ")
print(names(sample_names))


################################################################################
### 0. Define functions and colours ###
################################################################################

load_and_merge_heatmaps <- dget(
  paste0(func_dir, "load_and_merge_heatmaps.R")
)

create_matched_heatmap_annotations <- dget(
  paste0(func_dir, "create_matched_heatmap_annotations.R")
)

prepare_annotated_matched_heatmaps <- dget(
  paste0(func_dir, "prepare_annotated_matched_heatmaps.R")
)

subcluster_cols <- read.table(
  paste0(ref_dir, "CNV_colour_palette.txt"),
  header = F,
  comment.char = "",
  stringsAsFactors = F
)[,1]
colnames(subtype_cols) <- c("subcluster", "col")

all_cols <- read.table(
  paste0(ref_dir, "all_colour_palette.txt"),
  header = F,
  comment.char = "",
  stringsAsFactors = F
)[,1]

type_cols <- c(
  Primary = "#B066B2",
  Metastasis = "#F4D30B"
)


################################################################################
### 1. Load heatmap data and combine matched samples ###
################################################################################

if (!file.exists(
  paste0(Robject_dir, "matched_heatmaps_and_metadata.Rdata")
)) {

  all_data <- lapply(
    sample_names, 
    load_and_merge_heatmaps,
    in_path,
    coverage_filter,
    subcluster_method,
    subcluster_p,
    remove_artefacts,
    ref_dir, 
    subset_data = subset_data
  )
  
  saveRDS(
    all_data, 
    paste0(Robject_dir, "matched_heatmaps_and_metadata.Rdata")
  )

} else {

  all_data <- readRDS(
    paste0(Robject_dir, "matched_heatmaps_and_metadata.Rdata")
  )

}


################################################################################
### 2. Determine correlations between each primary subpop and main matched met 
# subpopulation ###
################################################################################

all_cors <- lapply(all_data, function(x) {

  met_meta <- x$metadata[x$metadata$type == "Metastasis",]
  split_met <- split(met_meta, met_meta$subcluster_id)
  main_met_index <- which.max(
    lapply(
      split_met,
      nrow
    )
  )
  main_met <- split_met[[main_met_index]]
  main_met_df <- x$heatmap[main_met$cell_ids, ]
  avg_main_met <- apply(main_met_df, 2, mean)
  avg_main_met <- avg_main_met[!is.na(avg_main_met)]

  prim_meta <- x$metadata[x$metadata$type == "Primary",]
  split_prim <- split(prim_meta, prim_meta$subcluster_id)

  cors <- lapply(split_prim, function(y) {

    prim_df <- x$heatmap[y$cell_ids, ]
    avg_prim <- apply(prim_df, 2, mean)
    avg_prim <- avg_prim[!is.na(avg_prim)]

    avg_prim <- avg_prim[names(avg_prim) %in% names(avg_main_met)]
    avg_main_met <- avg_main_met[names(avg_main_met) %in% names(avg_prim)]

    cor <- cor.test(
      as.numeric(avg_prim), 
      as.numeric(avg_main_met), 
      method = "pearson"
    )
    return(data.frame(cor$estimate, cor$p.value))
  })
  return(do.call("rbind", cors))
})

for (i in 1:length(all_cors)) {
  write.table(
    all_cors[[i]],
    paste0(
      table_dir, 
      names(all_cors)[i], 
      "_correlations_primaries_vs_main_met_subpop.txt"
    ),
    row.names = TRUE,
    col.names = TRUE,
    quote = FALSE
  )
}


################################################################################
### 3. Create heatmap annotations ###
################################################################################

if (!file.exists(
  paste0(Robject_dir, "matched_annotations.Rdata")
)) {

  all_annotations <- lapply(
    all_data,
    create_matched_heatmap_annotations,
    type_cols,
    subcluster_cols,
    func_dir,
    ref_dir
  )
  
  saveRDS(
    all_annotations, 
    paste0(Robject_dir, "matched_annotations.Rdata")
  )

} else {

  all_annotations <- readRDS(
    paste0(Robject_dir, "matched_annotations.Rdata")
  )

}


#############################################################################
### 3. Prepare heatmap and annotations ###
################################################################################

if (!file.exists(
  paste0(Robject_dir, "final_heatmap_data.Rdata")
)) {

  for (i in 1:length(all_data)) {
    if (i==1) {
      final_heatmap_data <- list(
        prepare_annotated_matched_heatmaps(
          all_data[[i]],
          all_annotations[[i]]
        )
      )
    } else {
      final_heatmap_data[[i]] <- prepare_annotated_matched_heatmaps(
        all_data[[i]],
        all_annotations[[i]]
      )
    }
  }
  
  saveRDS(
    final_heatmap_data, 
    paste0(Robject_dir, "final_heatmap_data.Rdata")
  )

} else {

  final_heatmap_data <- readRDS(
    paste0(Robject_dir, "final_heatmap_data.Rdata")
  )

}


###############################################################################
### 4. Plot heatmap and annotations ###
################################################################################

for (i in 1:length(sample_names)) {

  chr_labels <- all_annotations[[i]]$chr_data$lab_pos
  names(chr_labels) <- gsub("chr", "", names(chr_labels))
  names(chr_labels) <- gsub("^1$", "chr1", names(chr_labels))
  names(chr_labels) <- gsub("21", "\n21", names(chr_labels))

  names(all_annotations[[i]]$subcluster_cols) <- gsub(
    "_", " ",
    names(all_annotations[[i]]$subcluster_cols)
  )

  names(all_annotations[[i]]$subcluster_cols) <- gsub(
    "CNV", "CNA",
    names(all_annotations[[i]]$subcluster_cols)
  )

  if (subset_data) {
    out_dir <- paste0(out_path, names(sample_names)[i], "_sub/")
  } else {
    out_dir <- paste0(out_path, names(sample_names)[i], "/")
  }
  system(paste0("mkdir -p ", out_dir))

  pdf(
    paste0(
      out_dir, 
      names(sample_names)[i], 
      "_combined_and_rescaled.pdf"
    ), 
    height = 13, width = 24
  )

    grid.newpage()

      # plot heatmap:
      pushViewport(viewport(x = 0.155, y = 0.04, width = 0.82, height = 0.95, 
        just = c("left", "bottom")))
        #grid.rect()
        grid.draw(final_heatmap_data[[i]]$heatmap)
        decorate_heatmap_body("hm", {
          # draw chromosome lines and labels:
          for ( e in 1:length(all_annotations[[i]]$chr_data$end_pos) ) {
            grid.lines(
              x = c(
                all_annotations[[i]]$chr_data$end_pos[e], 
                all_annotations[[i]]$chr_data$end_pos[e]
              ), 
              y = c(0, 1), 
              gp = gpar(lwd = 2.5, col = "#383838")
            )
            grid.text(
              names(chr_labels)[e], chr_labels[e],
              unit(0, "npc") + unit(-3.1, "mm"), 
              gp=gpar(fontsize=18)
            )
          }
          # draw horizontal lines:
          for ( h in 1:length(final_heatmap_data[[i]]$hlines)) {
            grid.lines(
              x = c(0, 1), 
              y = c(final_heatmap_data[[i]]$hlines[h], final_heatmap_data[[i]]$hlines[h]), 
              gp = gpar(lwd = 2., col = "#383838")
            )
          }
        })
      popViewport()
  
      # plot heatmap legend:
      pushViewport(viewport(x = unit(3.5, "cm"), y = unit(28.3, "cm"), width = unit(4.5, "cm"), 
        height = unit(4, "cm")))
        draw(final_heatmap_data[[i]]$legend, x = unit(0.1, "cm"), y = unit(0.1, "cm"), just = c("left", "bottom"))
      popViewport()
  
      # plot type legend:
      pushViewport(viewport(x = unit(4, "cm"), y = unit(20, "cm"), width = unit(6, "cm"), 
        height = unit(10, "cm")))
        # add title:
        pushViewport(viewport(x = 0.06, y = 1, width = unit(2, "cm"), height = unit(0.5, "cm")))
          grid.text("Type", gp=gpar(fontsize=20), just = "left")
        popViewport()
        for (l in 1:length(type_cols)) {
          # add labels:
          pushViewport(viewport(
            x = 0.21, 
            y = 0.99-(0.08*l), 
            width = unit(2, "cm"), 
            height = unit(0.5, "cm")
          ))
            grid.text(names(type_cols)[l], gp=gpar(fontsize=18), just = "left")
          popViewport()
          # add dots:
          pushViewport(viewport(
            x = 0.1, 
            y = 0.99-(0.08*l), 
            width = unit(0.3, "cm"), 
            height = unit(0.3, "cm")
          ))
            grid.circle(gp=gpar(col = type_cols[l], fill = type_cols[l]))
          popViewport()
        }
      popViewport()
  
      # plot subcluster legend:
      pushViewport(viewport(x = unit(4, "cm"), y = unit(15.7, "cm"), width = unit(6, "cm"), 
        height = unit(10, "cm")))
        # add title:
        pushViewport(viewport(x = 0.06, y = 1, width = unit(2, "cm"), height = unit(0.5, "cm")))
          grid.text("CNA subclusters", gp=gpar(fontsize=20), just = "left")
        popViewport()
        for (l in 1:length(all_annotations[[i]]$subcluster_cols)) {
          # add labels:
          pushViewport(viewport(
            x = 0.21, 
            y = 0.99-(0.08*l), 
            width = unit(2, "cm"), 
            height = unit(0.5, "cm")
          ))
            #grid.rect()
            grid.text(
              names(all_annotations[[i]]$subcluster_cols)[l], 
              gp=gpar(fontsize=18), 
              just = "left"
            )
          popViewport()
          # add dots:
          pushViewport(viewport(
            x = 0.1, 
            y = 0.99-(0.08*l), 
            width = unit(0.3, "cm"), 
            height = unit(0.3, "cm")
          ))
            grid.circle(
              gp=gpar(col = all_annotations[[i]]$subcluster_cols[l], 
              fill = all_annotations[[i]]$subcluster_cols[l])
            )
          popViewport()
        }
      popViewport()
  
      # add GIN and QC annotation values:
      pushViewport(viewport(
        x=final_heatmap_data[[i]]$x_coord + 0.885, y=0.083, 
        width = 0.1, height = 0.1, 
        just = "top")
      )
        grid.text("GIN", rot=65, gp=gpar(fontsize=20))
      popViewport()
      pushViewport(viewport(
        x=final_heatmap_data[[i]]$x_coord + 0.905, y=0.083, 
        width = 0.1, height = 0.1, 
        just = "top")
      )
        grid.text("nUMI", rot=65, gp=gpar(fontsize=20))
      popViewport()
      pushViewport(viewport(
        x=final_heatmap_data[[i]]$x_coord + 0.925, y=0.081, 
        width = 0.1, height = 0.1, 
        just = "top")
      )
        grid.text("nGene", rot=65, gp=gpar(fontsize=20))
      popViewport()
  
  dev.off()

  #convert pdf to png:
  system(paste0("for p in ", out_dir, "*.pdf; do echo $p; f=$(basename $p); echo $f; ",
  "new=$(echo $f | sed 's/.pdf/.png/'); echo $new; ", 
  "convert -density 150 ", out_dir, "$f -quality 90 ", out_dir, "$new; done"))

  print(paste0("Heatmap created, output in ", out_dir))
}


