#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

project_name <- "thesis"
subproject_name <- "chapter_4"
args = commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
subcluster_method <- args[2]
subcluster_p <- args[3]
if (subcluster_p != "none") {
  subcluster_p <- as.numeric(subcluster_p)
}
coverage_filter <- args[4]
remove_artefacts <- args[5]
met <- args[6]
res <- args[7]
PC <- args[8]
epi_res <- args[9]
epi_PC <- args[10]
garnett_slot <- args[11]
remove_outliers <- args[12]
outlier_sd_multiplier <- as.numeric(args[13])
normal_markers <- strsplit( 
  args[14],
  "_"
)[[1]]
epi_markers <- strsplit( 
  args[15],
  "_"
)[[1]]
custom_epi_markers <- strsplit( 
  args[16],
  "_"
)[[1]]


project_name <- "thesis"
subproject_name <- "chapter_4"
sample_name <- "CID44971"
subcluster_method <- "random_trees"
subcluster_p <- "0.05"
if (subcluster_p != "none") {
  subcluster_p <- as.numeric(subcluster_p)
}
coverage_filter <- "filtered"
remove_artefacts <- "artefacts_not_removed"
met <- FALSE
res <- "PC_C_res.1"
PC <- "C"
epi_res <- "PC_C_res.1"
epi_PC <-"C"
garnett_slot <- "garnett_call_ext_major"
remove_outliers <- TRUE
outlier_sd_multiplier <- 3
normal_markers <- strsplit( 
  "EPCAM_KRT18_ESR1_KRT5_KRT14_PGR_ERBB2_MKI67",
  "_"
)[[1]]
epi_markers <- strsplit( 
  "EPCAM_KRT18_ESR1_KRT5_KRT14_PGR_ERBB2_MKI67",
  "_"
)[[1]]
custom_epi_markers <- strsplit(
  "MUCL1_NDRG1_IGFBP5",
  "_"
)[[1]]
#ELF5_GATA3_CD44_ZEB1

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("Subclustered by ", subcluster_method))
print(paste0("Subclustered by pval ", subcluster_p))
print(paste0("Filter by coverage? ", coverage_filter))
print(paste0("Remove artefacts? ", remove_artefacts))
print(paste0("Resolution for all cell plots = ", res))
print(paste0("PC for all cell plots = ", PC))
print(paste0("Resolution for epithelial cell plots = ", epi_res))
print(paste0("PC for epithelial cell plots = ", epi_PC))
print(paste0("Garnett slot = ", garnett_slot))
print(paste0("Remove outliers? ", remove_outliers))
print(paste0("No SDs within which to remove outliers ", outlier_sd_multiplier))
print(paste0("Normal markers ", normal_markers))
print(paste0("Epithelial markers ", epi_markers))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(naturalsort, lib.loc = lib_loc)
library(Seurat)
library(ggplot2)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/",
  subproject_name, "/")
ref_dir <- paste0(project_dir, "/refs/")
func_dir <- paste0(project_dir, "/scripts/functions/")
results_dir <- paste0(project_dir, "results/")
seurat_dir <- paste0(project_dir, "raw_files/seurat_objects/", 
  sample_name, "/")
in_path <- paste0(
  results_dir, "infercnv/", sample_name, "/",
  coverage_filter, "/", subcluster_method, "/p_", subcluster_p, "/", 
  remove_artefacts, "/"
)

Robject_dir <- paste0(in_path, "Rdata/")
system(paste0("mkdir -p ", Robject_dir))
plot_dir <- paste0(in_path, "plots/expression_plots/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(in_path, "tables/")
system(paste0("mkdir -p ", table_dir))

print(paste0("In path = ", in_path))
print(paste0("Reference directory = ", ref_dir))
print(paste0("R function directory = ", func_dir))
print(paste0("R object directory = ", Robject_dir))
print(paste0("Table directory = ", table_dir))
print(paste0("Plot directory = ", plot_dir))


################################################################################
### 0. Define colours and load functions###
################################################################################

expr_cols <- read.table(
  paste0(ref_dir, "expression_colour_palette.txt"),
  header = F,
  stringsAsFactors = F,
  comment.char = ""
)[,1]

normal_col_df <- read.table(
  paste0(ref_dir, "normal_colour_palette.txt"),
  header = F,
  stringsAsFactors = F,
  comment.char = ""
)
normal_cols <- normal_col_df[,2]
names(normal_cols) <- normal_col_df[,1]

subcluster_cols <- read.table(
  paste0(ref_dir, "CNV_colour_palette.txt"),
  header = F,
  stringsAsFactors = F,
  comment.char = ""
)[,1]

remove_UMAP_outliers <- dget(paste0(func_dir, "remove_UMAP_outliers.R"))


################################################################################
### 1. Load and format data ###
################################################################################

print("Loading metadata df...")

if (met) {
  epithelial_metadata <- readRDS(paste0(Robject_dir, 
    "/5b.final_epithelial_metadata_without_normals.Rdata"))
} else {
  epithelial_metadata <- readRDS(paste0(Robject_dir, 
    "/4b.final_epithelial_metadata_with_normals.Rdata"))
}

print("Loading seurat objects...")

# load all cell object:
seurat_10X <- readRDS(
  paste0(
    seurat_dir, "03_seurat_object_processed.Rdata"
  )
)
  
# choose resolution if needed:
if (res != "none") {
  Idents(seurat_10X) <- eval(
    parse(
      text = paste0("seurat_10X@meta.data$", res)
    )
  )
}

# update idents will cell type from Garnett:
annotated_idents <- paste0(
  eval(parse(text=paste0("seurat_10X@meta.data$", garnett_slot))), " ", 
  Idents(seurat_10X)
)
Idents(seurat_10X) <- factor(
  as.character(annotated_idents),
  levels = naturalsort(unique(as.character(annotated_idents)))
)

if (
  file.exists(
    paste0(
      seurat_dir, "05_seurat_object_epithelial_reclustered.Rdata"
    )
  )
) {
  # load reclustered epithelial object:
  seurat_epi <- readRDS(
    paste0(
      seurat_dir, "05_seurat_object_epithelial_reclustered.Rdata"
    )
  )
  
  # choose resolution if needed:
  if (epi_res != "none") {
    Idents(seurat_epi) <- eval(
      parse(
        text = paste0("seurat_epi@meta.data$", epi_res)
      )
    )
  }

  # only keep cells in metadata:
  seurat_epi <- subset(
    seurat_epi,
    cells = epithelial_metadata$cell_ids
  )

#  # label normal and unassigned cells:
#  temp_idents <- as.character(Idents(seurat_epi))
#  names(temp_idents) <- names(Idents(seurat_epi))
#
#  temp_idents[
#    names(temp_idents) %in% epithelial_metadata$cell_ids[
#      epithelial_metadata$normal_cell_call == "normal"
#    ]
#  ] <- "normal"
#
#  temp_idents[
#    names(temp_idents) %in% epithelial_metadata$cell_ids[
#      epithelial_metadata$normal_cell_call == "unassigned"
#    ]
#  ] <- "unassigned"
#
#  Idents(seurat_epi) <- factor(
#    temp_idents,
#    levels = naturalsort(unique(temp_idents))
#  )
#
#  # reassign levels:
#  if ("normal" %in% levels(Idents(seurat_epi)) & 
#  	"unassigned" %in% levels(Idents(seurat_epi))) {
#  	levels(Idents(seurat_epi)) <- c(
#      paste0(
#        "Expression_",
#        1:(length(levels(Idents(seurat_epi)))-2)
#      ),
#      "normal",
#      "unassigned"
#    )
#  } else if ("normal" %in% levels(Idents(seurat_epi))) {
#  	levels(Idents(seurat_epi)) <- c(
#      paste0(
#        "Expression_",
#        1:(length(levels(Idents(seurat_epi)))-1)
#      ),
#      "normal"
#    )
#  } else {
#  	levels(Idents(seurat_epi)) <- paste0(
#      "Expression_",
#      1:(length(levels(Idents(seurat_epi))))
#    )
#  }
  
}


################################################################################
### 2. Plot UMAP of all cells ###
################################################################################

# order idents so epithelial cells are printed on top:
ident_order <- levels(Idents(seurat_10X))
ident_order <- rev(
  c(
    ident_order[grep("Epithelial", ident_order)],
    ident_order[grep("Epithelial", ident_order, invert = T)]
  )
)

if (!file.exists(paste0(plot_dir, "all_cells_UMAP.png"))) {

  all_umap <- DimPlot(
    seurat_10X,
    cols = expr_cols,
    pt.size = 1.5,
    reduction = paste0("UMAP", PC),
    label = F,
    order = ident_order
  )
  
  png(
    file = paste0(plot_dir, "all_cells_UMAP.png"), 
    width = 15, 
    height = 8, 
    res = 300, 
    units = 'in'
  )
    print(all_umap)
  dev.off()

}

if (!file.exists(paste0(plot_dir, "all_cells_UMAP_epi_highlight.png"))) {

  all_umap_epi_highlight <- DimPlot(
    seurat_10X,
    cols = expr_cols,
    pt.size = 1.5,
    reduction = paste0("UMAP", PC),
    label = F,
    order = ident_order,
    cells.highlight = names(Idents(seurat_10X))[
      Idents(seurat_10X) %in% ident_order[grep("pithelial", ident_order)]
    ]
  )

  png(
    file = paste0(plot_dir, "all_cells_UMAP_epi_highlight.png"), 
    width = 15, 
    height = 8, 
    res = 300, 
    units = 'in'
  )
    print(all_umap_epi_highlight)
  dev.off()

}

if (
  file.exists(
    paste0(
      seurat_dir, "05_seurat_object_epithelial_reclustered.Rdata"
    )
  )
) {

  ################################################################################
  ### 3. Plot UMAP of reclustered epithelial cells ###
  ################################################################################
  
  if (
    !file.exists(paste0(
      plot_dir, "epithelial_cell_expr_clusters_UMAP_pre_filter.png")
    )
  ) {
  
    epi_umap <- DimPlot(
      seurat_epi_no_normal_annots,
      cols = expr_cols,
      pt.size = 1.5,
      reduction = paste0("UMAP", epi_PC),
      label = F
    ) + theme(
      axis.title.y = element_text(
        size=25,
        margin = margin(t = 0, r = 20, b = 0, l = 0)
      ),
      axis.text.y = element_text(size=20),
      axis.title.x = element_text(
        size=25,
        margin = margin(t = 20, r = 0, b = 0, l = 0)
      ),
      axis.text.x = element_text(size=20),
      legend.text = element_text(size=20)
    )
  
    png(
      file = paste0(plot_dir, "epithelial_cell_expr_clusters_UMAP_pre_filter.png"), 
      width = 15, 
      height = 8, 
      res = 300, 
      units = 'in'
    )
      print(epi_umap)
    dev.off()
  
  }

  # remove outliers - defined as those > 3 standard devs from mean:
  temp_seurat_epi <- remove_UMAP_outliers(seurat_epi)
  filtered_seurat_epi <- temp_seurat_epi[[1]]
  no_outlier_seurat_epi <- temp_seurat_epi[[2]]
  include_cells <- temp_seurat_epi[[3]]

#  if (!file.exists(
#    paste0(seurat_dir, "06_seurat_object_filtered_epithelial.Rdata"))
#  ) {
#    saveRDS(
#      filtered_seurat_epi,
#      paste0(seurat_dir, "06_seurat_object_filtered_epithelial.Rdata")
#    )
#  }
#
#  if (!file.exists(
#    paste0(seurat_dir, "07_seurat_object_no_outlier_epithelial.Rdata"))
#  ) {
#    saveRDS(
#      no_outlier_seurat_epi,
#      paste0(seurat_dir, "07_seurat_object_no_outlier_epithelial.Rdata")
#    )
#  }
  
  if (
    !file.exists(paste0(
      plot_dir, "epithelial_cell_expr_clusters_UMAP.png")
    )
  ) {
  
    epi_umap <- DimPlot(
      filtered_seurat_epi,
      cols = expr_cols,
      pt.size = 1.5,
      reduction = paste0("UMAP", epi_PC),
      label = F
    ) + theme(
      axis.title.y = element_text(
        size=25,
        margin = margin(t = 0, r = 20, b = 0, l = 0)
      ),
      axis.text.y = element_text(size=20),
      axis.title.x = element_text(
        size=25,
        margin = margin(t = 20, r = 0, b = 0, l = 0)
      ),
      axis.text.x = element_text(size=20),
      legend.text = element_text(size=20)
    )
  
    png(
      file = paste0(plot_dir, "epithelial_cell_expr_clusters_UMAP.png"), 
      width = 15, 
      height = 8, 
      res = 300, 
      units = 'in'
    )
      print(epi_umap)
    dev.off()
  
  }
  
  if (
    !file.exists(paste0(
      plot_dir, "epithelial_cell_expr_clusters_UMAP_no_outlier.png")
    )
  ) {
  
    epi_umap <- DimPlot(
      no_outlier_seurat_epi,
      cols = expr_cols,
      pt.size = 1.5,
      reduction = paste0("UMAP", epi_PC),
      label = F
    ) + theme(
      axis.title.y = element_text(
        size=25,
        margin = margin(t = 0, r = 20, b = 0, l = 0)
      ),
      axis.text.y = element_text(size=20),
      axis.title.x = element_text(
        size=25,
        margin = margin(t = 20, r = 0, b = 0, l = 0)
      ),
      axis.text.x = element_text(size=20),
      legend.text = element_text(size=20)
    )
  
    png(
      file = paste0(plot_dir, "epithelial_cell_expr_clusters_UMAP_no_outlier.png"), 
      width = 15, 
      height = 8, 
      res = 300, 
      units = 'in'
    )
      print(epi_umap)
    dev.off()
  
  }
  
  
  ################################################################################
  ### 4. Plot UMAP with normals labelled ###
  ################################################################################
  
  if (met) {

    filtered_seurats <- list(filtered_seurat_epi, no_outlier_seurat_epi)
  
    # update idents with subcluster ids:
    for (i in 1:length(filtered_seurats)) {

      epi_feature <- FeaturePlot(
        object = filtered_seurats[[i]],
        pt.size = 1.5,
        reduction = paste0("UMAP", epi_PC),
        features = c(epi_markers, custom_epi_markers),
        order = T
      )
      dev.off()
      
      if (i==1) {
        png(
          file = paste0(plot_dir, "epithelial_feature_UMAPs.png"), 
          width = 15, 
          height = 8, 
          res = 300, 
          units = 'in'
        )
      } else {
        png(
          file = paste0(plot_dir, "epithelial_feature_UMAPs_no_outlier.png"), 
          width = 15, 
          height = 8, 
          res = 300, 
          units = 'in'
        )
      }
        print(epi_feature)
      dev.off()

    }

  } else {

    filtered_seurats <- list(filtered_seurat_epi, no_outlier_seurat_epi)
  
    # update idents with subcluster ids:
    for (i in 1:length(filtered_seurats)) {
    
      normal_lab_seurat <- filtered_seurats[[i]]
      normal_lab_seurat@meta.data$normal_cell_call <- NA
      normal_lab_seurat@meta.data[epithelial_metadata$cell_ids,]$normal_cell_call <- 
        as.character(epithelial_metadata$normal_cell_call)
      
      Idents(normal_lab_seurat) <- factor(
        as.character(normal_lab_seurat@meta.data$normal_cell_call),
        levels = naturalsort(
          unique(as.character(normal_lab_seurat@meta.data$normal_cell_call))
        )
      )
      
      labs <- c("normal", "unassigned", "cancer")
      labs <- labs[
        labs %in% gsub(
          "_", " ", unique(
            Idents(normal_lab_seurat)[
              !is.na(Idents(normal_lab_seurat))
            ]
          )
        )
      ]
      
      normal_lab_umap <- DimPlot(
        normal_lab_seurat,
        pt.size = 1.5,
        reduction = paste0("UMAP", epi_PC),
        label = F,
        cols = normal_cols
      ) + theme(
          axis.title.y = element_text(
            size=25,
            margin = margin(t = 0, r = 20, b = 0, l = 0)
          ),
          axis.text.y = element_text(size=20),
          axis.title.x = element_text(
            size=25,
            margin = margin(t = 20, r = 0, b = 0, l = 0)
          ),
          axis.text.x = element_text(size=20),
          legend.text = element_text(size=20)
        )
      
      if (i==1) {
        png(
          file = paste0(plot_dir, "epithelial_cell_normal_call_UMAP.png"), 
          width = 15, 
          height = 8, 
          res = 300, 
          units = 'in'
        )
      } else {
        png(
          file = paste0(plot_dir, "epithelial_cell_normal_call_UMAP_no_outlier.png"), 
          width = 15, 
          height = 8, 
          res = 300, 
          units = 'in'
        )   
      }
        print(normal_lab_umap)
      dev.off()
    
    
      ################################################################################
      ### 5. Plot feature plots of epithelial marker expression ###
      ################################################################################
      
      normal_feature <- FeaturePlot(
        object = normal_lab_seurat,
        pt.size = 1.5,
        reduction = paste0("UMAP", epi_PC),
        features = normal_markers,
        order = T
      )
      dev.off()
      
      if (i==1) {
        png(
          file = paste0(plot_dir, "normal_feature_UMAPs.png"), 
          width = 15, 
          height = 8, 
          res = 300, 
          units = 'in'
        )
      } else {
        png(
          file = paste0(plot_dir, "normal_feature_UMAPs_no_outlier.png"), 
          width = 15, 
          height = 8, 
          res = 300, 
          units = 'in'
        )
      }
        print(normal_feature)
      dev.off()
  
      epi_feature <- FeaturePlot(
        object = normal_lab_seurat,
        pt.size = 1.5,
        reduction = paste0("UMAP", epi_PC),
        features = c(epi_markers, custom_epi_markers),
        order = T
      )
      dev.off()
      
      if (i==1) {
        png(
          file = paste0(plot_dir, "epithelial_feature_UMAPs.png"), 
          width = 15, 
          height = 8, 
          res = 300, 
          units = 'in'
        )
      } else {
        png(
          file = paste0(plot_dir, "epithelial_feature_UMAPs_no_outlier.png"), 
          width = 15, 
          height = 8, 
          res = 300, 
          units = 'in'
        )
      }
        print(epi_feature)
      dev.off()
    
    }
    
    
    ################################################################################
    ### 6. Plot UMAP of reclustered epithelial cells without normals ###
    ################################################################################
    
    # remove normals:
    if (remove_outliers) {
      no_normal <- include_cells[
        include_cells %in% epithelial_metadata$cell_ids[
          epithelial_metadata$normal_cell_call == "cancer"
        ]
      ]
    } else {
      no_normal <- no_outlier[
        no_outlier %in% epithelial_metadata$cell_ids[
          epithelial_metadata$normal_cell_call == "cancer"
        ]
      ]
    }
    seurat_epi_no_normal <- subset(
      no_outlier_seurat_epi,
      cells = no_normal
    )
  
    # remove outliers:
    seurat_epi_no_normal_filtered <- remove_UMAP_outliers(seurat_epi_no_normal)[[2]]
    
#    if (!file.exists(
#      paste0(seurat_dir, "08_seurat_object_no_normal_epithelial.Rdata"))
#    ) {
#      saveRDS(
#        seurat_epi_no_normal_filtered,
#        paste0(seurat_dir, "08_seurat_object_no_normal_epithelial.Rdata")
#      )
#    }
    
    if (
      !file.exists(paste0(
        plot_dir, "epithelial_cell_expr_clusters_UMAP_no_normal.png")
      )
    ) {
    
      epi_umap <- DimPlot(
        seurat_epi_no_normal_filtered,
        cols = expr_cols,
        pt.size = 1.5,
        reduction = paste0("UMAP", epi_PC),
        label = F
      )
      
      png(
        file = paste0(plot_dir, "epithelial_cell_expr_clusters_UMAP_no_normal.png"), 
        width = 15, 
        height = 8, 
        res = 300, 
        units = 'in'
      )
        print(epi_umap)
      dev.off()
    
    }
  
    
    ################################################################################
    ### 7. Plot feature plots of epithelial marker expression without normals ###
    ################################################################################
    
    epi_feature <- FeaturePlot(
      object = seurat_epi_no_normal_filtered,
      pt.size = 1.5,
      reduction = paste0("UMAP", epi_PC),
      features = c(epi_markers, custom_epi_markers),
      order = T
    )
    dev.off()
    
    png(
      file = paste0(plot_dir, "epithelial_feature_UMAPs_no_normal.png"), 
      width = 15, 
      height = 8, 
      res = 300, 
      units = 'in'
    )
      print(epi_feature)
    dev.off()

  }
  
  
  ################################################################################
  ### 7. Plot UMAP with CNV subclusters labelled ###
  ################################################################################
  
  if (met) {
    seurat_epi_no_normal_filtered <- no_outlier_seurat_epi
  }

  rownames(epithelial_metadata) <- epithelial_metadata$cell_ids
  # update idents with subcluster ids:
  seurat_subcluster <- seurat_epi_no_normal_filtered
  seurat_subcluster@meta.data$subcluster_id <- NA
  temp_meta <- epithelial_metadata[rownames(seurat_subcluster@meta.data),]

  seurat_subcluster@meta.data$subcluster_id <- 
    as.character(temp_meta$subcluster_id)
  
  Idents(seurat_subcluster) <- factor(
    as.character(seurat_subcluster@meta.data$subcluster_id),
    levels = naturalsort(unique(as.character(seurat_subcluster@meta.data$subcluster_id)))
  )
  
  labs <- gsub(
    "_", " ", levels(
      Idents(seurat_subcluster)[
        !is.na(Idents(seurat_subcluster))
      ]
    )
  )
  
  epi_umap <- DimPlot(
    seurat_subcluster,
    pt.size = 1.5,
    reduction = paste0("UMAP", epi_PC),
    label = F
  ) + 
  scale_color_manual(
    labels = labs, 
    values = subcluster_cols
  )
  
  png(
    file = paste0(plot_dir, "epithelial_cell_CNV_subclusters_UMAP.png"), 
    width = 15, 
    height = 8, 
    res = 300, 
    units = 'in'
  )
    print(epi_umap)
  dev.off()

}

