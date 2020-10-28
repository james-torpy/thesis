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
subset_samples <- as.logical(args[6])
pos_only <- as.logical(args[7])
external_list <- as.logical(args[8])
external_list_name <- as.logical(args[9])

project_name <- "thesis"
subproject_name <- "chapter_4"
args = commandArgs(trailingOnly=TRUE)
subcluster_method <- "random_trees"
subcluster_p <- "0.05"
if (subcluster_p != "none") {
  subcluster_p <- as.numeric(subcluster_p)
}
coverage_filter <- "filtered"
remove_artefacts <- "artefacts_not_removed"
subset_samples <- FALSE
pos_only <- FALSE
min_common_prop <- 0.2
lower_min_common_prop <- 0.1
external_list <- "none"
external_list_name <- "none"
#external_list <- c("ANGPT1", "CASP1", "FGFR1", "GM19589", "Il1R1",
#  "MGP", "MFSS1", "NEURL1B", "PHEX", "PLSCR4", "SERPINE2", "SLPI")
#external_list_name <- "Wagenblast"

up_met_sets <- c(
  "LIAO_METASTASIS", "NAKAMURA_METASTASIS",
  "PEDERSEN_METASTASIS_BY_ERBB2_ISOFORM_1",
  "PEDERSEN_METASTASIS_BY_ERBB2_ISOFORM_3",
  "PEDERSEN_METASTASIS_BY_ERBB2_ISOFORM_4",
  "PEDERSEN_METASTASIS_BY_ERBB2_ISOFORM_5",
  "PEDERSEN_METASTASIS_BY_ERBB2_ISOFORM_6",
  "PEDERSEN_METASTASIS_BY_ERBB2_ISOFORM_7",
  "TAVAZOIE_METASTASIS"
)
down_met_sets <- c(
  "GILDEA_METASTASIS",
  "WANG_METASTASIS_OF_BREAST_CANCER"
)

msig_brca_met_sets <- c(
  "PEDERSEN_METASTASIS_BY_ERBB2_ISOFORM_1",
  "PEDERSEN_METASTASIS_BY_ERBB2_ISOFORM_3",
  "PEDERSEN_METASTASIS_BY_ERBB2_ISOFORM_4",
  "PEDERSEN_METASTASIS_BY_ERBB2_ISOFORM_5",
  "PEDERSEN_METASTASIS_BY_ERBB2_ISOFORM_6",
  "PEDERSEN_METASTASIS_BY_ERBB2_ISOFORM_7",
  "TAVAZOIE_METASTASIS", "VANTVEER_BREAST_CANCER_METASTASIS_DN", 
  "VANTVEER_BREAST_CANCER_METASTASIS_UP",
  "WANG_METASTASIS_OF_BREAST_CANCER",
  "WANG_METASTASIS_OF_BREAST_CANCER_ESR1_UP",
  "WANG_METASTASIS_OF_BREAST_CANCER_ESR1_DN"
)

if (subset_samples) {

  # determine order of samples by subtype:
  ER <- c("CID4463")
  HER2 <- c("CID45171")
  TNBC <- c("CID44991")
  sample_names <- c(ER, HER2, TNBC)

  subtype_df <- data.frame(
    subtype = c(
      rep("ER", length(ER)),
      rep("HER2", length(HER2)),
      rep("TNBC", length(TNBC))
    ),
    sample = sample_names
  )

} else {

  # determine order of samples by subtype:
  ER <- c("CID3941", "CID4067", "CID4290A", "CID4461", 
    "CID4463", "CID4530N", "CID4535")
  HER2 <- c("CID3921", "CID3586", "CID3963", "CID4066", "CID45171")
  TNBC <- c("CID44971", "CID44991", "CID4513",
    "CID4515", "CID4523")
  sample_names <- c(ER, HER2, TNBC)

  subtype_df <- data.frame(
    subtype = c(
      rep("ER", length(ER)),
      rep("HER2", length(HER2)),
      rep("TNBC", length(TNBC))
    ),
    sample = sample_names
  )

}

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(ComplexHeatmap, lib.loc = lib_loc)
library(circlize, lib.loc = lib_loc)
library(dplyr)
library(textshape, lib.loc = lib_loc)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/",
  subproject_name, "/")
ref_dir <- paste0(project_dir, "/refs/")
func_dir <- paste0(project_dir, "/scripts/functions/")
results_dir <- paste0(project_dir, "results/")
in_path <- paste0(results_dir, "infercnv/")
out_dir <- paste0(results_dir, "infercnv/common_gene_expression/", 
  coverage_filter, "/", subcluster_method, "/p_", subcluster_p, "/", 
  remove_artefacts, "/")

if (subset_samples) {
  Robject_dir <- paste0(out_dir, "Rdata_sub/")
  plot_dir <- paste0(out_dir, "plots_sub/")
  table_dir <- paste0(out_dir, "tables_sub/")
} else {
  Robject_dir <- paste0(out_dir, "Rdata/")
  plot_dir <- paste0(out_dir, "plots/")
  table_dir <- out_dir
}

system(paste0("mkdir -p ", Robject_dir))
system(paste0("mkdir -p ", plot_dir))
system(paste0("mkdir -p ", table_dir))

print(paste0("In path = ", in_path))
print(paste0("Out dir = ", out_dir))
print(paste0("Reference directory = ", ref_dir))
print(paste0("R function directory = ", func_dir))
print(paste0("R object directory = ", Robject_dir))
print(paste0("Plot directory = ", plot_dir))


################################################################################
### 0. Load functions and colour palette ###
################################################################################

filter_common_DE_genes <- dget(paste0(func_dir, "filter_common_DE_genes.R"))

subtype_cols <- read.table(
  paste0(ref_dir, "subtype_colour_palette.txt"),
  header = F,
  comment.char = "",
  stringsAsFactors = F
)
m <- match(c("ER", "HER2", "TNBC"), subtype_cols$V1)
subtype_cols <- subtype_cols$V2[m]
names(subtype_cols) <- c("ER", "HER2", "TNBC")


################################################################################
### 1. Load DE data ###
################################################################################

if (!file.exists(paste0(Robject_dir, "common_subpop_DE_data.Rdata"))) {

  for (s in 1:length(sample_names)) {

    print(paste0("Loading ", sample_names[s], " CNV data..."))

    sample_table_dir <- paste0(
  	  in_path, sample_names[s], "/", coverage_filter, "/", subcluster_method, 
      "/p_", subcluster_p, "/", remove_artefacts, 
      "/tables/"
    )

    DE_data <- read.table(
      paste0(
        sample_table_dir, 
        "top_CNA_assoc_subpop_DE.txt"
      ),
      header = T
    )

    all_gene_DE_data <- read.table(
      paste0(
      	sample_table_dir, 
      	"all_CNA_assoc_subpop_DE.txt"
      ),
      header = T
    )
    all_gene_DE_data$sample_name <- sample_names[s] 

    # only keep genes with adjusted p value < 0.1:
    filtered_DE <- DE_data[DE_data$p_val_adj < 0.1,]

    if (s==1) {

      all_DE_data <- list(filtered_DE)
      names(all_DE_data) <- sample_names[s]

      all_gene_data <- list(all_gene_DE_data)
      names(all_gene_data) <- sample_names[s]

    } else {

      all_DE_data[[s]] <- filtered_DE
      names(all_DE_data)[s] <- sample_names[s]

      all_gene_data[[s]] <- all_gene_DE_data
      names(all_gene_data)[s] <- sample_names[s]

    }

  }

 saveRDS(
    all_DE_data, 
    paste0(Robject_dir, "common_subpop_DE_data.Rdata")
  )

  saveRDS(
    all_gene_data, 
    paste0(Robject_dir, "common_all_gene_subpop DE_data.Rdata")
  )

} else {

  all_DE_data <- readRDS(
    paste0(Robject_dir, "common_subpop_DE_data.Rdata")
  )
  all_gene_data <- readRDS(
    paste0(Robject_dir, "common_all_gene_subpop DE_data.Rdata")
  )
  
}


################################################################################
### 2. Identify DE genes present in a proportion of samples ###
################################################################################

# find genes present in all elements of all_DE_data list:
common_DE <- unique(
  unlist(
    lapply(all_DE_data, function(x) {
      gene_vec <- as.character(
        unlist(
          lapply(all_DE_data, function(y) {
            return(
              unique(
                x$gene[
                  which(x$gene %in% y$gene)
                ]
              )
            )
          })
        )
      )
      gene_tab <- table(gene_vec)
      return(
        names(gene_tab)[
          which(gene_tab >= ceiling(length(all_DE_data)*min_common_prop))
        ]
      )
    })
  )
)


################################################################################
### 3. Prepare DE genes present in a proportion of samples ###
################################################################################

# merge together met DE data:
if (pos_only) {

  DE_data <- lapply(
    all_gene_data, 
    filter_common_DE_genes, 
    common_DE, 
    "up",
    TRUE,
    func_dir,
    ref_dir,
    Robject_dir
  )

  empty_df <- data.frame(
    gene = sort(
      common_DE
    ),
    stringsAsFactors = FALSE
  )

  plot_data <- list(
    avg_logFC = empty_df,
    p_val_adj = empty_df
  )

} else {

  DE_data <- lapply(
    all_gene_data, 
    filter_common_DE_genes, 
    common_DE, 
    "bi",
    TRUE,
    func_dir,
    ref_dir,
    Robject_dir
  )

  # keep only common_DE genes:
  empty_df <- data.frame(
    gene = sort(
      c(
        paste0(common_DE, "_down"),
        paste0(common_DE, "_up")
      )
    ),
    stringsAsFactors = FALSE
  )

  plot_data <- list(
    avg_logFC = empty_df,
    p_val_adj = empty_df
  )

}

for (l in 1:length(plot_data)) {

  for (s in 1:length(DE_data)) {
    temp_df <- cbind(
      DE_data[[s]]$gene, 
      as.data.frame(
        eval(parse(text = paste0("DE_data[[s]]$", names(plot_data)[l])))
      )
    )
    colnames(temp_df) <- c("gene", names(plot_data)[l])
    plot_data[[l]] <- merge(
      plot_data[[l]],
      temp_df,
      by = "gene",
      all.x = TRUE
    )
    colnames(plot_data[[l]])[s+1] <- names(DE_data)[s]
  }

  # remove rows with all NAs:
  plot_data[[l]] <- plot_data[[l]][
    apply(plot_data[[l]], 1, function(x) {
      !(all(is.na(x)))
    }),
  ]

  plot_data[[l]] <- plot_data[[l]] %>%
    column_to_rownames(loc = 1)

}

# remove genes DE in < min proportion required:
plot_data$avg_logFC <- plot_data$avg_logFC[
  apply(plot_data$avg_logFC, 1, function(x) {
    length(which(!is.na(x))/length(x)) >= min_common_prop
  }),
]

# ensure p_vals rows match to logFC:
plot_data$p_val_adj <- plot_data$p_val_adj[rownames(plot_data$avg_logFC),]
print(paste0(
  "Are logFC rownames identical to p_vals rownames? ", 
  identical(rownames(plot_data$avg_logFC), rownames(plot_data$p_val_adj)))
)
# change significant values to the ComplexHeatmap code for dots, and all
# others to NA:
plot_data$p_val_dots <- plot_data$p_val_adj
plot_data$p_val_dots[plot_data$p_val_dots < 0.1] <- "*"
plot_data$p_val_dots[plot_data$p_val_dots != "*"] <- " "
plot_data$p_val_dots[is.na(plot_data$p_val_dots)] <- " "


################################################################################
### 3. Plot common DE genes ###
################################################################################

# split into up and down:
all_data <- list(
  up_data = lapply(plot_data, function(x) {
    return(x[grep("up", rownames(x)),])
  }),
  down_data = lapply(plot_data, function(x) {
    return(x[grep("down", rownames(x)),])
  })
)

for (i in 1:length(all_data)) {
  # define heatmap colours:
  na_less_vector <- unlist(all_data[[i]]$avg_logFC)
  na_less_vector <- na_less_vector[!is.na(na_less_vector)]

  if (i==1) {
    heatmap_cols <- colorRamp2(c(0, max(na_less_vector)), 
      c("white", "#870C0C"), space = "sRGB")
  } else {
    heatmap_cols <- colorRamp2(c(0, min(na_less_vector)), 
      c("white", "#0D3084"), space = "sRGB")
  }

  mat <- as.matrix(all_data[[i]]$avg_logFC)
  rownames(mat) <- gsub("_.*$", "", rownames(mat))
  pdots <- all_data[[i]]$p_val_dots

  if (i==1) {
    star_col <- "#1F0B99"
  } else {
    star_col <- "#870C0C"
  }

  # order subtype_df the same as matrix column:
  m <- match(colnames(mat), subtype_df$sample)
  subtype_df <- subtype_df[m,]

  # create subtype annotation:
  subtype_annot <- HeatmapAnnotation(
    subtype = as.character(subtype_df$subtype),
    col = list(
      subtype = c(
        subtype_cols["ER"], 
        subtype_cols["HER2"], 
        subtype_cols["TNBC"]
      )
    ),
    show_legend = F,
    show_annotation_name = F
  )

  # create heatmap legend:
  if (i==1) {

    lgd <- Legend(
      at = seq(
        round(min(na_less_vector), 1), 
        round(max(na_less_vector), 1),
        1
      ),
      col_fun = heatmap_cols, 
      title = "logFC", 
      direction = "vertical",
      grid_height = unit(1, "cm"),
      grid_width = unit(1, "cm"),
      #legend_height = unit(4.5, "cm"),
      #legend_width = unit(2, "cm"),
      labels_gp = gpar(fontsize = 20),
      title_gp = gpar(fontsize = 22, fontface = "plain")
    )

  } else {

   lgd <- Legend(
      at = seq(
        round(min(na_less_vector), 1), 
        round(max(na_less_vector), 1),
        1
      ),
      col_fun = heatmap_cols, 
      title = "logFC", 
      direction = "vertical",
      grid_height = unit(1, "cm"),
      grid_width = unit(1, "cm"),
      #legend_height = unit(4.5, "cm"),
      #legend_width = unit(2, "cm"),
      labels_gp = gpar(fontsize = 20),
      title_gp = gpar(fontsize = 22, fontface = "plain")
    )

  }
   
  # create main heatmap:
  final_heatmap <- Heatmap(
    mat, 
    na_col = "grey",
    name = paste0("hm"), 
    col = heatmap_cols,
    cluster_columns = F, cluster_rows = F,
    show_row_names = T, show_column_names = T,
    row_names_gp = gpar(fontsize = 17),
    column_names_gp = gpar(fontsize = 17),
    show_row_dend = T,
    show_heatmap_legend = F,
    use_raster = T, raster_device = c("png"),
    cell_fun = function(j, i, x, y, w, h, heatmap_cols) { # add text to each grid
        grid.text(pdots[i, j], x, y, gp=gpar(fontsize=30, col = star_col))
    },
    bottom_annotation = subtype_annot
  )

  # make grid object:
  ht_list <- final_heatmap
  annotated_heatmap <- grid.grabExpr(
    draw(ht_list
      #, gap = unit(6, "mm"), heatmap_legend_side = "left"
    )
  )
  dev.off()
  
  if (i==1) {
    pdf(paste0(plot_dir, "CNA_assoc_upreg_gene_heatmap.pdf"), 
      height = 25, width = 20) 
  } else {
    pdf(paste0(plot_dir, "CNA_assoc_downreg_gene_heatmap.pdf"), 
      height = 25, width = 20) 
  }
  
    grid.newpage()
  
      # plot heatmap:
      pushViewport(viewport(x = 0.08, y = 0.04, width = 0.9, height = 0.95, 
        just = c("left", "bottom")))
        grid.draw(annotated_heatmap)
      popViewport()

      # plot heatmap legend:
      pushViewport(viewport(x = unit(1.9, "cm"), y = unit(38, "cm"), width = unit(1, "cm"), 
        height = unit(2, "cm")))
        draw(lgd, x = unit(0.1, "cm"), y = unit(0.1, "cm"), just = c("left", "bottom"))
      popViewport()

      # print subtype labels:
      pushViewport(viewport(x = unit(13, "cm"), y = unit(1.6, "cm"), width = unit(3, "cm"), 
        height = unit(1, "cm")))
        grid.text("ER+", gp=gpar(fontsize=25, fontface = "bold", col = subtype_cols[1]))
      popViewport()

      pushViewport(viewport(x = unit(28.5, "cm"), y = unit(1.6, "cm"), width = unit(3, "cm"), 
        height = unit(1, "cm")))
        grid.text("HER2+", gp=gpar(fontsize=25, fontface = "bold", col = subtype_cols[2]))
      popViewport()

      pushViewport(viewport(x = unit(41.5, "cm"), y = unit(1.6, "cm"), width = unit(3, "cm"), 
        height = unit(1, "cm")))
        grid.text("TNBC", gp=gpar(fontsize=25, fontface = "bold", col = subtype_cols[3]))
      popViewport()

  dev.off()

}


################################################################################
### 4. Load metastasis genes ###
################################################################################

# load and format msig met gene sets :
if (!file.exists(paste0(ref_dir, "msigdbr_met_genes.Rdata"))) {
  library(msigdbr, lib.loc = lib_loc)
  # fetch all cancer associated GSEA gene set categories:
  m_df = msigdbr(species = "Homo sapiens")
  cancer_cats_temp <- m_df[m_df$gs_cat == "C2" & m_df$gs_subcat == "CGP",]
  cancer_cats <- rbind(
    cancer_cats_temp, 
    m_df[m_df$gs_cat == "C6",]
  )
  
  met_indices <- grep("metasta", cancer_cats$gs_name, ignore.case = T)
  remove_indices <- grep(
    "strom|lymphatic_vessel", 
    cancer_cats$gs_name, 
    ignore.case = TRUE
  )
  remove_indices <- unique(
    c(
      remove_indices,
      grep(
        "lee|landemaine|stk11", 
        cancer_cats$gs_name, 
        ignore.case = TRUE
      )
    )
  )
  met_indices <- met_indices[!(met_indices %in% remove_indices)]
  msig_met_genes <- data.frame(
    gene = cancer_cats$human_gene_symbol[met_indices],
    set = cancer_cats$gs_name[met_indices],
    cat = cancer_cats$gs_cat[met_indices],
    subcat = cancer_cats$gs_subcat[met_indices]
  )

  saveRDS(msig_met_genes, paste0(ref_dir, "msigdbr_met_genes.Rdata"))

} else {
  msig_met_genes <- readRDS(paste0(ref_dir, "msigdbr_met_genes.Rdata"))
}

# separate up- from down-regulated met assoc genes:
msig_met_up <- msig_met_genes[grep("up", msig_met_genes$set, ignore.case = T),]
msig_met_up <- rbind(
  msig_met_up,
  msig_met_genes[msig_met_genes$set %in% up_met_sets,]
)

msig_met_down <- msig_met_genes[grep("dn", msig_met_genes$set, ignore.case = T),]
msig_met_down <- rbind(
  msig_met_down,
  msig_met_genes[msig_met_genes$set %in% down_met_sets,]
)

# remove genes in both lists:
msig_met_up <- msig_met_up[which(!(msig_met_up$gene %in% msig_met_down$gene)),]
msig_met_down <- msig_met_down[which(!(msig_met_down$gene %in% msig_met_up$gene)),]

msig_brca_met_up <- msig_met_up[
  msig_met_up$set %in% msig_brca_met_sets,
]

msig_brca_met_down <- msig_met_down[
  msig_met_down$set %in% msig_brca_met_sets,
]

brca_met_up_genes <- c(as.character(unique(msig_brca_met_up$gene)), "MUCL1", "SLPI")
brca_met_down_genes <- unique(msig_brca_met_down$gene)


################################################################################
### 6. Identify less common DE genes ###
################################################################################

less_common_DE <- unique(
  unlist(
    lapply(all_DE_data, function(x) {
      gene_vec <- as.character(
        unlist(
          lapply(all_DE_data, function(y) {
            return(
              unique(
                x$gene[
                  which(x$gene %in% y$gene)
                ]
              )
            )
          })
        )
      )
      gene_tab <- table(gene_vec)
      return(
        names(gene_tab)[
          which(gene_tab >= ceiling(length(all_DE_data)*lower_min_common_prop))
        ]
      )
    })
  )
)

################################################################################
### 7. Filter less common genes for met genes only ###
################################################################################

common_met_up <- less_common_DE[less_common_DE %in% brca_met_up_genes]
common_met_down <- less_common_DE[less_common_DE %in% brca_met_down_genes]

common_met <- list(
  up = data.frame(
    gene = common_met_up,
    direction = "up"
  ),
  down = data.frame(
    gene = common_met_down,
    direction = "down"
  )
)

lapply(common_met, function(x) {

  met_data <- lapply(
    all_gene_data, 
    filter_common_DE_genes, 
    eval(parse(text = paste0("common_met_", unique(x$direction)))), 
    unique(x$direction),
    TRUE,
    func_dir,
    ref_dir,
    Robject_dir
  )
  empty_df <- data.frame(
    gene = sort(
      x$gene
    ),
    stringsAsFactors = FALSE
  )
  plot_data <- list(
    avg_logFC = empty_df,
    p_val_adj = empty_df
  )
  
  for (l in 1:length(plot_data)) {
  
    for (s in 1:length(met_data)) {
      temp_df <- cbind(
        met_data[[s]]$gene, 
        as.data.frame(
          eval(parse(text = paste0("met_data[[s]]$", names(plot_data)[l])))
        )
      )
      colnames(temp_df) <- c("gene", names(plot_data)[l])
      plot_data[[l]] <- merge(
        plot_data[[l]],
        temp_df,
        by = "gene",
        all.x = TRUE
      )
      colnames(plot_data[[l]])[s+1] <- names(met_data)[s]
    }
  
    # remove rows with all NAs:
    plot_data[[l]] <- plot_data[[l]][
      apply(plot_data[[l]], 1, function(x) {
        !(all(is.na(x)))
      }),
    ]
  
    plot_data[[l]] <- plot_data[[l]] %>%
      column_to_rownames(loc = 1)
  
  }
  
  # remove genes DE in < min proportion required:
  plot_data$avg_logFC <- plot_data$avg_logFC[
    apply(plot_data$avg_logFC, 1, function(x) {
      length(which(!is.na(x))/length(x)) >= min_common_prop
    }),
  ]
  
  # ensure p_vals rows match to logFC:
  plot_data$p_val_adj <- plot_data$p_val_adj[rownames(plot_data$avg_logFC),]
  print(paste0(
    "Are logFC rownames identical to p_vals rownames? ", 
    identical(rownames(plot_data$avg_logFC), rownames(plot_data$p_val_adj)))
  )
  # change significant values to the ComplexHeatmap code for dots, and all
  # others to NA:
  plot_data$p_val_dots <- plot_data$p_val_adj
  plot_data$p_val_dots[plot_data$p_val_dots < 0.1] <- "*"
  plot_data$p_val_dots[plot_data$p_val_dots != "*"] <- " "
  plot_data$p_val_dots[is.na(plot_data$p_val_dots)] <- " "


################################################################################
### 8. Plot common met genes ###
################################################################################

  # define heatmap colours:
  na_less_vector <- unlist(plot_data$avg_logFC)
  na_less_vector <- na_less_vector[!is.na(na_less_vector)]

  if (unique(x$direction) == "up") {
    heatmap_cols <- colorRamp2(c(0, max(na_less_vector)), 
      c("white", "#870C0C"), space = "sRGB")
  } else {
    heatmap_cols <- colorRamp2(c(min(na_less_vector), 0), 
      c("white", "#0D3084"), space = "sRGB")
  }

  mat <- as.matrix(plot_data$avg_logFC)
  rownames(mat) <- gsub("_.*$", "", rownames(mat))
  pdots <- plot_data$p_val_dots

  if (unique(x$direction) == "up") {
    star_col <- "#1F0B99"
  } else {
    star_col <- "#870C0C"
  }

  final_heatmap <- Heatmap(
    mat, 
    na_col = "grey",
    name = paste0("hm"), 
    col = heatmap_cols,
    cluster_columns = F, cluster_rows = F,
    show_row_names = T, show_column_names = T,
    show_row_dend = T,
    show_heatmap_legend = T,
    use_raster = T, raster_device = c("png"),
    cell_fun = function(j, i, x, y, w, h, heatmap_cols) { # add text to each grid
        grid.text(pdots[i, j], x, y, gp=gpar(fontsize=30, col = star_col))
    }
  )

  ht_list <- final_heatmap
  annotated_heatmap <- grid.grabExpr(
    draw(ht_list
      #, gap = unit(6, "mm"), heatmap_legend_side = "left"
    )
  )
  dev.off()
  
  if (unique(x$direction) == "up") {
    pdf(paste0(plot_dir, "CNA_assoc_upreg_brca_met_gene_heatmap.pdf"), 
      height = 13, width = 20) 
  } else {
    pdf(paste0(plot_dir, "CNA_assoc_downreg_brca_met_gene_heatmap.pdf"), 
      height = 13, width = 20) 
  }
  
    grid.newpage()
  
      # plot heatmap:
      pushViewport(viewport(x = 0.155, y = 0.04, width = 0.82, height = 0.95, 
        just = c("left", "bottom")))
        grid.draw(annotated_heatmap)
      popViewport()
  dev.off()

})




