#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

project_name <- "thesis"
subproject_name <- "Figure_2.7_identification_of_normals"

args = commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
normal_name <- args[2]
cancer_x_threshold_sd_multiplier <- as.numeric(args[3])
cancer_y_threshold_sd_multiplier <- as.numeric(args[4])
normal_x_threshold_sd_multiplier <- as.numeric(args[5])
normal_y_threshold_sd_multiplier <- as.numeric(args[6])

sample_name <- "CID4515"
normal_name <- "CID4520N"
cancer_x_threshold_sd_multiplier <- 2
cancer_y_threshold_sd_multiplier <- 1.5
normal_x_threshold_sd_multiplier <- 1
normal_y_threshold_sd_multiplier <- 1.5

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.6.0/"
library(cluster, lib.loc = lib_loc)
library(scales, lib.loc = lib_loc)
library(Seurat)
library(reshape2)
library(ggplot2)
library(cowplot)
library(fpc, lib.loc = lib_loc)
library(dplyr)
library(cluster)
#library(naturalsort, lib.loc = lib_loc)

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/",
  subproject_name, "/")
ref_dir <- paste0(project_dir, "/refs/")
func_dir <- paste0(project_dir, "/scripts/functions/")
results_dir <- paste0(project_dir, "results/")
raw_dir <- paste0(project_dir, "raw_files/")
#seurat_dir <- paste0(raw_dir, "seurat_objects/", sample_name, "/")

in_dir <- paste0(results_dir, "infercnv/", sample_name, "/")
input_dir <- paste0(results_dir, "infercnv/", sample_name, "/input_files/")
Robject_dir <- paste0(in_dir, "Rdata/")
system(paste0("mkdir -p ", Robject_dir))
plot_dir <- paste0(in_dir, "plots/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(in_dir, "tables/")
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


################################################################################
### 1. Load InferCNV output and create heatmap and metadata dfs ###
################################################################################

if (!file.exists(paste0(Robject_dir, "/1b.initial_epithelial_metadata.Rdata"))) {

  # load InferCNV output:
  print("Loading InferCNV output files...")
  infercnv_output <- as.data.frame(t(read.table(paste0(in_dir, 
    "infercnv.12_denoised.observations.txt"))))

  # load metadata df:
  metadata_df <- read.table(paste0(input_dir, "metadata.txt"), header = F,
    sep = "\t", as.is = TRUE)
  colnames(metadata_df) <- c("cell_ids", "cell_type")
  row.names(metadata_df) <- metadata_df$cell_ids

  # determine the epithelial cells and only include these in heatmap:
  print(paste0("Number of heatmap rows before non-epithelial thrown: ", 
    nrow(infercnv_output)))
  epithelial_ids <- metadata_df$cell_ids[grep("pithelial", metadata_df$cell_type)]
  epithelial_heatmap <- infercnv_output[rownames(infercnv_output) %in% epithelial_ids,]
  print(paste0("Number of heatmap rows after non-epithelial thrown: ", 
    nrow(epithelial_heatmap)))

  # create epithelial_metadata df and only include epithelial cells in epithelial_heatmap:
  print("Creating epithelial metadata df...")
  epithelial_metadata <- metadata_df[rownames(epithelial_heatmap),]

  saveRDS(epithelial_heatmap, paste0(Robject_dir, "/1a.initial_epithelial_heatmap.Rdata"))
  saveRDS(epithelial_metadata, paste0(Robject_dir, "/1b.initial_epithelial_metadata.Rdata"))

} else {

  print("Loading heatmap and metadata dfs...")
  epithelial_heatmap <- readRDS(paste0(Robject_dir, "/1a.initial_epithelial_heatmap.Rdata"))
  epithelial_metadata <- readRDS(paste0(Robject_dir, "/1b.initial_epithelial_metadata.Rdata"))

}


################################################################################
### 2. Calculate CNA and correlation values ###
################################################################################

if (!file.exists(paste0(Robject_dir, 
  "epithelial_metadata_with_cell_type_QC_CNA_and_correlation_values.Rdata"))) {
  
  # determine CNA values and add to epithelial_metadata:
  print("Determining CNA values and adding to epithelial metadata df...")
  # scale infercnv values to -1:1, square values and take the mean:
  scaled_df <- as.data.frame(rescale(as.matrix(epithelial_heatmap), c(-1,1)))
  CNA_values <- apply(scaled_df, 1, function(y) {
    return(mean(y^2))
  })
  CNA_value_df <- data.frame(
    row.names = names(CNA_values),
    CNA_value = CNA_values
  )
  epithelial_metadata <- cbind(epithelial_metadata, CNA_value_df)
  print(paste0(
    "Are epithelial_metadata rownames still in the same order as epithelial_heatmap?? ",
    identical(rownames(epithelial_heatmap), rownames(epithelial_metadata))
  ))

  # determine correlation with top 5% cancer values and add to epithelial_metadata:
  print(paste0(
    "Determining correlation with top 5% cancer values and adding to epithelial ", 
    "metadata df..."
  ))

  # plot distribution of CNV scores to determine correlation test to use:
  score_density <- density(unlist(epithelial_heatmap))
  png(
    paste0(plot_dir, "CNV_score_distribution.png"),
    height = 5,
    width = 7,
    res = 300,
    units = "in"
  )
    plot(score_density)
  dev.off()

  # determine top 5% cancer cells:
  CNA_order <- order(CNA_value_df$CNA_value, decreasing=T)
  ordered_CNA_values  <- data.frame(
    row.names = rownames(CNA_value_df)[CNA_order],
    CNA_value = CNA_value_df[CNA_order,]
  )
  top_cancer <- head(ordered_CNA_values, nrow(ordered_CNA_values)*0.05)
  
  # find average genome-wide CNV predictions across genome:
  top_cancer_CNV_average <- apply(epithelial_heatmap[rownames(top_cancer),], 2, mean)
  # find correlations of each cell's CNVs with top_GIN_CNV_average:
  cancer_correlations <- apply(epithelial_heatmap, 1, function(x) {
    if (length(unique(as.numeric(x))) == 1) {
      cor_result <- data.frame(cor.estimate="no_CNVs_recorded", 
        cor.p.value="no_CNVs_recorded")
    } else {
      cor <- cor.test(as.numeric(x), top_cancer_CNV_average, method = "pearson")
      cor_result <- data.frame(cor$estimate, cor$p.value)
    }
    return(cor_result)
  })
  correlation_df <- do.call("rbind", cancer_correlations)

  # add to epithelial_metadata:
  epithelial_metadata <- cbind(epithelial_metadata, correlation_df)
  print(paste0(
    "Are epithelial_metadata rownames still in the same order as epithelial_heatmap?? ",
    identical(rownames(epithelial_heatmap), rownames(epithelial_metadata))
  ))
  saveRDS(
    epithelial_metadata, paste0(Robject_dir, 
    "epithelial_metadata_with_cell_type_QC_CNA_and_correlation_values.Rdata")
  )

} else {

  epithelial_metadata <- readRDS(
    paste0(Robject_dir, 
    "epithelial_metadata_with_cell_type_QC_CNA_and_correlation_values.Rdata")
  )

}


################################################################################
### 3. Call normals and add to epithelial_metadata ###
################################################################################

# create and add normal cell annotations:
if (!file.exists(paste0(Robject_dir, 
  "epithelial_metadata_with_cell_type_QC_CNA_correlation_and_", 
  "normal_call_values.Rdata")) | !file.exists(paste0(plot_dir, 
  "normal_call_quad_plot_mean_of_scaled_squares.png"))) {

  print(paste0(
    "Determing normal epithelial cells and adding to epithelial ",
    "metadata df..."
  ))
  # create density plot of infercnv values:
  density_plot <- density(epithelial_metadata$CNA_value, bw="SJ")
  pdf(paste0(plot_dir, "CNA_density_plot.pdf"))
    plot(density_plot, main=NA, xlab = "Infercnv value")
  dev.off()

  # prepare df for quad plots:
  quad_df <- data.frame(
    row.names = rownames(epithelial_metadata),
    CNA_value = epithelial_metadata$CNA_value, 
    cor.estimate = epithelial_metadata$cor.estimate
  )


  ################################################################################
  ### 4. Call normals and add to epithelial_metadata ###
  ################################################################################

  # scale data:
  scaled_quad_df <- scale(quad_df) %>% as.data.frame()
  # run silhouette cluster analysis to determine clusters and thresholds:
#  increasePower  <- function(v) { v ^ (1:length(v)) }
#  test_vec <- increasePower(rep(0.1, 20))
#  for (n in test_vec[5:20]) {
#    print(n)
#    pamk_result <- pamk(scaled_quad_df, krange=1:4, alpha=n)
#    print(paste0("Number of cluster detected = ", pamk_result$nc))
#  }

  pamk_result <- pamk(scaled_quad_df, krange=1:4, alpha=1e-15)
  print(paste0("Number of cluster detected = ", pamk_result$nc))

  silhouette_result <- pam(scaled_quad_df, pamk_result$nc)
  saveRDS(silhouette_result, paste0(Robject_dir, "silhouette_result.Rdata"))

#  # order silhouette widths for plotting:
#  sil_df <- as.data.frame(silhouette_result$silinfo$widths)
#  split_sil <- split(sil_df, sil_df$cluster)
#  ordered_sils <- do.call(
#    "rbind",
#    lapply(split_sil, function(x) x[order(x$cluster),])
#  )
#  
#  ######
#  # plot silhouette widths:
#  sil_df <- sil_df[order(sil_df$cluster)]
#  p <- ggplot(, aes(x = ))
#  ######

  # if no. clusters estimated to be > 1, use cluster information to set 
  # normal vs cancer thresholds:
  if (pamk_result$nc > 1) {

    sil_values <- as.data.frame(silhouette_result$silinfo$widths)
    sil_result <- data.frame(row.names=names(silhouette_result$clustering),
      cluster=silhouette_result$clustering,
      sil_width=sil_values$sil_width)
    # add sil_result to epithelial_metadata:
    epithelial_metadata <- cbind(epithelial_metadata, sil_result)
    
    # determine normal and cancer clusters by determining the max CNA values and
    # correlation with top 5% cancer:
    cluster_split <- split(epithelial_metadata, epithelial_metadata$cluster)
    names(cluster_split) <- paste0("cluster_", names(cluster_split))
    # determine order of clusters by adding mean CNA and correlation values:
    cluster_means <- sort(
      unlist(
        lapply(cluster_split, function(x) {
          return(mean(x$CNA_value) + mean(x$cor.estimate))
        })
      )
    )
    # determine second cluster from axes as cancer cluster closest to axes:
    first_cancer_cluster <- names(cluster_means[2])
    first_cancer_df <- eval(parse(text=paste0("cluster_split$", first_cancer_cluster)))
    # make intercepts 1 std dev from mean towards axes:
    # define x-axis as 2 std devs left of mean:
    CNA_mean <- mean(first_cancer_df$CNA_value)
    CNA_std_dev <- sd(first_cancer_df$CNA_value)
    x_int <- round(CNA_mean - (cancer_x_threshold_sd_multiplier*CNA_std_dev), 3)
    # define x-axis as 2 std devs left of mean:
    cor_mean <- mean(first_cancer_df$cor.estimate)
    cor_std_dev <- sd(first_cancer_df$cor.estimate)
    y_int <- round(cor_mean - (cancer_y_threshold_sd_multiplier*cor_std_dev), 3)
    # define normal and cancer cells:
    epithelial_metadata$normal_cell_call <- "cancer"
    epithelial_metadata$normal_cell_call[
      epithelial_metadata$CNA_value < x_int & epithelial_metadata$cor.estimate < y_int
    ] <- "normal"
    epithelial_metadata$normal_cell_call[
      epithelial_metadata$CNA_value < x_int & epithelial_metadata$cor.estimate > y_int
    ] <- "unassigned"
    epithelial_metadata$normal_cell_call[
      epithelial_metadata$CNA_value > x_int & epithelial_metadata$cor.estimate < y_int
    ] <- "unassigned"
    
    # create quad plot:
    p <- ggplot(epithelial_metadata, 
                aes(x=CNA_value, y=cor.estimate, color=as.factor(normal_cell_call)))
    p <- p + geom_point()
    p <- p + scale_color_manual(values=c("black", "#74add1", "#b2182b"), 
                                  labels=c("Cancer", "Normal", "Unassigned"))
    p <- p + xlab("Infercnv level")
    p <- p + ylab("Corr. with top 5% cancer (p<0.05)")
    p <- p + theme(legend.title = element_blank())
    p <- p + geom_vline(xintercept = x_int)
    p <- p + geom_hline(yintercept = y_int)
    quad_plot <- p
    png(paste0(plot_dir, 
        "normal_call_quad_plot_mean_of_scaled_squares.png"), 
        width = 450, height = 270)
        print(quad_plot)
    dev.off()
    quad_plot
    ggsave(paste0(plot_dir, "normal_call_quad_plot_mean_of_scaled_squares.pdf"))
    dev.off()
    # create quad plot with clusters marked:
    p <- ggplot(epithelial_metadata, 
      aes(x=CNA_value, y=cor.estimate, color=as.factor(cluster)))
    p <- p + geom_point()
    p <- p + scale_color_manual(values=c("#9B59B6", "#b8e186", "#18ffff", "#ef6c00"), 
      labels=c("Cluster_1", "Cluster_2", "Cluster_3", "Cluster_4"))
    p <- p + xlab("CNA value")
    p <- p + ylab("Corr. with top 5% cancer")
    p <- p + theme(legend.title = element_blank())
    p <- p + geom_vline(xintercept = x_int)
    p <- p + geom_hline(yintercept = y_int)
    quad_plot_clusters <- p
    png(paste0(plot_dir, 
      "normal_call_quad_plot_mean_of_scaled_squares_clusters.png"), width = 430, height = 250)
      print(quad_plot_clusters)
    dev.off()
    quad_plot_clusters
    ggsave(paste0(plot_dir, "normal_call_quad_plot_mean_of_scaled_squares_clusters.pdf"))
    dev.off()

  } else {

    epithelial_metadata$cluster <- "Cluster_1"

    CNA_mean <- mean(epithelial_metadata$CNA_value)
    CNA_std_dev <- sd(epithelial_metadata$CNA_value)
    cor_mean <- mean(epithelial_metadata$cor.estimate)
    cor_std_dev <- sd(epithelial_metadata$cor.estimate)
    x_int1 <- CNA_mean - (cancer_x_threshold_sd_multiplier*CNA_std_dev)
    y_int1 <- cor_mean - (cancer_y_threshold_sd_multiplier*cor_std_dev)
    normal_outliers <- rownames(epithelial_metadata)[
      epithelial_metadata$cor.estimate < y_int1 & epithelial_metadata$CNA_value < CNA_mean
    ]
    x_int2 <- CNA_mean + (normal_x_threshold_sd_multiplier*cor_std_dev)
    y_int2 <- cor_mean + (normal_y_threshold_sd_multiplier*cor_std_dev)
    cancer_outliers <- rownames(epithelial_metadata)[
      epithelial_metadata$CNA_value > x_int2 & epithelial_metadata$cor.estimate > y_int2
    ]
    if (length(normal_outliers) >= length(cancer_outliers)) {
      x_int <- round(x_int1, 3)
      y_int <- round(y_int1, 3)
    } else if (length(normal_outliers) < length(cancer_outliers)) {
      x_int <- round(x_int2, 3)
      y_int <- round(y_int2, 3)
    }

    # define normal and cancer cells:
    epithelial_metadata$normal_cell_call <- "cancer"
    epithelial_metadata$normal_cell_call[
      epithelial_metadata$CNA_value < x_int & epithelial_metadata$cor.estimate < y_int
    ] <- "normal"
    epithelial_metadata$normal_cell_call[
      epithelial_metadata$CNA_value < x_int & epithelial_metadata$cor.estimate > y_int
    ] <- "unassigned"
    epithelial_metadata$normal_cell_call[
      epithelial_metadata$CNA_value > x_int & epithelial_metadata$cor.estimate < y_int
    ] <- "unassigned"
    # create quad plot:
    if (!("cancer" %in% unique(epithelial_metadata$normal_cell_call))) {
      p <- ggplot(epithelial_metadata, 
                  aes(x=CNA_value, y=cor.estimate, color=as.factor(normal_cell_call)))
      p <- p + geom_point()
      p <- p + scale_color_manual(values=c("#74add1", "#b2182b"), 
                                    labels=c("Normal", "Unassigned"))
      p <- p + xlab("Infercnv level")
      p <- p + ylab("Corr. with top 5% cancer (p<0.05)")
      p <- p + theme(legend.title = element_blank())
      p <- p + geom_vline(xintercept = x_int)
      p <- p + geom_hline(yintercept = y_int)
    } else if (!("normal" %in% unique(epithelial_metadata$normal_cell_call))) {
      p <- ggplot(epithelial_metadata, 
                  aes(x=CNA_value, y=cor.estimate, color=as.factor(normal_cell_call)))
      p <- p + geom_point()
      p <- p + scale_color_manual(values=c("black", "#b2182b"), 
                                    labels=c("Cancer", "Unassigned"))
      p <- p + xlab("Infercnv level")
      p <- p + ylab("Corr. with top 5% cancer (p<0.05)")
      p <- p + theme(legend.title = element_blank())
      p <- p + geom_vline(xintercept = x_int)
      p <- p + geom_hline(yintercept = y_int)
    } else if (!("unassigned" %in% unique(epithelial_metadata$normal_cell_call))) {
      p <- ggplot(epithelial_metadata, 
        aes(x=CNA_value, y=cor.estimate, color=as.factor(normal_cell_call)))
      p <- p + geom_point()
      p <- p + scale_color_manual(values=c("black", "#74add1"), 
                                    labels=c("Cancer", "Normal"))
      p <- p + xlab("Infercnv level")
      p <- p + ylab("Corr. with top 5% cancer (p<0.05)")
      p <- p + theme(legend.title = element_blank())
      p <- p + geom_vline(xintercept = x_int)
      p <- p + geom_hline(yintercept = y_int)
    } else if (length(unique(epithelial_metadata$normal_cell_call)) == 3) {
      p <- ggplot(epithelial_metadata, 
          aes(x=CNA_value, y=cor.estimate, color=as.factor(normal_cell_call)))
      p <- p + geom_point()
      p <- p + scale_color_manual(values=c("black", "#74add1", "#b2182b"), 
                                    labels=c("Cancer", "Normal", "Unassigned"))
      p <- p + xlab("Infercnv level")
      p <- p + ylab("Corr. with top 5% cancer (p<0.05)")
      p <- p + theme(legend.title = element_blank())
      p <- p + geom_vline(xintercept = x_int)
      p <- p + geom_hline(yintercept = y_int)
    }
    png(paste0(plot_dir, 
      "normal_call_quad_plot_mean_of_scaled_squares.png"), 
      width = 450, height = 270)
      print(p)
    dev.off()

    # create quad plot with clusters marked:
    p <- ggplot(epithelial_metadata, 
      aes(x=CNA_value, y=cor.estimate, color=as.factor(cluster)))
    p <- p + geom_point()
    p <- p + scale_color_manual(values=c("#9B59B6", "#b8e186", "#18ffff", "#ef6c00"), 
      labels=c("Cluster_1", "Cluster_2", "Cluster_3", "Cluster_4"))
    p <- p + xlab("CNA value")
    p <- p + ylab("Corr. with top 5% cancer")
    p <- p + theme(legend.title = element_blank())
    p <- p + geom_vline(xintercept = x_int)
    p <- p + geom_hline(yintercept = y_int)
    quad_plot_clusters <- p
    png(paste0(plot_dir, 
      "normal_call_quad_plot_mean_of_scaled_squares_clusters.png"), width = 430, height = 250)
      print(quad_plot_clusters)
    dev.off()
    quad_plot_clusters
    ggsave(paste0(plot_dir, "normal_call_quad_plot_mean_of_scaled_squares_clusters.pdf"))
    dev.off()
  }

  print(paste0(
    "Are epithelial_metadata rownames still in the same order as epithelial_heatmap?? ",
    identical(rownames(epithelial_heatmap), rownames(epithelial_metadata))
  ))

  saveRDS(epithelial_metadata, paste0(Robject_dir, 
    "epithelial_metadata_with_cell_type_QC_CNA_correlation_and_", 
    "normal_call_values.Rdata"))

} else {
  epithelial_metadata <- readRDS(paste0(Robject_dir, 
    "epithelial_metadata_with_cell_type_QC_CNA_correlation_and_", 
    "normal_call_values.Rdata"))
}


################################################################################
### 6. Add sc50 calls to metadata ###
################################################################################

if (!file.exists(paste0(Robject_dir, "epithelial_metadata_final.Rdata"))) {
  print(paste0(
    "Adding sc50 calls to epithelial metadata df..."
  ))
  # read in data and select data for current sample:
  sc50 <- read.table(paste0(ref_dir, "/sc50_calls.txt"), header = T, sep = "\t",
    stringsAsFactors = F)
  sc50$SampleID <- gsub("4290", "4290A", sc50$SampleID)
  sample_sc50 <- sc50[grep(sample_name, sc50$SampleID),]

  if (length(sample_sc50$Sept.SC50_Calls) > 0) {

    # select only cells contained in epithelial_metadata:
    print(paste0("Original sc50 call no. for ", sample_name, " = ", 
      nrow(sample_sc50)))
    rownames(sample_sc50) <- sample_sc50$SampleID
    sample_sc50 <- sample_sc50[rownames(epithelial_metadata),]
     print(paste0("sc50 call no. for ", sample_name, 
      " after selecting only cells in epithelial metadata df = ", nrow(sample_sc50)))

    # cbind to epithelial_metadata:
    epithelial_metadata <- cbind(epithelial_metadata, sample_sc50$Sept.SC50_Calls)
    colnames(epithelial_metadata)[ncol(epithelial_metadata)] <- "sc50"

  } else {

    print(paste0(
      "No sc50 calls found for sample ", sample_name, ", adding NAs for this column..."
    ))
    epithelial_metadata$sc50 <- NA
  }

  print(paste0(
    "Are epithelial_metadata rownames still in the same order as epithelial_heatmap?? ",
    identical(rownames(epithelial_heatmap), rownames(epithelial_metadata))
  ))

  saveRDS(epithelial_metadata, paste0(Robject_dir, "epithelial_metadata_final.Rdata"))
} else {
  epithelial_metadata <- readRDS(paste0(Robject_dir, "epithelial_metadata_final.Rdata"))
}


################################################################################
### 5. Order heatmap and metadata ###
################################################################################

# reorder cells starting with normals, unassigned and ending with cancer:
heatmap_metadata <- epithelial_metadata
heatmap_metadata_split <- split(heatmap_metadata, heatmap_metadata$normal_cell_call)
heatmap_metadata <- do.call("rbind",
  list(
    heatmap_metadata_split$normal,
    heatmap_metadata_split$unassigned,
    heatmap_metadata_split$cancer
  )
)
epithelial_heatmap <- epithelial_heatmap[rownames(heatmap_metadata),]
print(paste0(
  "Are heatmap_metadata rownames in the same order as epithelial_heatmap?? ",
  identical(rownames(epithelial_heatmap), rownames(heatmap_metadata))
))

saveRDS(epithelial_heatmap, paste0(Robject_dir, "epithelial_heatmap.Rdata"))
saveRDS(heatmap_metadata, paste0(Robject_dir, "heatmap_metadata.Rdata"))


################################################################################
### 7. Create non-labelled heatmap annotations ###
################################################################################

# define group annotation colours:
extra_colours <- c(brewer.pal(12, "Set3"), brewer.pal(12, "Paired"))
col_palette <- c(brewer.pal(8, "Dark2")[c(3:6,8)], brewer.pal(12, "Set3"), brewer.pal(8, "Accent"),
  "#660200", "#918940", "black", "#9ECAE1", "#3F007D", "#67000D", "#FDD0A2", "#08519C",
  "#DBB335", "#5AA050", "#807DBA", "#1A6000", "#F16913", "#FD8D3C", "#DEEBF7", 
  "#7F2704", "#DADAEB", "#FC9272", "#BCBDDC", extra_colours, "#b2182b", "#85929E", 
  "#9B59B6", "#74add1","#1b7837", "#b8e186", "#fed976","#e7298a", "#18ffff", "#ef6c00",
  "#A93226", "black","orange", "#b8bc53", "#5628ce", "#fa909c", "#8ff331","#270e26")
cluster_number <- length(unique(heatmap_metadata$cell_type))
cluster_cols <- col_palette[1:cluster_number]
names(cluster_cols) <- unique(heatmap_metadata$cell_type)
cluster_cols <- cluster_cols[levels(heatmap_metadata$cell_type)]

# create group annotations:
group_annotation_df <- subset(heatmap_metadata, select = cell_type)
group_annotation <- Heatmap(
  group_annotation_df, 
  col = cluster_cols, 
  name = "group_annotation", 
  width = unit(4, "mm"), 
  show_row_names = F, show_column_names = F,
  show_heatmap_legend = F
)

# create CNA annotation:
CNA_value_annotation <- rowAnnotation(
  correlation_annotation = anno_barplot(
    heatmap_metadata$CNA_value,
    gp = gpar(
      col = "#D95F02", 
      width = unit(4, "cm")
    ), 
    border = FALSE, 
    which = "row", 
    axis = F
  )
)
CNA_value_annotation@name <- "CNA_value"
# create QC annotations:
nUMI_annotation <- rowAnnotation(
  correlation_annotation = anno_barplot(
    heatmap_metadata$nUMI,
    gp = gpar(
      col = "#D8B72E", 
      width = unit(4, "cm")
    ), 
    border = FALSE, 
    which = "row", 
    axis = F
  )
)
nUMI_annotation@name <- "nUMI"
nGene_annotation <- rowAnnotation(
  correlation_annotation = anno_barplot(
    heatmap_metadata$nGene, name = "nGene",
    gp = gpar(
      col = "#9ECAE1", 
      width = unit(4, "cm")
    ), 
    border = FALSE, 
    which = "row", 
    axis = F
  )
)
nGene_annotation@name <- "nGene"
# create normal call annotation:
normal_call_annot_df <- subset(heatmap_metadata, select = normal_cell_call)
normal_call_annotation <- Heatmap(normal_call_annot_df, 
  col = c("unassigned" = "#E7E4D3", "normal" = "#1B7837", "cancer" = "#E7298A"), 
  name = "normal_call_annotation", width = unit(6, "mm"), 
  show_row_names = F, show_column_names = F, 
  show_heatmap_legend = F
)

# create array CNV annotation:
all_array_CNVs <- read.table(paste0(ref_dir, "all_array_CNVs.txt"))
colnames(all_array_CNVs) <- gsub("CID4499_1", "CID44991", colnames(all_array_CNVs))
if (any(colnames(all_array_CNVs) %in% sample_name)) {
  if (!file.exists(paste0(Robject_dir, "array_CNV_annotation.Rdata"))) {
    array_CNV_annotation <- create_array_CNV_annotation(epithelial_heatmap, all_array_CNVs)
    saveRDS(array_CNV_annotation, paste0(Robject_dir, "array_CNV_annotation.Rdata"))
    grid_array_heatmap <- grid.grabExpr(draw(array_CNV_annotation$array_CNV_heatmap, 
      heatmap_legend_side = "left"))
  } else {
    array_CNV_annotation <- readRDS(paste0(Robject_dir, "array_CNV_annotation.Rdata"))
    grid_array_heatmap <- grid.grabExpr(draw(array_CNV_annotation$array_CNV_heatmap, 
      heatmap_legend_side = "left"))
  }
}


################################################################################
### 8. Generate heatmap ###
################################################################################

# fetch chromosome boundary co-ordinates:
if (!file.exists(paste0(Robject_dir, "chromosome_data.Rdata"))) {
  chr_data <- fetch_chromosome_boundaries(epithelial_heatmap, ref_dir)
  saveRDS(chr_data, paste0(Robject_dir, "chromosome_data.Rdata"))
} else {
  chr_data <- readRDS(paste0(Robject_dir, "chromosome_data.Rdata"))
}
# prepare df for plotting:
plot_object <- epithelial_heatmap
colnames(plot_object) <- rep("la", ncol(plot_object))
# define heatmap colours:
na_less_vector <- unlist(plot_object)
na_less_vector <- na_less_vector[!is.na(na_less_vector)]
heatmap_cols <- colorRamp2(c(min(na_less_vector), 1, max(na_less_vector)), 
      c("#00106B", "white", "#680700"), space = "sRGB")

print("Generating final heatmap...")
# create main CNV heatmap:
final_heatmap <- Heatmap(
  plot_object, name = paste0("hm"), 
  col = heatmap_cols,
  cluster_columns = F, cluster_rows = F,
  show_row_names = F, show_column_names = T,
  column_names_gp = gpar(col = "white"),
  show_row_dend = F,
  show_heatmap_legend = F,
  use_raster = T, raster_device = c("png")
)

# determine where starting co-ordinates for heatmap are based upon longest cluster name
# (0.00604 units per character):
longest_cluster_name <- max(nchar(unique(as.character(heatmap_metadata$cell_type))))
x_coord <- longest_cluster_name*0.0037

# create non-labelled heatmap:
ht_list <- group_annotation + final_heatmap + normal_call_annotation + 
  CNA_value_annotation + nUMI_annotation + nGene_annotation
annotated_heatmap <- grid.grabExpr(
  draw(ht_list, gap = unit(6, "mm"), heatmap_legend_side = "left")
)

# plot final annotated heatmap:
pdf(paste0(plot_dir, "infercnv_plot.pdf"), height = 13, width = 18)   
  grid.newpage()
  pushViewport(viewport(x = 0.005, y = 0.065, width = 0.99, height = 0.78, 
  	just = c("left", "bottom")))
      grid.draw(annotated_heatmap)
      decorate_heatmap_body("hm", {
        for ( e in 1:length(chr_data$end_pos) ) {
        grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), 
        	gp = gpar(lwd = 1, col = "#383838"))
        }
      })
    popViewport()
    pushViewport(viewport(x=x_coord-0.05, y=0.98, width = 0.1, height = 0.1, 
      just = "bottom"))
      grid.text(as.character(sample_name), gp = gpar(fontsize = 14))
    popViewport()
    if (exists("grid_array_heatmap")) {
    	pushViewport(viewport(x = x_coord + 0.845, y = 0.86, 
     	  width = 0.8657, height = 0.13, just = c("right", "bottom")))
      grid.draw(grid_array_heatmap)
      popViewport()
    }
    
dev.off()

# write epithelial_metadata as table:
write.table(heatmap_metadata, paste0(table_dir, "epithelial_metadata.txt"), 
  sep="\t", quote=F, row.names=F, col.name=T)

print(paste0("Heatmap created, output in ", plot_dir))


################################################################################
### 9. Remove clusters without enough epithelial cells for supp plots ###
################################################################################

if (!file.exists(paste0(Robject_dir, 
    "heatmap_metadata_filtered.Rdata")) | !file.exists(paste0(Robject_dir, 
    "epithelial_heatmap_filtered.Rdata"))) {
  print("Removing non-epithelial cluster cells from heatmap and metadata dfs...")
  
  if (!exists("seurat_10X")) {
    seurat_10X <- readRDS(
      paste0(seurat_dir, "03_seurat_object_processed.Rdata")
    )
  }
  Idents(seurat_10X) <- seurat_10X@meta.data$PC_A_res.1
  cluster_ids <- as.character(unique(heatmap_metadata$cell_type))
  for (l in 1:length(cluster_ids)) {
    # fetch garnett calls for cluster id no:
    simple_id <- gsub("^.*_", "", cluster_ids[l])
    garnett_calls <- table(
      seurat_10X@meta.data$garnett_call_ext_major[
        Idents(seurat_10X) == simple_id
      ]
    )
    # if 1.5x or more non-epithelial than epithelial cells, add to remove list:
    non_epithelial <- garnett_calls[grep("pithelial", names(garnett_calls), 
      invert=T)]
    epithelial <- garnett_calls[grep("pithelial", names(garnett_calls))]
    if (any(non_epithelial > (1.5*epithelial))) {
      if (!exists("remove_clusters")) {
        remove_clusters <- c(cluster_ids[l])
      } else {
        remove_clusters <- c(remove_clusters, cluster_ids[l])
      }
    }
  }

  # remove cells attributed to clusters to be removed:
  if (exists("remove_clusters")) {
    remove_cells <- heatmap_metadata$cell_ids[
      heatmap_metadata$cell_type %in% remove_clusters
    ]
    print(paste0("No. epithelial cells to be removed = ", length(remove_cells)))
    print(paste0("No. rows in epithelial metadata df before removing ",
      "non-epithelial cluster cells = ", nrow(heatmap_metadata)))
    heatmap_metadata <- heatmap_metadata[
      !(heatmap_metadata$cell_ids %in% remove_cells),
    ]
    print(paste0("No. rows in epithelial metadata df after removing ",
    "non-epithelial cluster cells = ", nrow(heatmap_metadata)))
  
    epithelial_heatmap <- epithelial_heatmap[!(rownames(epithelial_heatmap) %in% remove_cells),]
    print(paste0(
      "Are heatmap_metadata rownames in the same order as epithelial_heatmap?? ",
      identical(rownames(epithelial_heatmap), rownames(heatmap_metadata))
    ))
  } else {
    print(paste0(
      "Are heatmap_metadata rownames in the same order as epithelial_heatmap?? ",
      identical(rownames(epithelial_heatmap), rownames(heatmap_metadata))
    ))
  }

  # adjust factor levels of heatmap_metadata$cell_type for missing cell types:
  heatmap_metadata$cell_type <- gsub("_", " ", heatmap_metadata$cell_type)
  heatmap_metadata$cell_type <- factor(
    heatmap_metadata$cell_type,
    levels = naturalsort(unique(heatmap_metadata$cell_type))
  )

  saveRDS(epithelial_heatmap,  paste0(Robject_dir, "epithelial_heatmap_filtered.Rdata"))
  saveRDS(heatmap_metadata,  paste0(Robject_dir, "heatmap_metadata_filtered.Rdata"))
} else {
  epithelial_heatmap <- readRDS(
    paste0(Robject_dir, "epithelial_heatmap_filtered.Rdata")
  )
  heatmap_metadata <- readRDS(
    paste0(Robject_dir, "heatmap_metadata_filtered.Rdata")
  )
}


################################################################################
### 10. Create non-labelled heatmap annotations ###
################################################################################

# redefine cluster and colour number:
cluster_number <- length(unique(heatmap_metadata$cell_type))
cluster_cols <- col_palette[1:cluster_number]
names(cluster_cols) <- unique(heatmap_metadata$cell_type)
cluster_cols <- cluster_cols[levels(heatmap_metadata$cell_type)]

# create group annotations:
group_annotation_labelled_df <- subset(heatmap_metadata, select = cell_type)
group_annotation_labelled <- Heatmap(
  group_annotation_labelled_df, 
  col = cluster_cols, 
  name = "group_annotation_labelled", 
  width = unit(4, "mm"), 
  show_row_names = F, show_column_names = F,
  heatmap_legend_param = list(
      title = "Cluster", title_gp = gpar(fontsize = 12, fontface = "bold"), 
      labels_gp = gpar(fontsize = 12), 
      at = as.character(levels(group_annotation_labelled_df$cell_type))
    )
)

# create CNA annotation:
CNA_value_annotation <- rowAnnotation(
  correlation_annotation = anno_barplot(
    heatmap_metadata$CNA_value,
    gp = gpar(
      col = "#D95F02", 
      width = unit(4, "cm")
    ), 
    border = FALSE, 
    which = "row", 
    axis = F
  )
)
CNA_value_annotation@name <- "CNA_value"
# create QC annotations:
nUMI_annotation <- rowAnnotation(
  correlation_annotation = anno_barplot(
    heatmap_metadata$nUMI,
    gp = gpar(
      col = "#D8B72E", 
      width = unit(4, "cm")
    ), 
    border = FALSE, 
    which = "row", 
    axis = F
  )
)
nUMI_annotation@name <- "nUMI"
nGene_annotation <- rowAnnotation(
  correlation_annotation = anno_barplot(
    heatmap_metadata$nGene, name = "nGene",
    gp = gpar(
      col = "#9ECAE1", 
      width = unit(4, "cm")
    ), 
    border = FALSE, 
    which = "row", 
    axis = F
  )
)
nGene_annotation@name <- "nGene"

# subset group annotation df to remove filtered epithelial cells:
normal_call_annot_labelled_df <- data.frame(
  row.names = rownames(normal_call_annot_df)[
    rownames(normal_call_annot_df) %in% rownames(heatmap_metadata)
  ],
  cell_type = normal_call_annot_df[rownames(heatmap_metadata),]
)
normal_call_annotation_labelled <- Heatmap(normal_call_annot_labelled_df, 
  col = c("unassigned" = "#E7E4D3", "normal" = "#1B7837", "cancer" = "#E7298A"), 
  name = "normal_call_annotation_labelled", width = unit(6, "mm"), 
  show_row_names = F, show_column_names = F, 
  heatmap_legend_param = list(title = "Normal calls", 
    title_gp = gpar(fontsize = 12, fontface = "bold"), 
    labels_gp = gpar(fontsize = 12))
)


################################################################################
### 11. Generate labelled heatmap ###
################################################################################

# prepare df for plotting:
plot_object <- epithelial_heatmap
colnames(plot_object) <- rep("la", ncol(plot_object))
# define heatmap colours:
na_less_vector <- unlist(plot_object)
na_less_vector <- na_less_vector[!is.na(na_less_vector)]
heatmap_cols <- colorRamp2(c(min(na_less_vector), 1, max(na_less_vector)), 
      c("#00106B", "white", "#680700"), space = "sRGB")

final_heatmap_labelled <- Heatmap(
  plot_object, name = paste0("hm"), 
  col = heatmap_cols,
  cluster_columns = F, cluster_rows = F,
  show_row_names = F, show_column_names = T,
  column_names_gp = gpar(col = "white"),
  show_row_dend = F,
  show_heatmap_legend = T,
  heatmap_legend_param = list(title = "Modified\nexpression", color_bar = "continuous", 
    grid_height = unit(1.5, "cm"), grid_width = unit(1.5, "cm"), legend_direction = "horizontal",
    title_gp = gpar(fontsize = 12, fontface = "bold"), labels_gp = gpar(fontsize = 6)),
  use_raster = T, raster_device = c("png")
)

# create labelled heatmap:
ht_list_labelled <- group_annotation_labelled + final_heatmap_labelled + 
  normal_call_annotation_labelled + CNA_value_annotation + nUMI_annotation + 
  nGene_annotation
annotated_heatmap_labelled <- grid.grabExpr(
  draw(ht_list_labelled, gap = unit(6, "mm"), heatmap_legend_side = "left")
)

# plot final annotated heatmap:
pdf(paste0(plot_dir, "infercnv_plot_labelled.pdf"), height = 13, width = 18)   
  grid.newpage()
  pushViewport(viewport(x = 0.005, y = 0.065, width = 0.99, height = 0.78, 
    just = c("left", "bottom")))
      grid.draw(annotated_heatmap_labelled)
      decorate_heatmap_body("hm", {
        for ( e in 1:length(chr_data$end_pos) ) {
          grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), 
            gp = gpar(lwd = 1, col = "#383838"))
          grid.text(gsub("chr", "", names(chr_data$lab_pos)[e]), chr_data$lab_pos[e], 
            unit(0, "npc") + unit(-2.1, "mm"), gp=gpar(fontsize=12))
        }
      })
    popViewport()
    pushViewport(viewport(x=x_coord-0.05, y=0.98, width = 0.1, height = 0.1, 
      just = "bottom"))
      grid.text(as.character(sample_name), gp = gpar(fontsize = 14))
    popViewport()
    if (exists("grid_array_heatmap")) {
      pushViewport(viewport(x = x_coord + 0.845, y = 0.86, 
        width = 0.801, height = 0.13, just = c("right", "bottom")))
      grid.draw(grid_array_heatmap)
      popViewport()
    }

    pushViewport(viewport(x=x_coord + 0.857, y=0.01, width = 0.1, height = 0.1, just = "bottom"))
      grid.text("type", rot=65)
    popViewport()
    pushViewport(viewport(x=x_coord + 0.882, y=0.01, width = 0.1, height = 0.1, just = "bottom"))
      grid.text("CNA", rot=65)
    popViewport()
    pushViewport(viewport(x=x_coord + 0.91, y=0.01, width = 0.1, height = 0.1, just = "bottom"))
      grid.text("nUMI", rot=65)
    popViewport()
    pushViewport(viewport(x=x_coord + 0.936, y=0.01, width = 0.1, height = 0.1, just = "bottom"))
      grid.text("nGene", rot=65)
    popViewport()
    
dev.off()

print(paste0("Heatmaps created, output in ", plot_dir))

#convert pdf to png:
system(paste0("for p in ", plot_dir, "*.pdf; do echo $p; f=$(basename $p); echo $f; ",
"new=$(echo $f | sed 's/.pdf/.png/'); echo $new; ", 
"convert -density 150 ", plot_dir, "$f -quality 90 ", plot_dir, "$new; done"))
