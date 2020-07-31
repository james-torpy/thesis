define_normals <- function(
  temp_heatmap, 
  temp_metadata, 
  x_thresh_multiplier,
  y_thresh_multiplier,
  plot_dir, 
  Robject_dir
) {

  # determine CNA values and add to temp_metadata:
  print("Determining CNA values and adding to epithelial metadata df...")
  # scale infercnv values to -1:1, square values and take the mean:
  scaled_df <- as.data.frame(rescale(as.matrix(temp_heatmap), c(-1,1)))
  CNA_values <- apply(scaled_df, 1, function(y) {
    return(mean(y^2))
  })
  CNA_value_df <- data.frame(
    cell_ids = names(CNA_values),
    CNA_value = CNA_values
  )
  temp_metadata <- merge(temp_metadata, CNA_value_df, by="cell_ids")

  # determine correlation with top 5% cancer values and add to temp_metadata:
  print(paste0(
    "Determining correlation with top 5% cancer values and adding to epithelial ", 
    "metadata df..."
  ))
  # determine top 5% cancer cells:
  CNA_order <- order(CNA_value_df$CNA_value, decreasing=T)
  ordered_CNA_values  <- data.frame(
    row.names = rownames(CNA_value_df)[CNA_order],
    CNA_value = CNA_value_df[CNA_order,]
  )
  top_cancer <- head(ordered_CNA_values, nrow(ordered_CNA_values)*0.05)
  
  # find average genome-wide CNV predictions across genome:
  top_cancer_CNV_average <- apply(temp_heatmap[rownames(top_cancer),], 2, mean)
  # find correlations of each cell's CNVs with top_GIN_CNV_average:
  cancer_correlations <- apply(temp_heatmap, 1, function(x) {
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

  # add to temp_metadata:
  correlation_df$cell_ids <- rownames(correlation_df)
  temp_metadata <- merge(temp_metadata, correlation_df, by = "cell_ids")
  temp_metadata <- temp_metadata[temp_metadata$cor.estimate != "no_CNVs_recorded",]
  
  # prepare df for quad plots:
  quad_df <- data.frame(
    row.names = temp_metadata$cell_ids,
    CNA_value = as.numeric(temp_metadata$CNA_value), 
    cor.estimate = as.numeric(temp_metadata$cor.estimate)
  )

  p <- ggplot(temp_metadata, 
                  aes(x=CNA_value, y=cor.estimate))
  p <- p + geom_point()
  p <- p + xlab("CNA value")
  p <- p + ylab("Corr. with top 5% cancer")
  p <- p + theme(legend.title = element_blank())
  png(paste0(plot_dir, 
      "undefined_normal_quad_plot.png"), 
      width = 450, height = 270)
      print(p)
  dev.off()

  if (!file.exists(paste0(Robject_dir, "quad clustering_result.Rdata"))) {
    # scale data:
    scaled_quad_df <- scale(quad_df) %>% as.data.frame()
    # run silhouette cluster analysis to determine clusters and thresholds:
    pamk_result <- pamk(scaled_quad_df, krange=1:6)
    pamk_result$nc
    silhouette_result <- pam(scaled_quad_df, pamk_result$nc)

    saveRDS(
      pamk_result, 
      paste0(Robject_dir, "quad_pamk_result.Rdata")
    )
    saveRDS(
      silhouette_result, 
      paste0(Robject_dir, "quad_clustering_result.Rdata")
    )
  
  } else {
    pamk_result <- readRDS(
      paste0(Robject_dir, "quad_pamk_result.Rdata")
    )
    silhouette_result <- readRDS(
      paste0(Robject_dir, "quad_clustering_result.Rdata")
    )
  }
  
  # if no. clusters estimated to be > 1, use cluster information to set 
  # normal vs cancer thresholds:
  if (pamk_result$nc > 1) {
  
    sil_values <- as.data.frame(silhouette_result$silinfo$widths)
    sil_result <- data.frame(
      cell_ids=names(silhouette_result$clustering),
      cluster=silhouette_result$clustering,
      sil_width=sil_values$sil_width
    )
    # add sil_result to temp_metadata:
    temp_metadata <- merge(temp_metadata, sil_result, by = "cell_ids")
    rownames(temp_metadata) <- temp_metadata$cell_ids
    
    # determine normal and cancer clusters by determining the max CNA values and
    # correlation with top 5% cancer:
    cluster_split <- split(temp_metadata, temp_metadata$cluster)
    names(cluster_split) <- paste0("cluster_", names(cluster_split))

    # identify the medoid cells:
    cluster_medoids <- rownames(silhouette_result$medoids)
    
    # determine order of clusters by identifying the medoid:
    cluster_medoids <- subset(
      temp_metadata[cluster_medoids,], select = c(CNA_value, cor.estimate, cluster)
    )

    # label clusters with medoid CNV value < 0.15 and correlation value < 0.25
    # as normal:
    cluster_medoids$type <- "cancer"
    cluster_medoids$type[
      cluster_medoids$CNA_value < 0.15 & cluster_medoids$cor.estimate < 0.25
    ] <- "normal"

    # label first cancer cluster:
    cancer_medoids <- cluster_medoids[
      cluster_medoids$type == "cancer",
    ]
    first_cancer_cluster <- rownames(cancer_medoids)[
      which.min(cancer_medoids$CNA_value + cancer_medoids$cor.estimate)
    ]
    cluster_medoids[first_cancer_cluster,]$type <- "first_cancer"
    first_cancer_df <- temp_metadata[
      temp_metadata$cluster == cluster_medoids$cluster[cluster_medoids$type == "first_cancer"],
    ]

    # define x-axis threshold as n std devs left of first cancer medoid:
    CNA_std_dev <- sd(first_cancer_df$CNA_value)
    x_int <- round(
      cluster_medoids$CNA_value[cluster_medoids$type == "first_cancer"] - 
      (CNA_std_dev*x_thresh_multiplier), 3
    )

    # define y-axis threshold as n std devs below of first cancer medoid:
    cor_std_dev <- sd(first_cancer_df$cor.est)
    y_int <- round(
      cluster_medoids$cor.est[cluster_medoids$type == "first_cancer"] - 
      (cor_std_dev*y_thresh_multiplier), 3
    )

    # define normal and cancer cells:
    temp_metadata$normal_cell_call <- "cancer"
    temp_metadata$normal_cell_call[
      temp_metadata$CNA_value < x_int & temp_metadata$cor.estimate < y_int
    ] <- "normal"
    temp_metadata$normal_cell_call[
      temp_metadata$CNA_value < x_int & temp_metadata$cor.estimate > y_int
    ] <- "unassigned"
    temp_metadata$normal_cell_call[
      temp_metadata$CNA_value > x_int & temp_metadata$cor.estimate < y_int
    ] <- "unassigned"
    
    # create quad plot:
    p <- ggplot(temp_metadata, 
                aes(x=CNA_value, y=cor.estimate, color=as.factor(normal_cell_call)))
    p <- p + geom_point()
    p <- p + scale_color_manual(values=c("black", "#74add1", "#b2182b"), 
                                  labels=c("Cancer", "Normal", "Unassigned"))
    p <- p + xlab("CNA value")
    p <- p + ylab("Corr. with top 5% cancer")
    p <- p + theme(legend.title = element_blank())
    p <- p + geom_vline(xintercept = x_int)
    p <- p + geom_hline(yintercept = y_int)
    quad_plot <- p
    png(paste0(plot_dir, 
        "normal_quad_plot.png"), 
        width = 450, height = 270)
        print(quad_plot)
    dev.off()
  
    quad_plot
    ggsave(paste0(plot_dir, "normal_quad_plot.pdf"),
    width = 18, height = 13, units = c("cm"))
    dev.off()
  
    # create quad plot with clusters marked:
    p <- ggplot(temp_metadata, 
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
      "normal_quad_plot_clusters.png"), width = 430, height = 250)
      print(quad_plot_clusters)
    dev.off()
  
    quad_plot_clusters
    ggsave(paste0(plot_dir, "normal_quad_plot_clusters.pdf"),
    width = 18, height = 13, units = c("cm"))
    dev.off()
  
    return(temp_metadata)

  } else {

    CNA_mean <- mean(epithelial_metadata$CNA_value)
    CNA_std_dev <- sd(epithelial_metadata$CNA_value)
    cor_mean <- mean(epithelial_metadata$cor.estimate)
    cor_std_dev <- sd(epithelial_metadata$cor.estimate)
    x_int1 <- CNA_mean - (x_thresh_multiplier*CNA_std_dev)
    y_int1 <- cor_mean - (y_thresh_multiplier*cor_std_dev)
    normal_outliers <- rownames(epithelial_metadata)[
      epithelial_metadata$cor.estimate < y_int1
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
      p <- p + xlab("CNA value")
      p <- p + ylab("Corr. with top 5% cancer")
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
    p
    ggsave(paste0(plot_dir, "normal_call_quad_plot_mean_of_scaled_squares.pdf"),
      width = 18, height = 13, units = c("cm"))
    dev.off()
  }
}