define_normals <- function(
  temp_heatmap, 
  temp_metadata, 
  x_outlier_multiplier,
  x_thresh_multiplier,
  y_outlier_multiplier,
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

    # define x outlier thresholds as n std devs left of first cancer medoid:
    CNA_std_dev1 <- sd(first_cancer_df$CNA_value)
    x_outlier_left <- round(
      cluster_medoids$CNA_value[cluster_medoids$type == "first_cancer"] - 
      (CNA_std_dev1*x_outlier_multiplier), 3
    )
    x_outlier_right <- round(
      cluster_medoids$CNA_value[cluster_medoids$type == "first_cancer"] + 
      (CNA_std_dev1*x_outlier_multiplier), 3
    )

    # define y outlier threshold as n std devs below of first cancer medoid:
    cor_std_dev1 <- sd(first_cancer_df$cor.est)
    y_outlier_bottom <- round(
      cluster_medoids$cor.est[cluster_medoids$type == "first_cancer"] - 
      (cor_std_dev1*y_outlier_multiplier), 3
    )
    y_outlier_top <- round(
      cluster_medoids$cor.est[cluster_medoids$type == "first_cancer"] +
      (cor_std_dev1*y_outlier_multiplier), 3
    )

    # remove outliers based on these thresholds:
    first_cancer_df$outlier <- TRUE
    first_cancer_df$outlier[
      first_cancer_df$CNA_value > x_outlier_left & 
      first_cancer_df$CNA_value < x_outlier_right &
      first_cancer_df$cor.est > y_outlier_bottom & 
      first_cancer_df$cor.est < y_outlier_top
    ] <- FALSE
    
    # define 2nd x-axis threshold as n std devs left of cluster medoid 
    # without outliers included:
    CNA_std_dev2 <- sd(first_cancer_df$CNA_value[!(first_cancer_df$outlier)])
    x_int <- round(
      cluster_medoids$CNA_value[cluster_medoids$type == "first_cancer"] - 
      (CNA_std_dev2*x_thresh_multiplier), 3
    )

    if ("normal" %in% cluster_medoids$type) {
      # define x and y-axis thresholds as halfway between first cancer 
      # and normal medoids:
      y_diff <- cluster_medoids$cor.est[
        cluster_medoids$type == "first_cancer"
      ] - cluster_medoids$cor.est[
        cluster_medoids$type == "normal"
      ]
      y_int <- cluster_medoids$cor.est[
        cluster_medoids$type == "normal"
      ] + (y_diff)/2
    } else {
      # define 1st y-axis threshold as n std devs below of cancer medoid
      # without outliers included:
      cor_std_dev2 <- sd(first_cancer_df$cor.est[!(first_cancer_df$outlier)])
      y_int <- round(
        cluster_medoids$cor.est[cluster_medoids$type == "first_cancer"] - 
        (cor_std_dev2*y_thresh_multiplier), 3
      )
    }

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

    png(paste0(plot_dir, 
        "normal_quad_plot.png"), 
        width = 450, height = 270)
        print(p)
    dev.off()
  
    p
    ggsave(paste0(plot_dir, "normal_quad_plot.pdf"),
    width = 18, height = 13, units = c("cm"))
    dev.off()

    # add outlier lines:
    p <- p + geom_vline(xintercept = x_outlier_left, colour="red")
    p <- p + geom_vline(xintercept = x_outlier_right, colour="red")
    p <- p + geom_hline(yintercept = y_outlier_bottom, colour="red")
    p <- p + geom_hline(yintercept = y_outlier_top, colour="red")
    png(paste0(plot_dir, 
        "normal_quad_plot_outlier_lines.png"), 
        width = 450, height = 270)
        print(p)
    dev.off()
  
    p
    ggsave(paste0(plot_dir, "normal_quad_plot_outlier_lines.pdf"),
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
    png(paste0(plot_dir, 
      "normal_quad_plot_clusters.png"), width = 430, height = 250)
      print(p)
    dev.off()
    
    p
    ggsave(paste0(plot_dir, "normal_quad_plot_clusters.pdf"),
    width = 18, height = 13, units = c("cm"))
    dev.off()
  
    return(temp_metadata)

  } else {

    # identify the medoid cell of the only cluster, assuming it to be 
    # cancer:
    cluster_medoid <- rownames(silhouette_result$medoids)

    # define 1st x-axis threshold as n std devs left of cluster medoid:
    CNA_std_dev1 <- sd(temp_metadata$CNA_value)
    x_outlier_left <- round(
      temp_metadata$CNA_value[temp_metadata$cell_ids %in% cluster_medoid] - 
      (CNA_std_dev1*x_thresh_multiplier), 3
    )
    x_outlier_right <- round(
      temp_metadata$CNA_value[temp_metadata$cell_ids %in% cluster_medoid] +
      (CNA_std_dev1*x_thresh_multiplier), 3
    )

    # define 1st y-axis threshold as n std devs below of cancer medoid:
    cor_std_dev1 <- sd(temp_metadata$cor.est)
    y_outlier_bottom <- round(
      temp_metadata$cor.est[temp_metadata$cell_ids %in% cluster_medoid] - 
      (cor_std_dev1*y_thresh_multiplier), 3
    )
    y_outlier_top <- round(
      temp_metadata$cor.est[temp_metadata$cell_ids %in% cluster_medoid] + 
      (cor_std_dev1*y_thresh_multiplier), 3
    )

    # remove outliers based on these thresholds:
    temp_metadata$outlier <- TRUE
    temp_metadata$outlier[
      temp_metadata$CNA_value > x_outlier_left & 
      temp_metadata$CNA_value < x_outlier_right &
      temp_metadata$cor.est > y_outlier_bottom & 
      temp_metadata$cor.est < y_outlier_top
    ] <- FALSE

    # define 2nd x-axis threshold as n std devs left of cluster medoid 
    # without outliers included:
    CNA_std_dev2 <- sd(temp_metadata$CNA_value[!(temp_metadata$outlier)])
    x_int <- round(
      temp_metadata$CNA_value[temp_metadata$cell_ids %in% cluster_medoid] - 
      (CNA_std_dev2*x_thresh_multiplier), 3
    )

    # define 1st y-axis threshold as n std devs below of cancer medoid
    # without outliers included:
    cor_std_dev2 <- sd(temp_metadata$cor.est[!(temp_metadata$outlier)])
    y_int <- round(
      temp_metadata$cor.est[temp_metadata$cell_ids %in% cluster_medoid] - 
      (cor_std_dev2*y_thresh_multiplier), 3
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
    
    # create quad plots:
    p <- ggplot(temp_metadata, 
                aes(x=CNA_value, y=cor.estimate, color=as.factor(normal_cell_call)))
    p <- p + geom_point()
    p <- p + scale_color_manual(values=c("black", "#74add1", "#b2182b"), 
                                  labels=c("Cancer", "Normal", "Unassigned"))
    p <- p + xlab("CNA value")
    p <- p + ylab("Corr. with top 5% cancer")
    p <- p + theme(legend.title = element_blank())
    p <- p + geom_hline(yintercept = y_int)
    p <- p + geom_vline(xintercept = x_int)

    png(paste0(plot_dir, 
        "normal_quad_plot.png"), 
        width = 450, height = 270)
        print(p)
    dev.off()
  
    p
    ggsave(paste0(plot_dir, "normal_quad_plot_outlier_lines.pdf"),
    width = 18, height = 13, units = c("cm"))
    dev.off()
  
    # add outlier lines:
    p <- p + geom_vline(xintercept = x_outlier_left, colour="red")
    p <- p + geom_vline(xintercept = x_outlier_right, colour="red")
    p <- p + geom_hline(yintercept = y_outlier_bottom, colour="red")
    p <- p + geom_hline(yintercept = y_outlier_top, colour="red")
    png(paste0(plot_dir, 
        "normal_quad_plot_outlier_lines.png"), 
        width = 450, height = 270)
        print(p)
    dev.off()
  
    p
    ggsave(paste0(plot_dir, "normal_quad_plot_outlier_lines.pdf"),
    width = 18, height = 13, units = c("cm"))
    dev.off()
  
    # create quad plot with clusters marked:
    temp_metadata$cluster <- "only_cluster"
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
    png(paste0(plot_dir, 
      "normal_quad_plot_clusters.png"), width = 430, height = 250)
      print(p)
    dev.off()
  
    p
    ggsave(paste0(plot_dir, "normal_quad_plot_clusters.pdf"),
    width = 18, height = 13, units = c("cm"))
    dev.off()
  
    return(temp_metadata)

  }
}