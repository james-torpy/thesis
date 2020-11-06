
merge_subclusters <- function(
    heatmap_df,
    temp_metadata,
    cor_thresh,
    coverage_prop
  ) {

  # adjust order of heatmap:
  heatmap_df <- heatmap_df[temp_metadata$cell_ids,]
  # find mean signal in each subcluster:
  subcluster_heatmaps <- split(
    heatmap_df , 
    as.character(temp_metadata$subcluster_id)
  )
  mean_signals <- lapply(subcluster_heatmaps, function(x) {
    apply(x, 2, mean)
  })
  # correlate signal from each subcluster with all others and merge any
  # with correlation > r:
  for (s in 1:length(mean_signals)) {
    cors <- lapply(mean_signals, function(x) {
      temp_cor <- cor.test(mean_signals[[s]], x, method = "pearson")
      cor_df <- data.frame(
        R2 = round(temp_cor$estimate, 2), 
        p_val = round(temp_cor$p.value, 3)
      )
      return(cor_df)
    })
    if (!exists("all_cors")) {
      all_cors <- data.frame(
        do.call("rbind", cors)
      )
      all_cors$query <- names(mean_signals)[s]
      all_cors$subject <- names(mean_signals)
    } else {
      current_cors <- data.frame(
        do.call("rbind", cors)
      )
      current_cors$query <- names(mean_signals)[s]
      current_cors$subject <- names(mean_signals)
      all_cors <- rbind(
        all_cors,
        current_cors
      )
    }
  }

  # reduce correlation df to only those that passed threshold:
  final_cors <- all_cors[
    all_cors$R2 >= cor_thresh & all_cors$p_val < 0.05 &
    all_cors$query != all_cors$subject,
  ]

  if (nrow(final_cors) > 0) {

    final_cors$keep <- TRUE
    for (f in 1:nrow(final_cors)) {
      if (f!=1) {
        if (final_cors$query[f] == final_cors$subject[f-1] &
          final_cors$subject[f] == final_cors$query[f-1]) {
          final_cors$keep[f] <- FALSE
        }
      }
    }
    final_cors <- final_cors[final_cors$keep,]
    # reduce correlation df to only those subclusters which differ
    # enough in mean coverage:
  
    for (c in 1:nrow(final_cors)) {
      if (final_cors$query[c] %in% temp_metadata$subcluster_id & 
        final_cors$subject[c] %in% temp_metadata$subcluster_id)
      # find query and subject nUMI:
      mean_UMI <- c(
        query = round(
          mean(
            temp_metadata$nUMI[
              temp_metadata$subcluster_id == final_cors$query[c]
            ]
          ), 0
        ),
        subject = round(
          mean(
            temp_metadata$nUMI[
              temp_metadata$subcluster_id == final_cors$subject[c]
            ]
          ), 0
        )
      )
      # test whether the two mean coverages are different, if so 
      # merge subclusters:
      diff_prop <- min(mean_UMI)/max(mean_UMI)
      if (diff_prop <= coverage_prop) {
        
        temp_metadata$subcluster_id[
          temp_metadata$subcluster_id == final_cors$subject[c]
        ] <- final_cors$query[c]
      }
    }
  }
  return(temp_metadata)
}