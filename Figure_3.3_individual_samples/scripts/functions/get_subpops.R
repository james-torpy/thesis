get_subpops <- function(infercnv_obj, temp_metadata) {

  subcluster_labels <- infercnv_obj@tumor_subclusters$subclusters$all_observations

  # rename each cluster:
  names(subcluster_labels) <- paste0("CNV_", 1:length(subcluster_labels))
  
  # match cell ids with subcluster labels and bind in df:
  for (s in 1:length(subcluster_labels)) {
    if (s==1) {
      subcluster_df <- data.frame(
        row.names = names(subcluster_labels[[s]]),
        subcluster_id = rep(
          names(subcluster_labels)[s], 
          length(subcluster_labels[[s]])
        ) 
      )
    } else {
      subcluster_df <- rbind(
        subcluster_df,
        data.frame(
          row.names = names(subcluster_labels[[s]]),
          subcluster_id = rep(
            names(subcluster_labels)[s], 
            length(subcluster_labels[[s]])
          ) 
        )
      )  
    }
  }
  # add to temp_metadata:
  rownames(temp_metadata) <- temp_metadata$cell_ids
  temp_metadata <- merge(temp_metadata, subcluster_df, by="row.names")
  rownames(temp_metadata) <- temp_metadata$Row.names
  temp_metadata <- subset(temp_metadata, select = -Row.names)

  return(temp_metadata)
}