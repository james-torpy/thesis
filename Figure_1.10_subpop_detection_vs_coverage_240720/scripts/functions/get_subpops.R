get_subpops <- function(infercnv_obj, temp_metadata) {

  subcluster_labels <- infercnv_obj@tumor_subclusters$subclusters$all_observations

  # rename each cluster:
  names(subcluster_labels) <- paste0(
    "CNV_", 1:length(names(subcluster_labels))
  )

  # add subcluster information to metadata:
  for (i in 1:length(subcluster_labels)) {
    temp_metadata$subcluster_ids[
      temp_metadata$cell_ids %in% subcluster_labels[[i]]
    ]
  }

#  # convert to df and split:
#  subcluster_df <- data.frame(
#    row.names = gsub("^.*\\.", "", names(unlist(subcluster_labels))),
#    cell_ids = gsub("^.*\\.", "", names(unlist(subcluster_labels))),
#    subcluster_ids = gsub("\\.CID.*$", "", names(unlist(subcluster_labels))),
#    stringsAsFactors = F
#  )
#
#  split_subcluster_ids <- split(subclust)
#
#  subcluster_df <- subcluster_df[naturalorder(subcluster_df$subcluster_ids),]
#
#  # add to temp_metadata:
#  temp_metadata <- merge(temp_metadata, subcluster_df, by="cell_ids")
#  rownames(temp_metadata) <- temp_metadata$cell_ids
#  # reorder temp_metadata:
#  temp_metadata <- temp_metadata[subcluster_df$cell_ids,]

  return(temp_metadata)
}