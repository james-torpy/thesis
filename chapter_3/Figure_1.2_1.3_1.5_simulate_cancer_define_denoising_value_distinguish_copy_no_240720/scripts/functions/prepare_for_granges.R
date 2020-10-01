prepare_for_granges <- function(indices_df) {
  for (r in 1:nrow(indices_df)) {
    if (indices_df$start_chr[r] != indices_df$end_chr[r]) {

      # add additional range for second chromosome:
      new_range <- indices_df[r,]
      new_range$start <- chr_data$ends[
        names(chr_data$ends) == new_range$start_chr
      ] + 1
      new_range$start_chr <- new_range$end_chr
      indices_df <- rbind(indices_df, new_range)

      # adjust current range:
      indices_df$end[r] <- new_range$start-1
      indices_df$end_chr[r] <- indices_df$start_chr[r]

    }
  }
  indices_df <- indices_df[order(indices_df$start),]
  indices_df <- indices_df[naturalorder(indices_df$start_chr),]

  indices_df$chr_start <- chr_positions$index[indices_df$start]
  indices_df$chr_end <- chr_positions$index[indices_df$end]
  return(indices_df)
}
