# define function to plot fold change modified from original values:
plot_median_FC <- function(
  orig_df,
  mod_df,
  log_orig_FC_df,
  CNV_indices
) {
  for (i in 1:nrow(CNV_indices)) {
    # generate mean original segment fold change vector with each value representing 
    # a gene:
    if (CNV_indices$start[i] != CNV_indices$end[i]){
      average_orig <- apply(
        orig_df[CNV_indices$start[i]:CNV_indices$end[i],], 1, mean
      )
    } else {
      average_orig <- mean(orig_df[CNV_indices$start[i],])
    }
    # add 0.1 to all values:
    #average_original_counts[average_original_counts == 0] <- 1e-3
    average_orig <- average_orig + 1e-4
    # determine median:
    median_average_orig <- median(average_orig)
    # divide by median to get fold change from median:
    orig_FC <- average_orig/median_average_orig
    # check median of original fold change = 1:
    median_orig_FC <- median(orig_FC)
    # generate mean modified fold change from original median vector with each value 
    # representing a gene:
    if (CNV_indices$start[i] != CNV_indices$end[i]){
      average_modified <- apply(
        mod_df[CNV_indices$start[i]:CNV_indices$end[i],], 1, mean
      )
    } else {
      average_modified <- mean(mod_df[CNV_indices$start[i],])
    }
    # add 0.1 to all values to account for zeros:
    average_modified <- average_modified + 1e-4
    # divide by median to get fold change from original median and add to df for plotting:
    modified_FC <- average_modified/median_average_orig
    # take log of the median of modified fold change to mark on plot:
    median_modified_FC <- median(modified_FC)    
    log_median_modified_FC <- log10(median_modified_FC)
    # add to CNV_record:
    CNV_indices$log_median_modified_FC[i] <- log_median_modified_FC
    # take the log10 of all fold changes:
    log_modified_FC <- log10(modified_FC)
   
    # add modified fold change values to log_original_FC_df:
    if (!exists("log_modified_FC_df")) {
      log_modified_FC_df <- log_orig_FC_df
    }
    log_modified_FC_df$count[
      CNV_indices$start[i]:CNV_indices$end[i]
    ] <- log_modified_FC
  }
  return(
    list(
      log_mod_FC_df = log_modified_FC_df,
      CNV_record = CNV_indices
    )
  )
}  