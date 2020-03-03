
### function that, for each length, introduces a gain or loss of that length ###
introduce_gaps <- function(gap_lengths, type = "gain", random_starts, epithelial_df, 
  extended_gap_record = NULL, CNV_record = NULL, modified_df = NULL, 
  log_modified_fold_change_df = NULL) {
  for (i in 1:length(feature_lengths)) {

    if (!is.null(CNV_record)) {
      total_CNV_length <- sum(CNV_record$end-CNV_record$start)
    } else {
      total_CNV_length <- 0
    }

    if ( total_CNV_length < (genome_length*0.4) ) {

      print(paste0(
        "Total CNV length still less than genome length x 0.4, adding another CNV..."
      ))
      writeLines("\n")

      for (r in 1:1000) {
        print(r)
        # choose start position at random:
        start_position <- random_starts[r]
      
        # define end position:
        end_position <- start_position+feature_lengths[i]
        gap_region <- start_position:end_position
        # check gap region is at least 1 gap CNV length span genes away from either end of genome:
        if (start_position > gap_CNV_length & end_position < (genome_length-gap_CNV_length) ) {
          # define extended gap region including gap_CNV_length number of genes flanking both sides of gap:
          extended_gap_start <- start_position-gap_CNV_length
          extended_gap_end <- end_position+gap_CNV_length
          extended_gap_region <- extended_gap_start:extended_gap_end
          # check this region does not any prior extended gap regions and if not, break loop:
          if (!is.null(extended_gap_record)) {
    
            print("extended_gap_record exists, checking region does not overlap with any prior extended gap regions...")
            
            no_break_indicator <- FALSE
    
            # check if proposed CNV overlaps within 150 genes either side of every previous CNV added:
            for (n in 1:nrow(extended_gap_record) ) {
              if (extended_gap_record$start[n] > 150 & extended_gap_record$end[n] < (genome_length - 150)) {
                prior_region <- (extended_gap_record$start[n]-150):(extended_gap_record$end[n]+150)
              } else if (extended_gap_record$start[n] <= 150 & extended_gap_record$end[n] < (genome_length - 150)) {
                prior_region <- 1:(extended_gap_record$end[n]+150)
              } else if (extended_gap_record$start[n] > 150 & extended_gap_record$end[n] >= (genome_length - 150)) {
                prior_region <- (extended_gap_record$start[n]-150):genome_length
              }
              
              if ( length(Reduce(intersect, list(prior_region, extended_gap_region))) > 0 ) {
                no_break_indicator <- TRUE
              }
            }
            if (no_break_indicator) {
              print("Overlaps found between prior extended gap region/s, trying different region...")
            } else {
              print("No overlaps found, keeping region...")
              writeLines("\n")
              break
            }
      
          } else {
            print("extended_gap_record does not exist, keeping region...")
            break
            writeLines("\n")
          }
    
        } else {
          print("Region overlaps with end of genome, trying different region...")
        }
      
      }
      
      # choose CNV multiplier depending on type:
      if (type == "gain") {
        CNV_multiplier <- 3
      } else {
        CNV_multiplier <- 0
      }
      # record region so no future extended gaps overlap:
      if (is.null(extended_gap_record)) {      
        extended_gap_record <- data.frame(
          start = extended_gap_start,
          end = extended_gap_end,
          multiplier = CNV_multiplier
        )
      } else {
        extended_gap_record <- rbind(
          extended_gap_record, 
          data.frame(
            start = extended_gap_start,
            end = extended_gap_end,
            multiplier = CNV_multiplier
          )
        )
      }
      # define gap CNVs (CNVs that flank gap):
      CNV1_start <- extended_gap_start
      CNV1_end <- extended_gap_start+99
      CNV2_start <- extended_gap_end-99
      CNV2_end <- extended_gap_end
      
      flanking_CNVs <- list(
        CNV1 = c(start=CNV1_start, end=CNV1_end),
        CNV2 = c(start=CNV2_start, end=CNV2_end)
      )
      for (j in 1:length(flanking_CNVs)) {
        # add CNV to modified_df if exists, or original epithelial_df if not:
        if (is.null(modified_df)) {
          modified_df <- epithelial_df
        }
        modified_df[flanking_CNVs[[j]]["start"]:flanking_CNVs[[j]]["end"],] <- 
          modified_df[flanking_CNVs[[j]]["start"]:flanking_CNVs[[j]]["end"],]*CNV_multiplier
  
        # generate mean original segment fold change vector with each value representing 
        # a gene:
        average_original_counts <- apply(
          epithelial_df[flanking_CNVs[[j]]["start"]:flanking_CNVs[[j]]["end"],], 1, mean
        )
        # add 0.1 to all values:
        #average_original_counts[average_original_counts == 0] <- 1e-3
        average_original_counts <- average_original_counts + 1e-4
        # determine median:
        median_average_original_counts <- median(average_original_counts)
        # divide by median to get fold change from median:
        original_fold_change <- average_original_counts/median_average_original_counts
        # check median of original fold change = 1:
        median_original_fold_change <- median(original_fold_change)
    
        # generate mean modified fold change from original median vector with each value 
        # representing a gene:
        average_modified_counts <- apply(
          modified_df[flanking_CNVs[[j]]["start"]:flanking_CNVs[[j]]["end"],], 1, mean
        )
        # add 0.1 to all values to account for zeros:
        average_modified_counts <- average_modified_counts + 1e-4
        # divide by median to get fold change from original median and add to df for plotting:
        modified_fold_change <- average_modified_counts/median_average_original_counts
        # take log of the median of modified fold change to mark on plot:
        median_modified_fold_change <- median(modified_fold_change)    
        log_median_modified_fold_change <- log10(median_modified_fold_change)
    
        # take the log10 of all fold changes:
        log_modified_fold_change <- log10(modified_fold_change)
    
        # add modified fold change values to log_original_fold_change_df:
        if (is.null(log_modified_fold_change_df)) {
          log_modified_fold_change_df <- log_original_fold_change_df
        }
        log_modified_fold_change_df$count[
          flanking_CNVs[[j]]["start"]:flanking_CNVs[[j]]["end"]
        ] <- log_modified_fold_change
        # record region so no future CNVs overlap:
        if (is.null(CNV_record)) {
          if (CNV_multiplier==0) {
            CNV_record <- data.frame(
              start = flanking_CNVs[[j]]["start"],
              end = flanking_CNVs[[j]]["end"],
              multiplier = CNV_multiplier,
              log_median_modified_FC = -3
            )
          } else {
            CNV_record <- data.frame(
              start = flanking_CNVs[[j]]["start"],
              end = flanking_CNVs[[j]]["end"],
              multiplier = CNV_multiplier,
              log_median_modified_FC = log_median_modified_fold_change
            )
          }
          
        } else {
          if (CNV_multiplier==0) {
            CNV_record <- rbind(
              CNV_record, 
              data.frame(
                start = flanking_CNVs[[j]]["start"],
                end = flanking_CNVs[[j]]["end"],
                multiplier = CNV_multiplier,
                log_median_modified_FC = -3
              )
            )
          } else {
            CNV_record <- rbind(
              CNV_record, 
              data.frame(
                start = flanking_CNVs[[j]]["start"],
                end = flanking_CNVs[[j]]["end"],
                multiplier = CNV_multiplier,
                log_median_modified_FC = log_median_modified_fold_change
              )
            )
          }
        }
      }
          
    } else {
      print(paste0(
        "Total CNV length is too close to genome length x 0.4, no CNVs added..."
      ))
      writeLines("\n")
    }
  }
  return(list(
    extended_gap_record = extended_gap_record,
    CNV_record = CNV_record,
    modified_df = modified_df,
    log_modified_fold_change_df = log_modified_fold_change_df
  ))
}