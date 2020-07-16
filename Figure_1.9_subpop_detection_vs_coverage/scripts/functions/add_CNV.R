add_CNV <- function(
  df, 
  log_orig_FC_df,
  no_CNV,
  seed_record,
  subpops,
  genome_length,
  chromosome_ends,
  ref_dir
) {

  # load function to plot fold change modified from original values:
  plot_median_FC <- dget(paste0(func_dir, "plot_median_FC.R"))

  # choose random indices for start positions:
  set.seed(seed_record["start_position",])
  gene_no <- nrow(epithelial_df)
  random_starts <- sample(1:gene_no, 1000)

  # choose random indices for lengths:
  set.seed(seed_record["length",])
  random_lengths <- CNV_lengths[
    sample(1:length(CNV_lengths), 1000, replace=T)
  ]

   # choose random indices for multipliers:
  set.seed(seed_record["multiplier",])
  random_multipliers <- CNV_multipliers[
    sample(1:length(CNV_multipliers), CNV_no, replace = T)
  ]

  for (i in 1:CNV_no) {

    if (exists("CNV_record")) {
      total_CNV_length <- sum(CNV_record$end-CNV_record$start)
    } else {
      total_CNV_length <- 0
    }


    if ( total_CNV_length < (genome_length*0.4) ) {

      print(paste0(
        "Total CNV length still less than genome length x 0.4, adding another CNV..."
      ))
      writeLines("\n")

      for (r in 1:length(random_starts)) {
        # choose start position at random:
        start_position <- random_starts[r]
      
        # choose CNV length at random:
        CNV_length <- random_lengths[r]
        end_position <- start_position+CNV_length
        CNV_region <- start_position:end_position
      
        # check this region does not overlap end of genome or any prior CNVs and if not, break loop:
        if (end_position < nrow(epithelial_df)) {
          if (exists("CNV_record")) {
    
            print("CNV_record exists, checking region does not overlap with any prior CNV regions...")
            
            break_indicator <- TRUE
    
            for (n in 1:nrow(CNV_record) ) {
              prior_region <- CNV_record$start[n]:CNV_record$end[n]
              if ( length(Reduce(intersect, list(prior_region, CNV_region))) > 0 ) {
                break_indicator <- FALSE
              }
            }

            # if region overlaps with chromosome midpoint, choose different region:
            if (any(chromosome_ends %in% CNV_region)) {
              break_indicator <- FALSE
            }

            if (break_indicator) {

              print("No overlaps found, keeping region...")
              writeLines("\n")

              # remove rth index from random_starts and random_lengths:
              random_starts <- random_starts[-r]
              random_lengths <- random_lengths[-r]

              break

            } else {

              print("Overlaps found between prior CNV region/s, or a chromosome end, trying different region...")

            }
      
          } else if (any(chromosome_ends %in% CNV_region)) {

            print("Overlaps found between a chromosome end, trying different region...")

          } else {

            print("CNV_record does not exist, keeping region...")
            writeLines("\n")

            # remove rth index from random_starts and random_lengths:
            random_starts <- random_starts[-r]
            random_lengths <- random_lengths[-r]

            break

          }
  
        } else {
          print("Region overlaps with end of genome, trying different region...")
        }

        # remove rth index from random_starts and random_lengths:
        random_starts <- random_starts[-r]
        random_lengths <- random_lengths[-r]

      }
    
      # choose CNV multiplier at random:
      CNV_multiplier <- random_multipliers[i]
  
      # record region so no future CNVs overlap:
      if (!(exists("CNV_record"))) {
        CNV_record <- data.frame(
          start = start_position,
          end = end_position,
          multiplier = CNV_multiplier,
          log_median_modified_FC = NA
        )
      } else {
        CNV_record <- rbind(
          CNV_record, 
          data.frame(
            start = start_position,
            end = end_position,
            multiplier = CNV_multiplier,
            log_median_modified_FC = NA
          )
        )
      }
    } else {
      print(paste0(
        "Total CNV length is too close to genome length x 0.4, no CNVs added..."
      ))
      writeLines("\n")
    }
  }

  if (subpops == "subpops") {

    # create list of subpopulations, one of which has the first CNV and
    # one which doesn't
    CNV_record$differentiating <- FALSE
    CNV_records <- list(
      sub1 = CNV_record,
      sub2 = CNV_record
    )
    CNV_records$sub2 <- CNV_records$sub2[-1,]

    # record removed CNV:
    CNV_records$sub1$differentiating[1] <- TRUE

     # if number of cells is odd, remove one for even split:
    if ((ncol(epithelial_df) %% 2) == 1) {
      epithelial_df <- epithelial_df[,-1]
    }
    # split epithelial_df in two:
    epithelial_dfs <- list(
      sub1 = epithelial_df[,1:(ncol(epithelial_df)/2)],
      sub2 = epithelial_df[,((ncol(epithelial_df)/2)+1):(ncol(epithelial_df))]
    )

  } else {
    CNV_records <- list(CNV_record)
    epithelial_dfs <- list(epithelial_df)
  }

  # add CNVs to modified_df:
  for (l in 1:length(CNV_records)) {

    modified_df <- epithelial_dfs[[l]]

    for (i in 1:nrow(CNV_records[[l]])) {

      modified_df[CNV_records[[l]]$start[i]:CNV_records[[l]]$end[i],] <- 
        modified_df[
          CNV_records[[l]]$start[i]:CNV_records[[l]]$end[i],
        ]*CNV_records[[l]]$multiplier[i]

    }

    mod_data <- plot_median_FC(
      epithelial_dfs[[l]],
      modified_df,
      log_orig_FC_df,
      CNV_records[[l]]
    )

    # add non-CNV regions to CNV records:
    CNV_record_gr <- GRanges(
      seqnames = Rle("genome"),
      ranges = IRanges(start = mod_data$CNV_record$start, end = mod_data$CNV_record$end),
      strand = Rle("*"),
      multiplier = mod_data$CNV_record$multiplier,
      log_median_modified_FC = mod_data$CNV_record$log_median_modified_FC,
      differentiating = mod_data$CNV_record$differentiating
    )
    # identify gaps between marked CNV regions as non-CNV regions and fill in values accordingly:
    non_CNV_record <- gaps(CNV_record_gr)
    non_CNV_record$multiplier = rep(1, length(non_CNV_record))
    non_CNV_record$log_median_modified_FC = rep(0, length(non_CNV_record))
    CNV_record_gr <- c(CNV_record_gr, non_CNV_record)
    # convert back to data frame and order:
    final_CNV_record <- data.frame(
      start = start(ranges(CNV_record_gr)),
      end = end(ranges(CNV_record_gr)),
      multiplier = CNV_record_gr$multiplier,
      log_median_modified_FC = CNV_record_gr$log_median_modified_FC,
      differentiating = CNV_record_gr$differentiating
    )
    final_CNV_record <- final_CNV_record[order(final_CNV_record$start),]
    # fill in last non-CNV segment:
    final_CNV_record <- rbind(
      final_CNV_record,
      data.frame(
        start = final_CNV_record$end[nrow(final_CNV_record)]+1,
        end = nrow(epithelial_dfs[[l]]),
        multiplier = 1,
        log_median_modified_FC = 0,
        differentiating = FALSE
      )
    )
    
    # add chromosome information:
    final_CNV_record$start_chr <- "chr1"
    final_CNV_record$end_chr <- "chr1"
    for (k in 1:length(chromosome_ends)) {
      if (k==1) {
  
        final_CNV_record$start_chr[
          final_CNV_record$start <= unlist(chromosome_ends[k])
        ] <- names(chromosome_ends)[k]
  
        final_CNV_record$end_chr[
          final_CNV_record$end <= unlist(chromosome_ends[k])
        ] <- names(chromosome_ends)[k]
  
      } else {
  
        final_CNV_record$start_chr[
          final_CNV_record$start <= unlist(chromosome_ends[k]) & 
          final_CNV_record$start > unlist(chromosome_ends[k-1])
        ] <- names(chromosome_ends)[k]
  
        final_CNV_record$end_chr[
          final_CNV_record$end <= unlist(chromosome_ends[k]) & 
          final_CNV_record$end > unlist(chromosome_ends[k-1])
        ] <- names(chromosome_ends)[k]
  
      }
    }

    final_CNV_record$differentiating[
      is.na(final_CNV_record$differentiating)
    ] <- FALSE

    if (l==1) {
      modified_dfs <- list(
        sub1 = modified_df
      )
      final_CNV_records <- list(
        sub1 = final_CNV_record
      )
      log_mod_FC_dfs <- list(
        sub1 = mod_data$log_mod_FC_df
      )
    } else {
      modified_dfs$sub2 <- modified_df
      final_CNV_records$sub2 <- final_CNV_record
      log_mod_FC_dfs$sub2 <- mod_data$log_mod_FC_df
    }

  }

  # split dfs and records and return:
  return(list(
    modified_dfs = modified_dfs,
    log_mod_FC_dfs = log_mod_FC_dfs,
    CNV_records = final_CNV_records
  ))

}