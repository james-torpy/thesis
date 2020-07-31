count_calls <- function(estimated_df, known_gr) {

  # add estimated CNVs to granges objects and check for overlaps:
  estimate_gr <- GRanges(
    seqnames = Rle(estimated_df$start_chr),
    ranges = IRanges(
      start = estimated_df$chr_start, 
      end = estimated_df$chr_end
    ),
    strand = Rle("*"),
    end_chr = estimated_df$end_chr,
    copy_no = estimated_df$copy_no
  )

  olaps <- findOverlaps(estimate_gr, known_gr)

  # annotate estimated CNVs with actual copy number:
  estimate_gr$actual_copy_number <- NA
  estimate_gr$actual_copy_number[queryHits(olaps)] <- 
    known_gr$copy_no[subjectHits(olaps)]

  # add chromosomal end indices to chr_data:
  chr_data$chr_ends <- chr_data$ends
  for (c in rev(2:length(chr_data$chr_ends))) {
    chr_data$chr_ends[c] <- chr_data$chr_ends[c] - chr_data$chr_ends[c-1]
  }

  # convert estimate_gr to df:
  estimate_df <- annoGR2DF(estimate_gr)
  estimate_df <- subset(estimate_df, select = -width)
  colnames(estimate_df) <- c(
    "start_chr", "start", "end", "end_chr", "copy_no", "actual_copy_number"
  )

  # resolve split CNVs:
  # if start of range == 1 and end of previous range == chromosomal end, 
  # merge:
  estimate_df$keep <- T
  for (r in 2:nrow(estimate_df)) {
    if (
      estimate_df$start[r] == 1 & 
      estimate_df$end[r-1] == chr_data$chr_ends[
        estimate_df$end_chr[r-1]
      ]
    ) {

      new_range <- estimate_df[r,]
      new_range$start_chr <- estimate_df$start_chr[r-1]
      new_range$start <- estimate_df$start[r-1]
      estimate_df <- rbind(
        estimate_df,
        new_range
      )

      estimate_df$keep[r] <- F
      estimate_df$keep[r-1] <- F

    }
  }
  estimate_df <- estimate_df[estimate_df$keep,]
  estimate_df <- estimate_df[order(estimate_df$start),]
  estimate_df <- estimate_df[order(estimate_df$start_chr),]
  estimate_df <- subset(estimate_df, select = -keep)
  # remove NA values, which mark artefacts:
  estimate_df <- estimate_df[!(is.na(estimate_df$actual_copy_number)),]

  # count correct estimates:
  estimate_df$correct <- 
    estimate_df$copy_no == estimate_df$actual_copy_number
  
  return(
    list(
      correct = length(which(estimate_df$correct)),
      total = nrow(estimate_df),
      estimated_vs_known = estimate_df
    )
  )

}