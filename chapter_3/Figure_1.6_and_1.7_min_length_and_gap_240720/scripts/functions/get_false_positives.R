
get_false_positives <- function(
  neutral_region, average_epithelial, min_CNV_length = 21, 
  gene_annotation
) {
  # set default accuracy and CNV_calls to be changed if necessary:
  neutral_region$accuracy <- "true_negative"
  neutral_region$CNV_call <- "neutral"
  neutral_region$slice_type <- "full"
  neutral_region$slice_type[
    neutral_region$end-neutral_region$start+1 < min_CNV_length
  ] <- "partial"

  # subset neutral_region:
  neutral_region <- subset(neutral_region, select = c(start, end, 
    start_chr, end_chr, type, CNV_call, accuracy, 
    slice_type))

  # fetch the mean signal at each gene:
  mean_segment_signal <- average_epithelial[
    neutral_region$start:neutral_region$end
  ]

  print(paste0("Segment coordinates are: ", neutral_region$start, ":",
    neutral_region$end))
  print(paste0(length(mean_segment_signal), " genes in segment"))

  index_types <- c("gain", "loss")

  for (i in 1:length(index_types)) {

    # determine genes with a CNV signal:
    if (index_types[i] == "gain") {
      signal_indices <- which(mean_segment_signal > 1.03)
    } else if (index_types[i] == "loss") {
      signal_indices <- which(mean_segment_signal < 0.97)
    }

    print(paste0(length(signal_indices), " genes with ", index_types[i]))

    # if total number of positive signal genes >= min_CNV_length, determine whether the 20 
    # genes around each gene has positive signal, and if so identify it as false 
    # positive:
    if (exists("signal_indices")) {
      if (length(signal_indices) >= min_CNV_length) {
        for (m in 1:length(signal_indices)) {

          # determine min_CNV_length of 20 genes around gene:
          if ( signal_indices[m] >= floor(min_CNV_length/2) & 
            signal_indices[m]+floor(min_CNV_length/2) < signal_indices[length(signal_indices)] ) {

            ind_range <- seq(
              signal_indices[m]-floor(min_CNV_length/2), signal_indices[m]+floor(min_CNV_length/2)
            )

          } else if ( !(signal_indices[m] >= floor(min_CNV_length/2)) & 
            signal_indices[m]+floor(min_CNV_length/2) < signal_indices[length(signal_indices)] ) {

            ind_range <- 1:min_CNV_length

          } else if ( signal_indices[m] >= floor(min_CNV_length/2) & 
            !(signal_indices[m]+floor(min_CNV_length/2) < signal_indices[length(signal_indices)]) ) {

            ind_range <- 
              (signal_indices[length(signal_indices)]-(min_CNV_length-1)):signal_indices[length(signal_indices)]

          }

          # determine whether all genes in that range have CNV signal:
          range_length <- length(which(ind_range %in% signal_indices))
    
          if (range_length == min_CNV_length) {
  
            # determine false CNV range:
            temp_coords <- data.frame(
              start = gene_annotation$new_number[
                which(gene_annotation$gene %in% names(mean_segment_signal)[ind_range[1]])
              ],
              end = gene_annotation$new_number[
                which(gene_annotation$gene %in% names(mean_segment_signal)[ind_range[length(ind_range)]])
              ]
            )
  
            temp_granges <- GRanges(
              seqnames = Rle(rep("all", nrow(temp_coords))),
              ranges = IRanges(
                start = temp_coords$start,
                end = temp_coords$end
              ),
              strand = Rle(rep("*", nrow(temp_coords))),
              CNV_call = index_types[i],
              accuracy = "false_positive"
            )
            
            if (exists("false_CNV_gr")) {
              false_CNV_gr <- c(false_CNV_gr, temp_granges)
            } else {
              false_CNV_gr <- temp_granges
            }
            
          }
        }
  
        if (exists("false_CNV_gr")) {
          # reduce overlapping coordinates to one:
          false_CNV_gr <- reduce(false_CNV_gr)

         # determine false CNV range chromosomes:
         false_CNV_gr$start_chr <- gene_annotation$chromosome[
           gene_annotation$new_number %in% start(false_CNV_gr)
         ]
         false_CNV_gr$end_chr <- gene_annotation$chromosome[
           gene_annotation$new_number %in% end(false_CNV_gr)
         ]
          # fill in other metadata:
          false_CNV_gr$accuracy <- "false_positive"
          false_CNV_gr$CNV_call <- index_types[i]
          false_CNV_lengths <- end(false_CNV_gr)-start(false_CNV_gr)+1
  
          print(paste0(length(false_CNV_gr), " false ", index_types[i], " identified ", 
            "with total length = ", sum(false_CNV_lengths)))

          # create all segment gr:
          if (index_types[i] == "gain") {
            segment_gr <- GRanges(
              seqnames = Rle(rep("all", nrow(neutral_region))),
              ranges = IRanges(
                start = neutral_region$start,
                end = neutral_region$end),
              strand = Rle(rep("*", nrow(neutral_region)))
            )
          } else if (index_types[i] == "loss") {
            if (exists("new_neutral_region")) {
              segment_gr <- GRanges(
                seqnames = Rle(rep("all", nrow(new_neutral_region))),
                ranges = IRanges(
                  start = new_neutral_region$start,
                  end = new_neutral_region$end),
                strand = Rle(rep("*", nrow(new_neutral_region)))
              )
            } else {
              segment_gr <- GRanges(
                seqnames = Rle(rep("all", nrow(neutral_region))),
                ranges = IRanges(
                  start = neutral_region$start,
                  end = neutral_region$end),
                strand = Rle(rep("*", nrow(neutral_region)))
              )
            }
          }

          # add false positive ranges back to neutral segment:
          segment_gr <- c(segment_gr, false_CNV_gr)
          segment_gr <- disjoin(segment_gr)
          # determine slice chromosomes:
          segment_gr$start_chr <- gene_annotation$chromosome[
            gene_annotation$new_number %in% start(segment_gr)
          ]
          segment_gr$end_chr <- gene_annotation$chromosome[
            gene_annotation$new_number %in% end(segment_gr)
          ]
          # add other metadata:
          segment_gr$accuracy <- "true_negative"
          segment_gr$accuracy[
            start(segment_gr) %in% start(false_CNV_gr)
          ] <- "false_positive"
          segment_gr$CNV_call <- neutral_region$CNV_call
          segment_gr$CNV_call[
            start(segment_gr) %in% start(false_CNV_gr)
          ] <- index_types[i]
            
          # define ranges as full (>=range) or partial (<range):
          segment_gr$slice_type <- "full"
          segment_gr$slice_type[
            end(segment_gr)-start(segment_gr)+1 < min_CNV_length
          ] <- "partial"

          if (index_types[i] == "gain") {

            # convert to data frame:
            false_CNV_segment <- data.frame(
              start = start(ranges(false_CNV_gr)),
              end = end(ranges(false_CNV_gr)),
              type = "neutral",
              start_chr = false_CNV_gr$start_chr,
              end_chr = false_CNV_gr$end_chr,
              CNV_call = false_CNV_gr$CNV_call,
              accuracy = false_CNV_gr$accuracy,
              slice_type = "full"
            )

            # keep new neutral segment to use for losses:
            new_neutral_gr <- segment_gr[
              segment_gr$CNV_call == "neutral"
            ]
            # define chromosomes of new neutral segment:
            new_neutral_gr$start_chr <- gene_annotation$chromosome[
              gene_annotation$new_number %in% start(new_neutral_gr)
            ]
            new_neutral_gr$end_chr <- gene_annotation$chromosome[
              gene_annotation$new_number %in% end(new_neutral_gr)
            ]
            new_neutral_region <- data.frame(
              start = start(ranges(new_neutral_gr)),
              end = end(ranges(new_neutral_gr)),
              type = "neutral",
              start_chr = new_neutral_gr$start_chr,
              end_chr = new_neutral_gr$end_chr,
              CNV_call = new_neutral_gr$CNV_call,
              accuracy = new_neutral_gr$accuracy,
              slice_type = new_neutral_gr$slice_type
            )

          } else if (index_types[i] == "loss") {

            # convert to data frame:
            false_CNV_segment <- data.frame(
              start = start(ranges(segment_gr)),
              end = end(ranges(segment_gr)),
              type = "neutral",
              start_chr = segment_gr$start_chr,
              end_chr = segment_gr$end_chr,
              CNV_call = segment_gr$CNV_call,
              accuracy = segment_gr$accuracy,
              slice_type = segment_gr$slice_type
            )
          }
        }
      }
    }

    # if no loss false_CNV_segment add neutral_region/new neutral entries
    # depending on whether new_neutral_gr exists:
    if ( index_types[i] == "loss" & !exists("false_CNV_segment")) {
      if (exists("new_neutral_region")) {
        false_CNV_segment <- data.frame(
          start = start(ranges(new_neutral_gr)),
          end = end(ranges(new_neutral_gr)),
          type = "neutral",
          start_chr = new_neutral_gr$start_chr,
          end_chr = new_neutral_gr$end_chr,
          CNV_call = new_neutral_gr$CNV_call,
          accuracy = new_neutral_gr$accuracy,
          slice_type = new_neutral_gr$slice_type
        )
      } else {
        false_CNV_segment <- neutral_region
      }
    }
    
    # add false segment to list if exists:
    if (exists("false_CNV_segment")) {
      if (!exists("false_positive_df")) {
        false_positive_df <- false_CNV_segment
      } else {
        false_positive_df <- rbind(false_positive_df, false_CNV_segment)
        false_positive_df <- false_positive_df[order(false_positive_df$start),]
      }
      rm(false_CNV_segment)
    }
    if (exists("false_CNV_gr")) {
      rm(false_CNV_gr)
    }
  }
  writeLines("\n")
  return(false_positive_df)
}