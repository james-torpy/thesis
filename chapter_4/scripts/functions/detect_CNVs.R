
detect_CNVs <- function(df, min_proportion, neutral_value) {
  
  # if at least 50% of cells have signal above or below neutral_value, 
  # label gene as loss/gain respectively:
  for (c in 1:ncol(df)) {
  
    # determine if loss:
    if ( length(which(
        round(df[,c], 6) < neutral_value
      )) >= min_proportion*length(df[,c]) ) {
  
      if (c==1) {
        call_vector <- c("loss")
      } else {
        call_vector[c] <- "loss"
      }
  
    } else if ( length(which(
        round(df[,c], 6) > neutral_value
      )) >= min_proportion*length(df[,c]) ) {
  
      if (c==1) {
        call_vector <- c("gain")
      } else {
        call_vector[c] <- "gain"
      }
  
    } else {
  
      if (c==1) {
        call_vector <- c("neutral")
      } else {
        call_vector[c] <- "neutral"
      }
  
    }
  
  }
  
  # label stretches of losses/gains as CNVs:
  # fetch indices of CNVs
  CNV_indices <- data.table(
    call = rle(call_vector)$values,
    length = rle(call_vector)$lengths
  )
  for (r in 1:nrow(CNV_indices)) {
  
    if (r==1) {
      CNV_indices$start[r] <- 1
      CNV_indices$end[r] <- CNV_indices$length[r]
    } else {
      CNV_indices$start[r] <- CNV_indices$end[r-1]+1
      CNV_indices$end[r] <- CNV_indices$start[r] + (CNV_indices$length[r]-1)
    }
  
  }
  
  # find chromosomes CNVs belong to:
  for (k in 1:length(chr_data$ends)) {
    if (k==1) {

      CNV_indices$start_chr[
        CNV_indices$start <= chr_data$ends[k]
      ] <- names(chr_data$ends)[k]

      CNV_indices$end_chr[
        CNV_indices$end <= chr_data$ends[k]
      ] <- names(chr_data$ends)[k]

    } else {

      CNV_indices$start_chr[
        CNV_indices$start <= chr_data$ends[k] & 
        CNV_indices$start > chr_data$ends[k-1]
      ] <- names(chr_data$ends)[k]

      CNV_indices$end_chr[
        CNV_indices$end <= chr_data$ends[k] & 
        CNV_indices$end > chr_data$ends[k-1]
      ] <- names(chr_data$ends)[k]

    }
  }

  # make CNVs overlapping 2 chromosomes two distinct ranges:
  CNV_indices$keep <- TRUE
  for (n in 1:nrow(CNV_indices)) { 
    if (CNV_indices$call[n] != "neutral" & 
        CNV_indices$start_chr[n] != CNV_indices$end_chr[n]) {
      # mark row for removal:
      CNV_indices$keep[n] <- FALSE
      # identify end of first chromosome:
      chr_end <- chr_data$ends[as.character(CNV_indices$start_chr[n])]
      new_rows <- data.frame(
        call = rep(CNV_indices$call[n], 2),
        length = rep(NA, 2),
        start = c(CNV_indices$start[n], chr_end+1),
        end = c(chr_end, CNV_indices$end[n]),
        start_chr = c(
          as.character(CNV_indices$start_chr[n]),
          as.character(CNV_indices$end_chr[n])
        ),
        end_chr = c(
          as.character(CNV_indices$start_chr[n]), 
          as.character(CNV_indices$end_chr[n])
        ),
        keep = rep(TRUE, 2)
      )
      new_rows$length <- new_rows$end-new_rows$start+1
      CNV_indices <- rbind(
        CNV_indices,
        new_rows
      )
    }
  }
  # remove CNV indices overlapping 2 chromosomes:
  CNV_indices <- CNV_indices[CNV_indices$keep,]
  CNV_indices <- subset(CNV_indices, select=-c(keep, end_chr))
  colnames(CNV_indices) <- gsub("start_chr", "chr", colnames(CNV_indices))

  # identify start and end genes of each CNV:
  CNV_indices$start_gene <- colnames(df)[CNV_indices$start]
  CNV_indices$end_gene <- colnames(df)[CNV_indices$end]
  
  # calculate genomic start of CNVs:
  merge_coords <- subset(gene_coords, select = c(gene_id, start))
  colnames(merge_coords) <- c("start_gene", "genomic_start")
  CNV_indices <- merge(
    CNV_indices,
    merge_coords,
    by = "start_gene"
  )
  # calculate genomic start of CNVs:
  merge_coords <- subset(gene_coords, select = c(gene_id, end))
  colnames(merge_coords) <- c("end_gene", "genomic_end")
  CNV_indices <- merge(
    CNV_indices,
    merge_coords,
    by = "end_gene"
  )
  # calculate genomic length of CNVs in mb:
  CNV_indices$genomic_length <- 
    round((CNV_indices$genomic_end-CNV_indices$genomic_start)/1000000, 1)

  # calculate midpoints for labelling:
  CNV_indices$midpoints <- 
    CNV_indices$start + floor((CNV_indices$end-CNV_indices$start)/2)
  
  return(CNV_indices[order(CNV_indices$start),])

}