filter_CNA_assoc_DE <- function(
  DE,
  CNA_dat
  ){
  # keep only downregulated genes in losses:
  DE$CNA_assoc <- FALSE
  for (l in 1:length(CNA_dat$CNA)) {

#    print("number is:")
#    print(l)

    # determine which genes are in losses, if any:
    temp_CNA_DE <- DE[DE$gene %in% CNA_dat$CNA[[l]],]

    if (nrow(temp_CNA_DE) > 0) {
      temp_CNA_DE$CNA_assoc <- TRUE
      # add CNA data for loss genes:
      to_bind <- subset(
        CNA_dat$CNA_coords[l,], 
        select = c(call, length, start, end, chr, genomic_start, genomic_end)
      )
      if (nrow(temp_CNA_DE) > 1) {
        for (j in 1:(nrow(temp_CNA_DE)-1)) {
          to_bind <- rbind(to_bind, to_bind[1,])
        }
      }
      temp_CNA_DE <- cbind(temp_CNA_DE, to_bind)

#      print("temp_CNA_DE is:")
#      print(temp_CNA_DE)

      if (!exists("filtered_DE")) {
        filtered_DE <- temp_CNA_DE
      } else {
        filtered_DE <- rbind(filtered_DE, temp_CNA_DE)
      }
    } else {

#      print("not enough rows")

    }
#    if (exists("filtered_DE")) {
#      print("filtered_DE is:")
#      print(head(filtered_DE))
#    }

  }
  if (exists("filtered_DE")) {
    return(filtered_DE)
  }
}