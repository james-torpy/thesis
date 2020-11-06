keep_top_common_DE <- function(
    gene_entry,
    direction
  ) {

    if (nrow(gene_entry) == 1) {
    
        return(
          subset(
            gene_entry,
            select = c(gene, p_val_adj, avg_logFC, sig)
          )
        )

    } else {

      if (!(all(gene_entry$sig == FALSE))) {
        gene_entry <- gene_entry[gene_entry$sig == TRUE,]
      }

      if (nrow(gene_entry) == 1) {
    
        return(
          subset(
            gene_entry,
            select = c(gene, p_val_adj, avg_logFC, sig)
          )
        )

      } else {

        if (direction == "up") {
          
          res <- subset(
            aggregate(
              avg_logFC~gene,
              gene_entry, 
              max
            ),
            select = avg_logFC,
          )

          return(
            subset(
              merge(res, gene_entry, by = "avg_logFC"),
              select = c(gene, p_val_adj, avg_logFC, sig)
            )
          )

        } else if (direction == "down") {

          res <- subset(
            aggregate(
              avg_logFC~gene,
              gene_entry, 
              min
            ),
            select = avg_logFC,
          )

          return(
            subset(
              merge(res, gene_entry, by = "avg_logFC"),
              select = c(gene, p_val_adj, avg_logFC, sig)
            )
          )

        }
      }
    }
  }