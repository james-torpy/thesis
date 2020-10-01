keep_top_common_DE <- function(
    gene_entry,
    direction = "up"
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

      if (direction == "up") {
        return(
          subset(
            aggregate(
              .~gene,
              gene_entry, 
              max
            ),
            select = c(gene, p_val_adj, avg_logFC, sig)
          )
        )
      } else if (direction == "bi") {
        if (names(split_df)[k] == "down") {
          return(
            subset(
              aggregate(
                .~gene,
                gene_entry, 
                min
              ),
              select = c(gene, p_val_adj, avg_logFC, sig)
            )
          )
        } else if (names(split_df)[k] == "up") {
          return(
            subset(
              aggregate(
                .~gene,
                gene_entry, 
                max
              ),
              select = c(gene, p_val_adj, avg_logFC, sig)
            )
          )
        }
        
      }
    }
  }