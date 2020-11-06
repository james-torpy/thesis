find_common_DE <- function(
  DE_list,
  min_common_prop,
  max_common_no
) {

  return(
    unique(
      unlist(
        lapply(DE_list, function(x) {
          gene_vec <- as.character(
            unlist(
              lapply(DE_list, function(y) {
                return(
                  unique(
                    x$gene[
                      x$gene %in% y$gene
                    ]
                  )
                )
              })
            )
          )
          gene_tab <- table(gene_vec)

          # only take genes DE in min proportion of samples, but no more than
          # max number specified
          res <- names(gene_tab)[
            which(gene_tab >= ceiling(length(DE_list)*min_common_prop))
          ]
          if (length(res) > max_common_no) {
            res <- res[1:max_common_no]
          }
          return(
            res
          )
        })
      )
    )
  )

}