prepare_common_plot_data <- function(
  DE_data_list
) {

  for (l in 1:length(DE_data_list$plot_data)) {

    for (s in 1:length(DE_data_list$merged_DE)) {
      temp_df <- cbind(
        DE_data_list$merged_DE[[s]]$gene, 
        as.data.frame(
          eval(
            parse(
              text = paste0(
                "DE_data_list$merged_DE[[s]]$", names(DE_data_list$plot_data)[l]
              )
            )
          )
        )
      )
      colnames(temp_df) <- c("gene", names(DE_data_list$plot_data)[l])
      DE_data_list$plot_data[[l]] <- merge(
        DE_data_list$plot_data[[l]],
        temp_df,
        by = "gene",
        all.x = TRUE
      )
      colnames(DE_data_list$plot_data[[l]])[s+1] <- names(DE_data_list$merged_DE)[s]
    }
  
    # remove rows with all NAs:
    DE_data_list$plot_data[[l]] <- DE_data_list$plot_data[[l]][
      apply(DE_data_list$plot_data[[l]], 1, function(x) {
        !(all(is.na(x)))
      }),
    ]
  
    DE_data_list$plot_data[[l]] <- DE_data_list$plot_data[[l]] %>%
      column_to_rownames(loc = 1)
  
  }
  
  # remove genes DE in < min proportion required:
  DE_data_list$plot_data$avg_logFC <- DE_data_list$plot_data$avg_logFC[
    apply(DE_data_list$plot_data$avg_logFC, 1, function(x) {
      length(which(!is.na(x))/length(x)) >= min_common_prop
    }),
  ]

  # ensure p_vals rows match to logFC:
  DE_data_list$plot_data$p_val_adj <- DE_data_list$plot_data$p_val_adj[
    rownames(DE_data_list$plot_data$avg_logFC),
  ]
  print(paste0(
    "Are logFC rownames identical to p_vals rownames? ", 
    identical(rownames(DE_data_list$plot_data$avg_logFC), rownames(DE_data_list$plot_data$p_val_adj)))
  )
  # change significant values to the ComplexHeatmap code for dots, and all
  # others to NA:
  DE_data_list$plot_data$p_val_dots <- DE_data_list$plot_data$p_val_adj
  DE_data_list$plot_data$p_val_dots[
    DE_data_list$plot_data$p_val_dots < 0.1
  ] <- "*"
  DE_data_list$plot_data$p_val_dots[
    DE_data_list$plot_data$p_val_dots != "*"
  ] <- " "
  DE_data_list$plot_data$p_val_dots[
    is.na(DE_data_list$plot_data$p_val_dots)
  ] <- " "

  # order matrices by least NAs, then most significant values:
  NA_sums <- apply(
    DE_data_list$plot_data$avg_logFC, 
    1, 
    function(x) sum(is.na(x))
  )

  sig_sums <- apply(
    DE_data_list$plot_data$p_val_dots, 
    1, 
    function(x) sum(x == "*")
  )

  DE_data_list$plot_data <- lapply(DE_data_list$plot_data, function(x) {
    x[order(NA_sums, -sig_sums), ]
  })

  return(DE_data_list)

}


