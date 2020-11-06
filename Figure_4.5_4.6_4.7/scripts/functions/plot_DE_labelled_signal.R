plot_DE_labelled_signal <- function(
  DE_results,
  signal_data,
  CNA_data,
  chr_data
) {

  # load functions:
  subtype_signal_plot <- dget(paste0(func_dir, "subtype_signal_plot.R"))
  print_subcluster_signal <- dget(
    paste0(func_dir, "print_subcluster_signal.R")
  )

  stagger_labels <- function(df, stag_no) {
    df$stagger <- 0
    for (j in 1:length(df$index)) {
      # for each gene, check all preceeding genes:
      if (j!=1) {
        if ( (df$index[j] - df$index[j-1]) < 400 ) {
          if (df$stagger[j-1] == 0) {
            df$stagger[j] <- 1
          } else {
            df$stagger[j] <- 0
          }
        }
      }
    }
    return(df)
  }

  # label direction of DE genes:
  DE_results$direction[DE_results$avg_logFC < 0] <- "down"
  DE_results$direction[DE_results$avg_logFC > 0] <- "up"

  for (i in 1:length(signal_data$subpop_matrices)) {
    print(i)
    # index genes to annotate:
    annot_gene <- data.frame(
      gene = as.character(
        DE_results$gene
      ),
      index = match(
        as.character(
          DE_results$gene
        ),
        colnames(signal_data$subpop_matrices[[i]])
      ),
      direction = DE_results$direction
    )
    # order annot gene and stagger labels:
    annot_gene <- annot_gene[order(annot_gene$index),]
    annot_gene <- stagger_labels(annot_gene)
  
    signal_plot <- subtype_signal_plot(
      plot_df = signal_data$subpop_matrices[[i]],
      neutral_score = signal_data$neutral_value,
      CNA_indices = CNA_data[[i]],
      annotate_genes = annot_gene
    )
  
    if (i==1) {
      signal_plots <- list(signal_plot)
    } else {
      signal_plots[[i]] <- signal_plot
    }

    signal_plot_with_artefacts <- subtype_signal_plot(
      plot_df = signal_data$subpop_matrices[[i]],
      neutral_score = signal_data$neutral_value,
      CNA_indices = CNA_data[[i]],
      annotate_genes = annot_gene,
      include_artefacts = TRUE
    )
  
    if (i==1) {
      signal_plots_with_artefacts <- list(signal_plot_with_artefacts)
    } else {
      signal_plots_with_artefacts[[i]] <- signal_plot_with_artefacts
    }
  
  }
  names(signal_plots) <- names(signal_data$subpop_matrices)
  names(signal_plots_with_artefacts) <- names(signal_data$subpop_matrices)
  
  print_subcluster_signal(
    s_plots = signal_plots,
    chr_coords = chr_data,
    plot_dir,
    annotate_genes = annot_gene,
    no_genes = ncol(signal_data$subpop_matrices[[i]]),
    with_artefacts = FALSE
  )

  print_subcluster_signal(
    s_plots = signal_plots_with_artefacts,
    chr_coords = chr_data,
    plot_dir,
    annotate_genes = annot_gene,
    no_genes = ncol(signal_data$subpop_matrices[[i]]),
    with_artefacts = TRUE
  )

}
