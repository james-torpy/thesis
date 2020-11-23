plot_DE <- function(
  seurat_object,
  min_pct,
  logfc_thresh,
  return_thresh, 
  only.pos = FALSE,
  plot_dir,
  filter_sig = TRUE,
  CNA_assoc_plot = TRUE,
  min_dist = diff_CNA_min_distance,
  raw_dir,
  table_dir,
  Robject_dir,
  ref_dir,
  func_dir,
  file_prefix,
  return_table = FALSE,
  merge_subclusters = "none",
  specific_DE = "none",
  group1 = "none",
  group2 = "none",
  specific_genes = "none",
  no_genes_returned = 50
) {

  # load functions:
  library(dplyr)

  # load functions:
  filter_CNA_assoc_DE <- dget(paste0(func_dir, "filter_CNA_assoc_DE.R"))
  DE_heatmap <- dget(paste0(func_dir, "DE_heatmap.R"))

  if (specific_DE[1] != "none") {
    # change cluster labels to 'group1' or 'group2':
    idents <- as.character(Idents(seurat_object))
    for (id in specific_DE[[1]]) {
      idents[idents == id] <- group1
    }
    idents[idents != group1] <- group2
    Idents(seurat_object) <- factor(idents)
  }

  plots_done <- FALSE
  if (CNA_assoc_plot) {
    if (file.exists(paste0(table_dir, file_prefix, ".txt")) & 
      file.exists(paste0(table_dir, file_prefix, "_CNA_assoc.txt"))) {
      plots_done <- TRUE
    }
  }

  if (!(file.exists(paste0(Robject_dir, file_prefix, "_initial_DE.Rdata")))) {

    # find DE genes using MAST:
    if (specific_DE[1] == "none") {
      DE_res <- FindAllMarkers(
        only.pos = only.pos,
        object = seurat_object,
        min.pct = min_pct, 
        logfc.threshold = logfc_thresh,
        test.use = 'MAST',
        return.thresh = return_thresh
      )
    } else {
      DE_res <- FindMarkers(
        only.pos = only.pos,
        object = seurat_object,
        ident.1 = group1,
        ident.2 = group2,
        min.pct = min_pct, 
        logfc.threshold = logfc_thresh, 
        test.use = 'MAST'
      )
      DE_res$gene <- rownames(DE_res)
    }

    saveRDS(DE_res, paste0(Robject_dir, file_prefix, "_initial_DE.Rdata"))

  } else {

    DE_res <- readRDS(paste0(Robject_dir, file_prefix, "_initial_DE.Rdata"))

  }

  if (!plots_done) {

    # write table:
    if (nrow(DE_res) > 0) {

      # sort with genes with most FC at top:
      if ("cluster" %in% colnames(DE_res)) {
        DE_sorted <- arrange(DE_res, cluster, desc(avg_logFC))
      } else {
        DE_sorted <- arrange(DE_res, desc(avg_logFC))
      }
      
      # filter for significant genes if needed:
      if (filter_sig) {
        DE_sorted <- DE_sorted[DE_sorted$p_val_adj < 0.1,]
      }

      # calculate and add variances of genes for iDEA:
      DE_sorted_counts <- seurat_object@assays$RNA[DE_sorted$gene]
      DE_sorted$variance <- apply(DE_sorted_counts, 1, function(x) (sd(x))^2)

      # round values and calculate average FC:
      DE_sorted$avg_logFC <- round(DE_sorted$avg_logFC, 3)
      DE_sorted$avg_FC <- round(exp(DE_sorted$avg_logFC), 1)

      write.table(
        DE_sorted, 
        paste0(table_dir, file_prefix, ".txt"),
        col.names = TRUE,
        row.names = FALSE,
        quote = FALSE
      )
  
      if (CNA_assoc_plot) {
  
        # load CNA and gene coordinate data:
        gene_coords <- read.table(paste0(ref_dir, "infercnv_gene_order.txt"))
        colnames(gene_coords) <- c("gene", "chr", "start", "end")
        CNA_data <- readRDS(paste0(Robject_dir, "CNV_indices_and_lengths.Rdata"))

        # merge CNA data if necessary:
        if (merge_subclusters[1] != "none") {
          to_merge <- CNA_data[merge_subclusters]
          temp_merged <- do.call("rbind", to_merge)
          new_CNA_data <- append(
            CNA_data[!(names(CNA_data) %in% merge_subclusters)],
            list(temp_merged)
          )
          names(new_CNA_data) <- c(
            names(CNA_data)[!(names(CNA_data) %in% merge_subclusters)],
            paste(merge_subclusters, collapse = "_")
          )
          CNA_data <- new_CNA_data
        }

        # remove artefacts:
        CNA_data_and_artefacts <- CNA_data
        CNA_data <- lapply(CNA_data, function(x) {
          return(x[grep("artefact", x$call, invert = T),])
        })

        # load CNV matrix to retrieve containing genes:
        CNV_matrix <- readRDS(paste0(Robject_dir, "5a.final_epithelial_heatmap_without_normals.Rdata")) 
        epithelial_genes <- colnames(CNV_matrix)
            
        # identify genes within CNA areas in at least one subpopulation, in the same direction as DE:
        # collate all gain and loss genes:
        collate_CNA_genes <- function(CNA_input) {
          
          CNAs <- list(
            losses = CNA_input[CNA_input$call == "loss",],
            gains = CNA_input[CNA_input$call == "gain",]
          )

          CNA_gene_vecs <- lapply(CNAs, function(y) {

            for (i in 1:nrow(y)) {
              gene_seg <- gene_coords$gene[
                (which(gene_coords$gene == y$start_gene[i])):(which(gene_coords$gene == y$end_gene[i]))
              ]
              gene_seg <- gene_seg[gene_seg %in% epithelial_genes]
              if (i==1) {
                exp_CNA_genes <- list(one = as.character(gene_seg))
              } else {
                exp_CNA_genes[[i]] <- as.character(gene_seg)
              }
            }

            names(exp_CNA_genes) <- seq(1, length(exp_CNA_genes), 1)
            return(exp_CNA_genes) 

          })

          CNA_res <- c(
            CNA_gene_vecs,
            CNAs
          )
          names(CNA_res)[3:4] <- c("loss_coords", "gain_coords") 

          return(CNA_res)

        }

        CNA_genes <- lapply(CNA_data, collate_CNA_genes)
        CNA_genes_and_artefacts <- lapply(CNA_data_and_artefacts, collate_CNA_genes)

        # fetch lists of all gain and loss associated genes for each subpopulation:
        all_loss_genes <- lapply(CNA_genes, function(x) {
          return(unique(unlist(x$losses)))
        })
        all_gain_genes <- lapply(CNA_genes, function(x) {
          return(unique(unlist(x$gains)))
        })

        all_loss_genes_and_artefacts <- lapply(CNA_genes_and_artefacts, function(x) {
          return(unique(unlist(x$losses)))
        })
        all_gain_genes_and_artefacts <- lapply(CNA_genes_and_artefacts, function(x) {
          return(unique(unlist(x$gains)))
        })

        # filter out genes not in CNAs:
        if ("cluster" %in% colnames(DE_sorted)) {

          # split by subcluster:
          split_DE <- split(DE_sorted, as.character(DE_sorted$cluster))

          for (i in 1:length(split_DE)) {

            # split into up and downregulated genes:
            temp_DE <- split(split_DE[[i]], split_DE[[i]]$avg_logFC > 0)
            if (length(temp_DE) > 1) {
               names(temp_DE) <- c("down", "up")
            } else if (temp_DE[[1]]$avg_logFC[1] < 0) {
              names(temp_DE) <- "down"
            } else if (temp_DE[[1]]$avg_logFC[1] > 0) {
              names(temp_DE) <- "up"
            }
           
            # find corresponding CNA genes and coords:
            spec_CNA_genes <- eval(parse(text = paste0("CNA_genes$", names(split_DE)[i])))
            spec_data <- list(
              loss_data = list(
                CNA = spec_CNA_genes$losses,
                CNA_coords = spec_CNA_genes$loss_coords
              ),
              gain_data = list(
                CNA = spec_CNA_genes$gains,
                CNA_coords = spec_CNA_genes$gain_coords
              )
            )
  
            # keep only CNA-associated genes DE in correct direction:
            if (!is.null(temp_DE$down)) {
              losses <- filter_CNA_assoc_DE(temp_DE$down, spec_data$loss_data)
            }
            if (!is.null(temp_DE$up)) {
              amps <- filter_CNA_assoc_DE(temp_DE$up, spec_data$gain_data)
            }

            if (!exists("losses")) {
              losses <- NULL
            }
            if (!exists("amps")) {
              amps <- NULL
            }

            # fetch lists of all gain and loss associated genes from all other
            # subpopulations:
            other_loss_genes <- all_loss_genes_and_artefacts[
              !(names(all_loss_genes_and_artefacts) %in% names(split_DE)[i])
            ]
            other_gain_genes <- all_gain_genes_and_artefacts[
              !(names(all_gain_genes_and_artefacts) %in% names(split_DE)[i])
            ]

            # keep only genes not in CNAs from at least one other subpopulation:
            if (!(is.null(losses))) {
              losses$in_diff_CNA <- FALSE
              for (g in 1:length(losses$gene)) {
                for (o in 1:length(other_loss_genes)) {
                  if (!(losses$gene[g] %in% other_loss_genes[[o]])) {
                    losses$in_diff_CNA[g] <- TRUE
                  }
                }
              }
            }

            if (!(is.null(amps))) {
              amps$in_diff_CNA <- FALSE
              for (g in 1:length(amps$gene)) {
                for (o in 1:length(other_loss_genes)) {
                  if (!(amps$gene[g] %in% other_gain_genes[[o]])) {
                    amps$in_diff_CNA[g] <- TRUE
                  }
                }
              }
            }

            if (!(is.null(losses)) & !(is.null(amps))) {
              CNA_assoc_DE_temp <- rbind(losses, amps)
            } else if (!(is.null(losses))) {
              CNA_assoc_DE_temp <- losses
            } else if (!(is.null(amps))) {
              CNA_assoc_DE_temp <- amps
            }
            
            if (exists("CNA_assoc_DE_temp")) {
              if (!exists("CNA_assoc_DE")) {
                CNA_assoc_DE <- list(CNA_assoc_DE_temp)
              } else {
                CNA_assoc_DE <- append(CNA_assoc_DE, list(CNA_assoc_DE_temp))
              }
            }
  
          }
          
          if (exists("CNA_assoc_DE")) {
            CNA_assoc_DE <- do.call("rbind", CNA_assoc_DE)
            # remove DE genes not in areas of differential CNAs:
            CNA_assoc_DE <- CNA_assoc_DE[CNA_assoc_DE$in_diff_CNA,]
          }



        } else {

          # isolate CNA genes in group 1:
          spec_CNA_genes <- CNA_genes[names(CNA_genes) %in% unlist(specific_DE[1])]

          group_CNA_genes <- list(
            loss = unique(
              unlist(
                lapply(spec_CNA_genes, function(x) {
                  return(x$losses)
                })
              )
            ),
            gain = unique(
              unlist(
                lapply(spec_CNA_genes, function(x) {
                  return(x$gains)
                })
              )
            )
          )

          # split down and up DE genes:
          split_DE <- split(DE_sorted, DE_sorted$avg_logFC > 0)
          names(split_DE) <- c("down", "up")

          loss_DE <- split_DE$down[
            split_DE$down$gene %in% group_CNA_genes$loss,
          ]

          gain_DE <- split_DE$up[
            split_DE$up$gene %in% group_CNA_genes$gain,
          ]

          # fetch lists of all gain and loss associated genes from all other
          # subpopulations:
          group1_subpops <- eval(
            parse(
              text = paste0("specific_DE$", group1)
            )
          )
          group2_subpops <- eval(
            parse(
              text = paste0("specific_DE$", group2)
            )
          )

          group1_loss_genes <- unique(
            unlist(
              all_loss_genes[
                names(all_loss_genes) %in% group1_subpops
              ]
            )
          )
          
          group1_gain_genes <- unique(
            unlist(
              all_gain_genes[
                names(all_loss_genes) %in% group1_subpops
              ]
            )
          )

          group2_loss_genes <- unique(
            unlist(
              all_loss_genes[
                names(all_loss_genes) %in% group2_subpops
              ]
            )
          )
          group2_gain_genes <- unique(
            unlist(
              all_gain_genes[
                names(all_loss_genes) %in% group2_subpops
              ]
            )
          )

          # keep downregulated genes in areas of loss in group1, but not group2:
          loss_DE$in_diff_CNA <- FALSE
          loss_DE$in_diff_CNA[
            loss_DE$gene %in% group1_loss_genes & !(loss_DE$gene %in% group2_loss_genes)
          ] <- TRUE

          # add downregulated genes in areas of gain in group2, but not group1:
          loss_DE$in_diff_CNA[
            loss_DE$gene %in% group2_gain_genes & !(loss_DE$gene %in% group1_gain_genes)
          ] <- TRUE

          # keep upregulated genes in areas of gain in group1, but not group2:
          gain_DE$in_diff_CNA <- FALSE
          gain_DE$in_diff_CNA[
            gain_DE$gene %in% group1_gain_genes & !(gain_DE$gene %in% group2_gain_genes)
          ] <- TRUE

          # add upregulated genes in areas of loss in group2, but not group1:
          gain_DE$in_diff_CNA[
            gain_DE$gene %in% group2_loss_genes & !(gain_DE$gene %in% group1_loss_genes)
          ] <- TRUE

          CNA_assoc_DE <- rbind(
            loss_DE[loss_DE$in_diff_CNA,],
            gain_DE[gain_DE$in_diff_CNA,]
          )
 
        }

        if (exists("CNA_assoc_DE")) {
          write.table(
            CNA_assoc_DE, 
            paste0(table_dir, file_prefix, "_CNA_assoc.txt"),
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE
          )
        }
       
      }

    } else {

      # write dummy file for Snakemake:
      write.table(
        "No DE genes detected", 
        paste0(table_dir, file_prefix, ".txt"),
        col.names = TRUE,
        row.names = FALSE,
        quote = FALSE
      )

    }

  } else {
  
    DE_sorted <- read.table(
      paste0(table_dir, file_prefix, ".txt"),
      header = TRUE
    )

    if (CNA_assoc_plot) {
      CNA_assoc_DE <- read.table(
        paste0(table_dir, file_prefix, "_CNA_assoc.txt"),
        header = TRUE
      )
    }

  }

  if (nrow(DE_sorted) > 0) {

    # plot DE results:
    DE_results <- DE_heatmap(
      DE = DE_sorted,
      seurat_object = seurat_object,
      file_prefix,
      Robject_dir,
      no_genes_returned
    )

    if (exists("CNA_assoc_DE")) {
      if (nrow(CNA_assoc_DE) > 0) {

        # plot DE results:
        CNA_assoc_DE_results <- DE_heatmap(
          DE = CNA_assoc_DE,
          seurat_object = seurat_object,
          file_prefix,
          Robject_dir,
          no_genes_returned,
          data_type = "CNA_assoc",
          specific_genes = DE_results$gene
        )
  
        if (return_table) {
          return(
            list(
              all = DE_sorted,
              CNA_assoc = CNA_assoc_DE[CNA_assoc_DE$gene %in% CNA_assoc_DE_results,]
            )
          )
        }
  
      } else {
  
        if (return_table) {
          return(
            list(
              all = DE_sorted
            )
          )
        }
  
      }
    }
    
  } else {

    print("No DE genes identified")

    # create dummy file for Snakemake:
    system(
      paste0(
        "touch ", plot_dir, file_prefix, ".png"
      )
    )
  }

}
