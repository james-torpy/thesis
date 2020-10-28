plot_DE <- function(
  seurat_object,
  min_pct,
  logfc_thresh,
  return_thresh, 
  only.pos = FALSE,
  plot_dir,
  filter_sig = TRUE,
  CNA_assoc_plot = TRUE,
  raw_dir,
  table_dir,
  Robject_dir,
  ref_dir,
  func_dir,
  file_prefix,
  specific_DE = "none",
  group1 = "none",
  group2 = "none",
  specific_genes = "none",
  return_table = FALSE,
  no_genes_returned = 50
) {

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

  if (!file.exists(paste0(table_dir, file_prefix, ".txt"))) {

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

        # remove artefacts:
        CNA_data <- lapply(CNA_data, function(x) {
          return(x[grep("artefact", x$call, invert = T),])
        })

        # load CNV matrix to retrieve containing genes:
        CNV_matrix <- readRDS(paste0(Robject_dir, "5a.final_epithelial_heatmap_without_normals.Rdata")) 
        epithelial_genes <- colnames(CNV_matrix)
            
        # identify genes within CNA areas in at least one subpopulation, in the same direction as DE:
        # collate all gain and loss genes:
        CNA_genes <- lapply(CNA_data, function(x) {
          
          CNAs <- list(
            losses = x[x$call == "loss"],
            gains = x[x$call == "gain"]
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

        })

        
 
        # filter out genes not in CNAs:
        if ("cluster" %in% colnames(DE_sorted)) {

          # split by subcluster:
          split_DE <- split(DE_sorted, DE_sorted$cluster) 

           for (i in 1:length(split_DE)) {

            # split into up and downregulated genes:
            temp_DE <- split(split_DE[[i]], split_DE[[i]]$avg_logFC > 0)
            names(temp_DE) <- c("down", "up")
  
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
  
            losses <- filter_CNA_assoc_DE(temp_DE$down, spec_data$loss_data)
            amps <- filter_CNA_assoc_DE(temp_DE$up, spec_data$gain_data)
  
            if (!(is.null(losses)) & !(is.null(amps))) {
              CNA_assoc_DE_temp <- rbind(losses, amps)
            } else if (!(is.null(losses))) {
              CNA_assoc_DE_temp <- losses
            } else if (!(is.null(amps))) {
              CNA_assoc_DE_temp <- amps
            }
  
            if (i==1) {
              CNA_assoc_DE <- list(CNA_assoc_DE_temp)
            } else {
              CNA_assoc_DE[[i]] <- CNA_assoc_DE_temp
            }
  
          }
  
          CNA_assoc_DE <- do.call("rbind", CNA_assoc_DE)

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
          
        }

        

        write.table(
          CNA_assoc_DE, 
          paste0(table_dir, file_prefix, "_CNA_assoc.txt"),
          col.names = TRUE,
          row.names = FALSE,
          quote = FALSE
        )
  
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
      DE_sorted,
      seurat_object,
      file_prefix,
      Robject_dir
    )

    if (nrow(CNA_assoc_DE) > 0) {

      # plot DE results:
      CNA_assoc_DE_results <- DE_heatmap(
        CNA_assoc_DE,
        seurat_object,
        file_prefix,
        Robject_dir,
        data_type = "CNA_assoc"
      )

      if (return_table) {
        return(CNA_assoc_DE)
      }

    } else {

      if (return_table) {
        return(DE_sorted)
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
