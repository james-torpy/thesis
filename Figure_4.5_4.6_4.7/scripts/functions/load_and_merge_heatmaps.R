# function to load and combine heatmap data:
load_and_merge_heatmaps <- function(
  samples, 
  in_path,
  coverage_filter,
  subcluster_method,
  subcluster_p,
  remove_artefacts,
  ref_dir, 
  subset_data = FALSE
) {

  for (s in 1:length(samples)) {

    print(paste0("Loading ", samples[s], " heatmap df..."))

    sample_Robject_dir <- paste0(
      in_path, samples[s], "/", coverage_filter, "/", subcluster_method, 
      "/p_", subcluster_p, "/", remove_artefacts, 
      "/Rdata/"
    )
    epi_heatmap <- readRDS(
      paste0(
        sample_Robject_dir, "5a.final_epithelial_heatmap_without_normals.Rdata"
      )
    )
    print(paste0("Loading ", samples[s], " metadata..."))
    epi_meta <- readRDS(
      paste0(
        sample_Robject_dir, "5b.final_epithelial_metadata_without_normals.Rdata"
      )
    )
    epi_meta <- subset(
      epi_meta, 
      select = c(cell_ids, cell_type, nUMI, nGene, CNA_value, subcluster_id)
    )
    rownames(epi_meta) <- epi_meta$cell_ids
    # order by CNV subcluster:
    epi_meta <- epi_meta[
      naturalorder(epi_meta$subcluster_id),
    ]
    epi_meta$subcluster_id <- factor(
      epi_meta$subcluster_id,
      levels = naturalsort(unique(epi_meta$subcluster_id))
    )
    epi_meta$sample <- samples[s]
    if (s==1) {
      epi_meta$type = "Primary"
    } else {
      epi_meta$type = "Metastasis"
    }
    epi_heatmap <- epi_heatmap[epi_meta$cell_ids,]
    print(paste0(
      "Are epi_meta rownames in the same order as epi_heatmap? ",
      identical(rownames(epi_heatmap), rownames(epi_meta))
    ))
  
    if (s==1) {
      epi_heatmaps <- list(epi_heatmap)
      names(epi_heatmaps) <- samples[s]
      gene_list <- colnames(epi_heatmap)
      epi_metas <- list(epi_meta)
      names(epi_metas) <- samples[s]
    } else {
      epi_heatmaps[[s]] <- epi_heatmap
      names(epi_heatmaps)[s] <- samples[s]
      gene_list <- c(gene_list, colnames(epi_heatmap))
      epi_metas[[s]] <- epi_meta
      names(epi_metas)[s] <- samples[s]
    }
  }

  # get complete gene list from samples and order:
  print("Adding missing genes and collating heatmap and metadata dfs...")
  gene_order <- read.table(paste0(ref_dir, "infercnv_gene_order.txt"))
  gene_list <- unique(gene_list)
  gene_list <- as.character(gene_order$V1[gene_order$V1 %in% gene_list])

  # keep only genes in both samples:
  for (k in 1:length(epi_heatmaps)) {
    if (k==1) {
      gene_present <- data.frame(
        row.names = gene_list,
        present = gene_list %in% colnames(epi_heatmaps[[k]])
      )
      colnames(gene_present)[k] <- names(epi_heatmaps)[k]
    } else {
      gene_present$present <- gene_list %in% colnames(epi_heatmaps[[k]])
      colnames(gene_present)[k] <- names(epi_heatmaps)[k]
    }
  }
  gene_freqs <- apply(gene_present, 1, function(x) length(which(x)))
  keep_genes <- names(gene_freqs[gene_freqs >= 1])
  
  # add all missing genes as columns to all heatmap dfs:
  complete_epi_heatmaps <- lapply(epi_heatmaps, function(x) {
    
    missing_genes <- keep_genes[!(keep_genes %in% colnames(x))]
    missing_genes_df <- data.frame(matrix(NA, nrow = nrow(x), 
      ncol = length(missing_genes)))
    colnames(missing_genes_df) <- missing_genes
  
    complete_df <- cbind(x, missing_genes_df)
    m <- match(keep_genes, colnames(complete_df))
    complete_df <- complete_df[,m]
  
    return(complete_df)
  
  })
  
  # collate all dfs and check rows line up:
  group_epi_heatmap <- do.call(rbind, complete_epi_heatmaps)
  rownames(group_epi_heatmap) <- gsub("^.*\\.C", "C", rownames(group_epi_heatmap))
  print("Are all rows present in group heatmap df?")
  identical(nrow(group_epi_heatmap), sum(unlist(lapply(complete_epi_heatmaps, nrow))))
  print("Are all genes present in group heatmap df?")
  identical(colnames(group_epi_heatmap), keep_genes)
  
  # subset group_epi_heatmap if needed:
  if (subset_data) {
    group_epi_heatmap <- group_epi_heatmap[,1:300]
  }
  
  # collate heatmap metadata:
  group_epi_meta <- do.call("rbind", epi_metas)
  rownames(group_epi_meta) <- gsub(
    "^.*\\.", "", rownames(group_epi_meta)
  )

  print(paste0(
    "Are epi_meta rownames still in the same order as epi_heatmap?? ",
    identical(rownames(group_epi_heatmap), rownames(group_epi_meta))
  ))

  # select only group_epi_meta rows in group_epi_heatmap and order:
  group_epi_meta <- group_epi_meta[rownames(group_epi_heatmap),]
  m <- match(rownames(group_epi_heatmap), rownames(group_epi_meta))
  group_epi_meta <- group_epi_meta[m,]
  
  print(paste0(
    "Are epi_meta rownames still in the same order as epi_heatmap?? ",
    identical(rownames(group_epi_heatmap), rownames(group_epi_meta))
  ))
  
  return(
    list(
      heatmap = group_epi_heatmap, 
      metadata = group_epi_meta
    )
  )

}