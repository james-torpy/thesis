annotate_METABRIC_CNV_freq <- function(
  observation_df, 
  cnv_frequencies, 
  temp_subtype, 
  chr_ends,
  chr_lengths,
  check_correlation = FALSE,
  func_dir
) {

  library(ggplot2)
  gg_color_hue <- dget(paste0(func_dir, 
    "gg_color_hue.R"))
  
  if (temp_subtype == "all") {

    cnv <- subset(
      cnv_frequencies, 
      select = c(FEATURE, Chr, start, stop, ATTRIBUTE, Loss, Gain)
    )

  } else {

    cnv <- cnv_frequencies[
      ,c(
        1:4, 
        6, 
        which(
          colnames(cnv_frequencies) %in% c(paste0("Loss", temp_subtype), 
          paste0("Gain", temp_subtype))
        )
      )
    ]

  }
  
  # keep only protein-coding genes:
  cnv <- cnv[grep("protein_coding", cnv$ATTRIBUTE),]
  
  # append meta_cnv info onto observation_df column to put it in order:
  temp <- as.data.frame(t(observation_df[1,]))
  temp$FEATURE <- rownames(temp)
  cnv_df <- merge(temp, cnv, by="FEATURE", all.x = T)
  
  # remove duplicates from cnv_df:
  cnv_df <- cnv_df[-which(duplicated(cnv_df$FEATURE)),]
  
  # order cnv_df by observation_df order:
  m <- match(colnames(observation_df), cnv_df$FEATURE)
  cnv_df <- cnv_df[m,]
  
  # remove unwanted columns:
  cnv_df <- cnv_df[,grep("CID|IND|PID|PDX|start|stop|ATTRIBUTE", colnames(cnv_df), invert=T)]
  
  # replace NAs with 0:
  ind <- grep("Loss", colnames(cnv_df))
  cnv_df[,ind][is.na(cnv_df[,ind])] <- 0
  ind <- grep("Gain", colnames(cnv_df))
  cnv_df[,ind][is.na(cnv_df[,ind])] <- 0
  
  for (j in 1:length(chr_lengths)) {
    if (j==1) {
      cnv_df$Chr[1:chr_lengths[j]] <- j
    } else {
      cnv_df$Chr[chr_lengths[j-1]:chr_lengths[j]] <- j
    }
  }
  
  # change '23' to 'X' in chromosome column:
  cnv_df$Chr[cnv_df$Chr==23] <- "X"
  
  # make loss values negative:
  ind <- grep("Loss", colnames(cnv_df))
  cnv_df[,ind] <- -cnv_df[,ind]
  
  # create area plot of CNVs:
  # adjust plot df:
  cnv_df$BPcum <- seq(1, nrow(cnv_df))
  
  m_df <- melt(cnv_df, id.vars = c("FEATURE", "Chr", "BPcum"))
  m_df <- m_df[order(m_df$BPcum),]
  
  # prepare df for labelling:
  colnames(m_df) <- c("Gene", "Chr", "Genomic_Region", "CNV_type", "Frequency")
  m_df$CNV_type <- gsub(temp_subtype, "", m_df$CNV_type)

  # correlate mean signal with frequency values:
  if (check_correlation) {

    # determine whether there is a gain or loss at each gene:
    

  }
  
  # define colours:
  cols <- c("#F8766D", "#00BFC4")
  
  # create area plot:
  p <- ggplot(m_df, aes(x=Genomic_Region, y=Frequency, fill = CNV_type))
  p <- p + geom_bar(stat = "identity")
  p <- p + scale_fill_manual(values = cols)
  p <- p + scale_x_continuous(label = names( c(chr_ends, "") ), breaks = c(0, chr_ends),
                              expand = c(0,0))
  p <- p + theme_bw()
  p <- p + scale_y_continuous(name=paste0(temp_subtype), breaks = c(-1, 0, 1), 
                              limits = c(-1, 1))
  if (temp_subtype == "all") {

    p <- p + theme(axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      text = element_text(size=22),
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      panel.grid.major.x = element_line(
       colour = "#383838", size = 0.2
      ),
      panel.grid.minor.x = element_blank(),
      legend.position="none"
    )

  } else {

    p <- p + theme(axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      text = element_text(size=22),
      axis.title.y = element_text(size=22, angle=0, vjust = 0.5),
      panel.border = element_blank(),
      panel.grid.major.x = element_line(
       colour = "#383838", size = 0.2
      ),
      panel.grid.minor.x = element_blank(),
      legend.position="none"
    )

  }
  return(grid.grabExpr({print(p)}))
}
