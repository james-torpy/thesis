plot_densities <- function(
  CNV_data,
  data_type,
  subtype_cols,
  subtype_df,
  plot_dir,
  leg_pos = "topright",
  max_val = "none"
) {

  # change names of lists to sample ids and split by subtype:
  names(CNV_data) <- gsub("\\..*$", "", names(CNV_data))
  split_vec <- subtype_df$subtype[match(names(CNV_data), subtype_df$sample)]
  spl <- split(CNV_data, split_vec)

  # plot distributions of CNV numbers:
  count_dens <- density(CNV_data)
  ER_dens <- density(spl$ER)
  HER2_dens <- density(spl$HER2)
  TNBC_dens <- density(spl$TNBC)

  # define x and y limits:
  if (max_val == "none") {
    x_limit <- range(
      count_dens$x,
      ER_dens$x,
      HER2_dens$x,
      TNBC_dens$xz
    )
  } else {
    x_limit <- c(
      min(
        count_dens$x,
        ER_dens$x,
        HER2_dens$x,
        TNBC_dens$xz
      ),
      max_val
    )
  }

  y_limit <- range(
    count_dens$y,
    ER_dens$y,
    HER2_dens$y,
    TNBC_dens$y
  )

  temp_xlab <- gsub("_", " ", paste0("CNV_", data_type))
  if (data_type == "lengths") {
    temp_xlab <- paste0(temp_xlab, " (no. genes)")
  } else if (data_type == "genomic_lengths") {
    temp_xlab <- paste0(temp_xlab, " (mb)")
  }

  pdf(paste0(plot_dir, data_type, "_density_plot.pdf"))
    plot(
      count_dens,
      main=NA,
      xlab = temp_xlab
    )
  dev.off()

  pdf(paste0(plot_dir, data_type, "_density_plot_subtypes.pdf"))
    plot(
      count_dens,
      main=NA,
      xlab = temp_xlab ,
      xlim = x_limit,
      ylim = y_limit
    )
    lines(
      ER_dens,
      main=NA,
      xlab = temp_xlab ,
      col = subtype_cols[1],
      lty = "dashed"
    )
    lines(
      HER2_dens,
      main=NA,
      xlab = temp_xlab ,
      col = subtype_cols[2],
      lty = "dashed"
    )
    lines(
      TNBC_dens,
      main=NA,
      xlab = temp_xlab ,
      col = subtype_cols[3],
      lty = "dashed"
    )
    legend(
      leg_pos,
      legend = c("all", "ER+", "HER2", "TNBC"),
      col=c("black", subtype_cols),
      pch=1
    )
  dev.off()

}