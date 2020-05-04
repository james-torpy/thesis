create_extended_vector <- function(accuracy_df, column) {
  for (j in 1:nrow(accuracy_df)) {
    vector_length <- length(seq(accuracy_df$start[j], accuracy_df$end[j]))
    if (j==1) {
      result_vector <- rep(
        eval(parse(text=paste0("accuracy_df$", column, "[j]"))), vector_length
      )
    } else {
      result_vector <- c(
        result_vector,
        rep(
          eval(parse(text=paste0("accuracy_df$", column, "[j]"))), vector_length
        )
      )
    }
  }
  return(result_vector)
}
