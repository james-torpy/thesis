function (expr_matrix, hclust_method, window_size, plot = FALSE)
{
    sm_expr_data = apply(expr_matrix, 2, caTools::runmean, k = window_size)
    sm_expr_data = scale(sm_expr_data, center = TRUE, scale = FALSE)
    d = dist(t(sm_expr_data))
    h_obs = hclust(d, method = hclust_method)
    permute_col_vals <- function(df) {
        num_cells = nrow(df)
        for (i in seq(ncol(df))) {
            df[, i] = df[sample(x = seq_len(num_cells), size = num_cells,
                replace = FALSE), i]
        }
        df
    }
    flog.info(sprintf("random trees, using %g parallel threads",
             infercnv.env$GLOBAL_NUM_THREADS))
    if (infercnv.env$GLOBAL_NUM_THREADS > future::availableCores()) {
        flog.warn(sprintf("not enough cores available, setting to num avail c
ores: %g",
            future::availableCores()))
        infercnv.env$GLOBAL_NUM_THREADS <- future::availableCores()
    }
    registerDoParallel(cores = infercnv.env$GLOBAL_NUM_THREADS)
    num_rand_iters = 100
    max_rand_heights <- foreach(i = seq_len(num_rand_iters)) %dopar%
        {
            rand.tumor.expr.data = t(permute_col_vals(t(expr_matrix)))
            sm.rand.tumor.expr.data = apply(rand.tumor.expr.data,
                2, caTools::runmean, k = window_size)
            sm.rand.tumor.expr.data = scale(sm.rand.tumor.expr.data,
                center = TRUE, scale = FALSE)
            rand.dist = dist(t(sm.rand.tumor.expr.data))
            h_rand <- hclust(rand.dist, method = hclust_method)
            max_rand_height <- max(h_rand$height)
            max_rand_height
        }
    max_rand_heights <- as.numeric(max_rand_heights)
    h = h_obs$height
    max_height = max(h)
    message(sprintf("Lengths for original tree branches (h): %s",
             paste(h, sep = ",", collapse = ",")))
    message(sprintf("Max height: %g", max_height))
    message(sprintf("Lengths for max heights: %s", paste(max_rand_heights,
        sep = ",", collapse = ",")))
    e = ecdf(max_rand_heights)
    pval = 1 - e(max_height)
    message(sprintf("pval: %g", pval))
    params_list <- list(h_obs = h_obs, max_h = max_height, rand_max_height_dist = max_rand_heights,
        ecdf = e)
    if (plot) {
        .plot_tree_height_dist(params_list)
    }
    return(params_list)
}