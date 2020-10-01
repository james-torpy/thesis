function (object, fit_noGS = FALSE, init_beta = NULL, init_tau = c(-2,
    0.5), min_degene = 5, em_iter = 15, mcmc_iter = 1000, fit.tol = 1e-05,
    modelVariant = F, verbose = TRUE, ...)
{
    num_core <- object@num_core
    if (num_core > 4) {
        if (num_core > detectCores()) {
            warning("iDEA:: the number of cores you're setting is larger than detected cores!")
            num_core = detectCores()
        }
    }
    if (.Platform$OS.type == "windows") {
        num_core <- 1
    }
    num_gene <- object@num_gene
    num_annot <- length(object@annotation)
    cat(paste("## ===== iDEA INPUT SUMMARY ==== ##\n"))
    cat(paste("## number of annotations: ", num_annot, "\n"))
    cat(paste("## number of genes: ", num_gene, "\n"))
    cat(paste("## number of cores: ", num_core, "\n"))
    if (is.null(init_beta)) {
        init_beta <- object@summary[, 1]
    }
    res_idea <- NULL
    if (fit_noGS) {
        if (verbose) {
            cat(paste("## fitting the model without gene set ... \n"))
        }
        Annot <- as.matrix(data.frame(rep(1, num_gene)))
        if (!modelVariant) {
            t1 <- system.time(model1 <- try(res <- EMMCMCStepSummary(object@summary[,
                1], object@summary[, 2], as.matrix(Annot), init_beta,
                init_tau[1], em_iter, mcmc_iter, min_degene)))
            if (class(model1) != "try-error") {
                rownames(res$pip) <- object@gene_id
                res$converged <- TRUE
                res$ctime <- t1[3]
            }
            else {
                res <- NULL
            }
            object@noGS <- res
        }
        else {
            t1 <- system.time(model1 <- try(res <- EMMCMCStepSummaryVariant(object@summary[,
                1], object@summary[, 2], as.matrix(Annot), init_beta,
                init_tau[1], em_iter, mcmc_iter, min_degene)))
            if (class(model1) != "try-error") {
                rownames(res$pip) <- object@gene_id
                res$converged <- TRUE
                res$ctime <- t1[3]
            }
            else {
                res <- NULL
            }
            object@noGS <- res
        }
    }
    else {
        if (!modelVariant) {
            if (verbose) {
                cat(paste("## fitting the model with gene sets information... \n"))
            }
            # problem lies here:
            res_idea <- pbmclapply(1:num_annot, FUN = function(x) {
                Annot <- rep(0, object@num_gene)
                Annot[object@annotation[[x]]] <- 1
                Annot <- Annot - mean(Annot)
                Annot <- as.matrix(data.frame(rep(1, num_gene),
                  Annot))
                t1 <- system.time(model1 <- try(res <- EMMCMCStepSummary(object@summary[,
                  1], object@summary[, 2], as.matrix(Annot),
                  init_beta, init_tau, em_iter, mcmc_iter, min_degene)))
                if (class(model1) != "try-error") {
                  rownames(res$pip) <- object@gene_id
                  colnames(res$pip) <- "PIP"
                  rownames(res$beta) <- object@gene_id
                  colnames(res$beta) <- "beta"
                  rownames(res$annot_coef) <- c("tau_1", "tau_2")
                  colnames(res$annot_coef) <- "annot_coef"
                  rownames(res$annot_var) <- c("tau_1", "tau_2")
                  colnames(res$annot_var) <- "annot_var"
                  res$converged <- TRUE
                  res$ctime <- t1[3]
                }
                else {
                  res <- NULL
                }
                return(res)
            }, mc.cores = getOption("mc.cores", num_core))
        }
        else {
            if (verbose) {
                cat(paste("## fitting the iDEA variant model with gene sets information... \n"))
            }
            res_idea <- pbmclapply(1:num_annot, FUN = function(x) {
                Annot <- rep(0, object@num_gene)
                Annot[object@annotation[[x]]] <- 1
                Annot <- Annot - mean(Annot)
                Annot <- as.matrix(data.frame(rep(1, num_gene),
                  Annot))
                t1 <- system.time(model1 <- try(res <- EMMCMCStepSummaryVariant(object@summary[,
                  1], object@summary[, 2], as.matrix(Annot),
                  init_beta, init_tau, em_iter, mcmc_iter, min_degene)))
                if (class(model1) != "try-error") {
                  rownames(res$pip) <- object@gene_id
                  colnames(res$pip) <- "PIP"
                  rownames(res$beta) <- object@gene_id
                  colnames(res$beta) <- "beta"
                  rownames(res$annot_coef) <- c("tau_1", "tau_2")
                  colnames(res$annot_coef) <- "annot_coef"
                  rownames(res$annot_var) <- c("tau_1", "tau_2")
                  colnames(res$annot_var) <- "annot_var"
                  res$converged <- TRUE
                  res$ctime <- t1[3]
                }
                else {
                  res <- NULL
                }
                return(res)
            }, mc.cores = getOption("mc.cores", num_core))
        }
    }
    names(res_idea) <- object@annot_id
    object@de <- res_idea
    rm(res_idea)
    return(object)
}