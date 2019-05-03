#' Fast Mann Whitney U Test and auROC
#' 
#' Compute auROC and p-value based on Gaussian approximation of 
#' Mann Whitney U statistic
#' 
#' @param Xt sample by feature matrix. 
#' @param cols vector of group labels. 
#' @param verbose boolean, TRUE for warnings and messages. 
#' 
#' @examples
#' 
#' data(wilcoxauc)
#' res <- fast_diff_exp(exprs, y)
#' ## Also works on sparse matrices
#' res <- fast_diff_exp(as(exprs, 'dgCMatrix'), y)
#' head(res)
#' 
#' @return table with the following columns: 
#' \itemize{
#' \item 
#' \item 
#' \item 
#' }
#' @export 
fast_diff_exp <- function(X, y, verbose=TRUE) {
    ## Check and possibly correct input values
    if (is(X, 'dgeMatrix')) X <- as.matrix(X)
    if (is(X, 'data.frame')) X <- as.matrix(X)
    if (is(X, 'DataFrame')) X <- as.matrix(X)
    if (is(X, 'data.table')) X <- as.matrix(X)
    if (is(X, 'dgTMatrix')) X <- as(X, 'dgCMatrix')
    if (is(X, 'TsparseMatrix')) X <- as(X, 'dgCMatrix')
    
    if (ncol(X) != length(y)) stop("ERROR: number of columns of X does not match length of y")
    y <- factor(y)
    idx_use <- which(!is.na(y))
    if (verbose & length(idx_use) < length(y)) {
        message('Removing NA values from labels')
    }
    y <- y[idx_use]
    X <- X[, idx_use]    
    
    ## Compute primary statistics
    group.size <- as.numeric(table(y))
    n1n2 <- group.size * (ncol(X) - group.size);
    if (is(X, 'dgCMatrix')) {
        rank_res <- rank_matrix(Matrix:::t(X))        
    } else {
        rank_res <- rank_matrix(t(X))        
    }

    ustat <- compute_ustat(rank_res$X_ranked, y, n1n2, group.size) 
    auc <- t(matrix(t(ustat / n1n2), ncol = ncol(ustat)))
    pvals <- compute_pval(ustat, rank_res$ties, ncol(X), n1n2) 
    fdr <- apply(pvals, 2, function(x) p.adjust(x, 'BH'))

    ### Auxiliary Statistics (AvgExpr, PctIn, LFC)
    group_sums <- sumGroups(X, y, 1)
    group_nnz <- nnzeroGroups(X, y, 1)
    group_pct <- sweep(group_nnz, 1, as.numeric(table(y)), "/") %>% t()
    group_pct_out <- -group_nnz %>% 
        sweep(2, colSums(group_nnz) , "+") %>% 
        sweep(1, as.numeric(length(y) - table(y)), "/") %>% t()
    group_means <- sweep(group_sums, 1, as.numeric(table(y)), "/") %>% t()
  
    cs <- colSums(group_sums)
    gs <- as.numeric(table(y))
    lfc <- Reduce(cbind, lapply(1:length(levels(y)), function(g) {
        group_means[, g] - ((cs - group_sums[g, ]) / (length(y) - gs[g]))
    }))

                 
    
    res_list <- list(auc = auc, 
                pval = pvals,
                padj = fdr, 
                pct_in = group_pct, 
                pct_out = group_pct_out,
                avgExpr = group_means, 
                logFC = lfc)
    return(tidy_results(res_list, row.names(X), levels(y)))
}


