# t <- Matrix:::t


fast_wilcox <- function(Xt, cols) {
    if (is(Xt, 'matrix')) {
        fast_wilcox_dense(Xt, cols)        
    } else if (is(Xt, 'dgeMatrix') | is(Xt, 'data.frame') | is(Xt, 'DataFrame') | is(Xt, 'data.table')) {
        fast_wilcox_dense(as.matrix(Xt), cols)
    } else if (is(Xt, 'dgCMatrix')) {
        fast_wilcox_sparse(Xt, cols)
    } else if (is(Xt, 'dgTMatrix') | is(Xt, 'TsparseMatrix')) {
        fast_wilcox_sparse(as(Xt, 'dgCMatrix'), cols)
    }
}




fast_wilcox_sparse <- function(Xt, cols) {
    if (any(Xt < 0)) {
        stop("ERROR: Wilcox function does not support negative values yet.")
    }
    Xr <- rank_matrix(Xt)  
    grs <- sumGroups(Xr, cols)

    # calculate number of non-zero entries per group
    gnzz <- nnzeroGroups(Xr, cols)
    group.size <- as.numeric(table(cols))

    # add contribution of zero entries to the grs
    gnz <- (group.size - gnzz)

    # rank of a 0 entry for each gene
    zero.ranks <- (nrow(Xr) - diff(Xr@p) + 1) / 2
    ustat <- t((t(gnz) * zero.ranks)) + grs - group.size * (group.size + 1 ) / 2
        
    # standardize to get Z-score and pval
    n1n2 <- group.size * (nrow(Xr) - group.size);
    .x <- t(matrix(rep(colSums(gnz) ^ 3 - colSums(gnz), nrow(gnz)), ncol(gnz), nrow(gnz)))
    usigma <- sqrt((nrow(Xr) + 1 - .x / (nrow(Xr) * (nrow(Xr) - 1))) * n1n2 / 12)

    ustat_norm <- t((ustat - (n1n2 / 2)) / usigma)
     pvals <- ustat_norm %>% abs %>% as.numeric %>% 
        pnorm(lower.tail = FALSE, log.p = TRUE) %>% 
        bh.adjust(log = TRUE) %>% 
#         qnorm(lower.tail = FALSE, log.p = TRUE) %>% 
        matrix(ncol = ncol(ustat_norm))
    pvals <- -pvals * sign(ustat_norm)
    rownames(pvals) <- colnames(Xr)
    colnames(pvals) <- levels(cols)[1:ncol(pvals)]
  
    ## compute auROC from ustat
    auc = matrix(t(ustat / n1n2), ncol = ncol(pvals))
    rownames(auc) <- colnames(Xr)
    colnames(auc) <- levels(cols)[1:ncol(auc)]
    
    return(list(fdr = pvals, auc = auc))
}



fast_wilcox_dense <- function(Xt, cols) {
    Xr <- rank_matrix(Xt)  
    grs <- sumGroups(Xr, cols)
    group.size <- as.numeric(table(cols))
    ustat <- grs - group.size * (group.size + 1 ) / 2
        
    # standardize to get Z-score and pval
    n1n2 <- group.size * (nrow(Xr) - group.size);
    usigma <- sqrt((nrow(Xr) + 1) * n1n2 / 12)
    ustat_norm <- t(sweep(ustat - (n1n2 / 2), 1, usigma, '/'))
     pvals <- ustat_norm %>% abs %>% as.numeric %>% 
        pnorm(lower.tail = FALSE, log.p = TRUE) %>% 
        bh.adjust(log = TRUE) %>% 
#         qnorm(lower.tail = FALSE, log.p = TRUE) %>% 
        matrix(ncol = ncol(ustat_norm))
    pvals <- -pvals * sign(ustat_norm)
    rownames(pvals) <- colnames(Xr)
    colnames(pvals) <- levels(cols)[1:ncol(pvals)]
  
    ## compute auROC from ustat
    auc = matrix(t(ustat / n1n2), ncol = ncol(pvals))
    rownames(auc) <- colnames(Xr)
    colnames(auc) <- levels(cols)[1:ncol(auc)]
    
    return(list(fdr = pvals, auc = auc))
}

