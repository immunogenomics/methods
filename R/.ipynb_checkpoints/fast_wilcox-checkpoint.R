fast_wilcox <- function(Xt, cols, correct=TRUE) {
#     if (is(Xt, 'matrix')) {
#         fast_wilcox_dense(Xt, cols)        
#     } else if (is(Xt, 'dgeMatrix') | is(Xt, 'data.frame') | is(Xt, 'DataFrame') | is(Xt, 'data.table')) {
#         fast_wilcox_dense(as.matrix(Xt), cols)
#     } else if (is(Xt, 'dgCMatrix')) {
#         fast_wilcox_sparse(Xt, cols)
#     } else if (is(Xt, 'dgTMatrix') | is(Xt, 'TsparseMatrix')) {
#         fast_wilcox_sparse(as(Xt, 'dgCMatrix'), cols)
#     }

    if (is(Xt, 'dgeMatrix')) Xt <- as.matrix(Xt)
    if (is(Xt, 'data.frame')) Xt <- as.matrix(Xt)
    if (is(Xt, 'DataFrame')) Xt <- as.matrix(Xt)
    if (is(Xt, 'data.table')) Xt <- as.matrix(Xt)
    if (is(Xt, 'dgTMatrix')) Xt <- as(Xt, 'dgCMatrix')
    if (is(Xt, 'TsparseMatrix')) Xt <- as(Xt, 'dgCMatrix')
    
    
    if (any(Xt < 0) & is(Xt, 'dgCMatrix')) {
        stop("ERROR: Wilcox function does not support negative values yet.")
    }
    rank_res <- rank_matrix(Xt)
    Xr <- rank_res$X_ranked
    grs <- sumGroups(Xr, cols)

    # calculate number of non-zero entries per group
    group.size <- as.numeric(table(cols))
    n1n2 <- group.size * (nrow(Xr) - group.size);

    if (is(Xt, 'dgCMatrix')) {
        gnz <- (group.size - nnzeroGroups(Xr, cols))
        zero.ranks <- (nrow(Xr) - diff(Xr@p) + 1) / 2
        ustat <- t((t(gnz) * zero.ranks)) + grs - group.size * (group.size + 1 ) / 2        
    } else {
        ustat <- grs - group.size * (group.size + 1 ) / 2
    }
    
    # standardize to get Z-score and pval
    z <- ustat - .5 * n1n2
    if (correct) {
        z <- z - sign(z) * .5
    }
    usigma <- compute_usigma(rank_res$ties, nrow(Xr), n1n2)
    z <- t(z / usigma)
    
    pvals <- matrix(2 * pnorm(-abs(as.numeric(z))), ncol = ncol(z))
    rownames(pvals) <- colnames(Xr)
    colnames(pvals) <- levels(cols)[1:ncol(pvals)]  
    fdr <- apply(pvals, 2, function(x) p.adjust(x, 'BH'))
    
    ## compute auROC from ustat
    auc = matrix(t(ustat / n1n2), ncol = ncol(pvals))
    rownames(auc) <- colnames(Xr)
    colnames(auc) <- levels(cols)[1:ncol(auc)]
    
    return(list(fdr = fdr, pvals = pvals, auc = auc, ustat = ustat))    
}




# fast_wilcox_sparse <- function(Xt, cols, correct=TRUE) {
#     if (any(Xt < 0)) {
#         stop("ERROR: Wilcox function does not support negative values yet.")
#     }
#     rank_res <- rank_matrix(Xt)
#     Xr <- rank_res$X_ranked
# #     ties <- rank_res$ties
#     grs <- sumGroups(Xr, cols)
    
#     # calculate number of non-zero entries per group
#     gnzz <- nnzeroGroups(Xr, cols)
#     group.size <- as.numeric(table(cols))
#     n1n2 <- group.size * (nrow(Xr) - group.size);

#     # add contribution of zero entries to the grs
#     gnz <- (group.size - gnzz)
#     zero.ranks <- (nrow(Xr) - diff(Xr@p) + 1) / 2
#     ustat <- t((t(gnz) * zero.ranks)) + grs - group.size * (group.size + 1 ) / 2
    
#     # standardize to get Z-score and pval
#     z <- ustat - .5 * n1n2
#     if (correct) {
#         z <- z - sign(z) * .5
#     }
#     usigma <- compute_usigma(rank_res$ties, nrow(Xr), n1n2)
#     z <- t(z / usigma)
    
#     pvals <- matrix(2 * pnorm(-abs(as.numeric(z))), ncol = ncol(z))
#     rownames(pvals) <- colnames(Xr)
#     colnames(pvals) <- levels(cols)[1:ncol(pvals)]  
#     fdr <- apply(pvals, 2, function(x) p.adjust(x, 'BH'))
  
#     ## compute auROC from ustat
#     auc = matrix(t(ustat / n1n2), ncol = ncol(pvals))
#     rownames(auc) <- colnames(Xr)
#     colnames(auc) <- levels(cols)[1:ncol(auc)]
    
#     return(list(fdr = pvals, auc = auc))
# }



# fast_wilcox_dense <- function(Xt, cols, correct=TRUE) {
#     rank_res <- rank_matrix(Xt)
#     Xr <- rank_res$X_ranked
    
#     grs <- sumGroups(Xr, cols)
#     group.size <- as.numeric(table(cols))
#     n1n2 <- group.size * (nrow(Xr) - group.size);
#     ustat <- grs - group.size * (group.size + 1 ) / 2

#     # standardize to get Z-score and pval
#     z <- ustat - .5 * n1n2
#     if (correct) {
#         z <- z - sign(z) * .5
#     }
#     usigma <- compute_usigma(rank_res$ties, nrow(Xr), n1n2)
#     z <- t(z / usigma)
    
#     pvals <- matrix(2 * pnorm(-abs(as.numeric(z))), ncol = ncol(z))
#     rownames(pvals) <- colnames(Xr)
#     colnames(pvals) <- levels(cols)[1:ncol(pvals)]  
#     fdr <- apply(pvals, 2, function(x) p.adjust(x, 'BH'))
    
#     ## compute auROC from ustat
#     auc = matrix(t(ustat / n1n2), ncol = ncol(pvals))
#     rownames(auc) <- colnames(Xr)
#     colnames(auc) <- levels(cols)[1:ncol(auc)]
    
#     return(list(fdr = fdr, pvals = pvals, auc = auc, ustat = ustat))
# }

