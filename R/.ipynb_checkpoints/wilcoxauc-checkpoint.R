#' Fast Wilcoxon rank sum test and auROC 
#' 
#' Computes auROC and Wilcoxon p-value based on Gaussian approximation 
#' 
#' @param X A feature-by-sample matrix, Seurat object, or SingleCellExperiment object
#' @param y vector of group labels. 
#' @param groups_use (optional) which groups from y vector to test. 
#' @param group_by (Seurat & SCE) name of groups variable ('e.g. Cluster').
#' @param assay (Seurat & SCE) name of feature matrix slot (e.g. 'data' or 'logcounts'). 
#' @param seurat_assay (Seurat) name of Seurat Assay (e.g. 'RNA'). 
#' @param verbose boolean, TRUE for warnings and messages. 
#' 
#' @examples
#' 
#' library(presto)
#' data(exprs)
#' data(y)
#' 
#' @return table with the following columns: 
#' \itemize{
#' \item 
#' \item 
#' \item 
#' }
#' @export 
wilcoxauc <- function(X, ...) {
    UseMethod('wilcoxauc')
}

wilcoxauc.seurat <- function(X, ...) {
    stop('wilcoxauc only implemented for Seurat Version 3, please upgrade to run.')
}

wilcoxauc.Seurat <- function(object, group_by=NULL, assay='data', seurat_assay='RNA', groups_use=NULL) {
    X_matrix <- Seurat::GetAssayData(object, assay=seurat_assay, slot=assay)
    if (is.null(group_by)) {
        y <- Seurat::Idents(object)
    } else {
        y <- Seurat::FetchData(object, group_by) %>% unlist %>% as.character()        
    }
#     return(y)
    wilcoxauc(X_matrix, y, groups_use)
}


# wilcoxauc.CellDataSet <- function(object, slot='exprs', group_by='Cluster') {
#     X <- assayData(object)[[slot]]
#     y <- phenoData(object)[[group_by]]
#     wilcoxauc(X, y)    
# }

wilcoxauc.SingleCellExperiment <- function(object, group_by=NULL, assay=NULL, groups_use=NULL) {
    if (is.null(group_by)) {
        stop('Must specify group_by with SingleCellExperiment')
    } else if (!group_by %in% names(colData(object))) {
        stop('group_by value is not defined in colData.')
    }
    y <- colData(object)[[group_by]]
    
    if (is.null(assay)) {
        standard_assays <- c('normcounts', 'logcounts', 'cpm', 'tpm', 'weights', 'counts')
        standard_assays <- factor(standard_assays, standard_assays)
        available_assays <- names(assays(object))    
        available_assays <- intersect(standard_assays, available_assays)
        if (length(available_assays) == 0) {
            stop('No assays in SingleCellExperiment object')
        } else {
            assay <- available_assays[1]
        }    
    } 
    
    X_matrix <- eval(call(assay, object))

    wilcoxauc(X_matrix, y, groups_use)
}

wilcoxauc.default <- function(X, y, groups_use=NULL, verbose=TRUE) {
    ## Check and possibly correct input values
    if (is(X, 'dgeMatrix')) X <- as.matrix(X)
    if (is(X, 'data.frame')) X <- as.matrix(X)
    if (is(X, 'DataFrame')) X <- as.matrix(X)
    if (is(X, 'data.table')) X <- as.matrix(X)
    if (is(X, 'dgTMatrix')) X <- as(X, 'dgCMatrix')
    if (is(X, 'TsparseMatrix')) X <- as(X, 'dgCMatrix')
    
    if (ncol(X) != length(y)) stop("ERROR: number of columns of X does not match length of y")
    if (!is.null(groups_use)) {
        idx_use <- which(y %in% groups_use)
        y <- y[idx_use]
        X <- X[, idx_use]
    }
    
    y <- factor(y)
    idx_use <- which(!is.na(y))
    if (length(idx_use) < length(y)) {
        y <- y[idx_use]
        X <- X[, idx_use]
        if (verbose) 
            message('Removing NA values from labels')        
    }
    
    
#     features_use <- which(apply(!is.na(X), 1, all))
#     if (verbose & length(features_use) < nrow(X)) {
#         message('Removing features with NA values')
#     }
#     X <- X[features_use, ]
    if (is.null(row.names(X))) {
        row.names(X) <- paste0('Feature', seq_len(nrow(X)))
    }
    
    ## Compute primary statistics
    group.size <- as.numeric(table(y))
    n1n2 <- group.size * (ncol(X) - group.size);
    if (is(X, 'dgCMatrix')) {
        rank_res <- rank_matrix(Matrix:::t(X))        
    } else {
        rank_res <- rank_matrix(X)
#         rank_res <- rank_matrix(t(X))
    }

    ustat <- compute_ustat(rank_res$X_ranked, y, n1n2, group.size) 
    auc <- t(matrix(t(ustat / n1n2), ncol = ncol(ustat)))
    pvals <- compute_pval(ustat, rank_res$ties, ncol(X), n1n2) 
    fdr <- apply(pvals, 2, function(x) p.adjust(x, 'BH'))

    ### Auxiliary Statistics (AvgExpr, PctIn, LFC, etc)
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
                statistic = t(ustat),
                logFC = lfc)
    return(tidy_results(res_list, row.names(X), levels(y)))
}


