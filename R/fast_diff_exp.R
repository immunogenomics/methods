
fast_diff_exp <- function(X, y) {
    if (is(X, 'dgeMatrix')) X <- as.matrix(X)
    if (is(X, 'data.frame')) X <- as.matrix(X)
    if (is(X, 'DataFrame')) X <- as.matrix(X)
    if (is(X, 'data.table')) X <- as.matrix(X)
    if (is(X, 'dgTMatrix')) X <- as(X, 'dgCMatrix')
    if (is(X, 'TsparseMatrix')) X <- as(X, 'dgCMatrix')
    
    
    if (ncol(X) != length(y)) stop("ERROR: number of columns of X does not match length of y")
    wilcox_res <- fast_wilcox(Matrix:::t(X), y)
  
    ### Auxiliary Statistics (AvgExpr, PctIn, LFC)
    group_sums <- sumGroups(X, y, 1)
    group_nnz <- nnzeroGroups(X, y, 1)
    group_pct <- sweep(group_nnz, 1, as.numeric(table(y)), "/")
    group_pct_out <- -group_nnz %>% 
    sweep(2, colSums(group_nnz) , "+") %>% 
    sweep(1, as.numeric(length(y) - table(y)), "/")      
    

    ## Group means can be over (1) all cells or (2) expressing cells
    group_means <- sweep(group_sums, 1, as.numeric(table(y)), "/")
    group_means_expressing <- group_sums / group_nnz
  
    cs <- colSums(group_sums)
    gs <- as.numeric(table(y))
    lfc <- Reduce(rbind, lapply(1:length(levels(y)), function(g) {
        group_means[g, ] - ((cs - group_sums[g, ]) / (length(y) - gs[g]))
    }))

    group_pct %<>% t %>% data.frame()
    group_pct_out %<>% t %>% data.frame()
    group_means %<>% t %>% data.frame()
    group_means_expressing %<>% t %>% data.frame()
    lfc %<>% t %>% data.frame()

    row.names(group_pct) <- row.names(X)
    row.names(group_pct_out) <- row.names(X)
    row.names(group_means) <- row.names(X)
    row.names(group_means_expressing) <- row.names(X)
    row.names(lfc) <- row.names(X)

    colnames(group_pct) <- levels(y)
    colnames(group_pct_out) <- levels(y)
    colnames(group_means) <- levels(y)
    colnames(group_means_expressing) <- levels(y)
    colnames(lfc) <- levels(y)
  
    return(list(auc = wilcox_res$auc, fdr = wilcox_res$fdr, 
                pct_in = group_pct, pct_out = group_pct_out,
                means = group_means, 
                means_nz = group_means_expressing, lfc = lfc))
}
