t <- Matrix:::t

## BH P-value adjustment with a log option
bh.adjust <- function(x, log = FALSE) {
  nai <- which(!is.na(x))
  ox <- x
  x <- x[nai]
  id <- order(x, decreasing = FALSE)
  if(log) {
    q <- x[id] + log(length(x)/seq_along(x))
  } else {
    q <- x[id]*length(x)/seq_along(x)
  }
  a <- rev(cummin(rev(q)))[order(id)]
  ox[nai] <- a
  ox
}

fast_wilcox <- function(Xt, cols) {
    Xr = Matrix(Xt, sparse = TRUE)
##    Xr@x <- as.numeric(rank_vec(Xr@x, Xr@p, nrow(Xr), ncol(Xr), "mean"))
    rank_vec(Xr@x, Xr@p, nrow(Xr), ncol(Xr), "mean")  
  
    grs <- sumGroups(Xr@x, Xr@p, Xr@i, ncol(Xr), as.integer(cols) - 1, length(unique(cols)))

    # calculate number of non-zero entries per group
    gnzz <- nnzeroGroups(Xr@p, Xr@i, ncol(Xr), as.integer(cols) - 1, length(unique(cols)))
    group.size <- as.numeric(table(cols))

    # add contribution of zero entries to the grs
    gnz <- (group.size-gnzz)

    # rank of a 0 entry for each gene
    zero.ranks <- (nrow(Xr) - diff(Xr@p) + 1) / 2
    ustat <- t((t(gnz) * zero.ranks)) + grs - group.size * (group.size + 1 ) / 2

    # standardize to get Z-score and pval
    n1n2 <- group.size * (nrow(Xr) - group.size);
##    usigma <- sqrt((nrow(Xr) + 1 - (colSums(gnz)^3 - colSums(gnz)) / (nrow(Xr) * (nrow(Xr) - 1))) * n1n2 / 12)  
    usigma <- sqrt((nrow(Xr) + 1 - (gnz ^ 3 - gnz) / (nrow(Xr) * (nrow(Xr) - 1))) * n1n2 / 12)
    ustat_norm <- t((ustat - (n1n2 / 2)) / usigma)    
    pvals <- ustat_norm %>% abs %>% as.numeric %>% 
        pnorm(lower.tail = FALSE, log.p = TRUE) %>%
        bh.adjust(log = TRUE) %>% 
#        qnorm(lower.tail = FALSE, log.p = TRUE) %>% 
        matrix(ncol = ncol(ustat_norm))
    pvals <- pvals * sign(ustat_norm)
    rownames(pvals) <- colnames(Xr)
    colnames(pvals) <- levels(cols)[1:ncol(pvals)]
    
    ## compute auROC from ustat
    auc = matrix(t(ustat / n1n2), ncol = ncol(pvals))
    rownames(auc) <- colnames(Xr)
    colnames(auc) <- levels(cols)[1:ncol(auc)]
    
    return(list(fdr = pvals, auc = auc))
}


fast_diff_exp <- function(X, y) {
  if (!class(X) == "dgCMatrix") stop(sprintf("ERROR: X must be a dgCMatrix"))
  if (ncol(X) != length(y)) stop("ERROR: number of columns of X does not match length of y")
  
  wilcox_res <- fast_wilcox(t(X), y)
  
  ### Auxiliary Statistics (AvgExpr, PctIn, LFC)
  group_sums <- sumGroupsT(X@x, X@p, X@i, ncol(X), nrow(X), as.integer(y) - 1, length(unique(y)))
  group_nnz <- nnzeroGroupsT(X@x, X@p, X@i, ncol(X), nrow(X), as.integer(y) - 1, length(unique(y))) 
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
