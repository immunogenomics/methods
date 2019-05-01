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

rank_matrix <- function(X) {
    UseMethod('rank_matrix')
}

rank_matrix.dgCMatrix <- function(X) {
    Xr <- Matrix(X, sparse = TRUE)
    cpp_rank_matrix_dgc(Xr@x, Xr@p, nrow(Xr), ncol(Xr))
    return(Xr)
}

rank_matrix.matrix <- function(X) {
    return(cpp_rank_matrix_dense(X))
}



#' X is samples by features
#' 
#' assumes that y is a factor with no missing values
sumGroups <- function(X, y, MARGIN=2) {
    if (MARGIN == 2 & nrow(X) != length(y)) {
        stop('wrong dims')
    } else if (MARGIN == 1 & ncol(X) != length(y)) {
        stop('wrong dims')        
    }
    UseMethod('sumGroups')
}

sumGroups.dgCMatrix <- function(X, y, MARGIN=2) {
    if (MARGIN == 1) {
        cpp_sumGroups_dgc_T(X@x, X@p, X@i, ncol(X), nrow(X), as.integer(y) - 1, length(unique(y)))        
    } else {
        cpp_sumGroups_dgc(X@x, X@p, X@i, ncol(X), as.integer(y) - 1, length(unique(y)))        
    }
}

sumGroups.matrix <- function(X, y, MARGIN=2) {
    if (MARGIN == 1) {
        cpp_sumGroups_dense_T(X, as.integer(y) - 1, length(unique(y)))        
    } else {
        cpp_sumGroups_dense(X, as.integer(y) - 1, length(unique(y)))
    }
}



#' X is samples by features
#' 
#' assumes that y is a factor with no missing values
nnzeroGroups <- function(X, y, MARGIN=2) {
    if (MARGIN == 2 & nrow(X) != length(y)) {
        stop('wrong dims')
    } else if (MARGIN == 1 & ncol(X) != length(y)) {
        stop('wrong dims')        
    }
    UseMethod('nnzeroGroups')
}

nnzeroGroups.dgCMatrix <- function(X, y, MARGIN=2) {
    if (MARGIN == 1) {
        cpp_nnzeroGroups_dgc_T(X@p, X@i, ncol(X), nrow(X), as.integer(y) - 1, length(unique(y)))        
    } else {
        cpp_nnzeroGroups_dgc(X@p, X@i, ncol(X), as.integer(y) - 1, length(unique(y)))
    }
}

nnzeroGroups.matrix <- function(X, y, MARGIN=2) {
    if (MARGIN == 1) {
        cpp_sumGroups_dense_T(X != 0, as.integer(y) - 1, length(unique(y)))        
    } else {
        cpp_sumGroups_dense(X != 0, as.integer(y) - 1, length(unique(y)))
    }
}
