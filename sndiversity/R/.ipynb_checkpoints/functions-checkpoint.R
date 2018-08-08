#' Function to compute Simpson's Index for each data point, given a kNN
#' @export
compute_simpson_from_knn <- function(dknn, batch_labels) {
  perplexity <- ncol(dknn$nn.idx) / 3
  batch_labels <- as.integer(factor(batch_labels)) - 1
  n_batches <- length(unique(batch_labels))
  lambda <- compute_diversity_Cpp(t(dknn$nn.dists), t(dknn$nn.idx) - 1, 
                                     batch_labels, n_batches, perplexity)
  ## compute maximal theoretical diversity
#  max_diversity <- 1 - (n_batches * ((1 / n_batches) ^ 2))
  ## scale diversity from 0 to 1
#  diversity <- diversity - min(diversity)
#  diversity <- diversity / max_diversity
  return(lambda)
}

#' Function to compute diversity score for each data point
#' @param X data matrix, samples (rows) by features (columns) (e.g. matrix(runif(500), 100, 5))
#' @param batch_labels 1D array of discrete batch assignment (e.g. c("A", "A", "B", "C", "A", "B"))
#' @param perplexity effective number of neighbors (defined in SNE/tSNE algorithms)
#' @param ratio_use ratio of data to use. 1 by default. Faster with smaller ratio. 
#' @param nn_eps noise parameter for kNN computation. 0 is exact, higher values are approximations.  
#' @export
compute_simpson_from_features <- function(X, batch_labels, perplexity = 30, ratio_use = 1, nn_eps = 0) {
  if (ratio_use < 1) {
    samples_use <- sample.int(nrow(X), nrow(X) * ratio_use)
    X <- X[samples_use, ]
    batch_labels <- batch_labels[samples_use]
  }

  ## TODO: extend to non-discrete batch with P * Phi
  dknn <- RANN::nn2(data = X, k = 3 * perplexity, eps = nn_eps)  
  lambda <- compute_simpson_from_knn(dknn, batch_labels)

  return(lambda)
}

