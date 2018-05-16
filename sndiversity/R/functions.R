#' Function to compute diversity score for each data point
#' @param X data matrix, samples (rows) by features (columns) (e.g. matrix(runif(500), 100, 5))
#' @param batch_labels 1D array of discrete batch assignment (e.g. c("A", "A", "B", "C", "A", "B"))
#' @param perplexity effective number of neighbors (defined in SNE/tSNE algorithms)
#' @param ratio_use ratio of data to use. 1 by default. Faster with smaller ratio. 
#' @export
compute_batch_diversity <- function(X, batch_labels, perplexity = 30, ratio_use = 1) {
  if (ratio_use < 1) {
    samples_use <- sample.int(nrow(X), nrow(X) * ratio_use)
    X <- X[samples_use, ]
    batch_labels <- batch_labels[samples_use]
  }

  ## TODO: extend to non-discrete batch with P * Phi
  dknn <- FNN::get.knn(X, k = 3 * perplexity)
  batch_labels <- as.integer(factor(batch_labels)) - 1
  n_batches <- length(unique(batch_labels))
  diversity <- compute_diversity_Cpp(t(dknn$nn.dist), t(dknn$nn.index) - 1, 
                                     batch_labels, n_batches, perplexity)

  ## compute maximal theoretical diversity
  max_diversity <- 1 - (n_batches * ((1/n_batches) ^ 2))

  ## scale diversity from 0 to 1
  diversity <- diversity - min(diversity)
  diversity <- diversity / max_diversity
  return(diversity)
}

