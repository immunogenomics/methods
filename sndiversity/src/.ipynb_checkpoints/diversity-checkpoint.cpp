#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

float Hbeta(arma::mat& D, float beta, arma::vec& P, int idx) {
  P = arma::exp(-D.col(idx) * beta);
  float sumP = sum(P);
  float H;
  if (sumP == 0){
      H = 0;
      P = D.col(idx) * 0;
  } else {
      H = log(sumP) + beta * sum(D.col(idx) % P) / sumP;
      P /= sumP;
  }
  return(H);
}

// [[Rcpp::export]]
arma::vec compute_diversity_Cpp(arma::mat& D, arma::umat& knn_idx, arma::vec& batch_labels, int n_batches,
                float perplexity = 15, float tol = 1e-5) {
  int n = D.n_cols;
  arma::vec P = zeros<arma::vec>(D.n_rows);
  arma::vec diversity = zeros<arma::vec>(n);
  float logU = log(perplexity);

  float hbeta, beta, betamin, betamax, H, Hdiff;
  int tries;
  for (int i = 0; i < n ; i++) {
    beta = 1;
    betamin = -datum::inf;
    betamax = datum::inf;
    H = Hbeta(D, beta, P, i);
    Hdiff = H - logU;
    tries = 0;
    // first get neighbor probabilities
    while(std::abs(Hdiff) > tol && tries < 50) {
      if (Hdiff > 0){
        betamin = beta;
        if (!is_finite(betamax)) beta *= 2;
        else beta = (beta + betamax) / 2;
      } else{
        betamax = beta;
        if (!is_finite(betamin)) beta /= 2;
        else beta = (beta + betamin) / 2;
      }

      H = Hbeta(D, beta, P, i);
      Hdiff = H - logU;
      tries++;
    }

    // then compute diversity
    for (int b = 0; b < n_batches; b++) {
      uvec q = find(batch_labels.elem(knn_idx.col(i)) == b); // indices of cells belonging to batch (b)
      if (q.n_elem > 0) {
        float sumP = sum(P.elem(q));
        diversity.row(i) += sumP * sumP;         
      }
    }
  }
  return(diversity);
}
