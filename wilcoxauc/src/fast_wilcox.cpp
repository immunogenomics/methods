#include <iostream>
#include <vector>
#include <algorithm>
#include <omp.h>
#include "wilcoxauc_types.h"
//#include <progress.hpp>
//#include <progress_bar.hpp>


using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace std;

// [[Rcpp::export]]
int inPlaceRankMean(vec& v_temp, int idx_begin, int idx_end) {
    if (idx_begin > idx_end) return 1;
    std::vector<pair<float, size_t> > v_sort(idx_end - idx_begin + 1);
    for (size_t i = idx_begin; i <= idx_end; i++) {
        v_sort[i - idx_begin] = make_pair(v_temp[i], i - idx_begin);
    }
    
    sort(v_sort.begin(), v_sort.end());

    float rank_sum = 0, n = 1;    
    size_t i;
    for (i = 1U; i < v_sort.size(); i++) {
        if (v_sort[i].first != v_sort[i - 1].first) {
            // if current val != prev val
            // set prev val to something
            for (unsigned j = 0; j < n; j++) {
                v_temp[v_sort[i - 1 - j].second + idx_begin] = (rank_sum / n) + 1;  
            }            
            // restart count ranks
            rank_sum = i;
            n = 1;
        } else {
            // if curr val is a tie, 
            // don't set anything yet, start computing mean rank
            rank_sum += i;
            n++;
        }
    }
    // set the last element(s)
    for (unsigned j = 0; j < n; j++)
        v_temp[v_sort[i - 1 - j].second + idx_begin] = (rank_sum / n) + 1;  
  return 0;
}

// [[Rcpp::export]]
int inPlaceRankMin(vec& v_temp, int idx_begin, int idx_end) {
    if (idx_begin > idx_end) return 1;
    std::vector<pair<float, size_t> > v_sort(idx_end - idx_begin + 1);
    for (size_t i = idx_begin; i <= idx_end; i++) {
        v_sort[i - idx_begin] = make_pair(v_temp[i], i - idx_begin);
    }
    
    sort(v_sort.begin(), v_sort.end());

    float rank_min = 0, n = 1;    
    size_t i;
    for (i = 1U; i < v_sort.size(); i++) {
        if (v_sort[i].first != v_sort[i - 1].first) {
            // if current val != prev val
            // set prev val to something
            for (unsigned j = 0; j < n; j++) {
                v_temp[v_sort[i - 1 - j].second + idx_begin] = rank_min + 1;  
            }            
            // restart count ranks
            rank_min = i;
            n = 1;
        } else {
            // if curr val is a tie, 
            // don't set anything yet
            n++;
        }
    }
    // set the last element(s)
    for (unsigned j = 0; j < n; j++)
        v_temp[v_sort[i - 1 - j].second + idx_begin] = rank_min + 1;
  return 0;
}


// [[Rcpp::export]]
int rank_vec(vec& x, const vec& p, int nrow, int ncol, std::string method) {
  omp_set_num_threads(8);
  #pragma omp parallel for
  for (int i = 0; i < ncol; i++) {
    if (p[i+1] == p[i]) continue;
    int n_zero = nrow - (p[i+1] - p[i]);
    if (method == "mean") {
        inPlaceRankMean(x, p[i], p[i + 1] - 1);        
    } else if (method == "min") {
        inPlaceRankMin(x, p[i], p[i + 1] - 1);        
    }
    x.rows(p[i], p[i + 1] - 1) += n_zero;
  }
//  return x;
  return 0;
}


// [[Rcpp::export]]
mat sumGroups(const vec& x, const vec& p, const vec& i, int ncol, const uvec& groups, int ngroups) {
    // Here, columns are genes
    mat res = arma::zeros<mat>(ngroups, ncol);
    for (int c = 0; c < ncol; c++) {
        for (int j = p[c]; j < p[c + 1]; j++) {
            // i[j] gives the row num
            // group_map gives the group num of that row
            res(groups[i[j]], c) += x[j];
        }
    }    
    return res;
}



// [[Rcpp::export]]
mat nnzeroGroups(const vec& p, const vec& i, int ncol, const uvec& groups, int ngroups) {
    mat res = arma::zeros<mat>(ngroups, ncol);
    for (int c = 0; c < ncol; c++) {
        for (int j = p[c]; j < p[c + 1]; j++) {
            res(groups[i[j]], c)++;
        }
    }    
    return res;
}


// [[Rcpp::export]]
mat sumGroupsT(const vec& x, const vec& p, const vec& i, int ncol, int nrow, const uvec& groups, int ngroups) {
    // Here, columns are samples
    mat res = arma::zeros<mat>(ngroups, nrow);
    for (int c = 0; c < ncol; c++) {
        for (int j = p[c]; j < p[c + 1]; j++) {
            // i[j] gives the row num
            // group_map gives the group num of that row
            res(groups[c], i[j]) += x[j];
        }
    }    
    return res;
}

// [[Rcpp::export]]
mat nnzeroGroupsT(const vec& x, const vec& p, const vec& i, int ncol, int nrow, const uvec& groups, int ngroups) {
    // Here, columns are samples
    mat res = arma::zeros<mat>(ngroups, nrow);
    for (int c = 0; c < ncol; c++) {
        for (int j = p[c]; j < p[c + 1]; j++) {
            // i[j] gives the row num
            // group_map gives the group num of that row
            res(groups[c], i[j])++;
        }
    }    
    return res;
}
