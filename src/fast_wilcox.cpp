#include <iostream>
#include <vector>
#include <algorithm>
#include <omp.h>
#include "wilcoxauc_types.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace std;



// [[Rcpp::export]]
mat cpp_sumGroups_dgc(const vec& x, const uvec& p, const vec& i, unsigned ncol, const uvec& groups, unsigned ngroups) {
    // Here, columns are genes
    mat res = arma::zeros<mat>(ngroups, ncol);
    for (unsigned c = 0; c < ncol; c++) {
        for (unsigned j = p[c]; j < p[c + 1]; j++) {
            // i[j] gives the row num
            // group_map gives the group num of that row
            res(groups[i[j]], c) += x[j];
        }
    }    
    return res;
}

// [[Rcpp::export]]
mat cpp_sumGroups_dgc_T(const vec& x, const vec& p, const vec& i, int ncol, int nrow, const uvec& groups, int ngroups) {
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
mat cpp_sumGroups_dense(const mat& X, const uvec& groups, unsigned ngroups) {
    // Here, columns are genes
    mat res = arma::zeros<mat>(ngroups, X.n_cols);
    for (unsigned r = 0; r < X.n_rows; r++) {
    // for (unsigned c = 0; c < X.n_cols; c++) {
        res.row(groups[r]) += sum(X.row(r), 0);
    //    }
    }    
    return res;
}

// [[Rcpp::export]]
mat cpp_sumGroups_dense_T(const mat& X, const uvec& groups, unsigned ngroups) {
    // Here, columns are genes
    mat res = arma::zeros<mat>(ngroups, X.n_rows);
    for (unsigned c = 0; c < X.n_cols; c++) {
    // for (unsigned c = 0; c < X.n_cols; c++) {
        res.row(groups[c]) += sum(X.col(c), 1).t();
    //    }
    }    
    return res;
}



// [[Rcpp::export]]
mat cpp_nnzeroGroups_dgc(const uvec& p, const vec& i, unsigned ncol, const uvec& groups, unsigned ngroups) {
    mat res = arma::zeros<mat>(ngroups, ncol);
    for (unsigned c = 0; c < ncol; c++) {
        for (unsigned j = p[c]; j < p[c + 1]; j++) {
            res(groups[i[j]], c)++;
        }
    }    
    return res;
}


// [[Rcpp::export]]
std::list<float> cpp_in_place_rank_mean(vec& v_temp, int idx_begin, int idx_end) {
    std::list<float> ties;
    
    if (idx_begin > idx_end) return ties;
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
            if (n > 1) ties.push_back(n);
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

    return ties;
}


// [[Rcpp::export]]
std::vector<std::list<float> > cpp_rank_matrix_dgc(vec& x, const vec& p, int nrow, int ncol) { //, std::string method) {
//   omp_set_num_threads(8);
//   #pragma omp parallel for
    vector<list<float> > ties(ncol);
    int n_zero;
    for (int i = 0; i < ncol; i++) {
        if (p[i+1] == p[i]) continue;
        n_zero = nrow - (p[i+1] - p[i]);
        ties[i] = cpp_in_place_rank_mean(x, p[i], p[i + 1] - 1);
        ties[i].push_back(n_zero);
        x.rows(p[i], p[i + 1] - 1) += n_zero;
    }
    return ties;
}


// [[Rcpp::export]]
Rcpp::List cpp_rank_matrix_dense(mat X) {
    // sizes of tied groups
    vector<list<float> > ties(X.n_cols);
    
    for (unsigned c = 0; c < X.n_cols; c++) {
        std::vector<pair<float, size_t> > v_sort(X.n_rows);
        for (size_t i = 0; i < X.n_rows; i++) {
            v_sort[i] = make_pair(X.col(c)[i], i);
        }

        sort(v_sort.begin(), v_sort.end());

        float rank_sum = 0, n = 1;
        size_t i;
        for (i = 1U; i < v_sort.size(); i++) {
            if (v_sort[i].first != v_sort[i - 1].first) {
                // if current val != prev val
                // set prev val to something
                for (unsigned j = 0; j < n; j++) {
                    X.col(c)[v_sort[i - 1 - j].second] = (rank_sum / n) + 1;  
                }            
                // restart count ranks
                rank_sum = i;
                if (n > 1) ties[c].push_back(n);
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
            X.col(c)[v_sort[i - 1 - j].second] = (rank_sum / n) + 1;
    }
    return Rcpp::List::create(Named("X_ranked") = X, Named("ties") = ties);
}




// [[Rcpp::export]]
mat cpp_nnzeroGroups_dgc_T(const vec& p, const vec& i, int ncol, int nrow, const uvec& groups, int ngroups) {
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




/*
// [[Rcpp::export]]
int cpp_in_place_rank_min(vec& v_temp, int idx_begin, int idx_end) {
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

*/
