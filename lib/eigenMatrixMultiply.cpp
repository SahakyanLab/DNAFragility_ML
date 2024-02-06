// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

#define ARMA_WARN_LEVEL 0
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
#pragma GCC diagnostic ignored "-Wunused-variable"

#include <RcppArmadillo.h>
#include <RcppEigen.h>

// [[Rcpp::export]]
SEXP eigenMatrixMultiply(Eigen::MatrixXd A, Eigen::MatrixXd B){
    Eigen::MatrixXd C = A * B;
    return Rcpp::wrap(C);
}