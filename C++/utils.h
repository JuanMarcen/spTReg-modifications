#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>
#include <random>
#include <chrono>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

// Aux functions
void reportProgress(int iter, int nBurnin, int nSims, int nReport, int save_idx = 0,
                    std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now());
arma::vec rig(const int N, const arma::vec& mu, const double lambda);
arma::vec RandomMultiNormalC(const arma::mat& Q, const arma::vec& b);
arma::mat inv_covariance_matrix(
    const double precision, 
    const double decay, 
    const double varsigma, 
    const double varphi, 
    const double cprec,
    const arma::mat& dmat,
    const arma::vec& dvec,
    const arma::mat& dmatc
);
arma::mat covariance_matrix(
    const double precision, //hp(0, m)
    const double decay, //hp(1, m)
    const double varsigma, //hp(2, m)
    const double varphi, //hp(3, m)
    const double cprec, //hp(4, m)
    const arma::mat& dmat,
    const arma::vec& dvec,
    const arma::mat& dmatc
);
arma::mat covariance_matrix2(
    const double precision, //hp(0, m)
    const double decay, //hp(1, m)
    const double varsigma, //hp(2, m)
    const double varphi, //hp(3, m)
    const double cprec, //hp(4, m)
    const arma::mat& dmat,
    const arma::vec& dvec1,
    const arma::vec& dvec2,
    const arma::mat& dmatc //nxn0
);
arma::mat inv_conv_covariance_matrix(
    const double precision, 
    const double decay, 
    const double varsigma, 
    const double varphi, 
    const double cprec, 
    const arma::mat& dmat,
    const arma::mat& dmatcoast,
    const arma::mat& dr,
    const double lencoast
);
arma::mat conv_covariance_matrix(
    const double precision, //hp(0, m)
    const double decay, //hp(1, m)
    const double varsigma, //hp(2, m)
    const double varphi, //hp(3, m)
    const double cprec, //hp(4, m)
    const arma::mat& dmat,
    const arma::mat& dmatcoast,
    const arma::mat& dr,
    const double lencoast
);
arma::mat conv_covariance_matrix2(
    const double precision, //hp(0, m)
    const double decay, //hp(1, m)
    const double varsigma, //hp(2, m)
    const double varphi, //hp(3, m)
    const double cprec, //hp(4, m)
    const arma::mat& dmat,
    const arma::mat& dmatcoast1,
    const arma::mat& dmatcoast2,
    const arma::mat& dr,
    const double lencoast
);
arma::mat dist_mat(const arma::mat& A, const arma::mat& B);
arma::mat inv_mat(
    const double precision, //hp(0, m)
    const double decay, //hp(1, m)
    const double varsigma, //hp(2, m)
    const double varphi, //hp(3, m)
    const double cprec, //hp(4, m)
    const arma::mat& dmat, //coastal
    const arma::vec& dvec,
    const arma::mat& dmatc,
    const arma::mat& dmatcoast,//conv
    const arma::mat& dr,
    const double lencoast,
    int model
);

#endif