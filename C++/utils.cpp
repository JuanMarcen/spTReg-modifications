#include "utils.h"

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

void reportProgress(
    int iter,
    int nBurnin,
    int nSims,
    int nReport,
    int save_idx,
    std::chrono::steady_clock::time_point start_time) {
  
  int iter_from_start = iter - (1 - nBurnin);
  int total_iters = nBurnin + nSims;
  
  if (nReport > 0 && (iter_from_start % nReport) == 0) {
    
    auto now = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - start_time).count() / 1000.0;
    
    double pct = 100.0 * double(iter_from_start) / double(std::max(1, total_iters));
    
    double est_total = elapsed / std::max(0.001, pct / 100.0); // total estimado
    double remaining = std::max(0.0, est_total - elapsed);
    
    if (iter <= 0) {
      int burned = iter + nBurnin;
      Rcpp::Rcout << "Burn-in: iteration " << burned << " / " << nBurnin
                  << " (" << std::floor(pct) << "%)"
                  << " | elapsed: " << std::round(elapsed) << "s"
                  << " | ETA: " << std::round(remaining) << "s"
                  << "\n";
    } else {
      Rcpp::Rcout << "Sampling: iteration " << iter << " / " << nSims
                  << ", saved: " << save_idx
                  << " (" << std::floor(pct) << "%)"
                  << " | elapsed: " << std::round(elapsed) << "s"
                  << " | ETA: " << std::round(remaining) << "s"
                  << "\n";
    }
  }
}

// [[Rcpp::export]]
arma::vec RandomMultiNormalC(
    const arma::mat& Q, 
    const arma::vec& b) {
  
  int n = Q.n_rows;
  arma::vec x(n);
  
  arma::mat L = arma::chol(Q, "lower");
  arma::vec y = arma::solve(arma::trimatl(L), b);
  arma::vec mu = arma::solve(arma::trimatu(L.t()), y);
  arma::vec z = arma::randn(n);
  arma::vec u = arma::solve(arma::trimatu(L.t()), z);
  x = mu + u;
  
  return x;
}
// [[Rcpp::export]]
arma::vec rig(
    const int N,
    const arma::vec& mu,
    const double lambda) {
  
  arma::vec X(N);
  arma::vec V = arma::square(arma::randn(N)); // chi^2(1)
  arma::vec U = arma::randu(N);               // uniform
  
  arma::vec C = mu / (2.0 * lambda);
  double mu_i, X_i, W, P1;
  
  for (int i = 0; i < N; ++i) {
    mu_i = mu(i);
    W = mu_i * V(i);
    X_i = mu_i + C(i) * (W - std::sqrt(W * (4.0 * lambda + W)));
    P1 = mu_i / (mu_i + X_i);
    if (U(i) > P1)
      X_i = mu_i * mu_i / X_i;
    X(i) = X_i;
  }
  
  return X;
}


// [[Rcpp::export]]
arma::mat inv_covariance_matrix(
    const double precision, //hp(0, m)
    const double decay, //hp(1, m)
    const double varsigma, //hp(2, m)
    const double varphi, //hp(3, m)
    const double cprec, //hp(4, m)
    const arma::mat& dmat,
    const arma::vec& dvec,
    const arma::mat& dmatc
){
  // new definition of the covariance matrices
  arma::vec expdc = exp(- varsigma* dvec); //nx1
  arma::mat Mcoast = exp(- varphi * dmatc); //nxn
  // column multiplication
  Mcoast.each_row() %= expdc.t();
  // row multiplication
  Mcoast.each_col() %= expdc;
  
  arma::mat R = arma::inv_sympd(1.0 / precision * exp(- decay * dmat) + 1.0 / cprec * Mcoast);
  
  return R;
}

// [[Rcpp::export]]
arma::mat covariance_matrix(
    const double precision, //hp(0, m)
    const double decay, //hp(1, m)
    const double varsigma, //hp(2, m)
    const double varphi, //hp(3, m)
    const double cprec, //hp(4, m)
    const arma::mat& dmat,
    const arma::vec& dvec,
    const arma::mat& dmatc
){
  // new definition of the covariance matrices
  arma::vec expdc = exp(- varsigma* dvec); //nx1
  arma::mat Mcoast = exp(- varphi * dmatc); //nxn
  // column multiplication
  Mcoast.each_row() %= expdc.t();
  // row multiplication
  Mcoast.each_col() %= expdc;
  
  arma::mat R = 1.0 / precision * exp(- decay * dmat) + 1.0 / cprec * Mcoast;
  
  return R;
}


// [[Rcpp::export]]
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
){
  // new definition of the covariance matrices
  arma::vec expdc1 = exp(- varsigma* dvec1); //nx1
  arma::vec expdc2 = exp(- varsigma* dvec2); //n0x1
  arma::mat Mcoast = exp(- varphi * dmatc); //nxn0
  // column multiplication
  Mcoast.each_row() %= expdc2.t();
  // row multiplication
  Mcoast.each_col() %= expdc1;
  
  arma::mat R = 1.0 / precision * exp(- decay * dmat) + 1.0 / cprec * Mcoast;
  
  return R;
}


// [[Rcpp::export]]
arma::mat inv_conv_covariance_matrix(
    const double precision, //hp(0, m)
    const double decay, //hp(1, m)
    const double varsigma, //hp(2, m)
    const double varphi, //hp(3, m)
    const double cprec, //hp(4, m)
    const arma::mat& dmat,
    const arma::mat& dmatcoast,
    const arma::mat& dr,
    const double lencoast
){
  arma::mat expdc = exp(- varsigma * dmatcoast); //nxM
  arma::mat covcoast = 1.0 / cprec * exp(- varphi * dr); //MxM
  
  double factor = lencoast * covcoast.n_rows;
  double factor2 = factor * factor;
  // no multiplicar por factor
  arma::mat Mcoast =  expdc * covcoast * expdc.t(); //nxn
  
  arma::mat R = arma::inv_sympd(1.0 / precision * exp(- decay * dmat) + Mcoast);
  
  return R;
}

// [[Rcpp::export]]
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
){
  arma::mat expdc = exp(- varsigma * dmatcoast); //nxM
  arma::mat covcoast = 1.0 / cprec * exp(- varphi * dr); //MxM
  
  double factor = lencoast * covcoast.n_rows;
  double factor2 = factor * factor;
  // no multiplicar por factor
  arma::mat Mcoast = expdc * covcoast * expdc.t(); //nxn
  
  arma::mat R = 1.0 / precision * exp(- decay * dmat) + Mcoast;
  
  return R;
}

// [[Rcpp::export]]
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
){
  arma::mat expdc1 = exp(- varsigma * dmatcoast1); //nxM
  arma::mat expdc2 = exp(- varsigma * dmatcoast2);
  arma::mat covcoast = 1.0 / cprec * exp(- varphi * dr); //MxM
  
  double factor = lencoast * covcoast.n_rows;
  double factor2 = factor * factor;
  
  // no multiplicar por factor
  arma::mat Mcoast = expdc1 * covcoast * expdc2.t(); //nxn
  
  arma::mat R = 1.0 / precision * exp(- decay * dmat) + Mcoast;
  
  return R;
}

// [[Rcpp::export]]
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
){
  if (model == 1) return(inv_covariance_matrix(precision, decay, varsigma, varphi, cprec,
                         dmat, dvec, dmatc));
  if (model == 2) return(inv_conv_covariance_matrix(precision, decay, varsigma, varphi, cprec, 
                         dmat, dmatcoast, dr, lencoast));
}

// [[Rcpp::export]]
arma::mat dist_mat(const arma::mat& A, const arma::mat& B) {
  int n1 = A.n_rows;
  int n2 = B.n_rows;
  arma::mat D(n1, n2);
  
  for (int i = 0; i < n1; ++i) {
    for (int j = 0; j < n2; ++j) {
      D(i,j) = arma::norm(A.row(i) - B.row(j));
    }
  }
  return D;
}
