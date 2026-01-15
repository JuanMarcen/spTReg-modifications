#include "utils.h"

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

// spatial quantile model
// [[Rcpp::export]]
Rcpp::List spQuantileRcpp(
    const double tau,      
    arma::vec Y,           
    const arma::mat& X,  
    const arma::mat& V,   
    const std::vector<arma::mat>& X_alpha, 
    const arma::mat& dist, 
    const arma::vec& dist_coast, //new
    const arma::mat& dist_coast_points, //new 
    const arma::mat& dmatcoast_conv, //convolution
    const arma::mat& drmat_conv,
    const double lencoast_conv,
    const arma::vec& M,    
    const arma::mat& P,
    const std::vector<arma::vec>& M_beta_alpha,
    const std::vector<arma::mat>& P_beta_alpha,
    const double da,    
    const double db,    // add priors resto hp
    const double ga,
    const double gb,
    const double ra,
    const double rb,
    const double na,    
    const double nb,    
    arma::vec beta,        
    arma::mat alpha,   
    double prec, 
    arma::mat hp, 
    std::vector<arma::vec>& beta_alpha,
    const int N,          
    const int n,        
    const int p,
    const int r,
    const arma::vec& p_alpha, 
    const arma::uvec& s,  
    const int nSims,      
    const int nThin,
    const int nBurnin,
    const int nReport, 
    const int model
    
) {
  
  int nKeep = nSims / nThin;
  const arma::uvec missing_idx = arma::find_nonfinite(Y);
  const arma::uword missing_n = missing_idx.n_elem;
  arma::vec mm(missing_n), ss(missing_n);
  arma::mat keep_Y(nKeep, missing_n);
  const double mean_Y = arma::mean(Y.elem(arma::find_finite(Y)));
  Y.elem(missing_idx).fill(mean_Y);
  
  //fitted values
  arma::mat fittedmat(nKeep, N);
  arma::vec fitted_val(N);
  
  // loglikelihoods and DIC
  arma::vec logLik(nKeep);
  double logLikmean;
  arma::mat E_aux(nKeep, N);
  arma::vec e_aux(N);
  arma::vec vecprec(nKeep);
  double DIC;
  
  const double c1 = (1 - 2 * tau) / (tau * (1 - tau));
  const double c2 = tau * (1 - tau) / 2;
  const double c3 = 2 + c1 * c1 * c2;
  const double c4 = sqrt(c3 / c2);
  const arma::vec PM = P * M;
  std::vector<arma::vec> PM_beta_alpha(r);
  for (int m = 0; m < r; ++m) {
    PM_beta_alpha[m] = P_beta_alpha[m] * M_beta_alpha[m];
  }
  
  
  arma::vec c2dxi(N);
  arma::vec Xaux(p);
  
  arma::vec xi(N, arma::fill::ones);
  
  
  arma::vec Xb = X * beta;
  arma::vec e  = Y - Xb - c1 * xi;
  arma::vec alpha_m(n);
  arma::vec V_m(N);
  arma::uvec s_missing(missing_n);
  for (int m = 0; m < r; ++m) {
    alpha_m = alpha.col(m);
    V_m = V.col(m);
    e -= V_m % alpha_m.elem(s);
  }
  
  std::vector<arma::vec> Xb_alpha(r);
  for (int m = 0; m < r; ++m) {
    Xb_alpha[m] = X_alpha[m] * beta_alpha[m];
  }
  
  
  int Ndn = N / n;
  
  std::vector<arma::uvec> s_group(n);
  for (int i = 0; i < n; ++i) {
    s_group[i] = arma::regspace<arma::uvec>(i*Ndn, i*Ndn + Ndn - 1);
  }
  
  arma::vec V_block(Ndn), e_block(Ndn);
  arma::vec c2dxi_block(Ndn), V_c2dxi_block(Ndn);
  
  
  
  // modificado para que ean matrices que se corresponden con hp 
  arma::mat accept(5, r, arma::fill::zeros);//new
  int total = 0;
  arma::mat ratio(5, r);//new
  arma::mat sd(5, r, arma::fill::ones);//new
  arma::mat lsd(5, r, arma::fill::zeros); //new
  
  arma::cube R(n, n, r);
  arma::vec Rlogdet(r);
  std::vector<arma::mat> xR(r);
  
  //new function in utils 
  for (int m = 0; m < r; ++m){
    // R.slice(m) = inv_covariance_matrix(hp(0, m), hp(1, m), hp(2, m), hp(3, m),
    //        hp(4, m), dist, dist_coast, dist_coast_points);
    // R.slice(m) = inv_conv_covariance_matrix(hp(0, m), hp(1, m), hp(2, m), hp(3, m),
    //         hp(4, m), dist, dmatcoast_conv, drmat_conv, lencoast_conv);
    R.slice(m) = inv_mat(hp(0, m), hp(1, m), hp(2, m), hp(3, m), hp(4, m),
            dist, dist_coast, dist_coast_points, dmatcoast_conv, drmat_conv, lencoast_conv, 
            model);
    Rlogdet(m) = arma::log_det_sympd(R.slice(m));
  }
  
  
  //new aux for all hp
  double precision_aux;
  double lprecision_aux = 0;
  arma::vec lprecision = log(hp.row(0).t());
  double decay_aux;
  double ldecay_aux = 0;
  arma::vec ldecay = log(hp.row(1).t());
  double varsigma_aux;
  double lvarsigma_aux = 0;
  arma::vec lvarsigma = log(hp.row(2).t());
  double varphi_aux;
  double lvarphi_aux = 0;
  arma::vec lvarphi = log(hp.row(3).t());
  double cprec_aux;
  double lcprec_aux = 0;
  arma::vec lcprec = log(hp.row(4).t());
  
  arma::mat R_aux(n, n);
  double Rlogdet_aux;
  
  arma::vec vn(n);
  double vtRv, vtRv_aux;
  double ALPHA = 0;
  
  arma::mat Qp(p, p);
  arma::vec bp(p);
  arma::mat Qn(n, n);
  arma::vec bn(n);
  std::vector<arma::mat> Qp_alpha(r);
  std::vector<arma::vec> bp_alpha(r);
  double A = 1.5 * N + ga;
  double B;
  
  //new ajuste de numero de cols 2
  int save_idx = 0;
  int nCols1 = p + 1;
  int nCols2 = r * (n + 5) + arma::accu(p_alpha);
  arma::mat keep(nKeep, nCols1);
  arma::mat keep_alpha(nKeep, nCols2);
  
  auto start_time = std::chrono::steady_clock::now();
  
  for (int iter = 1 - nBurnin; iter <= nSims; ++iter) {
    
    reportProgress(iter, nBurnin, nSims, nReport, 
                   (iter > 0 ? iter / nThin : 0), start_time);
    
    e += c1 * xi;
    xi = 1.0 / rig(N, c4 / arma::abs(e), prec * c3);
    c2dxi = c2 / xi;
    e -= c1 * xi;
    
    //new 
    if (p > 0){
      e += Xb;
      Qp = P;
      bp = PM;
      for (int i = 0; i < N; ++i) {
        Xaux = prec * c2dxi(i) * X.row(i).t();
        Qp += Xaux * X.row(i);
        bp += Xaux * e(i);
      }
      beta = RandomMultiNormalC(Qp, bp);
      Xb = X * beta;
      e -= Xb;
    }
    
    //new R slice m no requiere de multiplicaciones extras
    //new MH para todos hp
    if (r > 0) {
      for (int m = 0; m < r; ++m) {
        alpha_m = alpha.col(m);
        V_m = V.col(m);
        e += V_m % alpha_m.elem(s);
        Qn = R.slice(m);
        bn = Qn * Xb_alpha[m];
        for (int i = 0; i < n; ++i) {
          V_block = V_m.elem(s_group[i]);
          e_block = e.elem(s_group[i]);
          c2dxi_block = c2dxi.elem(s_group[i]);
          V_c2dxi_block = V_block % c2dxi_block;
          Qn(i,i) += prec * arma::accu(V_c2dxi_block % V_block);
          bn(i)   += prec * arma::accu(V_c2dxi_block % e_block);
        }
        alpha_m = RandomMultiNormalC(Qn, bn);
        alpha.col(m) = alpha_m;
        e -= V_m % alpha_m.elem(s);
        
        xR[m] = X_alpha[m].t() * R.slice(m);
        Qp_alpha[m] = xR[m] * X_alpha[m] + P_beta_alpha[m];
        bp_alpha[m] = xR[m] * alpha_m + PM_beta_alpha[m];
        beta_alpha[m] = RandomMultiNormalC(Qp_alpha[m], bp_alpha[m]);
        Xb_alpha[m] = X_alpha[m] * beta_alpha[m]; 
        
        // decay
        ldecay_aux  = R::rnorm(ldecay(m), sd(1, m));
        decay_aux   = exp(ldecay_aux);
        // R_aux       = inv_covariance_matrix(hp(0, m), decay_aux, hp(2, m), hp(3, m),
        //                                    hp(4, m), dist, dist_coast, dist_coast_points);
        // R_aux = inv_conv_covariance_matrix(hp(0, m), decay_aux, hp(2, m), hp(3, m),
        //                 hp(4, m), dist, dmatcoast_conv, drmat_conv, lencoast_conv);
        R_aux = inv_mat(hp(0, m), decay_aux, hp(2, m), hp(3, m), hp(4, m),
                        dist, dist_coast, dist_coast_points, dmatcoast_conv, drmat_conv, lencoast_conv, 
                        model);
        Rlogdet_aux = arma::log_det_sympd(R_aux);
        vn       = alpha_m - Xb_alpha[m];
        vtRv_aux = arma::as_scalar(vn.t() * R_aux * vn);
        vtRv     = arma::as_scalar(vn.t() * R.slice(m) * vn);
        ALPHA = 
          (Rlogdet_aux - vtRv_aux) / 2 + 
          da * ldecay_aux - db * decay_aux - 
          ((Rlogdet(m) - vtRv) / 2 +
          da * ldecay(m) - db * hp(1, m));
        if (log(R::runif(0, 1)) < ALPHA) {
          ++accept(1, m);
          hp(1, m) = decay_aux;
          ldecay(m) = ldecay_aux;
          R.slice(m) = R_aux;
          Rlogdet(m) = Rlogdet_aux;
          vtRv = vtRv_aux;
        }
        
        // precision
        lprecision_aux  = R::rnorm(lprecision(m), sd(0, m));
        precision_aux   = exp(lprecision_aux);
        // R_aux       = inv_covariance_matrix(precision_aux, hp(1, m), hp(2, m), hp(3, m),
        //                                     hp(4, m), dist, dist_coast, dist_coast_points);
        // R_aux = inv_conv_covariance_matrix(precision_aux, hp(1, m), hp(2, m), hp(3, m),
        //                 hp(4, m), dist, dmatcoast_conv, drmat_conv, lencoast_conv);
        R_aux = inv_mat(precision_aux, hp(1, m), hp(2, m), hp(3, m), hp(4, m),
                        dist, dist_coast, dist_coast_points, dmatcoast_conv, drmat_conv, lencoast_conv, 
                        model);
        Rlogdet_aux = arma::log_det_sympd(R_aux);
        vn       = alpha_m - Xb_alpha[m];
        vtRv_aux = arma::as_scalar(vn.t() * R_aux * vn);
        ALPHA = 
          (Rlogdet_aux - vtRv_aux) / 2 + 
          ga * lprecision_aux - gb * precision_aux - 
          ((Rlogdet(m) - vtRv) / 2 +
          ga * lprecision(m) - gb * hp(0, m));
        if (log(R::runif(0, 1)) < ALPHA) {
          ++accept(0, m);
          hp(0, m) = precision_aux;
          lprecision(m) = lprecision_aux;
          R.slice(m) = R_aux;
          Rlogdet(m) = Rlogdet_aux;
          vtRv = vtRv_aux;
        }
        
        // varsigma
        lvarsigma_aux  = R::rnorm(lvarsigma(m), sd(2, m));
        varsigma_aux   = exp(lvarsigma_aux);
        // R_aux       = inv_covariance_matrix(hp(0, m), hp(1, m), varsigma_aux, hp(3, m),
        //                                     hp(4, m), dist, dist_coast, dist_coast_points);
        // R_aux = inv_conv_covariance_matrix(hp(0, m), hp(1, m), varsigma_aux, hp(3, m),
        //                 hp(4, m), dist, dmatcoast_conv, drmat_conv, lencoast_conv);
        R_aux = inv_mat(hp(0, m), hp(1, m), varsigma_aux, hp(3, m), hp(4, m),
                        dist, dist_coast, dist_coast_points, dmatcoast_conv, drmat_conv, lencoast_conv, 
                        model);
        Rlogdet_aux = arma::log_det_sympd(R_aux);
        vn       = alpha_m - Xb_alpha[m];
        vtRv_aux = arma::as_scalar(vn.t() * R_aux * vn);
        ALPHA = 
          (Rlogdet_aux - vtRv_aux) / 2 + 
          ga * lvarsigma_aux - gb * varsigma_aux - 
          ((Rlogdet(m) - vtRv) / 2 +
          ga * lvarsigma(m) - gb * hp(2, m));
        if (log(R::runif(0, 1)) < ALPHA) {
          ++accept(2, m);
          hp(2, m) = varsigma_aux;
          lvarsigma(m) = lvarsigma_aux;
          R.slice(m) = R_aux;
          Rlogdet(m) = Rlogdet_aux;
          vtRv = vtRv_aux;
        }
        
        // varphi
        lvarphi_aux  = R::rnorm(lvarphi(m), sd(3, m));
        varphi_aux   = exp(lvarphi_aux);
        // R_aux       = inv_covariance_matrix(hp(0, m), hp(1, m), hp(2, m), varphi_aux,
        //                                     hp(4, m), dist, dist_coast, dist_coast_points);
        // R_aux = inv_conv_covariance_matrix(hp(0, m), hp(1, m), hp(2, m), varphi_aux,
        //                 hp(4, m), dist, dmatcoast_conv, drmat_conv, lencoast_conv);
        R_aux = inv_mat(hp(0, m), hp(1, m), hp(2, m), varphi_aux, hp(4, m),
                        dist, dist_coast, dist_coast_points, dmatcoast_conv, drmat_conv, lencoast_conv, 
                        model);
        Rlogdet_aux = arma::log_det_sympd(R_aux);
        vn       = alpha_m - Xb_alpha[m];
        vtRv_aux = arma::as_scalar(vn.t() * R_aux * vn);
        ALPHA = 
          (Rlogdet_aux - vtRv_aux) / 2 + 
          ra * lvarphi_aux - rb * varphi_aux - 
          ((Rlogdet(m) - vtRv) / 2 +
          ra * lvarphi(m) - rb * hp(3, m));
        if (log(R::runif(0, 1)) < ALPHA) {
          ++accept(3, m);
          hp(3, m) = varphi_aux;
          lvarphi(m) = lvarphi_aux;
          R.slice(m) = R_aux;
          Rlogdet(m) = Rlogdet_aux;
          vtRv = vtRv_aux;
        }
        
        // cprec
        lcprec_aux  = R::rnorm(lcprec(m), sd(4, m));
        cprec_aux   = exp(lcprec_aux);
        // R_aux       = inv_covariance_matrix(hp(0, m), hp(1, m), hp(2, m), hp(3, m),
        //                                     cprec_aux, dist, dist_coast, dist_coast_points);
        // R_aux = inv_conv_covariance_matrix(hp(0, m), hp(1, m), hp(2, m), hp(3, m),
        //                 cprec_aux, dist, dmatcoast_conv, drmat_conv, lencoast_conv);
        R_aux = inv_mat(hp(0, m), hp(1, m), hp(2, m), hp(3, m), cprec_aux,
                        dist, dist_coast, dist_coast_points, dmatcoast_conv, drmat_conv, lencoast_conv, 
                        model);
        Rlogdet_aux = arma::log_det_sympd(R_aux);
        vn       = alpha_m - Xb_alpha[m];
        vtRv_aux = arma::as_scalar(vn.t() * R_aux * vn);
        ALPHA = 
          (Rlogdet_aux - vtRv_aux) / 2 + 
          ga * lcprec_aux - gb * cprec_aux - 
          ((Rlogdet(m) - vtRv) / 2 +
          ga * lcprec(m) - gb * hp(4, m));
        if (log(R::runif(0, 1)) < ALPHA) {
          ++accept(4, m);
          hp(4, m) = cprec_aux;
          lcprec(m) = lcprec_aux;
          R.slice(m) = R_aux;
          Rlogdet(m) = Rlogdet_aux;
          vtRv = vtRv_aux;
        }
      }
    }
    
    // tuneo sd en cada hp
    if (iter == 0) {
      accept.zeros();
      total = 0;
    } else if ((iter < 1) && (++total % 25 == 0)) {
      ratio = accept / total;
      for (int m = 0; m < r; ++m) {
        for (int j = 0; j < 5; ++j){
          if (ratio(j, m) > 0.33) {
            lsd(j, m) += 1 / sqrt(total / 25);
          } else {
            lsd(j, m) -= 1 / sqrt(total / 25);
          }
          sd(j, m) = exp(lsd(j, m));
        }
      }
    }
    
    B = arma::accu(xi + 0.5 * c2dxi % arma::square(e)) + gb;
    prec = R::rgamma(A, 1.0 / B);
    
    if (missing_n > 0) {
      e.elem(missing_idx) -= Y.elem(missing_idx);
      mm = X.rows(missing_idx) * beta + c1 * xi.elem(missing_idx);
      ss = arma::sqrt(xi.elem(missing_idx) / (prec * c2));
      Y.elem(missing_idx) = mm + ss % arma::randn(missing_n);
      if (r > 0) {
        for (int m = 0; m < r; ++m) {
          alpha_m = alpha.col(m);
          V_m = V.col(m);
          s_missing = s.elem(missing_idx);
          Y.elem(missing_idx) += V_m.elem(missing_idx) % alpha_m.elem(s_missing);
        }
      }
      e.elem(missing_idx) += Y.elem(missing_idx);
      
      if (iter > 0 && iter % nThin == 0) {
        keep_Y.row(save_idx) = Y.elem(missing_idx).t();
      }
    }
    
    
    
    // small modifications for good saving
    if (iter > 0 && iter % nThin == 0) {
      // predictions 
      fittedmat.row(save_idx) = (Y - e).t();
      
      //logLiks
      vecprec(save_idx) = prec;
      e_aux = e % (1 / xi);
      E_aux.row(save_idx) = e_aux.t();
      logLik(save_idx) = - 0.5 * N * log(1 / sqrt(prec)) - 0.5 *  sqrt(c2 / prec) * arma::sum(e_aux);
      
      if (p > 0){
        keep(save_idx, arma::span(0, p - 1)) = beta.t();
        keep(save_idx, p) = prec;
      }else{
        keep(save_idx, 0) = prec;
      }
      if (r > 0) {
        nCols2 = 0;
        for (int m = 0; m < r; ++m) {
          keep_alpha(save_idx, arma::span(nCols2, n - 1 + nCols2)) = alpha.col(m).t();
          keep_alpha(save_idx, arma::span(n + nCols2, n + p_alpha[m] - 1 + nCols2)) = beta_alpha[m].t();
          keep_alpha(save_idx, arma::span(n + p_alpha[m] + nCols2, n + p_alpha[m] + 4 + nCols2)) = hp.col(m).t();
          nCols2 += n + p_alpha[m] + 5;
        }
      }
      ++save_idx;
    }
  }
  
  fitted_val = arma::mean(fittedmat, 0).t();
  
  double meanprec = arma::mean(1 / sqrt(vecprec));
  logLikmean = - 0.5 * log(meanprec) - 0.5 *  sqrt(c2 / meanprec) * arma::sum(mean(E_aux, 0).t());
  
  DIC = 2 * mean(logLik) - logLikmean;
  
  if (missing_n == 0) {
    return Rcpp::List::create(
      Rcpp::Named("params") = keep,
      Rcpp::Named("process") = keep_alpha,
      Rcpp::Named("fitted") = fitted_val,
      Rcpp::Named("DIC") = DIC
    );
  } else {
    return Rcpp::List::create(
      Rcpp::Named("params") = keep,
      Rcpp::Named("process") = keep_alpha,
      Rcpp::Named("fitted") = fitted_val,
      Rcpp::Named("DIC") = DIC,
      Rcpp::Named("missing") = keep_Y
    );
  }
}

// [[Rcpp::export]]
arma::mat krigeBayesRcpp(
    const arma::mat& w,
    const arma::mat& hp, //Bx5
    const arma::mat& coords,
    const arma::mat& newcoords,
    const arma::mat& dr, //200x200
    const arma::mat& dcoast, //nx200
    const arma::mat& newdcoast, //n0x200
    const double lencoast,
    const arma::vec& dvec, 
    const arma::vec& newdvec,
    const arma::mat& dmatc,
    const arma::mat& newdmatc,
    const arma::mat& combdmatc,
    const int model
) {
  
  int B  = w.n_rows;
  int n0 = newcoords.n_rows;
  
  arma::mat out(B, n0);
  
  arma::mat d22 = dist_mat(coords, coords);
  arma::mat d11 = dist_mat(newcoords, newcoords);
  arma::mat d21 = dist_mat(coords, newcoords);
  
  int n = coords.n_rows;
  double precision, decay, varsigma, varphi, cprec;
  arma::mat R22(n,n), R11(n0,n0), R21(n,n0); 
  arma::mat C(n0,n0), L(n0,n0), R12R22inv(n0,n); 
  arma::vec wb(n), z(n0);
  
  
  for (int b = 0; b < B; ++b) {
    precision = hp(b,0);
    decay = hp(b,1);
    varsigma = hp(b,2);
    varphi = hp(b, 3);
    cprec = hp(b, 4);
    
    wb = w.row(b).t();
    
    switch(model){
    case 1:
      R22 = covariance_matrix(precision, decay, varsigma, varphi, cprec,
                              d22, dvec, dmatc);
      R11 = covariance_matrix(precision, decay, varsigma, varphi, cprec,
                              d11, newdvec, newdmatc);
      R21 = covariance_matrix2(precision, decay, varsigma, varphi, cprec,
                               d21, dvec, newdvec, combdmatc);
      break;
    case 2:
      R22 = conv_covariance_matrix(precision, decay, varsigma, varphi, cprec,
                                   d22, dcoast, dr, lencoast);
      R11 = conv_covariance_matrix(precision, decay, varsigma, varphi, cprec,
                                   d11, newdcoast, dr, lencoast);
      R21 = conv_covariance_matrix2(precision, decay, varsigma, varphi, cprec,
                                    d21, dcoast, newdcoast, dr, lencoast);
      break;
      
    }
    // R22 = conv_covariance_matrix(precision, decay, varsigma, varphi, cprec,
    //                              d22, dcoast, dr, lencoast);
    // R22 = covariance_matrix(precision, decay, varsigma, varphi, cprec,
    //                         d22, dvec, dmatc);
    // R11 = conv_covariance_matrix(precision, decay, varsigma, varphi, cprec,
    //                              d11, newdcoast, dr, lencoast);
    // R11 = covariance_matrix(precision, decay, varsigma, varphi, cprec,
    //                         d11, newdvec, newdmatc);
    // R21 = conv_covariance_matrix2(precision, decay, varsigma, varphi, cprec,
    //                              d21, dcoast, newdcoast, dr, lencoast);
    // R21 = covariance_matrix2(precision, decay, varsigma, varphi, cprec,
    //                         d21, dvec, newdvec, combdmatc);
    
    
    
    R12R22inv = arma::solve(R22, R21, arma::solve_opts::fast).t();
    C = R11 - R12R22inv * R21;
    L = arma::chol(C, "lower");
    
    z = arma::randn(n0);
    
    out.row(b) = (R12R22inv * wb + (L * z)).t();
  }
  return out;
  
}