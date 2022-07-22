#define Check_Headers
#include <RcppArmadillo.h>  
// [[Rcpp::depends(RcppArmadillo)]]
#include <cmath>
#include <omp.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
using namespace std;
using namespace arma;
using namespace Rcpp;

vec rep1(vec &x, vec &each) { 
  vec x_rep(arma::sum(each));
  int ind = 0, m = x.n_elem; 
  for (int i = 0; i < m; i++) {
    x_rep.subvec(ind,ind+each(i)-1) = x(i) * ones(each(i));
    ind += each(i);
  }
  return x_rep;
}

void ind2uppsub(unsigned int index, unsigned int dim, unsigned int &row, unsigned int &col) {
  row = 0, col = dim-1;
  unsigned int n = dim*(dim-1)/2 - (dim-row)*(dim-row-1)/2 + col;
  while (index > n) {
    ++row;
    n = dim*(dim-1)/2 - (dim-row)*(dim-row-1)/2 + col;
  }
  while (index < n) {
    --col;
    --n;
  }
}

mat info_beta_omp(const mat &Z, const vec &pq) {
  int threads = omp_get_max_threads();
  omp_set_num_threads(threads);
  unsigned int p = Z.n_cols;
  unsigned int loops = p * (1 + p) / 2;
  mat output(p, p);
  #pragma omp parallel for schedule(static)
  for (unsigned int i = 0; i < loops; i++) {
    unsigned int r, c;
    ind2uppsub(i, p, r, c);
    output(r,c) = dot(Z.col(r), Z.col(c)%pq);
    output(c,r) = output(r,c);
  }
  return(output);
}

double Loglkd_1(const vec &Y, const vec &Z_beta, const vec &gamma_obs) {
  return sum((gamma_obs+Z_beta)%Y-log(1+exp(gamma_obs+Z_beta)));
}


// SerBIN used for finding residuals
// [[Rcpp::export]]
List SerBIN(vec &Y, mat &Z, vec &n_prov, vec gamma, vec beta) {
  int iter = 0, n = Z.n_rows, m = n_prov.n_elem, ind;
  vec gamma_obs(n);
  double crit = 100.0; 
  double loglkd = Loglkd_1(Y, Z * beta, rep1(gamma, n_prov)), d_loglkd, v, lambda, s = 0.01, t = 0.6;
  vec gamma_obs_tmp(n), gamma_tmp(m), beta_tmp(Z.n_cols);
  while (iter < 10000) {
    R_CheckUserInterrupt();
    if (crit < 1e-5) {
      break;
    }
    iter++;
    gamma_obs = rep1(gamma, n_prov);
    vec Z_beta = Z * beta;
    vec p = 1 / (1 + exp(-gamma_obs-Z_beta));
    vec Yp = Y - p, pq = p % (1-p);
    vec score_gamma(m), info_gamma_inv(m);
    mat info_betagamma(Z.n_cols,m);
    ind = 0;
    for (int i = 0; i < m; i++) {
      score_gamma(i) = sum(Yp(span(ind,ind+n_prov(i)-1)));
      info_gamma_inv(i) = 1 / sum(pq(span(ind,ind+n_prov(i)-1)));
      info_betagamma.col(i) = 
        sum(Z.rows(ind,ind+n_prov(i)-1).each_col()%(p.subvec(ind,ind+n_prov(i)-1)%(1-p.subvec(ind,ind+n_prov(i)-1)))).t();
      ind += n_prov(i);
    }
    vec score_beta = Z.t() * Yp;
    mat info_beta(Z.n_cols, Z.n_cols);
    info_beta = info_beta_omp(Z, pq); 
    mat mat_tmp1 = trans(info_betagamma.each_row()%info_gamma_inv.t()); 
    mat schur_inv = inv_sympd(info_beta-mat_tmp1.t()*info_betagamma.t());
    mat mat_tmp2 = mat_tmp1*schur_inv;
    vec d_gamma = info_gamma_inv%score_gamma + mat_tmp2*(mat_tmp1.t()*score_gamma-score_beta);
    vec d_beta = schur_inv*score_beta - mat_tmp2.t()*score_gamma;
    v = 1.0; 
    gamma_tmp = gamma + v * d_gamma;
    gamma_obs_tmp = rep1(gamma_tmp, n_prov);
    vec Z_beta_tmp = Z * (beta+v*d_beta);
    d_loglkd = Loglkd_1(Y, Z_beta_tmp, gamma_obs_tmp) - loglkd;
    lambda = dot(score_gamma, d_gamma) + dot(score_beta, d_beta);
    while (d_loglkd < s*v*lambda) {
      v = t*v;
      gamma_tmp = gamma + v * d_gamma;
      gamma_obs_tmp = rep1(gamma_tmp, n_prov);
      Z_beta_tmp = Z * (beta+v*d_beta);
      d_loglkd = Loglkd_1(Y, Z_beta_tmp, gamma_obs_tmp) - loglkd;
    }
    gamma += v * d_gamma;
    gamma = clamp(gamma, median(gamma) - 10, median(gamma) + 10);
    beta += v * d_beta;
    loglkd += d_loglkd;
    crit = norm(v*d_beta, "inf");
  }
  List ret = List::create(_["gamma"]=gamma, _["beta"]=beta);
  return ret;
}