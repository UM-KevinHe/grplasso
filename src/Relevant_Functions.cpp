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

arma::vec rep1(arma::vec &x, arma::vec &each) {
  arma::vec x_rep(arma::sum(each));
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

arma::mat info_beta_omp(const arma::mat &Z, const arma::vec &pq) {
  int threads = omp_get_max_threads();
  omp_set_num_threads(threads);
  unsigned int p = Z.n_cols;
  unsigned int loops = p * (1 + p) / 2;
  arma::mat output(p, p);
  #pragma omp parallel for schedule(static)
  for (unsigned int i = 0; i < loops; i++) {
    unsigned int r, c;
    ind2uppsub(i, p, r, c);
    output(r,c) = dot(Z.col(r), Z.col(c)%pq);
    output(c,r) = output(r,c);
  }
  return(output);
}

double Loglkd_1(const arma::vec &Y, const arma::vec &Z_beta, const arma::vec &gamma_obs) {
  return sum((gamma_obs+Z_beta)%Y-log(1+exp(gamma_obs+Z_beta)));
}


// SerBIN used for finding residuals
// [[Rcpp::export]]
List SerBIN(arma::vec &Y, arma::mat &Z, arma::vec &n_prov, arma::vec gamma, arma::vec beta) {
  int iter = 0, n = Z.n_rows, m = n_prov.n_elem, ind;
  arma::vec gamma_obs(n);
  double crit = 100.0;
  double loglkd = Loglkd_1(Y, Z * beta, rep1(gamma, n_prov)), d_loglkd, v, lambda, s = 0.01, t = 0.6;
  arma::vec gamma_obs_tmp(n), gamma_tmp(m), beta_tmp(Z.n_cols);
  while (iter < 10000) {
    R_CheckUserInterrupt();
    if (crit < 1e-5) {
      break;
    }
    iter++;
    gamma_obs = rep1(gamma, n_prov);
    arma::vec Z_beta = Z * beta;
    arma::vec p = 1 / (1 + exp(-gamma_obs-Z_beta));
    arma::vec Yp = Y - p, pq = p % (1-p);
    arma::vec score_gamma(m), info_gamma_inv(m);
    arma::mat info_betagamma(Z.n_cols,m);
    ind = 0;
    for (int i = 0; i < m; i++) {
      score_gamma(i) = sum(Yp(span(ind,ind+n_prov(i)-1)));
      info_gamma_inv(i) = 1 / sum(pq(span(ind,ind+n_prov(i)-1)));
      info_betagamma.col(i) =
        sum(Z.rows(ind,ind+n_prov(i)-1).each_col()%(p.subvec(ind,ind+n_prov(i)-1)%(1-p.subvec(ind,ind+n_prov(i)-1)))).t();
      ind += n_prov(i);
    }
    arma::vec score_beta = Z.t() * Yp;
    arma::mat info_beta(Z.n_cols, Z.n_cols);
    info_beta = info_beta_omp(Z, pq);
    arma::mat mat_tmp1 = trans(info_betagamma.each_row()%info_gamma_inv.t());
    arma::mat schur_inv = inv_sympd(info_beta-mat_tmp1.t()*info_betagamma.t());
    arma::mat mat_tmp2 = mat_tmp1*schur_inv;
    arma::vec d_gamma = info_gamma_inv%score_gamma + mat_tmp2*(mat_tmp1.t()*score_gamma-score_beta);
    arma::vec d_beta = schur_inv*score_beta - mat_tmp2.t()*score_gamma;
    v = 1.0;
    gamma_tmp = gamma + v * d_gamma;
    gamma_obs_tmp = rep1(gamma_tmp, n_prov);
    arma::vec Z_beta_tmp = Z * (beta+v*d_beta);
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
  List ret = List::create(_["gamma"] = gamma, _["beta"] = beta);
  return ret;
}

tuple<arma::vec, arma::vec, arma::vec, arma::mat, arma::mat> Update_logit(arma::vec t, arma::mat X, arma::vec delta_obs, arma::vec gamma, arma::vec beta, int max_t, int c, int r){
  arma::vec score_gamma(max_t, fill::zeros);
  arma::vec score_beta(c, fill::zeros);
  arma::vec info_gamma(max_t, fill::zeros);
  arma::mat info_beta(c, c, fill::zeros);
  arma::mat info_betagamma(max_t, c, fill::zeros);

  for (int i = 0 ; i < r ; i++){
    for (int s = 1 ; s <= t(i) ; s++){
        double lambda = 1/(1 + exp(-gamma(s-1) - dot(X.row(i), beta)));
        score_gamma(s-1) = score_gamma(s-1) - lambda;
        score_beta = score_beta - lambda * X.row(i).t();
        info_gamma(s-1) = info_gamma(s-1) + lambda * (1 - lambda);
        info_beta = info_beta + (lambda * (1 - lambda)) * (X.row(i).t() * X.row(i));
        info_betagamma.row(s-1) = info_betagamma.row(s-1) + (lambda * (1-lambda)) * X.row(i);
        if (t(i) == s && delta_obs(i) == 1){
          score_gamma(s-1) = score_gamma(s-1) + 1;
          score_beta = score_beta + X.row(i).t();
        }
    }
  }

  return make_tuple(score_gamma, score_beta, info_gamma, info_beta, info_betagamma);
}

// [[Rcpp::export]]
List NR_residuals(arma::vec t, arma::mat X, arma::vec delta_obs, arma::vec gamma, arma::vec beta, double tol, int max_iter){
  int r = X.n_rows;
  int c = X.n_cols;
  int max_t = t.max();

  arma::vec gamma2 = gamma;
  arma::vec beta2 = beta;
  arma::vec beta_change(c);
  arma::mat A_inv(max_t, max_t, fill::zeros);
  arma::vec A(max_t, fill::zeros);
  arma::mat B(max_t, c, fill::zeros);
  arma::mat C(c, c, fill::zeros);
  arma::mat schur(c, c, fill::zeros);
  arma::vec score_gamma(max_t, fill::zeros);
  arma::vec score_beta(c, fill::zeros);

  auto update = Update_logit(t, X, delta_obs, gamma, beta, max_t, c, r);
  int iter = 0;

  for (int i = 0 ; i <= max_iter; i++){
    iter++;
    R_CheckUserInterrupt();
    double MaxChange_beta = 0;
    tie(score_gamma, score_beta, A, C, B) = update;
    A = 1/A;
    A_inv = diagmat(A);
    schur = (C - (B.t()) * (A_inv) * (B)).i();
    gamma2 = gamma + (A_inv + A_inv * B * schur * (B.t()) * A_inv) * score_gamma - A_inv * B * schur * score_beta;
    beta2 = beta + schur * score_beta-schur * (B.t())* A_inv * score_gamma;
    beta_change = beta2 - beta;
    gamma = gamma2;
    beta = beta2;

    for (i = 0; i < c; i++){
      if (fabs(beta_change(i)) > MaxChange_beta) {
          MaxChange_beta = fabs(beta_change(i));
      }
    }
    if(MaxChange_beta < tol){
      break;
    }
    update = Update_logit(t, X, delta_obs, gamma, beta, max_t, c, r);
  }

  List result = List::create(_["beta"] = beta, _["gamma"] = gamma, _["iter"] = iter);
  return result;
}
