#include "utils.h"

// ============================================================
//  Vector / Matrix Utility Functions
// ============================================================

arma::vec rep(arma::vec &x, arma::vec &each) {
  arma::vec x_rep(arma::sum(each));
  int ind = 0, m = x.n_elem;
  for (int i = 0; i < m; i++) {
    x_rep.subvec(ind, ind + each(i) - 1) = x(i) * ones(each(i));
    ind += each(i);
  }
  return x_rep;
}

double vec_crossprod(arma::vec &w, arma::vec &r) {
  double vec_crossprod = dot(w, r);
  return(vec_crossprod);
}

double mean_crossprod(arma::mat &Z, arma::vec &r, int j, int n_obs) {
  double crossprod = dot(Z.col(j), r);
  return(crossprod / n_obs);
}

arma::vec inner_product(arma::mat &Z) {
  int p = Z.n_cols;
  arma::vec ip(p);
  for (int i = 0; i < p; i++) {
    ip(i) = dot(Z.col(i), Z.col(i));
  }
  return(ip);
}

double weighted_inner_product(arma::mat &Z, arma::vec &w, int j) {
  int n = Z.n_rows;
  double weighted_ip = 0;
  for (int i = 0; i < n; i++) {
    weighted_ip += w(i) * Z(i, j) * Z(i, j);
  }
  return(weighted_ip);
}

double crossprod(arma::mat &Z, arma::vec &r, int j) {
  double crossprod = dot(Z.col(j), r);
  return(crossprod);
}

double w_crossprod(arma::mat &Z, arma::vec &r, arma::vec &w, int j) {
  int n = Z.n_rows;
  double w_crossprod = 0;
  for (int i = 0; i < n; i++) {
    w_crossprod += w(i) * Z(i, j) * r(i);
  }
  return(w_crossprod);
}


// ============================================================
//  Probability / Link Functions
// ============================================================

double p_binomial(double &eta) {
  return(1 / (1 + exp(-eta)));
}

double p_binomial_Surv(double &gamma, double &eta) {
  return(1 / (1 + exp(-gamma - eta)));
}


// ============================================================
//  Soft-thresholding Operator
// ============================================================

double Soft_thres(double z, double l) {
  if (z > l) {
    return(z - l);
  } else if (z < -l) {
    return(z + l);
  } else {
    return(0);
  }
}


// ============================================================
//  Log-likelihood Functions
// ============================================================

double Loglkd(arma::vec &Y, arma::vec &eta) {
  return sum((eta) % Y - log(1 + exp(eta)));
}

double Loglkd_Surv(int n_obs, arma::vec &delta_obs, arma::vec &time, arma::vec &gamma, arma::vec &eta) {
  double loglkd = 0;
  for (int j = 0; j < n_obs; j++) {
    for (int i = 0; i < time(j); i++) {
      loglkd -= log(1 + exp(gamma(i) + eta(j)));
      if (i == (time(j) - 1) && delta_obs(j) == 1) {
        loglkd += gamma(i) + eta(j);
      }
    }
  }
  return(loglkd);
}

// *** MODIFIED: added FirthPenalty function ***
// Computes 0.5 * sum_i log(W_i) where W_i = sum_j mu_ij*(1-mu_ij)
double FirthPenalty(arma::vec &p, arma::vec &n_prov, arma::vec &ind_start, int n_gamma) {
  double penalty = 0.0;
  for (int i = 0; i < n_gamma; i++) {
    arma::vec p_i = p(span(ind_start(i), ind_start(i) + n_prov(i) - 1));
    double W_i = sum(p_i % (1 - p_i));
    if (W_i > 0) {
      penalty += log(W_i);
    }
  }
  return 0.5 * penalty;
}

// Solves Firth-corrected gamma when beta = 0 via Newton iteration.
// When beta = 0, all mu_ij within cluster i are equal to mu_i = expit(gamma_i),
// so the score simplifies to: n_i*(Ybar_i - mu_i) + 0.5*(1 - 2*mu_i) = 0,
// and the information simplifies to: n_i*w_i - 0.5*(d_i^2 - 2*w_i).
// [[Rcpp::export]]
arma::vec FirthGamma(arma::vec &Y_bar, arma::vec &n_prov, arma::vec &gamma_init, int max_iter, double tol) {
  int m = gamma_init.n_elem;
  arma::vec gamma = gamma_init;
  for (int iter = 0; iter < max_iter; iter++) {
    double max_step = 0.0;
    for (int i = 0; i < m; i++) {
      double mu = 1.0 / (1.0 + exp(-gamma(i)));
      double w = mu * (1.0 - mu);
      double d = 1.0 - 2.0 * mu;
      double U = n_prov(i) * (Y_bar(i) - mu) + 0.5 * d;
      double J = n_prov(i) * w - 0.5 * (d * d - 2.0 * w);
      double step = U / J;
      gamma(i) += step;
      if (fabs(step) > max_step) max_step = fabs(step);
    }
    if (max_step < tol) break;
  }
  return gamma;
}
// *** END MODIFIED ***


// ============================================================
//  Deviance & Max Lambda (Rcpp exported)
// ============================================================

// [[Rcpp::export]]
double Deviance(arma::vec &Y, arma::vec &p) {
  double Dev = 0;
  int n_obs = Y.n_elem;
  for (int i = 0; i < n_obs; i++) {
    if (p(i) != 0 && p(i) != 1) {
      Dev -= 2 * Y(i) * log(p(i)) + 2 * (1 - Y(i)) * log(1 - p(i));
    }
  }
  return(Dev);
}

// [[Rcpp::export]]
double Z_max_grLasso(arma::mat &x, arma::vec &r, arma::vec &K, arma::vec &m) {
  int J = K.n_elem - 1;
  double z_max = 0, z;
  for (int g = 0; g < J; g++) {
    int Kg = K(g + 1) - K(g);
    arma::vec Z(Kg);
    for (int j = K(g); j < K(g + 1); j++) {
      arma::vec x_tmp = x.col(j);
      Z(j - K(g)) = dot(x_tmp, r);
    }
    z = arma::norm(Z) / m(g);
    if (z > z_max) {
      z_max = z;
    }
  }
  return(z_max);
}


// ============================================================
//  Upper-triangular Index & Parallel Information Matrix
// ============================================================

void ind2uppsub(unsigned int index, unsigned int dim, unsigned int &row, unsigned int &col) {
  row = 0; col = dim - 1;
  unsigned int n = dim * (dim - 1) / 2 - (dim - row) * (dim - row - 1) / 2 + col;
  while (index > n) {
    ++row;
    n = dim * (dim - 1) / 2 - (dim - row) * (dim - row - 1) / 2 + col;
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
    output(r, c) = dot(Z.col(r), Z.col(c) % pq);
    output(c, r) = output(r, c);
  }
  return(output);
}
