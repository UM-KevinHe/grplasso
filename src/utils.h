#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <cmath>
#include <omp.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
using namespace std;
using namespace arma;
using namespace Rcpp;

// ---------- Vector/Matrix helpers ----------
arma::vec rep(arma::vec &x, arma::vec &each);
double vec_crossprod(arma::vec &w, arma::vec &r);
double mean_crossprod(arma::mat &Z, arma::vec &r, int j, int n_obs);
arma::vec inner_product(arma::mat &Z);
double weighted_inner_product(arma::mat &Z, arma::vec &w, int j);
double crossprod(arma::mat &Z, arma::vec &r, int j);
double w_crossprod(arma::mat &Z, arma::vec &r, arma::vec &w, int j);

// ---------- Probability / Link ----------
double p_binomial(double &eta);
double p_binomial_Surv(double &gamma, double &eta);

// ---------- Soft-thresholding ----------
double Soft_thres(double z, double l);

// ---------- Log-likelihood ----------
double Loglkd(arma::vec &Y, arma::vec &eta);
double Loglkd_Surv(int n_obs, arma::vec &delta_obs, arma::vec &time, arma::vec &gamma, arma::vec &eta);

// *** MODIFIED: added FirthPenalty and FirthGamma declarations ***
double FirthPenalty(arma::vec &p, arma::vec &n_prov, arma::vec &ind_start, int n_gamma);
arma::vec FirthGamma(arma::vec &Y_bar, arma::vec &n_prov, arma::vec &gamma_init, int max_iter, double tol);
// *** END MODIFIED ***

// ---------- Deviance ----------
double Deviance(arma::vec &Y, arma::vec &p);

// ---------- Upper-triangular index conversion (for parallel info matrix) ----------
void ind2uppsub(unsigned int index, unsigned int dim, unsigned int &row, unsigned int &col);
arma::mat info_beta_omp(const arma::mat &Z, const arma::vec &pq);

#endif // UTILS_H
