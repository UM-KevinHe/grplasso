#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include <cmath>

using namespace std;
using namespace arma;
using namespace Rcpp;

vec rep(vec &x,vec &each) { 
  vec x_rep(arma::sum(each));
  int ind = 0, m = x.n_elem; 
  for (int i = 0; i < m; i++) {
    x_rep.subvec(ind,ind+each(i)-1) = x(i) * ones(each(i));
    ind += each(i);
  }
  return x_rep;
}

double Loglkd(const vec &Y, const vec &Z_beta, const vec &gamma_obs) {
  return sum((gamma_obs+Z_beta)%Y-log(1+exp(gamma_obs+Z_beta)));
}

vec norm_matrix(mat x) { 
  int p = x.n_cols;
  vec norm(p);
  for (int i = 0; i < p; i++) {
    norm(i) = dot(x.col(i), x.col(i));
  }
  return(norm);
}

double crossprod(mat x, vec y, int j) {
  double crossprod = 0;
  crossprod = dot(x.col(j), y);
  return(crossprod);
}

double w_crossprod(mat x, vec y, vec w, int j) {
  int n = x.n_rows;
  double w_crossprod = 0;
  for (int i = 0; i < n; i++){
    w_crossprod += w(i) * x(i,j) * y(i);
  }
  return(w_crossprod);
}

double weighted_norm_matrix(mat x, vec w, int j) {
  int n = x.n_rows;
  double weighted_norm = 0;
  for (int i = 0; i < n; i++){
    weighted_norm += w(i) * x(i,j) * x(i,j);
  }
  return(weighted_norm);
}

// [[Rcpp::export]]
double Soft_thres(double z, double l) {
  if (z > l) {
    return(z - l);
  } else if (z < -l) {
    return(z + l);
  } else {
    return(0);
  }
}

// [[Rcpp::export]]
List pp_Lasso(vec &Y, mat &Z, vec &n_prov, vec &gamma, vec &beta, bool backtrack = true, bool MM = false,
               int max_iter = 10000, double bound = 10.0, double lam1 = 0, double tol = 1e-5) {
  int iter = 0, n_obs = Z.n_rows, n_beta = Z.n_cols, m = n_prov.n_elem, ind;
  vec gamma_obs(n_obs), old_beta(n_beta);
  vec beta_change(n_beta);  
  vec p(n_obs), Yp(n_obs), pq(n_obs);
  uvec update_order;
  vec r(n_obs);
  double crit = 100.0; 
  cout << "Implementing pp-LASSO (with Rcpp) ..." << endl;
  
  if (backtrack == true){ 
    double loglkd, d_loglkd, v, lambda, s = 0.01, t = 0.8;
    vec gamma_obs_tmp(n_obs), gamma_tmp(m);
    while (iter < max_iter) {
      R_CheckUserInterrupt();
      if (crit < tol) {
        break;
      }
      old_beta = beta;
      iter++;
      gamma_obs = rep(gamma, n_prov);
      vec Z_beta = Z * beta;
      p = 1 / (1 + exp(- gamma_obs - Z_beta));
      Yp = Y - p;
      pq = p % (1 - p);
      vec score_gamma(m), d_gamma(m);
      ind = 0;
      for (int i = 0; i < m; i++) {
        score_gamma(i) = sum(Yp(span(ind,ind+n_prov(i)-1))); 
        d_gamma(i) = score_gamma(i) / sum(pq(span(ind,ind+n_prov(i)-1)));  
        ind += n_prov(i);  
      }
      v = 1.0; 
      loglkd = Loglkd(Y, Z_beta, gamma_obs);  
      gamma_tmp = gamma + v * d_gamma; 
      gamma_obs_tmp = rep(gamma_tmp, n_prov);  
      d_loglkd = Loglkd(Y, Z_beta, gamma_obs_tmp) - loglkd;  
      lambda = dot(score_gamma, d_gamma);  
      while (d_loglkd < s*v*lambda) { 
        v = t*v;
        gamma_tmp = gamma + v * d_gamma;
        gamma_obs_tmp = rep(gamma_tmp, n_prov); 
        d_loglkd = Loglkd(Y, Z_beta, gamma_obs_tmp) - loglkd;
      }
      gamma += v * d_gamma;
      gamma = clamp(gamma, median(gamma)-bound, median(gamma)+bound);
      gamma_obs = rep(gamma, n_prov);
           
      update_order = randperm(n_beta);
      if (MM == true) {
        p = 1 / (1 + exp(- gamma_obs - Z_beta)); 
        r = (Y - p)/0.25; 
        vec innprod_Z = norm_matrix(Z); 
        for (int j = 0; j < n_beta; j++) { 
          double prod_Zr = crossprod(Z, r, update_order(j)); 
          double beta_initial = prod_Zr/innprod_Z(update_order(j)) + old_beta(update_order(j)); 
          double threshold = lam1/(2 * 0.25 * innprod_Z(update_order(j))); 
          beta(update_order(j)) = Soft_thres(beta_initial, threshold);
          beta_change(update_order(j)) = beta(update_order(j)) - old_beta(update_order(j)); 
          if (beta_change(update_order(j)) != 0){ 
            vec r_shift = Z.col(update_order(j)) * beta_change(update_order(j));  
            r -= r_shift;
          }
        }
      } else { 
        vec eta = gamma_obs + Z_beta;
        p = 1 / (1 + exp(- eta));
        r = (Y - p)/(p % (1 - p));
        for (int j = 0; j < n_beta; j++) {
          p = 1 / (1 + exp(- eta)); 
          vec w = p % (1 - p); 
          double weighted_innprod_Z = weighted_norm_matrix(Z, w, update_order(j)); 
          double weighted_prod_Zr = w_crossprod(Z, r, w, update_order(j));
          double beta_initial = weighted_prod_Zr/weighted_innprod_Z + old_beta(update_order(j));
          double threshold = lam1/(2 * weighted_innprod_Z);
          beta(update_order(j)) = Soft_thres(beta_initial, threshold);
          beta_change(update_order(j)) = beta(update_order(j)) - old_beta(update_order(j));
          if (beta_change(update_order(j)) != 0){
            vec r_shift = Z.col(update_order(j)) * beta_change(update_order(j));
            r -= r_shift;
            eta += r_shift;
          }
        }
      }
      crit = norm(beta_change, "inf");
      cout << "Iter " << iter << ": Maximum change of Beta is " << setprecision(3) << scientific << crit << ";" << endl;
    }
  } else {
    while (iter < max_iter) {
      R_CheckUserInterrupt();
      if (crit < tol) {
        break;
      }
      old_beta = beta;
      iter++;
      gamma_obs = rep(gamma, n_prov);
      vec Z_beta = Z * beta;
      p = 1 / (1 + exp(- gamma_obs - Z_beta));
      Yp = Y - p;
      pq = p % (1 - p);
      vec score_gamma(m), d_gamma(m);
      ind = 0;
      for (int i = 0; i < m; i++) {
        score_gamma(i) = sum(Yp(span(ind,ind+n_prov(i)-1))); 
        d_gamma(i) = score_gamma(i) / sum(pq(span(ind,ind+n_prov(i)-1)));
        gamma(i) += d_gamma(i);
        ind += n_prov(i);  
      }
      gamma = clamp(gamma, median(gamma)-bound, median(gamma)+bound);
      gamma_obs = rep(gamma, n_prov);

      // update beta
      update_order = randperm(n_beta);
      if (MM == true) {
        p = 1 / (1 + exp(- gamma_obs - Z_beta));
        r = (Y - p)/0.25; 
        vec innprod_Z = norm_matrix(Z);
        for (int j = 0; j < n_beta; j++) {
          double prod_Zr = crossprod(Z, r, update_order(j));
          double beta_initial = prod_Zr/innprod_Z(update_order(j)) + old_beta(update_order(j));
          double threshold = lam1/(2 * 0.25 * innprod_Z(update_order(j)));
          beta(update_order(j)) = Soft_thres(beta_initial, threshold);
          beta_change(update_order(j)) = beta(update_order(j)) - old_beta(update_order(j));
          if (beta_change(update_order(j)) != 0){
            vec r_shift = Z.col(update_order(j)) * beta_change(update_order(j));
            r -= r_shift;
          }
        }
      } else { 
        vec eta = gamma_obs + Z_beta; 
        p = 1 / (1 + exp(- eta));
        r = (Y - p)/(p % (1 - p)); 
        for (int j = 0; j < n_beta; j++) {
          p = 1 / (1 + exp(- eta));
          vec w = p % (1 - p);
          double weighted_innprod_Z = weighted_norm_matrix(Z, w, update_order(j)); 
          double weighted_prod_Zr = w_crossprod(Z, r, w, update_order(j));
          double beta_initial = weighted_prod_Zr/weighted_innprod_Z + old_beta(update_order(j));
          double threshold = lam1/(2 * weighted_innprod_Z);
          beta(update_order(j)) = Soft_thres(beta_initial, threshold);
          beta_change(update_order(j)) = beta(j) - old_beta(update_order(j));
          if (beta_change(update_order(j)) != 0){
            vec r_shift = Z.col(update_order(j)) * beta_change(update_order(j));
            r -= r_shift;
            eta += r_shift;
          }
        }
      }
      crit = norm(beta_change, "inf");
      cout << "Iter " << iter << ": Maximum change of Beta is " << setprecision(3) << scientific << crit << ";" << endl;
    }
  }
  cout << "pp-Lasso (with Rcpp) converged after " << iter << " iterations!" << endl;
  List ret = List::create(_["gamma"] = gamma, _["beta"] = beta); //Export beta and gamma
  return ret;
}

// [[Rcpp::export]]
List logis_fe_prov(vec &Y, mat &Z, vec &n_prov, vec gamma, vec beta, int backtrack=0,
                   int max_iter=10000, double bound=10.0, double tol=1e-5) {
  
  int iter = 0, n = Z.n_rows, m = n_prov.n_elem, ind;
  vec gamma_obs(n);
  double crit = 100.0; 
  cout << "Implementing Newton-Raphson algorithm (Rcpp) for fixed provider effects model ..." << endl;
  if (backtrack==1) {
    double loglkd, d_loglkd, v, lambda, s = 0.01, t = 0.6;
    vec gamma_obs_tmp(n), gamma_tmp(m), beta_tmp(Z.n_cols);
    while (iter < max_iter) {
      if (crit < tol) {
        break;
      }
      iter++;
      gamma_obs = rep(gamma, n_prov);
      vec Z_beta = Z * beta;
      vec p = 1 / (1 + exp(-gamma_obs-Z_beta));
      vec Yp = Y - p, pq = p % (1-p);
      vec score_gamma(m), d_gamma(m);
      ind = 0;
      for (int i = 0; i < m; i++) {
        score_gamma(i) = sum(Yp(span(ind,ind+n_prov(i)-1)));
        d_gamma(i) = score_gamma(i) / sum(pq(span(ind,ind+n_prov(i)-1)));
        ind += n_prov(i);
      }
      v = 1.0; 
      loglkd = Loglkd(Y, Z_beta, gamma_obs);
      gamma_tmp = gamma + v * d_gamma;
      gamma_obs_tmp = rep(gamma_tmp, n_prov);
      d_loglkd = Loglkd(Y, Z_beta, gamma_obs_tmp) - loglkd;
      lambda = dot(score_gamma, d_gamma);
      while (d_loglkd < s*v*lambda) {
        v = t*v;
        gamma_tmp = gamma + v * d_gamma;
        gamma_obs_tmp = rep(gamma_tmp, n_prov);
        d_loglkd = Loglkd(Y, Z_beta, gamma_obs_tmp) - loglkd;
      }
      gamma += v * d_gamma; 
      gamma = clamp(gamma, median(gamma)-bound, median(gamma)+bound); 
      gamma_obs = rep(gamma, n_prov);
      
      p = 1/(1+exp(-gamma_obs-Z_beta));
      pq = p % (1-p);
      vec score_beta = Z.t() * (Y-p);
      mat info_beta = Z.t() * (Z.each_col()%pq); 
      vec d_beta = solve(info_beta, score_beta, solve_opts::fast+solve_opts::likely_sympd);
      v = 1.0;
      loglkd = Loglkd(Y, Z_beta, gamma_obs);
      beta_tmp = beta + v * d_beta;
      d_loglkd = Loglkd(Y, Z*beta_tmp, gamma_obs) - loglkd;
      lambda = dot(score_beta, d_beta);
      while (d_loglkd < s*v*lambda) {
        v = t * v;
        beta_tmp = beta + v * d_beta;
        d_loglkd = Loglkd(Y, Z*beta_tmp, gamma_obs) - loglkd;
      }
      beta += v * d_beta;
      //t4 = clock();
      crit = norm(v*d_beta, "inf");
      cout << "Iter " << iter << ": Inf norm of running diff in est reg parm is " << setprecision(3) << scientific << crit << ";" << endl;
    }
  } else if (backtrack==0) {
    while (iter < max_iter) {
      if (crit < tol) {
        break;
      }
      iter++;
      gamma_obs = rep(gamma, n_prov);
      vec Z_beta = Z * beta;
      vec p = 1 / (1 + exp(-gamma_obs-Z_beta));
      vec Yp = Y - p, pq = p % (1-p);
      ind = 0;
      for (int i = 0; i < m; i++) {
        gamma(i) += sum(Yp(span(ind,ind+n_prov(i)-1))) / 
          sum(pq(span(ind,ind+n_prov(i)-1)));
        ind += n_prov(i);
      }
      gamma = clamp(gamma, median(gamma)-bound, median(gamma)+bound);
      gamma_obs = rep(gamma, n_prov);
      p = 1/(1+exp(-gamma_obs-Z_beta));
      pq = p % (1-p);
      vec score_beta = Z.t() * Yp;
      mat info_beta = Z.t() * (Z.each_col()%pq);
      vec d_beta = solve(info_beta, score_beta, solve_opts::fast+solve_opts::likely_sympd); 
      beta += d_beta;
      crit = norm(d_beta, "inf");
      cout << "Iter " << iter << ": Inf norm of running diff in est reg parm is " << setprecision(3) << scientific << crit << ";" << endl;
    }
  }
  
  cout << "Newton-Raphson algorithm (Rcpp) converged after " << iter << " iterations!" << endl;
  List ret = List::create(_["gamma"]=gamma, _["beta"]=beta); 
  return ret;
}
