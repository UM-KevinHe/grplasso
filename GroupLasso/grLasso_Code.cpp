#define Check_Headers
#include <RcppArmadillo.h>  
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
#include <cmath>
#include <omp.h>
#include <chrono>


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
using namespace Rcpp;
using namespace std;
using namespace arma;
using std::chrono::high_resolution_clock;
using std::chrono::duration;
using std::chrono::milliseconds;
using std::chrono::seconds;

vec rep(vec &x,vec &each) { 
  vec x_rep(arma::sum(each));
  int ind = 0, m = x.n_elem; 
  for (int i = 0; i < m; i++) {
    x_rep.subvec(ind,ind+each(i)-1) = x(i) * ones(each(i));
    ind += each(i);
  }
  return x_rep;
}

double mean_crossprod(mat &Z, vec &Y, int j, int n_obs) {
  //double crossprod = 0;
  //for (int i = 0; i < n_obs; i++){
  //  crossprod += Z(i, j) * Y(i);
  //}
  double crossprod = dot(Z.col(j), Y);
  return(crossprod/n_obs);
}

double Loglkd(vec &Y, vec &eta) {
  return sum((eta) % Y - log(1 + exp(eta)));
}

// Pr(y=1) for binomial
double p_binomial(double &eta) {
  if (eta > 10) {
    return(1);
  } else if (eta < -10) {
    return(0);
  } else {
    return(1/(1+exp(-eta)));
  }
}

//soft-thresholding operator
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
double Z_max_grLasso(mat &x, vec &y, vec &K, vec &m){ // "K": a vector contains the start index of each group; "m": group.multiplier
  int J = K.n_elem - 1;  //number of penalized group
  double z_max = 0, z;
  for (int g = 0; g < J; g++){
    int Kg = K(g + 1) - K(g); //number of features in group g
    vec Z(Kg);
    for (int j = K(g); j < K(g + 1); j++) {
      vec x_tmp = x.col(j);
      Z(j - K(g)) = dot(x_tmp, y);
    }
    z = arma::norm(Z)/m(g);
    if (z > z_max){
      z_max = z;
    }
  }
  return(z_max);
}

// [[Rcpp::export]]
double Deviance(vec &Y, vec &p){
  double Dev = 0;
  int n_obs = Y.n_elem;
  for (int i = 0; i < n_obs; i++){
    if (p(i) != 0 && p(i) != 1){
      Dev -= 2 * Y(i) * log(p(i)) + 2 * (1 - Y(i)) * log(1 - p(i));
    }
  }
  return(Dev);
}

// update penalized beta
void gd_glm(vec &beta, mat &Z, vec &r, vec &eta, vec &old_beta, int g, vec &K1, int n_obs, int n_beta, double lambda, double &df, double &MaxChange_beta){
  int K = K1(g + 1) - K1(g); //number of features in group g
  vec beta_initial(K);
  for (int j = K1(g); j < K1(g + 1); j++){
    beta_initial(j - K1(g)) = mean_crossprod(Z, r, j, n_obs) + old_beta(j);
  }
  double beta_initial_norm = arma::norm(beta_initial, 2); // ||Z_g|| in grpreg paper
  double len = Soft_thres(beta_initial_norm, lambda/0.25); 

  if (len != 0 || old_beta(K1(g)) != 0){ // old_beta(K1(g)) != 0: It's enough to consider only one beta in this group
    for (int j = K1(g); j < K1(g + 1); j++){
      beta(j) = len * beta_initial(j - K1(g)) / beta_initial_norm; 
      double beta_change = beta(j) - old_beta(j);
      if (fabs(beta_change) > MaxChange_beta) {
        MaxChange_beta = fabs(beta_change);
      }
      r -= beta_change * Z.col(j);
      eta += beta_change * Z.col(j);
    }
  }
  //update df: since beta has been shrunk, we cannot directly add K
  if (len > 0){
    df += K * len / beta_initial_norm;
  }
}


// estimataion for one given lambda
tuple<vec, vec, vec, double, double, int> grp_lasso_fit(vec &Y, mat &Z, vec &n_prov, vec gamma, vec beta, vec eta, int K0, vec &K1,  
                                                        double lambda, int &tol_iter, int max_iter, vec &group_multiplier, int max_n_prov, 
                                                        bool backtrack, double bound, double tol, vec &ind_start, vec &active_group, 
                                                        int n_beta, int n_gamma, int n_obs, int n_group, bool single_intercept, int threads){
  vec old_beta, old_gamma, p(n_obs), r(n_obs), r_shift, Z_tmp;
  double Dev, df, MaxChange_beta, shift, lambda_g, v = 0.25;
  int iter = 0, violations;

  while (tol_iter < max_iter) {
    while (tol_iter < max_iter) {
      R_CheckUserInterrupt();
      old_gamma = gamma;
      old_beta = beta;
      df = 0;
      tol_iter++;
      iter++;
      MaxChange_beta = 0;
      for (int i = 0; i < n_obs; i++){
        p(i) = p_binomial(eta(i)); 
      }
      
      //auto t1 = high_resolution_clock::now();
      if (single_intercept == true) {  // use coordinate descent to update intercept
        r = (Y - p)/v; //initial pseudo residual vector
        shift = sum(r)/n_obs;
        gamma += shift;
        eta += shift;
        r -= shift;
      } else {  //use Newton method to update gamma
        double omega_min = 1e-10, info_gamma;
        int nProcessors = threads;
        vec score_gamma(n_gamma), d_gamma(n_gamma), Yp(n_obs), pq(n_obs), gamma_obs(n_obs);
        if (backtrack == true) {
          vec gamma_shift_tmp, eta_tmp(n_obs);
          double loglkd, d_loglkd, u, k, s = 0.01, t = 0.8;
          Yp = Y - p;  //Score
          pq = p % (1 - p); //Info
          omp_set_num_threads(nProcessors);
          #pragma omp parallel for schedule(static)
          for (int i = 0; i < n_gamma; i++) {
            score_gamma(i) = sum(Yp(span(ind_start(i), ind_start(i) + n_prov(i) - 1)));
            info_gamma = sum(pq(span(ind_start(i), ind_start(i) + n_prov(i) - 1)));
            info_gamma = std::max(omega_min, std::min(info_gamma, 0.25 * max_n_prov));
            d_gamma(i) = score_gamma(i) / info_gamma;  
          }
          u = 1.0; //initial step size "u"
          loglkd = Loglkd(Y, eta);  
          gamma_shift_tmp = u * d_gamma;
          eta_tmp = eta + rep(gamma_shift_tmp, n_prov);
          d_loglkd = Loglkd(Y, eta_tmp) - loglkd;  
          k = dot(score_gamma, d_gamma);  
          while (d_loglkd < s * u * k) {
            u = t * u;
            gamma_shift_tmp = u * d_gamma;
            eta_tmp = eta + rep(gamma_shift_tmp, n_prov);
            d_loglkd = Loglkd(Y, eta_tmp) - loglkd;
          }
          gamma = gamma + u * d_gamma;
          gamma = clamp(gamma, median(gamma) - bound, median(gamma) + bound);
          vec gamma_shift = gamma - old_gamma;
          eta += rep(gamma_shift, n_prov);
          for (int i = 0; i < n_obs; i++){
            p(i) = p_binomial(eta(i)); 
          }
          r = (Y - p)/v;  //initial pseudo residual vector
        } else { 
          Yp = Y - p;
          pq = p % (1 - p);
          omp_set_num_threads(nProcessors);
          #pragma omp parallel for schedule(static)
          for (int i = 0; i < n_gamma; i++) {
            score_gamma(i) = sum(Yp(span(ind_start(i), ind_start(i) + n_prov(i) - 1)));
            info_gamma = sum(pq(span(ind_start(i), ind_start(i) + n_prov(i) - 1)));
            info_gamma = std::max(omega_min, std::min(info_gamma, 0.25 * max_n_prov));
           d_gamma(i) = score_gamma(i) / info_gamma;  
          }
          gamma = gamma + d_gamma;
          gamma = clamp(gamma, median(gamma) - bound, median(gamma) + bound);
          vec gamma_shift = gamma - old_gamma;
          eta += rep(gamma_shift, n_prov);
          for (int i = 0; i < n_obs; i++){
            p(i) = p_binomial(eta(i)); 
          }
          r = (Y - p)/v;  //initial pseudo residual vector
        }
      }
      //auto t2 = high_resolution_clock::now();
      //duration<double, std::milli> ms_double1 = t2 - t1; 
      //cout << "<#> Update Gamma:" << ms_double1.count() << "ms\n";

      //auto t3 = high_resolution_clock::now();
      uvec update_order_unpenalized = randperm(K0);  // Randomized coordinate descent for unpenalized covariate
      for (int j = 0; j < K0; j++){  // If K0 is zero, then the whole iteration will be skipped
        shift = mean_crossprod(Z, r, update_order_unpenalized(j), n_obs);
        if ((fabs)(shift) > MaxChange_beta) {
          MaxChange_beta = fabs(shift);
        }
        beta(update_order_unpenalized(j)) = old_beta(update_order_unpenalized(j)) + shift; 
        r -= Z.col(update_order_unpenalized(j)) * shift;
        eta += Z.col(update_order_unpenalized(j)) * shift;  
        df++;
      }
      //auto t4 = high_resolution_clock::now();
      //duration<double, std::milli> ms_double2 = t4 - t3; 
      //cout << "<<##>> Update Unpenalized Beta:" << ms_double2.count() << "ms\n";

      //auto t5 = high_resolution_clock::now();
      for (int g = 0; g < n_group; g++){
        if (active_group(g) == 1){
          lambda_g = lambda * group_multiplier(g);
          gd_glm(beta, Z, r, eta, old_beta, g, K1, n_obs, n_beta, lambda_g, df, MaxChange_beta);
        }
      }
      //auto t6 = high_resolution_clock::now();
      //duration<double, std::milli> ms_double3 = t6 - t5; 
      //cout << "<<<###>>> Update Penalized Beta:" << ms_double3.count() << "ms\n";

      if (MaxChange_beta < tol){
        break;
      }
    }
    //outer loop: update active set
    violations = 0;
    for (int g = 0; g < n_group; g++){
      if (active_group(g) == 0){ //only check non-active groups
        lambda_g = lambda * group_multiplier(g);
        gd_glm(beta, Z, r, eta, old_beta, g, K1, n_obs, n_beta, lambda_g, df, MaxChange_beta);
        if (beta(K1(g)) != 0){ //add group to active set
            active_group(g) = 1;
            violations++;
        }
      }
    }
    if (violations == 0){
      break;
    }
  }
  for (int i = 0; i < n_obs; i++){
    p(i) = p_binomial(eta(i)); 
  }
  Dev = Deviance(Y, p);
  return make_tuple(beta, gamma, eta, Dev, df, iter);
}



// [[Rcpp::export]]
List grp_lasso(vec &Y, mat &Z, vec &n_prov, vec &gamma, vec &beta, int K0, vec &K1, vec &lambda_seq, bool lambda_early_stop, 
               double stop_dev_ratio, vec &group_multiplier, int max_iter, double tol, double nullDev, bool backtrack, 
               double bound, int initial_active_group, double nvar_max, double group_max, bool trace_lambda, 
               bool single_intercept, int threads) {
  int n_obs = Z.n_rows, n_beta = Z.n_cols, n_gamma = n_prov.n_elem, n_lambda = lambda_seq.n_elem, max_n_prov = max(n_prov);;
  int n_group = K1.n_elem - 1; // n_group: number of penalized group
  int tol_iter = 0;
  
  mat beta_matrix(n_beta, n_lambda, fill::zeros), gamma_matrix(n_gamma, n_lambda, fill::zeros); // parameter estimation
  mat eta_matrix(n_obs, n_lambda, fill::zeros); // linear predictor
  vec Dev_vec(n_lambda, fill::zeros);  // "Deviance" is defined as "loss"
  vec iter_vec(n_lambda, fill::zeros);  // number of iterations for each lambda
  vec df_vec(n_lambda, fill::zeros); // number of non-zero beta

  vec active_group(n_group, fill::zeros); 
  if (K0 == 0){ // if there's no unpenalized beta, we need set the smallest group as the first active group
    active_group(initial_active_group) = 1;
  }

  vec ind_start(n_gamma); // index of the first observation of each provider
  ind_start(0) = 0;
  for (int i = 1; i < n_gamma; i++) {
    ind_start(i) = ind_start(i - 1) + n_prov(i - 1);
  }

  // initialize eta
  vec eta = rep(gamma, n_prov) + Z * beta;

  for (int l = 0; l < n_lambda; l++){
    R_CheckUserInterrupt();
    if (trace_lambda == true){
      cout << "processing lambda: " << l + 1 << " (total: " << l + 1 << "/" << n_lambda << ")..." << endl;
    }
    double lambda = lambda_seq(l);
    auto fit = grp_lasso_fit(Y, Z, n_prov, gamma, beta, eta, K0, K1, lambda, tol_iter, max_iter, group_multiplier, max_n_prov, backtrack, bound, tol, ind_start, active_group, n_beta, n_gamma, n_obs, n_group, single_intercept, threads);
    double Dev_l, df_l;
    int iter_l;
    tie(beta, gamma, eta, Dev_l, df_l, iter_l) = fit;
    beta_matrix.col(l) = beta;
    gamma_matrix.col(l) = gamma;
    eta_matrix.col(l) = eta;
    Dev_vec(l) = Dev_l;
    df_vec(l) = df_l;
    iter_vec(l) = iter_l;

    // check dfmax, gmax; (note, "nv" doesn't equal "df")
    int ng = 0, nv = 0;
    for (int g = 0; g < n_group; g++){
      if (beta(K1(g)) != 0){
         ng++;
         nv += (K1(g + 1) - K1(g));
      }
    }
    if (ng > group_max || nv > nvar_max || tol_iter == max_iter) {
      if (tol_iter == max_iter) {
        cout << "Reach maximum number of iterations, algorithm stops..." << endl;
      } else {
        cout << "Reach maximum number of groups or variables, algorithm stops..." << endl;
      }
      
      for (int ll = (l + 1); ll < n_lambda; ll++){
        iter_vec(ll) = NA_REAL;
      }
      break;
    }

    if (lambda_early_stop == true){
      // early stop, when deviance difference is small between two lambda (start from second lambda)
      if (l != 0){
        double Dev_ratio = (Dev_vec(l) - Dev_vec(l - 1))/(Dev_vec(l) - nullDev);
        if (Dev_ratio < stop_dev_ratio){
          for (int ll = (l + 1); ll < n_lambda; ll++){
            iter_vec(ll) = NA_REAL;
          }
        break;
        }
      }
    }
  }
  
  List result = List::create(_["gamma"] = gamma_matrix, _["beta"] = beta_matrix, _["Deviance"] = Dev_vec, _["Eta"] = eta_matrix, _["Df"] = df_vec, _["iter"] = iter_vec);
  return result;
}