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

vec rep(vec &x,vec &each) { 
  vec x_rep(arma::sum(each));
  int ind = 0, m = x.n_elem; 
  for (int i = 0; i < m; i++) {
    x_rep.subvec(ind,ind+each(i)-1) = x(i) * ones(each(i));
    ind += each(i);
  }
  return x_rep;
}

double Loglkd_Surv(int n_obs, vec &delta_obs, vec &time, vec &gamma, vec &eta) {
  double loglkd = 0;
  for (int j = 0; j < n_obs; j++){  //j: individual
    for (int i = 0; i < time(j); i++){ //i: time point
      loglkd -= log(1 + exp(gamma(i) + eta(j)));
      if (i == (time(j) - 1) && delta_obs(j) == 1) { //if jth individual failed at final time point, then likelihood should add "gamma(i) + eta(j)"
        loglkd += gamma(i) + eta(j);
      }
    }
  }
  return(loglkd);
}

double weighted_inner_product(mat &Z, vec &w, int j) {
  int n = Z.n_rows;
  double weighted_ip = 0;
  for (int i = 0; i < n; i++){
    weighted_ip += w(i) * Z(i,j) * Z(i,j);
  }
  return(weighted_ip);
}

// Weighted cross product of r with jth column of Z
double w_crossprod(mat &Z, vec &r, vec &w, int j) {
  int n = Z.n_rows;
  double w_crossprod = 0;
  for (int i = 0; i < n; i++){
    w_crossprod += w(i) * Z(i,j) * r(i);
  }
  return(w_crossprod);
}


// where eta = X\beta
double p_binomial_Surv(double &gamma, double &eta) {
  return(1/(1+exp(-gamma - eta)));
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
double Z_max_grLasso(mat &x, vec &r, vec &K, vec &m){ // "K": a vector contains the start index of each group; "m": group.multiplier
  int J = K.n_elem - 1;  //number of penalized group
  double z_max = 0, z;
  for (int g = 0; g < J; g++){
    int Kg = K(g + 1) - K(g); //number of features in group g
    vec Z(Kg);
    for (int j = K(g); j < K(g + 1); j++) {
      vec x_tmp = x.col(j);
      Z(j - K(g)) = dot(x_tmp, r);
    }
    z = arma::norm(Z)/m(g);
    if (z > z_max){
      z_max = z;
    }
  }
  return(z_max);
}

// [[Rcpp::export]]
vec DiscSurv_residuals(int n_obs, vec &delta_obs, vec &time, vec &gamma, vec &eta){
  vec residuals(n_obs);
  for (int j = 0; j < n_obs; j++){  //j: individual
    for (int i = 0; i < time(j); i++){ //i: time point
      residuals(j) += p_binomial_Surv(gamma(i), eta(j));
    }
  } 
  residuals = -residuals + delta_obs;
  return(residuals);
}

// [[Rcpp::export]]
mat predict_linear_predictor(int n_lambda, int n_obs, int expand_n_obs, vec &time, mat &gamma, mat &eta) {
  mat predict_linear_predictor(expand_n_obs, n_lambda);
  for (int l = 0; l < n_lambda; l++){
    int sum_index = 0;
    for (int j = 0; j < n_obs; j++){  //j: individual
      for (int i = 0; i < time(j); i++){ //i: time point
        predict_linear_predictor(sum_index, l) = p_binomial_Surv(gamma(i, l), eta(j, l));
        sum_index += 1;
      }
    } 
  }
  return(predict_linear_predictor);
}

// Function1: with logit-link
// update penalized beta
void gd_Surv(vec &beta, mat &Z, vec &r, vec &eta, vec &time, vec &old_beta, vec &w, int p, int n_obs, double lambda, double &df, double &MaxChange_beta){
  double beta_initial = w_crossprod(Z, r, w, p)/weighted_inner_product(Z, w, p) + old_beta(p);   
  double threshold = n_obs * lambda/weighted_inner_product(Z, w, p);
  double len = Soft_thres(beta_initial, threshold); 
  
  if (len != old_beta(p)){ // if beta has changed, then update eta and r
    beta(p) = len;
    double beta_change = beta(p) - old_beta(p);
    if (fabs(beta_change) > MaxChange_beta) {
      MaxChange_beta = fabs(beta_change);
    }
    if (beta_change != 0){
      r -= beta_change * Z.col(p);
      eta += beta_change * Z.col(p);
    }
  }
  if (len != 0){
    df += len / beta_initial;
  }
}


double gd_Surv_BetaChange(mat &Z, vec &r, vec &eta, vec &time, vec &w, int p, int n_obs, double lambda){
  double beta_initial = w_crossprod(Z, r, w, p)/weighted_inner_product(Z, w, p);   
  double threshold = n_obs * lambda/weighted_inner_product(Z, w, p);
  double len = Soft_thres(beta_initial, threshold); 
  
  if (len != 0){ 
    return(fabs(len)); 
  } else {
    return(0);
  }
}


tuple<vec, vec, vec, vec, double, int> pp_DiscSurv_fit(vec &delta_obs, int max_timepoint, mat &Z, vec &time, vec &n_prov, vec gamma, vec beta, vec alpha, vec eta,
                                                       int K0, vec &K1, vec &sum_failure, vec failure_each_center, double lambda, int &tol_iter, int max_total_iter, 
                                                       int max_each_iter, vec &penalized_multiplier, bool backtrack, bool MM, double bound, double tol, vec &ind_start, 
                                                       vec &active_var, int n_obs, int n_var, int threads, bool actSet, int actIter, int activeVarNum, bool actSetRemove){
  vec old_beta = beta, old_alpha = alpha, r(n_obs), r_shift, w(n_obs);
  int n_alpha = alpha.n_elem;
  int nProcessors = threads; 
  double df, MaxChange_beta, shift;
  double v = 0.25, omega_min = 1e-10;
  double p_gamma;
  int iter = 0; //"iter" counts the number of iterations for each lambda
  
  while (tol_iter < max_total_iter && iter < max_each_iter) { //"tol_iter" counts the number of iterations for the entire lambda sequence
    int inner_loop_iter = 0; // count the number of iterations for a new updated active set
    R_CheckUserInterrupt();
    // inner loop: update variables in the current active set
    while (tol_iter < max_total_iter && iter < max_each_iter && inner_loop_iter < actIter) {

      // the maximum number of inner iterations is "actIter", and after that we will update the current active set;
      R_CheckUserInterrupt();
      df = 0;
      iter++;
      inner_loop_iter++;
      MaxChange_beta = 0;

      // 1. update gamma (time effect)
      // initial score and fisher information matrix of gamma
      // "sum_failure": a vector contains the total number of failures at each time point (donnot need to consider which center one obs belongs to)
      vec score_gamma = - sum_failure, info_gamma(max_timepoint, fill::zeros), d_gamma(max_timepoint, fill::zeros); 
      if (backtrack == true) {
        vec gamma_tmp, eta_tmp(n_obs);
        double loglkd, d_loglkd, u = 1.0, k, s = 0.01, t = 0.8;
        for (int j = 0; j < n_obs; j++){ //iterated over everyone in the dataset
          for (int i = 0; i < time(j); i++){  //iterated over time
            p_gamma = p_binomial_Surv(gamma(i), eta(j));
            score_gamma(i) += p_gamma;
            info_gamma(i) -= p_gamma * (1 - p_gamma); //a diagonal matrix
          }
        }
        omp_set_num_threads(nProcessors);
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < max_timepoint; i++){
          double tmp_info_gamma = info_gamma(i);
          d_gamma(i) = score_gamma(i)/std::min(-omega_min, tmp_info_gamma);
        }
        loglkd = Loglkd_Surv(n_obs, delta_obs, time, gamma, eta); //eta will not change when only update time effect
        gamma_tmp = u * d_gamma + gamma;
        d_loglkd = Loglkd_Surv(n_obs, delta_obs, time, gamma_tmp, eta) - loglkd;  
        k = dot(score_gamma, d_gamma);  
        while (d_loglkd < s * u * k) {
          u = t * u;
          gamma_tmp = u * d_gamma + gamma;
          d_loglkd = Loglkd_Surv(n_obs, delta_obs, time, gamma_tmp, eta) - loglkd;
        }
        gamma = gamma + u * d_gamma;
        // for some timepoint with all failure or no failure (only have cencor), we should bound gamma, or it will become +Inf or -Inf (with baseline hazard = 0/1)
        gamma = clamp(gamma, median(gamma) - bound, median(gamma) + bound); 
      } else {
        for (int j = 0; j < n_obs; j++){ 
          for (int i = 0; i < time(j); i++){
            p_gamma = p_binomial_Surv(gamma(i), eta(j));
            score_gamma(i) += p_gamma;
            info_gamma(i) -= p_gamma * (1 - p_gamma);
          }
        }    

        int nProcessors = threads;
        omp_set_num_threads(nProcessors);
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < max_timepoint; i++){
          double tmp_info_gamma = info_gamma(i);
          d_gamma(i) = score_gamma(i)/std::min(-omega_min, tmp_info_gamma);
        }
        gamma = gamma + d_gamma;
        gamma = clamp(gamma, median(gamma) - bound, median(gamma) + bound); 
      } 
 
      // 2. update alpha
      double temp_p;
      vec p_obs(n_obs, fill::zeros); //each element of p_obs sum over all time points for that observation
      vec pq_obs(n_obs, fill::zeros);
      for (int j = 0; j < n_obs; j++){  //j: individual
        for (int i = 0; i < time(j); i++){ //i: time point
          temp_p = p_binomial_Surv(gamma(i), eta(j));
          p_obs(j) += temp_p;
          pq_obs(j) += temp_p * (1 - temp_p);
        }
      } 
      
      vec score_alpha(n_alpha, fill::zeros), info_alpha(n_alpha, fill::zeros), d_alpha(n_alpha, fill::zeros); 
      if (backtrack == true){
        vec alpha_shift_tmp, eta_tmp(n_obs);
        double loglkd, d_loglkd, u2 = 1.0, k2, s2 = 0.01, t2 = 0.8;
        omp_set_num_threads(nProcessors);
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < n_alpha; i++) {
          score_alpha(i) = sum(p_obs(span(ind_start(i), ind_start(i) + n_prov(i) - 1))); 
          double tmp_info_alpha = -sum(pq_obs(span(ind_start(i), ind_start(i) + n_prov(i) - 1)));
          info_alpha(i) = std::min(-omega_min, tmp_info_alpha);
        }
        score_alpha -= failure_each_center;
        for (int i = 0; i < n_alpha; i++) {
          d_alpha(i) = score_alpha(i)/info_alpha(i);
        }
        loglkd = Loglkd_Surv(n_obs, delta_obs, time, gamma, eta);
        alpha_shift_tmp = u2 * d_alpha;
        eta_tmp = eta + rep(alpha_shift_tmp, n_prov);
        d_loglkd = Loglkd_Surv(n_obs, delta_obs, time, gamma, eta_tmp) - loglkd;
        k2 = dot(score_alpha, d_alpha);
        while (d_loglkd < s2 * u2 * k2) {
          u2 = t2 * u2;
          alpha_shift_tmp = u2 * d_alpha;
          eta_tmp = eta + rep(alpha_shift_tmp, n_prov);
          d_loglkd = Loglkd_Surv(n_obs, delta_obs, time, gamma, eta_tmp) - loglkd;
        }
        alpha = alpha + u2 * d_alpha;
        alpha = clamp(alpha, median(alpha) - bound, median(alpha) + bound);
        vec alpha_shift = alpha - old_alpha;
        eta += rep(alpha_shift, n_prov);
      } else {
        omp_set_num_threads(nProcessors);
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < n_alpha; i++) {
          score_alpha(i) = sum(p_obs(span(ind_start(i), ind_start(i) + n_prov(i) - 1)));
          double tmp_info_alpha = -sum(pq_obs(span(ind_start(i), ind_start(i) + n_prov(i) - 1)));
          info_alpha(i) = std::min(-omega_min, tmp_info_alpha);
        } 
        score_alpha -= failure_each_center;
        for (int i = 0; i < n_alpha; i++) {
          d_alpha(i) = score_alpha(i)/info_alpha(i);
        }
        alpha = alpha + d_alpha;
        alpha = clamp(alpha, median(alpha) - bound, median(alpha) + bound);

        vec alpha_shift = alpha - old_alpha;
        eta += rep(alpha_shift, n_prov);
      }
        
      // 3. update beta
      // p_eta: 把每个人所有时间点的exp/1+exp加起来;
      // w: 把每个人所有时间点的exp/(1+exp)^2加起来;
      vec p_eta(n_obs, fill::zeros);
      
      //initial pseudo residual vector
      if (MM == true){ // surrogate function use W = 0.25 * diag{m1, m2, ..., mn}, which is a n*n matrix
        for (int j = 0; j < n_obs; j++){  //j: individual
          for (int i = 0; i < time(j); i++){ //i: time point
            p_eta(j) += p_binomial_Surv(gamma(i), eta(j));
          }
          w(j) = v * time(j);
        } 
      } else {
        w.zeros();
        for (int j = 0; j < n_obs; j++){  
          for (int i = 0; i < time(j); i++){ 
            double p_temp = p_binomial_Surv(gamma(i), eta(j));
            p_eta(j) += p_temp;
            w(j) += p_temp * (1 - p_temp);
          }
        } 
        if (any(w == 0)) {
          w.replace(0, omega_min); 
        }
      }
      vec score_eta = - delta_obs + p_eta;
      r = - score_eta / w; //initial pseudo outcome
           
      // 3.1 update unpenalized beta
      uvec update_order_unpenalized = randperm(K0);  // Randomized coordinate descent for unpenalized covariate 

      for (int p = 0; p < K0; p++){  // only need w (we don't need compute "p_eta")
        if (MM == true){
          for (int j = 0; j < n_obs; j++){
            w(j) = v * time(j);
          } 
        } else {
          w.zeros();
          for (int j = 0; j < n_obs; j++){  
            for (int i = 0; i < time(j); i++){ 
              double p_temp = p_binomial_Surv(gamma(i), eta(j));
              w(j) += p_temp * (1 - p_temp);
            }
          }
        }
        shift = w_crossprod(Z, r, w, update_order_unpenalized(p))/weighted_inner_product(Z, w, update_order_unpenalized(p)); 
        if (fabs(shift) > MaxChange_beta) {
          MaxChange_beta = fabs(shift);
        }
        beta(update_order_unpenalized(p)) = old_beta(update_order_unpenalized(p)) + shift; 
        r -= Z.col(update_order_unpenalized(p)) * shift;
        eta += Z.col(update_order_unpenalized(p)) * shift; 
        df++;
      }
      
      // 3.2 update penalized beta
      uvec update_order_penalized = randperm(n_var);
      for (int p = 0; p < n_var; p++) {
        int update_column_index = K1(update_order_penalized(p));
        if (active_var(update_order_penalized(p)) == 1){
          double lambda_m = lambda * penalized_multiplier(update_order_penalized(p));
          if (MM == true){
            for (int j = 0; j < n_obs; j++){
              w(j) = v * time(j);
            } 
          } else {
            w.zeros();
            for (int j = 0; j < n_obs; j++){  
              for (int i = 0; i < time(j); i++){ 
                double p_temp = p_binomial_Surv(gamma(i), eta(j));
                w(j) += p_temp * (1 - p_temp);
              }
            }
          }      
          gd_Surv(beta, Z, r, eta, time, old_beta, w, update_column_index, n_obs, lambda_m, df, MaxChange_beta);
        }
      }

      old_beta = beta;
      old_alpha = alpha;
      if (MaxChange_beta < tol){ 
        break;
      }
    }
    
    if (actSet == true){
      if (actSetRemove == true) {
        for (int p = 0; p < n_var; p++){
          if (active_var(p) == 1) {
            if (beta(p) == 0){
              active_var(p) = 0;
            }
          }
        }        
      }
      
      vec Current_Change_beta(n_var, fill::zeros); 
      for (int p = 0; p < n_var; p++) {
        if (active_var(p) == 0){
          double lambda_m = lambda * penalized_multiplier(p);
          if (MM == true){
            for (int j = 0; j < n_obs; j++){
              w(j) = v * time(j);
            } 
          } else {
            w.zeros();
            for (int j = 0; j < n_obs; j++){  
              for (int i = 0; i < time(j); i++){ 
                double p_temp = p_binomial_Surv(gamma(i), eta(j));
                w(j) += p_temp * (1 - p_temp);
              }
            }
          }            
          Current_Change_beta(p) = gd_Surv_BetaChange(Z, r, eta, time, w, p, n_obs, lambda_m);
        }
      }
      
      int if_add_new = 0;
      uvec descend_beta_change_index = sort_index(Current_Change_beta, "descend");
      vec descend_beta_change = sort(Current_Change_beta, "descend");
      
      for (int i = 0; i < activeVarNum; i++){
        if (descend_beta_change(i)!= 0){ 
          if_add_new++;
          active_var(descend_beta_change_index(i)) = 1; 
        } else {
          break;
        }
      } 
      
      
      if (if_add_new == 0){
        break;
      }
    } else {
      break; 
    }
    
  }
  
  return make_tuple(beta, gamma, alpha, eta, df, iter);
}

// [[Rcpp::export]]
List pp_DiscSurv_lasso(vec &delta_obs, int max_timepoint, mat &Z, vec &n_prov, vec &time, vec &gamma, vec &beta, vec &alpha, int K0, vec &K1, 
                    vec &sum_failure, vec failure_each_center, vec &lambda_seq, vec &penalized_multiplier, int max_total_iter, 
                    int max_each_iter, double tol, bool backtrack, bool MM, double bound, int initial_active_var, double nvar_max, 
                    bool trace_lambda, int threads, bool actSet, int actIter, int activeVarNum, bool actSetRemove) {
  int n_obs = Z.n_rows, n_beta = Z.n_cols, n_gamma = gamma.n_elem, n_alpha = alpha.n_elem, n_lambda = lambda_seq.n_elem;
  int n_var = K1.n_elem - 1; // n_var: number of penalized variables
  int tol_iter = 0;
  
  mat beta_matrix(n_beta, n_lambda, fill::zeros), gamma_matrix(n_gamma, n_lambda, fill::zeros), alpha_matrix(n_alpha, n_lambda, fill::zeros); // parameter estimation
  mat eta_matrix(n_obs, n_lambda, fill::zeros); // linear predictor
  vec iter_vec(n_lambda, fill::zeros);  // number of iterations for each lambda
  vec df_vec(n_lambda, fill::zeros); // number of non-zero beta
  vec active_var(n_var, fill::zeros); 
  
  if (actSet == true){
    if (K0 == 0){ // if there's no unpenalized beta, we need put the first variable into the active set
      active_var(initial_active_var) = 1;
    } 
  } else {
    active_var.ones(); 
  }
  
  vec ind_start(n_alpha); // index of the first observation of each center
  ind_start(0) = 0;
  for (int i = 1; i < n_alpha; i++) {
    ind_start(i) = ind_start(i - 1) + n_prov(i - 1);
  }
  
  // initialize eta: eta_{kj} = alpha_k + Z_{kj} * beta
  vec eta = rep(alpha, n_prov) + Z * beta;
  
  
  for (int l = 0; l < n_lambda; l++){
    R_CheckUserInterrupt();
    if (trace_lambda == true){
      cout << "processing lambda: " << l + 1 << " (total: " << l + 1 << "/" << n_lambda << ")..." << endl;
    }
    double lambda = lambda_seq(l);
    
    // compare with "pp_Surv": add "n_prov" and "alpha"
    auto fit = pp_DiscSurv_fit(delta_obs, max_timepoint, Z, time, n_prov, gamma, beta, alpha, eta, K0, K1, sum_failure, failure_each_center, lambda, tol_iter, max_total_iter, max_each_iter, penalized_multiplier, backtrack, MM, bound, tol, ind_start, active_var, n_obs, n_var, threads, actSet, actIter, activeVarNum, actSetRemove);
    
    double df_l;
    int iter_l;
    tie(beta, gamma, alpha, eta, df_l, iter_l) = fit;
    beta_matrix.col(l) = beta;
    gamma_matrix.col(l) = gamma;
    alpha_matrix.col(l) = alpha;
    eta_matrix.col(l) = eta;
    df_vec(l) = df_l;
    iter_vec(l) = iter_l;
    
    // check whether the iteration number for the current lambda has reached the maximum
    if (iter_l == max_each_iter) { 
      cout << "Warning: lambda " << l + 1 << "/" << n_lambda << " failed to converge within " << max_each_iter << " iterations!" << endl;
    }
    
    // check whether the entire lambda sequence should stop
    int nv = 0;
    for (int j = 0; j < n_var; j++){
      if (beta(K1(j)) != 0){
        nv++;
      }
    }
    
    //if the current number of penalized variables has already reached nvar_max, then the number of selected variables for the remaining lambda must >= nvar_max.
    if (nv > nvar_max || tol_iter == max_total_iter) { 
      if (tol_iter == max_total_iter) {
        cout << "Algorithm has reached the maximum number of total iterations, stops..." << endl;
      } else {
        cout << "Algorithm has selected the maximum number of penalized variables, stops..." << endl;
      }
      // the estimating process for the remaining lambda will be skipped
      for (int ll = (l + 1); ll < n_lambda; ll++){
        iter_vec(ll) = NA_REAL;
      }
      break;  //break lambda sequence
    }
  }
  List result = List::create(_["gamma"] = gamma_matrix, _["beta"] = beta_matrix, _["alpha"] = alpha_matrix, _["Eta"] = eta_matrix, _["Df"] = df_vec, _["iter"] = iter_vec);
  return result;
}