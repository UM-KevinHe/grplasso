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

double vec_crossprod(vec &w, vec &r) {
  double vec_crossprod = dot(w, r);
  return(vec_crossprod);
}

double mean_crossprod(mat &Z, vec &r, int j, int n_obs) {
  double crossprod = dot(Z.col(j), r);
  return(crossprod/n_obs);
}

double Loglkd(vec &Y, vec &eta) {
  return sum((eta) % Y - log(1 + exp(eta)));
}

vec inner_product(mat &Z) {
  int p = Z.n_cols;
  vec ip(p);
  for (int i = 0; i < p; i++) {
    ip(i) = dot(Z.col(i), Z.col(i));
  }
  return(ip);
}

double weighted_inner_product(mat &Z, vec &w, int j) {
  int n = Z.n_rows;
  double weighted_ip = 0;
  for (int i = 0; i < n; i++){
    weighted_ip += w(i) * Z(i,j) * Z(i,j);
  }
  return(weighted_ip);
}

// cross product of r with jth column of Z
double crossprod(mat &Z, vec &r, int j) {
  double crossprod = dot(Z.col(j), r);
  return(crossprod);
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


// Pr(y=1) for binomial
double p_binomial(double &eta) {
  /*
  if (eta > 10) {
    return(1);
  } else if (eta < -10) {
    return(0);
  } else {
    return(1/(1+exp(-eta)));
  }
  */
  return(1/(1+exp(-eta)));
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



// Function1: "ppLasso"
// "pplasso" is used for solving non-group lasso problems. The algorithms will be very different between that use MM and that do not.

// update penalized beta (simple lasso with MM algorithm)
void gd_glm_lasso_MM(vec &beta, mat &Z, vec &r, vec &eta, vec &old_beta, vec &inner_product_Z, int j, int n_obs, double lambda, double &df, double &MaxChange_beta){
  double beta_initial = crossprod(Z, r, j)/inner_product_Z(j) + old_beta(j);   
  double threshold = n_obs * lambda/(0.25 * inner_product_Z(j));
  double len = Soft_thres(beta_initial, threshold); 

  if (len != old_beta(j)){ // if beta has changed, then update eta and r
    beta(j) = len;
    double beta_change = beta(j) - old_beta(j);
    if (fabs(beta_change) > MaxChange_beta) {
      MaxChange_beta = fabs(beta_change);
    }
    if (beta_change != 0){
      r -= beta_change * Z.col(j);
      eta += beta_change * Z.col(j);
    }
  }
  if (len != 0){
    df += len / beta_initial;
  }
}


// update penalized beta (simple lasso without MM algorithm)
void gd_glm_lasso_noMM(vec &beta, mat &Z, vec &r, vec &eta, vec &old_beta, int j, int n_obs, double lambda, double &df, double &MaxChange_beta){
  vec p(n_obs);
  for (int i = 0; i < n_obs; i++){
    p(i) = p_binomial(eta(i)); 
  }
  vec w = p % (1 - p);  
  double beta_initial = w_crossprod(Z, r, w, j)/weighted_inner_product(Z, w, j) + old_beta(j);   
  double threshold = n_obs * lambda/weighted_inner_product(Z, w, j);
  double len = Soft_thres(beta_initial, threshold); 

  if (len != old_beta(j)){ // if beta has changed, then update eta and r
    beta(j) = len;
    double beta_change = beta(j) - old_beta(j);
    if (fabs(beta_change) > MaxChange_beta) {
      MaxChange_beta = fabs(beta_change);
    }
    if (beta_change != 0){
      r -= beta_change * Z.col(j);
      eta += beta_change * Z.col(j);
    }
  }
  if (len != 0){
    df += len / beta_initial;
  }
}

// "XXX_BetaChange" are functions to determine which variables should be added to the new active set
// We put those betas that can be updated the most currently into the next active set
double gd_glm_lasso_MM_BetaChange(mat &Z, vec &r, vec &inner_product_Z, int j, int n_obs, double lambda){
  double beta_initial = crossprod(Z, r, j)/inner_product_Z(j);   // old_beta = 0;
  double threshold = n_obs * lambda/(0.25 * inner_product_Z(j));
  double len = Soft_thres(beta_initial, threshold); 

  if (len != 0){ //any beta's that are not in the current active set should be zero!
    return(fabs(len)); //"len" equals beta_change
  } else {
    return(0); //if beta hasn't been updated, then 
  }
}


double gd_glm_lasso_noMM_BetaChange(mat &Z, vec &r, vec &eta, int j, int n_obs, double lambda){
  vec p(n_obs);
  for (int i = 0; i < n_obs; i++){
    p(i) = p_binomial(eta(i)); 
  }
  vec w = p % (1 - p);  
  double beta_initial = w_crossprod(Z, r, w, j)/weighted_inner_product(Z, w, j);   // old_beta = 0;
  double threshold = n_obs * lambda/weighted_inner_product(Z, w, j);
  double len = Soft_thres(beta_initial, threshold); 
  
  if (len != 0){ 
    return(fabs(len)); 
  } else {
    return(0);
  }
}


tuple<vec, vec, vec, double, double, int> pp_lasso_fit(vec &Y, mat &Z, vec &n_prov, vec gamma, vec beta, vec eta, int K0, vec &K1,  
                                                        double lambda, int &tol_iter, int max_total_iter, int max_each_iter, vec &penalized_multiplier, int max_n_prov, 
                                                        bool backtrack, bool MM, double bound, double tol, vec &ind_start, vec &active_var, 
                                                        int n_beta, int n_gamma, int n_obs, int n_var, bool single_intercept, int threads,
                                                        bool actSet, int actIter, int activeVarNum, bool actSetRemove){
  vec old_beta = beta, old_gamma = gamma, p(n_obs), r(n_obs), r_shift, Z_tmp;
  double Dev, df, MaxChange_beta, shift;
  double v = 0.25, omega_min = 1e-10;
  int iter = 0; //"iter" counts the number of iterations for each lambda
  vec inner_product_Z = inner_product(Z);

  while (tol_iter < max_total_iter && iter < max_each_iter) { //"tol_iter" counts the number of iterations for the entire lambda sequence
    int inner_loop_iter = 0; // count the number of iterations for a new updated active set
    R_CheckUserInterrupt();
    // inner loop: update variables in the current active set

    while (tol_iter < max_total_iter && iter < max_each_iter && inner_loop_iter < actIter) {
      // the maximum number of inner iterations is "actIter", and after that we will update the current active set;
      R_CheckUserInterrupt();
      df = 0;
      tol_iter++;
      iter++;
      inner_loop_iter++;
      MaxChange_beta = 0;

      for (int i = 0; i < n_obs; i++){
        p(i) = p_binomial(eta(i)); 
      }

      if (single_intercept == true) {  // use coordinate descent to update intercept
        if (MM == true){
          r = (Y - p)/v; //initial pseudo residual vector
          shift = sum(r)/n_obs;  //add a "1" column of Z as the intercept
          gamma += shift;
          eta += shift;
          r -= shift;
        } else {
          vec w = p % (1 - p);
          if (any(w == 0)) {
            w.replace(0, omega_min); 
          }
          r = (Y - p)/w;
          shift = vec_crossprod(w, r)/sum(w);
          gamma += shift;
          eta += shift;
          r -= shift;
        }
      } else {  //use Newton method to update gamma
        double info_gamma;
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
        }

        for (int i = 0; i < n_obs; i++){
          p(i) = p_binomial(eta(i)); 
        }
        //initial pseudo residual vector if using Newton Method to update gamma effect
        if (MM == true){
          r = (Y - p)/v; 
        } else {
          vec w = p % (1 - p);
          if (any(w == 0)) {
            w.replace(0, omega_min); 
          }
          r = (Y - p)/w;  
        }
      }

      // update beta
      uvec update_order_unpenalized = randperm(K0);  // Randomized coordinate descent for unpenalized covariate 
      uvec update_order_penalized = randperm(n_var);
      if (MM == true){ // Use MM algorithm, where w = 0.25. 
        // update unpenalized beta
        for (int j = 0; j < K0; j++){  // If K0 is zero, then the whole iteration will be skipped
          shift = crossprod(Z, r, update_order_unpenalized(j))/inner_product_Z(update_order_unpenalized(j));
          if (fabs(shift) > MaxChange_beta) {
            MaxChange_beta = fabs(shift);
          }
          beta(update_order_unpenalized(j)) = old_beta(update_order_unpenalized(j)) + shift; 
          r -= Z.col(update_order_unpenalized(j)) * shift;
          eta += Z.col(update_order_unpenalized(j)) * shift;  
          df++;
        }

        // update penalized beta
        for (int j = 0; j < n_var; j++) {
          int update_column_index = K1(update_order_penalized(j));
          if (active_var(update_order_penalized(j)) == 1){
            double lambda_m = lambda * penalized_multiplier(update_order_penalized(j));
            gd_glm_lasso_MM(beta, Z, r, eta, old_beta, inner_product_Z, update_column_index, n_obs, lambda_m, df, MaxChange_beta);
          }
        }
      } else {  //without using MM algorithm
        for (int j = 0; j < K0; j++){  // If K0 is zero, then the whole iteration will be skipped
          // update unpenalized beta
          for (int i = 0; i < n_obs; i++){
            p(i) = p_binomial(eta(i)); 
          }
          vec w = p % (1 - p);
          shift = w_crossprod(Z, r, w, update_order_unpenalized(j))/weighted_inner_product(Z, w, update_order_unpenalized(j));
          if (fabs(shift) > MaxChange_beta) {
            MaxChange_beta = fabs(shift);
          }
          beta(update_order_unpenalized(j)) = old_beta(update_order_unpenalized(j)) + shift; 
          r -= Z.col(update_order_unpenalized(j)) * shift;
          eta += Z.col(update_order_unpenalized(j)) * shift;  
          df++;
        }

        // update penalized beta
        for (int j = 0; j < n_var; j++) {
          int update_column_index = K1(update_order_penalized(j));
          if (active_var(update_order_penalized(j)) == 1){
            double lambda_m = lambda * penalized_multiplier(update_order_penalized(j));            
            gd_glm_lasso_noMM(beta, Z, r, eta, old_beta, update_column_index, n_obs, lambda_m, df, MaxChange_beta);
          }
        }
      }

      old_gamma = gamma;
      old_beta = beta;
      //cout << "MaxChange_beta (inner): " << MaxChange_beta << endl;

      if (MaxChange_beta < tol){ // if beta's in the current active set has converged, then jump out of inner loop.
        break;
      }
    }

    /* the following chunk is for checking my code...
    cout << "check active set..." << endl;
    for (int i = 0; i < n_var; i++){
      if (active_var(i) == 1){
        cout << i << " ";
      }
    }    
    cout << endl;
    cout << "#####---------------------------#####" << endl;
    */

    if (actSet == true){
      //outer loop: update active set
      if (actSetRemove == true) {
        // remove zero beta's from the current active set
        for (int j = 0; j < n_var; j++){
          if (active_var(j) == 1) {
            if (beta(j) == 0){
              active_var(j) = 0;
            }
          }
        }        
      }
    
     vec Current_Change_beta(n_var, fill::zeros);  // Store the one-step updates of the beta's that are not in the current active set
    
     if (MM == true){
        for (int j = 0; j < n_var; j++) {
          if (active_var(j) == 0) { //only check variables that are not in the current active set
            // We checked the current amount of updates for all betas (but the beta has not actually been updated)
            double lambda_m = lambda * penalized_multiplier(j);
           // "Current_Change_beta" is a vector that can contains "zero"
            Current_Change_beta(j) = gd_glm_lasso_MM_BetaChange(Z, r, inner_product_Z, j, n_obs, lambda_m); 
          } 
        // Note that the "Current_Change_beta" of those beta's have already in the current active set should be zero!
        }
      } else {
        for (int j = 0; j < n_var; j++) {
          if (active_var(j) == 0) { //only check variables that are not in the current active set
            double lambda_m = lambda * penalized_multiplier(j);
            Current_Change_beta(j) = gd_glm_lasso_noMM_BetaChange(Z, r, eta, j, n_obs, lambda_m); 
         } 
        }
      }

      int if_add_new = 0;
      uvec descend_beta_change_index = sort_index(Current_Change_beta, "descend"); // find index of largest "Current_Change_beta".
      vec descend_beta_change = sort(Current_Change_beta, "descend"); // find largest "Current_Change_beta"

      for (int i = 0; i < activeVarNum; i++){ // only select the largest "actNum" beta's into the new active set
        if (descend_beta_change(i)!= 0){ 
          if_add_new++;
          // the first "actNum" numbers of "descend_beta_change_index" can still contains zero,
          // we should not include them into our new active set.
          // e.g. The current "active_var" is (1, 0, 0, 1, 0), and "descend_beta_change" is (beta2 = 0.4, beta3 = 0.3, beta5 = 0, beta1 = 0, beta4 = 0), and we will choose maximum 3 beta's,
          // we should only add beta2 & beta3 in the current active set, and new "active_var" should be (1, 1, 1, 1, 0)
          active_var(descend_beta_change_index(i)) = 1; 
        } else {
          break;  // remaining "descend_beta_change" should be zero
        }
      }

      if (if_add_new == 0){
        break;
      }
    } else {
      break; // if we don't use "active set", then inner loop converges implies the whole algorithm has converged.
    }
  }  //outer loop ends


  for (int i = 0; i < n_obs; i++){
    p(i) = p_binomial(eta(i)); 
  }
  Dev = Deviance(Y, p);
  return make_tuple(beta, gamma, eta, Dev, df, iter);
}



// [[Rcpp::export]]
List pp_lasso(vec &Y, mat &Z, vec &n_prov, vec &gamma, vec &beta, int K0, vec &K1, vec &lambda_seq, bool lambda_early_stop,
              double stop_dev_ratio, vec &penalized_multiplier, int max_total_iter, int max_each_iter, double tol, double nullDev, 
              bool backtrack, bool MM, double bound, int initial_active_var, double nvar_max, bool trace_lambda, bool single_intercept, 
              int threads, bool actSet, int actIter, int activeVarNum, bool actSetRemove) {
  int n_obs = Z.n_rows, n_beta = Z.n_cols, n_gamma = n_prov.n_elem, n_lambda = lambda_seq.n_elem, max_n_prov = max(n_prov);
  int n_var = K1.n_elem - 1; // n_var: number of penalized variables
  int tol_iter = 0;
  
  mat beta_matrix(n_beta, n_lambda, fill::zeros), gamma_matrix(n_gamma, n_lambda, fill::zeros); // parameter estimation
  mat eta_matrix(n_obs, n_lambda, fill::zeros); // linear predictor
  vec Dev_vec(n_lambda, fill::zeros);  // "Deviance" is defined as "loss"
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
    auto fit = pp_lasso_fit(Y, Z, n_prov, gamma, beta, eta, K0, K1, lambda, tol_iter, max_total_iter, max_each_iter, penalized_multiplier, max_n_prov, backtrack, MM, bound, tol, ind_start, active_var, n_beta, n_gamma, n_obs, n_var, single_intercept, threads, actSet, actIter, activeVarNum, actSetRemove);
    double Dev_l, df_l;
    int iter_l;
    tie(beta, gamma, eta, Dev_l, df_l, iter_l) = fit;
    beta_matrix.col(l) = beta;
    gamma_matrix.col(l) = gamma;
    eta_matrix.col(l) = eta;
    Dev_vec(l) = Dev_l;
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
    //if the current number of variables has already reached nvar_max, then the number of selected variables for the remaining lambda must >= nvar_max.
    if (nv > nvar_max || tol_iter == max_total_iter) { 
      if (tol_iter == max_total_iter) {
        cout << "Algorithm has reached the maximum number of total iterations, stops..." << endl;
      } else {
        cout << "Algorithm has selected the maximum number of variables, stops..." << endl;
      }
      // the estimating process for the remaining lambda will be skipped
      for (int ll = (l + 1); ll < n_lambda; ll++){
        iter_vec(ll) = NA_REAL;
      }
      break;  //break lambda sequence
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



// Function2: "grplasso"
// "grplasso" is used for solving group lasso problems. The algorithm is based on MM algorithm.

void gd_glm_grplasso(vec &beta, mat &Z, vec &r, vec &eta, vec &old_beta, int g, vec &K1, int n_obs, double lambda, double &df, double &MaxChange_beta){
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



// "gd_glm_BetaChange" is used for determining which groups should be added to the new active set
double gd_glm_grplasso_BetaChange(mat &Z, vec &r, int g, vec &K1, int n_obs, double lambda){
  int K = K1(g + 1) - K1(g); //number of features in group g
  vec beta_initial(K);
  for (int j = K1(g); j < K1(g + 1); j++){
    beta_initial(j - K1(g)) = mean_crossprod(Z, r, j, n_obs); // "old_beta" should be zero
  }
  double beta_initial_norm = arma::norm(beta_initial, 2); 
  double len = Soft_thres(beta_initial_norm, lambda/0.25); // "len" is a non-negative number

  if (len != 0){ // which means all the beta's within this group will be updated (to non-zero) and added in the new active set.
    return(len); 
  } else {
    return(0);
  }
}

// estimataion for one given lambda
tuple<vec, vec, vec, double, double, int> grp_lasso_fit(vec &Y, mat &Z, vec &n_prov, vec gamma, vec beta, vec eta, int K0, vec &K1, double lambda, 
                                                        int &tol_iter, int max_total_iter, int max_each_iter, vec &group_multiplier, int max_n_prov, bool backtrack, 
                                                        double bound, double tol, vec &ind_start, vec &active_group, int n_beta, int n_gamma, 
                                                        int n_obs, int n_group, bool single_intercept, int threads, bool actSet, int actIter, 
                                                        int activeGroupNum, bool actSetRemove){
  vec old_beta = beta, old_gamma = gamma, p(n_obs), r(n_obs), r_shift, Z_tmp;
  double Dev, df, MaxChange_beta, shift, lambda_g, v = 0.25;
  int iter = 0;

  while (tol_iter < max_total_iter) {
    int inner_loop_iter = 0; 
    R_CheckUserInterrupt();
    while (tol_iter < max_total_iter && iter < max_each_iter) {
      R_CheckUserInterrupt();
      df = 0;
      tol_iter++;
      iter++;
      inner_loop_iter++;
      MaxChange_beta = 0;
      for (int i = 0; i < n_obs; i++){
        p(i) = p_binomial(eta(i)); 
      }
      
      //auto t1 = high_resolution_clock::now();
      if (single_intercept == true) {  // use coordinate descent to update intercept
        // Since we use the "MM algorithm", the "surrogate function" will be different for different iterations;
        // Consequently, we cannot directly use the "pseudo residual" from previous iteration!
        // In other words, updating the "pseudo residual" only works when using the same surrogate function.
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
        }
        for (int i = 0; i < n_obs; i++){
          p(i) = p_binomial(eta(i)); 
        }
        // note that "pseudo residual" is based on "group descent method", so we should re-define initial r within each iteration after Newton method has been performed
        r = (Y - p)/v;  //initial pseudo residual vector
      }

      //auto t2 = high_resolution_clock::now();
      //duration<double, std::milli> ms_double1 = t2 - t1; 
      //cout << "<#> Update Gamma:" << ms_double1.count() << "ms\n";

      //auto t3 = high_resolution_clock::now();
      uvec update_order_unpenalized = randperm(K0);  // Randomized coordinate descent for unpenalized covariate
      for (int j = 0; j < K0; j++){  // If K0 is zero, then the whole iteration will be skipped
        shift = mean_crossprod(Z, r, update_order_unpenalized(j), n_obs);
        if (fabs(shift) > MaxChange_beta) {
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
          gd_glm_grplasso(beta, Z, r, eta, old_beta, g, K1, n_obs, lambda_g, df, MaxChange_beta);
        }
      }
      //auto t6 = high_resolution_clock::now();
      //duration<double, std::milli> ms_double3 = t6 - t5; 
      //cout << "<<<###>>> Update Penalized Beta:" << ms_double3.count() << "ms\n";

      old_gamma = gamma;
      old_beta = beta;

      if (MaxChange_beta < tol){
        break;
      }
    }

    /* the following chunk is for checking my code...
    cout << "check active group..." << endl;
    for (int g = 0; g < n_group; g++){
      if (active_group(g) == 1){
        cout << g << " ";
      }
    }    
    cout << endl;
    cout << "#####---------------------------#####" << endl;
    */

    if (actSet == true){
      // outer loop: update active set
      if (actSetRemove == true){
      // remove zero groups from the current active set
        for (int g = 0; g < n_group; g++){
          if (active_group(g) == 1) {
            if (beta(K1(g)) == 0){
              active_group(g) = 0;
            }
          }
        }
      }
      
      vec Current_len_group(n_group, fill::zeros);  
    
      for (int g = 0; g < n_group; g++) {
        if (active_group(g) == 0) { 
          double lambda_g = lambda * group_multiplier(g);
          Current_len_group(g) = gd_glm_grplasso_BetaChange(Z, r, g, K1, n_obs, lambda_g); 
        }
      }

      int if_add_new = 0;
      uvec descend_len_index = sort_index(Current_len_group, "descend");
      vec descend_len_value = sort(Current_len_group, "descend");

      for (int i = 0; i < activeGroupNum; i++){ 
        if (descend_len_value(i)!= 0){ 
          if_add_new++;
          active_group(descend_len_index(i)) = 1; 
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

  for (int i = 0; i < n_obs; i++){
    p(i) = p_binomial(eta(i)); 
  }
  Dev = Deviance(Y, p);
  return make_tuple(beta, gamma, eta, Dev, df, iter);
}



// [[Rcpp::export]]
List grp_lasso(vec &Y, mat &Z, vec &n_prov, vec &gamma, vec &beta, int K0, vec &K1, vec &lambda_seq, bool lambda_early_stop, 
               double stop_dev_ratio, vec &group_multiplier, int max_total_iter, int max_each_iter, double tol, double nullDev, 
               bool backtrack, double bound, int initial_active_group, double nvar_max, double group_max, bool trace_lambda, 
               bool single_intercept, int threads, bool actSet, int actIter, int activeGroupNum, bool actSetRemove) {
  int n_obs = Z.n_rows, n_beta = Z.n_cols, n_gamma = n_prov.n_elem, n_lambda = lambda_seq.n_elem, max_n_prov = max(n_prov);
  int n_group = K1.n_elem - 1; // n_group: number of penalized group
  int tol_iter = 0;
  
  mat beta_matrix(n_beta, n_lambda, fill::zeros), gamma_matrix(n_gamma, n_lambda, fill::zeros); // parameter estimation
  mat eta_matrix(n_obs, n_lambda, fill::zeros); // linear predictor
  vec Dev_vec(n_lambda, fill::zeros);  // "Deviance" is defined as "loss"
  vec iter_vec(n_lambda, fill::zeros);  // number of iterations for each lambda
  vec df_vec(n_lambda, fill::zeros); // number of non-zero beta
  vec active_group(n_group, fill::zeros); 

  if (actSet == true){
    if (K0 == 0){ // if there's no unpenalized beta, we need put the first variable into the active set
      active_group(initial_active_group) = 1;
    } 
  } else {
    active_group.ones(); 
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
    auto fit = grp_lasso_fit(Y, Z, n_prov, gamma, beta, eta, K0, K1, lambda, tol_iter, max_total_iter, max_each_iter, group_multiplier, max_n_prov, backtrack, bound, tol, ind_start, active_group, n_beta, n_gamma, n_obs, n_group, single_intercept, threads, actSet, actIter, activeGroupNum, actSetRemove); 
    double Dev_l, df_l;
    int iter_l;
    tie(beta, gamma, eta, Dev_l, df_l, iter_l) = fit;
    beta_matrix.col(l) = beta;
    gamma_matrix.col(l) = gamma;
    eta_matrix.col(l) = eta;
    Dev_vec(l) = Dev_l;
    df_vec(l) = df_l;
    iter_vec(l) = iter_l;

    // check whether the iteration number for the current lambda has reached the maximum
    if (iter_l == max_each_iter) { 
      cout << "Warning: lambda " << l + 1 << "/" << n_lambda << " failed to converge within " << max_each_iter << " iterations!" << endl;
    }

    // check dfmax, gmax; (note, "nv" doesn't equal "df")
    int ng = 0, nv = 0;
    for (int g = 0; g < n_group; g++){
      if (beta(K1(g)) != 0){
         ng++;
         nv += (K1(g + 1) - K1(g));
      }
    }
    if (ng > group_max || nv > nvar_max || tol_iter == max_total_iter) {
      if (tol_iter == max_total_iter) {
        cout << "Algorithm has reached the maximum number of total iterations, stops..." << endl;
      } else if (ng > group_max) {
        cout << "Algorithm has selected the maximum number of groups, stops..." << endl;
      } else {
        cout << "Algorithm has selected the maximum number of variables, stops..." << endl;
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
