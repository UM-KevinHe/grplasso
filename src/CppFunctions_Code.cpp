#define Check_Headers
#include <RcppArmadillo.h>
#include <iostream>
#include <cmath>
#include <omp.h>
#include <chrono>
#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
using namespace Rcpp;
using namespace std;
using namespace arma;

arma::vec rep(arma::vec &x,arma::vec &each) {
  arma::vec x_rep(arma::sum(each));
  int ind = 0, m = x.n_elem;
  for (int i = 0; i < m; i++) {
    x_rep.subvec(ind,ind+each(i)-1) = x(i) * ones(each(i));
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
  return(crossprod/n_obs);
}

double Loglkd(arma::vec &Y, arma::vec &eta) {
  return sum((eta) % Y - log(1 + exp(eta)));
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
  for (int i = 0; i < n; i++){
    weighted_ip += w(i) * Z(i,j) * Z(i,j);
  }
  return(weighted_ip);
}

// cross product of r with jth column of Z
double crossprod(arma::mat &Z, arma::vec &r, int j) {
  double crossprod = dot(Z.col(j), r);
  return(crossprod);
}

// Weighted cross product of r with jth column of Z
double w_crossprod(arma::mat &Z, arma::vec &r, arma::vec &w, int j) {
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
double Z_max_grLasso(arma::mat &x, arma::vec &r, arma::vec &K, arma::vec &m){ // "K": a vector contains the start index of each group; "m": group.multiplier
  int J = K.n_elem - 1;  //number of penalized group
  double z_max = 0, z;
  for (int g = 0; g < J; g++){
    int Kg = K(g + 1) - K(g); //number of features in group g
    arma::vec Z(Kg);
    for (int j = K(g); j < K(g + 1); j++) {
      arma::vec x_tmp = x.col(j);
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
double Deviance(arma::vec &Y, arma::vec &p){
  double Dev = 0;
  int n_obs = Y.n_elem;
  for (int i = 0; i < n_obs; i++){
    if (p(i) != 0 && p(i) != 1){
      Dev -= 2 * Y(i) * log(p(i)) + 2 * (1 - Y(i)) * log(1 - p(i));
    }
  }
  return(Dev);
}

double Loglkd_Surv(int n_obs, arma::vec &delta_obs, arma::vec &time, arma::vec &gamma, arma::vec &eta) {
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

// where eta = X\beta
double p_binomial_Surv(double &gamma, double &eta) {
  return(1/(1+exp(-gamma - eta)));
}


// [[Rcpp::export]]
arma::vec DiscSurv_residuals(int n_obs, arma::vec &delta_obs, arma::vec &time, arma::vec &alpha, arma::vec &eta){
  arma::vec residuals(n_obs);
  for (int j = 0; j < n_obs; j++){  //j: individual
    for (int i = 0; i < time(j); i++){ //i: time point
      residuals(j) += p_binomial_Surv(alpha(i), eta(j));
    }
  }
  residuals = - residuals + delta_obs;
  return(residuals);
}

// [[Rcpp::export]]
arma::mat predict_linear_predictor(int n_lambda, int n_obs, int expand_n_obs, arma::vec &time, arma::mat &gamma, arma::mat &eta) {
  arma::mat predict_linear_predictor(expand_n_obs, n_lambda);
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


// Function1: "ppLasso"
// "pplasso" is used for solving non-group lasso problems. The algorithms will be very different between that use MM and that do not.

// update penalized beta (simple lasso with MM algorithm)
void gd_glm_lasso_MM(arma::vec &beta, arma::mat &Z, arma::vec &r, arma::vec &eta, arma::vec &old_beta, arma::vec &inner_product_Z, int j, int n_obs, double lambda, double &df, double &MaxChange_beta){
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

void gd_glm_lasso_noMM(arma::vec &beta, arma::vec &w, arma::mat &Z, arma::vec &r, arma::vec &eta, arma::vec &old_beta, int j, int n_obs, double lambda, double &df, double &MaxChange_beta){
  double w_z_square = weighted_inner_product(Z, w, j);
  double beta_initial = w_crossprod(Z, r, w, j)/w_z_square + old_beta(j);
  double threshold = n_obs * lambda/w_z_square;
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
double gd_glm_lasso_MM_BetaChange(arma::mat &Z, arma::vec &r, arma::vec &inner_product_Z, int j, int n_obs, double lambda){
  double beta_initial = crossprod(Z, r, j)/inner_product_Z(j);   // old_beta = 0;
  double threshold = n_obs * lambda/(0.25 * inner_product_Z(j));
  double len = Soft_thres(beta_initial, threshold);

  if (len != 0){ //any beta's that are not in the current active set should be zero!
    return(fabs(len)); //"len" equals beta_change
  } else {
    return(0); //if beta hasn't been updated, then
  }
}


double gd_glm_lasso_noMM_BetaChange(arma::mat &Z, arma::vec &r, arma::vec &eta, int j, int n_obs, double lambda){
  arma::vec p(n_obs);
  for (int i = 0; i < n_obs; i++){
    p(i) = p_binomial(eta(i));
  }
  arma::vec w = p % (1 - p);
  double beta_initial = w_crossprod(Z, r, w, j)/weighted_inner_product(Z, w, j);   // old_beta = 0;
  double threshold = n_obs * lambda/weighted_inner_product(Z, w, j);
  double len = Soft_thres(beta_initial, threshold);

  if (len != 0){
    return(fabs(len));
  } else {
    return(0);
  }
}


tuple<arma::vec, arma::vec, arma::vec, double, double, int> pp_lasso_fit(arma::vec &Y, arma::mat &Z, arma::vec &n_prov, arma::vec gamma, arma::vec beta, arma::vec eta, 
                                                                         int K0, arma::vec &K1,  double lambda, int &tol_iter, int max_total_iter, int max_each_iter, 
                                                                         arma::vec &penalized_multiplier, int max_n_prov, bool backtrack, bool MM, double bound, double tol, 
                                                                         arma::vec &ind_start, arma::vec &active_var, int n_obs, int n_var, bool single_intercept, int threads, 
                                                                         bool actSet, int actIter, int activeVarNum, bool actSetRemove){
  arma::vec old_beta = beta, old_gamma = gamma, p(n_obs), r(n_obs), r_shift;
  int n_gamma = gamma.n_elem;
  double Dev, df, MaxChange_beta, shift;
  double v = 0.25, omega_min = 1e-20;
  int iter = 0; //"iter" counts the number of iterations for each lambda
  arma::vec inner_product_Z = inner_product(Z), w(n_obs);

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
          w = p % (1 - p); //initial "w" for current iteration; "w" here doesn't change during within current iteration!
          if (any(w == 0)) {
            w.replace(0, omega_min);
          }
          r = (Y - p)/w;  //initial pseudo residual
          shift = vec_crossprod(w, r)/sum(w);
          gamma += shift;
          eta += shift;
          r -= shift;
        }
      } else {  //use Newton method to update gamma
        double info_gamma;
        int nProcessors = threads;
        arma::vec score_gamma(n_gamma), d_gamma(n_gamma), Yp(n_obs), pq(n_obs), gamma_obs(n_obs);
        if (backtrack == true) {
          arma::vec gamma_shift_tmp, eta_tmp(n_obs);
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
          arma::vec gamma_shift = gamma - old_gamma;
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
          arma::vec gamma_shift = gamma - old_gamma;
          eta += rep(gamma_shift, n_prov);
        }

        //initialize "w" and pseudo residual "r"  if using Newton Method to update gamma
        for (int i = 0; i < n_obs; i++){
          p(i) = p_binomial(eta(i));
        }

        if (MM == true){
          r = (Y - p)/v;
        } else {
          w = p % (1 - p); //again, "w" here doesn't change during within current iteration!
          if (any(w == 0)) {
            w.replace(0, omega_min);
          }
          r = (Y - p)/w;
        }
      }

      // update beta
      arma::uvec update_order_unpenalized = randperm(K0);  // Randomized coordinate descent for unpenalized covariate
      arma::uvec update_order_penalized = randperm(n_var);
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
            gd_glm_lasso_noMM(beta, w, Z, r, eta, old_beta, update_column_index, n_obs, lambda_m, df, MaxChange_beta);
          }
        }

      }

      old_gamma = gamma;
      old_beta = beta;

      if (MaxChange_beta < tol){ // if beta's in the current active set has converged, then jump out of inner loop.
        break;
      }
    }

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

     arma::vec Current_Change_beta(n_var, fill::zeros);  // Store the one-step updates of the beta's that are not in the current active set

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
      arma::uvec descend_beta_change_index = sort_index(Current_Change_beta, "descend"); // find index of largest "Current_Change_beta".
      arma::vec descend_beta_change = sort(Current_Change_beta, "descend"); // find largest "Current_Change_beta"

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
List pp_lasso(arma::vec &Y, arma::mat &Z, arma::vec &n_prov, arma::vec &gamma, arma::vec &beta, int K0, arma::vec &K1, arma::vec &lambda_seq, bool lambda_early_stop,
              double stop_dev_ratio, arma::vec &penalized_multiplier, int max_total_iter, int max_each_iter, double tol, double nullDev,
              bool backtrack, bool MM, double bound, int initial_active_var, double nvar_max, bool trace_lambda, bool single_intercept,
              int threads, bool actSet, int actIter, int activeVarNum, bool actSetRemove) {
  int n_obs = Z.n_rows, n_beta = Z.n_cols, n_gamma = n_prov.n_elem, n_lambda = lambda_seq.n_elem, max_n_prov = max(n_prov);
  int n_var = K1.n_elem - 1; // n_var: number of penalized variables
  int tol_iter = 0;

  arma::mat beta_matrix(n_beta, n_lambda, fill::zeros), gamma_matrix(n_gamma, n_lambda, fill::zeros); // parameter estimation
  arma::mat eta_matrix(n_obs, n_lambda, fill::zeros); // linear predictor
  arma::vec Dev_vec(n_lambda, fill::zeros);  // "Deviance" is defined as "loss"
  arma::vec iter_vec(n_lambda, fill::zeros);  // number of iterations for each lambda
  arma::vec df_vec(n_lambda, fill::zeros); // number of non-zero beta
  arma::vec active_var(n_var, fill::zeros);

  if (actSet == true){
    if (K0 == 0){ // if there's no unpenalized beta, we need put the first variable into the active set
      active_var(initial_active_var) = 1;
    }
  } else {
    active_var.ones();
  }


  arma::vec ind_start(n_gamma); // index of the first observation of each provider
  ind_start(0) = 0;
  for (int i = 1; i < n_gamma; i++) {
    ind_start(i) = ind_start(i - 1) + n_prov(i - 1);
  }

  // initialize eta
  arma::vec eta = rep(gamma, n_prov) + Z * beta;

  for (int l = 0; l < n_lambda; l++){
    R_CheckUserInterrupt();
    if (trace_lambda == true){
      cout << "processing lambda: " << l + 1 << " (total: " << l + 1 << "/" << n_lambda << ")..." << endl;
    }
    double lambda = lambda_seq(l);
    auto fit = pp_lasso_fit(Y, Z, n_prov, gamma, beta, eta, K0, K1, lambda, tol_iter, max_total_iter, max_each_iter, penalized_multiplier, max_n_prov, backtrack, MM, bound, tol, ind_start, active_var, n_obs, n_var, single_intercept, threads, actSet, actIter, activeVarNum, actSetRemove);
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
        cout << "Algorithm has selected the maximum number of penalized variables, stops..." << endl;
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

void gd_glm_grplasso(arma::vec &beta, arma::mat &Z, arma::vec &r, arma::vec &eta, arma::vec &old_beta, int g, arma::vec &K1, int n_obs, double lambda, double &df, double &MaxChange_beta){
  int K = K1(g + 1) - K1(g); //number of features in group g
  arma::vec beta_initial(K);
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
double gd_glm_grplasso_BetaChange(arma::mat &Z, arma::vec &r, int g, arma::vec &K1, int n_obs, double lambda){
  int K = K1(g + 1) - K1(g); //number of features in group g
  arma::vec beta_initial(K);
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
tuple<arma::vec, arma::vec, arma::vec, double, double, int> grp_lasso_fit(arma::vec &Y, arma::mat &Z, arma::vec &n_prov, arma::vec gamma, arma::vec beta, arma::vec eta, 
                                                                          int K0, arma::vec &K1, double lambda, int &tol_iter, int max_total_iter, int max_each_iter, 
                                                                          arma::vec &group_multiplier, int max_n_prov, bool backtrack, double bound, double tol, 
                                                                          arma::vec &ind_start, arma::vec &active_group, int n_obs, int n_group, bool single_intercept, 
                                                                          int threads, bool actSet, int actIter, int activeGroupNum, bool actSetRemove){
  arma::vec old_beta = beta, old_gamma = gamma, p(n_obs), r(n_obs), r_shift;
  int n_gamma = gamma.n_elem;
  double Dev, df, MaxChange_beta, shift, lambda_g, v = 0.25;
  int iter = 0;

  while (tol_iter < max_total_iter) {
    int inner_loop_iter = -1;
    inner_loop_iter = inner_loop_iter +1; //just for removing unused warning.
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
        arma::vec score_gamma(n_gamma), d_gamma(n_gamma), Yp(n_obs), pq(n_obs), gamma_obs(n_obs);
        if (backtrack == true) {
          arma::vec gamma_shift_tmp, eta_tmp(n_obs);
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
          arma::vec gamma_shift = gamma - old_gamma;
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
          arma::vec gamma_shift = gamma - old_gamma;
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
      arma::uvec update_order_unpenalized = randperm(K0);  // Randomized coordinate descent for unpenalized covariate
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

      // update penalized beta
      // note that all groups are iterated without "active set method"
      for (int g = 0; g < n_group; g++){
        if (active_group(g) == 1){
          lambda_g = lambda * group_multiplier(g);
          gd_glm_grplasso(beta, Z, r, eta, old_beta, g, K1, n_obs, lambda_g, df, MaxChange_beta);
        }
      }
      old_gamma = gamma;
      old_beta = beta;

      if (MaxChange_beta < tol){
        break;
      }
    }

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

      // check whether "non-active" groups can still be updated; if not, algorithm ends
      arma::vec Current_len_group(n_group, fill::zeros);
      for (int g = 0; g < n_group; g++) {
        if (active_group(g) == 0) {
          double lambda_g = lambda * group_multiplier(g);
          Current_len_group(g) = gd_glm_grplasso_BetaChange(Z, r, g, K1, n_obs, lambda_g); //beta's are not truly updated
        }
      }

      int if_add_new = 0; // count how many "non-active" groups can still be further updated
      arma::uvec descend_len_index = sort_index(Current_len_group, "descend");
      arma::vec descend_len_value = sort(Current_len_group, "descend");

      for (int i = 0; i < activeGroupNum; i++){
        if (descend_len_value(i)!= 0){
          if_add_new++;
          active_group(descend_len_index(i)) = 1;
        } else {  // all remaining groups must stay at zero as well
          break;
        }
      }

      if (if_add_new == 0){ //if no new groups can be updated, the algorithm ends
        break;
      }
    } else {  //while without "active set", all the parameters have already been updated to global optimal
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
List grp_lasso(arma::vec &Y, arma::mat &Z, arma::vec &n_prov, arma::vec &gamma, arma::vec &beta, int K0, arma::vec &K1, arma::vec &lambda_seq, 
               bool lambda_early_stop, double stop_dev_ratio, arma::vec &group_multiplier, int max_total_iter, int max_each_iter, double tol, double nullDev,
               bool backtrack, double bound, int initial_active_group, double nvar_max, double group_max, bool trace_lambda,
               bool single_intercept, int threads, bool actSet, int actIter, int activeGroupNum, bool actSetRemove) {
  int n_obs = Z.n_rows, n_beta = Z.n_cols, n_gamma = n_prov.n_elem, n_lambda = lambda_seq.n_elem, max_n_prov = max(n_prov);
  int n_group = K1.n_elem - 1; // n_group: number of penalized group
  int tol_iter = 0;

  arma::mat beta_matrix(n_beta, n_lambda, fill::zeros), gamma_matrix(n_gamma, n_lambda, fill::zeros); // parameter estimation
  arma::mat eta_matrix(n_obs, n_lambda, fill::zeros); // linear predictor
  arma::vec Dev_vec(n_lambda, fill::zeros);  // "Deviance" is defined as "loss"
  arma::vec iter_vec(n_lambda, fill::zeros);  // number of iterations for each lambda
  arma::vec df_vec(n_lambda, fill::zeros); // number of non-zero beta
  arma::vec active_group(n_group, fill::zeros);

  if (actSet == true){
    if (K0 == 0){ // if there's no unpenalized beta, we need put the first variable into the active set
      active_group(initial_active_group) = 1;
    }
  } else {
    active_group.ones();
  }

  arma::vec ind_start(n_gamma); // index of the first observation of each provider
  ind_start(0) = 0;
  for (int i = 1; i < n_gamma; i++) {
    ind_start(i) = ind_start(i - 1) + n_prov(i - 1);
  }

  // initialize eta
  arma::vec eta = rep(gamma, n_prov) + Z * beta;

  for (int l = 0; l < n_lambda; l++){
    R_CheckUserInterrupt();
    if (trace_lambda == true){
      cout << "processing lambda: " << l + 1 << " (total: " << l + 1 << "/" << n_lambda << ")..." << endl;
    }
    double lambda = lambda_seq(l);
    auto fit = grp_lasso_fit(Y, Z, n_prov, gamma, beta, eta, K0, K1, lambda, tol_iter, max_total_iter, max_each_iter, group_multiplier, max_n_prov, backtrack, bound, tol, ind_start, active_group, n_obs, n_group, single_intercept, threads, actSet, actIter, activeGroupNum, actSetRemove);
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


// Function3: discrete survival model with logit-link
void gd_Surv(arma::vec &beta, arma::mat &Z, arma::vec &r, arma::vec &eta, arma::vec &time, arma::vec &old_beta, arma::vec &w, int p, int n_obs, double lambda, double &df, double &MaxChange_beta){
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


double gd_Surv_BetaChange(arma::mat &Z, arma::vec &r, arma::vec &eta, arma::vec &time, arma::vec &w, int p, int n_obs, double lambda){
  double beta_initial = w_crossprod(Z, r, w, p)/weighted_inner_product(Z, w, p);
  double threshold = n_obs * lambda/weighted_inner_product(Z, w, p);
  double len = Soft_thres(beta_initial, threshold);

  if (len != 0){
    return(fabs(len));
  } else {
    return(0);
  }
}


tuple<arma::vec, arma::vec, arma::vec, arma::vec, double, int> pp_DiscSurv_fit(arma::vec &delta_obs, int max_timepoint, arma::mat &Z, arma::vec &time, arma::vec &n_prov, arma::vec gamma,
                                                                               arma::vec beta, arma::vec alpha, arma::vec eta, int K0, arma::vec &K1, arma::vec &sum_failure, arma::vec failure_each_center,
                                                                               double lambda, int &tol_iter, int max_total_iter, int max_each_iter, arma::vec &penalized_multiplier, bool backtrack,
                                                                               bool MM, double bound, double tol, arma::vec &ind_start, arma::vec &active_var, int n_obs, int n_var, int threads,
                                                                               bool actSet, int actIter, int activeVarNum, bool actSetRemove){
  // note: In the current ".cpp" function, "gamma" denotes the time effect, while "alpha" denotes the center effect
  // which are different from the notations wrote in the paper
  
  arma::vec old_beta = beta, old_alpha = alpha, r(n_obs), r_shift, w(n_obs);
  int n_alpha = alpha.n_elem;
  double df, MaxChange_beta, shift;
  double v = 0.25, omega_min = 1e-10;
  double p_gamma;
  int iter = 0;

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
      arma::vec score_gamma = - sum_failure, info_gamma(max_timepoint, fill::zeros), d_gamma(max_timepoint, fill::zeros);
      if (backtrack == true) {
        arma::vec gamma_tmp, eta_tmp(n_obs);
        double loglkd, d_loglkd, u = 1.0, k, s = 0.01, t = 0.8;
        for (int j = 0; j < n_obs; j++){ //iterated over everyone in the dataset
          for (int i = 0; i < time(j); i++){  //iterated over time
            p_gamma = p_binomial_Surv(gamma(i), eta(j));
            score_gamma(i) += p_gamma;
            info_gamma(i) -= p_gamma * (1 - p_gamma); //a diagonal matrix
          }
        }
        int nProcessors = threads;
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

      // 2. update alpha (center effect)
      double temp_p;
      arma::vec p_obs(n_obs, fill::zeros); // each element of p_obs sum over all time points for that observation
      arma::vec pq_obs(n_obs, fill::zeros);
      for (int j = 0; j < n_obs; j++){  //j: individual
        for (int i = 0; i < time(j); i++){ //i: time point
          temp_p = p_binomial_Surv(gamma(i), eta(j));
          p_obs(j) += temp_p;
          pq_obs(j) += temp_p * (1 - temp_p);
        }
      }

      arma::vec score_alpha(n_alpha, fill::zeros), info_alpha(n_alpha, fill::zeros), d_alpha(n_alpha, fill::zeros);
      if (backtrack == true){
        arma::vec alpha_shift_tmp, eta_tmp(n_obs);
        double loglkd, d_loglkd, u2 = 1.0, k2, s2 = 0.01, t2 = 0.8;
        int nProcessors = threads;
        omp_set_num_threads(nProcessors);
        #pragma omp parallel for schedule(static)
        for (int i = 1; i < n_alpha; i++) { // start from the second center
          score_alpha(i) = sum(p_obs(span(ind_start(i), ind_start(i) + n_prov(i) - 1)));
          double tmp_info_alpha = -sum(pq_obs(span(ind_start(i), ind_start(i) + n_prov(i) - 1)));
          info_alpha(i) = std::min(-omega_min, tmp_info_alpha);
        }
        score_alpha -= failure_each_center; // double sum of delta_ij
        for (int i = 1; i < n_alpha; i++) {
          d_alpha(i) = score_alpha(i)/info_alpha(i); // d_alpha of the first center will always be zero (alpha_1 is zero and never update)
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
        arma::vec alpha_shift = alpha - old_alpha;
        eta += rep(alpha_shift, n_prov);
      } else {
        int nProcessors = threads;
        omp_set_num_threads(nProcessors);
        #pragma omp parallel for schedule(static)
        for (int i = 1; i < n_alpha; i++) {
          score_alpha(i) = sum(p_obs(span(ind_start(i), ind_start(i) + n_prov(i) - 1)));
          double tmp_info_alpha = -sum(pq_obs(span(ind_start(i), ind_start(i) + n_prov(i) - 1)));
          info_alpha(i) = std::min(-omega_min, tmp_info_alpha);
        }
        score_alpha -= failure_each_center;
        for (int i = 1; i < n_alpha; i++) {
          d_alpha(i) = score_alpha(i)/info_alpha(i);
        }
        alpha = alpha + d_alpha;
        alpha = clamp(alpha, median(alpha) - bound, median(alpha) + bound);

        arma::vec alpha_shift = alpha - old_alpha;
        eta += rep(alpha_shift, n_prov);
      }

      // 3. update beta
      arma::vec p_eta(n_obs, fill::zeros);
      //initialize pseudo residual
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
      arma::vec score_eta = - delta_obs + p_eta;
      r = - score_eta / w; //initial pseudo outcome

      // 3.1 update unpenalized beta
      arma::uvec update_order_unpenalized = randperm(K0);  // Randomized coordinate descent for unpenalized covariate
      for (int p = 0; p < K0; p++){
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
      arma::uvec update_order_penalized = randperm(n_var);
      for (int p = 0; p < n_var; p++) {
        int update_column_index = K1(update_order_penalized(p));
        if (active_var(update_order_penalized(p)) == 1){
          double lambda_m = lambda * penalized_multiplier(update_order_penalized(p));
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

      arma::vec Current_Change_beta(n_var, fill::zeros);
      for (int p = 0; p < n_var; p++) {
        if (active_var(p) == 0){
          double lambda_m = lambda * penalized_multiplier(p);
          Current_Change_beta(p) = gd_Surv_BetaChange(Z, r, eta, time, w, p, n_obs, lambda_m);
        }
      }

      int if_add_new = 0;
      arma::uvec descend_beta_change_index = sort_index(Current_Change_beta, "descend");
      arma::vec descend_beta_change = sort(Current_Change_beta, "descend");

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
List pp_DiscSurv_lasso(arma::vec &delta_obs, int max_timepoint, arma::mat &Z, arma::vec &n_prov, arma::vec &time, arma::vec &gamma, arma::vec &beta,
                       arma::vec &alpha, int K0, arma::vec &K1, arma::vec &sum_failure, arma::vec failure_each_center, arma::vec &lambda_seq, arma::vec &penalized_multiplier, 
                       int max_total_iter, int max_each_iter, double tol, bool backtrack, bool MM, double bound, int initial_active_var, double nvar_max,
                       bool trace_lambda, int threads, bool actSet, int actIter, int activeVarNum, bool actSetRemove) {
  int n_obs = Z.n_rows, n_beta = Z.n_cols, n_gamma = gamma.n_elem, n_alpha = alpha.n_elem, n_lambda = lambda_seq.n_elem;
  int n_var = K1.n_elem - 1; // n_var: number of penalized variables
  int tol_iter = 0;

  arma::mat beta_matrix(n_beta, n_lambda, fill::zeros), gamma_matrix(n_gamma, n_lambda, fill::zeros), alpha_matrix(n_alpha, n_lambda, fill::zeros); // parameter estimation
  arma::mat eta_matrix(n_obs, n_lambda, fill::zeros); // linear predictor
  arma::vec iter_vec(n_lambda, fill::zeros);  // number of iterations for each lambda
  arma::vec df_vec(n_lambda, fill::zeros); // number of non-zero beta
  arma::vec active_var(n_var, fill::zeros);

  if (actSet == true){
    if (K0 == 0){ // if there's no unpenalized beta, we need put the first variable into the active set
      active_var(initial_active_var) = 1;
    }
  } else {
    active_var.ones();
  }

  arma::vec ind_start(n_alpha); // index of the first observation of each center
  ind_start(0) = 0;
  for (int i = 1; i < n_alpha; i++) {
    ind_start(i) = ind_start(i - 1) + n_prov(i - 1);
  }

  // initialize eta: eta_{kj} = alpha_k + Z_{kj} * beta
  // alpha_1 (effect of the first center) is set to be 0, and will never be updated.
  arma::vec eta = rep(alpha, n_prov) + Z * beta;


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


tuple<vec, vec, vec, double, int> DiscSurv_fit(arma::vec &delta_obs, int max_timepoint, arma::mat &Z, arma::vec &time, arma::vec gamma, arma::vec beta, arma::vec eta, int K0, arma::vec &K1,
                                               arma::vec &sum_failure, double lambda, int &tol_iter, int max_total_iter, int max_each_iter, arma::vec &penalized_multiplier, bool backtrack, 
                                               bool MM, double bound, double tol, arma::vec &active_var, int n_obs, int n_var, int threads, bool actSet, int actIter, int activeVarNum, 
                                               bool actSetRemove){
  // note: In the current ".cpp" function, "gamma" denotes the time effect
  arma::vec old_beta = beta, r(n_obs), r_shift, w(n_obs);
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
      arma::vec score_gamma = - sum_failure, info_gamma(max_timepoint, fill::zeros), d_gamma(max_timepoint, fill::zeros); 
      if (backtrack == true) {
        arma::vec gamma_tmp, eta_tmp(n_obs);
        double loglkd, d_loglkd, u = 1.0, k, s = 0.01, t = 0.8;
        for (int j = 0; j < n_obs; j++){  //j: individual
          for (int i = 0; i < time(j); i++){ //i: time point
            p_gamma = p_binomial_Surv(gamma(i), eta(j));
            score_gamma(i) += p_gamma;
            info_gamma(i) -= p_gamma * (1 - p_gamma); //a diagonal matrix
          }
        }

        int nProcessors = threads; 
        omp_set_num_threads(nProcessors);
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < max_timepoint; i++){
          double tmp_info_gamma = info_gamma(i);
          d_gamma(i) = score_gamma(i)/std::min(-omega_min, tmp_info_gamma);
        }
        loglkd = Loglkd_Surv(n_obs, delta_obs, time, gamma, eta); //eta = Z * beta. It will not change when only update gamma
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
      
      // 2. update beta
      arma::vec p_eta(n_obs, fill::zeros);
      //initialize pseudo residual
      if (MM == true){ // surrogate function use w = 0.25 * diag{k1, k2, ..., kn}, which is a n*n matrix
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
      arma::vec score_eta = p_eta - delta_obs;
      r = - score_eta / w; 
      

      
      // 2.1 update unpenalized beta
      arma::uvec update_order_unpenalized = randperm(K0);  
      for (int p = 0; p < K0; p++){  
        shift = w_crossprod(Z, r, w, update_order_unpenalized(p))/weighted_inner_product(Z, w, update_order_unpenalized(p)); 
        if (fabs(shift) > MaxChange_beta) {
          MaxChange_beta = fabs(shift);
        }
        beta(update_order_unpenalized(p)) = old_beta(update_order_unpenalized(p)) + shift; 
        r -= Z.col(update_order_unpenalized(p)) * shift;
        eta += Z.col(update_order_unpenalized(p)) * shift;  //eta = Z\beta
        df++;
      }
      
      // 2.2 update penalized beta
      arma::uvec update_order_penalized = randperm(n_var);
      for (int p = 0; p < n_var; p++) {
        int update_column_index = K1(update_order_penalized(p));
        if (active_var(update_order_penalized(p)) == 1){
          double lambda_m = lambda * penalized_multiplier(update_order_penalized(p));
          gd_Surv(beta, Z, r, eta, time, old_beta, w, update_column_index, n_obs, lambda_m, df, MaxChange_beta);
        }
      }
      
      old_beta = beta;
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
      
      arma::vec Current_Change_beta(n_var, fill::zeros); 
      for (int p = 0; p < n_var; p++) {
        if (active_var(p) == 0){
          double lambda_m = lambda * penalized_multiplier(p); 
          Current_Change_beta(p) = gd_Surv_BetaChange(Z, r, eta, time, w, p, n_obs, lambda_m);
        }
      }
      
      int if_add_new = 0;
      arma::uvec descend_beta_change_index = sort_index(Current_Change_beta, "descend");
      arma::vec descend_beta_change = sort(Current_Change_beta, "descend");
      
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
  
  return make_tuple(beta, gamma, eta, df, iter);
}

// [[Rcpp::export]]
List DiscSurv_lasso(arma::vec &delta_obs, int max_timepoint, arma::mat &Z, arma::vec &time, arma::vec &gamma, arma::vec &beta, 
                    int K0, arma::vec &K1, arma::vec &sum_failure, arma::vec &lambda_seq, arma::vec &penalized_multiplier, int max_total_iter, 
                    int max_each_iter, double tol, bool backtrack, bool MM, double bound, int initial_active_var, double nvar_max, 
                    bool trace_lambda, int threads, bool actSet, int actIter, int activeVarNum, bool actSetRemove) {
  int n_obs = Z.n_rows, n_beta = Z.n_cols, n_gamma = gamma.n_elem, n_lambda = lambda_seq.n_elem;
  int n_var = K1.n_elem - 1; 
  int tol_iter = 0;
  
  arma::mat beta_matrix(n_beta, n_lambda, fill::zeros), gamma_matrix(n_gamma, n_lambda, fill::zeros); 
  arma::mat eta_matrix(n_obs, n_lambda, fill::zeros); // linear predictor
  arma::vec iter_vec(n_lambda, fill::zeros);  // number of iterations for each lambda
  arma::vec df_vec(n_lambda, fill::zeros); 
  arma::vec active_var(n_var, fill::zeros); 
  
  if (actSet == true){
    if (K0 == 0){ // if there's no unpenalized beta, we need put the first variable into the active set
      active_var(initial_active_var) = 1;
    } 
  } else {
    active_var.ones(); 
  }
  
  // initialize eta
  arma::vec eta = Z * beta;
  
  for (int l = 0; l < n_lambda; l++){
    R_CheckUserInterrupt();
    if (trace_lambda == true){
      cout << "processing lambda: " << l + 1 << " (total: " << l + 1 << "/" << n_lambda << ")..." << endl;
    }
    double lambda = lambda_seq(l);
    
    auto fit = DiscSurv_fit(delta_obs, max_timepoint, Z, time, gamma, beta, eta, K0, K1, sum_failure, lambda, tol_iter, max_total_iter, max_each_iter, penalized_multiplier, backtrack, MM, bound, tol, active_var, n_obs, n_var, threads, actSet, actIter, activeVarNum, actSetRemove);
    
    double df_l;
    int iter_l;
    tie(beta, gamma, eta, df_l, iter_l) = fit;
    beta_matrix.col(l) = beta;
    gamma_matrix.col(l) = gamma;
    eta_matrix.col(l) = eta;
    df_vec(l) = df_l;
    iter_vec(l) = iter_l;
    
    if (iter_l == max_each_iter) { 
      cout << "Warning: lambda " << l + 1 << "/" << n_lambda << " failed to converge within " << max_each_iter << " iterations!" << endl;
    }
    
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
      for (int ll = (l + 1); ll < n_lambda; ll++){
        iter_vec(ll) = NA_REAL;
      }
      break; 
    }
    
  }
  
  List result = List::create(_["gamma"] = gamma_matrix, _["beta"] = beta_matrix, _["Eta"] = eta_matrix, _["Df"] = df_vec, _["iter"] = iter_vec);
  return result;
}





// stratified cox model

void gd_stratCox(arma::vec &beta, arma::mat &Z, arma::vec &r, arma::vec &eta, arma::vec &old_beta, 
                 int g, arma::vec &K1, int n_obs, double lambda, double &df, double &MaxChange_beta){
  int K = K1(g + 1) - K1(g); 
  arma::vec beta_initial(K);
  for (int j = K1(g); j < K1(g + 1); j++){
    beta_initial(j - K1(g)) = mean_crossprod(Z, r, j, n_obs) + old_beta(j);
  }
  double beta_initial_norm = arma::norm(beta_initial, 2);
  double len = Soft_thres(beta_initial_norm, lambda/1); //v = 1

  if (len != 0 || old_beta(K1(g)) != 0){ 
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
  if (len > 0){
    df += K * len / beta_initial_norm;
  }
}

double gd_stratCox_BetaChange(arma::mat &Z, arma::vec &r, int g, arma::vec &K1, int n_obs, double lambda){
  int K = K1(g + 1) - K1(g); //number of features in group g
  arma::vec beta_initial(K);
  for (int j = K1(g); j < K1(g + 1); j++){
    beta_initial(j - K1(g)) = mean_crossprod(Z, r, j, n_obs); // "old_beta" should be zero
  }
  double beta_initial_norm = arma::norm(beta_initial, 2);
  double len = Soft_thres(beta_initial_norm, lambda/1);

  if (len != 0){ 
    return(len);
  } else {
    return(0);
  }
}


tuple<arma::vec, arma::vec, double, double, int> StratCox_lasso_fit(arma::vec &delta_obs, arma::mat &Z, arma::vec &n_each_prov, arma::vec beta, arma::vec eta, int K0, arma::vec &K1,
                                                                    double lambda, int &tol_iter, int max_total_iter, int max_each_iter, arma::vec &group_multiplier, int count_stratum,
                                                                    double tol, arma::vec &ind_start, arma::vec &active_group, int n_obs, int n_group,
                                                                    bool actSet, int actIter, int activeGroupNum, bool actSetRemove){

  arma::vec old_beta = beta, r(n_obs), r_shift;
  arma::vec haz(n_obs), rsk(n_obs), h(n_obs); // (1) haz: exp(eta); (2) rsk: sum(eta); (3) h: delta/sum(eta);
  double loss, df, MaxChange_beta, shift, lambda_g;
  double s, v = 1; // (1) s: l'(eta); (2) v: l''(eta) qpproximate to 1
  int iter = 0;

  while (tol_iter < max_total_iter) {
    int inner_loop_iter = -1;
    inner_loop_iter = inner_loop_iter +1; //just for removing unused warning.
    R_CheckUserInterrupt();
    while (tol_iter < max_total_iter && iter < max_each_iter) {
      R_CheckUserInterrupt();
      df = 0;
      tol_iter++;
      iter++;
      inner_loop_iter++;
      MaxChange_beta = 0;


      // calculate haz and rsk
      haz = exp(eta);
      for (int j = 0; j < count_stratum; j++){
        rsk(ind_start(j) + n_each_prov(j) - 1) = haz(ind_start(j) + n_each_prov(j) - 1); //the last record (with max time) in j^th stratum
        for (int i = ind_start(j) + n_each_prov(j) - 2; i >= ind_start(j); i--){
          rsk(i) = rsk(i + 1) + haz(i);
        }
      }

      // calculate l'(eta)
      for (int j = 0; j < count_stratum; j++){
        h(ind_start(j)) = delta_obs(ind_start(j))/rsk(ind_start(j));
        for (int i = ind_start(j) + 1; i < ind_start(j) + n_each_prov(j); i++){
          h(i) = h(i - 1) + delta_obs(i)/rsk(i);
        }
      }

      double a; // a: exp(eta) * [sum(delta/sum(exp(eta)))]

      for (int i = 0; i < n_obs; i++){ //compute the initial pseudo residual r = l'(eta)
        a = haz(i) * h(i);
        s = delta_obs(i) - a; //l'(eta)
        if (a == 0){
          r(i) = 0;
        } else {
          r(i) = s/v;
        }
      }

      loss = 0; //loss: current likelihood function (before update beta) without penalty
      for (int i = 0; i < n_obs; i++){
         loss += delta_obs(i) *(haz(i) + log(rsk(i)));
      }
      
      //update unpenalized groups
      arma::uvec update_order_unpenalized = randperm(K0); //randomly update
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

      // Update penalized groups
      // note that all groups are iterated if user choose not using the "active set method"
      for (int g = 0; g < n_group; g++){
        if (active_group(g) == 1){
          lambda_g = lambda * group_multiplier(g);
          gd_stratCox(beta, Z, r, eta, old_beta, g, K1, n_obs, lambda_g, df, MaxChange_beta);
        }
      }
      old_beta = beta;

      if (MaxChange_beta < tol){
        break;
      }
    }

    if (actSet == true){
    // update active set
      if (actSetRemove == true){
        for (int g = 0; g < n_group; g++){
          if (active_group(g) == 1) {
            if (beta(K1(g)) == 0){
              active_group(g) = 0;
            }
          }
        }
      }
      // check whether "non-active" groups can still be updated; if not, algorithm ends
      arma::vec Current_len_group(n_group, fill::zeros);
      for (int g = 0; g < n_group; g++) {
        if (active_group(g) == 0) {
          double lambda_g = lambda * group_multiplier(g);
          Current_len_group(g) = gd_stratCox_BetaChange(Z, r, g, K1, n_obs, lambda_g);
        }
      }

      int if_add_new = 0; 
      arma::uvec descend_len_index = sort_index(Current_len_group, "descend");
      arma::vec descend_len_value = sort(Current_len_group, "descend");

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

  return make_tuple(beta, eta, loss, df, iter);

}

// [[Rcpp::export]]
List StratCox_lasso(arma::vec &delta_obs, arma::mat &Z, arma::vec &n_each_prov, arma::vec &beta, int K0, arma::vec &K1, 
                    arma::vec &lambda_seq, bool lambda_early_stop, double stop_loss_ratio, arma::vec &group_multiplier, 
                    int max_total_iter, int max_each_iter, double tol, int initial_active_group, double nvar_max, 
                    double group_max, bool trace_lambda, bool actSet, int actIter, int activeGroupNum, bool actSetRemove) {
  int n_obs = Z.n_rows, n_beta = Z.n_cols, n_lambda = lambda_seq.n_elem, n_group = K1.n_elem - 1;
  int tol_iter = 0;
  int count_stratum = n_each_prov.n_elem;

  arma::mat beta_matrix(n_beta, n_lambda, fill::zeros); 
  arma::mat eta_matrix(n_obs, n_lambda, fill::zeros); // linear predictor
  arma::vec iter_vec(n_lambda, fill::zeros);  // number of iterations for each lambda
  arma::vec df_vec(n_lambda, fill::zeros); 
  arma::vec loss_vec(n_lambda, fill::zeros); //loss: likelihood function (without penalty)
  arma::vec active_group(n_group, fill::zeros);

  if (actSet == true){
    if (K0 == 0){ 
      active_group(initial_active_group) = 1;
    }
  } else {
    active_group.ones();
  }

  arma::vec ind_start(count_stratum); // index of the first observation within each provider
  ind_start(0) = 0;
  for (int i = 1; i < count_stratum; i++) {
    ind_start(i) = ind_start(i - 1) + n_each_prov(i - 1);
  }
 
  // initialize eta
  arma::vec eta = Z * beta;

  for (int l = 0; l < n_lambda; l++){
    R_CheckUserInterrupt();
    if (trace_lambda == true){
      cout << "processing lambda: " << l + 1 << " (total: " << l + 1 << "/" << n_lambda << ")..." << endl;
    }
    double lambda = lambda_seq(l);
    
    auto fit = StratCox_lasso_fit(delta_obs, Z, n_each_prov, beta, eta, K0, K1, lambda, tol_iter, max_total_iter, max_each_iter, group_multiplier, count_stratum, tol, ind_start, active_group, n_obs, n_group, actSet, actIter, activeGroupNum, actSetRemove);
    double loss_l, df_l;
    int iter_l;
    tie(beta, eta, loss_l, df_l, iter_l) = fit;
    beta_matrix.col(l) = beta;
    eta_matrix.col(l) = eta;
    loss_vec(l) = loss_l;
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
      if (l != 0){
        double null_lkd = loss_vec(0);
        double loss_ratio = fabs((loss_vec(l) - loss_vec(l - 1))/(loss_vec(l) - null_lkd));
        if (loss_ratio < stop_loss_ratio){
          for (int ll = (l + 1); ll < n_lambda; ll++){
            iter_vec(ll) = NA_REAL;
          }
        break;
        }
      }
    }
  }

  List result = List::create(_["beta"] = beta_matrix, _["loss"] = loss_vec, _["Eta"] = eta_matrix, _["Df"] = df_vec, _["iter"] = iter_vec);
  return result;
}