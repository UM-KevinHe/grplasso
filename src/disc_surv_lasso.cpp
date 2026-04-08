#include "utils.h"

// ============================================================
//  Discrete Survival Model with Logit Link + LASSO
//  — Two variants:
//    (A) pp_DiscSurv: with provider/center effects (alpha)
//    (B) DiscSurv:    without provider effects
//  — Shared coordinate descent updater (gd_Surv)
// ============================================================


// ---------- Coordinate descent update for survival beta ----------
void gd_Surv(arma::vec &beta, arma::mat &Z, arma::vec &r, arma::vec &eta, arma::vec &time,
             arma::vec &old_beta, arma::vec &w, int p, int n_obs, double lambda, double &df, double &MaxChange_beta) {
  double beta_initial = w_crossprod(Z, r, w, p) / weighted_inner_product(Z, w, p) + old_beta(p);
  double threshold = n_obs * lambda / weighted_inner_product(Z, w, p);
  double len = Soft_thres(beta_initial, threshold);

  if (len != old_beta(p)) {
    beta(p) = len;
    double beta_change = beta(p) - old_beta(p);
    if (fabs(beta_change) > MaxChange_beta) {
      MaxChange_beta = fabs(beta_change);
    }
    if (beta_change != 0) {
      r -= beta_change * Z.col(p);
      eta += beta_change * Z.col(p);
    }
  }
  if (len != 0) {
    df += len / beta_initial;
  }
}

// ---------- BetaChange for active-set screening ----------
double gd_Surv_BetaChange(arma::mat &Z, arma::vec &r, arma::vec &eta, arma::vec &time,
                           arma::vec &w, int p, int n_obs, double lambda) {
  double beta_initial = w_crossprod(Z, r, w, p) / weighted_inner_product(Z, w, p);
  double threshold = n_obs * lambda / weighted_inner_product(Z, w, p);
  double len = Soft_thres(beta_initial, threshold);

  if (len != 0) {
    return(fabs(len));
  } else {
    return(0);
  }
}


// ============================================================
//  (A) pp_DiscSurv — with provider/center effects
// ============================================================

tuple<arma::vec, arma::vec, arma::vec, arma::vec, double, int> pp_DiscSurv_fit(
    arma::vec &delta_obs, int max_timepoint, arma::mat &Z, arma::vec &time, arma::vec &n_prov,
    arma::vec gamma, arma::vec beta, arma::vec alpha, arma::vec eta,
    int K0, arma::vec &K1, arma::vec &sum_failure, arma::vec failure_each_center,
    double lambda, int &tol_iter, int max_total_iter, int max_each_iter,
    arma::vec &penalized_multiplier, bool backtrack, bool MM, double bound, double tol,
    arma::vec &ind_start, arma::vec &active_var, int n_obs, int n_var, int threads,
    bool actSet, int actIter, int activeVarNum, bool actSetRemove) {

  arma::vec old_beta = beta, old_alpha = alpha, r(n_obs), r_shift, w(n_obs);
  int n_alpha = alpha.n_elem;
  double df, MaxChange_beta, shift;
  double v = 0.25, omega_min = 1e-10;
  double p_gamma;
  int iter = 0;

  while (tol_iter < max_total_iter && iter < max_each_iter) {
    int inner_loop_iter = 0;
    R_CheckUserInterrupt();
    while (tol_iter < max_total_iter && iter < max_each_iter && inner_loop_iter < actIter) {
      R_CheckUserInterrupt();
      df = 0;
      iter++;
      inner_loop_iter++;
      MaxChange_beta = 0;

      // 1. update gamma (time effect)
      arma::vec score_gamma = -sum_failure, info_gamma(max_timepoint, fill::zeros), d_gamma(max_timepoint, fill::zeros);
      if (backtrack == true) {
        arma::vec gamma_tmp, eta_tmp(n_obs);
        double loglkd, d_loglkd, u = 1.0, k, s = 0.01, t = 0.8;
        for (int j = 0; j < n_obs; j++) {
          for (int i = 0; i < time(j); i++) {
            p_gamma = p_binomial_Surv(gamma(i), eta(j));
            score_gamma(i) += p_gamma;
            info_gamma(i) -= p_gamma * (1 - p_gamma);
          }
        }
        int nProcessors = threads;
        omp_set_num_threads(nProcessors);
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < max_timepoint; i++) {
          double tmp_info_gamma = info_gamma(i);
          d_gamma(i) = score_gamma(i) / std::min(-omega_min, tmp_info_gamma);
        }
        loglkd = Loglkd_Surv(n_obs, delta_obs, time, gamma, eta);
        gamma_tmp = u * d_gamma + gamma;
        d_loglkd = Loglkd_Surv(n_obs, delta_obs, time, gamma_tmp, eta) - loglkd;
        k = dot(score_gamma, d_gamma);
        while (d_loglkd < s * u * k) {
          u = t * u;
          gamma_tmp = u * d_gamma + gamma;
          d_loglkd = Loglkd_Surv(n_obs, delta_obs, time, gamma_tmp, eta) - loglkd;
        }
        gamma = gamma + u * d_gamma;
        gamma = clamp(gamma, median(gamma) - bound, median(gamma) + bound);
      } else {
        for (int j = 0; j < n_obs; j++) {
          for (int i = 0; i < time(j); i++) {
            p_gamma = p_binomial_Surv(gamma(i), eta(j));
            score_gamma(i) += p_gamma;
            info_gamma(i) -= p_gamma * (1 - p_gamma);
          }
        }
        int nProcessors = threads;
        omp_set_num_threads(nProcessors);
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < max_timepoint; i++) {
          double tmp_info_gamma = info_gamma(i);
          d_gamma(i) = score_gamma(i) / std::min(-omega_min, tmp_info_gamma);
        }
        gamma = gamma + d_gamma;
        gamma = clamp(gamma, median(gamma) - bound, median(gamma) + bound);
      }

      // 2. update alpha (center effect)
      double temp_p;
      arma::vec p_obs(n_obs, fill::zeros);
      arma::vec pq_obs(n_obs, fill::zeros);
      for (int j = 0; j < n_obs; j++) {
        for (int i = 0; i < time(j); i++) {
          temp_p = p_binomial_Surv(gamma(i), eta(j));
          p_obs(j) += temp_p;
          pq_obs(j) += temp_p * (1 - temp_p);
        }
      }

      arma::vec score_alpha(n_alpha, fill::zeros), info_alpha(n_alpha, fill::zeros), d_alpha(n_alpha, fill::zeros);
      if (backtrack == true) {
        arma::vec alpha_shift_tmp, eta_tmp(n_obs);
        double loglkd, d_loglkd, u2 = 1.0, k2, s2 = 0.01, t2 = 0.8;
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
          d_alpha(i) = score_alpha(i) / info_alpha(i);
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
          d_alpha(i) = score_alpha(i) / info_alpha(i);
        }
        alpha = alpha + d_alpha;
        alpha = clamp(alpha, median(alpha) - bound, median(alpha) + bound);
        arma::vec alpha_shift = alpha - old_alpha;
        eta += rep(alpha_shift, n_prov);
      }

      // 3. update beta
      arma::vec p_eta(n_obs, fill::zeros);
      if (MM == true) {
        for (int j = 0; j < n_obs; j++) {
          for (int i = 0; i < time(j); i++) {
            p_eta(j) += p_binomial_Surv(gamma(i), eta(j));
          }
          w(j) = v * time(j);
        }
      } else {
        w.zeros();
        for (int j = 0; j < n_obs; j++) {
          for (int i = 0; i < time(j); i++) {
            double p_temp = p_binomial_Surv(gamma(i), eta(j));
            p_eta(j) += p_temp;
            w(j) += p_temp * (1 - p_temp);
          }
        }
        if (any(w == 0)) { w.replace(0, omega_min); }
      }
      arma::vec score_eta = -delta_obs + p_eta;
      r = -score_eta / w;

      // 3.1 update unpenalized beta
      arma::uvec update_order_unpenalized = randperm(K0);
      for (int p = 0; p < K0; p++) {
        shift = w_crossprod(Z, r, w, update_order_unpenalized(p)) / weighted_inner_product(Z, w, update_order_unpenalized(p));
        if (fabs(shift) > MaxChange_beta) { MaxChange_beta = fabs(shift); }
        beta(update_order_unpenalized(p)) = old_beta(update_order_unpenalized(p)) + shift;
        r -= Z.col(update_order_unpenalized(p)) * shift;
        eta += Z.col(update_order_unpenalized(p)) * shift;
        df++;
      }

      // 3.2 update penalized beta
      arma::uvec update_order_penalized = randperm(n_var);
      for (int p = 0; p < n_var; p++) {
        int update_column_index = K1(update_order_penalized(p));
        if (active_var(update_order_penalized(p)) == 1) {
          double lambda_m = lambda * penalized_multiplier(update_order_penalized(p));
          gd_Surv(beta, Z, r, eta, time, old_beta, w, update_column_index, n_obs, lambda_m, df, MaxChange_beta);
        }
      }

      old_beta = beta;
      old_alpha = alpha;
      if (MaxChange_beta < tol) { break; }
    }

    if (actSet == true) {
      if (actSetRemove == true) {
        for (int p = 0; p < n_var; p++) {
          if (active_var(p) == 1) {
            if (beta(p) == 0) { active_var(p) = 0; }
          }
        }
      }

      arma::vec Current_Change_beta(n_var, fill::zeros);
      for (int p = 0; p < n_var; p++) {
        if (active_var(p) == 0) {
          double lambda_m = lambda * penalized_multiplier(p);
          Current_Change_beta(p) = gd_Surv_BetaChange(Z, r, eta, time, w, p, n_obs, lambda_m);
        }
      }

      int if_add_new = 0;
      arma::uvec descend_beta_change_index = sort_index(Current_Change_beta, "descend");
      arma::vec descend_beta_change = sort(Current_Change_beta, "descend");

      for (int i = 0; i < activeVarNum; i++) {
        if (descend_beta_change(i) != 0) {
          if_add_new++;
          active_var(descend_beta_change_index(i)) = 1;
        } else {
          break;
        }
      }

      if (if_add_new == 0) { break; }
    } else {
      break;
    }
  }

  return make_tuple(beta, gamma, alpha, eta, df, iter);
}


// [[Rcpp::export]]
List pp_DiscSurv_lasso(arma::vec &delta_obs, int max_timepoint, arma::mat &Z, arma::vec &n_prov, arma::vec &time,
                       arma::vec &gamma, arma::vec &beta, arma::vec &alpha,
                       int K0, arma::vec &K1, arma::vec &sum_failure, arma::vec failure_each_center,
                       arma::vec &lambda_seq, arma::vec &penalized_multiplier,
                       int max_total_iter, int max_each_iter, double tol, bool backtrack, bool MM, double bound,
                       int initial_active_var, double nvar_max, bool trace_lambda, int threads,
                       bool actSet, int actIter, int activeVarNum, bool actSetRemove) {

  int n_obs = Z.n_rows, n_beta = Z.n_cols, n_gamma = gamma.n_elem, n_alpha = alpha.n_elem, n_lambda = lambda_seq.n_elem;
  int n_var = K1.n_elem - 1;
  int tol_iter = 0;

  arma::mat beta_matrix(n_beta, n_lambda, fill::zeros), gamma_matrix(n_gamma, n_lambda, fill::zeros), alpha_matrix(n_alpha, n_lambda, fill::zeros);
  arma::mat eta_matrix(n_obs, n_lambda, fill::zeros);
  arma::vec iter_vec(n_lambda, fill::zeros);
  arma::vec df_vec(n_lambda, fill::zeros);
  arma::vec active_var(n_var, fill::zeros);

  if (actSet == true) {
    if (K0 == 0) { active_var(initial_active_var) = 1; }
  } else {
    active_var.ones();
  }

  arma::vec ind_start(n_alpha);
  ind_start(0) = 0;
  for (int i = 1; i < n_alpha; i++) {
    ind_start(i) = ind_start(i - 1) + n_prov(i - 1);
  }

  arma::vec eta = rep(alpha, n_prov) + Z * beta;

  for (int l = 0; l < n_lambda; l++) {
    R_CheckUserInterrupt();
    if (trace_lambda == true) {
      cout << "processing lambda: " << l + 1 << " (total: " << l + 1 << "/" << n_lambda << ")..." << endl;
    }
    double lambda = lambda_seq(l);

    auto fit = pp_DiscSurv_fit(delta_obs, max_timepoint, Z, time, n_prov, gamma, beta, alpha, eta, K0, K1,
                               sum_failure, failure_each_center, lambda, tol_iter, max_total_iter, max_each_iter,
                               penalized_multiplier, backtrack, MM, bound, tol, ind_start, active_var,
                               n_obs, n_var, threads, actSet, actIter, activeVarNum, actSetRemove);

    double df_l;
    int iter_l;
    tie(beta, gamma, alpha, eta, df_l, iter_l) = fit;
    beta_matrix.col(l) = beta;
    gamma_matrix.col(l) = gamma;
    alpha_matrix.col(l) = alpha;
    eta_matrix.col(l) = eta;
    df_vec(l) = df_l;
    iter_vec(l) = iter_l;

    if (iter_l == max_each_iter) {
      cout << "Warning: lambda " << l + 1 << "/" << n_lambda << " failed to converge within " << max_each_iter << " iterations!" << endl;
    }

    int nv = 0;
    for (int j = 0; j < n_var; j++) {
      if (beta(K1(j)) != 0) { nv++; }
    }

    if (nv > nvar_max || tol_iter == max_total_iter) {
      if (tol_iter == max_total_iter) {
        cout << "Algorithm has reached the maximum number of total iterations, stops..." << endl;
      } else {
        cout << "Algorithm has selected the maximum number of penalized variables, stops..." << endl;
      }
      for (int ll = (l + 1); ll < n_lambda; ll++) { iter_vec(ll) = NA_REAL; }
      break;
    }
  }

  List result = List::create(_["gamma"] = gamma_matrix, _["beta"] = beta_matrix, _["alpha"] = alpha_matrix,
                             _["Eta"] = eta_matrix, _["Df"] = df_vec, _["iter"] = iter_vec);
  return result;
}


// ============================================================
//  (B) DiscSurv — without provider effects
// ============================================================

tuple<vec, vec, vec, double, int> DiscSurv_fit(
    arma::vec &delta_obs, int max_timepoint, arma::mat &Z, arma::vec &time,
    arma::vec gamma, arma::vec beta, arma::vec eta,
    int K0, arma::vec &K1, arma::vec &sum_failure,
    double lambda, int &tol_iter, int max_total_iter, int max_each_iter,
    arma::vec &penalized_multiplier, bool backtrack, bool MM, double bound, double tol,
    arma::vec &active_var, int n_obs, int n_var, int threads,
    bool actSet, int actIter, int activeVarNum, bool actSetRemove) {

  arma::vec old_beta = beta, r(n_obs), r_shift, w(n_obs);
  double df, MaxChange_beta, shift;
  double v = 0.25, omega_min = 1e-10;
  double p_gamma;
  int iter = 0;

  while (tol_iter < max_total_iter && iter < max_each_iter) {
    int inner_loop_iter = 0;
    R_CheckUserInterrupt();
    while (tol_iter < max_total_iter && iter < max_each_iter && inner_loop_iter < actIter) {
      R_CheckUserInterrupt();
      df = 0;
      iter++;
      inner_loop_iter++;
      MaxChange_beta = 0;

      // 1. update gamma (time effect)
      arma::vec score_gamma = -sum_failure, info_gamma(max_timepoint, fill::zeros), d_gamma(max_timepoint, fill::zeros);
      if (backtrack == true) {
        arma::vec gamma_tmp, eta_tmp(n_obs);
        double loglkd, d_loglkd, u = 1.0, k, s = 0.01, t = 0.8;
        for (int j = 0; j < n_obs; j++) {
          for (int i = 0; i < time(j); i++) {
            p_gamma = p_binomial_Surv(gamma(i), eta(j));
            score_gamma(i) += p_gamma;
            info_gamma(i) -= p_gamma * (1 - p_gamma);
          }
        }
        int nProcessors = threads;
        omp_set_num_threads(nProcessors);
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < max_timepoint; i++) {
          double tmp_info_gamma = info_gamma(i);
          d_gamma(i) = score_gamma(i) / std::min(-omega_min, tmp_info_gamma);
        }
        loglkd = Loglkd_Surv(n_obs, delta_obs, time, gamma, eta);
        gamma_tmp = u * d_gamma + gamma;
        d_loglkd = Loglkd_Surv(n_obs, delta_obs, time, gamma_tmp, eta) - loglkd;
        k = dot(score_gamma, d_gamma);
        while (d_loglkd < s * u * k) {
          u = t * u;
          gamma_tmp = u * d_gamma + gamma;
          d_loglkd = Loglkd_Surv(n_obs, delta_obs, time, gamma_tmp, eta) - loglkd;
        }
        gamma = gamma + u * d_gamma;
        gamma = clamp(gamma, median(gamma) - bound, median(gamma) + bound);
      } else {
        for (int j = 0; j < n_obs; j++) {
          for (int i = 0; i < time(j); i++) {
            p_gamma = p_binomial_Surv(gamma(i), eta(j));
            score_gamma(i) += p_gamma;
            info_gamma(i) -= p_gamma * (1 - p_gamma);
          }
        }
        int nProcessors = threads;
        omp_set_num_threads(nProcessors);
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < max_timepoint; i++) {
          double tmp_info_gamma = info_gamma(i);
          d_gamma(i) = score_gamma(i) / std::min(-omega_min, tmp_info_gamma);
        }
        gamma = gamma + d_gamma;
        gamma = clamp(gamma, median(gamma) - bound, median(gamma) + bound);
      }

      // 2. update beta
      arma::vec p_eta(n_obs, fill::zeros);
      if (MM == true) {
        for (int j = 0; j < n_obs; j++) {
          for (int i = 0; i < time(j); i++) {
            p_eta(j) += p_binomial_Surv(gamma(i), eta(j));
          }
          w(j) = v * time(j);
        }
      } else {
        w.zeros();
        for (int j = 0; j < n_obs; j++) {
          for (int i = 0; i < time(j); i++) {
            double p_temp = p_binomial_Surv(gamma(i), eta(j));
            p_eta(j) += p_temp;
            w(j) += p_temp * (1 - p_temp);
          }
        }
        if (any(w == 0)) { w.replace(0, omega_min); }
      }
      arma::vec score_eta = p_eta - delta_obs;
      r = -score_eta / w;

      // 2.1 update unpenalized beta
      arma::uvec update_order_unpenalized = randperm(K0);
      for (int p = 0; p < K0; p++) {
        shift = w_crossprod(Z, r, w, update_order_unpenalized(p)) / weighted_inner_product(Z, w, update_order_unpenalized(p));
        if (fabs(shift) > MaxChange_beta) { MaxChange_beta = fabs(shift); }
        beta(update_order_unpenalized(p)) = old_beta(update_order_unpenalized(p)) + shift;
        r -= Z.col(update_order_unpenalized(p)) * shift;
        eta += Z.col(update_order_unpenalized(p)) * shift;
        df++;
      }

      // 2.2 update penalized beta
      arma::uvec update_order_penalized = randperm(n_var);
      for (int p = 0; p < n_var; p++) {
        int update_column_index = K1(update_order_penalized(p));
        if (active_var(update_order_penalized(p)) == 1) {
          double lambda_m = lambda * penalized_multiplier(update_order_penalized(p));
          gd_Surv(beta, Z, r, eta, time, old_beta, w, update_column_index, n_obs, lambda_m, df, MaxChange_beta);
        }
      }

      old_beta = beta;
      if (MaxChange_beta < tol) { break; }
    }

    if (actSet == true) {
      if (actSetRemove == true) {
        for (int p = 0; p < n_var; p++) {
          if (active_var(p) == 1) {
            if (beta(p) == 0) { active_var(p) = 0; }
          }
        }
      }

      arma::vec Current_Change_beta(n_var, fill::zeros);
      for (int p = 0; p < n_var; p++) {
        if (active_var(p) == 0) {
          double lambda_m = lambda * penalized_multiplier(p);
          Current_Change_beta(p) = gd_Surv_BetaChange(Z, r, eta, time, w, p, n_obs, lambda_m);
        }
      }

      int if_add_new = 0;
      arma::uvec descend_beta_change_index = sort_index(Current_Change_beta, "descend");
      arma::vec descend_beta_change = sort(Current_Change_beta, "descend");

      for (int i = 0; i < activeVarNum; i++) {
        if (descend_beta_change(i) != 0) {
          if_add_new++;
          active_var(descend_beta_change_index(i)) = 1;
        } else {
          break;
        }
      }

      if (if_add_new == 0) { break; }
    } else {
      break;
    }
  }

  return make_tuple(beta, gamma, eta, df, iter);
}


// [[Rcpp::export]]
List DiscSurv_lasso(arma::vec &delta_obs, int max_timepoint, arma::mat &Z, arma::vec &time,
                    arma::vec &gamma, arma::vec &beta,
                    int K0, arma::vec &K1, arma::vec &sum_failure, arma::vec &lambda_seq,
                    arma::vec &penalized_multiplier, int max_total_iter, int max_each_iter, double tol,
                    bool backtrack, bool MM, double bound, int initial_active_var, double nvar_max,
                    bool trace_lambda, int threads, bool actSet, int actIter, int activeVarNum, bool actSetRemove) {

  int n_obs = Z.n_rows, n_beta = Z.n_cols, n_gamma = gamma.n_elem, n_lambda = lambda_seq.n_elem;
  int n_var = K1.n_elem - 1;
  int tol_iter = 0;

  arma::mat beta_matrix(n_beta, n_lambda, fill::zeros), gamma_matrix(n_gamma, n_lambda, fill::zeros);
  arma::mat eta_matrix(n_obs, n_lambda, fill::zeros);
  arma::vec iter_vec(n_lambda, fill::zeros);
  arma::vec df_vec(n_lambda, fill::zeros);
  arma::vec active_var(n_var, fill::zeros);

  if (actSet == true) {
    if (K0 == 0) { active_var(initial_active_var) = 1; }
  } else {
    active_var.ones();
  }

  arma::vec eta = Z * beta;

  for (int l = 0; l < n_lambda; l++) {
    R_CheckUserInterrupt();
    if (trace_lambda == true) {
      cout << "processing lambda: " << l + 1 << " (total: " << l + 1 << "/" << n_lambda << ")..." << endl;
    }
    double lambda = lambda_seq(l);

    auto fit = DiscSurv_fit(delta_obs, max_timepoint, Z, time, gamma, beta, eta, K0, K1, sum_failure, lambda,
                            tol_iter, max_total_iter, max_each_iter, penalized_multiplier, backtrack, MM, bound, tol,
                            active_var, n_obs, n_var, threads, actSet, actIter, activeVarNum, actSetRemove);

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
    for (int j = 0; j < n_var; j++) {
      if (beta(K1(j)) != 0) { nv++; }
    }

    if (nv > nvar_max || tol_iter == max_total_iter) {
      if (tol_iter == max_total_iter) {
        cout << "Algorithm has reached the maximum number of total iterations, stops..." << endl;
      } else {
        cout << "Algorithm has selected the maximum number of penalized variables, stops..." << endl;
      }
      for (int ll = (l + 1); ll < n_lambda; ll++) { iter_vec(ll) = NA_REAL; }
      break;
    }
  }

  List result = List::create(_["gamma"] = gamma_matrix, _["beta"] = beta_matrix,
                             _["Eta"] = eta_matrix, _["Df"] = df_vec, _["iter"] = iter_vec);
  return result;
}


// ============================================================
//  Residuals & Prediction helpers (Rcpp exported)
// ============================================================

// [[Rcpp::export]]
arma::vec DiscSurv_residuals(int n_obs, arma::vec &delta_obs, arma::vec &time, arma::vec &alpha, arma::vec &eta) {
  arma::vec residuals(n_obs);
  for (int j = 0; j < n_obs; j++) {
    for (int i = 0; i < time(j); i++) {
      residuals(j) += p_binomial_Surv(alpha(i), eta(j));
    }
  }
  residuals = -residuals + delta_obs;
  return(residuals);
}

// [[Rcpp::export]]
arma::mat predict_linear_predictor(int n_lambda, int n_obs, int expand_n_obs, arma::vec &time, arma::mat &gamma, arma::mat &eta) {
  arma::mat predict_linear_predictor(expand_n_obs, n_lambda);
  for (int l = 0; l < n_lambda; l++) {
    int sum_index = 0;
    for (int j = 0; j < n_obs; j++) {
      for (int i = 0; i < time(j); i++) {
        predict_linear_predictor(sum_index, l) = p_binomial_Surv(gamma(i, l), eta(j, l));
        sum_index += 1;
      }
    }
  }
  return(predict_linear_predictor);
}
