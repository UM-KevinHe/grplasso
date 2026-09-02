#include "utils.h"

// ============================================================
//  Penalized Provider Lasso (pp_lasso)
//  — Binomial GLM with provider-specific intercepts
//  — Individual-level LASSO penalty (non-group)
//  — Supports MM and Newton coordinate descent
// ============================================================


// ---------- Coordinate descent update: MM algorithm ----------
void gd_glm_lasso_MM(arma::vec &beta, arma::mat &Z, arma::vec &r, arma::vec &eta, arma::vec &old_beta,
                     arma::vec &inner_product_Z, int j, int n_obs, double lambda, double &df, double &MaxChange_beta) {
  double beta_initial = crossprod(Z, r, j) / inner_product_Z(j) + old_beta(j);
  double threshold = n_obs * lambda / (0.25 * inner_product_Z(j));
  double len = Soft_thres(beta_initial, threshold);

  if (len != old_beta(j)) {
    beta(j) = len;
    double beta_change = beta(j) - old_beta(j);
    if (fabs(beta_change) > MaxChange_beta) {
      MaxChange_beta = fabs(beta_change);
    }
    if (beta_change != 0) {
      r -= beta_change * Z.col(j);
      eta += beta_change * Z.col(j);
    }
  }
  if (len != 0) {
    df += len / beta_initial;
  }
}

// ---------- Coordinate descent update: without MM ----------
void gd_glm_lasso_noMM(arma::vec &beta, arma::vec &w, arma::mat &Z, arma::vec &r, arma::vec &eta, arma::vec &old_beta,
                        int j, int n_obs, double lambda, double &df, double &MaxChange_beta) {
  double w_z_square = weighted_inner_product(Z, w, j);
  double beta_initial = w_crossprod(Z, r, w, j) / w_z_square + old_beta(j);
  double threshold = n_obs * lambda / w_z_square;
  double len = Soft_thres(beta_initial, threshold);

  if (len != old_beta(j)) {
    beta(j) = len;
    double beta_change = beta(j) - old_beta(j);
    if (fabs(beta_change) > MaxChange_beta) {
      MaxChange_beta = fabs(beta_change);
    }
    if (beta_change != 0) {
      r -= beta_change * Z.col(j);
      eta += beta_change * Z.col(j);
    }
  }
  if (len != 0) {
    df += len / beta_initial;
  }
}

// ---------- BetaChange for active-set screening: MM ----------
double gd_glm_lasso_MM_BetaChange(arma::mat &Z, arma::vec &r, arma::vec &inner_product_Z, int j, int n_obs, double lambda) {
  double beta_initial = crossprod(Z, r, j) / inner_product_Z(j);
  double threshold = n_obs * lambda / (0.25 * inner_product_Z(j));
  double len = Soft_thres(beta_initial, threshold);

  if (len != 0) {
    return(fabs(len));
  } else {
    return(0);
  }
}

// ---------- BetaChange for active-set screening: without MM ----------
double gd_glm_lasso_noMM_BetaChange(arma::mat &Z, arma::vec &r, arma::vec &eta, int j, int n_obs, double lambda) {
  arma::vec p(n_obs);
  for (int i = 0; i < n_obs; i++) {
    p(i) = p_binomial(eta(i));
  }
  arma::vec w = p % (1 - p);
  double beta_initial = w_crossprod(Z, r, w, j) / weighted_inner_product(Z, w, j);
  double threshold = n_obs * lambda / weighted_inner_product(Z, w, j);
  double len = Soft_thres(beta_initial, threshold);

  if (len != 0) {
    return(fabs(len));
  } else {
    return(0);
  }
}


// ---------- Single-lambda fit ----------
tuple<arma::vec, arma::vec, arma::vec, double, double, int> pp_lasso_fit(
    arma::vec &Y, arma::mat &Z, arma::vec &n_prov, arma::vec gamma, arma::vec beta, arma::vec eta,
    int K0, arma::vec &K1, double lambda, int &tol_iter, int max_total_iter, int max_each_iter,
    arma::vec &penalized_multiplier, int max_n_prov, bool backtrack, bool MM, double bound, double tol,
    arma::vec &ind_start, arma::vec &active_var, int n_obs, int n_var, bool single_intercept, int threads,
    bool actSet, int actIter, int activeVarNum, bool actSetRemove) {

  arma::vec old_beta = beta, old_gamma = gamma, p(n_obs), r(n_obs), r_shift;
  int n_gamma = gamma.n_elem;
  double Dev, df, MaxChange_beta, shift;
  double v = 0.25, omega_min = 1e-20;
  int iter = 0;
  arma::vec inner_product_Z = inner_product(Z), w(n_obs);

  while (tol_iter < max_total_iter && iter < max_each_iter) {
    int inner_loop_iter = 0;
    R_CheckUserInterrupt();

    while (tol_iter < max_total_iter && iter < max_each_iter && inner_loop_iter < actIter) {
      R_CheckUserInterrupt();
      df = 0;
      tol_iter++;
      iter++;
      inner_loop_iter++;
      MaxChange_beta = 0;

      for (int i = 0; i < n_obs; i++) {
        p(i) = p_binomial(eta(i));
      }

      if (single_intercept == true) {
        if (MM == true) {
          r = (Y - p) / v;
          shift = sum(r) / n_obs;
          gamma += shift;
          eta += shift;
          r -= shift;
        } else {
          w = p % (1 - p);
          if (any(w == 0)) { w.replace(0, omega_min); }
          r = (Y - p) / w;
          shift = vec_crossprod(w, r) / sum(w);
          gamma += shift;
          eta += shift;
          r -= shift;
        }
      } else {
        double info_gamma;
        int nProcessors = threads;
        arma::vec score_gamma(n_gamma), d_gamma(n_gamma), Yp(n_obs), pq(n_obs), gamma_obs(n_obs);
        if (backtrack == true) {
          arma::vec gamma_shift_tmp, eta_tmp(n_obs);
          double loglkd, d_loglkd, u, k, s = 0.01, t = 0.8;
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
          u = 1.0;
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

        for (int i = 0; i < n_obs; i++) {
          p(i) = p_binomial(eta(i));
        }

        if (MM == true) {
          r = (Y - p) / v;
        } else {
          w = p % (1 - p);
          if (any(w == 0)) { w.replace(0, omega_min); }
          r = (Y - p) / w;
        }
      }

      // update beta
      arma::uvec update_order_unpenalized = randperm(K0);
      arma::uvec update_order_penalized = randperm(n_var);
      if (MM == true) {
        for (int j = 0; j < K0; j++) {
          shift = crossprod(Z, r, update_order_unpenalized(j)) / inner_product_Z(update_order_unpenalized(j));
          if (fabs(shift) > MaxChange_beta) { MaxChange_beta = fabs(shift); }
          beta(update_order_unpenalized(j)) = old_beta(update_order_unpenalized(j)) + shift;
          r -= Z.col(update_order_unpenalized(j)) * shift;
          eta += Z.col(update_order_unpenalized(j)) * shift;
          df++;
        }
        for (int j = 0; j < n_var; j++) {
          int update_column_index = K1(update_order_penalized(j));
          if (active_var(update_order_penalized(j)) == 1) {
            double lambda_m = lambda * penalized_multiplier(update_order_penalized(j));
            gd_glm_lasso_MM(beta, Z, r, eta, old_beta, inner_product_Z, update_column_index, n_obs, lambda_m, df, MaxChange_beta);
          }
        }
      } else {
        for (int j = 0; j < K0; j++) {
          shift = w_crossprod(Z, r, w, update_order_unpenalized(j)) / weighted_inner_product(Z, w, update_order_unpenalized(j));
          if (fabs(shift) > MaxChange_beta) { MaxChange_beta = fabs(shift); }
          beta(update_order_unpenalized(j)) = old_beta(update_order_unpenalized(j)) + shift;
          r -= Z.col(update_order_unpenalized(j)) * shift;
          eta += Z.col(update_order_unpenalized(j)) * shift;
          df++;
        }
        for (int j = 0; j < n_var; j++) {
          int update_column_index = K1(update_order_penalized(j));
          if (active_var(update_order_penalized(j)) == 1) {
            double lambda_m = lambda * penalized_multiplier(update_order_penalized(j));
            gd_glm_lasso_noMM(beta, w, Z, r, eta, old_beta, update_column_index, n_obs, lambda_m, df, MaxChange_beta);
          }
        }
      }

      old_gamma = gamma;
      old_beta = beta;

      if (MaxChange_beta < tol) {
        break;
      }
    }

    if (actSet == true) {
      if (actSetRemove == true) {
        for (int j = 0; j < n_var; j++) {
          if (active_var(j) == 1) {
            if (beta(j) == 0) { active_var(j) = 0; }
          }
        }
      }

      arma::vec Current_Change_beta(n_var, fill::zeros);

      if (MM == true) {
        for (int j = 0; j < n_var; j++) {
          if (active_var(j) == 0) {
            double lambda_m = lambda * penalized_multiplier(j);
            Current_Change_beta(j) = gd_glm_lasso_MM_BetaChange(Z, r, inner_product_Z, j, n_obs, lambda_m);
          }
        }
      } else {
        for (int j = 0; j < n_var; j++) {
          if (active_var(j) == 0) {
            double lambda_m = lambda * penalized_multiplier(j);
            Current_Change_beta(j) = gd_glm_lasso_noMM_BetaChange(Z, r, eta, j, n_obs, lambda_m);
          }
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

      if (if_add_new == 0) {
        break;
      }
    } else {
      break;
    }
  }

  for (int i = 0; i < n_obs; i++) {
    p(i) = p_binomial(eta(i));
  }
  Dev = Deviance(Y, p);
  return make_tuple(beta, gamma, eta, Dev, df, iter);
}


// ---------- Full lambda-path wrapper (Rcpp exported) ----------
// [[Rcpp::export]]
List pp_lasso(arma::vec &Y, arma::mat &Z, arma::vec &n_prov, arma::vec &gamma, arma::vec &beta, int K0, arma::vec &K1,
              arma::vec &lambda_seq, bool lambda_early_stop, double stop_dev_ratio, arma::vec &penalized_multiplier,
              int max_total_iter, int max_each_iter, double tol, double nullDev,
              bool backtrack, bool MM, double bound, int initial_active_var, double nvar_max, bool trace_lambda,
              bool single_intercept, int threads, bool actSet, int actIter, int activeVarNum, bool actSetRemove) {

  int n_obs = Z.n_rows, n_beta = Z.n_cols, n_gamma = n_prov.n_elem, n_lambda = lambda_seq.n_elem, max_n_prov = max(n_prov);
  int n_var = K1.n_elem - 1;
  int tol_iter = 0;

  arma::mat beta_matrix(n_beta, n_lambda, fill::zeros), gamma_matrix(n_gamma, n_lambda, fill::zeros);
  arma::mat eta_matrix(n_obs, n_lambda, fill::zeros);
  arma::vec Dev_vec(n_lambda, fill::zeros);
  arma::vec iter_vec(n_lambda, fill::zeros);
  arma::vec df_vec(n_lambda, fill::zeros);
  arma::vec active_var(n_var, fill::zeros);

  if (actSet == true) {
    if (K0 == 0) { active_var(initial_active_var) = 1; }
  } else {
    active_var.ones();
  }

  arma::vec ind_start(n_gamma);
  ind_start(0) = 0;
  for (int i = 1; i < n_gamma; i++) {
    ind_start(i) = ind_start(i - 1) + n_prov(i - 1);
  }

  arma::vec eta = rep(gamma, n_prov) + Z * beta;

  for (int l = 0; l < n_lambda; l++) {
    R_CheckUserInterrupt();
    if (trace_lambda == true) {
      Rcpp::Rcout << "processing lambda: " << l + 1 << " (total: " << l + 1 << "/" << n_lambda << ")..." << endl;
    }
    double lambda = lambda_seq(l);
    auto fit = pp_lasso_fit(Y, Z, n_prov, gamma, beta, eta, K0, K1, lambda, tol_iter, max_total_iter, max_each_iter,
                            penalized_multiplier, max_n_prov, backtrack, MM, bound, tol, ind_start, active_var,
                            n_obs, n_var, single_intercept, threads, actSet, actIter, activeVarNum, actSetRemove);
    double Dev_l, df_l;
    int iter_l;
    tie(beta, gamma, eta, Dev_l, df_l, iter_l) = fit;
    beta_matrix.col(l) = beta;
    gamma_matrix.col(l) = gamma;
    eta_matrix.col(l) = eta;
    Dev_vec(l) = Dev_l;
    df_vec(l) = df_l;
    iter_vec(l) = iter_l;

    if (iter_l == max_each_iter) {
      Rcpp::Rcout << "Warning: lambda " << l + 1 << "/" << n_lambda << " failed to converge within " << max_each_iter << " iterations!" << endl;
    }

    int nv = 0;
    for (int j = 0; j < n_var; j++) {
      if (beta(K1(j)) != 0) { nv++; }
    }
    if (nv > nvar_max || tol_iter == max_total_iter) {
      if (tol_iter == max_total_iter) {
        Rcpp::Rcout << "Algorithm has reached the maximum number of total iterations, stops..." << endl;
      } else {
        Rcpp::Rcout << "Algorithm has selected the maximum number of penalized variables, stops..." << endl;
      }
      for (int ll = (l + 1); ll < n_lambda; ll++) { iter_vec(ll) = NA_REAL; }
      break;
    }

    if (lambda_early_stop == true) {
      if (l != 0) {
        double Dev_ratio = (Dev_vec(l) - Dev_vec(l - 1)) / (Dev_vec(l) - nullDev);
        if (Dev_ratio < stop_dev_ratio) {
          for (int ll = (l + 1); ll < n_lambda; ll++) { iter_vec(ll) = NA_REAL; }
          break;
        }
      }
    }
  }

  List result = List::create(_["gamma"] = gamma_matrix, _["beta"] = beta_matrix, _["Deviance"] = Dev_vec,
                             _["Eta"] = eta_matrix, _["Df"] = df_vec, _["iter"] = iter_vec);
  return result;
}
