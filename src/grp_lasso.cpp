#include "utils.h"

// ============================================================
//  Group Lasso (grp_lasso)
//  — Binomial GLM with provider-specific intercepts
//  — Group-level LASSO penalty (MM algorithm)
// ============================================================


// ---------- Group descent update ----------
void gd_glm_grplasso(arma::vec &beta, arma::mat &Z, arma::vec &r, arma::vec &eta, arma::vec &old_beta,
                     int g, arma::vec &K1, int n_obs, double lambda, double &df, double &MaxChange_beta) {
  int K = K1(g + 1) - K1(g);
  arma::vec beta_initial(K);
  for (int j = K1(g); j < K1(g + 1); j++) {
    beta_initial(j - K1(g)) = mean_crossprod(Z, r, j, n_obs) + old_beta(j);
  }
  double beta_initial_norm = arma::norm(beta_initial, 2);
  double len = Soft_thres(beta_initial_norm, lambda / 0.25);

  if (len != 0 || old_beta(K1(g)) != 0) {
    for (int j = K1(g); j < K1(g + 1); j++) {
      beta(j) = len * beta_initial(j - K1(g)) / beta_initial_norm;
      double beta_change = beta(j) - old_beta(j);
      if (fabs(beta_change) > MaxChange_beta) {
        MaxChange_beta = fabs(beta_change);
      }
      r -= beta_change * Z.col(j);
      eta += beta_change * Z.col(j);
    }
  }
  if (len > 0) {
    df += K * len / beta_initial_norm;
  }
}

// ---------- BetaChange for active-set screening ----------
double gd_glm_grplasso_BetaChange(arma::mat &Z, arma::vec &r, int g, arma::vec &K1, int n_obs, double lambda) {
  int K = K1(g + 1) - K1(g);
  arma::vec beta_initial(K);
  for (int j = K1(g); j < K1(g + 1); j++) {
    beta_initial(j - K1(g)) = mean_crossprod(Z, r, j, n_obs);
  }
  double beta_initial_norm = arma::norm(beta_initial, 2);
  double len = Soft_thres(beta_initial_norm, lambda / 0.25);

  if (len != 0) {
    return(len);
  } else {
    return(0);
  }
}


// ---------- Single-lambda fit ----------
// *** MODIFIED: added bool firth parameter ***
tuple<arma::vec, arma::vec, arma::vec, double, double, int> grp_lasso_fit(
    arma::vec &Y, arma::mat &Z, arma::vec &n_prov, arma::vec gamma, arma::vec beta, arma::vec eta,
    int K0, arma::vec &K1, double lambda, int &tol_iter, int max_total_iter, int max_each_iter,
    arma::vec &group_multiplier, int max_n_prov, bool backtrack, double bound, double tol,
    arma::vec &ind_start, arma::vec &active_group, int n_obs, int n_group, bool single_intercept,
    int threads, bool actSet, int actIter, int activeGroupNum, bool actSetRemove,
    bool firth) {  // *** MODIFIED: added firth ***

  arma::vec old_beta = beta, old_gamma = gamma, p(n_obs), r(n_obs), r_shift;
  int n_gamma = gamma.n_elem;
  double Dev, df, MaxChange_beta, shift, lambda_g, v = 0.25;
  int iter = 0;

  while (tol_iter < max_total_iter && iter < max_each_iter) {
    int inner_loop_iter = -1;
    inner_loop_iter = inner_loop_iter + 1;
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
        r = (Y - p) / v;
        shift = sum(r) / n_obs;
        gamma += shift;
        eta += shift;
        r -= shift;
      } else {
        double omega_min = 1e-10, info_gamma;
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
            // *** MODIFIED: Firth correction for score and information ***
            if (firth) {
              arma::vec pq_i = pq(span(ind_start(i), ind_start(i) + n_prov(i) - 1));
              arma::vec p_i = p(span(ind_start(i), ind_start(i) + n_prov(i) - 1));
              double W_i = info_gamma;  // sum of w_ij = sum of pq_i
              arma::vec d_i = 1 - 2 * p_i;  // d_ij = 1 - 2*mu_ij
              double A_i = dot(pq_i, d_i);  // sum of w_ij * d_ij
              // U_i* = score + 0.5 * A_i / W_i
              score_gamma(i) += 0.5 * A_i / W_i;
              // J_i* = W_i - 0.5 * (sum(w_ij*(d_ij^2 - 2*w_ij))/W_i - A_i^2/W_i^2)
              double B_i = dot(pq_i, d_i % d_i - 2 * pq_i);  // sum of w_ij*(d_ij^2 - 2*w_ij)
              info_gamma = W_i - 0.5 * (B_i / W_i - A_i * A_i / (W_i * W_i));
            }
            // *** END MODIFIED ***
            info_gamma = std::max(omega_min, std::min(info_gamma, 0.25 * max_n_prov));
            d_gamma(i) = score_gamma(i) / info_gamma;
          }
          u = 1.0;
          // *** MODIFIED: use penalized log-likelihood for backtracking when firth ***
          loglkd = Loglkd(Y, eta);
          if (firth) {
            loglkd += FirthPenalty(p, n_prov, ind_start, n_gamma);
          }
          // *** END MODIFIED ***
          gamma_shift_tmp = u * d_gamma;
          eta_tmp = eta + rep(gamma_shift_tmp, n_prov);
          // *** MODIFIED: use penalized log-likelihood for backtracking when firth ***
          d_loglkd = Loglkd(Y, eta_tmp);
          if (firth) {
            arma::vec p_tmp(n_obs);
            for (int ii = 0; ii < n_obs; ii++) { p_tmp(ii) = p_binomial(eta_tmp(ii)); }
            d_loglkd += FirthPenalty(p_tmp, n_prov, ind_start, n_gamma);
          }
          d_loglkd -= loglkd;
          // *** END MODIFIED ***
          k = dot(score_gamma, d_gamma);
          while (d_loglkd < s * u * k) {
            u = t * u;
            gamma_shift_tmp = u * d_gamma;
            eta_tmp = eta + rep(gamma_shift_tmp, n_prov);
            // *** MODIFIED: use penalized log-likelihood for backtracking when firth ***
            d_loglkd = Loglkd(Y, eta_tmp);
            if (firth) {
              arma::vec p_tmp(n_obs);
              for (int ii = 0; ii < n_obs; ii++) { p_tmp(ii) = p_binomial(eta_tmp(ii)); }
              d_loglkd += FirthPenalty(p_tmp, n_prov, ind_start, n_gamma);
            }
            d_loglkd -= loglkd;
            // *** END MODIFIED ***
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
            // *** MODIFIED: Firth correction for score and information ***
            if (firth) {
              arma::vec pq_i = pq(span(ind_start(i), ind_start(i) + n_prov(i) - 1));
              arma::vec p_i = p(span(ind_start(i), ind_start(i) + n_prov(i) - 1));
              double W_i = info_gamma;
              arma::vec d_i = 1 - 2 * p_i;
              double A_i = dot(pq_i, d_i);
              score_gamma(i) += 0.5 * A_i / W_i;
              double B_i = dot(pq_i, d_i % d_i - 2 * pq_i);
              info_gamma = W_i - 0.5 * (B_i / W_i - A_i * A_i / (W_i * W_i));
            }
            // *** END MODIFIED ***
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
        r = (Y - p) / v;
      }

      arma::uvec update_order_unpenalized = randperm(K0);
      for (int j = 0; j < K0; j++) {
        shift = mean_crossprod(Z, r, update_order_unpenalized(j), n_obs);
        if (fabs(shift) > MaxChange_beta) { MaxChange_beta = fabs(shift); }
        beta(update_order_unpenalized(j)) = old_beta(update_order_unpenalized(j)) + shift;
        r -= Z.col(update_order_unpenalized(j)) * shift;
        eta += Z.col(update_order_unpenalized(j)) * shift;
        df++;
      }

      // update penalized groups
      for (int g = 0; g < n_group; g++) {
        if (active_group(g) == 1) {
          lambda_g = lambda * group_multiplier(g);
          gd_glm_grplasso(beta, Z, r, eta, old_beta, g, K1, n_obs, lambda_g, df, MaxChange_beta);
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
        for (int g = 0; g < n_group; g++) {
          if (active_group(g) == 1) {
            if (beta(K1(g)) == 0) { active_group(g) = 0; }
          }
        }
      }

      arma::vec Current_len_group(n_group, fill::zeros);
      for (int g = 0; g < n_group; g++) {
        if (active_group(g) == 0) {
          double lambda_g = lambda * group_multiplier(g);
          Current_len_group(g) = gd_glm_grplasso_BetaChange(Z, r, g, K1, n_obs, lambda_g);
        }
      }

      int if_add_new = 0;
      arma::uvec descend_len_index = sort_index(Current_len_group, "descend");
      arma::vec descend_len_value = sort(Current_len_group, "descend");

      for (int i = 0; i < activeGroupNum; i++) {
        if (descend_len_value(i) != 0) {
          if_add_new++;
          active_group(descend_len_index(i)) = 1;
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
// *** MODIFIED: added bool firth parameter ***
// [[Rcpp::export]]
List grp_lasso(arma::vec &Y, arma::mat &Z, arma::vec &n_prov, arma::vec &gamma, arma::vec &beta, int K0, arma::vec &K1,
               arma::vec &lambda_seq, bool lambda_early_stop, double stop_dev_ratio, arma::vec &group_multiplier,
               int max_total_iter, int max_each_iter, double tol, double nullDev,
               bool backtrack, double bound, int initial_active_group, double nvar_max, double group_max,
               bool trace_lambda, bool single_intercept, int threads, bool actSet, int actIter, int activeGroupNum, bool actSetRemove,
               bool firth) {  // *** MODIFIED: added firth ***

  int n_obs = Z.n_rows, n_beta = Z.n_cols, n_gamma = n_prov.n_elem, n_lambda = lambda_seq.n_elem, max_n_prov = max(n_prov);
  int n_group = K1.n_elem - 1;
  int tol_iter = 0;

  arma::mat beta_matrix(n_beta, n_lambda, fill::zeros), gamma_matrix(n_gamma, n_lambda, fill::zeros);
  arma::mat eta_matrix(n_obs, n_lambda, fill::zeros);
  arma::vec Dev_vec(n_lambda, fill::zeros);
  arma::vec iter_vec(n_lambda, fill::zeros);
  arma::vec df_vec(n_lambda, fill::zeros);
  arma::vec active_group(n_group, fill::zeros);

  if (actSet == true) {
    if (K0 == 0) { active_group(initial_active_group) = 1; }
  } else {
    active_group.ones();
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
      cout << "processing lambda: " << l + 1 << " (total: " << l + 1 << "/" << n_lambda << ")..." << endl;
    }
    double lambda = lambda_seq(l);
    // *** MODIFIED: pass firth to grp_lasso_fit ***
    auto fit = grp_lasso_fit(Y, Z, n_prov, gamma, beta, eta, K0, K1, lambda, tol_iter, max_total_iter, max_each_iter,
                             group_multiplier, max_n_prov, backtrack, bound, tol, ind_start, active_group,
                             n_obs, n_group, single_intercept, threads, actSet, actIter, activeGroupNum, actSetRemove,
                             firth);  // *** MODIFIED: added firth ***
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
      cout << "Warning: lambda " << l + 1 << "/" << n_lambda << " failed to converge within " << max_each_iter << " iterations!" << endl;
    }

    int ng = 0, nv = 0;
    for (int g = 0; g < n_group; g++) {
      if (beta(K1(g)) != 0) {
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
