#include "utils.h"

// ============================================================
//  Stratified Cox Model with Group LASSO
//  — Cox proportional hazards, stratified by provider
//  — Group-level LASSO penalty
// ============================================================


// ---------- Group descent update ----------
void gd_stratCox(arma::vec &beta, arma::mat &Z, arma::vec &r, arma::vec &eta, arma::vec &old_beta,
                 int g, arma::vec &K1, int n_obs, double lambda, double &df, double &MaxChange_beta) {
  int K = K1(g + 1) - K1(g);
  arma::vec beta_initial(K);
  for (int j = K1(g); j < K1(g + 1); j++) {
    beta_initial(j - K1(g)) = mean_crossprod(Z, r, j, n_obs) + old_beta(j);
  }

  double beta_initial_norm = arma::norm(beta_initial, 2);
  double len = Soft_thres(beta_initial_norm, lambda / 1); // v = 1

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
double gd_stratCox_BetaChange(arma::mat &Z, arma::vec &r, int g, arma::vec &K1, int n_obs, double lambda) {
  int K = K1(g + 1) - K1(g);
  arma::vec beta_initial(K);
  for (int j = K1(g); j < K1(g + 1); j++) {
    beta_initial(j - K1(g)) = mean_crossprod(Z, r, j, n_obs);
  }
  double beta_initial_norm = arma::norm(beta_initial, 2);
  double len = Soft_thres(beta_initial_norm, lambda / 1);

  if (len != 0) {
    return(len);
  } else {
    return(0);
  }
}


// ---------- Single-lambda fit ----------
tuple<arma::vec, arma::vec, double, double, int> StratCox_lasso_fit(
    arma::vec &delta_obs, arma::mat &Z, arma::vec &weight, arma::vec &n_each_prov,
    arma::vec beta, arma::vec eta, int K0, arma::vec &K1,
    double lambda, int &tol_iter, int max_total_iter, int max_each_iter,
    arma::vec &group_multiplier, int count_stratum, double tol, arma::vec &ind_start,
    arma::vec &active_group, int n_obs, int n_group,
    bool actSet, int actIter, int activeGroupNum, bool actSetRemove) {

  arma::vec old_beta = beta, r(n_obs), r_shift;
  arma::vec haz(n_obs), rsk(n_obs), h(n_obs);
  double loss = 0, df, MaxChange_beta, shift, lambda_g;
  double s, v = 1.0;
  int iter = 0;

  while (tol_iter < max_total_iter) {
    int inner_loop_iter = -1;
    inner_loop_iter = inner_loop_iter + 1;
    R_CheckUserInterrupt();
    while (tol_iter < max_total_iter && iter < max_each_iter) {
      R_CheckUserInterrupt();
      df = 0;
      tol_iter++;
      iter++;
      inner_loop_iter++;
      MaxChange_beta = 0;

      // calculate haz, weighted haz, and weighted risk set
      haz = arma::exp(eta);
      arma::vec whaz = weight % haz;  // w_i * exp(eta_i)
      for (int j = 0; j < count_stratum; j++) {
        rsk(ind_start(j) + n_each_prov(j) - 1) = whaz(ind_start(j) + n_each_prov(j) - 1);
        for (int i = ind_start(j) + n_each_prov(j) - 2; i >= ind_start(j); i--) {
          rsk(i) = rsk(i + 1) + whaz(i);
        }
      }

      // calculate l'(eta): h(i) = sum_{k<=i} w_k * delta_k / rsk_k
      for (int j = 0; j < count_stratum; j++) {
        h(ind_start(j)) = weight(ind_start(j)) * delta_obs(ind_start(j)) / rsk(ind_start(j));
        for (int i = ind_start(j) + 1; i < ind_start(j) + n_each_prov(j); i++) {
          h(i) = h(i - 1) + weight(i) * delta_obs(i) / rsk(i);
        }
      }

      double a;
      for (int i = 0; i < n_obs; i++) {
        a = haz(i) * h(i);
        if (a == 0) {
          r(i) = 0;
        } else {
          s = weight(i) * (delta_obs(i) - a);
          r(i) = s / v;
        }
      }

      // weighted partial log-likelihood: sum_i w_i * delta_i * (eta_i - log(rsk_i))
      loss = 0;
      for (int i = 0; i < n_obs; i++) {
        loss += weight(i) * delta_obs(i) * (eta(i) - log(rsk(i)));
      }

      // update unpenalized groups
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
          gd_stratCox(beta, Z, r, eta, old_beta, g, K1, n_obs, lambda_g, df, MaxChange_beta);
        }
      }

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
          Current_len_group(g) = gd_stratCox_BetaChange(Z, r, g, K1, n_obs, lambda_g);
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

  return make_tuple(beta, eta, loss, df, iter);
}


// ---------- Full lambda-path wrapper (Rcpp exported) ----------
// [[Rcpp::export]]
List StratCox_lasso(arma::vec &delta_obs, arma::mat &Z, arma::vec &weight, arma::vec &n_each_prov, arma::vec &beta,
                    int K0, arma::vec &K1, arma::vec &lambda_seq, bool lambda_early_stop, double stop_loss_ratio,
                    arma::vec &group_multiplier, int max_total_iter, int max_each_iter, double tol,
                    int initial_active_group, double nvar_max, double group_max, bool trace_lambda,
                    bool actSet, int actIter, int activeGroupNum, bool actSetRemove) {

  int n_obs = Z.n_rows, n_beta = Z.n_cols, n_lambda = lambda_seq.n_elem, n_group = K1.n_elem - 1;
  int tol_iter = 0;
  int count_stratum = n_each_prov.n_elem;

  arma::mat beta_matrix(n_beta, n_lambda, fill::zeros);
  arma::mat eta_matrix(n_obs, n_lambda, fill::zeros);
  arma::vec iter_vec(n_lambda, fill::zeros);
  arma::vec df_vec(n_lambda, fill::zeros);
  arma::vec loss_vec(n_lambda, fill::zeros);
  arma::vec active_group(n_group, fill::zeros);

  if (actSet == true) {
    if (K0 == 0) { active_group(initial_active_group) = 1; }
  } else {
    active_group.ones();
  }

  arma::vec ind_start(count_stratum);
  ind_start(0) = 0;
  for (int i = 1; i < count_stratum; i++) {
    ind_start(i) = ind_start(i - 1) + n_each_prov(i - 1);
  }

  arma::vec eta = Z * beta;

  for (int l = 0; l < n_lambda; l++) {
    R_CheckUserInterrupt();
    if (trace_lambda == true) {
      Rcpp::Rcout << "processing lambda: " << l + 1 << " (total: " << l + 1 << "/" << n_lambda << ")..." << endl;
    }
    double lambda = lambda_seq(l);

    auto fit = StratCox_lasso_fit(delta_obs, Z, weight, n_each_prov, beta, eta, K0, K1, lambda, tol_iter,
                                  max_total_iter, max_each_iter, group_multiplier, count_stratum, tol, ind_start,
                                  active_group, n_obs, n_group, actSet, actIter, activeGroupNum, actSetRemove);
    double loss_l, df_l;
    int iter_l;
    tie(beta, eta, loss_l, df_l, iter_l) = fit;
    beta_matrix.col(l) = beta;
    eta_matrix.col(l) = eta;
    loss_vec(l) = loss_l;
    df_vec(l) = df_l;
    iter_vec(l) = iter_l;

    if (iter_l == max_each_iter) {
      Rcpp::Rcout << "Warning: lambda " << l + 1 << "/" << n_lambda << " failed to converge within " << max_each_iter << " iterations!" << endl;
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
        Rcpp::Rcout << "Algorithm has reached the maximum number of total iterations, stops..." << endl;
      } else if (ng > group_max) {
        Rcpp::Rcout << "Algorithm has selected the maximum number of groups, stops..." << endl;
      } else {
        Rcpp::Rcout << "Algorithm has selected the maximum number of variables, stops..." << endl;
      }
      for (int ll = (l + 1); ll < n_lambda; ll++) { iter_vec(ll) = NA_REAL; }
      break;
    }

    if (lambda_early_stop == true) {
      if (l != 0) {
        double null_lkd = loss_vec(0);
        double loss_ratio = fabs((loss_vec(l) - loss_vec(l - 1)) / (loss_vec(l) - null_lkd));
        if (loss_ratio < stop_loss_ratio) {
          for (int ll = (l + 1); ll < n_lambda; ll++) { iter_vec(ll) = NA_REAL; }
          break;
        }
      }
    }
  }

  List result = List::create(_["beta"] = beta_matrix, _["loss"] = loss_vec,
                             _["Eta"] = eta_matrix, _["Df"] = df_vec, _["iter"] = iter_vec);
  return result;
}
