// Weight solvers. These run on worker threads, so nothing here may call R.
// Interrupt checks go through Hook, which only the main thread sets.
// Ports fw.step, sc.weight.fw(.covariates), collapsed.form, sparsify_function,
// sum_normalize and synthdid_estimate's weights step by step, so results match
// the R code up to rounding.

#ifndef SYNTHDID_CORE_H
#define SYNTHDID_CORE_H

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace synthdid {

// Periodic callback (interrupt check). Main thread only.
typedef void (*Hook)();

// R's which.min: first minimum, NaN skipped.
inline arma::uword first_argmin(const arma::vec& v) {
  arma::uword best = 0;
  double best_val = std::numeric_limits<double>::infinity();
  bool found = false;
  for (arma::uword j = 0; j < v.n_elem; ++j) {
    if (std::isnan(v[j])) continue;
    if (!found || v[j] < best_val) {
      best = j;
      best_val = v[j];
      found = true;
    }
  }
  return best;
}

// One exact-line-search FW step on ||Ax - b||^2 + eta ||x||^2, x in simplex.
// Requires Ax == A * x on entry; updates both in place.
inline void fw_step(const arma::mat& A, const arma::vec& b, double eta,
                    arma::vec& x, arma::vec& Ax) {
  const arma::vec half_grad = A.t() * (Ax - b) + eta * x;
  const arma::uword i = first_argmin(half_grad);

  arma::vec d_x = -x;
  d_x[i] = 1.0 - x[i];
  if (arma::all(d_x == 0.0)) return;

  const arma::vec d_err = A.col(i) - Ax;
  const double num = -arma::dot(half_grad, d_x);
  const double den = arma::dot(d_err, d_err) + eta * arma::dot(d_x, d_x);

  // Matches R's clamp of +-Inf; 0/0 (NaN in R) means no move.
  double step;
  if (den > 0.0) {
    step = num / den;
  } else {
    step = (num > 0.0) ? 1.0 : 0.0;
  }
  step = std::min(1.0, std::max(0.0, step));
  if (step == 0.0) return;

  x += step * d_x;
  Ax += step * d_err;
}

inline void center_columns(arma::mat& M) { M.each_row() -= arma::mean(M, 0); }

// How the solver gets its gradient. Direct recomputes it every step, as the
// original R code did. Gram precomputes A'A once, which makes each step much
// cheaper on long solves. Auto starts direct and switches if the matrix fits.
enum SolverMode { SOLVER_AUTO = 0, SOLVER_DIRECT = 1, SOLVER_GRAM = 2 };

struct SolverConfig {
  SolverMode mode;
  double gram_max_bytes;  // largest allowed G, in bytes
  SolverConfig() : mode(SOLVER_AUTO), gram_max_bytes(64.0 * 1024.0 * 1024.0) {}
};

// Same step as fw_step, but using the precomputed A'A, so it never multiplies
// by the full matrix. Ax, Gx and the gradient g must be current on entry and
// are updated in place.
inline void fw_step_gram(const arma::mat& A, const arma::mat& G, double eta,
                         arma::vec& x, arma::vec& Ax, arma::vec& Gx,
                         arma::vec& g, arma::vec* Ax_orig = nullptr) {
  const arma::vec half_grad = g + eta * x;
  const arma::uword i = first_argmin(half_grad);

  arma::vec d_x = -x;
  d_x[i] = 1.0 - x[i];
  if (arma::all(d_x == 0.0)) return;

  const arma::vec d_err = A.col(i) - Ax;
  const double num = -arma::dot(half_grad, d_x);
  const double den = arma::dot(d_err, d_err) + eta * arma::dot(d_x, d_x);

  double step;
  if (den > 0.0) {
    step = num / den;
  } else {
    step = (num > 0.0) ? 1.0 : 0.0;
  }
  step = std::min(1.0, std::max(0.0, step));
  if (step == 0.0) return;

  // A' d_err = G e_i - G x.
  const arma::vec dG = G.col(i) - Gx;
  g += step * dG;
  Gx += step * dG;
  x += step * d_x;
  Ax += step * d_err;
  if (Ax_orig) *Ax_orig += step * d_err;
}

// C++ version of sc.weight.fw from synthdid. The last column of Y is the target.
// The stopping rule always uses the objective computed from Ax - b, not from
// A'A: with large outcome values, rounding at the switch to A'A can be bigger
// than the tolerance and stop the solver too early.
inline arma::vec sc_weight_fw(const arma::mat& Y, double zeta2, bool intercept,
                              arma::vec lambda, double min_decrease,
                              int max_iter, const SolverConfig& cfg,
                              std::vector<double>* vals = nullptr,
                              Hook hook = nullptr, int* gram_start = nullptr) {
  const arma::uword N0 = Y.n_rows;
  const arma::uword T0 = Y.n_cols - 1;

  arma::mat Yc = Y;
  if (intercept) center_columns(Yc);
  const arma::mat A = Yc.cols(0, T0 - 1);
  const arma::vec b = Yc.col(T0);
  const double eta = static_cast<double>(N0) * zeta2;
  const double tol = min_decrease * min_decrease;

  // Building G costs ~p direct iterations, so auto only switches after p.
  const arma::uword p = A.n_cols;
  const double gram_bytes = static_cast<double>(p) * static_cast<double>(p) * 8.0;
  const bool gram_allowed = cfg.mode != SOLVER_DIRECT && gram_bytes <= cfg.gram_max_bytes;
  const int switch_at = (cfg.mode == SOLVER_GRAM) ? 0 : static_cast<int>(p);
  // Periodic exact recompute of Gx and g bounds drift; amortized O(n).
  const int refresh_every = std::max<int>(static_cast<int>(p), 16);

  if (gram_start) *gram_start = -1;
  bool gram = false;
  arma::mat As, G;    // shifted regressors and their Gram matrix
  arma::vec bs, Asx;  // shifted target and As * lambda
  arma::vec Gx, g;
  int since_refresh = 0;

  arma::vec Ax = A * lambda;
  double prev = 0.0, curr = 0.0;
  int t = 0;
  while (t < max_iter && (t < 2 || prev - curr > tol)) {
    if (!gram && gram_allowed && t >= switch_at) {
      // Weights sum to one, so shifting A's columns and b by v leaves Ax - b and
      // the iterates unchanged. Row means strip shared levels that would
      // otherwise dominate G and cost precision.
      const arma::vec v = arma::mean(A, 1);
      As = A.each_col() - v;
      bs = b - v;
      Asx = Ax - v;
      G = As.t() * As;
      Gx = As.t() * Asx;
      g = As.t() * (Asx - bs);
      gram = true;
      since_refresh = 0;
      if (gram_start) *gram_start = t;
    }
    ++t;
    if (gram) {
      fw_step_gram(As, G, eta, lambda, Asx, Gx, g, &Ax);
      if (++since_refresh >= refresh_every) {
        Gx = As.t() * Asx;
        g = As.t() * (Asx - bs);
        since_refresh = 0;
      }
    } else {
      fw_step(A, b, eta, lambda, Ax);
    }
    const arma::vec err = Ax - b;
    prev = curr;
    curr = zeta2 * arma::dot(lambda, lambda) +
           arma::dot(err, err) / static_cast<double>(N0);
    if (vals) vals->push_back(curr);
    if (hook && t % 1000 == 0) hook();
  }
  return lambda;
}

// Mirrors update.weights in R.
struct WeightState {
  arma::vec lambda;
  arma::vec omega;
  arma::vec err_lambda;
  arma::vec err_omega;
  double val;
};

inline WeightState update_weights(const arma::mat& Y, const arma::vec& lambda,
                                  const arma::vec& omega, double zeta_lambda2,
                                  double zeta_omega2, bool lambda_intercept,
                                  bool omega_intercept, bool update_lambda,
                                  bool update_omega) {
  const arma::uword N0 = Y.n_rows - 1;
  const arma::uword T0 = Y.n_cols - 1;
  WeightState s;

  // lambda: control rows.
  arma::mat Y_lambda = Y.rows(0, N0 - 1);
  if (lambda_intercept) center_columns(Y_lambda);
  const arma::mat A_l = Y_lambda.cols(0, T0 - 1);
  const arma::vec b_l = Y_lambda.col(T0);
  s.lambda = lambda;
  arma::vec Ax_l = A_l * s.lambda;
  if (update_lambda) {
    fw_step(A_l, b_l, static_cast<double>(N0) * zeta_lambda2, s.lambda, Ax_l);
  }
  s.err_lambda = Ax_l - b_l;

  // omega: pre-treatment columns.
  arma::mat Y_omega = Y.cols(0, T0 - 1).t();
  if (omega_intercept) center_columns(Y_omega);
  const arma::mat A_o = Y_omega.cols(0, N0 - 1);
  const arma::vec b_o = Y_omega.col(N0);
  s.omega = omega;
  arma::vec Ax_o = A_o * s.omega;
  if (update_omega) {
    fw_step(A_o, b_o, static_cast<double>(T0) * zeta_omega2, s.omega, Ax_o);
  }
  s.err_omega = Ax_o - b_o;

  s.val = zeta_omega2 * arma::dot(s.omega, s.omega) +
          zeta_lambda2 * arma::dot(s.lambda, s.lambda) +
          arma::dot(s.err_omega, s.err_omega) / static_cast<double>(T0) +
          arma::dot(s.err_lambda, s.err_lambda) / static_cast<double>(N0);
  return s;
}

inline arma::mat contract(const arma::cube& X, const arma::vec& beta) {
  arma::mat out(X.n_rows, X.n_cols, arma::fill::zeros);
  for (arma::uword k = 0; k < X.n_slices; ++k) out += beta[k] * X.slice(k);
  return out;
}

struct CovariateResult {
  arma::vec lambda;
  arma::vec omega;
  arma::vec beta;
};

// Port of sc.weight.fw.covariates. Y and each X slice in collapsed form.
inline CovariateResult sc_weight_fw_covariates(
    const arma::mat& Y, const arma::cube& X, double zeta_lambda2,
    double zeta_omega2, bool lambda_intercept, bool omega_intercept,
    double min_decrease, int max_iter, const arma::vec& lambda,
    const arma::vec& omega, arma::vec beta, bool update_lambda,
    bool update_omega, std::vector<double>* vals = nullptr,
    Hook hook = nullptr) {
  const arma::uword N0 = Y.n_rows - 1;
  const arma::uword T0 = Y.n_cols - 1;
  const arma::uword K = X.n_slices;
  const double tol = min_decrease * min_decrease;

  arma::vec lambda_ext(T0 + 1), omega_ext(N0 + 1), grad(K);

  WeightState w = update_weights(Y - contract(X, beta), lambda, omega,
                                 zeta_lambda2, zeta_omega2, lambda_intercept,
                                 omega_intercept, update_lambda, update_omega);
  double prev = 0.0, curr = 0.0;
  int t = 0;
  while (t < max_iter && (t < 2 || std::abs(prev - curr) > tol)) {
    ++t;
    lambda_ext.head(T0) = w.lambda;
    lambda_ext[T0] = -1.0;
    omega_ext.head(N0) = w.omega;
    omega_ext[N0] = -1.0;
    for (arma::uword k = 0; k < K; ++k) {
      const arma::mat& Xk = X.slice(k);
      const double g_l =
          arma::dot(w.err_lambda, Xk.rows(0, N0 - 1) * lambda_ext);
      const double g_o =
          arma::dot(w.err_omega, Xk.cols(0, T0 - 1).t() * omega_ext);
      grad[k] = -(g_l / static_cast<double>(N0) + g_o / static_cast<double>(T0));
    }
    beta -= grad / static_cast<double>(t);
    w = update_weights(Y - contract(X, beta), w.lambda, w.omega, zeta_lambda2,
                       zeta_omega2, lambda_intercept, omega_intercept,
                       update_lambda, update_omega);
    prev = curr;
    curr = w.val;
    // Beta steps ignore covariate scale and can diverge; fail loudly.
    if (!std::isfinite(curr) || !beta.is_finite()) {
      throw std::runtime_error(
          "The covariate solver diverged: the coefficients on the covariates became non-finite. "
          "Standardizing the covariates (standardize = TRUE) or using covariate.method = 'projected' avoids this.");
    }
    if (vals) vals->push_back(curr);
    if (hook && t % 1000 == 0) hook();
  }
  CovariateResult out;
  out.lambda = w.lambda;
  out.omega = w.omega;
  out.beta = beta;
  return out;
}

// Projected adjustment (Kranz 2022): two-way FE regression of Y on X over the
// in_sample cells. FE removed by alternating demeaning (one pass if balanced).
// rcond flags covariates collinear with the FE; singular -> pinv.
inline arma::vec projected_beta_masked(const arma::mat& Y, const arma::cube& X,
                                       const arma::umat& in_sample,
                                       double* rcond = nullptr) {
  const arma::uword N = Y.n_rows, T = Y.n_cols, K = X.n_slices;
  const arma::uvec cells = arma::find(in_sample);

  arma::vec row_count = arma::conv_to<arma::vec>::from(arma::sum(in_sample, 1));
  arma::rowvec col_count = arma::conv_to<arma::rowvec>::from(arma::sum(in_sample, 0));
  const arma::mat mask = arma::conv_to<arma::mat>::from(in_sample);

  auto demean = [&](arma::mat M) {
    M %= mask;
    const double scale = std::max(1.0, arma::abs(M).max());
    for (int pass = 0; pass < 10000; ++pass) {
      arma::vec row_means = arma::sum(M, 1);
      for (arma::uword i = 0; i < N; ++i) row_means[i] = row_count[i] > 0 ? row_means[i] / row_count[i] : 0.0;
      M.each_col() -= row_means;
      M %= mask;
      arma::rowvec col_means = arma::sum(M, 0);
      for (arma::uword t = 0; t < T; ++t) col_means[t] = col_count[t] > 0 ? col_means[t] / col_count[t] : 0.0;
      M.each_row() -= col_means;
      M %= mask;
      const double change = std::max(arma::abs(row_means).max(), arma::abs(col_means).max());
      if (pass > 0 && change <= 1e-13 * scale) break;
    }
    return arma::vec(M.elem(cells));
  };

  const arma::vec y = demean(Y);
  arma::mat Z(cells.n_elem, K);
  arma::vec raw_ss(K);
  for (arma::uword k = 0; k < K; ++k) {
    const arma::vec raw = X.slice(k).elem(cells);
    raw_ss[k] = arma::dot(raw, raw);
    Z.col(k) = demean(X.slice(k));
  }
  // Correlation-scale so rcond and solve are unit-free. No variation left
  // after demeaning means collinear with the FE.
  const arma::mat ZZ = Z.t() * Z;
  const arma::vec Zy = Z.t() * y;
  arma::vec d = ZZ.diag();
  bool degenerate = false;
  for (arma::uword k = 0; k < K; ++k) {
    if (!(d[k] > 1e-20 * std::max(raw_ss[k], 1e-300))) { degenerate = true; d[k] = 1.0; }
  }
  const arma::vec root = arma::sqrt(d);
  const arma::mat scaled = ZZ / (root * root.t());
  const double rc = K == 0 ? 1.0 : (degenerate ? 0.0 : arma::rcond(scaled));
  if (rcond) *rcond = rc;
  const arma::vec rhs = Zy / root;
  const arma::vec scaled_beta = rc > 1e-12 ? arma::vec(arma::solve(scaled, rhs, arma::solve_opts::likely_sympd))
                                           : arma::vec(arma::pinv(scaled) * rhs);
  return scaled_beta / root;
}

// Block design: controls, plus treated pre-periods if untreated.
inline arma::vec projected_beta(const arma::mat& Y, const arma::cube& X,
                                arma::uword N0, arma::uword T0, bool untreated,
                                double* rcond = nullptr) {
  const arma::uword N = Y.n_rows, T = Y.n_cols;
  arma::umat in_sample(N, T, arma::fill::zeros);
  in_sample.rows(0, N0 - 1).fill(1);
  if (untreated && N > N0) in_sample.submat(N0, 0, N - 1, T0 - 1).fill(1);
  return projected_beta_masked(Y, X, in_sample, rcond);
}

// Port of collapsed.form: (N0 + 1) x (T0 + 1).
inline arma::mat collapsed_form(const arma::mat& Y, arma::uword N0,
                                arma::uword T0) {
  const arma::uword N = Y.n_rows, T = Y.n_cols;
  arma::mat out(N0 + 1, T0 + 1);
  out.submat(0, 0, N0 - 1, T0 - 1) = Y.submat(0, 0, N0 - 1, T0 - 1);
  out.submat(0, T0, N0 - 1, T0) = arma::mean(Y.submat(0, T0, N0 - 1, T - 1), 1);
  out.submat(N0, 0, N0, T0 - 1) = arma::mean(Y.submat(N0, 0, N - 1, T0 - 1), 0);
  out(N0, T0) = arma::accu(Y.submat(N0, T0, N - 1, T - 1)) /
                static_cast<double>((N - N0) * (T - T0));
  return out;
}

// Port of the default sparsify_function.
inline arma::vec sparsify(arma::vec v) {
  const double cut = v.max() / 4.0;
  v.elem(arma::find(v <= cut)).zeros();
  return v / arma::accu(v);
}

// Port of sum_normalize; uniform if sum is zero.
inline arma::vec sum_normalize(const arma::vec& x) {
  const double s = arma::accu(x);
  if (s != 0.0) return x / s;
  arma::vec u(x.n_elem);
  u.fill(1.0 / x.n_elem);
  return u;
}

// Mirrors attr(estimate, 'opts').
struct Options {
  double zeta_lambda2;
  double zeta_omega2;
  bool lambda_intercept;
  bool omega_intercept;
  double min_decrease;
  int max_iter;
  bool sparsify;              // false when sparsify = NULL
  int max_iter_pre_sparsify;  // first-round iterations when sparsify is on
  SolverConfig solver;        // only affects the no-covariate path
};

// synthdid_estimate's weights and point estimate, for resampling. Default
// sparsify between two rounds, or one round if off; custom sparsify stays in R.
inline double estimate(const arma::mat& Y, const arma::cube& X, arma::uword N0,
                       arma::uword T0, const Options& o, arma::vec lambda,
                       arma::vec omega, bool update_lambda, bool update_omega) {
  const arma::uword N1 = Y.n_rows - N0;
  const arma::uword T1 = Y.n_cols - T0;
  const arma::uword K = X.n_slices;
  arma::vec beta(K, arma::fill::zeros);

  if (K == 0) {
    const arma::mat Yc = collapsed_form(Y, N0, T0);
    const int first_round = o.sparsify ? o.max_iter_pre_sparsify : o.max_iter;
    if (update_lambda) {
      const arma::mat A = Yc.rows(0, N0 - 1);
      lambda = sc_weight_fw(A, o.zeta_lambda2, o.lambda_intercept, lambda,
                            o.min_decrease, first_round, o.solver);
      if (o.sparsify) {
        lambda = sc_weight_fw(A, o.zeta_lambda2, o.lambda_intercept,
                              sparsify(lambda), o.min_decrease, o.max_iter, o.solver);
      }
    }
    if (update_omega) {
      const arma::mat A = Yc.cols(0, T0 - 1).t();
      omega = sc_weight_fw(A, o.zeta_omega2, o.omega_intercept, omega,
                           o.min_decrease, first_round, o.solver);
      if (o.sparsify) {
        omega = sc_weight_fw(A, o.zeta_omega2, o.omega_intercept,
                             sparsify(omega), o.min_decrease, o.max_iter, o.solver);
      }
    }
  } else {
    const arma::mat Yc = collapsed_form(Y, N0, T0);
    arma::cube Xc(N0 + 1, T0 + 1, K);
    for (arma::uword k = 0; k < K; ++k) {
      Xc.slice(k) = collapsed_form(X.slice(k), N0, T0);
    }
    const CovariateResult r = sc_weight_fw_covariates(
        Yc, Xc, o.zeta_lambda2, o.zeta_omega2, o.lambda_intercept,
        o.omega_intercept, o.min_decrease, o.max_iter, lambda, omega, beta,
        update_lambda, update_omega);
    lambda = r.lambda;
    omega = r.omega;
    beta = r.beta;
  }

  arma::vec row_w(N0 + N1), col_w(T0 + T1);
  row_w.head(N0) = -omega;
  row_w.tail(N1).fill(1.0 / N1);
  col_w.head(T0) = -lambda;
  col_w.tail(T1).fill(1.0 / T1);
  const arma::mat Yb = (K == 0) ? Y : arma::mat(Y - contract(X, beta));
  return arma::dot(row_w, Yb * col_w);
}

}  // namespace synthdid

#endif  // SYNTHDID_CORE_H
