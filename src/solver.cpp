// R entry points for the single-estimate weight solvers.
// The numerical work lives in synthdid_core.h.

#include "synthdid_core.h"

// [[Rcpp::depends(RcppArmadillo)]]

namespace {

void check_interrupt() { Rcpp::checkUserInterrupt(); }

Rcpp::NumericVector as_r_vector(const arma::vec& v) {
  return Rcpp::NumericVector(v.begin(), v.end());
}

// Pads the objective trace with NA to max_iter entries, as the R solvers did.
Rcpp::NumericVector padded_vals(const std::vector<double>& vals, int max_iter) {
  Rcpp::NumericVector out(max_iter, NA_REAL);
  const std::size_t n = std::min(vals.size(), static_cast<std::size_t>(out.size()));
  for (std::size_t i = 0; i < n; ++i) out[i] = vals[i];
  return out;
}

}  // namespace

// Frank-Wolfe solver for synthetic control weights with exact line search.
// zeta2 is Re(zeta^2) as computed on the R side. solver_mode is 0 for auto,
// 1 for direct and 2 for Gram, see SolverMode in synthdid_core.h.
// [[Rcpp::export(name = ".sc_weight_fw_cpp")]]
Rcpp::List sc_weight_fw_cpp(const arma::mat& Y, double zeta2, bool intercept,
                            const arma::vec& lambda, double min_decrease,
                            int max_iter, int solver_mode,
                            double gram_max_bytes) {
  synthdid::SolverConfig cfg;
  cfg.mode = static_cast<synthdid::SolverMode>(solver_mode);
  cfg.gram_max_bytes = gram_max_bytes;
  std::vector<double> vals;
  int gram_start = -1;
  const arma::vec out = synthdid::sc_weight_fw(Y, zeta2, intercept, lambda,
                                               min_decrease, max_iter, cfg,
                                               &vals, &check_interrupt,
                                               &gram_start);
  return Rcpp::List::create(
      Rcpp::Named("lambda") = as_r_vector(out),
      Rcpp::Named("vals") = padded_vals(vals, max_iter),
      Rcpp::Named("gram.start") = gram_start < 0 ? NA_INTEGER : gram_start);
}

// Frank-Wolfe plus gradient solver for lambda, omega and beta with covariates.
// zeta_lambda2 and zeta_omega2 are Re(zeta^2) as computed on the R side.
// [[Rcpp::export(name = ".sc_weight_fw_covariates_cpp")]]
Rcpp::List sc_weight_fw_covariates_cpp(
    const arma::mat& Y, const arma::cube& X, double zeta_lambda2,
    double zeta_omega2, bool lambda_intercept, bool omega_intercept,
    double min_decrease, int max_iter, const arma::vec& lambda,
    const arma::vec& omega, const arma::vec& beta, bool update_lambda,
    bool update_omega) {
  std::vector<double> vals;
  const synthdid::CovariateResult r = synthdid::sc_weight_fw_covariates(
      Y, X, zeta_lambda2, zeta_omega2, lambda_intercept, omega_intercept,
      min_decrease, max_iter, lambda, omega, beta, update_lambda, update_omega,
      &vals, &check_interrupt);
  return Rcpp::List::create(Rcpp::Named("lambda") = as_r_vector(r.lambda),
                            Rcpp::Named("omega") = as_r_vector(r.omega),
                            Rcpp::Named("beta") = as_r_vector(r.beta),
                            Rcpp::Named("vals") = padded_vals(vals, max_iter));
}

// Coefficients for the projected covariate adjustment, with the reciprocal
// condition number of the demeaned covariates' cross-product. See
// projected_beta in synthdid_core.h.
// [[Rcpp::export(name = ".projected_beta_cpp")]]
Rcpp::List projected_beta_cpp(const arma::mat& Y, const arma::cube& X, int N0, int T0, bool untreated) {
  double rcond = 0.0;
  const arma::vec beta = synthdid::projected_beta(Y, X, static_cast<arma::uword>(N0), static_cast<arma::uword>(T0),
                                                  untreated, &rcond);
  return Rcpp::List::create(Rcpp::Named("beta") = as_r_vector(beta), Rcpp::Named("rcond") = rcond);
}

// The projected regression over an arbitrary sample of cells, marked by
// nonzero entries of in_sample, for designs whose sample is not a block, such
// as staggered adoption.
// [[Rcpp::export(name = ".projected_beta_masked_cpp")]]
Rcpp::List projected_beta_masked_cpp(const arma::mat& Y, const arma::cube& X, const arma::imat& in_sample) {
  double rcond = 0.0;
  const arma::umat mask = arma::conv_to<arma::umat>::from(in_sample != 0);
  const arma::vec beta = synthdid::projected_beta_masked(Y, X, mask, &rcond);
  return Rcpp::List::create(Rcpp::Named("beta") = as_r_vector(beta), Rcpp::Named("rcond") = rcond);
}
