// Threaded bootstrap/placebo/jackknife. R draws all index sets up front, so
// results are identical across thread counts and match the R path per seed.
// Don't call R from worker threads. Interrupts, progress and copying results
// back to R all happen on the main thread.

#include "synthdid_core.h"

#include <atomic>
#include <functional>
#include <chrono>
#include <cstdio>
#include <condition_variable>
#include <exception>
#include <mutex>
#include <string>
#include <thread>

// [[Rcpp::depends(RcppArmadillo)]]

namespace {

void interrupt_probe(void*) { R_CheckUserInterrupt(); }

// Main thread only. R_ToplevelExec stops the longjmp unwinding C++ frames.
bool interrupt_pending() { return R_ToplevelExec(interrupt_probe, nullptr) == FALSE; }

std::string format_seconds(double s) {
  char buf[32];
  if (s < 60.0) {
    std::snprintf(buf, sizeof(buf), "%.1fs", s);
  } else {
    const long total = static_cast<long>(s + 0.5);
    std::snprintf(buf, sizeof(buf), "%ldm%02lds", total / 60, total % 60);
  }
  return buf;
}

// e.g.
//   Bootstrap  [===========>            ]  112/200  56%  1.2s elapsed  ~0.9s left  18 threads
std::string progress_line(const std::string& label, std::size_t done,
                          std::size_t total, double elapsed, int threads,
                          bool interrupted) {
  const int width = 24;
  const int filled = total == 0 ? width : static_cast<int>((done * width) / total);
  std::string bar(width, ' ');
  for (int i = 0; i < filled; ++i) bar[i] = '=';
  if (filled < width && done > 0) bar[filled] = '>';
  const int pct = total == 0 ? 100 : static_cast<int>((done * 100) / total);

  std::string status;
  if (interrupted) {
    status = "interrupted";
  } else if (done >= total) {
    status = "done";
  } else if (done == 0) {
    status = "";
  } else {
    status = "~" + format_seconds(elapsed * (total - done) / done) + " left";
  }

  char head[160];
  // Fixed-width count so redraws don't shift.
  const int count_width = static_cast<int>(std::to_string(total).size());
  std::snprintf(head, sizeof(head), "%s  [%s]  %*zu/%zu  %3d%%  %s elapsed", label.c_str(),
                bar.c_str(), count_width, done, total, pct, format_seconds(elapsed).c_str());
  std::string line = head;
  if (!status.empty()) line += "  " + status;
  line += "  " + std::to_string(threads) + (threads == 1 ? " thread" : " threads");
  return line;
}

// Progress on stderr, main thread only. Style 0: redraw in place, <= every
// 200 ms. Style 1 (non-interactive): a new line per 10%.
class Progress {
 public:
  Progress(bool enabled, int style, const std::string& label, std::size_t total)
      : enabled_(enabled), style_(style), label_(label), total_(total),
        start_(std::chrono::steady_clock::now()), last_draw_(start_) {}

  void set_threads(int threads) { threads_ = threads; }

  void update(std::size_t done, bool force = false) {
    if (!enabled_) return;
    const auto now = std::chrono::steady_clock::now();
    if (style_ == 0) {
      if (!force && (done == last_done_ ||
                     now - last_draw_ < std::chrono::milliseconds(200))) {
        return;
      }
      draw_in_place(done, false);
      last_draw_ = now;
    } else {
      const std::size_t step = total_ == 0 ? 10 : (done * 10) / total_;
      if (started_ && step < next_step_ && !force) return;
      print_line(done, false);
      started_ = true;
      next_step_ = step + 1;
    }
    last_done_ = done;
  }

  void finish(std::size_t done, bool interrupted) {
    if (!enabled_) return;
    if (style_ == 0) {
      draw_in_place(done, interrupted);
      REprintf("\n");
    } else if (interrupted || done != last_done_ || !started_) {
      print_line(done, interrupted);
    }
    R_FlushConsole();
  }

 private:
  double elapsed() const {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start_).count();
  }

  void draw_in_place(std::size_t done, bool interrupted) {
    const std::string line = progress_line(label_, done, total_, elapsed(), threads_, interrupted);
    // Pad to overwrite a longer previous line.
    const std::size_t pad = line.size() < last_len_ ? last_len_ - line.size() : 0;
    REprintf("\r%s%s", line.c_str(), std::string(pad, ' ').c_str());
    last_len_ = line.size();
    R_FlushConsole();
  }

  void print_line(std::size_t done, bool interrupted) {
    REprintf("%s\n", progress_line(label_, done, total_, elapsed(), threads_, interrupted).c_str());
    R_FlushConsole();
  }

  bool enabled_;
  int style_;
  std::string label_;
  std::size_t total_;
  int threads_ = 1;
  std::chrono::steady_clock::time_point start_, last_draw_;
  std::size_t last_done_ = static_cast<std::size_t>(-1);
  std::size_t last_len_ = 0;
  std::size_t next_step_ = 0;
  bool started_ = false;
};

// Runs body(r) for r in [0, R) on num_threads threads (0 = all cores).
// Exceptions go to errors[r]. Main thread polls interrupts: on one, cancel,
// drain in-flight work, then throw. Returns threads used.
int run_replications(std::size_t R, int num_threads, Progress& progress,
                     std::vector<std::string>& errors,
                     const std::function<void(std::size_t)>& body) {
  std::atomic<std::size_t> next(0);
  std::atomic<std::size_t> completed(0);
  std::atomic<bool> cancel(false);
  std::mutex mtx;
  std::condition_variable cv;
  std::size_t active = 0;

  auto worker = [&]() {
    for (;;) {
      if (cancel.load()) break;
      const std::size_t r = next.fetch_add(1);
      if (r >= R) break;
      try {
        body(r);
      } catch (const std::exception& e) {
        errors[r] = e.what();
      } catch (...) {
        errors[r] = "unknown error";
      }
      completed.fetch_add(1);
    }
    std::lock_guard<std::mutex> lock(mtx);
    --active;
    cv.notify_one();
  };

  std::size_t n_threads = num_threads > 0
                              ? static_cast<std::size_t>(num_threads)
                              : static_cast<std::size_t>(std::thread::hardware_concurrency());
  n_threads = std::max<std::size_t>(1, std::min(n_threads, R));

  std::vector<std::thread> pool;
  pool.reserve(n_threads);
  for (std::size_t i = 0; i < n_threads; ++i) {
    {
      std::lock_guard<std::mutex> lock(mtx);
      ++active;
    }
    try {
      pool.emplace_back(worker);
    } catch (...) {
      std::lock_guard<std::mutex> lock(mtx);
      --active;
      break;
    }
  }

  progress.set_threads(static_cast<int>(pool.empty() ? 1 : pool.size()));
  progress.update(0, true);

  bool interrupted = false;
  if (pool.empty()) {
    // No threads started; run inline.
    for (std::size_t r = 0; r < R; ++r) {
      if (interrupt_pending()) {
        interrupted = true;
        break;
      }
      try {
        body(r);
      } catch (const std::exception& e) {
        errors[r] = e.what();
      }
      completed.fetch_add(1);
      progress.update(completed.load());
    }
  } else {
    std::unique_lock<std::mutex> lock(mtx);
    while (active > 0) {
      cv.wait_for(lock, std::chrono::milliseconds(50));
      if (active > 0 && !interrupted) {
        lock.unlock();
        if (interrupt_pending()) {
          interrupted = true;
          cancel.store(true);
        } else {
          progress.update(completed.load());
        }
        lock.lock();
      }
    }
    lock.unlock();
    for (std::thread& th : pool) th.join();
  }
  progress.finish(completed.load(), interrupted);
  if (interrupted) throw Rcpp::internal::InterruptedException();
  return static_cast<int>(pool.empty() ? 1 : pool.size());
}

}  // namespace

// Lets the sequential R path print the same progress format.
// [[Rcpp::export(name = ".progress_line_cpp")]]
std::string progress_line_cpp(std::string label, double done, double total,
                              double elapsed, int threads, bool interrupted) {
  return progress_line(label, static_cast<std::size_t>(done), static_cast<std::size_t>(total),
                       elapsed, threads, interrupted);
}

// Re-estimates the effect on each resampled panel. Each column of ind holds the
// rows of Y used in one replication (1-based), with the controls first. See
// resample_estimates() in R/vcov.R, which calls it.
// [[Rcpp::export(name = ".synthdid_resample_cpp")]]
Rcpp::List synthdid_resample_cpp(
    const arma::mat& Y, const arma::cube& X, int T0, const Rcpp::IntegerMatrix& ind,
    const Rcpp::IntegerVector& n0, const arma::vec& omega,
    const arma::vec& lambda, bool update_lambda, bool update_omega,
    double zeta_lambda2, double zeta_omega2, bool lambda_intercept,
    bool omega_intercept, double min_decrease, int max_iter, int num_threads,
    int solver_mode, double gram_max_bytes, bool sparsify,
    int max_iter_pre_sparsify, bool verbose, int progress_style,
    std::string label, bool projected, bool projection_untreated) {
  const std::size_t R = ind.ncol();
  const arma::uword n_rows = ind.nrow();
  const arma::uword N = Y.n_rows;
  const arma::uword T = Y.n_cols;

  if (static_cast<std::size_t>(n0.size()) != R) Rcpp::stop("length(n0) must equal ncol(ind)");
  if (T0 < 1 || static_cast<arma::uword>(T0) >= T) Rcpp::stop("T0 must be in [1, ncol(Y) - 1]");
  if (X.n_rows != N || X.n_cols != T) Rcpp::stop("X must have the same first two dimensions as Y");
  if (max_iter_pre_sparsify < 0) Rcpp::stop("max_iter_pre_sparsify must be non-negative");
  if (lambda.n_elem != 0 && lambda.n_elem != static_cast<arma::uword>(T0)) {
    Rcpp::stop("lambda must have length T0 or 0");
  }

  // Copy out of R memory before threads start.
  std::vector<arma::uvec> rows(R);
  std::vector<arma::uword> n0v(R);
  for (std::size_t r = 0; r < R; ++r) {
    const int n0r = n0[r];
    n0v[r] = static_cast<arma::uword>(n0r);
    if (n0r < 1 || static_cast<arma::uword>(n0r) >= n_rows) {
      Rcpp::stop("each replication needs at least one control and one treated unit");
    }
    rows[r].set_size(n_rows);
    for (arma::uword i = 0; i < n_rows; ++i) {
      const int v = ind(i, r);
      if (v == NA_INTEGER || v < 1 || static_cast<arma::uword>(v) > N) {
        Rcpp::stop("row indices out of range");
      }
      if (omega.n_elem > 0 && i < static_cast<arma::uword>(n0r) &&
          static_cast<arma::uword>(v) > omega.n_elem) {
        Rcpp::stop("control indices must refer to rows with an omega weight");
      }
      rows[r][i] = static_cast<arma::uword>(v - 1);
    }
  }

  synthdid::Options opts;
  opts.zeta_lambda2 = zeta_lambda2;
  opts.zeta_omega2 = zeta_omega2;
  opts.lambda_intercept = lambda_intercept;
  opts.omega_intercept = omega_intercept;
  opts.min_decrease = min_decrease;
  opts.max_iter = max_iter;
  opts.sparsify = sparsify;
  opts.max_iter_pre_sparsify = max_iter_pre_sparsify;
  opts.solver.mode = static_cast<synthdid::SolverMode>(solver_mode);
  opts.solver.gram_max_bytes = gram_max_bytes;

  arma::vec lambda_init = lambda;
  if (lambda_init.n_elem == 0) {
    lambda_init.set_size(T0);
    lambda_init.fill(1.0 / T0);
  }

  std::vector<double> estimates(R, NA_REAL);
  std::vector<std::string> errors(R);

  auto run_one = [&](std::size_t r) {
    const arma::uvec& idx = rows[r];
    const arma::uword N0r = n0v[r];
    const arma::mat Yr = Y.rows(idx);
    arma::cube Xr(n_rows, T, X.n_slices);
    for (arma::uword k = 0; k < X.n_slices; ++k) Xr.slice(k) = X.slice(k).rows(idx);

    arma::vec omega_init;
    if (omega.n_elem > 0) {
      omega_init = synthdid::sum_normalize(omega.elem(idx.head(N0r)));
    } else {
      omega_init.set_size(N0r);
      omega_init.fill(1.0 / N0r);
    }
    if (projected && Xr.n_slices > 0) {
      const arma::vec beta = synthdid::projected_beta(Yr, Xr, N0r, static_cast<arma::uword>(T0), projection_untreated);
      const arma::mat Yadj = Yr - synthdid::contract(Xr, beta);
      const arma::cube none(n_rows, T, 0);
      return synthdid::estimate(Yadj, none, N0r, static_cast<arma::uword>(T0), opts,
                                lambda_init, omega_init, update_lambda, update_omega);
    }
    return synthdid::estimate(Yr, Xr, N0r, static_cast<arma::uword>(T0), opts,
                              lambda_init, omega_init, update_lambda,
                              update_omega);
  };

  Progress progress(verbose, progress_style, label, R);
  const int threads = run_replications(R, num_threads, progress, errors,
                                       [&](std::size_t r) { estimates[r] = run_one(r); });

  Rcpp::CharacterVector err(R);
  for (std::size_t r = 0; r < R; ++r) err[r] = errors[r];
  return Rcpp::List::create(
      Rcpp::Named("estimates") = Rcpp::NumericVector(estimates.begin(), estimates.end()),
      Rcpp::Named("errors") = err,
      Rcpp::Named("threads") = threads);
}

// Staggered bootstrap and placebo, Algorithms 2 and 4 of Clarke et al. (2024).
// For each replication, cohort says which cohort each sampled row belongs to
// (0 for never-treated). A cohort missing from a replication gets NA.
// Called from staggered_resample() in R/staggered.R.
// [[Rcpp::export(name = ".staggered_resample_cpp")]]
Rcpp::List staggered_resample_cpp(
    const arma::mat& Y, const arma::cube& X, const Rcpp::IntegerMatrix& ind,
    const Rcpp::IntegerMatrix& cohort, const Rcpp::IntegerVector& cohort_period,
    const arma::mat& omega, const Rcpp::List& lambda, bool update_lambda, bool update_omega,
    const arma::vec& zeta_lambda2, const arma::vec& zeta_omega2, bool lambda_intercept,
    bool omega_intercept, const arma::vec& min_decrease, int max_iter, int num_threads,
    int solver_mode, double gram_max_bytes, bool sparsify,
    int max_iter_pre_sparsify, bool verbose, int progress_style,
    std::string label, bool projected, bool projection_untreated) {
  const std::size_t R = ind.ncol();
  const arma::uword n_rows = ind.nrow();
  const arma::uword N = Y.n_rows;
  const arma::uword T = Y.n_cols;
  const arma::uword A = cohort_period.size();
  const arma::uword N0 = omega.n_rows;

  if (cohort.nrow() != ind.nrow() || cohort.ncol() != ind.ncol()) Rcpp::stop("cohort must have the dimensions of ind");
  if (X.n_rows != N || X.n_cols != T) Rcpp::stop("X must have the same first two dimensions as Y");
  if (omega.n_cols != A || static_cast<arma::uword>(lambda.size()) != A || zeta_lambda2.n_elem != A ||
      zeta_omega2.n_elem != A || min_decrease.n_elem != A) {
    Rcpp::stop("each cohort needs unit weights, time weights and settings");
  }
  std::vector<arma::uword> T0(A);
  std::vector<arma::vec> lambda_init(A);
  for (arma::uword a = 0; a < A; ++a) {
    const int p = cohort_period[a];
    if (p < 2 || static_cast<arma::uword>(p) > T) Rcpp::stop("cohort periods must be in [2, ncol(Y)]");
    T0[a] = static_cast<arma::uword>(p - 1);
    lambda_init[a] = Rcpp::as<arma::vec>(lambda[a]);
    if (lambda_init[a].n_elem != T0[a]) Rcpp::stop("each cohort's time weights must have one entry per pre-treatment period");
  }

  // Copy out of R memory before threads start.
  std::vector<arma::uvec> rows(R);
  std::vector<arma::uvec> codes(R);
  for (std::size_t r = 0; r < R; ++r) {
    rows[r].set_size(n_rows);
    codes[r].set_size(n_rows);
    bool any_control = false, any_treated = false;
    for (arma::uword i = 0; i < n_rows; ++i) {
      const int v = ind(i, r), c = cohort(i, r);
      if (v == NA_INTEGER || v < 1 || static_cast<arma::uword>(v) > N) Rcpp::stop("row indices out of range");
      if (c == NA_INTEGER || c < 0 || static_cast<arma::uword>(c) > A) Rcpp::stop("cohort codes out of range");
      if (c == 0 && static_cast<arma::uword>(v) > N0) Rcpp::stop("controls must refer to rows with unit weights");
      rows[r][i] = static_cast<arma::uword>(v - 1);
      codes[r][i] = static_cast<arma::uword>(c);
      any_control = any_control || c == 0;
      any_treated = any_treated || c > 0;
    }
    if (!any_control || !any_treated) Rcpp::stop("each replication needs at least one control and one treated unit");
  }

  std::vector<synthdid::Options> opts(A);
  for (arma::uword a = 0; a < A; ++a) {
    opts[a].zeta_lambda2 = zeta_lambda2[a];
    opts[a].zeta_omega2 = zeta_omega2[a];
    opts[a].lambda_intercept = lambda_intercept;
    opts[a].omega_intercept = omega_intercept;
    opts[a].min_decrease = min_decrease[a];
    opts[a].max_iter = max_iter;
    opts[a].sparsify = sparsify;
    opts[a].max_iter_pre_sparsify = max_iter_pre_sparsify;
    opts[a].solver.mode = static_cast<synthdid::SolverMode>(solver_mode);
    opts[a].solver.gram_max_bytes = gram_max_bytes;
  }

  std::vector<double> estimates(R, NA_REAL);
  std::vector<double> taus(R * A, NA_REAL);
  std::vector<std::string> errors(R);

  auto run_one = [&](std::size_t r) {
    const arma::uvec& idx = rows[r];
    const arma::uvec& code = codes[r];
    arma::mat Yr = Y.rows(idx);
    arma::cube Xr(n_rows, T, X.n_slices);
    for (arma::uword k = 0; k < X.n_slices; ++k) Xr.slice(k) = X.slice(k).rows(idx);
    const arma::uvec controls = arma::find(code == 0);

    if (projected && Xr.n_slices > 0) {
      arma::umat in_sample(n_rows, T, arma::fill::zeros);
      for (arma::uword i = 0; i < n_rows; ++i) {
        if (code[i] == 0) {
          in_sample.row(i).fill(1);
        } else if (projection_untreated) {
          in_sample.row(i).head(T0[code[i] - 1]).fill(1);
        }
      }
      const arma::vec beta = synthdid::projected_beta_masked(Yr, Xr, in_sample);
      Yr -= synthdid::contract(Xr, beta);
      Xr.set_size(n_rows, T, 0);
    }

    double weighted = 0.0, total = 0.0;
    for (arma::uword a = 0; a < A; ++a) {
      const arma::uvec members = arma::find(code == a + 1);
      if (members.n_elem == 0) continue;
      const arma::uvec sub = arma::join_cols(controls, members);
      const arma::mat Ya = Yr.rows(sub);
      arma::cube Xa(sub.n_elem, T, Xr.n_slices);
      for (arma::uword k = 0; k < Xr.n_slices; ++k) Xa.slice(k) = Xr.slice(k).rows(sub);
      const arma::vec omega_a = omega.col(a);
      const arma::vec omega_init = synthdid::sum_normalize(omega_a.elem(idx.elem(controls)));
      const double tau = synthdid::estimate(Ya, Xa, controls.n_elem, T0[a], opts[a], lambda_init[a],
                                            omega_init, update_lambda, update_omega);
      taus[r + R * a] = tau;
      const double share = static_cast<double>(members.n_elem) * static_cast<double>(T - T0[a]);
      weighted += share * tau;
      total += share;
    }
    estimates[r] = weighted / total;
  };

  Progress progress(verbose, progress_style, label, R);
  const int threads = run_replications(R, num_threads, progress, errors, run_one);

  Rcpp::CharacterVector err(R);
  for (std::size_t r = 0; r < R; ++r) err[r] = errors[r];
  Rcpp::NumericMatrix tau_matrix(R, A);
  std::copy(taus.begin(), taus.end(), tau_matrix.begin());
  return Rcpp::List::create(
      Rcpp::Named("estimates") = Rcpp::NumericVector(estimates.begin(), estimates.end()),
      Rcpp::Named("cohort.estimates") = tau_matrix,
      Rcpp::Named("errors") = err,
      Rcpp::Named("threads") = threads);
}
