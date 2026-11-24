## =====================================================================
## MCMC MA fit for ARFIMA(p,d,q), p,q in {0,1}, under MCAR
## missingness using adaptMCMC (Vihola 2012 robust adaptive Metropolis,
## RAM sampler)
##
## Samples directly in (d, [phi], [theta], sigma2)
## MODEL / SIGN CONVENTION (matches stats::arima AND nsarfima::arfima.sim):
##   (1 - phi B)(1 - B)^d y_t = (1 + theta B) eps_t,  eps_t ~ N(0, sigma2)
## NOTE: fracdiff.sim uses the OPPOSITE MA sign (pass ma = -theta there).
##
## LIKELIHOOD: APPROXIMATE Gaussian likelihood (Palma & Chan, 1997) from
##   the m-truncated MA state-space form -- the same likelihood as
##   ml_palma_chan_arfima_pdq.R. The Kalman filter skips the update step
##   at missing times
##
## SIGMA2 PRIOR: set with SIGMA2_PRIOR below ("unif" or "invgamma").
## =====================================================================

suppressPackageStartupMessages({
  library(adaptMCMC)
  library(coda)
  library(nsarfima)
  library(Rcpp)
})

## ---------------------------------------------------------------------
## CONFIGURATION
## ---------------------------------------------------------------------
ROOT_DIR  <- 'ROOT_DIR'
DONE_FILE <- file.path(ROOT_DIR, 'done_files.txt')   # reps already in combined.sqlite
SIGMA2_PRIOR <- "invgamma"   # "unif" or "invgamma"

## ---------------------------------------------------------------------
## PRIOR / PARAMETER SUPPORT
## ---------------------------------------------------------------------
D_LOWER_BOUND <- 0     # d bounded BELOW at 0 (long-memory / d > 0 regime only)
D_UPPER_BOUND <- 0.5   # stationarity boundary

PHI_LOWER_BOUND   <- -1  # causality
PHI_UPPER_BOUND   <-  1
THETA_LOWER_BOUND <- -1  # invertibility
THETA_UPPER_BOUND <-  1

SIGMA2_UPPER <- 10   # "unif":     sigma2 ~ U(0, SIGMA2_UPPER)
SIGMA2_SHAPE <- 1    # "invgamma": sigma2 ~ IG(SIGMA2_SHAPE, SIGMA2_RATE)
SIGMA2_RATE  <- 1

if (!SIGMA2_PRIOR %in% c("unif", "invgamma"))
  stop("SIGMA2_PRIOR must be \"unif\" or \"invgamma\"")

## =====================================================================
## Simulator
## =====================================================================
## nsarfima::arfima.sim follows the stats::arima convention for BOTH ar
## and ma 
sim_arfima <- function(n, d, phi = NULL, theta = NULL, sigma2 = 1) {
  nsarfima::arfima.sim(n, d = d,
                       ar   = if (is.null(phi))   numeric() else phi,
                       ma   = if (is.null(theta)) numeric() else theta,
                       sig2 = sqrt(sigma2))
}

## =====================================================================
## Autocovariances
## =====================================================================

## ---- Closed-form ACF of ARFIMA(0,d,0), gamma(0..n-1) (Hosking, 1981) ----
arfima00_acf <- function(d, n, sigma2 = 1) {
  g <- numeric(n)
  g[1] <- sigma2 * gamma(1 - 2 * d) / gamma(1 - d)^2
  if (n > 1) {
    for (k in 1:(n - 1)) g[k + 1] <- g[k] * (k - 1 + d) / (k - d)
  }
  g
}

## ---- Autocovariance gamma(0..n-1) of ARFIMA(1,d,1) ----

arfima_acf <- function(d, phi, theta, n, sigma2 = 1, tol = 1e-14) {
  K  <- if (phi == 0) 0L else as.integer(ceiling(log(tol) / log(abs(phi))))
  L  <- n + K + 1
  gu <- arfima00_acf(d, L + 1, sigma2)                          # gamma_u(0..L)
  gv <- (1 + theta^2) * gu[1:L] +
    theta * (c(gu[2], gu[1:(L - 1)]) + gu[2:(L + 1)])           # gamma_v(0..L-1)
  if (phi == 0) return(gv[1:n])
  ks <- -K:K
  w  <- phi^abs(ks)
  vapply(0:(n - 1), function(h) sum(w * gv[abs(h + ks) + 1]), numeric(1)) /
    (1 - phi^2)
}

## =====================================================================
## psi-weights and Omega_1
## =====================================================================

## Weights of (1-B)^{-d}: pi_0 = 1, pi_1 = d, ...
arfima00_psi <- function(d, n) {
  psi <- numeric(n)
  psi[1] <- 1
  if (n > 1) {
    for (j in 2:n) psi[j] <- psi[j - 1] * (d + (j - 2)) / (j - 1)
  }
  psi
}

## psi-weights of ARFIMA(1,d,1): psi_j = phi psi_{j-1} + pi_j + theta pi_{j-1}
arfima_psi <- function(d, phi, theta, n) {
  pi_d <- arfima00_psi(d, n)
  psi  <- numeric(n)
  psi[1] <- 1
  if (n > 1) {
    for (j in 2:n) psi[j] <- phi * psi[j - 1] + pi_d[j] + theta * pi_d[j - 1]
  }
  psi
}

## Initial state covariance Omega_1 (m x m)
##   Omega_1[i, i+h] = gamma(h) - sum_{k=1}^{i} psi_k psi_{k+h}
Omega1_from <- function(psi, gam) {
  m <- length(gam)
  Omega <- matrix(0, m, m)
  for (h in 0:(m - 1)) {
    i   <- 1:(m - h)
    val <- gam[h + 1] - cumsum(psi[i] * psi[i + h])
    Omega[cbind(i, i + h)] <- val
    Omega[cbind(i + h, i)] <- val
  }
  Omega
}

## =====================================================================
## Truncated Kalman filter (Palma & Chan, 1997)
## =====================================================================
kalman_trunc <- function(psi, Omega1, y, m) {
  n <- length(y)
  H <- psi[2:(m + 1)]                 # psi_1, ..., psi_m
  Q <- outer(H, H)

  Omega <- Omega1
  x <- numeric(m)
  e <- rep(NA_real_, n)
  D <- rep(NA_real_, n)

  for (t in 1:n) {
    D[t]  <- Omega[1, 1] + 1
    Theta <- c(Omega[2:m, 1], 0) + H

    Omega_shift <- matrix(0, m, m)
    Omega_shift[1:(m - 1), 1:(m - 1)] <- Omega[2:m, 2:m]

    if (is.na(y[t])) {
      x     <- c(x[2:m], 0)
      Omega <- Omega_shift + Q
    } else {
      e[t]  <- y[t] - x[1]
      x     <- c(x[2:m], 0) + Theta * e[t] / D[t]
      Omega <- Omega_shift + Q - outer(Theta, Theta) / D[t]
    }
  }
  obs <- !is.na(y)
  list(e = e[obs], D = D[obs])
}

## ---- C++ implementation 
Rcpp::sourceCpp(code = '
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List kalman_trunc_cpp(NumericVector psi, NumericMatrix Omega1,
                      NumericVector y, int m) {
  int n = y.size();
  std::vector<double> H(m), x(m, 0.0), Th(m);
  for (int i = 0; i < m; i++) H[i] = psi[i + 1];      // psi_1, ..., psi_m

  NumericMatrix Om = clone(Omega1);
  std::vector<double> e, D;
  e.reserve(n); D.reserve(n);

  for (int t = 0; t < n; t++) {
    double Dt = Om(0, 0) + 1.0;
    for (int i = 0; i < m - 1; i++) Th[i] = Om(i + 1, 0) + H[i];
    Th[m - 1] = H[m - 1];

    bool obs = !ISNAN(y[t]);
    double et = 0.0, g = 0.0;
    if (obs) { et = y[t] - x[0]; g = et / Dt; }

    // x <- c(x[2:m], 0) + Theta * e / D   (gain term only if observed)
    for (int i = 0; i < m - 1; i++) x[i] = x[i + 1] + Th[i] * g;
    x[m - 1] = Th[m - 1] * g;

    double c = obs ? 1.0 / Dt : 0.0;
    for (int j = 0; j < m; j++) {
      for (int i = 0; i <= j; i++) {
        double v = (j < m - 1) ? Om(i + 1, j + 1) : 0.0;
        v += H[i] * H[j] - c * Th[i] * Th[j];
        Om(i, j) = v;
        Om(j, i) = v;
      }
    }
    if (obs) { e.push_back(et); D.push_back(Dt); }
  }
  return List::create(_["e"] = wrap(e), _["D"] = wrap(D));
}
')

## e_t, D_t at (d, phi, theta); phi = 0 / theta = 0 when the term is absent
kalman_run <- function(d, phi, theta, y, m) {
  psi <- arfima_psi(d, phi, theta, m + 1)
  gam <- arfima_acf(d, phi, theta, m)
  kalman_trunc_cpp(psi, Omega1_from(psi[1:m], gam), y, m)
}

## =====================================================================
## Approximate (truncated-filter) likelihood of the observed values
## =====================================================================

## ---- Gaussian log-likelihood ----
## Returns -Inf if the filter fails or gives non-finite / non-positive D_t.
arfima_trunc_loglik <- function(d, phi, theta, sigma2, y, m) {
  r <- tryCatch(kalman_run(d, phi, theta, y, m), error = function(e) NULL)
  if (is.null(r) || !all(is.finite(c(r$e, r$D))) || any(r$D <= 0)) return(-Inf)
  n_obs <- length(r$e)
  -0.5 * (n_obs * log(2 * pi * sigma2) + sum(log(r$D)) + sum(r$e^2 / r$D) / sigma2)
}

## =====================================================================
## MCMC
## =====================================================================

## ---- Parameter vector layout: (d, [phi], [theta], sigma2) ----
par_names <- function(p, q) c("d", if (p == 1) "phi", if (q == 1) "theta", "sigma2")

unpack_par <- function(par, p, q) {
  list(d      = par[1],
       phi    = if (p == 1) par[2] else 0,
       theta  = if (q == 1) par[2 + p] else 0,
       sigma2 = par[2 + p + q])
}

## ---- Inverse-gamma log density ----
dinvgamma_log <- function(x, shape, rate) {
  shape * log(rate) - lgamma(shape) - (shape + 1) * log(x) - rate / x
}

## ---- sigma2 log prior, chosen by SIGMA2_PRIOR ----
log_prior_sigma2 <- function(sigma2) {
  if (SIGMA2_PRIOR == "unif") {
    dunif(sigma2, 0, SIGMA2_UPPER, log = TRUE)
  } else {
    dinvgamma_log(sigma2, SIGMA2_SHAPE, SIGMA2_RATE)
  }
}

log_prior <- function(pp, p, q) {
  lp <- dunif(pp$d, D_LOWER_BOUND, D_UPPER_BOUND, log = TRUE) +
    log_prior_sigma2(pp$sigma2)
  if (p == 1) lp <- lp + dunif(pp$phi,   PHI_LOWER_BOUND,   PHI_UPPER_BOUND,   log = TRUE)
  if (q == 1) lp <- lp + dunif(pp$theta, THETA_LOWER_BOUND, THETA_UPPER_BOUND, log = TRUE)
  lp
}

## ---- Log-posterior target ----
make_log_posterior <- function(y, p, q, m) {
  k <- 2L + p + q
  function(par) {
    if (length(par) != k || !all(is.finite(par))) return(-Inf)
    pp <- unpack_par(par, p, q)
    if (pp$sigma2 <= 0) return(-Inf)
    if (pp$d <= D_LOWER_BOUND + 0.01 || pp$d >= D_UPPER_BOUND - 0.01) return(-Inf)
    if (p == 1 && (pp$phi   <= PHI_LOWER_BOUND   + 0.01 ||
                   pp$phi   >= PHI_UPPER_BOUND   - 0.01)) return(-Inf)
    if (q == 1 && (pp$theta <= THETA_LOWER_BOUND + 0.01 ||
                   pp$theta >= THETA_UPPER_BOUND - 0.01)) return(-Inf)

    ll <- arfima_trunc_loglik(pp$d, pp$phi, pp$theta, pp$sigma2, y, m)
    if (!is.finite(ll)) return(-Inf)

    ll + log_prior(pp, p, q)
  }
}

## ---- Starting values ----
get_mle <- function(target, init, max_iter = 20) {
  neg_target <- function(x) {
    v <- target(x)
    if (is.finite(v)) -v else 1e10
  }
  for (i in 1:max_iter) {
    opt  <- optim(init, neg_target)
    init <- opt$par
    if (opt$convergence == 0) break
  }
  init
}

## ---- Single-chain run for one replicate -----------------------------
run_mcmc_one_rep <- function(y, p = 0L, q = 0L, m,
                             n_burnin = 1000, n_sample = 20000,
                             seed = NULL, show_progress = TRUE) {
  if (!is.null(seed)) set.seed(seed)   # RNG discipline: seed per replicate

  target <- make_log_posterior(y, p, q, m)

  ## start the chain at the posterior mode (MLE / MAP, see get_mle)
  start <- c(0.25, if (p == 1) 0, if (q == 1) 0, 1)
  init  <- get_mle(target, start)

  ## initial proposal VARIANCES (diagonal); adaptation corrects these.
  init_scale <- c(0.005, if (p == 1) 0.01, if (q == 1) 0.01, 0.05)

  n_total <- n_burnin + n_sample

  fit <- adaptMCMC::MCMC(
    p               = target,
    n               = n_total,
    init            = init,
    scale           = init_scale,
    adapt           = n_burnin,   # freeze covariance after burn-in
    acc.rate        = 0.234,
    showProgressBar = show_progress
  )

  post_idx      <- (n_burnin + 1):n_total
  theta_samples <- fit$samples[post_idx, , drop = FALSE]
  colnames(theta_samples) <- par_names(p, q)

  list(
    theta_samples   = theta_samples,
    acceptance.rate = fit$acceptance.rate,
    cov.jump        = fit$cov.jump,
    ess             = coda::effectiveSize(coda::mcmc(theta_samples)),
    init            = setNames(init, par_names(p, q)),
    init_scale      = init_scale,
    raw_fit         = fit
  )
}

## ---- Summarise posterior draws: median, 90% interval, interval width ----
summarise_mcmc_rep <- function(theta_samples, method_label = "MCMC") {
  out <- data.frame(method = method_label)
  for (nm in colnames(theta_samples)) {
    s  <- theta_samples[, nm]
    qq <- unname(quantile(s, c(0.05, 0.95)))
    out[[paste0(nm, "_median")]] <- median(s)
    out[[paste0(nm, "_lower")]]  <- qq[1]
    out[[paste0(nm, "_upper")]]  <- qq[2]
    out[[paste0(nm, "_ciw")]]    <- qq[2] - qq[1]
  }
  out
}

## =====================================================================
## CLI entry point:
##   rep_start rep_end d_true n miss_pct m [phi_true] [theta_true]
## Omitted or "NA" phi_true/theta_true => that term is absent, so:
##   ...                -> ARFIMA(0,d,0)
##   ... 0.5            -> ARFIMA(1,d,0)
##   ... NA 0.3         -> ARFIMA(0,d,1)
##   ... 0.5 0.3        -> ARFIMA(1,d,1)
## =====================================================================
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  usage <- paste0(
    "Usage: Rscript mcmc_benchmark_arfima_pdq_approx.R rep_start rep_end d_true n miss_pct m [phi_true] [theta_true]\n",
    "       (omit or pass NA for phi_true / theta_true to drop that term)")
  if (length(args) < 6) stop(usage)

  rep_start <- as.integer(args[1])
  rep_end   <- as.integer(args[2])
  d_true    <- as.numeric(args[3])
  n         <- as.integer(args[4])
  miss_pct  <- as.numeric(args[5])
  m         <- as.integer(args[6])
  phi_true   <- if (length(args) >= 7) suppressWarnings(as.numeric(args[7])) else NA_real_
  theta_true <- if (length(args) >= 8) suppressWarnings(as.numeric(args[8])) else NA_real_
  p <- as.integer(!is.na(phi_true))
  q <- as.integer(!is.na(theta_true))
  if (m < 2) stop("m must be >= 2 (truncated filter only)")
  if (p == 1 && abs(phi_true)   >= 1) stop("phi_true must satisfy |phi| < 1 (causality)")
  if (q == 1 && abs(theta_true) >= 1) stop("theta_true must satisfy |theta| < 1 (invertibility)")

  s2_true <- 1
  dir.create(ROOT_DIR, recursive = TRUE, showWarnings = FALSE)

  ## reps already moved into the combined DB (their .rds may no longer be on disk)
  done_files <- if (file.exists(DONE_FILE)) readLines(DONE_FILE) else character(0)

  arma_tag  <- paste0(if (p == 1) sprintf("_phi%.2f", phi_true) else "",
                      if (q == 1) sprintf("_theta%.2f", theta_true) else "")
  prior_tag <- if (SIGMA2_PRIOR == "invgamma") "_s2ig" else ""

  for (rep_id in rep_start:rep_end) {
    out_file <- file.path(ROOT_DIR, sprintf("mcmc_rep%04d_d%.2f%s_pctmis%.2f%s_APPROXm%d.rds",
                                            rep_id, d_true, arma_tag, miss_pct, prior_tag, m))
    if (file.exists(out_file) || basename(out_file) %in% done_files) next   # resume-safe skip

    seed <- rep_id
    set.seed(seed)
    t.start <- Sys.time()

    ## simulate, then apply MCAR mask
    y_full  <- sim_arfima(n, d = d_true,
                          phi   = if (p == 1) phi_true   else NULL,
                          theta = if (q == 1) theta_true else NULL,
                          sigma2 = s2_true)
    present <- rep(TRUE, n)
    present[sample.int(n, round(miss_pct * n))] <- FALSE
    y_obs   <- y_full
    y_obs[!present] <- NA

    fit <- run_mcmc_one_rep(y_obs, p = p, q = q, m = m, seed = rep_id)
    fit$meta <- list(rep_id = rep_id, seed = seed, p = p, q = q,
                     d_true = d_true,
                     phi_true   = if (p == 1) phi_true   else NA_real_,
                     theta_true = if (q == 1) theta_true else NA_real_,
                     sigma2_true = s2_true, n = n, miss_pct = miss_pct,
                     likelihood = "approximate", m = m,
                     sigma2_prior = SIGMA2_PRIOR,
                     sigma2_prior_pars = if (SIGMA2_PRIOR == "unif")
                       c(lower = 0, upper = SIGMA2_UPPER) else
                       c(shape = SIGMA2_SHAPE, rate = SIGMA2_RATE),
                     runtime_secs = as.numeric(Sys.time() - t.start, units = "secs"))

    ## atomic write: a killed job can't leave a truncated .rds that the
    ## resume-skip would later trust
    tmp_file <- paste0(out_file, ".tmp", Sys.getpid())
    saveRDS(fit, tmp_file)
    file.rename(tmp_file, out_file)

    cat(fit$meta$runtime_secs, "seconds\n")
  }
}

## Run only via Rscript, so the file can be source()'d for testing.
if (sys.nframe() == 0) main()