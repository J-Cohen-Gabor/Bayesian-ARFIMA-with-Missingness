## =====================================================================
## MCMC exact fit for ARFIMA(p,d,q), p,q in {0,1}, under MCAR
## missingness using adaptMCMC (Vihola 2012 robust adaptive Metropolis,
## RAM sampler)
##
## Samples directly in (d, [phi], [theta], sigma2)
## MODEL / SIGN CONVENTION (matches stats::arima AND nsarfima::arfima.sim):
##   (1 - phi B)(1 - B)^d y_t = (1 + theta B) eps_t,  eps_t ~ N(0, sigma2)
## NOTE: fracdiff.sim uses the OPPOSITE MA sign (pass ma = -theta there).
##
## LIKELIHOOD: EXACT Gaussian likelihood of the observed values.
##   Under MCAR the observed values y_o = (y_{t_1}, ..., y_{t_no}) are
##   jointly Gaussian:
##     y_o ~ N(0, sigma2 * Gamma_oo),  Gamma_oo[i, j] = gamma(|t_i - t_j|)
##   where gamma() is the ARFIMA autocovariance at sigma2 = 1 (arfima_acf).
##   
## SIGMA2 PRIOR: set with SIGMA2_PRIOR below ("unif" or "invgamma").
## =====================================================================

suppressPackageStartupMessages({
  library(adaptMCMC)
  library(coda)
  library(nsarfima)
})

## ---------------------------------------------------------------------
## CONFIGURATION
## ---------------------------------------------------------------------
ROOT_DIR     <- 'ROOT_DIR'
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
## Exact likelihood of the observed values
## =====================================================================

prep_obs <- function(y) {
  t_obs <- which(!is.na(y))
  list(y_o     = y[t_obs],
       n_o     = length(t_obs),
       max_lag = max(t_obs) - min(t_obs),
       lag_idx = abs(outer(t_obs, t_obs, "-")) + 1)
}

## ---- Exact Gaussian log-likelihood of y_o ----
## Returns -Inf if Gamma_oo is not numerically positive definite.
arfima_exact_loglik <- function(d, phi, theta, sigma2, obs) {
  gam <- arfima_acf(d, phi, theta, obs$max_lag + 1)        # gamma(0..max_lag)
  G   <- matrix(gam[obs$lag_idx], obs$n_o, obs$n_o)        # Gamma_oo

  U <- tryCatch(chol(G), error = function(e) NULL)         # G = t(U) %*% U
  if (is.null(U)) return(-Inf)

  z <- backsolve(U, obs$y_o, transpose = TRUE)             # solves t(U) %*% z = y_o

  -obs$n_o / 2 * log(2 * pi * sigma2) - sum(log(diag(U))) - sum(z^2) / (2 * sigma2)
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
make_log_posterior <- function(y, p, q) {
  obs <- prep_obs(y)          # done once, not at every evaluation
  k   <- 2L + p + q
  function(par) {
    if (length(par) != k || !all(is.finite(par))) return(-Inf)
    pp <- unpack_par(par, p, q)
    if (pp$sigma2 <= 0) return(-Inf)
    if (pp$d <= D_LOWER_BOUND + 0.01 || pp$d >= D_UPPER_BOUND - 0.01) return(-Inf)
    if (p == 1 && (pp$phi   <= PHI_LOWER_BOUND   + 0.01 ||
                   pp$phi   >= PHI_UPPER_BOUND   - 0.01)) return(-Inf)
    if (q == 1 && (pp$theta <= THETA_LOWER_BOUND + 0.01 ||
                   pp$theta >= THETA_UPPER_BOUND - 0.01)) return(-Inf)

    ll <- arfima_exact_loglik(pp$d, pp$phi, pp$theta, pp$sigma2, obs)
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
run_mcmc_one_rep <- function(y, p = 0L, q = 0L,
                             n_burnin = 1000, n_sample = 20000,
                             seed = NULL, show_progress = TRUE) {
  if (!is.null(seed)) set.seed(seed)   # RNG discipline: seed per replicate

  target <- make_log_posterior(y, p, q)

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
##   rep_start rep_end d_true n miss_pct [phi_true] [theta_true]
## Omitted or "NA" phi_true/theta_true => that term is absent, so:
##   ...                -> ARFIMA(0,d,0)
##   ... 0.5            -> ARFIMA(1,d,0)
##   ... NA 0.3         -> ARFIMA(0,d,1)
##   ... 0.5 0.3        -> ARFIMA(1,d,1)
## =====================================================================
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  usage <- paste0(
    "Usage: Rscript mcmc_benchmark_arfima_pdq_exact.R rep_start rep_end d_true n miss_pct [phi_true] [theta_true]\n",
    "       (omit or pass NA for phi_true / theta_true to drop that term)")
  if (length(args) < 5) stop(usage)

  rep_start <- as.integer(args[1])
  rep_end   <- as.integer(args[2])
  d_true    <- as.numeric(args[3])
  n         <- as.integer(args[4])
  miss_pct  <- as.numeric(args[5])
  phi_true   <- if (length(args) >= 6) suppressWarnings(as.numeric(args[6])) else NA_real_
  theta_true <- if (length(args) >= 7) suppressWarnings(as.numeric(args[7])) else NA_real_
  p <- as.integer(!is.na(phi_true))
  q <- as.integer(!is.na(theta_true))
  if (p == 1 && abs(phi_true)   >= 1) stop("phi_true must satisfy |phi| < 1 (causality)")
  if (q == 1 && abs(theta_true) >= 1) stop("theta_true must satisfy |theta| < 1 (invertibility)")

  s2_true <- 1
  dir.create(ROOT_DIR, recursive = TRUE, showWarnings = FALSE)

  arma_tag  <- paste0(if (p == 1) sprintf("_phi%.2f", phi_true) else "",
                      if (q == 1) sprintf("_theta%.2f", theta_true) else "")
  prior_tag <- if (SIGMA2_PRIOR == "invgamma") "_s2ig" else ""

  for (rep_id in rep_start:rep_end) {
    out_file <- file.path(ROOT_DIR, sprintf("mcmc_rep%04d_d%.2f%s_pctmis%.2f%s_exact.rds",
                                            rep_id, d_true, arma_tag, miss_pct, prior_tag))
    if (file.exists(out_file)) next   # resume-safe skip

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

    fit <- run_mcmc_one_rep(y_obs, p = p, q = q, seed = rep_id)
    fit$meta <- list(rep_id = rep_id, seed = seed, p = p, q = q,
                     d_true = d_true,
                     phi_true   = if (p == 1) phi_true   else NA_real_,
                     theta_true = if (q == 1) theta_true else NA_real_,
                     sigma2_true = s2_true, n = n, miss_pct = miss_pct,
                     likelihood = "exact",
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