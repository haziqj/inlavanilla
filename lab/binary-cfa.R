################################################################################
#
# Binary CFA with the factors integrated out (October 2026): does
# inla_vanilla() recover the true values? Probit factor models in the theta
# parameterisation, the factors integrated out by Gauss-Hermite quadrature
# inside every logpost call. The models, priors, log-posterior and Stan
# reference are in lab/binary-cfa-model.R: one factor with 6 items (m = 12),
# and two correlated factors with 4 items each (m = 17, a 2-D inner integral).
# lab/binary-cfa-delta.R repeats the study in the delta parameterisation.
#
# Experiment 1 fits one simulated data set per model and compares the
# marginals with Stan (NUTS) on the full joint, factor scores sampled, so the
# reference depends on neither the quadrature nor the approximation.
# Experiment 2 repeats the fit on many simulated data sets and reports bias,
# RMSE and the coverage of the 95% intervals, skew-normal and Gaussian. Needs
# rstan. Takes about 15 minutes on 12 cores. Run from the package root:
#   Rscript lab/binary-cfa.R
#
# Results (October 2026, 500 replicates per setting): the 95% intervals cover
# the true values at the nominal rate. Mean coverage over the parameters was
# 0.944 to 0.947 in the three settings, every parameter within 0.93 to 0.97
# (Monte Carlo SE 0.01). The large loadings are biased upwards (lambda5 = 1.5:
# +0.12 at n = 500), the right skew of their posteriors, which Stan shares.
# Against Stan, one factor: mean gaps under 0.03 SDs, SD ratios 0.95 to 1.01.
# Two factors: the same except item 5, whose lambda and tau have SDs 20% too
# small and 97.5% quantiles 0.7 SDs short, both flagged (heavy tails, ridge
# bend 0.47). The ridge bend flag fires in 37% to 93% of the replicates,
# mostly on tau5, with no loss of coverage.
#
################################################################################

## ----- Configuration ---------------------------------------------------------

smoke <- FALSE
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("lab/binary-cfa-model.R")
n_reps <- if (smoke) 24 else 500
n_cores <- if (.Platform$OS.type == "windows") {
  1L
} else {
  min(4L, parallel::detectCores() - 1L) # each fork holds its own memory
}
n_draws <- if (smoke) 1000 else 5000 # per Stan chain, 4 chains

## ----- Checks of the quadrature and the gradient -----------------------------

set.seed(10)
for (nm in names(models)) {
  model <- models[[nm]]
  pat <- collapse_patterns(sim_cfa(500, model))
  k <- max(model$item_factor)
  theta <- true_theta(model) + rnorm(length(true_theta(model)), sd = 0.1)
  lp_at <- function(n_gh) {
    cfa_logpost(theta, pat$y, pat$count, model$item_factor, gh_grid(n_gh, k))
  }
  nodes <- gh_grid(model$n_gh, k)
  fd <- vapply(
    seq_along(theta),
    function(j) {
      h <- replace(numeric(length(theta)), j, 1e-5)
      f <- function(t) {
        cfa_logpost(t, pat$y, pat$count, model$item_factor, nodes)
      }
      (f(theta + h) - f(theta - h)) / 2e-5
    },
    numeric(1)
  )
  gr <- cfa_grad(theta, pat$y, pat$count, model$item_factor, nodes)
  cat(sprintf(
    "%s: %d patterns | logpost error of %d nodes vs %d: %.1e | grad vs FD: %.1e\n",
    nm,
    nrow(pat$y),
    model$n_gh^k,
    (3 * model$n_gh)^k,
    lp_at(model$n_gh) - lp_at(3 * model$n_gh),
    max(abs(gr - fd))
  ))
}

## =============================================================================
## Experiment 1 --- one data set against Stan
## =============================================================================

## ----- Comparison ------------------------------------------------------------

# Each parameter on the theta scale: the gaps in the mean and the quantiles in
# Stan posterior SDs, and the ratio of SDs, with the true value and whether
# each 95% interval holds it.
compare_stan <- function(fit, ref, truth) {
  s <- fit$summary
  ref_sd <- apply(ref, 2, sd)
  q_ref <- apply(ref, 2, quantile, c(0.025, 0.5, 0.975))
  data.frame(
    true = truth,
    iv_mean = s[, "Mean"],
    stan_mean = colMeans(ref),
    mean_gap = (s[, "Mean"] - colMeans(ref)) / ref_sd,
    sd_ratio = s[, "SD"] / ref_sd,
    q025_gap = (s[, "2.5%"] - q_ref[1, ]) / ref_sd,
    q975_gap = (s[, "97.5%"] - q_ref[3, ]) / ref_sd,
    iv_covers = s[, "2.5%"] <= truth & truth <= s[, "97.5%"],
    stan_covers = q_ref[1, ] <= truth & truth <= q_ref[3, ]
  )
}

## ----- Run -------------------------------------------------------------------

stan_model <- rstan::stan_model(model_code = cfa_stan)
set.seed(2026)
for (nm in names(models)) {
  model <- models[[nm]]
  y <- sim_cfa(500, model)
  cat("\n=====", nm, "(n = 500) =====\n")
  fit <- fit_cfa(y, model)
  ref <- run_stan(y, model, stan_model, n_draws)
  cat(sprintf(
    "inla_vanilla: %.2f s, %d logpost + %d grad calls | Stan: %.0f s sampling, max Rhat %.3f, %d divergences\n",
    fit$timing[["total"]],
    fit$counts[["logpost"]],
    fit$counts[["grad"]],
    ref$time,
    ref$rhat,
    ref$n_div
  ))
  print(summary(fit))
  cmp <- compare_stan(fit, ref$draws, true_theta(model))
  num <- vapply(cmp, is.numeric, logical(1))
  cmp[num] <- lapply(cmp[num], round, 3)
  print(cmp)
}

## =============================================================================
## Experiment 2 --- repeated sampling
## =============================================================================

## ----- Runner ----------------------------------------------------------------

# One replicate: simulate, fit, and keep each parameter's posterior mean and
# SD, the skew-normal 95% interval, the Gaussian one, N(mean_vb, Sigma_jj),
# which marginals = "gaussian" would report, and the checks summary() flags.
run_rep <- function(rep, model, n) {
  set.seed(rep)
  y <- sim_cfa(n, model)
  fit <- tryCatch(fit_cfa(y, model, silent = TRUE), error = identity)
  if (inherits(fit, "error")) {
    return(list(ok = FALSE))
  }
  s <- fit$summary
  sd_g <- sqrt(diag(fit$cov))
  dg <- diagnostics_internal(
    fit$diagnostics$optim,
    fit$diagnostics$vb,
    fit$diagnostics$marginals
  )
  list(
    ok = TRUE,
    est = cbind(
      mean = s[, "Mean"],
      sd = s[, "SD"],
      lo = s[, "2.5%"],
      hi = s[, "97.5%"],
      lo_g = fit$mean_vb - qnorm(0.975) * sd_g,
      hi_g = fit$mean_vb + qnorm(0.975) * sd_g
    ),
    warned = length(fit$warnings) > 0,
    checks = check_flags(dg, FALSE),
    param_flagged = apply(param_flags(dg, FALSE), 1, any),
    time = fit$timing[["total"]]
  )
}

summarise_reps <- function(reps, truth) {
  ok <- vapply(reps, `[[`, logical(1), "ok")
  reps <- reps[ok]
  est <- simplify2array(lapply(reps, `[[`, "est"))
  covers <- function(lo, hi) {
    rowMeans(est[, lo, ] <= truth & truth <= est[, hi, ])
  }
  checks <- vapply(reps, `[[`, logical(9), "checks")
  tab <- data.frame(
    true = truth,
    bias = rowMeans(est[, "mean", ]) - truth,
    rmse = sqrt(rowMeans((est[, "mean", ] - truth)^2)),
    avg_sd = rowMeans(est[, "sd", ]),
    emp_sd = apply(est[, "mean", ], 1, sd),
    cover_sn = covers("lo", "hi"),
    cover_gauss = covers("lo_g", "hi_g"),
    flag_rate = rowMeans(vapply(
      reps,
      `[[`,
      logical(length(truth)),
      "param_flagged"
    ))
  )
  list(
    tab = tab,
    n_ok = sum(ok),
    n_warned = sum(vapply(reps, `[[`, logical(1), "warned")),
    n_flagged = sum(colSums(checks, na.rm = TRUE) > 0),
    checks = rowSums(checks, na.rm = TRUE),
    time = median(vapply(reps, `[[`, numeric(1), "time"))
  )
}

## ----- Run -------------------------------------------------------------------

settings <- list(
  list(model = "one factor", n = 200),
  list(model = "one factor", n = 500),
  list(model = "two factors", n = 500)
)

for (st in settings) {
  model <- models[[st$model]]
  reps <- parallel::mclapply(
    seq_len(n_reps),
    run_rep,
    model = model,
    n = st$n,
    mc.cores = n_cores
  )
  res <- summarise_reps(reps, true_theta(model))
  mc_se <- sqrt(0.95 * 0.05 / res$n_ok)
  cat(sprintf(
    "\n===== %s, n = %d: %d replicates =====\n",
    st$model,
    st$n,
    n_reps
  ))
  cat(sprintf(
    "%d fitted, %d with a warning, %d with a summary() flag | median fit %.2f s\n",
    res$n_ok,
    res$n_warned,
    res$n_flagged,
    res$time
  ))
  fired <- res$checks[res$checks > 0]
  if (length(fired) > 0) {
    cat(
      "Fits flagged by check:",
      paste(names(fired), fired, collapse = ", "),
      "\n"
    )
  }
  cat(sprintf(
    "Mean coverage of the 95%% intervals: skew-normal %.3f, Gaussian %.3f (Monte Carlo SE per parameter %.3f)\n",
    mean(res$tab$cover_sn),
    mean(res$tab$cover_gauss),
    mc_se
  ))
  print(round(res$tab, 3))
}
