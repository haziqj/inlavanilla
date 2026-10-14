################################################################################
#
# Binary CFA in the delta parameterisation (October 2026), the sibling of
# lab/binary-cfa.R: does it give inla_vanilla() a better-behaved posterior?
# The delta parameterisation fixes Var(y*) = 1 instead of the residual
# variance, so its loadings and thresholds are standardised,
# lambda*_j = lambda_j / sqrt(1 + lambda_j^2) and
# tau*_j = tau_j / sqrt(1 + lambda_j^2). Since tau*_j = -qnorm(P(y_j = 1)),
# the threshold depends on the item's proportion alone, where the theta
# threshold scales with the loading.
#
# Same models, data and priors as lab/binary-cfa.R (lab/binary-cfa-model.R).
# The delta coordinates are atanh(lambda*_j), which equals asinh(lambda_j), and
# tau*_j. The log-posterior is the theta one at the mapped point plus the
# log-Jacobian, 2 sum_j log cosh(atanh(lambda*_j)). The posterior is therefore
# the same, and only the coordinates in which inla_vanilla() approximates it
# change, so one Stan run is the reference for both.
#
# Experiment 1 fits each model's data set in both coordinates and compares
# each fit with Stan on both scales, through sampling(). Experiment 2 fits both
# to the same simulated data sets and compares the diagnostics, the checks
# summary() flags and the coverage of the 95% intervals, each fit on its own
# scale (the loadings map one to one, so theirs are the same estimand). Needs
# rstan. Run from the package root:
#   Rscript lab/binary-cfa-delta.R
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
n_samp <- 20000 # copula draws per fit in Experiment 1

## =============================================================================
## Delta coordinates
## =============================================================================

delta_names <- function(model) {
  sub("^lambda", "atanh_lambda_d", sub("^tau", "tau_d", theta_names(model)))
}

true_delta <- function(model) {
  s <- sqrt(1 + model$lambda^2)
  setNames(
    c(
      asinh(model$lambda),
      model$tau / s,
      if (!is.null(model$rho)) atanh(model$rho)
    ),
    delta_names(model)
  )
}

# The same starting point as the theta fit
delta_init <- function(model) {
  init <- theta_init(model)
  p <- length(model$lambda)
  init[seq_len(p)] <- asinh(init[seq_len(p)])
  setNames(init, delta_names(model))
}

# lambda_j = sinh(u_j) and tau_j = tau*_j cosh(u_j), with
# u_j = atanh(lambda*_j). atanh_rho is the same in both.
delta_to_theta <- function(x, p) {
  u <- x[seq_len(p)]
  x[seq_len(p)] <- sinh(u)
  x[p + seq_len(p)] <- x[p + seq_len(p)] * cosh(u)
  x
}

log_cosh <- function(u) abs(u) + log1p(exp(-2 * abs(u))) - log(2)

delta_logpost <- function(theta, y, count, item_factor, nodes) {
  p <- ncol(y)
  cfa_logpost(delta_to_theta(theta, p), y, count, item_factor, nodes) +
    2 * sum(log_cosh(theta[seq_len(p)]))
}

# The chain rule through delta_to_theta(), plus the log-Jacobian's gradient
delta_grad <- function(theta, y, count, item_factor, nodes) {
  p <- ncol(y)
  u <- theta[seq_len(p)]
  tau_d <- theta[p + seq_len(p)]
  g <- cfa_grad(delta_to_theta(theta, p), y, count, item_factor, nodes)
  g_lambda <- g[seq_len(p)]
  g_tau <- g[p + seq_len(p)]
  g[seq_len(p)] <- g_lambda * cosh(u) + g_tau * tau_d * sinh(u) + 2 * tanh(u)
  g[p + seq_len(p)] <- g_tau * cosh(u)
  g
}

fit_coords <- function(y, model, coords, ...) {
  if (coords == "theta") {
    fit_cfa(y, model, ...)
  } else {
    fit_cfa(y, model, delta_logpost, delta_grad, delta_init(model), ...)
  }
}

## ----- Check of the gradient -------------------------------------------------

set.seed(10)
for (nm in names(models)) {
  model <- models[[nm]]
  pat <- collapse_patterns(sim_cfa(500, model))
  nodes <- gh_grid(model$n_gh, max(model$item_factor))
  theta <- true_delta(model) + rnorm(length(true_delta(model)), sd = 0.1)
  f <- function(t) delta_logpost(t, pat$y, pat$count, model$item_factor, nodes)
  fd <- vapply(
    seq_along(theta),
    function(j) {
      h <- replace(numeric(length(theta)), j, 1e-5)
      (f(theta + h) - f(theta - h)) / 2e-5
    },
    numeric(1)
  )
  gr <- delta_grad(theta, pat$y, pat$count, model$item_factor, nodes)
  cat(sprintf("%s: delta grad vs FD: %.1e\n", nm, max(abs(gr - fd))))
}

## =============================================================================
## Experiment 1 --- one data set, both coordinates, against Stan
## =============================================================================

## ----- Comparison ------------------------------------------------------------

# Theta-scale draws, with atanh_rho, to both scales, with rho
both_scales <- function(th, model) {
  p <- length(model$lambda)
  colnames(th) <- theta_names(model)
  lambda <- th[, seq_len(p), drop = FALSE]
  tau <- th[, p + seq_len(p), drop = FALSE]
  s <- sqrt(1 + lambda^2)
  rho <- if (!is.null(model$rho)) cbind(rho = tanh(th[, 2 * p + 1]))
  delta <- cbind(lambda / s, tau / s)
  colnames(delta) <- c(
    paste0("lambda_d", seq_len(p)),
    paste0("tau_d", seq_len(p))
  )
  list(theta = cbind(lambda, tau, rho), delta = cbind(delta, rho))
}

fit_draws <- function(fit, model, coords) {
  d <- sampling(fit, nsamp = n_samp)
  if (coords == "delta") {
    d <- t(apply(d, 1, delta_to_theta, p = length(model$lambda)))
  }
  both_scales(d, model)
}

# Gaps from the reference, in reference SDs: of the mean, and of the quantile
# furthest off among the 2.5%, 50% and 97.5%. Also the ratio of SDs.
gaps <- function(x, ref) {
  q <- c(0.025, 0.5, 0.975)
  ref_sd <- apply(ref, 2, sd)
  qg <- (apply(x, 2, quantile, q) - apply(ref, 2, quantile, q)) /
    rep(ref_sd, each = 3)
  cbind(
    mean = (colMeans(x) - colMeans(ref)) / ref_sd,
    sd = apply(x, 2, sd) / ref_sd,
    q = qg[cbind(apply(abs(qg), 2, which.max), seq_len(ncol(qg)))]
  )
}

side_by_side <- function(th, de) {
  out <- cbind(th[, 1], de[, 1], th[, 2], de[, 2], th[, 3], de[, 3])
  colnames(out) <- c("mean_th", "mean_de", "sd_th", "sd_de", "q_th", "q_de")
  round(out, 3)
}

fit_diag <- function(fit) {
  mg <- fit$diagnostics$marginals
  dg <- diagnostics_internal(fit$diagnostics$optim, fit$diagnostics$vb, mg)
  flags <- check_flags(dg, FALSE)
  data.frame(
    secs = round(fit$timing[["total"]], 2),
    calls = sum(fit$counts),
    hess_cond = signif(fit$diagnostics$optim$hess_cond, 3),
    vb_shift = round(max(abs(mg$vb_shift_sigma)), 3),
    nmad = signif(max(mg$nmad), 3),
    alpha = signif(max(abs(mg$alpha)), 3),
    ridge_bend = signif(max(mg$ridge_bend), 3),
    heavy_tails = sum(!mg$tails),
    flags = paste(names(flags)[flags %in% TRUE], collapse = " ")
  )
}

## ----- Run -------------------------------------------------------------------

# The data sets of lab/binary-cfa.R, simulated before sampling() draws from
# the same stream
stan_model <- rstan::stan_model(model_code = cfa_stan)
set.seed(2026)
data_sets <- lapply(models, sim_cfa, n = 500)
for (nm in names(models)) {
  model <- models[[nm]]
  y <- data_sets[[nm]]
  cat("\n=====", nm, "(n = 500) =====\n")
  fits <- list(
    theta = fit_coords(y, model, "theta", silent = TRUE),
    delta = fit_coords(y, model, "delta", silent = TRUE)
  )
  ref <- run_stan(y, model, stan_model, n_draws)
  cat(sprintf(
    "Stan: %.0f s sampling, max Rhat %.3f, %d divergences\n",
    ref$time,
    ref$rhat,
    ref$n_div
  ))
  print(do.call(rbind, lapply(fits, fit_diag)))
  ref_scales <- both_scales(ref$draws, model)
  draws <- list(
    theta = fit_draws(fits$theta, model, "theta"),
    delta = fit_draws(fits$delta, model, "delta")
  )
  for (sc in c("theta", "delta")) {
    cat(sprintf(
      "\nGaps from Stan on the %s scale (_th: theta fit, _de: delta fit)\n",
      sc
    ))
    print(side_by_side(
      gaps(draws$theta[[sc]], ref_scales[[sc]]),
      gaps(draws$delta[[sc]], ref_scales[[sc]])
    ))
  }
}

## =============================================================================
## Experiment 2 --- repeated sampling, both coordinates
## =============================================================================

## ----- Runner ----------------------------------------------------------------

rep_fit <- function(y, model, coords) {
  fit <- tryCatch(fit_coords(y, model, coords, silent = TRUE), error = identity)
  if (inherits(fit, "error")) {
    return(list(ok = FALSE))
  }
  truth <- if (coords == "theta") true_theta(model) else true_delta(model)
  s <- fit$summary
  mg <- fit$diagnostics$marginals
  dg <- diagnostics_internal(fit$diagnostics$optim, fit$diagnostics$vb, mg)
  list(
    ok = TRUE,
    covers = s[, "2.5%"] <= truth & truth <= s[, "97.5%"],
    checks = check_flags(dg, FALSE),
    param_flagged = apply(param_flags(dg, FALSE), 1, any),
    hess_cond = fit$diagnostics$optim$hess_cond,
    nmad = max(mg$nmad),
    ridge_bend = max(mg$ridge_bend),
    warned = length(fit$warnings) > 0,
    time = fit$timing[["total"]]
  )
}

# Seeded as in lab/binary-cfa.R, so the theta fits repeat its replicates
run_rep <- function(rep, model, n) {
  set.seed(rep)
  y <- sim_cfa(n, model)
  list(
    theta = rep_fit(y, model, "theta"),
    delta = rep_fit(y, model, "delta")
  )
}

summarise_coords <- function(reps) {
  ok <- vapply(reps, `[[`, logical(1), "ok")
  reps <- reps[ok]
  get <- function(key, type) vapply(reps, `[[`, type, key)
  checks <- vapply(reps, `[[`, logical(9), "checks")
  covers <- simplify2array(lapply(reps, `[[`, "covers"))
  flagged <- simplify2array(lapply(reps, `[[`, "param_flagged"))
  fired <- rowSums(checks, na.rm = TRUE)
  list(
    overall = data.frame(
      fitted = sum(ok),
      warned = sum(get("warned", logical(1))),
      flagged = sum(colSums(checks, na.rm = TRUE) > 0),
      cover_mean = round(mean(covers), 3),
      cover_min = round(min(rowMeans(covers)), 3),
      hess_cond = signif(median(get("hess_cond", numeric(1))), 3),
      nmad = signif(median(get("nmad", numeric(1))), 3),
      ridge_bend = signif(median(get("ridge_bend", numeric(1))), 3),
      secs = round(median(get("time", numeric(1))), 2)
    ),
    fired = fired[fired > 0],
    cover = rowMeans(covers),
    flag_rate = rowMeans(flagged)
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
  res <- lapply(
    c(theta = "theta", delta = "delta"),
    function(co) summarise_coords(lapply(reps, `[[`, co))
  )
  cat(sprintf(
    "\n===== %s, n = %d: %d replicates =====\n",
    st$model,
    st$n,
    n_reps
  ))
  cat("Medians of each fit's largest value: nmad, ridge_bend\n")
  print(do.call(rbind, lapply(res, `[[`, "overall")))
  for (co in names(res)) {
    fired <- res[[co]]$fired
    cat(sprintf(
      "Checks fired, %s: %s\n",
      co,
      if (length(fired) > 0) {
        paste(names(fired), fired, collapse = ", ")
      } else {
        "none"
      }
    ))
  }
  cat("Per parameter (delta rows: atanh(lambda*_j), tau*_j)\n")
  print(round(
    data.frame(
      cover_th = res$theta$cover,
      cover_de = res$delta$cover,
      flag_th = res$theta$flag_rate,
      flag_de = res$delta$flag_rate,
      row.names = theta_names(model)
    ),
    3
  ))
}
