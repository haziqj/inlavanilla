################################################################################
#
# Binary CFA with the factors integrated out (October 2026): does
# inla_vanilla() recover the true values? Probit (normal-ogive) factor models
# for binary items, y_ij = 1 when lambda_j eta_i + e_ij > tau_j, with standard
# normal factors and errors (the theta parameterisation, as in lavaan). Each
# logpost integrates the factors out by Gauss-Hermite quadrature, so
# inla_vanilla() sees only the loadings, the thresholds and the factor
# correlation, and the inner integral runs inside every logpost call.
#
# Models: one factor with 6 items (m = 12), and two correlated factors with 4
# items each (m = 17, a 2-D inner integral). Priors: lambda_j, tau_j ~ N(0, 2^2)
# and rho ~ U(-1, 1), on atanh(rho). The loadings are left unconstrained: the
# fit starts at positive loadings and stays on that mode, far from the mirror
# image with every sign flipped.
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
n_reps <- if (smoke) 24 else 500
n_cores <- if (.Platform$OS.type == "windows") {
  1L
} else {
  max(1L, parallel::detectCores() - 2L)
}
n_draws <- if (smoke) 1000 else 5000 # per Stan chain, 4 chains
prior_sd <- 2

models <- list(
  "one factor" = list(
    item_factor = rep(1, 6),
    lambda = c(0.6, 0.8, 1.0, 1.2, 1.5, 0.7),
    tau = c(-1, -0.5, 0, 0.3, 0.8, 1.2),
    rho = NULL,
    n_gh = 21
  ),
  "two factors" = list(
    item_factor = rep(1:2, each = 4),
    lambda = c(0.6, 0.8, 1.0, 1.2, 1.5, 0.7, 0.9, 1.1),
    tau = c(-1, -0.5, 0, 0.3, 0.8, 1.2, -0.2, 0.5),
    rho = 0.4,
    n_gh = 15
  )
)

## =============================================================================
## Model
## =============================================================================

## ----- Gauss-Hermite quadrature ----------------------------------------------

# Nodes and weights of the n-point Gauss-Hermite rule for N(0, 1), by
# Golub-Welsch, and their product grid in k dimensions.
gh_normal <- function(n) {
  b <- sqrt(seq_len(n - 1))
  jac <- matrix(0, n, n)
  jac[cbind(seq_len(n - 1), 2:n)] <- b
  jac[cbind(2:n, seq_len(n - 1))] <- b
  e <- eigen(jac, symmetric = TRUE)
  list(z = e$values, w = e$vectors[1, ]^2)
}

gh_grid <- function(n, k) {
  g <- gh_normal(n)
  z <- as.matrix(expand.grid(rep(list(g$z), k)))
  w <- as.matrix(expand.grid(rep(list(g$w), k)))
  list(z = unname(z), log_w = rowSums(log(w)))
}

## ----- Log-posterior and gradient --------------------------------------------

theta_names <- function(model) {
  p <- length(model$lambda)
  c(
    paste0("lambda", seq_len(p)),
    paste0("tau", seq_len(p)),
    if (!is.null(model$rho)) "atanh_rho"
  )
}

true_theta <- function(model) {
  setNames(
    c(model$lambda, model$tau, if (!is.null(model$rho)) atanh(model$rho)),
    theta_names(model)
  )
}

# Shared by the log-posterior and its gradient. At quadrature node q the
# factors are eta_q = L z_q, with L the Cholesky factor of the factor
# correlation matrix, and item j has linear predictor a_qj = lambda_j eta_q,f(j)
# - tau_j. The log-likelihood of response pattern i at node q is then a matrix
# product, y log Phi(a) + (1 - y) log Phi(-a), and the pattern's marginal
# log-likelihood is a log-sum-exp over the nodes.
cfa_parts <- function(theta, y, item_factor, nodes) {
  p <- ncol(y)
  lambda <- theta[seq_len(p)]
  tau <- theta[p + seq_len(p)]
  z <- nodes$z
  rho <- if (ncol(z) == 2) tanh(theta[[2 * p + 1]])
  eta <- z
  if (!is.null(rho)) {
    eta[, 2] <- rho * z[, 1] + sqrt(1 - rho^2) * z[, 2]
  }
  n_q <- nrow(z)
  a <- eta[, item_factor, drop = FALSE] *
    matrix(lambda, n_q, p, byrow = TRUE) -
    matrix(tau, n_q, p, byrow = TRUE)
  log_p1 <- pnorm(a, log.p = TRUE)
  log_p0 <- pnorm(-a, log.p = TRUE)
  ell <- tcrossprod(y, log_p1) +
    tcrossprod(1 - y, log_p0) +
    matrix(nodes$log_w, nrow(y), n_q, byrow = TRUE)
  ell_max <- apply(ell, 1, max)
  list(
    lambda = lambda,
    tau = tau,
    rho = rho,
    eta = eta,
    a = a,
    log_p1 = log_p1,
    log_p0 = log_p0,
    ell = ell,
    ll = ell_max + log(rowSums(exp(ell - ell_max)))
  )
}

cfa_logpost <- function(theta, y, count, item_factor, nodes) {
  cp <- cfa_parts(theta, y, item_factor, nodes)
  lp <- sum(count * cp$ll) +
    sum(dnorm(cp$lambda, 0, prior_sd, log = TRUE)) +
    sum(dnorm(cp$tau, 0, prior_sd, log = TRUE))
  # rho ~ U(-1, 1), with the Jacobian of rho = tanh(atanh_rho)
  if (!is.null(cp$rho)) {
    lp <- lp + log1p(-cp$rho^2) - log(2)
  }
  lp
}

# The score of each pattern is the average over its posterior node weights of
# the complete-data score, whose derivative in a_qj is the signed Mills ratio
# of the probit.
cfa_grad <- function(theta, y, count, item_factor, nodes) {
  cp <- cfa_parts(theta, y, item_factor, nodes)
  post <- count * exp(cp$ell - cp$ll)
  log_phi <- dnorm(cp$a, log = TRUE)
  d_a <- crossprod(post, y) *
    exp(log_phi - cp$log_p1) -
    crossprod(post, 1 - y) * exp(log_phi - cp$log_p0)
  eta_j <- cp$eta[, item_factor, drop = FALSE]
  out <- c(
    colSums(d_a * eta_j) - cp$lambda / prior_sd^2,
    -colSums(d_a) - cp$tau / prior_sd^2
  )
  if (!is.null(cp$rho)) {
    z <- nodes$z
    d_eta2 <- z[, 1] - cp$rho / sqrt(1 - cp$rho^2) * z[, 2]
    on_2 <- item_factor == 2
    d_rho <- sum(d_a[, on_2, drop = FALSE] * outer(d_eta2, cp$lambda[on_2]))
    out <- c(out, d_rho * (1 - cp$rho^2) - 2 * cp$rho)
  }
  out
}

## ----- Simulation ------------------------------------------------------------

sim_cfa <- function(n, model) {
  p <- length(model$lambda)
  eta <- matrix(rnorm(n * max(model$item_factor)), n)
  if (!is.null(model$rho)) {
    eta[, 2] <- model$rho * eta[, 1] + sqrt(1 - model$rho^2) * eta[, 2]
  }
  y_star <- eta[, model$item_factor, drop = FALSE] *
    matrix(model$lambda, n, p, byrow = TRUE) +
    matrix(rnorm(n * p), n)
  y <- 1 * (y_star > matrix(model$tau, n, p, byrow = TRUE))
  colnames(y) <- paste0("y", seq_len(p))
  y
}

# Unique response patterns and their counts: at most 2^p rows whatever n
collapse_patterns <- function(y) {
  key <- apply(y, 1, paste, collapse = "")
  first <- !duplicated(key)
  list(
    y = y[first, , drop = FALSE],
    count = as.vector(table(key)[key[first]])
  )
}

fit_cfa <- function(y, model, ...) {
  pat <- collapse_patterns(y)
  nodes <- gh_grid(model$n_gh, max(model$item_factor))
  init <- true_theta(model)
  init[] <- 0
  init[startsWith(names(init), "lambda")] <- 1
  inla_vanilla(
    cfa_logpost,
    init,
    cfa_grad,
    y = pat$y,
    count = pat$count,
    item_factor = model$item_factor,
    nodes = nodes,
    ...
  )
}

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

## ----- Stan model ------------------------------------------------------------

# The full joint, factor scores sampled, with the same priors. The loadings
# are constrained positive, which picks the same mode as inla_vanilla()'s
# start: the true loadings are all over 4 posterior SDs from zero. For one
# factor, rho is sampled from its prior and ignored.
cfa_stan <- "
data {
  int<lower=1> N;
  int<lower=1> J;
  int<lower=1, upper=2> K;
  array[N, J] int<lower=0, upper=1> y;
  array[J] int<lower=1, upper=K> f;
  real<lower=0> prior_sd;
}
parameters {
  vector<lower=0>[J] lambda;
  vector[J] tau;
  real<lower=-1, upper=1> rho;
  matrix[N, K] z;
}
model {
  matrix[N, K] eta = z;
  if (K == 2) {
    eta[, 2] = rho * z[, 1] + sqrt(1 - square(rho)) * z[, 2];
  }
  to_vector(z) ~ std_normal();
  lambda ~ normal(0, prior_sd);
  tau ~ normal(0, prior_sd);
  for (j in 1:J) {
    y[, j] ~ bernoulli(Phi(lambda[j] * eta[, f[j]] - tau[j]));
  }
}
"

run_stan <- function(y, model, stan_model) {
  fit <- rstan::sampling(
    stan_model,
    data = list(
      N = nrow(y),
      J = ncol(y),
      K = max(model$item_factor),
      y = matrix(as.integer(y), nrow(y)),
      f = model$item_factor,
      prior_sd = prior_sd
    ),
    chains = 4,
    cores = 4,
    iter = 2 * n_draws,
    warmup = n_draws,
    seed = 1,
    refresh = 0,
    init = function() {
      list(lambda = rep(1, ncol(y)), tau = rep(0, ncol(y)), rho = 0)
    }
  )
  pars <- c("lambda", "tau", if (!is.null(model$rho)) "rho")
  draws <- as.matrix(fit, pars = pars)
  if (!is.null(model$rho)) {
    draws[, "rho"] <- atanh(draws[, "rho"])
  }
  colnames(draws) <- theta_names(model)
  summ <- rstan::summary(fit, pars = pars)$summary
  elapsed <- rstan::get_elapsed_time(fit)
  list(
    draws = draws,
    time = max(rowSums(elapsed)),
    rhat = max(summ[, "Rhat"]),
    n_div = sum(rstan::get_divergent_iterations(fit))
  )
}

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
  ref <- run_stan(y, model, stan_model)
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
