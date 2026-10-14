################################################################################
#
# Comparison of inla_vanilla() with brms (Stan, NUTS) on multilevel models
# (October 2026). The random effects are integrated out analytically in each
# logpost, so inla_vanilla() sees only the fixed effects and the variance
# components. brms samples the full joint, random effects included, in its
# non-centred form, with the same priors. Long chains serve as the reference.
#
# Models: lme4's sleepstudy with correlated random intercepts and slopes (all
# 18 subjects, and the first 6 only, where the variance components are weakly
# identified), and the eight-schools model (Rubin, 1981), whose tau is
# notoriously skewed towards zero. Each parameter is compared on its natural
# scale: the gap in means and in quantiles, in reference posterior SDs, and
# the ratio of SDs. Needs brms (rstan backend) and lme4. Takes a few minutes,
# most of it compiling Stan. Run from the package root:
#   Rscript lab/multilevel-brms.R
#
################################################################################

## ----- Configuration ----------------------------------------------------------

suppressMessages(devtools::load_all(".", quiet = TRUE))
library(brms)
n_draws <- 20000 # per method: 4 chains x 5000 for brms
probs <- c(0.025, 0.5, 0.975)

## ----- Random intercept and slope model ---------------------------------------

# Marginal model per subject: y_j ~ N(X_j b, X_j G X_j' + sigma^2 I). Priors:
# b ~ N(0, 500^2), tau_0, tau_1, sigma ~ half-N(0, 50^2), rho ~ U(-1, 1), which
# is LKJ(1) for a 2 x 2 correlation matrix.
lmm_logpost <- function(theta, ys, xs) {
  b <- c(theta[["b0"]], theta[["b1"]])
  tau <- exp(c(theta[["log_tau0"]], theta[["log_tau1"]]))
  rho <- tanh(theta[["atanh_rho"]])
  sigma <- exp(theta[["log_sigma"]])
  g <- diag(tau^2)
  g[1, 2] <- g[2, 1] <- rho * tau[1] * tau[2]
  ll <- 0
  for (j in seq_along(ys)) {
    x <- xs[[j]]
    v <- x %*% g %*% t(x) + diag(sigma^2, nrow(x))
    r <- chol(v)
    z <- backsolve(r, ys[[j]] - x %*% b, transpose = TRUE)
    ll <- ll - sum(log(diag(r))) - 0.5 * sum(z^2) - nrow(x) / 2 * log(2 * pi)
  }
  lp_sd <- function(s) dnorm(s, 0, 50, log = TRUE) + log(2) + log(s)
  ll +
    sum(dnorm(b, 0, 500, log = TRUE)) +
    sum(lp_sd(tau)) +
    lp_sd(sigma) +
    log(0.5) +
    log1p(-rho^2)
}

lmm_natural <- function(theta) {
  c(
    b0 = theta[["b0"]],
    b1 = theta[["b1"]],
    tau0 = exp(theta[["log_tau0"]]),
    tau1 = exp(theta[["log_tau1"]]),
    rho = tanh(theta[["atanh_rho"]]),
    sigma = exp(theta[["log_sigma"]])
  )
}

lmm_case <- function(n_subjects) {
  dat <- lme4::sleepstudy
  keep <- levels(dat$Subject)[seq_len(n_subjects)]
  dat <- droplevels(dat[dat$Subject %in% keep, ])
  split_dat <- split(dat, dat$Subject)
  init_lm <- coef(lm(Reaction ~ Days, dat))
  list(
    logpost = lmm_logpost,
    init = c(
      b0 = init_lm[[1]],
      b1 = init_lm[[2]],
      log_tau0 = log(20),
      log_tau1 = log(5),
      atanh_rho = 0,
      log_sigma = log(sd(dat$Reaction) / 2)
    ),
    args = list(
      ys = lapply(split_dat, `[[`, "Reaction"),
      xs = lapply(split_dat, function(d) cbind(1, d$Days))
    ),
    natural = lmm_natural,
    formula = bf(Reaction ~ 0 + Intercept + Days + (Days | Subject)),
    data = dat,
    prior = c(
      prior(normal(0, 500), class = b),
      prior(normal(0, 50), class = sd),
      prior(normal(0, 50), class = sigma),
      prior(lkj(1), class = cor)
    ),
    brms_pars = c(
      "b_Intercept",
      "b_Days",
      "sd_Subject__Intercept",
      "sd_Subject__Days",
      "cor_Subject__Intercept__Days",
      "sigma"
    ),
    adapt_delta = 0.95
  )
}

## ----- Eight schools ----------------------------------------------------------

# Marginal model: y_j ~ N(mu, s_j^2 + tau^2). Priors: mu ~ N(0, 5^2),
# tau ~ half-Cauchy(0, 5).
schools_logpost <- function(theta, y, s) {
  mu <- theta[["mu"]]
  tau <- exp(theta[["log_tau"]])
  sum(dnorm(y, mu, sqrt(s^2 + tau^2), log = TRUE)) +
    dnorm(mu, 0, 5, log = TRUE) +
    log(2) +
    dcauchy(tau, 0, 5, log = TRUE) +
    log(tau)
}

schools <- data.frame(
  school = LETTERS[1:8],
  y = c(28, 8, -3, 7, -1, 1, 18, 12),
  s = c(15, 10, 16, 11, 9, 11, 10, 18)
)

schools_case <- list(
  logpost = schools_logpost,
  init = c(mu = mean(schools$y), log_tau = log(5)),
  args = list(y = schools$y, s = schools$s),
  natural = function(theta) {
    c(mu = theta[["mu"]], tau = exp(theta[["log_tau"]]))
  },
  formula = bf(y | se(s) ~ 1 + (1 | school)),
  data = schools,
  prior = c(
    prior(normal(0, 5), class = Intercept),
    prior(cauchy(0, 5), class = sd)
  ),
  brms_pars = c("b_Intercept", "sd_school__Intercept"),
  adapt_delta = 0.99
)

cases <- list(
  "sleepstudy, 18 subjects" = lmm_case(18),
  "sleepstudy, 6 subjects" = lmm_case(6),
  "eight schools" = schools_case
)

## ----- Runners ----------------------------------------------------------------

run_vanilla <- function(case) {
  fit <- do.call(
    inla_vanilla,
    c(list(logpost = case$logpost, init = case$init), case$args, silent = TRUE)
  )
  set.seed(1)
  draws <- sampling(fit, nsamp = n_draws, fun = case$natural)
  list(fit = fit, draws = draws)
}

run_brms <- function(case) {
  fit <- brm(
    case$formula,
    data = case$data,
    prior = case$prior,
    chains = 4,
    cores = 4,
    iter = 7500,
    warmup = 2500,
    seed = 1,
    refresh = 0,
    silent = 2,
    control = list(adapt_delta = case$adapt_delta, max_treedepth = 12)
  )
  # Wall time of sampling: the chains run in parallel, so the slowest chain
  elapsed <- rstan::get_elapsed_time(fit$fit)
  np <- nuts_params(fit)
  list(
    draws = as.matrix(fit, variable = case$brms_pars),
    time = max(rowSums(elapsed)),
    rhat = max(rhat(fit, pars = case$brms_pars)),
    n_div = sum(np$Value[np$Parameter == "divergent__"])
  )
}

compare_draws <- function(iv, ref) {
  colnames(ref) <- colnames(iv)
  ref_sd <- apply(ref, 2, sd)
  q_iv <- apply(iv, 2, quantile, probs)
  q_ref <- apply(ref, 2, quantile, probs)
  data.frame(
    brms_mean = colMeans(ref),
    iv_mean = colMeans(iv),
    brms_sd = ref_sd,
    iv_sd = apply(iv, 2, sd),
    mean_gap = (colMeans(iv) - colMeans(ref)) / ref_sd,
    sd_ratio = apply(iv, 2, sd) / ref_sd,
    q025_gap = (q_iv[1, ] - q_ref[1, ]) / ref_sd,
    q50_gap = (q_iv[2, ] - q_ref[2, ]) / ref_sd,
    q975_gap = (q_iv[3, ] - q_ref[3, ]) / ref_sd
  )
}

## ----- Run --------------------------------------------------------------------

for (nm in names(cases)) {
  case <- cases[[nm]]
  iv <- run_vanilla(case)
  br <- suppressWarnings(run_brms(case))
  cat("\n=====", nm, "=====\n")
  cat(sprintf(
    "inla_vanilla: %.2f s, %d logpost calls | brms: %.1f s sampling, max Rhat %.3f, %d divergences\n",
    iv$fit$timing[["total"]],
    iv$fit$counts[["logpost"]],
    br$time,
    br$rhat,
    br$n_div
  ))
  print(iv$fit)
  print(round(compare_draws(iv$draws, br$draws), 3))
}
