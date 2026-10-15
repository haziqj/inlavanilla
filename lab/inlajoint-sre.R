################################################################################
#
# A joint model of a longitudinal marker and a survival time with shared random
# effects, from Rustand's INLA book, section 4.3.1 (Chapter 4,
# https://rustand.fr/INLA_book/C4/chapter4.html#sec-SREind), fitted with
# INLAjoint and with inla_vanilla() (October 2026). The data are the book's
# simulation: 500 subjects, 2,006 measurements and 219 events. The book prints
# no seed, so this script sets its own.
#
# The model, for subject i at time t:
#   y_i(t) = beta0 + b_i0 + (beta1 + b_i1) t + beta2 trt_i + beta3 trt_i t + e,
#   b_i ~ N(0, Sigma_b), e ~ N(0, sigma_e^2),
#   h_i(t) = h_0(t) exp(gamma0 b_i0 + gamma1 b_i1),
# with log h_0 piecewise constant on 15 intervals (16 nodes, an rw1).
# INLAjoint's defaults (version 26.6.2) are the priors, matched exactly:
# beta ~ N(0, precision 0.01), 1 / sigma_e^2 ~ Gamma(1, 5e-5), the precision of
# b_i ~ Wishart(3, I), gamma ~ N(0, 1), and the 16 log-hazards
# f ~ N(0, (tau R + 0.01 I)^-1), with R the scaled rw1 structure, no
# constraint, no separate intercept, and a PC prior on tau with
# P(1 / sqrt(tau) > 0.5) = 0.01. The book's printed fit used Wishart(10, I),
# from an older INLAjoint.
#
# The 1,000 random effects integrate out inside logpost. Given the parameters,
# b_i | y_i is normal in closed form, and the survival term sees b_i only
# through u_i = gamma' b_i, which is normal too. So each subject needs one
# integral over u_i, by adaptive Gauss-Hermite quadrature. Two versions treat
# the log-hazards differently.
#
# Version A gives inla_vanilla() all 27 parameters, with f non-centred,
# f = U diag(1 / sqrt(tau lambda + 0.01)) eta and eta ~ N(0, I). Centred, the
# joint density grows like tau^7 as tau goes to infinity, so it has no mode.
#
# Version B integrates f out by a Laplace approximation inside logpost, as
# R-INLA does for its latent field, with Newton steps from analytic
# derivatives. That leaves 11 parameters, and sampling() draws f from its
# Laplace approximation given each draw.
#
# Reference: Stan (NUTS) on all 1,027 parameters, and importance sampling for
# the log marginal likelihood.
#
# Needs INLA, INLAjoint, mvtnorm, rstan and loo. Run from the package root:
#   Rscript lab/inlajoint-sre.R
#
# Results (October 2026; Stan 4 chains of 2,000 kept draws, R-hat 1.003, no
# divergences). Version B matches Stan, as closely as INLAjoint or closer. The
# 10 parameters other than the baseline variance are within 0.02 posterior SDs
# in mean and 0.10 in the 2.5% and 97.5% quantiles, with SDs within 3.1%
# (INLAjoint: 0.08, 0.17 and 3.1%, part of it Monte Carlo noise in its summary
# of the random-effect SDs). The 15 log-hazards are within 0.07 SDs in mean
# and 0.15 in the quantiles, with SDs within 3% (INLAjoint: 0.10, 0.37, and 3%
# to 8% too wide). The baseline variance 1 / tau has mean 0.0092 against
# Stan's 0.0089, and log tau a 97.5% quantile of 11.4 against Stan's 12.4
# (INLAjoint: 10.0 to 10.8). INLAjoint's marginal of log tau changes from run
# to run on several threads: its 2.5% quantile ranged from 2.0 to 3.0 over
# six fits, against Stan's 2.9. The log marginal likelihood is -3336.13
# against -3336.11 by importance sampling (Pareto k 0.33). INLAjoint's is
# -3336.08 from its Gaussian approximation and -3341.07 from its integration,
# once sum(log E) = -472.38 and the rw1 offset of -4.76 are removed. Version B
# took 3.2 s (2,671 calls of logpost, 1.2 ms each), and its draws 1.6 ms each
# for the inner Laplace (16 s for 10,000). INLAjoint took 2.9 s on 14 threads
# and 8.9 s on one. Every check in summary() passes except heavy tails for log
# tau, which is right: its right tail is the PC prior's.
#
# Version A gets the same 10 parameters right (within 0.02 SDs in mean and
# 0.11 in the quantiles) in 4.9 s (19,385 calls), but not the baseline. Log
# tau has an exponential right tail and the eta's a funnel, which skew-normal
# marginals cannot follow. Its draws put 1 / tau at 0, and the log-hazard SDs
# come out 15% to 32% too small. Its warning names log tau and three eta's.
#
# Two findings about the comparison itself. With sdcor = TRUE, INLAjoint's row
# "Baseline risk (variance)" holds 1 / sqrt(tau), so the book's 0.1015 is an
# SD. R-INLA normalises the rw1 with an intrinsic model's constant even with
# diagonal = 0.01: on a normal model with every hyperparameter fixed, its log
# marginal likelihood is 4.76 below the exact one, at any tau.
#
################################################################################

## ----- Configuration ---------------------------------------------------------

smoke <- FALSE
suppressMessages({
  devtools::load_all(".", quiet = TRUE)
  library(INLA)
  library(INLAjoint)
})
n_iter <- if (smoke) 1000 else 4000 # Stan iterations per chain, half warm-up
n_draws <- if (smoke) 2000 else 10000 # copula draws
n_is <- if (smoke) 2000 else 20000 # importance draws for the log marginal lik
n_nodes <- 15 # adaptive Gauss-Hermite nodes per subject
n_newton <- 8 # Newton steps for each subject's integrand mode
options(mc.cores = 4)

## =============================================================================
## Data
## =============================================================================

# The book's simulation, with its parameter values and the same random draws
# in the same order
set.seed(2026)
n <- 500
fw <- 5
mestime <- seq(0, fw, 1)
n_i <- length(mestime)
time <- rep(mestime, n)
id <- rep(1:n, each = n_i)
sigma_true <- matrix(c(1, 0.42, 0.42, 0.49), 2) # SDs 1 and 0.7, correlation 0.6
re <- mvtnorm::rmvnorm(n, rep(0, 2), sigma_true)
trt <- rbinom(n, 1, 0.5)
trt_y <- rep(trt, each = n_i)
lin_pred <- 1 +
  rep(re[, 1], each = n_i) +
  (-1 + rep(re[, 2], each = n_i)) * time +
  trt_y +
  0.5 * trt_y * time
y_all <- rnorm(n_i * n, lin_pred, 0.5)
u <- runif(n)
event_times <- -log(u) / (0.1 * exp(0.5 * re[, 1] + 0.3 * re[, 2]))
event_indic <- as.numeric(event_times < fw)
event_times <- pmin(event_times, fw)
dat_long <- data.frame(id, time, Y = y_all, trt = trt_y)
dat_long <- dat_long[dat_long$time < event_times[dat_long$id], ]
dat_surv <- data.frame(id = 1:n, event_times, event_indic, trt)

## =============================================================================
## INLAjoint
## =============================================================================

fit_inlajoint <- function() {
  joint(
    formSurv = inla.surv(event_times, event_indic) ~ -1,
    formLong = Y ~ time * trt + (1 + time | id),
    dataSurv = dat_surv,
    dataLong = dat_long,
    id = "id",
    timeVar = "time",
    assoc = "SRE_ind"
  )
}
t_ij <- system.time(fit_ij <- fit_inlajoint())[["elapsed"]]
inla.setOption(num.threads = "1:1")
t_ij_1 <- system.time(fit_inlajoint())[["elapsed"]]
inla.setOption(num.threads = NULL)
# With sdcor = TRUE, INLAjoint's row "Baseline risk (variance)" holds
# 1 / sqrt(tau), the SD, so the variance comes from the default summary
sum_ij <- summary(fit_ij, sdcor = TRUE)
sum_ij_var <- summary(fit_ij)

## =============================================================================
## Model for inla_vanilla()
## =============================================================================

## ----- Longitudinal sufficient statistics ------------------------------------

# Per subject: n_i, sum t, sum t^2, sum y, sum t y and sum y^2. The fixed part
# is a_i + c_i t, with a_i = beta0 + beta2 trt_i and c_i = beta1 + beta3 trt_i,
# so these give every residual sum the likelihood needs.
by_id <- function(x) rowsum(x, dat_long$id)[, 1]
st <- list(
  n = by_id(rep(1, nrow(dat_long))),
  t = by_id(dat_long$time),
  tt = by_id(dat_long$time^2),
  y = by_id(dat_long$Y),
  ty = by_id(dat_long$time * dat_long$Y),
  yy = by_id(dat_long$Y^2)
)

## ----- Baseline hazard -------------------------------------------------------

# INLAjoint's 15 intervals on [0, 5] and their 16 nodes. The hazard on
# [cut_k, cut_k+1) is exp(f_k), and node 16 carries no exposure.
cut <- seq(0, fw, length.out = 16)
n_cut <- length(cut)
expo <- outer(event_times, cut, function(t, s) {
  pmax(0, pmin(t, s + fw / 15) - s)
})
expo[, n_cut] <- 0
event_node <- findInterval(event_times, cut)

# The scaled rw1 structure, as INLA builds it, and its eigenvectors, with the
# constant vector (eigenvalue 0) first and positive
rw_r <- inla.scale.model(
  crossprod(diff(diag(n_cut))),
  constr = list(A = matrix(1, 1, n_cut), e = 0)
)
rw_eig <- eigen(as.matrix(rw_r), symmetric = TRUE)
rw_order <- rev(seq_len(n_cut))
rw_u <- rw_eig$vectors[, rw_order]
rw_lam <- c(0, rw_eig$values[rw_order][-1])
rw_u[, 1] <- abs(rw_u[, 1])
rw_diag <- 0.01
pc_lambda <- -log(0.01) / 0.5

# The log-hazards at the 16 nodes from eta and log tau
log_hazard <- function(eta, log_tau) {
  drop(rw_u %*% (eta / sqrt(exp(log_tau) * rw_lam + rw_diag)))
}

## ----- Gauss-Hermite rule ----------------------------------------------------

# Nodes and weights for the standard normal, from the Golub-Welsch algorithm
gh_normal <- function(k) {
  j <- seq_len(k - 1)
  jacobi <- matrix(0, k, k)
  jacobi[cbind(j, j + 1)] <- sqrt(j)
  jacobi[cbind(j + 1, j)] <- sqrt(j)
  e <- eigen(jacobi, symmetric = TRUE)
  list(z = e$values, w = e$vectors[1, ]^2)
}
set_rule <- function(k) {
  gh <- gh_normal(k)
  list(z = gh$z, log_w = log(gh$w) - dnorm(gh$z, log = TRUE), k = k)
}
rule <- set_rule(n_nodes)

## ----- Log-posterior pieces --------------------------------------------------

# The longitudinal log-likelihood with b_i integrated out, the priors of the
# longitudinal and association parameters, and the moments of u_i = gamma' b_i
# given y_i. Given the parameters, y_i is normal, and so is b_i | y_i, with
# precision P_i = Z_i' Z_i / sigma_e^2 + Sigma_b^-1 and mean m_i. The Wishart
# density of W = Sigma_b^-1 with r = 3 and R = I is exp(-tr(W) / 2) / (4 pi),
# and the Jacobian from (log sd0, log sd1, atanh rho) to W is
# 4 / ((sd0 sd1)^3 (1 - rho^2)^2).
long_part <- function(theta) {
  sig2 <- exp(2 * theta[["log_sigma_e"]])
  sd0 <- exp(theta[["log_sd_b0"]])
  sd1 <- exp(theta[["log_sd_b1"]])
  rho <- tanh(theta[["atanh_rho"]])
  a <- theta[["beta0"]] + theta[["beta2"]] * trt
  c <- theta[["beta1"]] + theta[["beta3"]] * trt
  s_r <- st$y - a * st$n - c * st$t
  s_tr <- st$ty - a * st$t - c * st$tt
  s_rr <- st$yy -
    2 * a * st$y -
    2 * c * st$ty +
    a^2 * st$n +
    2 * a * c * st$t +
    c^2 * st$tt
  one_rho <- 1 - rho^2
  w11 <- 1 / (sd0^2 * one_rho)
  w22 <- 1 / (sd1^2 * one_rho)
  w12 <- -rho / (sd0 * sd1 * one_rho)
  p11 <- st$n / sig2 + w11
  p12 <- st$t / sig2 + w12
  p22 <- st$tt / sig2 + w22
  det_p <- p11 * p22 - p12^2
  g1 <- s_r / sig2
  g2 <- s_tr / sig2
  m1 <- (p22 * g1 - p12 * g2) / det_p
  m2 <- (p11 * g2 - p12 * g1) / det_p
  ll <- sum(
    -0.5 *
      st$n *
      log(2 * pi * sig2) -
      0.5 * log(det_p) -
      0.5 * s_rr / sig2 +
      0.5 * (g1 * m1 + g2 * m2)
  ) -
    0.5 * n * (2 * log(sd0 * sd1) + log(one_rho))
  ga <- theta[["gamma0"]]
  gb <- theta[["gamma1"]]
  prior <- sum(dnorm(theta[1:4], 0, 10, log = TRUE)) +
    dgamma(1 / sig2, 1, 5e-5, log = TRUE) +
    log(2 / sig2) -
    (w11 + w22) / 2 -
    log(4 * pi) +
    log(4) -
    3 * log(sd0 * sd1) -
    2 * log(one_rho) +
    sum(dnorm(c(ga, gb), 0, 1, log = TRUE))
  list(
    ll = ll,
    prior = prior,
    mu = ga * m1 + gb * m2,
    s = sqrt((ga^2 * p22 - 2 * ga * gb * p12 + gb^2 * p11) / det_p)
  )
}

# The survival terms given the log-hazards f. With u_i ~ N(mu_i, s_i^2),
# subject i's term is the integral of exp(delta (f_k + u) - H_i e^u) over u,
# with H_i the baseline cumulative hazard at the subject's time. With
# u = mu + s z, the integrand's mode in z comes from Newton steps (it is concave
# with curvature below -1), then adaptive Gauss-Hermite around it. Returns the
# log-integrals and the mean and variance of e^u under each normalised
# integrand, which give the derivatives in f.
surv_terms <- function(f, mu, s) {
  h0 <- drop(expo %*% exp(f))
  f_ev <- event_indic * f[event_node]
  z <- numeric(n)
  for (k in seq_len(n_newton)) {
    e <- h0 * exp(mu + s * z)
    z <- z + (event_indic * s - s * e - z) / (s^2 * e + 1)
  }
  e <- h0 * exp(mu + s * z)
  omega <- 1 / sqrt(s^2 * e + 1)
  top <- event_indic * (f_ev + mu + s * z) - e - z^2 / 2
  zk <- z + outer(omega, rule$z)
  eu <- exp(mu + s * zk)
  wt <- exp(
    event_indic *
      (f_ev + mu + s * zk) -
      h0 * eu -
      zk^2 / 2 -
      top +
      rep(rule$log_w, each = n)
  )
  tot <- rowSums(wt)
  m_eu <- rowSums(wt * eu) / tot
  list(
    log_int = log(omega) + top - 0.5 * log(2 * pi) + log(tot),
    m_eu = m_eu,
    v_eu = pmax(rowSums(wt * eu^2) / tot - m_eu^2, 0)
  )
}

log_pc_prior <- function(log_tau) {
  log(pc_lambda / 2) - 0.5 * log_tau - pc_lambda * exp(-0.5 * log_tau)
}

## ----- Version A: all 27 parameters ------------------------------------------

# The log-hazards non-centred, f = U diag(1 / sqrt(tau lambda + 0.01)) eta with
# eta ~ N(0, I)
eta_idx <- 10 + seq_len(n_cut)
init_a <- c(
  beta0 = 1,
  beta1 = -1,
  beta2 = 0,
  beta3 = 0,
  log_sigma_e = log(0.5),
  log_sd_b0 = 0,
  log_sd_b1 = 0,
  atanh_rho = 0,
  gamma0 = 0,
  gamma1 = 0,
  setNames(numeric(n_cut), sprintf("eta%02d", seq_len(n_cut))),
  log_tau = log(pc_lambda^2)
)
# The level of the log-hazards, mean(f) = eta01 / sqrt(16 rw_diag), at the
# crude event rate
log_rate <- log(sum(event_indic) / sum(event_times))
init_a[["eta01"]] <- log_rate * sqrt(n_cut * rw_diag)

logpost_a <- function(theta) {
  lp <- long_part(theta)
  f <- log_hazard(theta[eta_idx], theta[["log_tau"]])
  lp$ll +
    sum(surv_terms(f, lp$mu, lp$s)$log_int) +
    lp$prior +
    sum(dnorm(theta[eta_idx], log = TRUE)) +
    log_pc_prior(theta[["log_tau"]])
}

## ----- Version B: the log-hazards integrated out -----------------------------

# A Laplace approximation over the 16 log-hazards f given the other 11
# parameters, as R-INLA does for its latent field. The survival terms are
# log-concave in f, so Newton's method from the last call's optimum converges
# in a few steps. With v_ik = E_ik exp(f_k), the gradient of the survival
# log-likelihood is the event count per interval minus sum_i v_ik E[e^u_i],
# and its Hessian is -diag(sum_i v_ik E[e^u_i]) + sum_i Var[e^u_i] v_i v_i'.
init_b <- init_a[c(1:10, length(init_a))]
event_count <- tabulate(event_node[event_indic == 1], n_cut)
inner <- new.env()
inner$f <- rep(log_rate, n_cut)

inner_laplace <- function(theta) {
  lp <- long_part(theta)
  q_eig <- exp(theta[["log_tau"]]) * rw_lam + rw_diag
  q <- rw_u %*% (q_eig * t(rw_u))
  f <- inner$f
  for (it in 1:50) {
    sv <- surv_terms(f, lp$mu, lp$s)
    v <- expo * rep(exp(f), each = n)
    rate <- colSums(v * sv$m_eu)
    prec <- diag(rate) - crossprod(v * sqrt(sv$v_eu)) + q
    step <- solve(prec, event_count - rate - drop(q %*% f))
    f <- f + step
    if (max(abs(step)) < 1e-9) break
  }
  sv <- surv_terms(f, lp$mu, lp$s)
  v <- expo * rep(exp(f), each = n)
  prec <- diag(colSums(v * sv$m_eu)) - crossprod(v * sqrt(sv$v_eu)) + q
  inner$f <- f
  list(lp = lp, sv = sv, f = f, q = q, q_eig = q_eig, prec = prec)
}

logpost_b <- function(theta) {
  il <- inner_laplace(theta)
  il$lp$ll +
    sum(il$sv$log_int) -
    0.5 * sum(il$f * (il$q %*% il$f)) +
    0.5 * sum(log(il$q_eig)) -
    0.5 * as.numeric(determinant(il$prec)$modulus) +
    il$lp$prior +
    log_pc_prior(theta[["log_tau"]])
}

## ----- Natural-scale quantities ----------------------------------------------

# Named as in INLAjoint's summary, then the log-hazards of the 15 intervals.
# In version B, f is drawn from its Laplace approximation given the draw.
natural <- function(theta, f) {
  c(
    beta0 = theta[["beta0"]],
    beta1 = theta[["beta1"]],
    beta2 = theta[["beta2"]],
    beta3 = theta[["beta3"]],
    sigma_e = exp(theta[["log_sigma_e"]]),
    sd_b0 = exp(theta[["log_sd_b0"]]),
    sd_b1 = exp(theta[["log_sd_b1"]]),
    rho = tanh(theta[["atanh_rho"]]),
    var_h0 = exp(-theta[["log_tau"]]),
    gamma0 = theta[["gamma0"]],
    gamma1 = theta[["gamma1"]],
    setNames(f[-n_cut], sprintf("f%02d", seq_len(n_cut - 1)))
  )
}
natural_a <- function(theta) {
  natural(theta, log_hazard(theta[eta_idx], theta[["log_tau"]]))
}
natural_b <- function(theta) {
  il <- inner_laplace(theta)
  natural(theta, il$f + backsolve(chol(il$prec), rnorm(n_cut)))
}

## =============================================================================
## inla_vanilla()
## =============================================================================

t_a <- system.time(fit_a <- inla_vanilla(logpost_a, init_a))[["elapsed"]]
set.seed(1)
t_a_draws <- system.time(
  draws_a <- sampling(fit_a, nsamp = n_draws, fun = natural_a)
)[["elapsed"]]

t_b <- system.time(fit_b <- inla_vanilla(logpost_b, init_b))[["elapsed"]]
set.seed(1)
t_b_draws <- system.time(
  draws_b <- sampling(fit_b, nsamp = n_draws, fun = natural_b)
)[["elapsed"]]

# Quadrature check at the mode: 15 against 41 nodes
rule <- set_rule(41)
lp_41 <- logpost_a(fit_a$mode)
rule <- set_rule(n_nodes)
quad_error <- lp_41 - logpost_a(fit_a$mode)

## =============================================================================
## Reference
## =============================================================================

## ----- Stan ------------------------------------------------------------------

stan_code <- "
data {
  int<lower=1> N;
  int<lower=1> n;
  int<lower=1> K;
  array[N] int<lower=1, upper=n> id;
  vector[N] time;
  vector[N] y;
  vector[n] trt;
  vector[n] delta;
  array[n] int<lower=1, upper=K> event_node;
  matrix[n, K] expo;
  matrix[K, K] rw_u;
  vector[K] rw_lam;
  real rw_diag;
  real pc_lambda;
}
parameters {
  vector[4] beta;
  real log_sigma_e;
  cov_matrix[2] sigma_b;
  matrix[n, 2] b;
  vector[2] gam;
  vector[K] eta;
  real log_tau;
}
transformed parameters {
  vector[K] f = rw_u * (eta ./ sqrt(exp(log_tau) * rw_lam + rw_diag));
}
model {
  vector[n] u = b * gam;
  y ~ normal(
    beta[1] + b[id, 1] + (beta[2] + b[id, 2]) .* time +
      (beta[3] + beta[4] * time) .* trt[id],
    exp(log_sigma_e)
  );
  target += dot_product(delta, f[event_node] + u) -
    dot_product(expo * exp(f), exp(u));
  target += -0.5 * n * log_determinant(sigma_b) -
    0.5 * sum(columns_dot_product(b', inverse_spd(sigma_b) * b'));
  beta ~ normal(0, 10);
  target += gamma_lpdf(exp(-2 * log_sigma_e) | 1, 5e-5) + log(2) -
    2 * log_sigma_e;
  sigma_b ~ inv_wishart(3, identity_matrix(2));
  gam ~ std_normal();
  eta ~ std_normal();
  target += log(pc_lambda / 2) - 0.5 * log_tau -
    pc_lambda * exp(-0.5 * log_tau);
}
generated quantities {
  real sigma_e = exp(log_sigma_e);
  real sd_b0 = sqrt(sigma_b[1, 1]);
  real sd_b1 = sqrt(sigma_b[2, 2]);
  real rho = sigma_b[1, 2] / (sd_b0 * sd_b1);
  real var_h0 = exp(-log_tau);
}
"
stan_mod <- rstan::stan_model(model_code = stan_code)
t_stan <- system.time(
  stan_fit <- rstan::sampling(
    stan_mod,
    data = list(
      N = nrow(dat_long),
      n = n,
      K = n_cut,
      id = dat_long$id,
      time = dat_long$time,
      y = dat_long$Y,
      trt = trt,
      delta = event_indic,
      event_node = event_node,
      expo = expo,
      rw_u = rw_u,
      rw_lam = rw_lam,
      rw_diag = rw_diag,
      pc_lambda = pc_lambda
    ),
    pars = c(
      "beta",
      "log_sigma_e",
      "sigma_b",
      "gam",
      "eta",
      "log_tau",
      "f",
      "sigma_e",
      "sd_b0",
      "sd_b1",
      "rho",
      "var_h0"
    ),
    chains = 4,
    iter = n_iter,
    seed = 2026,
    refresh = 0,
    control = list(adapt_delta = 0.95)
  )
)[["elapsed"]]
sd_stan <- rstan::extract(stan_fit)
draws_stan <- cbind(
  beta0 = sd_stan$beta[, 1],
  beta1 = sd_stan$beta[, 2],
  beta2 = sd_stan$beta[, 3],
  beta3 = sd_stan$beta[, 4],
  sigma_e = sd_stan$sigma_e,
  sd_b0 = sd_stan$sd_b0,
  sd_b1 = sd_stan$sd_b1,
  rho = sd_stan$rho,
  var_h0 = sd_stan$var_h0,
  gamma0 = sd_stan$gam[, 1],
  gamma1 = sd_stan$gam[, 2],
  setNames(
    as.data.frame(sd_stan$f[, -n_cut]),
    sprintf("f%02d", seq_len(n_cut - 1))
  )
)
draws_stan <- as.matrix(draws_stan)

## ----- Log marginal likelihood by importance sampling ------------------------

# A multivariate t (5 df) with the Stan draws' mean and covariance on the scale
# of version A's 27 parameters, Pareto-smoothed. A draw where logpost is not
# finite (a correlation of 1, say) gets weight 0. INLA's coxph is a Poisson
# likelihood on the expanded data, which adds sum(log E) over the event rows,
# the exposure of each event's interval, to its log marginal likelihood.
theta_stan <- cbind(
  sd_stan$beta,
  sd_stan$log_sigma_e,
  log(sd_stan$sd_b0),
  log(sd_stan$sd_b1),
  atanh(sd_stan$rho),
  sd_stan$gam,
  sd_stan$eta,
  sd_stan$log_tau
)
colnames(theta_stan) <- names(init_a)
is_df <- 5
is_mean <- colMeans(theta_stan)
is_chol <- chol(cov(theta_stan) * (is_df - 2) / is_df)
set.seed(3)
is_z <- matrix(rnorm(n_is * length(init_a)), n_is) /
  sqrt(rchisq(n_is, is_df) / is_df)
is_theta <- sweep(is_z %*% is_chol, 2, is_mean, "+")
colnames(is_theta) <- names(init_a)
is_logq <- mvtnorm::dmvt(
  is_theta,
  delta = is_mean,
  sigma = crossprod(is_chol),
  df = is_df,
  log = TRUE
)
is_lp <- suppressWarnings(apply(is_theta, 1, logpost_a))
is_lp[!is.finite(is_lp)] <- -Inf
is_lw <- is_lp - is_logq
is_psis <- suppressWarnings(loo::psis(is_lw, r_eff = NA))
logml_ref <- max(is_lw) + log(mean(exp(is_lw - max(is_lw))))
poisson_const <- sum(log(expo[cbind(1:n, event_node)])[event_indic == 1])

# R-INLA normalises the rw1 with the constant of an intrinsic model, even with
# diagonal > 0, which makes it proper. The offset comes from a normal model
# with the same rw1 and every hyperparameter fixed, whose log marginal
# likelihood is a normal density. It barely moves with tau (-4.78 at tau = 5,
# -4.76 at 5,000), and the longitudinal submodel on its own needs no offset.
rw_test_y <- rnorm(n_cut, log_rate, 1)
rw_test <- inla(
  y ~ -1 +
    f(
      idx,
      model = "rw1",
      constr = FALSE,
      diagonal = rw_diag,
      scale.model = TRUE,
      hyper = list(prec = list(initial = 6, fixed = TRUE))
    ),
  data = list(y = rw_test_y, idx = seq_len(n_cut)),
  control.family = list(hyper = list(prec = list(initial = 0, fixed = TRUE)))
)
rw_offset <- rw_test$mlik[1, 1] -
  mvtnorm::dmvnorm(
    rw_test_y,
    rep(0, n_cut),
    solve(exp(6) * as.matrix(rw_r) + rw_diag * diag(n_cut)) + diag(n_cut),
    log = TRUE
  )

## =============================================================================
## Comparison
## =============================================================================

## ----- Posterior summaries ---------------------------------------------------

par_main <- colnames(draws_stan)[1:11]
par_f <- colnames(draws_stan)[-(1:11)]
summarise <- function(draws) {
  cbind(
    mean = colMeans(draws),
    sd = apply(draws, 2, sd),
    t(apply(draws, 2, quantile, c(0.025, 0.975)))
  )
}
ref <- summarise(draws_stan)
s_a <- summarise(draws_a)
s_b <- summarise(draws_b)

# INLAjoint's summaries, in the same layout. The SDs and the correlation of
# the random effects, and the baseline variance, come from its summary(), the
# log-hazards from INLA's latent-field marginals. Each table has the columns
# mean, sd and the 2.5%, 50% and 97.5% quantiles.
ij_cols <- function(x) as.matrix(x)[, c(1, 2, 3, 5), drop = FALSE]
s_ij <- rbind(
  ij_cols(sum_ij$FixedEff[[1]]),
  ij_cols(sum_ij$ReffList[[1]]),
  ij_cols(sum_ij_var$SurvEff[[1]]),
  ij_cols(sum_ij$AssocLS),
  ij_cols(fit_ij$summary.random$baseline1.hazard[-n_cut, -1])
)
dimnames(s_ij) <- dimnames(ref)

gap <- function(s) {
  cbind(
    mean = (s[, 1] - ref[, 1]) / ref[, 2],
    sd_ratio = s[, 2] / ref[, 2],
    q025 = (s[, 3] - ref[, 3]) / ref[, 2],
    q975 = (s[, 4] - ref[, 4]) / ref[, 2]
  )
}
worst <- function(g) {
  c(
    max_mean_gap = max(abs(g[par_f, "mean"])),
    sd_ratio_min = min(g[par_f, "sd_ratio"]),
    sd_ratio_max = max(g[par_f, "sd_ratio"]),
    max_q_gap = max(abs(g[par_f, c("q025", "q975")]))
  )
}
cat("\n== Stan reference ==\n")
print(signif(ref[par_main, ], 3))
cat("\n== INLAjoint against Stan, in Stan's posterior SDs ==\n")
print(round(gap(s_ij)[par_main, ], 3))
cat(
  "\n== inla_vanilla(), version B, against Stan, in Stan's posterior SDs ==\n"
)
print(round(gap(s_b)[par_main, ], 3))
cat(
  "\n== inla_vanilla(), version A, against Stan, in Stan's posterior SDs ==\n"
)
print(round(gap(s_a)[par_main, ], 3))
cat("\n== Baseline log-hazards f01 to f15: worst over the 15 intervals ==\n")
print(round(
  rbind(
    INLAjoint = worst(gap(s_ij)),
    version_b = worst(gap(s_b)),
    version_a = worst(gap(s_a))
  ),
  3
))

# The smoothing precision on its own scale, where INLA reports it
cat("\n== log tau: mean and quantiles ==\n")
log_tau_ij <- fit_ij$internal.marginals.hyperpar[[
  "Log precision for baseline1.hazard"
]]
print(round(
  rbind(
    stan = c(
      mean(sd_stan$log_tau),
      quantile(sd_stan$log_tau, c(0.025, 0.5, 0.975))
    ),
    INLAjoint = c(
      inla.emarginal(identity, log_tau_ij),
      inla.qmarginal(c(0.025, 0.5, 0.975), log_tau_ij)
    ),
    version_b = fit_b$summary["log_tau", c("Mean", "2.5%", "50%", "97.5%")]
  ),
  2
))

## ----- Log marginal likelihood -----------------------------------------------

cat("\n== Log marginal likelihood ==\n")
cat(
  "INLAjoint's less sum(log E) =",
  round(poisson_const, 2),
  "and its rw1 offset",
  round(rw_offset, 2),
  "\n"
)
print(round(
  c(
    importance = logml_ref,
    version_b = fit_b$logml,
    version_a = fit_a$logml,
    inlajoint_integration = fit_ij$mlik[1, 1] - poisson_const - rw_offset,
    inlajoint_gaussian = fit_ij$mlik[2, 1] - poisson_const - rw_offset
  ),
  2
))
cat(
  "Importance sampling: Pareto k",
  round(is_psis$diagnostics$pareto_k, 2),
  ", ESS",
  round(is_psis$diagnostics$n_eff),
  "of",
  n_is,
  "\n"
)

## ----- Cost and checks -------------------------------------------------------

cat("\n== Cost ==\n")
cat(
  "INLAjoint:",
  round(t_ij, 2),
  "s with",
  inla.getOption("num.threads"),
  "threads,",
  round(t_ij_1, 2),
  "s on one\n"
)
cost <- function(label, fit, t_fit, t_draws) {
  cat(
    label,
    fit$counts[["logpost"]],
    "calls,",
    round(t_fit, 2),
    "s, then",
    round(t_draws, 2),
    "s for",
    n_draws,
    "draws\n"
  )
}
cost("inla_vanilla(), version B:", fit_b, t_b, t_b_draws)
cost("inla_vanilla(), version A:", fit_a, t_a, t_a_draws)
cat("Stan:", round(t_stan, 1), "s of sampling\n")
cat(
  "Stan: max R-hat",
  round(max(rstan::summary(stan_fit)$summary[, "Rhat"], na.rm = TRUE), 4)
)
cat(", divergences", rstan::get_num_divergent(stan_fit), "\n")
cat(
  "Quadrature: logpost at the mode, 41 nodes minus 15:",
  signif(quad_error, 2),
  "\n\n"
)
summary(fit_b)
