################################################################################
#
# A nonlinear mixed model with the random effects integrated out by RTMB's
# Laplace approximation (October 2026): can TMB + inla_vanilla() give a fast
# Bayesian fit of a model R-INLA cannot fit? The theophylline data (Theoph, in
# base R) hold the plasma concentrations of 12 subjects over 25 hours after an
# oral dose. The one-compartment model with first-order absorption of
# Pinheiro and Bates (2000, ch. 8) has random effects on the log absorption
# rate and the log clearance, which enter the mean nonlinearly.
#
# The population parameters theta are the three fixed effects (log ke, log
# ka, log Cl) and the log SDs of the two random effects and of the residual,
# with N(0, 10^2) priors on the fixed effects and Exponential(1) priors on the
# SDs. Four fits of the same posterior:
#   1. RTMB integrates out the 24 subject effects by Laplace, and
#      inla_vanilla() fits the 6 population parameters (the TMB route).
#   2. inla_vanilla() fits all 30 parameters, with RTMB's gradient.
#   3. Stan (NUTS, non-centred) as the reference.
#   4. nlme's ML fit, as a frequentist anchor.
# The subject effects of fit 1 come from sampling(): at each draw of theta,
# the inner Laplace mode and Hessian give a Gaussian draw of the effects.
#
# Needs RTMB and rstan. Run from the package root:
#   Rscript lab/theoph-rtmb.R
#
# Results (October 2026; Stan 4 chains of 3,000 kept draws, R-hat 1.005, no
# divergences): the TMB route works. It took 0.4 s (252 objective and 365
# gradient calls, each an inner Laplace), against 6 s of sampling for Stan,
# after a minute of compiling. The population means are within 0.05 posterior
# SDs of Stan's, and the SDs of the three log SDs within 4%. The fixed effects
# log ka and log Cl are 13% too narrow (97.5% quantiles 0.35 and 0.24 SDs
# low). summary() flags the ridge bend of log ka (0.38) but not that of log Cl
# (0.19). The subject effects, drawn from the inner Gaussian at each copula
# draw of theta (4 s for 10,000 draws), are within 3% of Stan's SDs for 22 of
# the 24. Log ka of subject 9, whose absorption is barely identified (Stan
# skewness 1.8), is 27% too narrow and 0.3 SDs off in mean, beyond what a
# Gaussian can follow. Putting all 30 parameters in theta is worse: every
# subject effect is about 12% too narrow, and the fit warns of a VB shift of
# 1.4 SDs in log sigma. The joint mode of a hierarchical model has its SDs
# too small (VB moves the random-effect log SDs up by 0.5 SDs), and the scans
# start there. Integrating the effects out at each theta, as INLA does,
# avoids that. The two log marginal likelihoods agree to 0.19 (-199.24 and
# -199.43), and the inner optimisation added no noise: obj$fn at the mode was
# the same before and after a call elsewhere.
#
################################################################################

## ----- Configuration ---------------------------------------------------------

smoke <- FALSE
suppressMessages(devtools::load_all(".", quiet = TRUE))
library(nlme)
n_iter <- if (smoke) 1000 else 6000 # Stan iterations per chain, half warm-up
n_draws <- if (smoke) 1000 else 10000 # copula draws
options(mc.cores = 4)

## =============================================================================
## Model and data
## =============================================================================

id <- as.integer(as.character(Theoph$Subject))
n_sub <- max(id)
time <- Theoph$Time
dose <- Theoph$Dose
conc <- Theoph$conc

pop_names <- c("b_ke", "b_ka", "b_cl", "log_sd_ka", "log_sd_cl", "log_sigma")

## ----- Log joint density -----------------------------------------------------

# A function of the parameter list, for RTMB. The concentration is SSfol()'s
# Dose ke ka / Cl (exp(-ke t) - exp(-ka t)) / (ka - ke). An Exponential(1)
# prior on an SD s = exp(l) gives l - exp(l) on the log scale, with the
# Jacobian.
log_joint <- function(p) {
  ke <- exp(p$b_ke)
  ka <- exp(p$b_ka + p$u_ka)[id]
  cl <- exp(p$b_cl + p$u_cl)[id]
  mu <- dose * ke * ka / cl * (exp(-ke * time) - exp(-ka * time)) / (ka - ke)
  sum(RTMB::dnorm(conc, mu, exp(p$log_sigma), log = TRUE)) +
    sum(RTMB::dnorm(p$u_ka, 0, exp(p$log_sd_ka), log = TRUE)) +
    sum(RTMB::dnorm(p$u_cl, 0, exp(p$log_sd_cl), log = TRUE)) +
    RTMB::dnorm(p$b_ke, 0, 10, log = TRUE) +
    RTMB::dnorm(p$b_ka, 0, 10, log = TRUE) +
    RTMB::dnorm(p$b_cl, 0, 10, log = TRUE) +
    p$log_sd_ka -
    exp(p$log_sd_ka) +
    p$log_sd_cl -
    exp(p$log_sd_cl) +
    p$log_sigma -
    exp(p$log_sigma)
}

par0 <- list(
  b_ke = -2.5,
  b_ka = 0.5,
  b_cl = -3,
  log_sd_ka = log(0.5),
  log_sd_cl = log(0.3),
  log_sigma = log(0.7),
  u_ka = numeric(n_sub),
  u_cl = numeric(n_sub)
)
init_pop <- unlist(par0[pop_names])
u_names <- c(paste0("u_ka[", 1:n_sub, "]"), paste0("u_cl[", 1:n_sub, "]"))
init_full <- c(init_pop, setNames(numeric(2 * n_sub), u_names))

## =============================================================================
## Fits
## =============================================================================

## ----- 1. RTMB Laplace over the subject effects ------------------------------

# TMB minimises over the random effects, so the objective is the negative log
# joint, and obj$fn the negative Laplace-approximate log posterior of theta.
obj_lap <- RTMB::MakeADFun(
  function(p) -log_joint(p),
  par0,
  random = c("u_ka", "u_cl"),
  silent = TRUE
)
t_lap <- system.time(
  fit_lap <- inla_vanilla(
    function(theta) -obj_lap$fn(theta),
    init_pop,
    grad = function(theta) -as.vector(obj_lap$gr(theta))
  )
)[["elapsed"]]

# The inner optimisation starts from the last inner mode, so a value can
# depend on the calls before it
f_before <- obj_lap$fn(fit_lap$mode)
invisible(obj_lap$fn(fit_lap$mode + 2))
inner_noise <- obj_lap$fn(fit_lap$mode) - f_before

## ----- 2. All 30 parameters --------------------------------------------------

obj_full <- RTMB::MakeADFun(log_joint, par0, silent = TRUE)
t_full <- system.time(
  fit_full <- inla_vanilla(obj_full$fn, init_full, grad = obj_full$gr)
)[["elapsed"]]

## ----- 3. Stan ---------------------------------------------------------------

stan_code <- "
data {
  int<lower=1> N;
  int<lower=1> J;
  array[N] int<lower=1, upper=J> id;
  vector[N] time;
  vector[N] dose;
  vector[N] conc;
}
parameters {
  real b_ke;
  real b_ka;
  real b_cl;
  real<lower=0> sd_ka;
  real<lower=0> sd_cl;
  real<lower=0> sigma;
  vector[J] z_ka;
  vector[J] z_cl;
}
transformed parameters {
  vector[J] u_ka = sd_ka * z_ka;
  vector[J] u_cl = sd_cl * z_cl;
}
model {
  real ke = exp(b_ke);
  vector[N] ka = exp(b_ka + u_ka[id]);
  vector[N] cl = exp(b_cl + u_cl[id]);
  vector[N] mu = ke * dose .* ka ./ cl .* (exp(-ke * time) - exp(-ka .* time))
    ./ (ka - ke);
  conc ~ normal(mu, sigma);
  b_ke ~ normal(0, 10);
  b_ka ~ normal(0, 10);
  b_cl ~ normal(0, 10);
  sd_ka ~ exponential(1);
  sd_cl ~ exponential(1);
  sigma ~ exponential(1);
  z_ka ~ std_normal();
  z_cl ~ std_normal();
}
generated quantities {
  real log_sd_ka = log(sd_ka);
  real log_sd_cl = log(sd_cl);
  real log_sigma = log(sigma);
}
"
stan_mod <- rstan::stan_model(model_code = stan_code)
t_stan <- system.time(
  stan_fit <- rstan::sampling(
    stan_mod,
    data = list(
      N = length(conc),
      J = n_sub,
      id = id,
      time = time,
      dose = dose,
      conc = conc
    ),
    chains = 4,
    iter = n_iter,
    seed = 2026,
    refresh = 0,
    control = list(adapt_delta = 0.95)
  )
)[["elapsed"]]
stan_draws <- as.matrix(stan_fit)

## ----- 4. nlme ---------------------------------------------------------------

fit_nlme <- nlme(
  conc ~ SSfol(Dose, Time, lKe, lKa, lCl),
  data = Theoph,
  fixed = lKe + lKa + lCl ~ 1,
  random = pdDiag(lKa + lCl ~ 1),
  start = c(lKe = -2.5, lKa = 0.5, lCl = -3),
  method = "ML"
)

## =============================================================================
## Comparison
## =============================================================================

## ----- Population parameters -------------------------------------------------

# Each fit against Stan, in Stan's posterior SDs: the mean and the 2.5% and
# 97.5% quantiles, and the SD as a ratio
probs <- c(0.025, 0.975)
ref <- cbind(
  mean = colMeans(stan_draws[, pop_names]),
  sd = apply(stan_draws[, pop_names], 2, sd),
  t(apply(stan_draws[, pop_names], 2, quantile, probs))
)
against_stan <- function(s) {
  cbind(
    mean = (s[, "Mean"] - ref[, "mean"]) / ref[, "sd"],
    sd_ratio = s[, "SD"] / ref[, "sd"],
    q025 = (s[, "2.5%"] - ref[, "2.5%"]) / ref[, "sd"],
    q975 = (s[, "97.5%"] - ref[, "97.5%"]) / ref[, "sd"]
  )
}

cat("\n== Stan reference ==\n")
print(signif(ref, 3))
cat("\n== RTMB Laplace + inla_vanilla(), against Stan ==\n")
print(round(against_stan(fit_lap$summary[pop_names, ]), 3))
cat("\n== inla_vanilla() on all 30 parameters, against Stan ==\n")
print(round(against_stan(fit_full$summary[pop_names, ]), 3))

cat("\n== nlme (ML) ==\n")
print(round(fixef(fit_nlme), 3))
print(round(log(as.numeric(VarCorr(fit_nlme)[, "StdDev"])), 3))

## ----- Subject effects -------------------------------------------------------

# Draws of theta from the copula, then of the subject effects from the inner
# Laplace Gaussian at that theta: N(mode, H^-1), with H the inner Hessian of
# the negative log joint
draw_subjects <- function(theta) {
  obj_lap$fn(theta)
  par <- obj_lap$env$last.par
  r <- obj_lap$env$random
  h <- as.matrix(obj_lap$env$spHess(par, random = TRUE))
  u <- par[r] + backsolve(chol(h), rnorm(length(r)))
  c(theta[["b_ka"]] + u[1:n_sub], theta[["b_cl"]] + u[n_sub + 1:n_sub])
}
subject_names <- c(paste0("lka[", 1:n_sub, "]"), paste0("lcl[", 1:n_sub, "]"))
set.seed(1)
t_draws <- system.time(
  sub_lap <- sampling(fit_lap, nsamp = n_draws, fun = draw_subjects)
)[["elapsed"]]
sub_full <- sampling(fit_full, nsamp = n_draws, fun = function(theta) {
  c(
    theta[["b_ka"]] + theta[paste0("u_ka[", 1:n_sub, "]")],
    theta[["b_cl"]] + theta[paste0("u_cl[", 1:n_sub, "]")]
  )
})
sub_stan <- cbind(
  stan_draws[, "b_ka"] + stan_draws[, paste0("u_ka[", 1:n_sub, "]")],
  stan_draws[, "b_cl"] + stan_draws[, paste0("u_cl[", 1:n_sub, "]")]
)
colnames(sub_lap) <- colnames(sub_full) <- colnames(sub_stan) <- subject_names

sub_error <- function(draws) {
  sd_ref <- apply(sub_stan, 2, sd)
  c(
    mean = max(abs(colMeans(draws) - colMeans(sub_stan)) / sd_ref),
    sd_ratio_min = min(apply(draws, 2, sd) / sd_ref),
    sd_ratio_max = max(apply(draws, 2, sd) / sd_ref)
  )
}
cat("\n== Subject log ka and log Cl, worst over the 24, against Stan ==\n")
print(round(rbind(laplace = sub_error(sub_lap), full = sub_error(sub_full)), 3))

sd_ratio_by_subject <- function(draws) {
  ratio <- apply(draws, 2, sd) / apply(sub_stan, 2, sd)
  matrix(ratio, 2, byrow = TRUE, dimnames = list(c("lka", "lcl"), 1:n_sub))
}
cat("\n== SD ratio to Stan by subject, Laplace route ==\n")
print(round(sd_ratio_by_subject(sub_lap), 2))
cat("\n== SD ratio to Stan by subject, all 30 parameters ==\n")
print(round(sd_ratio_by_subject(sub_full), 2))
cat("\n== Stan skewness of log ka by subject ==\n")
skew <- function(x) mean((x - mean(x))^3) / sd(x)^3
print(round(apply(sub_stan[, 1:n_sub], 2, skew), 2))

## ----- Cost and checks -------------------------------------------------------

cat("\n== Cost ==\n")
print(rbind(
  laplace = c(fit_lap$counts, seconds = round(t_lap, 2)),
  full = c(fit_full$counts, seconds = round(t_full, 2)),
  stan = c(logpost = NA, grad = NA, seconds = round(t_stan, 2))
))
cat("Subject draws for the Laplace fit:", round(t_draws, 2), "s\n")
cat("Inner warm-start noise in obj$fn:", signif(inner_noise, 2), "\n")
cat(
  "Log marginal likelihood: Laplace route",
  round(fit_lap$logml, 3),
  "| full",
  round(fit_full$logml, 3),
  "\n"
)
cat(
  "\nStan: max R-hat",
  round(max(rstan::summary(stan_fit)$summary[, "Rhat"]), 4)
)
cat(", divergences", rstan::get_num_divergent(stan_fit), "\n")
cat("\nWarnings, Laplace route:", length(fit_lap$warnings), "| full:")
cat(length(fit_full$warnings), "\n")
print(fit_lap)
