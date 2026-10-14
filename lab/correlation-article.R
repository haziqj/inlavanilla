################################################################################
#
# Precomputes the Stan reference for the article "Correlation matrices"
# (vignettes/articles/correlation.qmd), so that the site build needs no Stan.
# Fits a multivariate normal to the nine Holzinger and Swineford (1939) tests
# in lavaan, with N(0, 10^2) priors on the means, Exp(1) priors on the SDs and
# an LKJ(2) prior on the correlation matrix, the same model and priors as the
# article's logpost. Saves each parameter's posterior summary, on the natural
# scale, kernel density estimates of the 36 correlations and a summary of
# Cronbach's alpha for each three-test scale, computed in Stan, to
# vignettes/articles/correlation-stan.rds. Takes under a minute, most of it
# compiling Stan. Run from the package root:
#   Rscript lab/correlation-article.R
#
################################################################################

## ----- Configuration ---------------------------------------------------------

out_file <- "vignettes/articles/correlation-stan.rds"
probs <- c(0.025, 0.5, 0.975)
eta <- 2

data(HolzingerSwineford1939, package = "lavaan")
x <- as.matrix(HolzingerSwineford1939[paste0("x", 1:9)])
d <- ncol(x)
low <- lower.tri(diag(d))
pair_row <- row(diag(d))[low]
pair_col <- col(diag(d))[low]

lkj_stan <- "
data {
  int<lower=1> N;
  int<lower=1> D;
  array[N] vector[D] x;
  real<lower=0> eta;
}
parameters {
  vector[D] mu;
  vector<lower=0>[D] sigma;
  cholesky_factor_corr[D] L;
}
model {
  mu ~ normal(0, 10);
  sigma ~ exponential(1);
  L ~ lkj_corr_cholesky(eta);
  x ~ multi_normal_cholesky(mu, diag_pre_multiply(sigma, L));
}
generated quantities {
  corr_matrix[D] Omega = multiply_lower_tri_self_transpose(L);
  vector[3] alpha;
  for (s in 1:3) {
    matrix[3, 3] S = quad_form_diag(
      Omega[(3 * s - 2):(3 * s), (3 * s - 2):(3 * s)],
      sigma[(3 * s - 2):(3 * s)]
    );
    alpha[s] = 1.5 * (1 - trace(S) / sum(S));
  }
}
"

## ----- Fit and summarise ------------------------------------------------------

fit <- rstan::stan(
  model_code = lkj_stan,
  data = list(N = nrow(x), D = d, x = x, eta = eta),
  chains = 4,
  iter = 6000,
  warmup = 1000,
  seed = 1,
  refresh = 0
)

stan_pars <- c(
  sprintf("mu[%d]", seq_len(d)),
  sprintf("sigma[%d]", seq_len(d)),
  sprintf("Omega[%d,%d]", pair_row, pair_col)
)
draws <- as.matrix(fit, pars = c("mu", "sigma", "Omega"))[, stan_pars]
colnames(draws) <- c(
  paste0("mu_x", seq_len(d)),
  paste0("sigma_x", seq_len(d)),
  paste0("r", pair_col, pair_row)
)
alpha <- as.matrix(fit, pars = "alpha")
colnames(alpha) <- c("visual", "textual", "speed")
r_names <- grep("^r", colnames(draws), value = TRUE)
elapsed <- rstan::get_elapsed_time(fit)

summarise <- function(draws) {
  data.frame(
    mean = colMeans(draws),
    sd = apply(draws, 2, sd),
    t(apply(draws, 2, quantile, probs)),
    check.names = FALSE
  )
}
ref <- list(
  summary = summarise(draws),
  density = lapply(
    as.data.frame(draws[, r_names]),
    function(r) lapply(density(r, n = 128)[c("x", "y")], signif, 4)
  ),
  alpha = summarise(alpha),
  draws = nrow(draws),
  seconds = max(rowSums(elapsed)),
  rhat = max(rstan::summary(fit)$summary[, "Rhat"], na.rm = TRUE),
  divergences = sum(rstan::get_divergent_iterations(fit)),
  rstan_version = as.character(packageVersion("rstan"))
)
saveRDS(ref, out_file, compress = "xz")
str(ref, max.level = 1)
