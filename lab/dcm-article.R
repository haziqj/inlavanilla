################################################################################
#
# Precomputes the Stan reference for the article "Diagnostic classification
# models" (vignettes/articles/dcm.qmd), so that the site build needs no Stan.
# Fits the DINA model to the ECPE data (edmdata) with the latent classes summed
# out, the same model and priors as the article's logpost: uniform guessing and
# slipping probabilities, and a Dirichlet(1) prior on the eight class
# proportions. Saves each parameter's posterior summary, on the natural scale,
# and kernel density estimates of the class proportions to
# vignettes/articles/dcm-stan.rds. Takes a few minutes. Run from the package
# root:
#   Rscript lab/dcm-article.R
#
################################################################################

## ----- Configuration ----------------------------------------------------------

out_file <- "vignettes/articles/dcm-stan.rds"
probs <- c(0.025, 0.5, 0.975)

data(items_ecpe, qmatrix_ecpe, package = "edmdata")
y <- unclass(items_ecpe)
q <- unclass(qmatrix_ecpe)
alpha <- as.matrix(expand.grid(rep(list(0:1), ncol(q))))
eta <- (alpha %*% t(q)) ==
  matrix(rowSums(q), nrow(alpha), nrow(q), byrow = TRUE)
class_names <- apply(alpha, 1, paste, collapse = "")

dina_stan <- "
data {
  int<lower=1> N;
  int<lower=1> J;
  int<lower=1> C;
  matrix[N, J] y;
  matrix[C, J] eta;
}
parameters {
  vector<lower=0, upper=1>[J] g;
  vector<lower=0, upper=1>[J] s;
  simplex[C] pi;
}
model {
  matrix[J, C] p = eta' .* rep_matrix(1 - s, C)
    + (1 - eta') .* rep_matrix(g, C);
  row_vector[C] l0 = log(pi)' + rep_row_vector(1, J) * log1m(p);
  matrix[N, C] l = y * logit(p) + rep_matrix(l0, N);
  pi ~ dirichlet(rep_vector(1, C));
  for (i in 1:N) {
    target += log_sum_exp(l[i]);
  }
}
"

## ----- Fit and summarise ------------------------------------------------------

model <- rstan::stan_model(model_code = dina_stan)
init_one <- list(
  g = rep(0.2, ncol(y)),
  s = rep(0.1, ncol(y)),
  pi = rep(1 / nrow(alpha), nrow(alpha))
)
fit <- rstan::sampling(
  model,
  data = list(
    N = nrow(y),
    J = ncol(y),
    C = nrow(alpha),
    y = y,
    eta = eta * 1
  ),
  chains = 4,
  cores = 4,
  iter = 3500,
  warmup = 1000,
  seed = 1,
  refresh = 0,
  init = rep(list(init_one), 4)
)

draws <- as.matrix(fit, pars = c("g", "s", "pi"))
colnames(draws) <- c(
  sprintf("g%02d", seq_len(ncol(y))),
  sprintf("s%02d", seq_len(ncol(y))),
  paste0("pi", class_names)
)
pi_draws <- draws[, paste0("pi", class_names)]
elapsed <- rstan::get_elapsed_time(fit)

ref <- list(
  summary = data.frame(
    mean = colMeans(draws),
    sd = apply(draws, 2, sd),
    t(apply(draws, 2, quantile, probs)),
    check.names = FALSE
  ),
  density = lapply(
    as.data.frame(pi_draws),
    function(d) density(d, n = 256, from = 0)[c("x", "y")]
  ),
  draws = nrow(draws),
  seconds = max(rowSums(elapsed)),
  rhat = max(rstan::summary(fit)$summary[, "Rhat"], na.rm = TRUE),
  divergences = sum(rstan::get_divergent_iterations(fit)),
  rstan_version = as.character(packageVersion("rstan"))
)
saveRDS(ref, out_file, compress = "xz")
str(ref, max.level = 1)
