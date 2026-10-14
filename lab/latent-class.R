################################################################################
#
# Example: a latent class model of political participation in ten European
# countries. The data are ESS_datLCA from the latent package (Jimenez,
# Garnier-Villarreal & Franco, github.com/Marcosjnez/latent), eight yes/no
# items from the European Social Survey: in the last 12 months, did the
# respondent contact a politician, wear a campaign badge, sign a petition,
# boycott a product, donate to a party, take part in a demonstration, volunteer
# for an organisation, or post about politics online? Three latent classes
# share their item profiles across the countries, and each country has its own
# class sizes: a multigroup latent class model with invariant measurement, the
# same model as the one-step latent class regression on country that
# latent::lca() fits. logpost() sums the classes out, and sampling() gives the
# class profiles, the class sizes in each country and the class membership of a
# response pattern. If latent is installed, its estimates are a check on the
# fit, and if rstan is, Stan checks the marginals. Step through it in an R
# session.
#
# Results (October 2026): the fit took 11 s, and all ten random starts found
# the same mode. Against Stan (4 chains of 2000 draws, 8 minutes, R-hat 1.002),
# every mean is within 0.1 posterior SDs and every SD within 5%. summary()
# flags only BG_civic (ridge bend 0.59, mean gap 0.12), Bulgaria's nearly empty
# civic class, and Stan agrees there too. latent's estimates are within 0.2
# SDs of the posterior means except in Bulgaria and Hungary, where the civic
# class is nearly empty and the priors differ (BG_civic 0.4% against 1.0%).
# The activists are 2% of Hungarians and 14% of the French and the Finns. The
# civic class is 69% of Finns and 53% of the French but 1% of Bulgarians, so a
# respondent who only signed a petition is civic with probability 0.90 in
# Finland and inactive with probability 0.92 in Bulgaria.
#
################################################################################

## ----- Configuration ---------------------------------------------------------

library(inlavanilla)
set.seed(1)
n_class <- 3
n_starts <- 10 # random starts for the mode search
nsamp <- 4000 # draws from the fitted posterior
run_latent <- requireNamespace("latent", quietly = TRUE)
run_stan <- requireNamespace("rstan", quietly = TRUE)

items <- c(
  "contplt", # contacted a politician
  "badge", # wore a campaign badge
  "sgnptit", # signed a petition
  "bctprd", # boycotted a product
  "donprty", # donated to a party
  "pbldmna", # took part in a demonstration
  "volunfp", # volunteered for an organisation
  "pstplonl" # posted about politics online
)

# latent is on GitHub only, and needs a C++ toolchain with OpenMP to install.
# Without it, the data come from the package's repository, at a fixed commit.
ess_url <- paste0(
  "https://raw.githubusercontent.com/Marcosjnez/latent/",
  "38d2479265735acea91235c279342f7640b3623d/data/ESS_datLCA.rda"
)

## =============================================================================
## Data
## =============================================================================

if (run_latent) {
  ess <- latent::ESS_datLCA
} else {
  rda <- tempfile(fileext = ".rda")
  download.file(ess_url, rda, mode = "wb", quiet = TRUE)
  load(rda) # creates ESS_datLCA
  ess <- ESS_datLCA
}

# Items are coded 1 = yes and 2 = no. The 6% of respondents with a missing
# item are dropped, as in latent's workshop on these data.
ess <- ess[complete.cases(ess[items]), c(items, "cntry")]
ess[items] <- lapply(ess[items], function(v) as.integer(v == 1))

# The likelihood depends on the data only through the counts of each response
# pattern in each country: 1045 rows instead of 17057.
pat <- aggregate(list(count = rep(1L, nrow(ess))), ess, sum)
y <- as.matrix(pat[items])
count <- pat$count
country <- factor(pat$cntry)
n_country <- nlevels(country)

# The share answering yes to each item, by country
round(sapply(split(ess[items], ess$cntry), colMeans), 3)

## =============================================================================
## Model
## =============================================================================

## ----- Parameters ------------------------------------------------------------

# theta holds the class logits of each country, log P(class k) - log P(class
# 1), for k = 2, ..., K (class 1 is the reference), then the item logits
# log P(yes | class k) - log P(no | class k), item by item within each class.
unpack <- function(theta, n_country, n_class, n_item) {
  n_beta <- n_country * (n_class - 1)
  list(
    beta = cbind(0, matrix(theta[seq_len(n_beta)], n_country)),
    lambda = matrix(theta[-seq_len(n_beta)], n_item, n_class)
  )
}

row_logsumexp <- function(a) {
  a_max <- apply(a, 1, max)
  a_max + log(rowSums(exp(a - a_max)))
}

# The pieces that logpost() and its gradient share. Given its class, a
# respondent answers the items independently, so log P(y | class k) is
# sum_j y_j lambda_jk + log(1 - p_jk), and a[i, k] is the log of the joint
# probability of pattern i and class k in the pattern's country.
lca_terms <- function(theta, y, country, n_class) {
  p <- unpack(theta, nlevels(country), n_class, ncol(y))
  log_class <- p$beta - row_logsumexp(p$beta) # countries x classes
  log_no <- colSums(plogis(-p$lambda, log.p = TRUE))
  a <- log_class[as.integer(country), , drop = FALSE] +
    y %*% p$lambda +
    rep(log_no, each = nrow(y))
  c(p, list(log_class = log_class, a = a, lse = row_logsumexp(a)))
}

## ----- Log-posterior and gradient --------------------------------------------

# Priors: a flat Dirichlet(1) on the class sizes of each country, and a uniform
# prior on each P(yes | class). On the logit scale, with the Jacobians, these
# are prod_k P(class k) times (K - 1)! for each country, and a standard
# logistic density for each item logit, so fit$logml is a proper log marginal
# likelihood.
logpost <- function(theta, y, count, country, n_class) {
  tm <- lca_terms(theta, y, country, n_class)
  sum(count * tm$lse) +
    sum(tm$log_class) +
    nlevels(country) * lfactorial(n_class - 1) +
    sum(dlogis(tm$lambda, log = TRUE))
}

# The gradient needs only the posterior class probabilities of each pattern,
# as in the E step of EM: w[i, k] is the expected number of respondents with
# pattern i in class k.
logpost_grad <- function(theta, y, count, country, n_class) {
  tm <- lca_terms(theta, y, country, n_class)
  w <- count * exp(tm$a - tm$lse)
  n_ck <- rowsum(w, country) # expected class counts in each country
  p_class <- exp(tm$log_class)
  p_yes <- plogis(tm$lambda)
  grad_beta <- n_ck - rowSums(n_ck) * p_class
  grad_beta <- grad_beta + 1 - n_class * p_class # Dirichlet prior
  grad_lambda <- crossprod(y, w) - sweep(p_yes, 2, colSums(w), "*")
  grad_lambda <- grad_lambda + 1 - 2 * p_yes # logistic prior
  c(grad_beta[, -1], grad_lambda)
}

## ----- Check the gradient ----------------------------------------------------

m <- n_country * (n_class - 1) + length(items) * n_class
theta <- rnorm(m)
h <- 1e-5
grad_fd <- vapply(
  seq_len(m),
  function(j) {
    e <- h * (seq_len(m) == j)
    (logpost(theta + e, y, count, country, n_class) -
      logpost(theta - e, y, count, country, n_class)) /
      (2 * h)
  },
  numeric(1)
)
stopifnot(all.equal(
  logpost_grad(theta, y, count, country, n_class),
  grad_fd,
  tolerance = 1e-6
))

## =============================================================================
## Fit
## =============================================================================

## ----- Mode search from random starts ----------------------------------------

# A latent class likelihood has local maxima, and inla_vanilla() climbs to the
# one nearest init. The start comes from the best of a few nlminb() runs from
# random points, as latent::lca() uses random starts.
starts <- lapply(seq_len(n_starts), function(s) {
  nlminb(
    rnorm(m),
    function(t) -logpost(t, y, count, country, n_class),
    function(t) -logpost_grad(t, y, count, country, n_class)
  )
})
round(sort(vapply(starts, `[[`, numeric(1), "objective")), 2)
best <- starts[[which.min(vapply(starts, `[[`, numeric(1), "objective"))]]

## ----- Name the classes ------------------------------------------------------

# The labels of the classes are arbitrary: the posterior has K! copies of each
# mode, one per labelling, and inla_vanilla() approximates one of them. The
# classes are put in order of their mean P(yes), and the first, the inactive
# class, becomes the reference.
p_best <- unpack(best$par, n_country, n_class, length(items))
ord <- order(colMeans(plogis(p_best$lambda)))
class_names <- c("inactive", "civic", "activist")
beta <- p_best$beta[, ord] - p_best$beta[, ord[1]]
init <- c(beta[, -1], p_best$lambda[, ord])
names(init) <- c(
  paste(levels(country), rep(class_names[-1], each = n_country), sep = "_"),
  paste(items, rep(class_names, each = length(items)), sep = "_")
)

## ----- Fit on the logit scale ------------------------------------------------

fit <- inla_vanilla(
  logpost,
  init,
  grad = logpost_grad,
  y = y,
  count = count,
  country = country,
  n_class = n_class
)
summary(fit)

## =============================================================================
## Results on the probability scale
## =============================================================================

## ----- Class profiles and class sizes ----------------------------------------

# Each draw on the probability scale: P(yes | class) for each item and class,
# then P(class | country), named item_class and country_class.
natural <- function(theta, n_country, n_class, n_item) {
  p <- unpack(theta, n_country, n_class, n_item)
  p_class <- exp(p$beta - row_logsumexp(p$beta))
  c(plogis(p$lambda), p_class)
}
draws <- sampling(
  fit,
  nsamp = nsamp,
  fun = natural,
  n_country = n_country,
  n_class = n_class,
  n_item = length(items)
)
colnames(draws) <- c(
  paste(items, rep(class_names, each = length(items)), sep = "_"),
  paste(levels(country), rep(class_names, each = n_country), sep = "_")
)

post <- data.frame(
  mean = colMeans(draws),
  sd = apply(draws, 2, sd),
  t(apply(draws, 2, quantile, c(0.025, 0.975))),
  check.names = FALSE
)
profile_rows <- seq_len(length(items) * n_class)

# P(yes | class), each item a row and each class a column. The activists are
# the only class that demonstrates, wears badges or donates. The civic class
# signs, boycotts, volunteers and contacts politicians, and the inactive class
# does almost nothing.
profile <- matrix(post$mean[profile_rows], length(items), n_class)
dimnames(profile) <- list(items, class_names)
round(profile, 3)

# Class sizes in each country, with 95% intervals
size <- post[-profile_rows, ]
round(size, 3)
size_mean <- matrix(size$mean, n_country, n_class)

## ----- Class membership of a petition signer ---------------------------------

# P(class | y, country) for a respondent who signed a petition and did nothing
# else. The same answers point to different classes in different countries,
# because the class sizes differ.
signer <- as.integer(items == "sgnptit")
membership <- function(theta, y_new, n_country, n_class) {
  tm <- lca_terms(
    theta,
    matrix(y_new, n_country, length(y_new), byrow = TRUE),
    factor(seq_len(n_country)),
    n_class
  )
  exp(tm$a - tm$lse)
}
signer_draws <- sampling(
  fit,
  nsamp = nsamp,
  fun = membership,
  y_new = signer,
  n_country = n_country,
  n_class = n_class
)
signer_post <- matrix(colMeans(signer_draws), n_country, n_class)
dimnames(signer_post) <- list(levels(country), class_names)
round(signer_post, 3)

## =============================================================================
## Check against latent
## =============================================================================

# latent::lca() fits the same model by penalised maximum likelihood. Its
# default penalties add 1 / (K * 10) log P(class k | country) for each class
# and country, and (s / K) log P(answer | class k) for each answer and class,
# where s is the answer's share in the sample. On the logit scale these are
# Dirichlet priors, much weaker than ours, so its estimates should be within
# a posterior SD or so of the posterior means, except where a class is nearly
# empty and the prior decides how empty.
if (run_latent) {
  latent_fit <- latent::lca(
    data = ess,
    nclasses = n_class,
    multinomial = items,
    covariates = "cntry",
    adjustment = "none",
    verbose = FALSE
  )

  # Item profiles: the items are 0/1, so the second row of each item's table
  # is P(yes). latent's class order is its own, so the classes are matched to
  # ours by their mean P(yes).
  latent_profile <- sapply(
    latent::latInspect(latent_fit, what = "profile")$indicators,
    function(tab) tab[2, ]
  )
  latent_ord <- rownames(latent_profile)[order(rowMeans(latent_profile))]
  latent_profile <- t(latent_profile[latent_ord, ])

  # Class sizes: each country's class logits are the intercept plus its
  # coefficient (the first country, BG, is the baseline).
  b <- latent::latInspect(latent_fit, what = "coefs")[, latent_ord]
  stopifnot(rownames(b)[-1] == paste0("cntry_", levels(country)[-1]))
  eta <- sweep(rbind(0, b[-1, , drop = FALSE]), 2, b[1, ], "+")
  latent_size <- exp(eta - row_logsumexp(eta))

  # Gaps in posterior SDs
  latent_est <- c(latent_profile, latent_size)
  compare_latent <- data.frame(
    mean = post$mean,
    latent = latent_est,
    gap = (latent_est - post$mean) / post$sd,
    row.names = rownames(post)
  )
  print(round(compare_latent, 3))
}

## =============================================================================
## Check against Stan
## =============================================================================

if (run_stan) {
  lca_stan <- "
  data {
    int<lower=1> N; // response patterns
    int<lower=1> J; // items
    int<lower=1> G; // countries
    int<lower=1> K; // classes
    array[N, J] int<lower=0, upper=1> y;
    array[N] int<lower=1, upper=G> country;
    vector[N] count;
  }
  parameters {
    array[G] simplex[K] class_prop; // class sizes, flat Dirichlet
    matrix<lower=0, upper=1>[J, K] p; // P(yes | class), uniform
  }
  model {
    matrix[J, K] log_p = log(p);
    matrix[J, K] log1m_p = log1m(p);
    for (i in 1:N) {
      vector[K] lps = log(class_prop[country[i]]);
      for (j in 1:J) {
        lps += y[i, j] == 1 ? log_p[j]' : log1m_p[j]';
      }
      target += count[i] * log_sum_exp(lps);
    }
  }
  "
  # Every chain starts at the posterior mean, in our order of the classes
  init_stan <- function() {
    list(
      class_prop = lapply(seq_len(n_country), function(g) size_mean[g, ]),
      p = profile
    )
  }
  stan_fit <- rstan::stan(
    model_code = lca_stan,
    data = list(
      N = nrow(y),
      J = length(items),
      G = n_country,
      K = n_class,
      y = y,
      country = as.integer(country),
      count = count
    ),
    chains = 4,
    iter = 3000,
    warmup = 1000,
    seed = 1,
    refresh = 0,
    init = init_stan
  )
  stan_pars <- c(
    sprintf(
      "p[%d,%d]",
      rep(seq_along(items), n_class),
      rep(seq_len(n_class), each = length(items))
    ),
    sprintf(
      "class_prop[%d,%d]",
      rep(seq_len(n_country), n_class),
      rep(seq_len(n_class), each = n_country)
    )
  )
  stan_draws <- as.matrix(stan_fit, pars = c("p", "class_prop"))[, stan_pars]

  # A chain can still wander to another labelling of the classes, so each draw
  # is put in order of mean P(yes), the rule that named the classes.
  relabel <- function(d) {
    p <- matrix(d[profile_rows], length(items))
    s <- matrix(d[-profile_rows], n_country)
    ord <- order(colMeans(p))
    c(p[, ord], s[, ord])
  }
  stan_draws <- t(apply(stan_draws, 1, relabel))
  colnames(stan_draws) <- colnames(draws)

  # Gaps in posterior SDs of the Stan draws: under 0.1 is a close match
  compare <- data.frame(
    mean = post$mean,
    stan_mean = colMeans(stan_draws),
    sd = post$sd,
    stan_sd = apply(stan_draws, 2, sd),
    row.names = colnames(draws)
  )
  compare$mean_gap <- (compare$mean - compare$stan_mean) / compare$stan_sd
  compare$sd_ratio <- compare$sd / compare$stan_sd
  print(round(compare, 3))
}

## =============================================================================
## Plot the class sizes
## =============================================================================

# Posterior means and 95% intervals of the class sizes in each country, the
# countries in order of their share of activists, with Stan's intervals in
# grey beside them when it ran.
lower <- matrix(size[["2.5%"]], n_country)
upper <- matrix(size[["97.5%"]], n_country)
stan_lower <- lower
stan_upper <- upper
if (run_stan) {
  stan_size <- stan_draws[, -profile_rows]
  stan_lower <- matrix(apply(stan_size, 2, quantile, 0.025), n_country)
  stan_upper <- matrix(apply(stan_size, 2, quantile, 0.975), n_country)
}
by_activist <- order(size_mean[, n_class])
at <- seq_len(n_country)
op <- par(mfrow = c(1, n_class), mar = c(4, 3, 2, 1))
for (k in seq_len(n_class)) {
  plot(
    size_mean[by_activist, k],
    at,
    xlim = range(lower[, k], upper[, k], stan_lower[, k], stan_upper[, k]),
    yaxt = "n",
    pch = 19,
    xlab = "Class size",
    ylab = "",
    main = class_names[k]
  )
  axis(2, at, levels(country)[by_activist], las = 1)
  segments(lower[by_activist, k], at, upper[by_activist, k], at)
  if (run_stan) {
    segments(
      stan_lower[by_activist, k],
      at + 0.25,
      stan_upper[by_activist, k],
      at + 0.25,
      col = "grey"
    )
  }
}
par(op)
