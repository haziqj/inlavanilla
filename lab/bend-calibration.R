################################################################################
#
# Cross-model check of the scan-path bend level (October 2026). For each
# parameter of 12 small models and of the 39-parameter Political Democracy
# SEM, compares bend with the error of the skew-normal marginal: how much too
# small its SD is, and how much exact posterior mass its 95% interval holds.
# The references are a grid over +-8 Laplace SDs for the small models (m <= 3)
# and HMC (lab/sem-hmc.R, 4 chains of 10,000 draws) for the SEM. Also gives
# the shortfall that the normal-mean reference family of bend-threshold.md
# predicts from bend, 1 - sqrt(1 - 3 bend / 4).
#
# Result: bend predicts the scan-path narrowing to within about a factor of
# 1.5 in most models, overstating it for GLMs and understating it for
# Puromycin, and the level 0.25 (a predicted 10%) flags every narrowing of 10%
# or more. See bend-threshold.md. The SEM needs lavaan, Rcpp, RcppEigen and
# knitr, and compiles the C++ of vignettes/articles/rcpp.qmd. Takes about 30
# seconds. Run from the package root:
#   Rscript lab/bend-calibration.R
#
################################################################################

## -----  ---------------------------------------------------------------------

suppressMessages(devtools::load_all(".", quiet = TRUE))

## -----  ---------------------------------------------------------------------

# Exact marginals on a grid of +-w Laplace SDs around the mode
grid_marginals <- function(logpost, fit, n = 81, w = 8) {
  m <- length(fit$mode)
  sds <- sqrt(diag(fit$cov))
  axes <- lapply(seq_len(m), function(j) {
    fit$mode[j] + seq(-w, w, length.out = n) * sds[j]
  })
  g <- as.matrix(expand.grid(axes))
  lp <- apply(g, 1, function(x) logpost(setNames(x, names(fit$mode))))
  lp[!is.finite(lp)] <- -Inf
  wt <- exp(lp - max(lp))
  wt <- wt / sum(wt)
  lapply(seq_len(m), function(j) {
    list(x = axes[[j]], p = as.numeric(tapply(wt, g[, j], sum)))
  })
}

# The predicted SD shortfall of the normal-mean reference family
predicted_shortfall <- function(bend) {
  ifelse(0.75 * bend < 1, 1 - sqrt(pmax(0, 1 - 0.75 * bend)), NA)
}

# One row per parameter: bend, the predicted and actual SD shortfall, the
# mean error in exact SDs, and the exact mass inside the 95% interval.
# `mean_sd_cover` gives the exact mean, SD and mass in [lo, hi] of
# parameter j.
compare <- function(case, fit, mean_sd_cover) {
  s <- summary(fit, probs = c(0.025, 0.975))
  bend <- s$diagnostics$bend
  rows <- lapply(seq_along(fit$mode), function(j) {
    e <- mean_sd_cover(j, s$table[j, "2.5%"], s$table[j, "97.5%"])
    data.frame(
      case = case,
      par = names(fit$mode)[j],
      bend = bend[j],
      predicted = predicted_shortfall(bend[j]),
      sd_short = 1 - fit$summary[j, "SD"] / e[["sd"]],
      mean_err = (fit$summary[j, "Mean"] - e[["mean"]]) / e[["sd"]],
      cover95 = e[["cover"]]
    )
  })
  do.call(rbind, rows)
}

compare_grid <- function(case, logpost, init) {
  fit <- inla_vanilla(logpost, init, silent = TRUE)
  ex <- grid_marginals(logpost, fit)
  compare(case, fit, function(j, lo, hi) {
    e <- ex[[j]]
    mu <- sum(e$p * e$x)
    cdf <- approxfun(e$x, cumsum(e$p) - e$p / 2, rule = 2)
    c(mean = mu, sd = sqrt(sum(e$p * (e$x - mu)^2)), cover = cdf(hi) - cdf(lo))
  })
}

## =============================================================================
## Small models --- grid references
## =============================================================================

res <- list()

## -----  ---------------------------------------------------------------------

X <- cbind(1, scale(mtcars$wt), scale(mtcars$hp))
y <- mtcars$am
lp_logit <- function(theta) {
  eta <- as.numeric(X %*% theta)
  sum(dbinom(y, 1, plogis(eta), log = TRUE)) +
    sum(dnorm(theta, 0, 2.5, log = TRUE))
}
res$mtcars <- compare_grid(
  "logistic, mtcars",
  lp_logit,
  c(b0 = 0, b_wt = 0, b_hp = 0)
)

## -----  ---------------------------------------------------------------------

reg_x <- c(-1.75, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, 1.75)
reg_y <- c(1, 0, 2, 1, 3, 4, 2, 6)
lp_reg <- function(theta) {
  eta <- theta[["intercept"]] + theta[["slope"]] * reg_x
  sum(dpois(reg_y, exp(eta), log = TRUE)) +
    dnorm(theta[["intercept"]], 0, 2, log = TRUE) +
    dnorm(theta[["slope"]], 0, 2, log = TRUE)
}
res$pois <- compare_grid(
  "Poisson, n = 8",
  lp_reg,
  c(intercept = 0.2, slope = 0.1)
)

## -----  ---------------------------------------------------------------------

for (n in c(12, 20, 30, 50)) {
  set.seed(4)
  y_n <- rnorm(n, 10, 2)
  lp_norm <- local({
    y_n <- y_n
    function(theta) {
      sum(dnorm(y_n, theta[["mu"]], exp(theta[["log_sigma"]]), log = TRUE)) +
        dnorm(theta[["mu"]], 0, 100, log = TRUE) +
        dnorm(theta[["log_sigma"]], 0, 10, log = TRUE)
    }
  })
  res[[paste0("normal", n)]] <- compare_grid(
    paste0("normal, n = ", n),
    lp_norm,
    c(mu = mean(y_n), log_sigma = log(sd(y_n)))
  )
}

## -----  ---------------------------------------------------------------------

d <- subset(Puromycin, state == "treated")
lp_puro <- function(theta) {
  Vm <- exp(theta[[1]])
  K <- exp(theta[[2]])
  sigma <- exp(theta[[3]])
  mu <- Vm * d$conc / (K + d$conc)
  sum(dnorm(d$rate, mu, sigma, log = TRUE)) +
    dnorm(theta[[1]], log(200), 0.5, log = TRUE) +
    dnorm(theta[[2]], log(0.1), 1, log = TRUE) +
    dnorm(sigma, 0, 20, log = TRUE) +
    log(2) +
    theta[[3]]
}
res$puro <- compare_grid(
  "Puromycin, n = 12",
  lp_puro,
  c(log_Vm = log(200), log_K = log(0.1), log_sigma = log(10))
)

## -----  ---------------------------------------------------------------------

set.seed(5)
y_w <- rweibull(15, shape = 1.5, scale = 2)
lp_weib <- function(theta) {
  k <- exp(theta[["log_shape"]])
  sum(dweibull(y_w, shape = k, scale = exp(theta[["log_scale"]]), log = TRUE)) +
    dnorm(theta[["log_shape"]], 0, 1, log = TRUE) +
    dnorm(theta[["log_scale"]], 0, 2, log = TRUE)
}
res$weib <- compare_grid(
  "Weibull, n = 15",
  lp_weib,
  c(log_shape = 0, log_scale = 0)
)

## -----  ---------------------------------------------------------------------

for (n in c(30, 60, 150)) {
  set.seed(7)
  X_n <- cbind(1, rnorm(n), rnorm(n))
  y_n <- rbinom(n, 1, plogis(X_n %*% c(-0.5, 1.2, -0.8)))
  lp_sim <- local({
    X_n <- X_n
    y_n <- y_n
    function(theta) {
      eta <- as.numeric(X_n %*% theta)
      sum(dbinom(y_n, 1, plogis(eta), log = TRUE)) +
        sum(dnorm(theta, 0, 2.5, log = TRUE))
    }
  })
  res[[paste0("logistic", n)]] <- compare_grid(
    paste0("logistic, n = ", n),
    lp_sim,
    c(b0 = 0, b1 = 0, b2 = 0)
  )
}

## =============================================================================
## Political Democracy SEM --- HMC reference
## =============================================================================

## -----  ---------------------------------------------------------------------

# The R chunks of the article up to `init`, and its C++ file, which the
# article writes in two `cat` chunks
article <- "vignettes/articles/rcpp.qmd"
purled <- tempfile(fileext = ".R")
invisible(knitr::purl(article, purled, documentation = 0, quiet = TRUE))
code <- readLines(purled)
code <- code[seq_len(grep("^fit_inlavaan <- ", code)[1] - 1)]
code <- code[!grepl("^library\\(inlavanilla\\)", code)]
eval(parse(text = code))
qmd <- readLines(article)
fences <- grep("^```", qmd)
cpp_starts <- grep("^```\\{cat, engine.opts", qmd)
cpp <- unlist(lapply(cpp_starts, function(s) {
  qmd[(s + 1):(fences[fences > s][1] - 1)]
}))
cpp_file <- tempfile(fileext = ".cpp")
writeLines(cpp, cpp_file)
Rcpp::sourceCpp(cpp_file)
fit <- inla_vanilla(logpost_cpp, init, grad = grad_cpp, Y = Y, silent = TRUE)

## -----  ---------------------------------------------------------------------

hmc <- readLines("lab/sem-hmc.R")
eval(parse(text = hmc[seq_len(grep("^chains <- lapply", hmc))]))
draws <- do.call(rbind, chains)
res$sem <- compare("PolDem SEM", fit, function(j, lo, hi) {
  x <- draws[, j]
  c(mean = mean(x), sd = sd(x), cover = mean(x >= lo & x <= hi))
})

## =============================================================================
## Results
## =============================================================================

out <- do.call(rbind, res)
rownames(out) <- NULL
num <- vapply(out, is.numeric, logical(1))
out[num] <- lapply(out[num], signif, 3)
small <- out$case != "PolDem SEM"
cat("Small models, by bend:\n")
print(out[small, ][order(out$bend[small]), ], row.names = FALSE)
cat("\nPolDem SEM, by bend:\n")
print(out[!small, -1][order(out$bend[!small]), ], row.names = FALSE)

# How each level sorts the parameters, by actual SD shortfall
for (level in c(0.1, 0.25)) {
  flagged <- out$bend > level
  cat(
    "\nLevel ",
    level,
    ": ",
    sum(flagged),
    " flagged, of which ",
    sum(flagged & out$sd_short < 0.05),
    " under 5% short; ",
    sum(!flagged & out$sd_short >= 0.1),
    " unflagged at 10% or more short (",
    paste(out$par[!flagged & out$sd_short >= 0.1], collapse = ", "),
    ")\n",
    sep = ""
  )
}
