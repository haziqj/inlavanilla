################################################################################
#
# Benchmark of an idea that was tried and rejected (November 2026): replacing
# the exact skew-normal quantile map in the copula draws by a closed-form
# adjustment of the normal scores by the marginal's skewness.
#
# sampling() takes z ~ N(0, R) and maps each column to its marginal through
# z -> Q(Phi(z)), with Q the skew-normal quantile function, solved at 141 knots
# per parameter and splined. The candidate is the quadratic with the marginal's
# mean, SD and skewness, with its slope held between the exact tail slopes. It
# lives only in this script. The outcome is in lab/DESIGN.md, under "The copula
# draws keep the exact quantile map".
#
# The one performance measure is the error of a marginal quantile of the draws
# in posterior SDs, largest between the 0.5% and 99.5% quantiles. For scale, the
# Monte Carlo SE of a 2.5% quantile from 4000 draws is about 0.04 SDs. Takes
# about ten seconds. Run from the package root:
#   Rscript lab/ideas/2026-11-copula-skew/benchmark.R
#
################################################################################

## ----- Configuration ---------------------------------------------------------

suppressMessages(devtools::load_all(".", quiet = TRUE))
ns <- asNamespace("inlavanilla")
qsnorm <- get("qsnorm", ns)
psnorm <- get("psnorm", ns)
sn_score_map <- get("sn_score_map", ns)

alphas <- c(0.25, 0.5, 1, 1.5, 2, 3, 3.5, 5, 10, 20)
out_dir <- "lab/ideas/2026-11-copula-skew"

# Mean, SD and skewness of SN(xi, omega, alpha)
sn_moments <- function(xi, omega, alpha) {
  delta <- alpha / sqrt(1 + alpha^2)
  b2 <- 2 * delta^2 / pi
  list(
    mean = xi + omega * delta * sqrt(2 / pi),
    sd = omega * sqrt(1 - b2),
    skew = sign(alpha) * (4 - pi) / 2 * b2^1.5 / (1 - b2)^1.5
  )
}

## =============================================================================
## Maps from a normal score to SN(xi, omega, alpha)
## =============================================================================

# The package's map, exact to within 1e-6 posterior SDs
map_exact <- function(z, xi, omega, alpha) {
  if (alpha == 0) {
    return(xi + omega * z)
  }
  sn_score_map(xi, omega, alpha)(z)
}

# For the record: Cornish-Fisher, z + (gamma / 6) (z^2 - 1), with its variance
# and skewness only approximately right
map_cornish_fisher <- function(z, xi, omega, alpha) {
  m <- sn_moments(xi, omega, alpha)
  m$mean + m$sd * (z + m$skew / 6 * (z^2 - 1))
}

# For the record: the quadratic with the marginal's three moments, without the
# slope limit, held constant below its turning point
map_plain_quadratic <- function(z, xi, omega, alpha) {
  m <- sn_moments(xi, omega, alpha)
  c2 <- sqrt(2) * cos(acos(-m$skew / (2 * sqrt(2))) / 3 - 2 * pi / 3)
  c1 <- sqrt(1 - 2 * c2^2)
  if (c2 != 0) {
    turn <- -c1 / (2 * c2)
    z <- if (c2 > 0) pmax(z, turn) else pmin(z, turn)
  }
  m$mean + m$sd * (c1 * z + c2 * (z^2 - 1))
}

# The rejected candidate: the quadratic above, with its slope held between the
# slopes of the exact map in the far tails
map_skew <- function(z, xi, omega, alpha) {
  if (alpha == 0) {
    return(xi + omega * z)
  }
  mom <- sn_moments(xi, omega, alpha)
  c2 <- sqrt(2) * cos(acos(-mom$skew / (2 * sqrt(2))) / 3 - 2 * pi / 3)
  c1 <- sqrt(1 - 2 * c2^2)
  heavy <- omega / mom$sd
  light <- heavy / sqrt(1 + alpha^2)
  limits <- sort((c(light, heavy) - c1) / (2 * c2))
  z0 <- pmin(pmax(z, limits[1]), limits[2])
  mom$mean +
    mom$sd * (c1 * z0 + c2 * (z0^2 - 1) + (c1 + 2 * c2 * z0) * (z - z0))
}

## =============================================================================
## Experiment 1 --- the map against the exact quantiles, no Monte Carlo
## =============================================================================

# One row per alpha, all errors in posterior SDs. `max_err` is the performance
# measure and `ks` the largest gap between the distribution functions.
score_map <- function(map, alpha_grid = alphas) {
  probs <- seq(0.005, 0.995, by = 0.0025)
  rows <- lapply(alpha_grid, function(a) {
    s <- sn_moments(0, 1, a)$sd
    err <- (map(stats::qnorm(probs), 0, 1, a) - qsnorm(probs, 0, 1, a)) / s
    tails <- c(0.001, 0.025, 0.975, 0.999)
    tail_err <- (map(stats::qnorm(tails), 0, 1, a) - qsnorm(tails, 0, 1, a)) / s
    zg <- seq(-8, 8, by = 0.001)
    x <- map(zg, 0, 1, a)
    data.frame(
      alpha = a,
      skew = sn_moments(0, 1, a)$skew,
      max_err = max(abs(err)),
      q0.1 = tail_err[1],
      q2.5 = tail_err[2],
      q97.5 = tail_err[3],
      q99.9 = tail_err[4],
      ks = max(abs(psnorm(x, 0, 1, a) - stats::pnorm(zg))),
      monotone = all(diff(x) >= 0)
    )
  })
  do.call(rbind, rows)
}

## =============================================================================
## Experiment 2 --- closed-form posteriors, end to end
## =============================================================================

# One-parameter posteriors with exact quantiles on the theta scale, and the
# natural-scale functional plogis() or exp() with an exact mean.
closed_form <- list(
  binom_0_of_10 = list(a = 1, b = 11, kind = "logit"),
  binom_2_of_30 = list(a = 3, b = 29, kind = "logit"),
  gammapois_sum0_n5 = list(a = 1, b = 6, kind = "log"),
  gammapois_sum4_n5 = list(a = 5, b = 6, kind = "log")
)

cf_logpost <- function(spec) {
  if (spec$kind == "logit") {
    function(t) {
      p <- stats::plogis(t[[1]])
      stats::dbeta(p, spec$a, spec$b, log = TRUE) + log(p) + log1p(-p)
    }
  } else {
    function(t) spec$a * t[[1]] - spec$b * exp(t[[1]])
  }
}

cf_quantile <- function(spec, p) {
  if (spec$kind == "logit") {
    stats::qlogis(stats::qbeta(p, spec$a, spec$b))
  } else {
    log(stats::qgamma(p, spec$a, spec$b))
  }
}

cf_density <- function(spec, x) {
  if (spec$kind == "logit") {
    p <- stats::plogis(x)
    stats::dbeta(p, spec$a, spec$b) * p * (1 - p)
  } else {
    stats::dgamma(exp(x), spec$a, spec$b) * exp(x)
  }
}

cf_fun <- function(spec) if (spec$kind == "logit") stats::plogis else exp
cf_mean <- function(spec) {
  if (spec$kind == "logit") spec$a / (spec$a + spec$b) else spec$a / spec$b
}

cf_fit <- function(spec) {
  inla_vanilla(cf_logpost(spec), c(th = 0), silent = TRUE)
}

# The largest quantile error against the truth, for the exact map and the
# skewness map, and the mean of the natural-scale parameter from each
experiment_closed_form <- function(seed = 1, nsamp = 2e6L) {
  rows <- lapply(names(closed_form), function(nm) {
    spec <- closed_form[[nm]]
    sn <- cf_fit(spec)$sn[1, ]
    probs <- c(0.01, 0.025, 0.5, 0.975, 0.99)
    z <- stats::qnorm(probs)
    truth <- cf_quantile(spec, probs)
    sd_truth <- diff(cf_quantile(spec, stats::pnorm(c(-1, 1)))) / 2
    old <- map_exact(z, sn[["xi"]], sn[["omega"]], sn[["alpha"]])
    new <- map_skew(z, sn[["xi"]], sn[["omega"]], sn[["alpha"]])
    set.seed(seed)
    zz <- stats::rnorm(nsamp)
    f <- cf_fun(spec)
    fold <- f(map_exact(zz, sn[["xi"]], sn[["omega"]], sn[["alpha"]]))
    fnew <- f(map_skew(zz, sn[["xi"]], sn[["omega"]], sn[["alpha"]]))
    data.frame(
      model = nm,
      alpha = sn[["alpha"]],
      err_exact_map = max(abs(old - truth)) / sd_truth,
      err_skew_map = max(abs(new - truth)) / sd_truth,
      map_gap = max(abs(new - old)) / sd_truth,
      mean_exact = cf_mean(spec),
      mean_exact_map = mean(fold),
      mean_skew_map = mean(fnew),
      mc_se = stats::sd(fold) / sqrt(nsamp)
    )
  })
  do.call(rbind, rows)
}

## =============================================================================
## Experiment 3 --- the copula, under common random numbers
## =============================================================================

# The same normals through both maps. Ranks are unchanged by an increasing map,
# so Spearman's rho must agree exactly, and Pearson's may shift a little.
experiment_copula <- function(seed = 1, nsamp = 1e5L) {
  x <- c(0.2, 0.5, 0.9, 1.3, 1.6, 2.0, 2.4, 2.9, 3.3, 3.8)
  d <- c(1, 0, 1, 1, 0, 0, 1, 0, 1, 1)
  y <- c(0, 1, 0, 2, 1, 3, 2, 4, 6, 5)
  X <- cbind(1, x, d, deparse.level = 0)
  lp <- function(t) {
    eta <- drop(X %*% t)
    sum(stats::dpois(y, exp(eta), log = TRUE)) +
      sum(stats::dnorm(t, 0, 2, log = TRUE))
  }
  fit <- inla_vanilla(lp, c(b0 = 0, b1 = 0, b2 = 0), silent = TRUE)
  sn <- fit$sn
  m <- nrow(sn)
  set.seed(seed)
  z <- matrix(stats::rnorm(nsamp * m), nsamp, m) %*%
    chol(stats::cov2cor(fit$cov))
  draw <- function(map) {
    vapply(
      seq_len(m),
      function(j) map(z[, j], sn[j, "xi"], sn[j, "omega"], sn[j, "alpha"]),
      numeric(nsamp)
    )
  }
  old <- draw(map_exact)
  new <- draw(map_skew)
  lower <- lower.tri(diag(m))
  list(
    alpha = sn[, "alpha"],
    spearman_gap = max(abs(
      stats::cor(old, method = "spearman") -
        stats::cor(new, method = "spearman")
    )),
    pearson_exact_map = stats::cor(old)[lower],
    pearson_skew_map = stats::cor(new)[lower]
  )
}

## =============================================================================
## Experiment 4 --- exp() of a wide parameter, no Monte Carlo
## =============================================================================

# E exp(theta) has a closed form under SN(0, omega, alpha). A map with too much
# growth in the heavy tail shows up here first.
experiment_exp <- function() {
  e_exact <- function(omega, a) {
    d <- a / sqrt(1 + a^2)
    2 * exp(omega^2 / 2) * stats::pnorm(d * omega)
  }
  e_map <- function(map, omega, a) {
    stats::integrate(
      function(z) exp(map(z, 0, omega, a)) * stats::dnorm(z),
      -9,
      9,
      subdivisions = 3000,
      rel.tol = 1e-9
    )$value
  }
  g <- expand.grid(omega = c(1, 2, 3), alpha = c(3, -3, 10))
  ratio <- function(map) {
    mapply(function(o, a) e_map(map, o, a) / e_exact(o, a), g$omega, g$alpha)
  }
  g$plain_quadratic <- ratio(map_plain_quadratic)
  g$skew_map <- ratio(map_skew)
  g
}

## =============================================================================
## Experiment 5 --- cost
## =============================================================================

time_map <- function(map, alpha, z, repeats = 5L) {
  elapsed <- system.time(
    for (i in seq_len(repeats)) {
      for (a in alpha) {
        map(z, 0, 1, a)
      }
    }
  )[["elapsed"]]
  1000 * elapsed / (repeats * length(alpha))
}

experiment_cost <- function() {
  set.seed(1)
  alpha <- stats::rnorm(400, 0, 2)
  grid <- data.frame(
    nsamp = c(1000L, 100000L),
    ms_exact_map = NA,
    ms_skew_map = NA
  )
  for (i in seq_len(nrow(grid))) {
    z <- stats::rnorm(grid$nsamp[i])
    grid$ms_exact_map[i] <- time_map(map_exact, alpha, z, 3L)
    grid$ms_skew_map[i] <- time_map(map_skew, alpha, z, 3L)
  }
  grid
}

## =============================================================================
## Figure --- the draws against the skew-normal and the exact posterior
## =============================================================================

figure_draws <- function(file, nsamp = 4e5L) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(file, width = 2400, height = 1500, res = 240)
  on.exit(grDevices::dev.off())
  graphics::par(mfrow = c(2, 3), mar = c(4, 4, 2.5, 0.5), mgp = c(2.2, 0.7, 0))
  cols <- c(exact = "grey70", skew = "#b8751f")
  set.seed(1)
  z <- stats::rnorm(nsamp)

  draw_panel <- function(x, target, old, new, main) {
    d_old <- stats::density(old, adjust = 0.6, n = 1024)
    d_new <- stats::density(new, adjust = 0.6, n = 1024)
    graphics::plot(
      x,
      target,
      type = "l",
      lty = 2,
      main = main,
      xlab = "Parameter",
      ylab = "Density",
      ylim = c(0, max(target, d_old$y) * 1.1)
    )
    graphics::lines(d_old, col = cols[["exact"]], lwd = 3)
    graphics::lines(d_new, col = cols[["skew"]], lwd = 1.5)
    graphics::lines(x, target, lty = 2)
  }

  for (a in c(1, 3, 10)) {
    x <- seq(-1.5, 3.5, length.out = 500)
    draw_panel(
      x,
      2 * stats::dnorm(x) * stats::pnorm(a * x),
      map_exact(z, 0, 1, a),
      map_skew(z, 0, 1, a),
      bquote("Skew-normal, " * alpha == .(a))
    )
    if (a == 1) {
      graphics::legend(
        "topright",
        c("Exact map", "Skewness map", "Skew-normal"),
        col = c(cols[["exact"]], cols[["skew"]], "black"),
        lwd = c(3, 1.5, 1),
        lty = c(1, 1, 2),
        bty = "n",
        cex = 0.8
      )
    }
  }
  titles <- c(
    binom_0_of_10 = "0 of 10, logit scale",
    binom_2_of_30 = "2 of 30, logit scale",
    gammapois_sum4_n5 = "Gamma-Poisson, log scale"
  )
  for (nm in names(titles)) {
    spec <- closed_form[[nm]]
    sn <- cf_fit(spec)$sn[1, ]
    old <- map_exact(z, sn[["xi"]], sn[["omega"]], sn[["alpha"]])
    new <- map_skew(z, sn[["xi"]], sn[["omega"]], sn[["alpha"]])
    x <- seq(
      cf_quantile(spec, 1e-4),
      cf_quantile(spec, 1 - 1e-5),
      length.out = 500
    )
    draw_panel(x, cf_density(spec, x), old, new, titles[[nm]])
    if (nm == "binom_0_of_10") {
      graphics::legend(
        "topleft",
        c("Exact map", "Skewness map", "Exact posterior"),
        col = c(cols[["exact"]], cols[["skew"]], "black"),
        lwd = c(3, 1.5, 1),
        lty = c(1, 1, 2),
        bty = "n",
        cex = 0.8
      )
    }
  }
}

## =============================================================================
## Results
## =============================================================================

if (sys.nframe() == 0L) {
  cat("\nExperiment 1: quantile error in posterior SDs\n")
  for (nm in c("map_cornish_fisher", "map_plain_quadratic", "map_skew")) {
    cat("\n", nm, "\n", sep = "")
    print(format(score_map(get(nm)), digits = 2), row.names = FALSE)
  }

  cat("\nExperiment 2: closed-form posteriors\n")
  print(format(experiment_closed_form(), digits = 3), row.names = FALSE)

  cat("\nExperiment 3: copula under common random numbers\n")
  str(experiment_copula())

  cat("\nExperiment 4: E exp(theta), ratio to the exact value\n")
  print(signif(experiment_exp(), 4), row.names = FALSE)

  cat("\nExperiment 5: cost in ms per parameter\n")
  print(experiment_cost(), row.names = FALSE)

  figure_draws(file.path(out_dir, "draws-vs-exact.png"))
}
