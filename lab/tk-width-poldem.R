################################################################################
#
# The width term of the skew-normal scans on the Political Democracy SEM (39
# parameters, n = 75): the linear trace (width = "trace", INLAvaan's), the
# exact width along the conditional mode path (width = "cmp") and the
# Tierney-Kadane Laplace marginal (width = "tk"), each against an HMC
# reference. Trace against cmp isolates the error of the linearisation, cmp
# against tk the error of the straight path. Writes the figures and the
# numbers behind lab/tk-width-poldem.md.
#
# Run from the package root: Rscript lab/tk-width-poldem.R
#
################################################################################

## ----- Configuration ----------------------------------------------------------

# Short HMC chains, for a quick check of the script: Rscript ... --smoke
smoke <- "--smoke" %in% commandArgs(trailingOnly = TRUE)

devtools::load_all(".", quiet = TRUE)
library(Rcpp)
library(RcppEigen)
library(posterior)

lab_dir <- "lab"
cache_dir <- file.path(lab_dir, "tk-width-poldem_cache") # gitignored
fig_dir <- file.path(lab_dir, "tk-width-poldem-figures")
dir.create(cache_dir, showWarnings = FALSE)
dir.create(fig_dir, showWarnings = FALSE)
hmc_file <- file.path(
  cache_dir,
  if (smoke) "hmc-smoke.rds" else "hmc.rds"
)

# HMC: chains, warm-up and kept draws per chain
n_chains <- 4
n_warm <- if (smoke) 500 else 2000
n_keep <- if (smoke) 2000 else 25000
n_time <- if (smoke) 1 else 5 # repeats of each fit, for the timings

# Colours by method (fixed order), with line types as a second encoding
methods <- c("trace", "cmp", "tk")
method_col <- c(trace = "#2a78d6", cmp = "#eb6834", tk = "#1baf7a")
method_lty <- c(trace = 2, cmp = 3, tk = 1)
mcmc_col <- "#52514e"

## =============================================================================
## Model
## =============================================================================

## ----- Data and parameters ----------------------------------------------------

# As in vignettes/articles/rcpp.qmd: INLAvaan's parameterisation and default
# priors, so the marginals are comparable with INLAvaan::asem(model,
# PoliticalDemocracy, meanstructure = TRUE).
data(PoliticalDemocracy, package = "lavaan")
ov <- c("x1", "x2", "x3", paste0("y", 1:8))
Y <- as.matrix(PoliticalDemocracy[, ov])
cov_pairs <- rbind(
  c("y1", "y5"),
  c("y2", "y4"),
  c("y2", "y6"),
  c("y3", "y7"),
  c("y4", "y8"),
  c("y6", "y8")
)
parnames <- c(
  paste0("nu_", ov),
  c("lam_x2", "lam_x3", "a", "b", "c"),
  c("beta_dem60_ind60", "beta_dem65_ind60", "beta_dem65_dem60"),
  paste0("logv_", c("ind60", "dem60", "dem65")),
  paste0("logv_", ov),
  paste0("z_", cov_pairs[, 1], "_", cov_pairs[, 2])
)
init <- c(
  colMeans(Y),
  rep(1, 5),
  rep(0, 3),
  rep(0, 3),
  log(apply(Y, 2, var) / 2),
  rep(0, 6)
)
names(init) <- parnames
m <- length(init)

# Parameter classes, as in INLAvaan's ablation, and the map to the natural
# scale it reports (variances, and residual correlations)
par_class <- rep(
  c(
    "Intercepts",
    "Loadings",
    "Regressions",
    "Latent variances",
    "Residual variances",
    "Residual correlations"
  ),
  c(11, 5, 3, 3, 11, 6)
)
to_natural <- function(x, j) {
  if (j %in% 20:33) {
    exp(x)
  } else if (j %in% 34:39) {
    tanh(x)
  } else {
    x
  }
}

## ----- Log-posterior and gradient in C++ --------------------------------------

# The C++ of the rcpp article, taken from its `cat` chunks, so there is one
# copy of the model.
qmd <- readLines("vignettes/articles/rcpp.qmd")
starts <- grep("^```\\{cat", qmd)
cpp <- unlist(lapply(starts, function(s) {
  e <- s + which(qmd[(s + 1):length(qmd)] == "```")[1]
  qmd[(s + 1):(e - 1)]
}))
Rcpp::sourceCpp(code = paste(cpp, collapse = "\n"))
lp <- function(theta) logpost_cpp(theta, Y)
gr <- function(theta) grad_cpp(theta, Y)

## =============================================================================
## Fits
## =============================================================================

## ----- The three widths and the fallback --------------------------------------

# Each fit is repeated `n_time` times, and the median times are kept: one
# fit of the marginal stage takes under a second here
fit_one <- function(...) {
  reps <- lapply(seq_len(n_time), function(r) {
    t0 <- proc.time()[["elapsed"]]
    fit <- inla_vanilla(lp, init, gr, silent = TRUE, ...)
    fit$wall <- proc.time()[["elapsed"]] - t0
    fit
  })
  fit <- reps[[1]]
  fit$wall <- median(vapply(reps, `[[`, 0, "wall"))
  fit$timing[["marginals"]] <- median(vapply(
    reps,
    function(f) f$timing[["marginals"]],
    0
  ))
  fit
}
# The same mode and VB correction with no scans, to count the calls of the
# marginal stage alone
fit_base <- fit_one(marginals = "gaussian")
fits <- list(
  trace = fit_one(width = "trace"),
  cmp = fit_one(width = "cmp"),
  tk = fit_one(width = "tk"),
  `trace + fallback` = fit_one(width = "trace", tk_fallback = TRUE)
)
cost <- t(vapply(
  fits,
  function(f) {
    c(
      marginal_s = f$timing[["marginals"]],
      total_s = f$wall,
      marginal_grad = f$counts[["grad"]] - fit_base$counts[["grad"]],
      marginal_logpost = f$counts[["logpost"]] - fit_base$counts[["logpost"]]
    )
  },
  numeric(4)
))
print(round(cost, 2))

bend <- fits$trace$diagnostics$marginals$ridge_bend
names(bend) <- parnames
flagged <- parnames[bend > summary_levels$ridge_bend]
cat("Ridge bend flagged:", flagged, "\n")
cat("Fell back:", fits$`trace + fallback`$fallback, "\n")
# Scan points where the conditional optimisation or its Hessian failed
cat("TK points lost:", sum(!is.finite(fits$tk$scan$width)), "\n")

## =============================================================================
## MCMC reference
## =============================================================================

## ----- HMC in whitened coordinates --------------------------------------------

# lab/sem-hmc.R: plain HMC with an identity mass matrix in the whitened
# coordinates of the Laplace approximation, a jittered step size adapted in
# warm-up to an acceptance of 0.8, and 1 to 40 leapfrog steps.
run_hmc <- function() {
  L <- t(chol(fits$trace$cov))
  centre <- fits$trace$mode
  lp_u <- function(u) {
    v <- lp(centre + drop(L %*% u))
    if (is.finite(v)) v else -Inf
  }
  gr_u <- function(u) drop(crossprod(L, gr(centre + drop(L %*% u))))
  chain <- function(seed, l_max = 40) {
    set.seed(seed)
    u <- rnorm(m) * 0.5
    lp0 <- lp_u(u)
    g <- gr_u(u)
    log_eps <- log(0.2)
    out <- matrix(NA_real_, n_keep, m, dimnames = list(NULL, parnames))
    n_div <- 0L
    for (it in seq_len(n_warm + n_keep)) {
      e <- exp(log_eps) * runif(1, 0.8, 1.2)
      p <- rnorm(m)
      h0 <- lp0 - 0.5 * sum(p^2)
      u1 <- u
      g1 <- g
      p1 <- p + 0.5 * e * g1
      ok <- TRUE
      n_steps <- sample.int(l_max, 1)
      for (s in seq_len(n_steps)) {
        u1 <- u1 + e * p1
        g1 <- gr_u(u1)
        if (!all(is.finite(g1))) {
          ok <- FALSE
          break
        }
        if (s < n_steps) p1 <- p1 + e * g1
      }
      a <- 0
      if (ok) {
        p1 <- p1 + 0.5 * e * g1
        lp1 <- lp_u(u1)
        h1 <- lp1 - 0.5 * sum(p1^2)
        if (is.finite(h1)) {
          a <- min(1, exp(h1 - h0))
        }
        if (!is.finite(h1) || h0 - h1 > 1000) n_div <- n_div + 1L
      } else {
        n_div <- n_div + 1L
      }
      if (runif(1) < a) {
        u <- u1
        lp0 <- lp1
        g <- g1
      }
      if (it <= n_warm) {
        log_eps <- log_eps + (a - 0.8) / sqrt(it + 10)
      } else {
        out[it - n_warm, ] <- centre + drop(L %*% u)
      }
    }
    attr(out, "divergences") <- n_div
    out
  }
  t0 <- proc.time()[["elapsed"]]
  chains <- lapply(seq_len(n_chains), chain)
  list(
    chains = chains,
    secs = proc.time()[["elapsed"]] - t0,
    divergences = vapply(chains, attr, 0L, "divergences")
  )
}
hmc <- if (file.exists(hmc_file)) readRDS(hmc_file) else run_hmc()
saveRDS(hmc, hmc_file)
draws_arr <- as_draws_array(aperm(simplify2array(hmc$chains), c(1, 3, 2)))
conv <- summarise_draws(draws_arr, "rhat", "ess_bulk", "ess_tail")
draws <- do.call(rbind, hmc$chains)
cat(
  "HMC:",
  nrow(draws),
  "draws in",
  round(hmc$secs),
  "s; max R-hat",
  round(max(conv$rhat), 4),
  "; min bulk ESS",
  round(min(conv$ess_bulk)),
  "; min tail ESS",
  round(min(conv$ess_tail)),
  "; divergences",
  sum(hmc$divergences),
  "\n"
)

## =============================================================================
## Comparison
## =============================================================================

## ----- Metrics ----------------------------------------------------------------

trapz <- function(x, y) sum((head(y, -1) + tail(y, -1)) * diff(x) / 2)

# Jensen-Shannon discrepancy in % (100 JS / log 2), as INLAvaan's
# compare_mcmc() computes it: a KDE of the draws against the approximate
# density on a common grid, both normalised there. Unlike compare_mcmc(),
# whose grid is the approximation's own, the grid covers both densities, so
# that the mass of an MCMC tail beyond a too-narrow skew-normal counts.
js_pct <- function(x, p, q) {
  eps <- 1e-12
  p <- p / trapz(x, p)
  q <- q / trapz(x, q)
  mix <- 0.5 * (p + q)
  kl <- function(a) trapz(x, a * (log(pmax(a, eps)) - log(pmax(mix, eps))))
  100 * (0.5 * kl(p) + 0.5 * kl(q)) / log(2)
}

sn_of <- function(fit, j) fit$sn[j, c("xi", "omega", "alpha")]
dens <- function(fit, j, x) {
  s <- sn_of(fit, j)
  dsnorm(x, s[["xi"]], s[["omega"]], s[["alpha"]])
}

# Moments and quantiles of a skew-normal marginal, on the scale of logpost
# or, through the monotone map to the natural scale, on that scale
sn_stats <- function(fit, j, natural = FALSE) {
  s <- sn_of(fit, j)
  q <- qsnorm(c(1e-7, 1 - 1e-7), s[["xi"]], s[["omega"]], s[["alpha"]])
  x <- seq(q[1], q[2], length.out = 20001)
  f <- dsnorm(x, s[["xi"]], s[["omega"]], s[["alpha"]])
  f <- f / trapz(x, f)
  g <- if (natural) to_natural(x, j) else x
  mu <- trapz(x, g * f)
  qq <- qsnorm(c(0.025, 0.975), s[["xi"]], s[["omega"]], s[["alpha"]])
  if (natural) {
    qq <- to_natural(qq, j)
  }
  c(
    mean = mu,
    sd = sqrt(trapz(x, (g - mu)^2 * f)),
    q025 = qq[1],
    q975 = qq[2]
  )
}
mc_stats <- function(j, natural = FALSE) {
  d <- if (natural) to_natural(draws[, j], j) else draws[, j]
  q <- quantile(d, c(0.025, 0.975), names = FALSE)
  c(mean = mean(d), sd = sd(d), q025 = q[1], q975 = q[2])
}

# The Kolmogorov-Smirnov distance from the draws
ks_dist <- function(fit, j) {
  s <- sn_of(fit, j)
  d <- sort(draws[, j])
  n <- length(d)
  f <- psnorm(d, s[["xi"]], s[["omega"]], s[["alpha"]])
  max(pmax(abs(f - seq_len(n) / n), abs(f - (seq_len(n) - 1) / n)))
}

grid_for <- function(j) {
  q_mc <- quantile(draws[, j], c(1e-4, 1 - 1e-4), names = FALSE)
  q_sn <- unlist(lapply(fits[methods], function(f) {
    s <- sn_of(f, j)
    qsnorm(c(1e-4, 1 - 1e-4), s[["xi"]], s[["omega"]], s[["alpha"]])
  }))
  r <- range(q_mc, q_sn)
  seq(r[1], r[2], length.out = 1001)
}

## ----- Per parameter ----------------------------------------------------------

fit_names <- names(fits)
rows <- list()
pair_js <- matrix(
  NA_real_,
  m,
  2,
  dimnames = list(parnames, c("trace_vs_cmp", "cmp_vs_tk"))
)
for (j in seq_len(m)) {
  x <- grid_for(j)
  kde <- stats::density(draws[, j], from = min(x), to = max(x), n = length(x))
  f_mc <- stats::approx(kde$x, kde$y, xout = x, rule = 2)$y
  mc_u <- mc_stats(j)
  mc_n <- mc_stats(j, natural = TRUE)
  for (k in fit_names) {
    sn_u <- sn_stats(fits[[k]], j)
    sn_n <- sn_stats(fits[[k]], j, natural = TRUE)
    rows[[length(rows) + 1]] <- data.frame(
      par = parnames[j],
      class = par_class[j],
      method = k,
      width = fits[[k]]$diagnostics$marginals$width[j],
      js = js_pct(x, dens(fits[[k]], j, x), f_mc),
      ks = ks_dist(fits[[k]], j),
      d_mean = abs(sn_u[["mean"]] - mc_u[["mean"]]) / mc_u[["sd"]],
      sd_ratio = sn_u[["sd"]] / mc_u[["sd"]],
      d_q025 = abs(sn_u[["q025"]] - mc_u[["q025"]]) / mc_u[["sd"]],
      d_q975 = abs(sn_u[["q975"]] - mc_u[["q975"]]) / mc_u[["sd"]],
      d_mean_nat = abs(sn_n[["mean"]] - mc_n[["mean"]]) / mc_n[["sd"]],
      sd_ratio_nat = sn_n[["sd"]] / mc_n[["sd"]]
    )
  }
  # Between the methods themselves: the linearisation, then the path
  pair_js[j, ] <- c(
    js_pct(x, dens(fits$trace, j, x), dens(fits$cmp, j, x)),
    js_pct(x, dens(fits$cmp, j, x), dens(fits$tk, j, x))
  )
}
res <- do.call(rbind, rows)

# The floor of the JS discrepancy from MCMC error and the KDE: chains 1-2
# against chains 3-4, on the same grids
half <- rep(1:2, each = 2 * n_keep)
js_floor <- vapply(
  seq_len(m),
  function(j) {
    x <- grid_for(j)
    k1 <- stats::density(
      draws[half == 1, j],
      from = min(x),
      to = max(x),
      n = length(x)
    )
    k2 <- stats::density(
      draws[half == 2, j],
      from = min(x),
      to = max(x),
      n = length(x)
    )
    js_pct(x, k1$y, k2$y)
  },
  numeric(1)
)
names(js_floor) <- parnames
cat(
  "JS floor between halves of the draws (%): mean",
  signif(mean(js_floor), 2),
  ", max",
  signif(max(js_floor), 2),
  "\n"
)
res$method <- factor(res$method, levels = fit_names)
res$bend <- bend[res$par]
saveRDS(
  list(
    res = res,
    pair_js = pair_js,
    cost = cost,
    conv = conv,
    bend = bend,
    js_floor = js_floor,
    hmc = hmc[c("secs", "divergences")]
  ),
  file.path(cache_dir, "results.rds")
)

## ----- Aggregates -------------------------------------------------------------

agg <- function(d) {
  c(
    js = mean(d$js),
    js_median = median(d$js),
    ks = mean(d$ks),
    d_mean = mean(d$d_mean),
    sd_err = mean(abs(d$sd_ratio - 1)),
    d_q025 = mean(d$d_q025),
    d_q975 = mean(d$d_q975),
    d_mean_nat = mean(d$d_mean_nat),
    sd_err_nat = mean(abs(d$sd_ratio_nat - 1))
  )
}
overall <- t(sapply(split(res, res$method), agg))
bulk <- res[!res$par %in% c("logv_x2", "logv_dem65"), ]
others <- t(sapply(split(bulk, bulk$method), agg))
by_class <- t(sapply(
  split(res, list(res$method, res$class), sep = " | "),
  function(d) c(js = mean(d$js), d_mean = mean(d$d_mean))
))
options(width = 160)
cat("\nOverall, all 39 parameters:\n")
print(round(cbind(overall, cost[rownames(overall), ]), 4))
cat("\nWithout logv_x2 and logv_dem65 (the boundary tails):\n")
print(round(others, 4))
cat("\nBy class:\n")
print(round(by_class, 4))
cat("\nJS between methods (%), mean and max:\n")
print(round(rbind(mean = colMeans(pair_js), max = apply(pair_js, 2, max)), 4))
wide <- reshape(
  res[, c("par", "method", "js", "ks", "d_mean", "sd_ratio")],
  idvar = "par",
  timevar = "method",
  direction = "wide"
)
wide$bend <- bend[wide$par]
cat("\nPer parameter:\n")
print(wide, digits = 3, row.names = FALSE)

## =============================================================================
## Figures
## =============================================================================

## ----- Width and height curves ------------------------------------------------

z <- fits$trace$scan$z
lp_max <- fits$trace$scan$logpost[z == 0, 1]
curve_panels <- function(file, curves, ylab, title) {
  grDevices::png(file, width = 2400, height = 2400, res = 200)
  op <- par(
    mfrow = c(7, 6),
    mar = c(1.8, 2.2, 1.6, 0.4),
    oma = c(2, 2, 3.5, 0),
    mgp = c(1.2, 0.35, 0),
    tcl = -0.2,
    cex.axis = 0.7,
    las = 1
  )
  for (j in seq_len(m)) {
    ys <- lapply(curves, function(cv) cv[, j])
    ylim <- range(unlist(ys), finite = TRUE)
    plot(
      range(z),
      ylim,
      type = "n",
      bty = "n",
      xlab = "",
      ylab = "",
      main = parnames[j],
      font.main = 1,
      cex.main = 0.8,
      col.main = if (parnames[j] %in% flagged) "#e34948" else "black"
    )
    abline(h = 0, v = 0, col = "grey88")
    for (k in names(curves)) {
      lines(z, ys[[k]], col = method_col[[k]], lty = method_lty[[k]], lwd = 1.6)
    }
  }
  plot.new()
  legend(
    "center",
    legend = names(curves),
    col = method_col[names(curves)],
    lty = method_lty[names(curves)],
    lwd = 1.6,
    bty = "n",
    cex = 0.9
  )
  mtext("t (Laplace posterior SDs from the mode)", side = 1, outer = TRUE)
  mtext(ylab, side = 2, outer = TRUE, las = 0)
  mtext(title, side = 3, outer = TRUE, line = 1.5, font = 2)
  par(op)
  grDevices::dev.off()
}
curve_panels(
  file.path(fig_dir, "width-curves.png"),
  list(
    trace = fits$trace$scan$width,
    cmp = fits$cmp$scan$width,
    tk = fits$tk$scan$width
  ),
  "Width term gamma_j(t) - gamma_j(0)",
  "Width term: linear trace, exact along the CMP, and at the conditional modes (TK)"
)
curve_panels(
  file.path(fig_dir, "height-curves.png"),
  list(
    cmp = fits$cmp$scan$logpost - lp_max,
    tk = fits$tk$scan$logpost - lp_max
  ),
  "Height: log-posterior minus its maximum",
  "Height: along the straight CMP and at the conditional modes (TK)"
)

## ----- Marginals against HMC --------------------------------------------------

grDevices::png(
  file.path(fig_dir, "marginals-vs-hmc.png"),
  width = 2400,
  height = 2400,
  res = 200
)
op <- par(
  mfrow = c(7, 6),
  mar = c(1.8, 0.6, 1.6, 0.4),
  oma = c(1, 0, 3.5, 0),
  mgp = c(1.2, 0.35, 0),
  tcl = -0.2,
  cex.axis = 0.7
)
for (j in seq_len(m)) {
  x <- grid_for(j)
  kde <- stats::density(draws[, j], from = min(x), to = max(x), n = 512)
  fs <- lapply(methods, function(k) dens(fits[[k]], j, x))
  plot(
    range(x),
    c(0, max(kde$y, unlist(fs)) * 1.05),
    type = "n",
    bty = "n",
    yaxt = "n",
    xlab = "",
    ylab = "",
    main = parnames[j],
    font.main = 1,
    cex.main = 0.8,
    col.main = if (parnames[j] %in% flagged) "#e34948" else "black"
  )
  polygon(
    c(kde$x, rev(kde$x)),
    c(kde$y, numeric(length(kde$y))),
    col = grDevices::adjustcolor(mcmc_col, 0.25),
    border = NA
  )
  for (k in seq_along(methods)) {
    lines(
      x,
      fs[[k]],
      col = method_col[[methods[k]]],
      lty = method_lty[[methods[k]]],
      lwd = 1.4
    )
  }
}
plot.new()
legend(
  "center",
  legend = c("HMC", methods),
  fill = c(grDevices::adjustcolor(mcmc_col, 0.25), NA, NA, NA),
  border = NA,
  col = c(NA, method_col),
  lty = c(NA, method_lty),
  lwd = 1.4,
  bty = "n",
  cex = 0.9
)
mtext(
  "Skew-normal marginals by width method against HMC (red title: ridge bend flagged)",
  side = 3,
  outer = TRUE,
  line = 1.5,
  font = 2
)
par(op)
grDevices::dev.off()

## ----- JS by parameter --------------------------------------------------------

grDevices::png(
  file.path(fig_dir, "js-by-parameter.png"),
  width = 1600,
  height = 2200,
  res = 200
)
js_mat <- sapply(methods, function(k) res$js[res$method == k])
rownames(js_mat) <- parnames
ord <- order(js_mat[, "trace"])
op <- par(mar = c(3, 9, 2.5, 1), mgp = c(1.8, 0.5, 0), las = 1, cex.axis = 0.7)
plot(
  NA,
  xlim = c(0.01, max(js_mat) * 1.2),
  ylim = c(1, m),
  log = "x",
  yaxt = "n",
  bty = "n",
  xlab = "Jensen-Shannon discrepancy against HMC (%)",
  ylab = ""
)
for (k in seq_len(m)) {
  axis(
    2,
    at = k,
    labels = parnames[ord][k],
    lwd = 0,
    col.axis = if (parnames[ord][k] %in% flagged) "#e34948" else "black"
  )
}
abline(h = seq_len(m), col = "grey92")
# The MCMC noise floor of each parameter (chains 1-2 against 3-4)
points(
  js_floor[ord],
  seq_len(m),
  pch = "|",
  col = "grey60",
  cex = 0.9
)
pch <- c(trace = 1, cmp = 2, tk = 16)
for (k in methods) {
  points(
    pmax(js_mat[ord, k], 0.01),
    seq_len(m),
    col = method_col[[k]],
    pch = pch[[k]],
    cex = 0.9
  )
}
legend(
  "bottomright",
  legend = c(methods, "MCMC noise floor"),
  col = c(method_col, "grey60"),
  pch = c(pch, 124),
  bty = "n",
  cex = 0.85
)
title(
  "JS by parameter (sorted by trace; values under 0.01 drawn at 0.01)",
  font.main = 1,
  cex.main = 0.9
)
par(op)
grDevices::dev.off()
