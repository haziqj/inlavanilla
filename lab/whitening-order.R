################################################################################
#
# Does the fit depend on the order of the parameters? (November 2026). The
# route without a gradient whitened with lap$L, the triangular factor from the
# Cholesky of the precision, which depends on the order: the volume correction
# (gamma1_lp(), whitened_lp()), the direction of the ridge bend (bend_along())
# and the coordinates of the VB fallback. The other stages use the principal
# axes, which do not. This script fits the 28 models of lab/bend-screen.R
# (logpost only, no gradient) in their own order and in three random orders,
# by default and with the VB fallback forced (vb = "full", vb_maxit = 1), and
# reports the largest change across the orders: in posterior SDs for the
# means, quantiles and VB means, as a ratio for the SDs, in nats for the ridge
# bend and the log marginal likelihood.
#
# Result: see "One whitening factor" in lab/DESIGN.md.
#
# Takes a few minutes, and needs lavaan, Rcpp, RcppEigen and knitr for the
# SEM. Run from the package root, with SMOKE=1 for three models. The optional
# argument saves the fits in their own order, to compare two versions:
#   Rscript lab/whitening-order.R [out.rds]
#
################################################################################

## ----- Configuration ----------------------------------------------------------

smoke <- nzchar(Sys.getenv("SMOKE"))
out_file <- commandArgs(trailingOnly = TRUE)[1]
n_orders <- 3L
seed <- 2026L
settings <- list(
  default = list(),
  fallback = list(vb = "full", control = list(vb_maxit = 1L))
)

## ----- The models of lab/bend-screen.R ----------------------------------------

# Everything before its fits: its configuration and the 28 models
script <- readLines("lab/bend-screen.R")
eval(parse(text = script[seq_len(grep("^## Fits", script) - 2L)]))
if (smoke) {
  cases <- cases[1:3]
}

## =============================================================================
## Fits
## =============================================================================

# A fit with the parameters in the order `perm`, reported in their own order
fit_in_order <- function(case, perm, setting) {
  own <- names(case$init)
  logpost <- function(theta) case$logpost(theta[order(perm)])
  args <- c(
    list(logpost, case$init[perm], transform = FALSE, silent = TRUE),
    setting
  )
  fit <- tryCatch(do.call(inla_vanilla, args), error = function(e) NULL)
  if (is.null(fit)) {
    return(NULL)
  }
  bend <- fit$diagnostics$marginals$ridge_bend
  names(bend) <- names(fit$mode)
  list(
    summary = fit$summary[own, , drop = FALSE],
    mean_vb = fit$mean_vb[own],
    bend = bend[own],
    logml = fit$logml,
    fallback = isTRUE(fit$diagnostics$vb$fallback)
  )
}

# The largest change from `ref` to `alt`
compare <- function(ref, alt) {
  sd0 <- ref$summary[, "SD"]
  q <- c("2.5%", "97.5%")
  c(
    mean = max(abs(alt$summary[, "Mean"] - ref$summary[, "Mean"]) / sd0),
    sd = max(abs(alt$summary[, "SD"] / sd0 - 1)),
    quantile = max(abs(alt$summary[, q] - ref$summary[, q]) / sd0),
    vb_mean = max(abs(alt$mean_vb - ref$mean_vb) / sd0),
    bend = max(abs(alt$bend - ref$bend)),
    logml = abs(alt$logml - ref$logml)
  )
}

set.seed(seed)
orders <- lapply(cases, function(case) {
  replicate(n_orders, sample(length(case$init)), simplify = FALSE)
})

refs <- list()
rows <- list()
for (i in seq_along(cases)) {
  case <- cases[[i]]
  m <- length(case$init)
  for (s in names(settings)) {
    ref <- fit_in_order(case, seq_len(m), settings[[s]])
    refs[[paste(case$name, s, sep = " | ")]] <- ref
    alts <- lapply(orders[[i]], function(p) {
      fit_in_order(case, p, settings[[s]])
    })
    ok <- !vapply(alts, is.null, NA)
    diffs <- if (!is.null(ref) && any(ok)) {
      apply(
        vapply(alts[ok], function(a) compare(ref, a), numeric(6)),
        1,
        max
      )
    } else {
      stats::setNames(
        rep(NA_real_, 6),
        c("mean", "sd", "quantile", "vb_mean", "bend", "logml")
      )
    }
    rows[[length(rows) + 1L]] <- data.frame(
      case = case$name,
      setting = s,
      m = m,
      failed = sum(!ok) + is.null(ref),
      fallback = isTRUE(ref$fallback),
      t(diffs)
    )
  }
}
d <- do.call(rbind, rows)
names(d)[6:11] <- c("mean", "sd", "quantile", "vb_mean", "bend", "logml")

## =============================================================================
## Results
## =============================================================================

cat("Models:", length(cases), " fits that failed:", sum(d$failed), "\n\n")
shown <- d
shown[6:11] <- lapply(shown[6:11], signif, 2)
print(shown, row.names = FALSE)
cat("\nLargest change over all models, by setting:\n")
print(signif(
  sapply(split(d[6:11], d$setting), function(x) {
    apply(x, 2, max, na.rm = TRUE)
  }),
  2
))
if (!is.na(out_file)) {
  saveRDS(list(d = d, refs = refs), out_file)
}
