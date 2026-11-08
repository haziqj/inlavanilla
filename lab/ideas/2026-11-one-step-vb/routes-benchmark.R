################################################################################
#
# Cheap routes to the VB location on 24 posteriors, one dataset each, against
# the true posterior mean (quadrature, or HMC for three or more parameters).
# Reports the worst-parameter error in posterior SDs and the gradient calls.
# Takes about five minutes. Run from the package root:
#   Rscript lab/ideas/2026-11-one-step-vb/routes-benchmark.R [--smoke]
#
################################################################################

## ----- Configuration ---------------------------------------------------------

source("lab/ideas/2026-11-one-step-vb/common.R")

## ----- Run -------------------------------------------------------------------

if (smoke) {
  models <- models[c("beta_3_11", "banana_b1", "mtcars", "logit_m5")]
}
rows <- list()
meta <- list()
for (name in names(models)) {
  model <- models[[name]]
  t0 <- proc.time()[[3]]
  w <- counting(model)
  lap <- tryCatch(
    laplace_fit(w$lp, w$gr, model$init, ctrl),
    error = function(e) {
      message(name, " laplace failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(lap)) {
    next
  }
  lap_calls <- sum(w$calls())
  m <- length(lap$mode)
  sds <- sqrt(diag(lap$Sigma))
  truth <- if (!is.null(model$draws)) {
    draws_truth(model$draws, lap)
  } else if (m <= 2) {
    grid_truth(w$lp, lap)
  } else {
    hmc_truth(w$lp, w$gr, lap)
  }
  res <- tryCatch(run_routes(model, lap, w), error = function(e) {
    message(name, " routes failed: ", conditionMessage(e))
    NULL
  })
  if (is.null(res)) {
    next
  }
  target <- truth$mean - lap$mode
  cur <- res$shifts$current
  H_quad <- function(s) 0.5 * sum((lap$R %*% s)^2)
  for (route in names(res$shifts)) {
    s <- res$shifts[[route]]
    rows[[length(rows) + 1]] <- data.frame(
      model = name,
      m = m,
      route = route,
      worst = max(abs(s - target) / sds),
      rms = sqrt(mean(((s - target) / sds)^2)),
      vs_current = max(abs(s - cur) / sds),
      logml_gap = H_quad(s) - H_quad(cur),
      calls = res$calls[[route]]
    )
  }
  meta[[name]] <- data.frame(
    model = name,
    m = m,
    lap_calls = lap_calls,
    sweeps = res$sweeps,
    fallback = res$fallback,
    biggest_shift = max(abs(target) / sds),
    truth_se = truth$se,
    secs = round(proc.time()[[3]] - t0)
  )
  message(name, " done in ", round(proc.time()[[3]] - t0), " s")
}
out <- do.call(rbind, rows)
meta <- do.call(rbind, meta)
saveRDS(
  list(out = out, meta = meta),
  file.path(cache_dir, if (smoke) "routes-smoke.rds" else "routes.rds")
)
print(meta, row.names = FALSE)
wide <- reshape(
  out[, c("model", "m", "route", "worst")],
  idvar = c("model", "m"),
  timevar = "route",
  direction = "wide"
)
names(wide) <- sub("worst\\.", "", names(wide))
print(format(wide, digits = 2), row.names = FALSE)
