################################################################################
#
# Two checks behind the one-step VB correction. First, the guard: the uncapped
# first step against the current solver, and which guard levels fire. Second,
# the cost: the time and calls of each VB option against the mode and Hessian
# stage, with and without a gradient. Needs the cache of routes-simulation.R
# only for the last table. Run from the package root:
#   Rscript lab/ideas/2026-11-one-step-vb/guard-and-budget.R
#
################################################################################

## ----- Configuration ---------------------------------------------------------

source("lab/ideas/2026-11-one-step-vb/common.R")
n_warm <- 1000
n_keep <- 8000
options(width = 200)

# The first step of the axis rule, with no cap on its length
axis_step_uncapped <- function(w, lap) {
  m <- length(lap$mode)
  e <- eigen(lap$Sigma, symmetric = TRUE)
  A <- sweep(e$vectors, 2, sqrt(pmax(e$values, 0)), "*")
  r <- sqrt(3)
  g0 <- w$gr(lap$mode)
  s <- g0
  for (k in seq_len(m)) {
    s <- s +
      (w$gr(lap$mode + r * A[, k]) + w$gr(lap$mode - r * A[, k]) - 2 * g0) / 6
  }
  as.numeric(lap$Sigma %*% s)
}

# Seconds per run of f(), warm
time_it <- function(f) {
  f()
  n <- 0
  t0 <- proc.time()[[3]]
  repeat {
    f()
    n <- n + 1
    dt <- proc.time()[[3]] - t0
    if (dt > 0.4) break
  }
  dt / n
}

## =============================================================================
## Guard --- the uncapped step against the current solver
## =============================================================================

guard_rows <- function(name, model, r) {
  w <- counting(model)
  lap <- tryCatch(
    suppressWarnings(laplace_fit(w$lp, w$gr, model$init, ctrl)),
    error = function(e) NULL
  )
  if (is.null(lap)) {
    return(NULL)
  }
  sds <- sqrt(diag(lap$Sigma))
  truth <- if (!is.null(model$draws)) {
    draws_truth(model$draws, lap)
  } else if (length(lap$mode) <= 2) {
    grid_truth(w$lp, lap)
  } else {
    hmc_truth(w$lp, w$gr, lap)
  }
  current <- suppressWarnings(
    vb_correct(w$lp, w$gr, lap$mode, lap$Sigma, lap$R, lap$L, ctrl)
  )$shift
  first <- axis_step_uncapped(w, lap)
  move <- max(abs(first) / sds)
  target <- truth$mean - lap$mode
  err <- function(s) {
    c(
      worst = max(abs(s - target) / sds),
      rms = sqrt(mean(((s - target) / sds)^2))
    )
  }
  guarded <- function(level) {
    if (!is.finite(move) || move > level) current else first
  }
  res <- rbind(
    current = err(current),
    axis_uncapped = err(first),
    guard2 = err(guarded(2)),
    guard2.5 = err(guarded(2.5)),
    guard3 = err(guarded(3))
  )
  data.frame(
    family = name,
    rep = r,
    m = length(lap$mode),
    route = rownames(res),
    res,
    move = move,
    row.names = NULL
  )
}

fams <- list(
  logit_m20_n150 = gen_logit(20, 150),
  logit_m30_n250 = gen_logit(30, 250),
  glmm_J10_n6 = gen_glmm(10, 6),
  pois_m15_n60 = gen_pois(15, 60)
)
singles <- c(
  "boundary",
  "skewnorm_a5",
  "skewnorm_a10",
  "eight_schools",
  "poldem",
  "loggamma_1",
  "normal_n5"
)
reps <- if (smoke) 2L else 15L
rows <- list()
for (f in names(fams)) {
  for (r in seq_len(reps)) {
    rows[[length(rows) + 1]] <- guard_rows(f, fams[[f]](1000 + r), r)
  }
  message(f, " done")
}
for (nm in intersect(singles, names(models))) {
  rows[[length(rows) + 1]] <- guard_rows(nm, models[[nm]], 1)
  message(nm, " done")
}
out <- do.call(rbind, rows)
saveRDS(out, file.path(cache_dir, "guard.rds"))

w <- reshape(
  aggregate(worst ~ family + route, out, median),
  idvar = "family",
  timevar = "route",
  direction = "wide"
)
names(w) <- sub("worst\\.", "", names(w))
w[, -1] <- round(w[, -1], 3)
cat("\nMedian worst error, by family\n")
print(w, row.names = FALSE)

## =============================================================================
## Budget --- VB against the mode and Hessian stage
## =============================================================================

budget_row <- function(name, model) {
  w <- counting(model)
  w$reset()
  lap <- suppressWarnings(laplace_fit(w$lp, w$gr, model$init, ctrl))
  n_lap <- sum(w$calls())
  w$reset()
  lap0 <- suppressWarnings(laplace_fit(w$lp, NULL, model$init, ctrl))
  n_lap0 <- sum(w$calls())
  t_lap <- time_it(function() {
    suppressWarnings(laplace_fit(w$lp, w$gr, model$init, ctrl))
  })
  t_lap0 <- time_it(function() {
    suppressWarnings(laplace_fit(w$lp, NULL, model$init, ctrl))
  })
  run_vb <- function(gr) {
    suppressWarnings(
      vb_correct(w$lp, gr, lap$mode, lap$Sigma, lap$R, lap$L, ctrl)
    )
  }
  w$reset()
  run_vb(w$gr)
  n_vb <- sum(w$calls())
  w$reset()
  run_vb(NULL)
  n_vb0 <- sum(w$calls())
  t_vb <- time_it(function() run_vb(w$gr))
  t_vb0 <- time_it(function() run_vb(NULL))
  t_ax <- time_it(function() axis_step_uncapped(w, lap))
  rule <- vb_nodes(lap$Sigma)
  m <- length(lap$mode)
  data.frame(
    model = name,
    m = m,
    grad_budget_calls = n_lap,
    grad_vb_calls = n_vb,
    axis_calls = 2 * m + 1,
    t_vb_over_budget = round(t_vb / t_lap, 2),
    t_axis_over_budget = round(t_ax / t_lap, 2),
    nograd_budget_calls = n_lap0,
    nograd_vb_calls = n_vb0,
    stein_calls = nrow(rule$nodes) - 1,
    t_vb0_over_budget = round(t_vb0 / t_lap0, 2)
  )
}

ms <- list(
  mtcars = models$mtcars,
  logit_m10 = models$logit_m10,
  logit_m20 = models$logit_m20,
  logit_m30 = models$logit_m30,
  logit_m40 = models$logit_m40,
  pois_m15 = models$pois_m15,
  cfa_n60 = gen_cfa(60)(1001),
  glmm = gen_glmm(10, 6)(1001)
)
if (!is.null(models$poldem)) {
  ms$poldem <- models$poldem
}
budget <- do.call(rbind, Map(budget_row, names(ms), ms))
cat("\nCost against the mode and Hessian stage\n")
print(budget, row.names = FALSE)
