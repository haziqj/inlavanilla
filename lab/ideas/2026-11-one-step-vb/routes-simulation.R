################################################################################
#
# Replicates routes-benchmark.R over nine model families, 15 datasets each, and
# checks the log marginal likelihood against importance sampling. Run the two
# parts side by side, from the package root, about ten minutes each:
#   Rscript lab/ideas/2026-11-one-step-vb/routes-simulation.R A
#   Rscript lab/ideas/2026-11-one-step-vb/routes-simulation.R B
#
################################################################################

## ----- Configuration ---------------------------------------------------------

source("lab/ideas/2026-11-one-step-vb/common.R")
args <- commandArgs(trailingOnly = TRUE)
part <- args[1]
reps <- if (length(args) > 1) as.integer(args[2]) else 15L
n_warm <- 1000
n_keep <- 8000

## ----- Families --------------------------------------------------------------

families <- if (part == "A") {
  list(
    logit_m10_n80 = gen_logit(10, 80),
    logit_m20_n150 = gen_logit(20, 150),
    logit_m30_n250 = gen_logit(30, 250),
    logit_m20_n600 = gen_logit(20, 600)
  )
} else {
  list(
    pois_m15_n60 = gen_pois(15, 60),
    glmm_J10_n6 = gen_glmm(10, 6),
    treg_m5_n30 = gen_treg(5, 30),
    cfa_p6_n60 = gen_cfa(60),
    cfa_p6_n200 = gen_cfa(200)
  )
}

## ----- Run -------------------------------------------------------------------

rows <- list()
meta <- list()
for (fam in names(families)) {
  for (r in seq_len(reps)) {
    model <- families[[fam]](1000 + r)
    w <- counting(model)
    lap <- tryCatch(
      suppressWarnings(laplace_fit(w$lp, w$gr, model$init, ctrl)),
      error = function(e) NULL
    )
    if (is.null(lap)) {
      message(fam, " ", r, ": laplace failed")
      next
    }
    m <- length(lap$mode)
    sds <- sqrt(diag(lap$Sigma))
    truth <- tryCatch(hmc_truth(w$lp, w$gr, lap), error = function(e) NULL)
    res <- tryCatch(
      suppressWarnings(run_routes(model, lap, w)),
      error = function(e) {
        message(fam, " ", r, ": routes failed: ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(truth) || is.null(res)) {
      next
    }
    set.seed(r)
    iz <- is_logz(w$lp, truth$mean, truth$cov)

    # The cheap pass with the iterative solver as its guard
    a1 <- res$shifts$axis1
    guard <- !all(is.finite(a1)) || max(abs(a1) / sds) >= 1.5 - 1e-9
    res$shifts$axis1_g <- if (guard) res$shifts$current else a1
    res$calls <- c(
      res$calls,
      axis1_g = (2 * m + 1) + if (guard) res$calls[["current"]] else 0
    )

    target <- truth$mean - lap$mode
    cur <- res$shifts$current
    H_quad <- function(s) 0.5 * sum((lap$R %*% s)^2)
    for (route in names(res$shifts)) {
      s <- res$shifts[[route]]
      rows[[length(rows) + 1]] <- data.frame(
        family = fam,
        rep = r,
        m = m,
        route = route,
        worst = max(abs(s - target) / sds),
        rms = sqrt(mean(((s - target) / sds)^2)),
        vs_current = max(abs(s - cur) / sds),
        logml_err = lap$logml + H_quad(s) - iz$logz,
        calls = res$calls[[route]]
      )
    }
    meta[[length(meta) + 1]] <- data.frame(
      family = fam,
      rep = r,
      m = m,
      sweeps = res$sweeps,
      fallback = res$fallback,
      guard = guard,
      shift = max(abs(target) / sds),
      truth_se = truth$se,
      ess = iz$ess,
      lap_err = lap$logml - iz$logz
    )
  }
  message(fam, " done")
  saveRDS(
    list(out = do.call(rbind, rows), meta = do.call(rbind, meta)),
    file.path(cache_dir, paste0("simulation-", part, ".rds"))
  )
}
message("finished ", part)
