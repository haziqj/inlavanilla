################################################################################
#
# Why does one-step VB fail on some models, and can the fit tell? For each
# dataset: the error of the one-step mean and evidence against the truth, the
# gain ratio of the package, and rival diagnostics from the first pass. Run from
# the package root, parts A and B side by side, then C, about 25 minutes in all:
#   Rscript lab/ideas/2026-11-one-step-vb/why-one-step-fails.R A
#   Rscript lab/ideas/2026-11-one-step-vb/why-one-step-fails.R B
#   Rscript lab/ideas/2026-11-one-step-vb/why-one-step-fails.R C
#
################################################################################

## ----- Configuration ---------------------------------------------------------

source("lab/ideas/2026-11-one-step-vb/common.R")
args <- commandArgs(trailingOnly = TRUE)
part <- args[1]
reps <- if (length(args) > 1 && part != "C") as.integer(args[2]) else 15L
n_warm <- 1000
n_keep <- 6000

## ----- One dataset -----------------------------------------------------------

grid_truth_z <- function(lp, lap) {
  m <- length(lap$mode)
  L <- lap$L
  z1 <- if (m == 1) seq(-40, 40, by = 0.005) else seq(-30, 30, by = 0.1)
  z <- as.matrix(expand.grid(rep(list(z1), m)))
  theta <- sweep(z %*% t(L), 2, lap$mode, "+")
  v <- vapply(seq_len(nrow(theta)), function(i) lp(theta[i, ]), 0)
  v[!is.finite(v)] <- -Inf
  top <- max(v)
  w <- exp(v - top)
  cell <- prod(diff(z1)[1]^m) * abs(det(L))
  list(mean = colSums(theta * w) / sum(w), logz = top + log(sum(w) * cell))
}

study <- function(name, r, model, truth = NULL, logz = NULL) {
  w <- counting(model)
  lap <- tryCatch(
    suppressWarnings(laplace_fit(w$lp, w$gr, model$init, ctrl)),
    error = function(e) NULL
  )
  if (is.null(lap)) {
    return(NULL)
  }
  m <- length(lap$mode)
  sds <- sqrt(diag(lap$Sigma))
  A <- principal_axes(lap$Sigma)

  # The truth
  if (is.null(truth)) {
    if (!is.null(model$draws)) {
      tr <- draws_truth(model$draws, lap)
      set.seed(r)
      logz <- is_logz(w$lp, tr$mean, tr$cov)$logz
      truth <- tr$mean
    } else if (m <= 2) {
      tr <- grid_truth_z(w$lp, lap)
      truth <- tr$mean
      logz <- tr$logz
    } else {
      tr <- hmc_truth(w$lp, w$gr, lap)
      truth <- tr$mean
      set.seed(r)
      logz <- is_logz(w$lp, tr$mean, tr$cov)$logz
    }
  }

  # The first pass, in the whitened coordinates of the principal axes
  rad <- sqrt(3)
  white <- function(g) as.numeric(crossprod(A, g))
  pass <- function(d) {
    mu <- lap$mode + as.numeric(A %*% d)
    g0 <- white(w$gr(mu))
    gp <- gm <- matrix(NA_real_, m, m)
    for (k in seq_len(m)) {
      gp[, k] <- white(w$gr(mu + rad * A[, k]))
      gm[, k] <- white(w$gr(mu - rad * A[, k]))
    }
    list(
      f = g0 + rowSums(gp + gm - 2 * g0) / 6,
      J = sqrt(3) / 6 * (gp - gm)
    )
  }
  obj <- function(d) {
    mu <- lap$mode + as.numeric(A %*% d)
    v0 <- w$lp(mu)
    vals <- vapply(
      seq_len(m),
      function(k) w$lp(mu + rad * A[, k]) + w$lp(mu - rad * A[, k]) - 2 * v0,
      0
    )
    v0 + sum(vals) / 6
  }

  p0 <- pass(numeric(m))
  f <- p0$f
  if (!all(is.finite(f)) || !all(is.finite(p0$J))) {
    return(NULL)
  }
  Phi <- 0.5 * (p0$J + t(p0$J)) + diag(m)
  ev <- eigen(Phi, symmetric = TRUE)$values
  d_newton <- tryCatch(
    as.numeric(solve(diag(m) - Phi, f)),
    error = function(e) rep(NA_real_, m)
  )
  shift_on <- as.numeric(A %*% f)
  move <- max(abs(shift_on) / sds)
  G0 <- obj(numeric(m))
  G1 <- obj(f)
  pred <- 0.5 * sum(f^2)
  p1 <- pass(f)
  step2 <- as.numeric(A %*% (p1$f))
  shift_nw <- as.numeric(A %*% d_newton)

  cur <- suppressWarnings(
    vb_correct(w$lp, w$gr, lap$mode, lap$Sigma, lap$R, lap$L, ctrl)
  )
  real <- function(gr, method) {
    tryCatch(
      suppressWarnings(vb_correct(
        w$lp,
        gr,
        lap$mode,
        lap$Sigma,
        lap$R,
        lap$L,
        ctrl,
        method = method,
        budget = Inf
      )),
      error = function(e) NULL
    )
  }
  os_g <- real(w$gr, "onestep")
  os_s <- real(NULL, "onestep")
  fu_s <- real(NULL, "full")
  fit0 <- suppressWarnings(
    inla_vanilla(w$lp, model$init, w$gr, vb = "off", silent = TRUE)
  )
  mom <- sn_moments(fit0$sn[, "xi"], fit0$sn[, "omega"], fit0$sn[, "alpha"])
  shift_sn <- as.numeric(mom$mean) - lap$mode

  target <- truth - lap$mode
  err <- function(s) (s - target) / sds
  worst <- function(s) max(abs(err(s)))
  rms <- function(s) sqrt(mean(err(s)^2))
  L <- lap$logml
  data.frame(
    family = name,
    rep = r,
    m = m,
    true_shift = max(abs(target) / sds),
    worst_mode = worst(numeric(m)),
    worst_on = worst(shift_on),
    worst_full = worst(cur$shift),
    worst_newton = worst(shift_nw),
    worst_sn = worst(shift_sn),
    rms_on = rms(shift_on),
    rms_full = rms(cur$shift),
    rms_newton = rms(shift_nw),
    le_laplace = L - logz,
    le_on = L + pred - logz,
    le_full = L + cur$kld_global - logz,
    le_actual = L + (G1 - G0) - logz,
    move = move,
    curv_max = max(abs(ev)),
    curv_dir = sum(f * (Phi %*% f)) / sum(f^2),
    ratio = (G1 - G0) / pred,
    step2 = max(abs(step2) / sds),
    mean_gap = max(abs(shift_on - shift_sn) / sds),
    pred = pred,
    gain_real = if (is.null(os_g)) NA else os_g$gain_ratio,
    one_pass_g = if (is.null(os_g)) NA else os_g$one_pass,
    worst_g = if (is.null(os_g)) NA else worst(os_g$shift),
    gain_s = if (is.null(os_s)) NA else os_s$gain_ratio,
    one_pass_s = if (is.null(os_s)) NA else os_s$one_pass,
    worst_s = if (is.null(os_s)) NA else worst(os_s$shift),
    worst_full_s = if (is.null(fu_s)) NA else worst(fu_s$shift),
    gap_s = if (is.null(os_s)) NA else max(abs(os_s$shift - shift_sn) / sds),
    le_s = if (is.null(os_s)) NA else L + os_s$kld_global - logz
  )
}

## ----- Run -------------------------------------------------------------------

families <- switch(
  part,
  A = list(
    logit_m10_n80 = gen_logit(10, 80),
    logit_m20_n150 = gen_logit(20, 150),
    logit_m30_n250 = gen_logit(30, 250),
    logit_m20_n600 = gen_logit(20, 600)
  ),
  B = list(
    pois_m15_n60 = gen_pois(15, 60),
    glmm_J10_n6 = gen_glmm(10, 6),
    treg_m5_n30 = gen_treg(5, 30),
    cfa_p6_n60 = gen_cfa(60),
    cfa_p6_n200 = gen_cfa(200)
  ),
  C = list()
)
rows <- list()
for (fam in names(families)) {
  for (r in seq_len(reps)) {
    rows[[length(rows) + 1]] <- tryCatch(
      study(fam, r, families[[fam]](1000 + r)),
      error = function(e) {
        message(fam, " ", r, ": ", conditionMessage(e))
        NULL
      }
    )
  }
  message(fam, " done")
  saveRDS(
    do.call(rbind, rows),
    file.path(cache_dir, paste0("why-", part, ".rds"))
  )
}

if (part == "C") {
  for (nm in names(models)) {
    rows[[length(rows) + 1]] <- tryCatch(
      study(nm, 1, models[[nm]]),
      error = function(e) {
        message(nm, ": ", conditionMessage(e))
        NULL
      }
    )
    message(nm, " done")
  }
  # The diagnostic-test models, with a numerical gradient. Their truth is a
  # long random-walk chain for the means and importance sampling for log Z,
  # cached after the first run.
  qmd <- readLines("vignettes/articles/diagnostic-tests.qmd")
  starts <- grep("^```\\{r\\}", qmd)
  chunks <- lapply(starts, function(i) {
    e <- i + which(qmd[(i + 1):length(qmd)] == "```")[1]
    qmd[(i + 1):(e - 1)]
  })
  code <- unlist(chunks[seq_len(
    which(grepl(
      "^start <- c\\(",
      sapply(
        chunks,
        function(x) paste(x[!startsWith(x, "#|")], collapse = "\n")
      )
    ))[1] -
      1
  )])
  eval(parse(text = code[!grepl("library\\(inlavanilla\\)", code)]))
  start <- c(
    logit_prev = qlogis(0.2),
    setNames(rep(qlogis(0.7), 5), paste0("logit_Se_", tests)),
    setNames(rep(qlogis(0.7), 5), paste0("logit_Sp_", tests))
  )
  start_fm <- c(
    start[1],
    logit_eta1 = qlogis(0.2),
    logit_eta0 = qlogis(0.1),
    setNames(rep(qlogis(0.7), 5), paste0("logit_w1_", tests)),
    setNames(rep(qlogis(0.7), 5), paste0("logit_w0_", tests))
  )
  dts <- list(
    LC = list(f = function(t) lp_lc(t, y, count), init = start),
    LCRE = list(
      f = function(t) lp_lcre(t, y, count, gh),
      init = c(start, log_sigma1 = log(0.3), log_sigma0 = log(0.3))
    ),
    FM = list(f = function(t) lp_fm(t, y, count), init = start_fm)
  )
  rwm <- function(f, mode, Sigma, n, seed, thin = 20) {
    set.seed(seed)
    L <- t(chol(Sigma))
    m <- length(mode)
    scale <- 2.38 / sqrt(m) * 0.9
    x <- mode
    fx <- f(x)
    keep <- matrix(NA_real_, floor(n / thin), m)
    k <- 0
    for (i in seq_len(n)) {
      y1 <- x + scale * as.numeric(L %*% rnorm(m))
      fy <- f(y1)
      if (is.finite(fy) && log(runif(1)) < fy - fx) {
        x <- y1
        fx <- fy
      }
      if (i %% thin == 0) {
        k <- k + 1
        keep[k, ] <- x
      }
    }
    keep[-seq_len(floor(k / 5)), , drop = FALSE]
  }
  for (nm in names(dts)) {
    d <- dts[[nm]]
    fit <- suppressWarnings(inla_vanilla(
      d$f,
      d$init,
      vb = "off",
      silent = TRUE
    ))
    file <- file.path(cache_dir, paste0("truth-", nm, ".rds"))
    if (!file.exists(file)) {
      draws <- do.call(
        rbind,
        lapply(1:4, function(s) rwm(d$f, fit$mode, fit$cov, 400000, s))
      )
      iz <- is_logz(d$f, colMeans(draws), cov(draws) * 1.4, 150000, df = 5)
      saveRDS(list(mean = colMeans(draws), logz = iz$logz), file)
    }
    tr <- readRDS(file)
    sds0 <- sqrt(diag(fit$cov))
    model <- list(
      init = d$init,
      lp = d$f,
      gr = local({
        f0 <- d$f
        h <- 1e-4 * sds0
        function(x) fd_grad(f0, x, h)
      })
    )
    rows[[length(rows) + 1]] <- tryCatch(
      study(paste0("dt_", nm), 1, model, truth = tr$mean, logz = tr$logz),
      error = function(e) {
        message(nm, ": ", conditionMessage(e))
        NULL
      }
    )
    message(nm, " done")
  }
  saveRDS(do.call(rbind, rows), file.path(cache_dir, "why-C.rds"))
}
message("finished ", part)
