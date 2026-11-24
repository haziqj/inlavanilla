################################################################################
#
# The Sobol VB mean correction of INLAvaan 0.3.5, as research code
#
# INLAvaan 0.4.0 takes its VB correction from inlavanilla and drops its own
# Sobol rule. This is that rule, kept outside any package as a baseline for a
# study of the one-step VB correction: the same mode, Hessian and whitening
# factor, with only the rule changed. It is a copy of R/sobol.R and of the VB
# block of inlavaan() in INLAvaan 0.3.5 (tag v0.3.5, R/inlavaan.R, the Sobol
# branch), with the fit's closures as arguments: the fixed-point iteration with
# Anderson(1) acceleration, its nlminb() fallback, and the split-half Monte
# Carlo standard error. sobol-owen.rds is the table SobolOwen of that
# version's R/sysdata.rda (128 Owen-scrambled Sobol points in 300
# dimensions). check-v0.3.5.R checks that it reproduces the shift of
# INLAvaan 0.3.5 from the same inputs.
#
################################################################################

## ----- Configuration ---------------------------------------------------------

sobol_dir <- if (exists("sobol_dir")) sobol_dir else "lab/sobol-vb"

## ----- The rule --------------------------------------------------------------

sobol_owen_table <- function() {
  readRDS(file.path(sobol_dir, "sobol-owen.rds"))
}

# `n` scrambled Sobol points mapped through the Gaussian quantile function,
# centred, rotated by `L` (L L' = Sigma), with the mode itself prepended. The
# centring makes the node mean exactly zero, which the expectation over a
# mean-zero Gaussian needs.
sobol_vb_nodes <- function(n, L, table = sobol_owen_table()) {
  d <- ncol(L)
  if (n > nrow(table) || d > ncol(table)) {
    stop("The stored table has ", nrow(table), " points in ", ncol(table),
         " dimensions.")
  }
  z_std <- stats::qnorm(table[seq_len(n), seq_len(d), drop = FALSE])
  z_std <- sweep(z_std, 2, colMeans(z_std))
  rbind(0, z_std %*% t(L))
}

# The VB shift by INLAvaan 0.3.5's Sobol rule. `lp` and `gr` are the
# log-posterior and its gradient, `mode` the mode, `Sigma` the Laplace
# covariance, and `L` a factor with L L' = Sigma: INLAvaan 0.3.5 took the
# Cholesky factor in the canonical order of the names, and inlavanilla's
# cov_root() gives one that does not depend on the order. `fixed` holds the
# parameters whose shift is held at zero (INLAvaan's saturated means). Returns
# the shift, its Monte Carlo standard error from the two half-sets of nodes,
# the iterations, whether nlminb() took over, and the VB term of the log
# marginal likelihood, 0.5 shift' Sigma^-1 shift.
sobol_vb <- function(
  lp,
  gr,
  mode,
  Sigma,
  L,
  n_qmc = 64L,
  maxit = 25L,
  tol = 1e-3,
  maxstep = 1,
  fixed = integer(),
  table = sobol_owen_table()
) {
  m <- length(mode)
  zs <- sobol_vb_nodes(n_qmc, L, table)
  sds <- sqrt(diag(Sigma))

  ob_shift <- function(shift) {
    mu <- mode + shift
    total <- 0
    for (b in seq_len(nrow(zs))) {
      total <- total + lp(mu + zs[b, , drop = TRUE])
    }
    -total / nrow(zs)
  }
  sweep_nodes <- function(shift) {
    mu <- mode + shift
    ns <- nrow(zs)
    nhalf <- floor(ns / 2)
    gA <- gB <- numeric(m)
    for (b in seq_len(ns)) {
      g <- gr(mu + zs[b, , drop = TRUE])
      if (b <= nhalf) {
        gA <- gA + g
      } else {
        gB <- gB + g
      }
    }
    list(score = (gA + gB) / ns, gA = gA / nhalf, gB = gB / (ns - nhalf))
  }
  gA <- gB <- numeric(m)
  ob <- function(delta) ob_shift(as.numeric(L %*% delta))
  ob_gr <- function(delta) {
    sw <- sweep_nodes(as.numeric(L %*% delta))
    gA <<- sw$gA
    gB <<- sw$gB
    as.numeric(t(L) %*% (-1 * sw$score))
  }

  # The fixed-point iteration with Anderson(1) acceleration, which falls back
  # to nlminb() on a step over `maxstep` SDs, a score that is not finite, or no
  # convergence in `maxit` steps
  shift <- numeric(m)
  iter <- 0L
  move <- Inf
  fallback <- FALSE
  step_prev <- shift_prev <- NULL
  for (it in seq_len(maxit)) {
    sw <- sweep_nodes(shift)
    if (!all(is.finite(sw$score))) {
      fallback <- TRUE
      break
    }
    gA <- sw$gA
    gB <- sw$gB
    step <- as.numeric(Sigma %*% sw$score)
    step[fixed] <- 0
    move <- max(abs(step) / sds)
    if (move > maxstep) {
      fallback <- TRUE
      break
    }
    if (move < tol) {
      shift <- shift + step
      iter <- it
      break
    }
    update <- step
    if (!is.null(step_prev)) {
      df <- step - step_prev
      dx <- shift - shift_prev
      denom <- sum(df^2)
      if (denom > 0) {
        gam <- sum(step * df) / denom
        cand <- step - gam * (dx + df)
        cand[fixed] <- 0
        if (
          all(is.finite(cand)) &&
            max(abs(cand) / sds) <= min(maxstep, 2 * move)
        ) {
          update <- cand
        }
      }
    }
    step_prev <- step
    shift_prev <- shift
    shift <- shift + update
    iter <- it
  }
  if (move >= tol) {
    fallback <- TRUE
  }
  if (fallback) {
    nl <- stats::nlminb(
      start = numeric(m),
      objective = ob,
      gradient = ob_gr,
      control = list(rel.tol = 1e-8)
    )
    shift <- as.numeric(L %*% nl$par)
    shift[fixed] <- 0
    iter <- nl$iterations
  }
  mcse <- abs(as.numeric(Sigma %*% (gA - gB))) / 2
  mcse[fixed] <- 0
  names(shift) <- names(mcse) <- names(mode)
  list(
    shift = shift,
    mcse = mcse,
    iterations = iter,
    fallback = fallback,
    kld_global = 0.5 * sum(shift * solve(Sigma, shift))
  )
}
