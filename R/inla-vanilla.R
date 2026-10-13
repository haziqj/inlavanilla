#' INLA-style posterior approximation for any log-posterior
#'
#' @description
#' Approximates the posterior of a parameter vector from its log density
#' alone, using the approximation pipeline of INLAvaan with no latent field
#' and nothing model-specific. The pipeline has four stages:
#'
#' 1. **Joint Gaussian (Laplace).** The posterior mode \eqn{\vartheta^*} and
#'    the inverse Hessian \eqn{\Omega} there give a Gaussian approximation
#'    \eqn{N(\vartheta^*, \Omega)} and the Laplace estimate of the log marginal
#'    likelihood.
#' 2. **VB mean correction.** The Gaussian is recentred at
#'    \eqn{\tilde\vartheta = \vartheta^* + \hat\delta}, the mean that maximises
#'    its expected log-posterior, with the covariance unchanged.
#' 3. **Skew-normal marginals.** For each parameter, the log-posterior is
#'    scanned along the conditional mode path of the Gaussian, corrected for the
#'    changing volume of the other parameters, and fitted by a skew-normal
#'    \eqn{SN(\hat\xi_j, \hat\omega_j, \hat\alpha_j)}.
#' 4. **Gaussian copula draws.** [sampling()] couples the skew-normal marginals
#'    through the correlations of \eqn{\Omega}.
#'
#' @section Writing `logpost`:
#' The approximation works on an unconstrained scale: every parameter must be
#' free to take any real value. Reparameterise anything constrained, for
#' example a standard deviation \eqn{\sigma > 0} as \eqn{\log\sigma}, a
#' probability as its logit, or a correlation as its Fisher \eqn{z}. Means and
#' regression coefficients need no change.
#'
#' Return `-Inf` wherever the density is zero or cannot be evaluated, such as
#' outside the support or where a matrix is not positive definite (`NA` and
#' `NaN` count the same). The scans and the VB nodes reach several posterior
#' SDs from the mode, and an error raised there stops the fit.
#'
#' `logpost` is then the log density **of the parameters you pass**. When a
#' prior is stated on the natural scale, include the log-Jacobian of the
#' transformation: with \eqn{\vartheta = \log\sigma} and a prior
#' \eqn{p(\sigma)}, the prior term is \eqn{\log p(e^\vartheta) + \vartheta}.
#' Leaving it out targets a different posterior, not a worse approximation of
#' the same one, and nothing downstream can detect it.
#'
#' The choice of parameterisation also sets the quality of the approximation,
#' since each parameter gets its own skew-normal on the scale you choose.
#' Scales on which the posterior is close to Gaussian (log standard
#' deviations, Fisher-\eqn{z} correlations, logits) work best. Natural-scale or
#' derived quantities are available through `sampling(fit, fun = ...)`.
#'
#' The log marginal likelihood is an estimate of \eqn{\log p(y)} only when
#' `logpost` returns the full \eqn{\log p(y, \vartheta)}, including every
#' normalising constant of the likelihood and the prior. Otherwise it is off by
#' the missing constants.
#'
#' @section Accuracy:
#' The skew-normal marginals capture skewness that the Laplace approximation
#' misses, and they are accurate to the extent that a skew-normal fits each
#' marginal (see the diagnostics in [summary.inla_vanilla()]). Like INLAvaan,
#' each marginal is scanned along a straight line through the mode, and the
#' changing volume of the other parameters enters only to first order. When a
#' parameter's spread depends strongly on another, poorly determined
#' parameter, its marginal comes out too narrow, by a factor that vanishes as
#' O(1/n). Examples are a regression coefficient or a mean when the noise SD
#' is uncertain. For a normal mean with unknown variance the SD is about 11\%
#' too small at n = 12, 3\% at n = 50, and under 1\% at n = 200.
#'
#' @section Gradients:
#' `grad` is optional. Without it, derivatives come from finite differences
#' of `logpost`, and the number of calls grows quickly with the number of
#' parameters m. The mode and Hessian take about m^2 + 6m calls. The marginals
#' take about m(`n_grid` + 2m + 5). The VB correction takes about
#' (m + 1)(2m^2 + 1) calls per sweep, over a few sweeps, so it dominates
#' beyond a handful of parameters: at m = 20 a fit makes on the order of
#' 50,000 calls. With `grad`, the VB correction takes 2m^2 + 1 gradient calls
#' per sweep instead, and the fit is faster and more accurate. For more than
#' about 10 parameters, supply `grad`, or set `vb = FALSE` to skip the VB
#' correction (and its term in the log marginal likelihood).
#'
#' A supplied `grad` is checked against finite differences at `init`, and
#' its Jacobian at the mode is checked for symmetry, with a warning if either
#' check fails. Any R function works for `logpost` and `grad`, including
#' compiled code (Rcpp, TMB, RTMB). For a TMB object, which minimises, pass
#' `function(x) -obj$fn(x)` and `function(x) -as.vector(obj$gr(x))`, with
#' `init = obj$par`; its repeated names (one per element of a parameter
#' vector) are numbered, as `beta[1]`, `beta[2]`, and so on.
#'
#' @section Control:
#' `control` is a list with any of:
#' * `n_grid`: number of grid points per skew-normal scan, spread over
#'   \eqn{\pm 4} posterior SDs (default `21`).
#' * `log_threshold`: grid points whose log-posterior falls more than this
#'   below the maximum are ignored in the skew-normal fit (default `-6`).
#' * `vb_maxit`, `vb_tol`: iteration limit (default `25`) and convergence
#'   tolerance in posterior SDs (default `1e-3`) for the VB correction.
#' * `nlminb`: a list of [stats::nlminb()] controls for the mode search.
#'   `iter.max` and `eval.max` default to `1000` and `2000`.
#'
#' @param logpost A function `logpost(theta, ...)` returning the log-posterior
#'   density (up to a constant) at the named numeric vector `theta`, on the
#'   unconstrained scale. See the section on writing `logpost`.
#' @param init A named numeric vector of starting values. Its length sets the
#'   number of parameters and its names label them; unnamed values are called
#'   `theta[1]`, `theta[2]`, ..., and a repeated name is numbered in the same
#'   way.
#' @param grad Optional function `grad(theta, ...)` returning the gradient of
#'   `logpost` as a numeric vector. When `NULL` (the default), finite
#'   differences are used.
#' @param ... Further arguments passed to `logpost` and `grad`, such as the
#'   data. As in [stats::optim()] and [stats::nlminb()], they come after the
#'   gradient, and the arguments after them (`vb`, `verbose`, `control`) must
#'   be given by their full names, so data arguments cannot partially match
#'   them.
#' @param vb Logical. Apply the VB mean correction? Default `TRUE`. It enters
#'   the log marginal likelihood, and the skew-normal means are checked against
#'   it (see [summary.inla_vanilla()]).
#' @param verbose Logical. Show progress? Default `TRUE`.
#' @param control A list of tuning parameters; see the section on control.
#'
#' @returns An object of class `inla_vanilla`, with methods for [print()],
#'   [summary()], [plot()] and [sampling()]. It is a list holding, among
#'   other things:
#'   * `summary`: the posterior summary table (mean, SD, quantiles, mode).
#'   * `mode`, `cov`: the Laplace mode and covariance.
#'   * `mean_vb`: the VB-corrected mean (the mode when `vb = FALSE`).
#'   * `sn`: the skew-normal parameters `xi`, `omega` and `alpha` per
#'     parameter.
#'   * `logml`: the log marginal likelihood.
#'
#' @references
#' Jamil, H. and Rue, H. (2026). Approximate Bayesian inference for structural
#' equation models using integrated nested Laplace approximations.
#' \doi{10.48550/arXiv.2603.25690}
#'
#' Jamil, H. and Rue, H. (2026). Implementation and workflows for INLA-based
#' approximate Bayesian structural equation modelling.
#' \doi{10.48550/arXiv.2604.00671}
#'
#' @seealso [sampling()] for posterior draws and derived quantities,
#'   [summary.inla_vanilla()] for diagnostics.
#'
#' @examples
#' # Normal data with unknown mean and standard deviation, on the scale
#' # (mu, log sigma)
#' set.seed(1)
#' y <- rnorm(30, mean = 5, sd = 2)
#'
#' logpost <- function(theta, y) {
#'   mu <- theta[["mu"]]
#'   sigma <- exp(theta[["log_sigma"]])
#'   sum(dnorm(y, mu, sigma, log = TRUE)) +
#'     dnorm(mu, 0, 10, log = TRUE) + # prior on mu
#'     dexp(sigma, rate = 1, log = TRUE) + # prior on sigma ...
#'     theta[["log_sigma"]] # ... and its Jacobian
#' }
#'
#' fit <- inla_vanilla(logpost, init = c(mu = 0, log_sigma = 0), y = y)
#' fit
#'
#' # The posterior of sigma itself, through the draws
#' sigma <- sampling(fit, fun = function(theta) exp(theta[["log_sigma"]]))
#' quantile(sigma, c(0.025, 0.5, 0.975))
#' @export
inla_vanilla <- function(
  logpost,
  init,
  grad = NULL,
  ...,
  vb = TRUE,
  verbose = TRUE,
  control = list()
) {
  t_start <- elapsed()
  call <- match.call()

  check_function(logpost, "logpost")
  if (!is.null(grad)) {
    check_function(grad, "grad")
  }
  init <- check_init(init, verbose)
  check_flag(vb, "vb")
  check_flag(verbose, "verbose")
  ctrl <- vanilla_control(control)
  m <- length(init)
  parnames <- names(init)

  prog <- progress_new(verbose)

  # The engine calls these wrappers, which name the parameters, pass on `...`
  # and count the calls.
  n_lp <- 0L
  n_gr <- 0L
  last_x <- init
  # The user's functions with `...` bound once here: forwarding `...` to
  # internal helpers would let a data argument such as `g` partially match
  # one of their formals.
  logpost_x <- function(x) logpost(x, ...)
  grad_x <- if (is.null(grad)) NULL else function(x) grad(x, ...)
  lp <- function(x) {
    names(x) <- parnames
    last_x <<- x
    n_lp <<- n_lp + 1L
    if (prog$calls && n_lp %% 16L == 0L) {
      progress_tick(prog, set = n_lp)
    }
    # [[1L]] also accepts a 1 x 1 matrix, such as a quadratic form
    v <- logpost_x(x)[[1L]]
    if (is.na(v)) -Inf else v
  }
  gr <- if (is.null(grad)) {
    NULL
  } else {
    function(x) {
      names(x) <- parnames
      last_x <<- x
      n_gr <<- n_gr + 1L
      as.numeric(grad_x(x))
    }
  }
  # An error inside logpost or grad stops the fit; say where it happened
  # and at which parameter value.
  with_context <- function(stage, expr) {
    tryCatch(expr, error = function(e) {
      if (inherits(e, "inlavanilla_error")) {
        stop(e)
      }
      at <- paste(parnames, "=", signif(last_x, 4), collapse = ", ")
      abort(
        c(
          "The fit stopped with an error while {stage}.",
          "i" = "The last parameter value tried was {at}.",
          "i" = "If {.arg logpost} cannot be evaluated somewhere, return
                 {.code -Inf} there instead of failing."
        ),
        parent = e,
        call = quote(inla_vanilla())
      )
    })
  }

  with_context(
    "evaluating logpost at init",
    check_logpost_at_init(logpost_x, init)
  )
  if (!is.null(grad)) {
    with_context(
      "evaluating grad at init",
      check_grad_at_init(grad_x, lp, init)
    )
  }

  # 1a. Joint Gaussian -------------------------------------------------------
  t0 <- elapsed()
  progress_start(prog, "Finding the posterior mode")
  lap <- with_context(
    "finding the posterior mode",
    laplace_fit(lp, gr, init, ctrl)
  )
  if (!is.null(grad) && lap$asymmetry > 1e-3) {
    warn(c(
      "{.arg grad} does not look like the gradient of {.arg logpost}.",
      "x" = "Its Jacobian at the mode is not symmetric (relative asymmetry
             {signif(lap$asymmetry, 2)}), as the Hessian of a function must
             be.",
      "i" = "A wrong gradient gives a wrong mode and wrong marginals."
    ))
  }
  # The whitened gradient at the mode (how far, in posterior SDs, a Newton
  # step would still move it), not nlminb()'s return code, decides whether to
  # warn: without a gradient nlminb() can report "false convergence" at an
  # accurate mode, and it can report success short of one.
  grad_mode <- lap$grad_mode
  mode_ok <- is.finite(grad_mode) && grad_mode <= 1e-2
  if (!mode_ok) {
    hint <- if (grepl("limit", lap$opt$message)) {
      "Raise {.code iter.max} or {.code eval.max} in {.code control$nlminb}."
    } else {
      "Try other starting values, or supply {.arg grad}."
    }
    warn(c(
      "The posterior mode may be inaccurate.",
      "x" = "The gradient there is {signif(grad_mode, 2)} posterior SDs
             (whitened); it should be close to zero.",
      "i" = "{.fn nlminb} says: {.val {lap$opt$message}}.",
      "i" = hint
    ))
  }
  t_laplace <- elapsed() - t0
  progress_done(
    prog,
    if (mode_ok) {
      paste0(
        "Found the posterior mode in ",
        hl(lap$opt$iterations),
        " iteration",
        if (lap$opt$iterations == 1L) "" else "s"
      )
    } else {
      paste0(
        "Stopped near the posterior mode (gradient ",
        hl(paste0(format(signif(grad_mode, 2)), sd_unit())),
        ")"
      )
    },
    t_laplace,
    ok = mode_ok
  )

  # Gradients for the VB correction: the user's, or one-sided differences
  # stepping a thousandth of a posterior SD. A one-sided difference is biased
  # by half a step times the curvature; the bias is removed with the Laplace
  # curvature diag(H), which is close to the curvature at every node. This
  # takes m + 1 calls per gradient rather than the 2m of central differences,
  # and the VB stage is where most calls go without a gradient. In tests the
  # shift moved by about 1e-4 posterior SDs against central differences, far
  # below the error of the cubature rule itself.
  sds <- sqrt(diag(lap$Sigma))
  gr_vb <- if (is.null(gr)) {
    h_vb <- 1e-3 * sds
    bias <- 0.5 * h_vb * diag(lap$H)
    function(x) fd_grad_forward(lp, x, h = h_vb) + bias
  } else {
    gr
  }

  # 1b. VB mean correction ---------------------------------------------------
  t0 <- elapsed()
  vbc <- NULL
  if (vb) {
    progress_start(prog, "Applying the VB mean correction")
    vbc <- with_context(
      "applying the VB mean correction",
      vb_correct(lp, gr_vb, lap$mode, lap$Sigma, lap$R, lap$L, ctrl)
    )
    vb_max <- max(abs(vbc$shift) / sds)
    progress_done(
      prog,
      paste0(
        "VB mean correction, max shift ",
        hl(paste0(format(signif(vb_max, 2)), sd_unit())),
        if (vbc$fallback) " (fallback optimiser)" else ""
      ),
      elapsed() - t0,
      ok = !vbc$fallback && vb_max <= 1
    )
  }
  t_vb <- elapsed() - t0
  shift <- if (vb) vbc$shift else stats::setNames(numeric(m), parnames)

  # 2. Skew-normal marginals -------------------------------------------------
  t0 <- elapsed()
  progress_start(prog, "Fitting skew-normal marginals", total = m)
  sn <- with_context(
    "fitting the skew-normal marginals",
    fit_marginals(lp, gr, lap, ctrl, tick = function() progress_tick(prog))
  )
  t_marginals <- elapsed() - t0
  progress_done(
    prog,
    paste0(
      "Fitted ",
      hl(m),
      " skew-normal marginal",
      if (m == 1L) "" else "s"
    ),
    t_marginals
  )

  mean_vb <- lap$mode + shift
  logml <- lap$logml + if (vb) vbc$kld else 0
  mom <- sn_moments(sn[, "xi"], sn[, "omega"], sn[, "alpha"])

  diagnostics <- list(
    optim = list(
      convergence = lap$opt$convergence,
      message = lap$opt$message,
      iterations = lap$opt$iterations,
      grad_mode = grad_mode
    ),
    vb = if (vb) {
      list(
        iterations = vbc$iterations,
        fallback = vbc$fallback,
        kld = vbc$kld
      )
    },
    marginals = data.frame(
      vb_shift = shift / sds,
      sn_vs_vb = (mom$mean - mean_vb) / sds,
      alpha = sn[, "alpha"],
      gamma1 = sn[, "gamma1"],
      nmad = sn[, "nmad"],
      tails = as.logical(sn[, "tails"]),
      bend = sn[, "bend"],
      row.names = parnames
    )
  )

  out <- new_inla_vanilla(
    sn = sn[, c("xi", "omega", "alpha"), drop = FALSE],
    mode = lap$mode,
    cov = lap$Sigma,
    mean_vb = mean_vb,
    logml = logml,
    logml_laplace = lap$logml,
    diagnostics = diagnostics,
    counts = c(logpost = n_lp, grad = n_gr),
    timing = c(
      laplace = t_laplace,
      vb = t_vb,
      marginals = t_marginals,
      total = elapsed() - t_start
    ),
    vb = vb,
    has_grad = !is.null(grad),
    call = call
  )
  out
}

new_inla_vanilla <- function(
  sn,
  mode,
  cov,
  mean_vb,
  logml,
  logml_laplace,
  diagnostics,
  counts,
  timing,
  vb,
  has_grad,
  call
) {
  structure(
    list(
      summary = sn_summary(sn),
      sn = sn,
      mode = mode,
      cov = cov,
      mean_vb = mean_vb,
      logml = logml,
      logml_laplace = logml_laplace,
      diagnostics = diagnostics,
      counts = counts,
      timing = timing,
      vb = vb,
      has_grad = has_grad,
      call = call
    ),
    class = "inla_vanilla"
  )
}

# Posterior summary table from the skew-normal marginals: moments in closed
# form, quantiles by inversion, and the mode numerically.
sn_summary <- function(sn, probs = c(0.025, 0.5, 0.975)) {
  m <- nrow(sn)
  mom <- sn_moments(sn[, "xi"], sn[, "omega"], sn[, "alpha"])
  qq <- matrix(NA_real_, m, length(probs))
  md <- numeric(m)
  for (j in seq_len(m)) {
    qq[j, ] <- qsnorm(probs, sn[j, "xi"], sn[j, "omega"], sn[j, "alpha"])
    md[j] <- sn_mode(sn[j, "xi"], sn[j, "omega"], sn[j, "alpha"])
  }
  tab <- cbind(mom$mean, mom$sd, qq, md)
  dimnames(tab) <- list(
    rownames(sn),
    c("Mean", "SD", format_probs(probs), "Mode")
  )
  tab
}

format_probs <- function(probs) {
  # as stats::quantile() names them, with no padding
  paste0(formatC(100 * probs, format = "fg", width = 1, digits = 7), "%")
}

vanilla_control <- function(control) {
  defaults <- list(
    n_grid = 21L,
    log_threshold = -6,
    vb_maxit = 25L,
    vb_tol = 1e-3,
    nlminb = list()
  )
  if (!is.list(control)) {
    abort("{.arg control} must be a list.")
  }
  nms <- names(control)
  if (length(control) > 0L && (is.null(nms) || any(is.na(nms) | nms == ""))) {
    abort(c(
      "Every {.arg control} entry must be named.",
      "i" = "Valid entries are {.val {names(defaults)}}."
    ))
  }
  unknown <- setdiff(nms, names(defaults))
  if (length(unknown) > 0L) {
    abort(c(
      "Unknown {.arg control} entr{?y/ies}: {.val {unknown}}.",
      "i" = "Valid entries are {.val {names(defaults)}}."
    ))
  }
  ctrl <- utils::modifyList(defaults, control)
  is_num <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!is_num(ctrl$n_grid) || ctrl$n_grid < 7) {
    abort("{.code control$n_grid} must be a single number of at least 7.")
  }
  ctrl$n_grid <- as.integer(round(ctrl$n_grid))
  if (!is_num(ctrl$log_threshold) || ctrl$log_threshold >= 0) {
    abort("{.code control$log_threshold} must be a single negative number.")
  }
  if (!is_num(ctrl$vb_maxit) || ctrl$vb_maxit < 1) {
    abort("{.code control$vb_maxit} must be a single number of at least 1.")
  }
  ctrl$vb_maxit <- as.integer(round(ctrl$vb_maxit))
  if (!is_num(ctrl$vb_tol) || ctrl$vb_tol <= 0) {
    abort("{.code control$vb_tol} must be a single positive number.")
  }
  if (!is.list(ctrl$nlminb)) {
    abort(
      "{.code control$nlminb} must be a list of {.fn stats::nlminb} controls."
    )
  }
  ctrl$nlminb <- utils::modifyList(
    list(iter.max = 1000L, eval.max = 2000L),
    ctrl$nlminb
  )
  ctrl
}

check_function <- function(f, arg) {
  if (!is.function(f)) {
    abort(
      "{.arg {arg}} must be a function, not {.obj_type_friendly {f}}."
    )
  }
}

check_flag <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    abort("{.arg {arg}} must be {.code TRUE} or {.code FALSE}.")
  }
}

check_init <- function(init, verbose = TRUE) {
  if (!is.numeric(init) || length(init) == 0L) {
    abort(
      "{.arg init} must be a non-empty numeric vector, not {.obj_type_friendly {init}}."
    )
  }
  if (!all(is.finite(init))) {
    abort("{.arg init} must contain only finite values.")
  }
  nms <- names(init)
  if (is.null(nms)) {
    nms <- paste0("theta[", seq_along(init), "]")
  } else if (any(is.na(nms) | nms == "")) {
    abort("{.arg init} must be fully named, or not named at all.")
  } else if (anyDuplicated(nms)) {
    # As TMB names obj$par: one repeated name per element of a vector
    dup <- unique(nms[duplicated(nms)])
    renamed <- nms %in% dup
    for (d in dup) {
      i <- which(nms == d)
      nms[i] <- paste0(d, "[", seq_along(i), "]")
    }
    if (anyDuplicated(nms)) {
      abort("{.arg init} has names that cannot be made unique.")
    }
    if (verbose) {
      message(cli::format_message(c(
        "i" = "Numbered the repeated names in {.arg init}: {.val {nms[renamed]}}."
      )))
    }
  }
  stats::setNames(as.numeric(init), nms)
}

check_logpost_at_init <- function(logpost, init) {
  v <- logpost(init)
  if (!is.numeric(v) || length(v) != 1L) {
    abort(
      c(
        "{.arg logpost} must return a single number.",
        "x" = "{.code logpost(init)} returned {.obj_type_friendly {v}} of length {length(v)}."
      )
    )
  }
  if (!is.finite(v)) {
    abort(
      c(
        "{.arg logpost} is not finite at {.arg init}.",
        "x" = "{.code logpost(init)} returned {.val {v}}.",
        "i" = "Start inside the support of the posterior, and check that every
             parameter is on an unconstrained scale."
      )
    )
  }
  invisible(v)
}

check_grad_at_init <- function(grad, lp, init) {
  g <- grad(init)
  if (!is.numeric(g) || length(g) != length(init)) {
    abort(
      c(
        "{.arg grad} must return a numeric vector of length {length(init)}.",
        "x" = "{.code grad(init)} returned {.obj_type_friendly {g}} of length {length(g)}."
      )
    )
  }
  if (!all(is.finite(g))) {
    abort(
      "{.arg grad} is not finite at {.arg init}.",
    )
  }
  # Steps from the local curvature scale, so the check holds for posterior
  # SDs far from 1; a fixed step where the curvature search fails
  scales <- fd_scales(lp, init, lp(init))
  h <- ifelse(attr(scales, "ok"), 1e-3 * scales, 1e-5 * pmax(1, abs(init)))
  g_fd <- fd_grad(lp, init, h = h)
  bad <- abs(g - g_fd) > 1e-3 * (1 + abs(g_fd))
  if (any(bad)) {
    worst <- which.max(abs(g - g_fd) / (1 + abs(g_fd)))
    warn(c(
      "{.arg grad} does not match finite differences of {.arg logpost} at
       {.arg init}.",
      "x" = "For {.field {names(init)[worst]}}: {.arg grad} gives
             {.val {signif(g[worst], 4)}}, finite differences give
             {.val {signif(g_fd[worst], 4)}}.",
      "i" = "A wrong gradient gives a wrong mode and wrong marginals."
    ))
  }
  invisible(g)
}
