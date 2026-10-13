#' INLA-style posterior approximation for any log-posterior
#'
#' @description
#' Approximates the posterior of a parameter vector from its log density
#' alone, using the approximation pipeline of INLAvaan (Jamil & Rue, 2026a,
#' 2026b) with no latent field and nothing model-specific. The pipeline has
#' four stages:
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
#' Stages 2 and 3 are optional and independent of each other. `vb = FALSE`
#' skips the VB correction, which leaves the marginals unchanged (they are
#' scanned around the mode) but gives the Laplace log marginal likelihood.
#' `marginals = "gaussian"` skips the scans and uses the Gaussian's own
#' marginals \eqn{N(\tilde\vartheta_j, \Omega_{jj})}, centred at the VB mean
#' (or at the mode when `vb = FALSE`); [sampling()] then draws from the joint
#' Gaussian. So `vb = FALSE, marginals = "gaussian"` is the Laplace
#' approximation alone, and `marginals = "gaussian"` alone is the Laplace
#' approximation with the VB correction.
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
#' is uncertain. For a normal mean with unknown variance the SD is about 11%
#' too small at n = 12, 3% at n = 50, and under 1% at n = 200.
#'
#' A skew-normal also cannot follow a tail heavier than a Gaussian's. In a
#' structural equation model with 39 parameters and n = 75, the data could not
#' keep two residual variances away from zero, so on the log scale their
#' posteriors had a long, exponential left tail, set by the prior. Checked
#' against MCMC, those two skew-normal marginals were less than half as wide
#' as the posterior, with the 2.5% quantile about 2 posterior SDs too high;
#' the other 37 were within 7% in SD. On the variance scale, where the tail
#' near zero is compressed, the two were much closer: means 0.2 posterior SDs
#' too high and SDs 10 to 15% too small. INLAvaan's marginals of this model
#' are the same as inlavanilla's: Jamil and Rue (2026a, Section 4.4) find
#' the same miss in the lower tail of a variance near zero, the boundary
#' pile-up that a unimodal skew-normal cannot capture. The `bend` diagnostic
#' of [summary.inla_vanilla()] flagged both, but its size does not say how
#' far off a marginal is, and the `tails` check, which looks only 4 posterior
#' SDs either side of the mode (the Laplace SDs), missed them.
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
#' correction (and its term in the log marginal likelihood): the skew-normal
#' marginals do not depend on it. `marginals = "gaussian"` saves the calls of
#' the scans.
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
#'   `iter.max` and `eval.max` default to `1000` and `2000`, and `trace = 1`
#'   prints its iterations.
#'
#' @section Output:
#' By default each stage prints one line as it finishes, with the time it
#' took, and shows a progress bar while it runs for longer than half a
#' second, or than the option `cli.progress_show_after` when it is set. In
#' knitr or Quarto documents, where a bar cannot redraw in place, a larger
#' `cli.progress_show_after` in the setup chunk keeps bars out of the output.
#' Warnings, such as a gradient that does not match `logpost`, are signalled
#' as usual. A fit whose checks look off, such as a mode that was not
#' reached, ends with one warning that lists them, as INLAvaan's does (see
#' [summary.inla_vanilla()] and [diagnostics()]).
#'
#' `verbose = TRUE` adds a full diagnostic of the fit, in grey lines under
#' each stage. It reports what the fit computes anyway, with no further calls
#' of `logpost` or `grad`:
#' * Before the mode search: the parameters, where the derivatives come from,
#'   `logpost` at `init`, and for a supplied `grad` its largest relative
#'   discrepancy from finite differences there.
#' * The mode: the message, iterations and evaluations of [stats::nlminb()],
#'   and whether it was restarted with the parameters rescaled; the
#'   log-posterior at the mode, the whitened gradient there (which decides
#'   the restart) and the largest Newton step in posterior SDs (which the
#'   fit-time warning checks); how the Hessian was taken (and, for a
#'   supplied `grad`, the asymmetry of its Jacobian); the largest posterior
#'   correlation and the smallest eigenvalue of the correlation matrix, which
#'   is near zero when some combination of the parameters is poorly
#'   identified, and the condition number of the Hessian in posterior-SD
#'   units; and the Laplace log marginal likelihood.
#' * The VB correction: the largest step, in posterior SDs, at each
#'   fixed-point iteration, whether Anderson acceleration or the fallback
#'   optimiser was used, the largest shift, and the term it adds to the log
#'   marginal likelihood.
#' * The marginals: the shape, NMAD, bend, tails, scan-end mass and volume
#'   correction `gamma1` of each parameter, with the values that
#'   [summary.inla_vanilla()] flags marked.
#' * The calls of `logpost` and `grad` in each stage, and in all, with the
#'   total time.
#'
#' For the iteration trace of [stats::nlminb()] itself, pass
#' `control = list(nlminb = list(trace = 1))`.
#'
#' `silent = TRUE` prints nothing, for simulations and other repeated fits:
#' no progress bars and no stage lines, and every warning and message raised
#' during the fit, including those of `logpost` and `grad`, is muffled.
#' Errors are still raised. Output printed with [cat()] or [print()], by
#' `logpost`, `grad` or the trace of [stats::nlminb()], is not captured.
#'
#' Whatever the level, the warnings raised during the fit are recorded in
#' `fit$warnings`, so that a simulation can check them afterwards.
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
#'   gradient, and the arguments after them (`vb`, `marginals`, `verbose`,
#'   `silent`, `control`) must be given by their full names, so data
#'   arguments cannot partially match them. The arguments before them can be
#'   matched, as in `optim()`: a data argument named `g` or `gr` would be
#'   taken as `grad`, so give it a longer name.
#' @param vb Logical. Apply the VB mean correction? Default `TRUE`. It enters
#'   the log marginal likelihood, and the skew-normal means are checked against
#'   it (see [summary.inla_vanilla()]).
#' @param marginals The marginal of each parameter: `"skewnormal"` (the
#'   default) fits a skew-normal to a scan of the log-posterior; `"gaussian"`
#'   uses the Gaussian marginal of the Laplace approximation, centred at the
#'   VB mean when `vb = TRUE`, with no scans. Abbreviations such as `"skew"`
#'   and `"gauss"` work too.
#' @param verbose Logical. Print a full diagnostic of each stage under its
#'   progress line? Default `FALSE`. See the section on output.
#' @param silent Logical. Print nothing at all, and muffle the warnings and
#'   messages raised during the fit? Default `FALSE`. The warnings are still
#'   recorded in `fit$warnings`. Cannot be combined with `verbose = TRUE`.
#' @param control A list of tuning parameters; see the section on control.
#'
#' @returns An object of class `inla_vanilla`, with methods for [print()],
#'   [summary()], [plot()], [sampling()] and [diagnostics()]. It is a list
#'   holding, among other things:
#'   * `summary`: the posterior summary table (mean, SD, quantiles, mode).
#'   * `mode`, `cov`: the Laplace mode and covariance.
#'   * `mean_vb`: the VB-corrected mean (the mode when `vb = FALSE`).
#'   * `sn`: the skew-normal parameters `xi`, `omega` and `alpha` per
#'     parameter (with `alpha = 0` for Gaussian marginals).
#'   * `vb`, `marginals`, `verbose`, `silent`: the settings the fit used.
#'   * `logml`: the log marginal likelihood, and `logml_laplace` its Laplace
#'     part.
#'   * `diagnostics`: the checks of each stage, which [diagnostics()]
#'     returns under INLAvaan's names and [summary.inla_vanilla()] explains.
#'     `optim` holds the return code, message and iterations of
#'     [stats::nlminb()], the gradient of minus the log-posterior at the mode
#'     (`grad`), the whitened gradient there (`grad_mode`), and the condition
#'     number of the Hessian in posterior-SD units (`hess_cond`). `vb` holds
#'     the iterations of the VB correction, whether it fell back to
#'     [stats::nlminb()] (`fallback`), and the Kullback-Leibler divergence of
#'     the shift (`kld_global`); it is `NULL` when `vb = FALSE`. `marginals`
#'     is a data frame with a row per parameter: `mode_shift_sigma`, `kld`,
#'     `vb_shift`, `vb_shift_sigma`, `nmad`, `alpha`, `scan_end_mass`,
#'     `sn_vs_vb`, `gamma1`, `tails` and `bend`, with `NA` for the checks the
#'     fit did not run.
#'   * `scan`: the scans the skew-normals were fitted to, for
#'     `plot(fit, type = "sn_fit")`: the grid `z` in posterior SDs from the
#'     mode, the log-posterior at each point (`logpost`, a column per
#'     parameter), and the log scale factor `logC` of each fit on that grid.
#'     It is `NULL` for Gaussian marginals.
#'   * `warnings`: the messages of the warnings raised during the fit, as a
#'     character vector of plain text, each message once (empty when there
#'     were none). They are recorded whether or not `silent` is `TRUE`.
#'
#' @references
#' Jamil, H., & Rue, H. (2026a). *Approximate Bayesian inference for
#' structural equation models using integrated nested Laplace
#' approximations* (arXiv:2603.25690 \[stat.ME\]). arXiv.
#' <https://doi.org/10.48550/arXiv.2603.25690>
#'
#' Jamil, H., & Rue, H. (2026b). *Implementation and workflows for INLA-based
#' approximate Bayesian structural equation modelling* (arXiv:2604.00671
#' \[stat.CO\]). arXiv. <https://doi.org/10.48550/arXiv.2604.00671>
#'
#' @seealso [sampling()] for posterior draws and derived quantities,
#'   [summary.inla_vanilla()] and [diagnostics()] for diagnostics.
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
#'
#' # A full diagnostic of each stage
#' fit <- inla_vanilla(
#'   logpost,
#'   init = c(mu = 0, log_sigma = 0),
#'   y = y,
#'   verbose = TRUE
#' )
#'
#' # Nothing printed, as in a simulation; any warnings are kept in the fit
#' fit <- inla_vanilla(
#'   logpost,
#'   init = c(mu = 0, log_sigma = 0),
#'   y = y,
#'   silent = TRUE
#' )
#' fit$warnings
#' @export
inla_vanilla <- function(
  logpost,
  init,
  grad = NULL,
  ...,
  vb = TRUE,
  marginals = c("skewnormal", "gaussian"),
  verbose = FALSE,
  silent = FALSE,
  control = list()
) {
  t_start <- elapsed()
  call <- match.call()

  check_function(logpost, "logpost")
  if (!is.null(grad)) {
    check_function(grad, "grad")
  }
  check_flag(vb, "vb")
  marginals <- check_choice(marginals, c("skewnormal", "gaussian"), "marginals")
  check_flag(verbose, "verbose")
  check_flag(silent, "silent")
  if (verbose && silent) {
    abort(c(
      "{.arg verbose} and {.arg silent} cannot both be {.code TRUE}.",
      "i" = "{.code verbose = TRUE} prints a full diagnostic of the fit, and
             {.code silent = TRUE} prints nothing."
    ))
  }
  init <- check_init(init, show = !silent)
  ctrl <- vanilla_control(control)

  # The user's functions with `...` bound once here: forwarding `...` to
  # internal helpers would let a data argument such as `g` partially match
  # one of their formals.
  logpost_x <- function(x) logpost(x, ...)
  grad_x <- if (is.null(grad)) NULL else function(x) grad(x, ...)

  # Every warning raised during the fit, by the package or by logpost and
  # grad, is recorded once in fit$warnings. With silent = TRUE it is then
  # muffled, and so is every message; otherwise it is signalled as usual.
  # Errors pass through either way.
  prog <- progress_new(show = !silent, detail = verbose)
  withCallingHandlers(
    vanilla_fit(
      logpost_x,
      grad_x,
      init,
      vb = vb,
      marginals = marginals,
      ctrl = ctrl,
      prog = prog,
      call = call,
      t_start = t_start
    ),
    warning = function(w) {
      msg <- warning_text(w)
      if (!msg %in% prog$warnings) {
        prog$warnings <- c(prog$warnings, msg)
      }
      if (silent) {
        # tryInvokeRestart(), as suppressWarnings() does: a condition
        # signalled with signalCondition() has no muffle restart
        tryInvokeRestart("muffleWarning")
      }
    },
    message = function(m) {
      if (silent) {
        tryInvokeRestart("muffleMessage")
      }
    }
  )
}

# The stages of the fit, given the user's functions with `...` already bound,
# the checked arguments, and the progress reporter `prog`, which also holds
# the record of warnings.
vanilla_fit <- function(
  logpost_x,
  grad_x,
  init,
  vb,
  marginals,
  ctrl,
  prog,
  call,
  t_start
) {
  m <- length(init)
  parnames <- names(init)
  has_grad <- !is.null(grad_x)

  # The engine calls these wrappers, which name the parameters, count the
  # calls and tick the progress bar. The spinner counts logpost and grad
  # calls together, so stages driven by grad tick too.
  n_lp <- 0L
  n_gr <- 0L
  last_x <- init
  tick <- function() {
    n <- n_lp + n_gr
    if (prog$calls && n %% 16L == 0L) {
      progress_tick(prog, set = n)
    }
  }
  lp <- function(x) {
    names(x) <- parnames
    last_x <<- x
    n_lp <<- n_lp + 1L
    tick()
    # [[1L]] also accepts a 1 x 1 matrix, such as a quadratic form
    v <- logpost_x(x)[[1L]]
    if (is.na(v)) -Inf else v
  }
  gr <- if (!has_grad) {
    NULL
  } else {
    function(x) {
      names(x) <- parnames
      last_x <<- x
      n_gr <<- n_gr + 1L
      tick()
      as.numeric(grad_x(x))
    }
  }
  # The calls since the last stage, for verbose output
  calls_seen <- c(0L, 0L)
  calls_since <- function() {
    now <- c(n_lp, n_gr)
    out <- now - calls_seen
    calls_seen <<- now
    out
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

  lp_init <- with_context(
    "evaluating logpost at init",
    check_logpost_at_init(logpost_x, init)
  )
  grad_check <- if (has_grad) {
    with_context(
      "evaluating grad at init",
      check_grad_at_init(grad_x, lp, init)
    )
  }
  if (prog$detail) {
    report_start(prog, init, lp_init, grad_check, calls_since())
  }

  # 1a. Joint Gaussian -------------------------------------------------------
  t0 <- elapsed()
  progress_start(prog, "Finding the posterior mode")
  lap <- with_context(
    "finding the posterior mode",
    laplace_fit(lp, gr, init, ctrl)
  )
  if (has_grad && lap$asymmetry > 1e-3) {
    warn(c(
      "{.arg grad} does not look like the gradient of {.arg logpost}.",
      "x" = "Its Jacobian at the mode is not symmetric (relative asymmetry
             {signif(lap$asymmetry, 2)}), as the Hessian of a function must
             be.",
      "i" = "A wrong gradient gives a wrong mode and wrong marginals."
    ))
  }
  # INLAvaan's Newton step at the mode, in posterior SDs, decides whether the
  # mode is marked here and reported in the fit-time warning, not nlminb()'s
  # return code: without a gradient nlminb() can report "false convergence"
  # at an accurate mode, and it can report success short of one.
  mode_shift_max <- max(lap$mode_shift)
  mode_ok <- is.finite(mode_shift_max) &&
    mode_shift_max <= warn_levels$mode_shift
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
        "Stopped near the posterior mode (Newton step ",
        hl(with_sd(mode_shift_max)),
        ")"
      )
    },
    t_laplace,
    ok = mode_ok
  )
  if (prog$detail) {
    report_mode(prog, lap, calls_since())
  }

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
    # The fallback is a route of the solver, reported but not marked, as in
    # INLAvaan
    vb_max <- max(abs(vbc$shift) / sds)
    progress_done(
      prog,
      paste0(
        "VB mean correction, max shift ",
        hl(paste0(format(signif(vb_max, 2)), sd_unit())),
        if (vbc$fallback) " (fallback optimiser)" else ""
      ),
      elapsed() - t0,
      ok = is.finite(vb_max) && vb_max <= warn_levels$vb_shift
    )
    if (prog$detail) {
      report_vb(prog, vbc, sds, lap$logml, calls_since())
    }
  } else {
    progress_info(prog, "No VB correction (vb = FALSE)")
  }
  t_vb <- elapsed() - t0
  shift <- if (vb) vbc$shift else stats::setNames(numeric(m), parnames)
  mean_vb <- lap$mode + shift

  # 2. Skew-normal marginals -------------------------------------------------
  # With Gaussian marginals there are no scans: each marginal is the
  # Gaussian's own, N(mean_vb_j, Sigma_jj), a skew-normal with shape 0.
  t0 <- elapsed()
  if (marginals == "skewnormal") {
    progress_start(prog, "Fitting skew-normal marginals", total = m)
    mar <- with_context(
      "fitting the skew-normal marginals",
      fit_marginals(lp, gr, lap, ctrl, tick = function() progress_tick(prog))
    )
    sn <- cbind(
      mar$sn,
      scan_end_mass = scan_end_mass(mar$sn, lap$mode, sds)
    )
    # The scans, kept for plot(type = "sn_fit"): no further calls
    scan <- list(z = mar$z, logpost = mar$logpost, logC = sn[, "logC"])
    progress_done(
      prog,
      paste0(
        "Fitted ",
        hl(m),
        " skew-normal marginal",
        if (m == 1L) "" else "s"
      ),
      elapsed() - t0
    )
    if (prog$detail) {
      report_marginals(prog, sn, calls_since())
    }
  } else {
    na <- rep(NA_real_, m)
    sn <- cbind(
      xi = mean_vb,
      omega = sds,
      alpha = 0,
      gamma1 = na,
      nmad = na,
      tails = na,
      bend = na,
      logC = na,
      scan_end_mass = na
    )
    rownames(sn) <- parnames
    scan <- NULL
    progress_info(
      prog,
      "Gaussian marginals (marginals = \"gaussian\"): no scans"
    )
  }
  t_marginals <- elapsed() - t0

  logml <- lap$logml + if (vb) vbc$kld_global else 0
  mom <- sn_moments(sn[, "xi"], sn[, "omega"], sn[, "alpha"])

  # The diagnostics, in the names and order of INLAvaan's diagnostics() where
  # it has the quantity. Checks that did not run are NA: the VB columns with
  # vb = FALSE, the scan columns with Gaussian marginals.
  skew <- marginals == "skewnormal"
  na <- rep(NA_real_, m)
  diagnostics <- list(
    optim = list(
      convergence = lap$opt$convergence,
      message = lap$opt$message,
      iterations = lap$opt$iterations,
      # INLAvaan's sign: the gradient of nlminb()'s objective, -logpost
      grad = -lap$grad,
      grad_mode = lap$grad_mode,
      hess_cond = lap$hess_cond
    ),
    vb = if (vb) {
      list(
        iterations = vbc$iterations,
        fallback = vbc$fallback,
        kld_global = vbc$kld_global
      )
    },
    marginals = data.frame(
      mode_shift_sigma = lap$mode_shift,
      kld = if (vb) 0.5 * (shift / sds)^2 else na,
      vb_shift = if (vb) shift else na,
      vb_shift_sigma = if (vb) shift / sds else na,
      nmad = sn[, "nmad"],
      alpha = if (skew) sn[, "alpha"] else na,
      scan_end_mass = sn[, "scan_end_mass"],
      sn_vs_vb = if (vb && skew) (mom$mean - mean_vb) / sds else na,
      gamma1 = sn[, "gamma1"],
      tails = as.logical(sn[, "tails"]),
      bend = sn[, "bend"],
      row.names = parnames
    )
  )

  # One warning for every check that looks off, before the fit is built so
  # that fit$warnings records it
  warn_diagnostics(
    diagnostics_internal(
      diagnostics$optim,
      diagnostics$vb,
      diagnostics$marginals,
      pars = mean_vb
    ),
    gaussian = !skew
  )

  t_total <- elapsed() - t_start
  if (prog$detail) {
    report_end(prog, n_lp, n_gr, t_total)
  }
  new_inla_vanilla(
    sn = sn[, c("xi", "omega", "alpha"), drop = FALSE],
    mode = lap$mode,
    cov = lap$Sigma,
    mean_vb = mean_vb,
    logml = logml,
    logml_laplace = lap$logml,
    diagnostics = diagnostics,
    scan = scan,
    counts = c(logpost = n_lp, grad = n_gr),
    timing = c(
      laplace = t_laplace,
      vb = t_vb,
      marginals = t_marginals,
      total = t_total
    ),
    vb = vb,
    marginals = marginals,
    has_grad = has_grad,
    verbose = prog$detail,
    silent = !prog$show,
    warnings = prog$warnings,
    call = call
  )
}

new_inla_vanilla <- function(
  sn,
  mode,
  cov,
  mean_vb,
  logml,
  logml_laplace,
  diagnostics,
  scan,
  counts,
  timing,
  vb,
  marginals,
  has_grad,
  verbose,
  silent,
  warnings,
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
      scan = scan,
      counts = counts,
      timing = timing,
      vb = vb,
      marginals = marginals,
      has_grad = has_grad,
      verbose = verbose,
      silent = silent,
      warnings = warnings,
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
  sprintf("%s%%", formatC(100 * probs, format = "fg", width = 1, digits = 7))
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

# One of `choices` (lower case), or an unambiguous abbreviation of one in any
# case, as with match.arg(): "skew" and "Gauss" both work. The whole vector
# (the default) means the first. The package's own error otherwise.
check_choice <- function(x, choices, arg) {
  if (identical(x, choices)) {
    return(choices[[1L]])
  }
  i <- NA_integer_
  if (is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)) {
    i <- pmatch(tolower(x), choices)
  }
  if (is.na(i)) {
    abort("{.arg {arg}} must be one of {.or {.val {choices}}}.")
  }
  choices[[i]]
}

check_init <- function(init, show = TRUE) {
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
    if (show) {
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
  # Compare per curvature scale, the change in log-posterior over about one
  # posterior SD, so a gradient that is wrong by a factor is caught whatever
  # the scale of the parameter
  u <- ifelse(attr(scales, "ok"), scales, 1)
  score <- function(g, g_fd) abs(g - g_fd) * u / (1 + abs(g_fd) * u)
  err <- score(g, g_fd)
  where <- "at {.arg init}"
  # Near the mode both gradients are close to zero, so a gradient that is
  # wrong by a factor can still pass there. Check again one curvature scale
  # away, where the gradient is about one per scale; skip the point if
  # logpost or grad cannot be evaluated there.
  x1 <- init + u
  if (is.finite(lp(x1))) {
    g1 <- tryCatch(as.numeric(grad(x1)), error = function(e) NULL)
    if (length(g1) == length(init) && all(is.finite(g1))) {
      g1_fd <- fd_grad(lp, x1, h = h)
      err1 <- score(g1, g1_fd)
      if (all(is.finite(err1)) && max(err1) > max(err)) {
        err <- err1
        g <- g1
        g_fd <- g1_fd
        where <- "one curvature scale from {.arg init}"
      }
    }
  }
  if (any(err > 1e-3)) {
    worst <- which.max(err)
    warn(c(
      paste0(
        "{.arg grad} does not match finite differences of {.arg logpost} ",
        where,
        "."
      ),
      "x" = "For {.field {names(init)[worst]}}: {.arg grad} gives
             {.val {signif(g[worst], 4)}}, finite differences give
             {.val {signif(g_fd[worst], 4)}}.",
      "i" = "A wrong gradient gives a wrong mode and wrong marginals."
    ))
  }
  # The largest discrepancy, for verbose output
  invisible(max(err))
}
