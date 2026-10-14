# Diagnostics of a fit: diagnostics(), and the fit-time warning.
#
# Both read diagnostics_internal(), as INLAvaan's diagnostics() and its
# warn_fit_diagnostics() do, and the flags of summary() and print() come from
# check_flags() (R/methods.R) on the same numbers, so none of them can
# disagree. The warning checks, in INLAvaan's order, the mode shift, NMAD,
# the unscanned mass, the VB shift and the Hessian condition number, at
# INLAvaan's levels. Two checks are inlavanilla's own: whether the VB
# correction settled, after the mode shift, since the two together are what
# print() calls converged, and the heavy tails, after the unscanned mass.
# INLAvaan's bullet for a non-zero return code of nlminb() is left out:
# without a gradient nlminb() reports false convergence at accurate modes.
# The VB fallback is recorded only, as in INLAvaan: it is a route of the
# solver, and the settling check judges where it lands. summary() flags every warned check at the same level, and also
# the ridge bend and the mean gap, which never warn. A check that did
# not run (VB with vb = FALSE, the scans with Gaussian marginals) is NA and
# never flagged; one that ran and gave a value that is not finite is.
warn_levels <- list(
  mode_shift = 0.1,
  nmad = 0.1,
  unscanned_mass = 0.05,
  vb_shift = 1,
  hess_cond = 1e8,
  # inlavanilla's own: the VB correction settles as the mode search does
  vb_step = 0.1
)
# The ridge bend level is a marginal SD about 10% too small,
# 1 - sqrt(1 - 3 b / 4), in the normal mean with an unknown SD
# (lab/bend-threshold.md).
summary_levels <- list(ridge_bend = 0.25, mean_gap = 0.1)

# TRUE where `x` is above `level` or not finite.
off_level <- function(x, level) {
  !(is.finite(x) & x <= level)
}

# The global and per-parameter diagnostics, from the parts of a fit:
# fit$diagnostics$optim, $vb (NULL when vb = FALSE) and $marginals, and the
# VB mean (`pars`, for grad_inf_rel and grad_rel; NA when not given). It runs
# before the fit object exists, for the warning. The definitions are
# INLAvaan's (R/method-diagnostics.R), with `grad` the gradient of minus the
# log-posterior; the unscanned mass, the VB shifts and the KL divergences are
# computed in vanilla_fit(). The per-parameter frame has no row names, as
# INLAvaan's. `message` is nlminb()'s, for the warning.
diagnostics_internal <- function(optim, vb, marginals, pars = NULL) {
  d <- marginals
  m <- nrow(d)
  g <- unname(optim$grad)
  if (is.null(pars)) {
    pars <- rep(NA_real_, m)
  }
  pars <- unname(pars)
  vb_applied <- !is.null(vb)
  max_na <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
  mean_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

  global <- c(
    npar = m,
    converged = as.numeric(optim$convergence == 0L),
    iterations = as.integer(optim$iterations),
    grad_inf = max(abs(g)),
    grad_inf_rel = max(abs(g) / (abs(pars) + 1e-06)),
    grad_l2 = sqrt(sum(g^2)),
    mode_shift_max = max(d$mode_shift_sigma),
    hess_cond = optim$hess_cond,
    vb_applied = as.numeric(vb_applied),
    vb_shift_max = max_na(abs(d$vb_shift_sigma)),
    vb_kld_global = if (vb_applied) vb$kld_global else NA_real_,
    kld_max = max_na(d$kld),
    kld_mean = mean_na(d$kld),
    nmad_max = max_na(d$nmad),
    nmad_mean = mean_na(d$nmad),
    unscanned_mass_max = max_na(d$unscanned_mass),
    # inlavanilla's own
    vb_fallback = if (vb_applied) as.numeric(vb$fallback) else NA_real_,
    vb_step_max = max_na(d$vb_step_sigma),
    mean_gap_max = max_na(abs(d$mean_gap)),
    ridge_bend_max = max_na(d$ridge_bend),
    n_heavy_tails = if (all(is.na(d$tails))) {
      NA_real_
    } else {
      sum(!d$tails, na.rm = TRUE)
    }
  )
  class(global) <- c("inla_vanilla_diagnostics", "numeric")

  param <- data.frame(
    names = rownames(d),
    grad = g,
    grad_abs = abs(g),
    grad_rel = abs(g) / (abs(pars) + 1e-06),
    mode_shift_sigma = d$mode_shift_sigma,
    kld = d$kld,
    vb_shift = d$vb_shift,
    vb_shift_sigma = d$vb_shift_sigma,
    nmad = d$nmad,
    alpha = d$alpha,
    unscanned_mass = d$unscanned_mass,
    # inlavanilla's own
    vb_step_sigma = d$vb_step_sigma,
    mean_gap = d$mean_gap,
    gamma1 = d$gamma1,
    tails = d$tails,
    ridge_bend = d$ridge_bend,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  class(param) <- c("inla_vanilla_diagnostics_param", "data.frame")

  list(global = global, param = param, message = optim$message)
}

# The order in which to list parameters: the largest values first, and any
# that are not finite before them.
worst_first <- function(x) {
  order(ifelse(is.finite(x), x, Inf), decreasing = TRUE)
}

# A cli template listing parameters as INLAvaan's warning does,
# "`name` (value), ..." for at most `n` of them, then " and N other(s)".
# The names, and the values formatted to two decimals, are put in `env` as
# `<key>_names` and `<key>_values`, and the template refers to them by
# index. So a name is never read as cli markup (names come from `init` and
# can hold braces), and each list keeps variables of its own: the templates
# are filled in only when the message is shown.
list_template <- function(env, key, names, values = NULL, n = 3L) {
  i <- seq_len(min(n, length(names)))
  names_var <- paste0(key, "_names")
  assign(names_var, names[i], envir = env)
  items <- paste0("{.code {", names_var, "[", i, "]}}")
  if (!is.null(values)) {
    values_var <- paste0(key, "_values")
    assign(
      values_var,
      formatC(values[i], digits = 2, format = "f"),
      envir = env
    )
    items <- paste0(items, " ({", values_var, "[", i, "]})")
  }
  more <- length(names) - length(i)
  paste0(
    paste(items, collapse = ", "),
    if (more > 0L) paste0(" and ", more, " other", if (more > 1L) "s")
  )
}

# The parameters flagged by a per-parameter check, worst first, as a list
# template (see list_template()).
flagged_list <- function(env, key, dg, column, level, n) {
  x <- dg$param[[column]]
  if (column %in% c("vb_shift_sigma", "mean_gap")) {
    x <- abs(x)
  }
  bad <- which(off_level(x, level))
  bad <- bad[worst_first(x[bad])]
  list_template(env, key, dg$param$names[bad], x[bad], n = n)
}

# A single value, as INLAvaan's warning formats it.
fmt_g3 <- function(x) {
  trimws(formatC(x, digits = 3, format = "g"))
}

# The sentence for a flagged check, as a cli template whose values are put
# in `env`: INLAvaan's words for the checks it shares, with the SD unit of
# sd_unit(), which summary() also prints with every parameter listed
# (`n = Inf`), and inlavanilla's own for the others. The mode and VB
# bullets name the worst parameter only, as INLAvaan's do; summary() has its
# own lines for those two.
problem_text <- function(check, dg, env, n = 3L) {
  p <- dg$param
  glob <- dg$global
  switch(
    check,
    mode = {
      env$mode_worst <- p$names[worst_first(p$mode_shift_sigma)[1L]]
      env$mode_grad <- fmt_g3(glob[["grad_inf"]])
      env$mode_shift <- fmt_g3(glob[["mode_shift_max"]])
      paste0(
        "The gradient at the posterior mode is not zero (max |grad| = ",
        "{mode_grad}): a Newton step would move {.code {mode_worst}} by ",
        "{mode_shift}",
        sd_unit(),
        "."
      )
    },
    vb_step = {
      env$vbs_worst <- p$names[worst_first(p$vb_step_sigma)[1L]]
      env$vbs_step <- fmt_g3(glob[["vb_step_max"]])
      paste0(
        "The VB correction did not settle: one more step would move ",
        "{.code {vbs_worst}} by {vbs_step}",
        sd_unit(),
        "."
      )
    },
    nmad = paste0(
      "The fitted marginal deviates from the scanned posterior (NMAD > ",
      warn_levels$nmad,
      ") for ",
      flagged_list(env, "nmad", dg, "nmad", warn_levels$nmad, n),
      "."
    ),
    unscanned_mass = paste0(
      "The fitted marginal puts more than ",
      warn_levels$unscanned_mass,
      " of its mass beyond the scanned window (4",
      sd_unit(),
      " either side of the mode) for ",
      flagged_list(
        env,
        "sem",
        dg,
        "unscanned_mass",
        warn_levels$unscanned_mass,
        n
      ),
      "; its credible limits rely on extrapolation."
    ),
    tails = paste0(
      "Heavy tails for ",
      list_template(env, "tails", p$names[!p$tails %in% TRUE], n = n),
      ": the density is still over 5% of its peak at the ends of the ",
      "scanned window (4",
      sd_unit(),
      " either side of the mode), and the skew-normal extrapolates beyond it."
    ),
    vb = {
      shift <- abs(p$vb_shift_sigma)
      env$vb_worst <- p$names[worst_first(shift)[1L]]
      env$vb_shift <- fmt_g3(glob[["vb_shift_max"]])
      paste0(
        "The VB correction shifted {.code {vb_worst}} by {vb_shift}",
        sd_unit(),
        "; the Gaussian approximation at the mode may be inaccurate."
      )
    },
    hess_cond = {
      env$hess_cond <- fmt_g3(glob[["hess_cond"]])
      paste0(
        "The Hessian at the mode is near-singular (condition number = ",
        "{hess_cond}); some parameters may be weakly identified."
      )
    },
    ridge_bend = paste0(
      "Ridge bend > ",
      summary_levels$ridge_bend,
      " for ",
      flagged_list(
        env,
        "ridge_bend",
        dg,
        "ridge_bend",
        summary_levels$ridge_bend,
        n
      ),
      ": the marginal SD may be >10% too small."
    ),
    mean_gap = paste0(
      "Mean gap > ",
      with_sd(summary_levels$mean_gap),
      " for ",
      flagged_list(env, "gap", dg, "mean_gap", summary_levels$mean_gap, n),
      ": the posterior is far from Gaussian there, and either mean may be off."
    )
  )
}

# The bullets of the fit-time warning, in INLAvaan's order with tails after
# the unscanned mass, as cli templates with the environment that fills them
# in attribute "envir". A flagged mode adds nlminb()'s message and a hint.
diagnostic_problems <- function(dg, gaussian) {
  flags <- check_flags(dg, gaussian)
  env <- new.env(parent = baseenv())
  out <- character()
  warned <- c(
    "mode",
    "vb_step",
    "nmad",
    "unscanned_mass",
    "tails",
    "vb",
    "hess_cond"
  )
  for (check in warned[flags[warned]]) {
    out <- c(out, x = problem_text(check, dg, env))
    if (check == "vb_step") {
      out <- c(
        out,
        i = paste0(
          "The skew-normal marginals do not depend on it; ",
          "{.code vb = FALSE} leaves it out."
        )
      )
    }
    if (check == "mode") {
      env$mode_message <- dg$message
      out <- c(
        out,
        i = "{.fn nlminb} says: {.val {mode_message}}.",
        i = if (isTRUE(grepl("limit", dg$message))) {
          "Raise {.code iter.max} or {.code eval.max} in {.code control$nlminb}."
        } else {
          "Try other starting values, or supply {.arg grad}."
        }
      )
    }
  }
  attr(out, "envir") <- env
  out
}

# One warning, of class inlavanilla_diagnostics_warning, when any check is
# off; FALSE (invisibly) when none is.
warn_diagnostics <- function(dg, gaussian) {
  problems <- diagnostic_problems(dg, gaussian)
  if (length(problems) == 0L) {
    return(invisible(FALSE))
  }
  env <- attr(problems, "envir")
  env$n_issues <- sum(names(problems) == "x")
  warn(
    c(
      "Fit diagnostics flagged {n_issues} potential issue{?s}:",
      problems,
      "i" = paste0(
        "Inspect with {.code diagnostics(fit)} and ",
        "{.code diagnostics(fit, type = \"param\")}; ",
        "{.code vignette(\"diagnostics\", package = \"inlavanilla\")} ",
        "explains each check."
      )
    ),
    class = "inlavanilla_diagnostics_warning",
    .envir = env
  )
  invisible(TRUE)
}

#' Diagnostics of an INLAvanilla fit
#'
#' @description
#' The checks of each stage of the approximation, as numbers: a summary of
#' the whole fit, or a table with one row per parameter. They are the checks
#' of INLAvaan's `diagnostics()`, under the same names and in the same order,
#' followed by inlavanilla's own. The fit-time warning and the flags of
#' [summary.inla_vanilla()] are computed from the same numbers, and
#' [summary.inla_vanilla()] explains each check and the level at which it is
#' flagged.
#'
#' Jamil and Rue (2026) use the VB shift and NMAD as fit diagnostics
#' (Sections 3.1 and 4.4), and the gap between the skew-normal and VB means
#' (`mean_gap`) as an internal check (Section 5.1). The mode shift, the
#' Hessian condition number and the unscanned mass come from INLAvaan's
#' `diagnostics()`, not from the paper. `tails` and `ridge_bend` are
#' inlavanilla's own.
#'
#' @details
#' A posterior SD here is the Laplace SD at the mode,
#' \eqn{s_j = \sqrt{\Sigma_{jj}}}{s_j = sqrt(Sigma_jj)}, as in INLAvaan; it
#' is not the SD column of the summary table, which is the skew-normal's.
#' Values that a fit did not compute are `NA`: the VB quantities with
#' `vb = FALSE`, and the scan quantities with `marginals = "gaussian"`.
#'
#' **Global diagnostics** (`type = "global"`), a named numeric vector:
#' * `npar`: the number of parameters.
#' * `converged`: 1 if [stats::nlminb()] reported convergence (return code
#'   0), 0 otherwise, as INLAvaan defines it. No check uses it: without a
#'   gradient `nlminb()` can report false convergence at an accurate mode.
#'   `mode_shift_max` checks the mode instead.
#' * `iterations`: the iterations of `nlminb()`, in the run that was kept if
#'   it was restarted.
#' * `grad_inf`, `grad_inf_rel`, `grad_l2`: the largest absolute value, the
#'   largest relative value and the Euclidean norm of `grad` (see below).
#'   They depend on the scale of each parameter.
#' * `mode_shift_max`: the largest `mode_shift_sigma`. Unlike the gradient
#'   norms it is scale-free. It should be close to 0.
#' * `hess_cond`: the condition number of the Hessian at the mode, with each
#'   parameter in its own posterior-SD units: the ratio of the largest to
#'   the smallest eigenvalue of the posterior correlation matrix. Large
#'   values mean that some combination of the parameters is weakly
#'   identified. INLAvaan takes it on the raw scale of the parameters, where
#'   a mere change of units makes it large, so the two agree only when the
#'   posterior SDs are similar.
#' * `vb_applied`: 1 if the VB correction ran, 0 otherwise.
#' * `vb_shift_max`: the largest absolute `vb_shift_sigma`.
#' * `vb_kld_global`: the Kullback-Leibler divergence between the Laplace
#'   Gaussian and its VB-shifted copy,
#'   \eqn{\frac12 \delta^\top H \delta}{delta' H delta / 2}, which is also
#'   the VB term of the log marginal likelihood.
#' * `kld_max`, `kld_mean`: the largest and the mean `kld`.
#' * `nmad_max`, `nmad_mean`: the largest and the mean `nmad`.
#' * `unscanned_mass_max`: the largest `unscanned_mass`.
#'
#' and then inlavanilla's own:
#' * `vb_fallback`: 1 if the VB fixed-point iteration handed over to the
#'   fallback optimiser, 0 otherwise. Like INLAvaan, inlavanilla records it
#'   but flags nothing: the fallback is a route of the solver, not a measure
#'   of the approximation.
#' * `vb_step_max`: the largest `vb_step_sigma`. Like `mode_shift_max` it
#'   should be close to 0.
#' * `mean_gap_max`: the largest absolute `mean_gap`.
#' * `ridge_bend_max`: the largest `ridge_bend`.
#' * `n_heavy_tails`: the number of parameters whose `tails` is `FALSE`.
#'
#' **Per-parameter diagnostics** (`type = "param"`), a data frame:
#' * `names`: the parameter.
#' * `grad`: the gradient of minus the log-posterior at the mode, as INLAvaan
#'   reports it: from `grad` when the fit has one, and otherwise from central
#'   differences of `logpost` a thousandth of a posterior SD apart. It
#'   should be close to 0.
#' * `grad_abs`: its absolute value.
#' * `grad_rel`: `grad_abs / (abs(par) + 1e-6)`, with `par` the VB mean (the
#'   mode when `vb = FALSE`), as in INLAvaan.
#' * `mode_shift_sigma`: the Newton step at the mode, in posterior SDs,
#'   \eqn{|(\Sigma g)_j| / s_j} with \eqn{g} the gradient. It should be
#'   close to 0.
#' * `kld`: the Kullback-Leibler divergence between the Laplace Gaussian
#'   marginal and its VB-shifted copy,
#'   \eqn{\delta_j^2 / (2\Sigma_{jj})}{delta_j^2 / (2 Sigma_jj)} =
#'   `vb_shift_sigma^2 / 2`. It is a function of the VB shift alone, and
#'   says nothing about the skew-normal fit.
#' * `vb_shift`: the VB shift of the mean, \eqn{\delta_j}, the VB mean minus
#'   the mode, on the scale of `logpost`.
#' * `vb_shift_sigma`: the same in posterior SDs.
#' * `nmad`: the NMAD of the skew-normal fit.
#' * `alpha`: the shape of the fitted skew-normal marginal.
#' * `unscanned_mass`: the mass the fitted skew-normal puts beyond the
#'   scanned window.
#'
#' and then inlavanilla's own:
#' * `vb_step_sigma`: the step still to go at the VB mean, in posterior SDs:
#'   the fixed-point step there, the Newton step of the VB stationarity
#'   condition, as `mode_shift_sigma` is for the mode. After the fallback
#'   optimiser it is taken where the optimiser landed. On the fixed-point
#'   path it is the last step taken, which bounds the next, unless that step
#'   was 0.1 or more (a `vb_tol` above 0.1), when it is taken afresh. It
#'   should be close to 0 (`Inf` where the step is not finite).
#' * `mean_gap`, `gamma1`, `tails` and `ridge_bend`, which
#'   [summary.inla_vanilla()] describes.
#'
#' Some of INLAvaan's diagnostics are left out. `nsamp` counts posterior
#' draws, which [sampling()] makes on demand rather than during the fit.
#' `grad_num` and `grad_diff` compare an analytic gradient with finite
#' differences at the mode, which would take further calls; a supplied
#' `grad` is checked at `init` and through the symmetry of its Jacobian at
#' the mode instead, with a warning if either check fails. `vb_mcse_max`,
#' `vb_mcse_mean` and `vb_mcse_sigma` estimate the quadrature error of
#' INLAvaan's Sobol points, and the cubature rule of inlavanilla is
#' deterministic (INLAvaan gives `NA` for its own deterministic rule).
#'
#' INLAvaan also exports a `diagnostics()` generic, an S4 one. When both
#' packages are attached, the one attached last masks the other. With
#' inlavanilla attached last, its `diagnostics()` passes an INLAvaan fit on
#' to INLAvaan's. With INLAvaan attached last, call
#' `inlavanilla::diagnostics(fit)` for an inlavanilla fit.
#'
#' @param object An object of class `inla_vanilla`.
#' @param type `"global"` (the default) for the diagnostics of the whole fit,
#'   or `"param"` for those of each parameter.
#' @param ... Unused.
#'
#' @returns For `type = "global"`, a named numeric vector of class
#'   `inla_vanilla_diagnostics`. For `type = "param"`, a data frame of class
#'   `inla_vanilla_diagnostics_param`. Their `print()` methods use
#'   INLAvaan's formats.
#'
#' @references
#' Jamil, H., & Rue, H. (2026). *Approximate Bayesian inference for
#' structural equation models using integrated nested Laplace
#' approximations* (arXiv:2603.25690 \[stat.ME\]). arXiv.
#' <https://doi.org/10.48550/arXiv.2603.25690>
#'
#' @seealso [summary.inla_vanilla()] for what each check means and when it is
#'   flagged, and [plot.inla_vanilla()] with `type = "sn_fit"` for the
#'   skew-normal fits themselves.
#'
#' @examples
#' set.seed(1)
#' y <- rnorm(30, mean = 5, sd = 2)
#' logpost <- function(theta, y) {
#'   sum(dnorm(y, theta[["mu"]], exp(theta[["log_sigma"]]), log = TRUE)) +
#'     dnorm(theta[["mu"]], 0, 10, log = TRUE) +
#'     dnorm(theta[["log_sigma"]], 0, 2, log = TRUE)
#' }
#' fit <- inla_vanilla(logpost, c(mu = 0, log_sigma = 0), y = y, silent = TRUE)
#'
#' # The whole fit
#' diagnostics(fit)
#'
#' # One row per parameter
#' diagnostics(fit, type = "param")
#' @export
diagnostics <- function(object, ...) {
  UseMethod("diagnostics")
}

#' @rdname diagnostics
#' @export
diagnostics.default <- function(object, ...) {
  # INLAvaan's diagnostics() is an S4 generic: pass its fits on to it when it
  # is loaded, so that attaching inlavanilla after INLAvaan does not break it
  inlavaan <- if (isS4(object)) inlavaan_diagnostics()
  if (!is.null(inlavaan)) {
    return(inlavaan(object, ...))
  }
  abort(
    "There is no {.fn diagnostics} method for {.obj_type_friendly {object}}."
  )
}

# INLAvaan's diagnostics() generic when INLAvaan is loaded, or NULL.
inlavaan_diagnostics <- function() {
  if (isNamespaceLoaded("INLAvaan")) {
    getExportedValue("INLAvaan", "diagnostics")
  }
}

#' @rdname diagnostics
#' @export
diagnostics.inla_vanilla <- function(object, type = c("global", "param"), ...) {
  type <- check_choice(type, c("global", "param"), "type")
  dg <- diagnostics_internal(
    object$diagnostics$optim,
    object$diagnostics$vb,
    object$diagnostics$marginals,
    pars = object$mean_vb
  )
  dg[[type]]
}

# Printed as INLAvaan prints its diagnostics: counts as integers, the
# gradient norms and the quantities that span orders of magnitude in
# scientific notation, and the rest to four decimals.
#' @export
print.inla_vanilla_diagnostics <- function(x, ...) {
  nm <- names(x)
  int_names <- c(
    "npar",
    "converged",
    "iterations",
    "vb_applied",
    "vb_fallback",
    "n_heavy_tails"
  )
  sci_names <- c(
    "hess_cond",
    "mode_shift_max",
    "unscanned_mass_max",
    "vb_step_max"
  )
  formatted <- vapply(
    seq_along(x),
    function(i) {
      val <- unclass(x)[[i]]
      name <- nm[i]
      if (is.na(val)) {
        "NA"
      } else if (name %in% int_names) {
        as.character(as.integer(round(val)))
      } else if (startsWith(name, "grad_") || name %in% sci_names) {
        formatC(val, digits = 2, format = "e")
      } else {
        formatC(val, digits = 4, format = "f", drop0trailing = FALSE)
      }
    },
    character(1)
  )
  names(formatted) <- nm
  print(formatted, quote = FALSE, right = TRUE)
  invisible(x)
}

# A healthy unscanned mass is below the last of four decimals, so it keeps
# three significant figures instead, as in INLAvaan.
#' @export
print.inla_vanilla_diagnostics_param <- function(x, digits = 4, ...) {
  out <- x
  for (nm in names(out)[vapply(out, is.numeric, logical(1))]) {
    out[[nm]] <- if (nm == "unscanned_mass") {
      signif(out[[nm]], 3)
    } else {
      round(out[[nm]], digits)
    }
  }
  print.data.frame(out, ...)
  invisible(x)
}
