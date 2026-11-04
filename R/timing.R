#' Time taken by each stage of an INLAvanilla fit
#'
#' The elapsed (wall-clock) time of each stage of [inla_vanilla()], in
#' seconds, as INLAvaan's `timing()` gives them.
#'
#' The segments are:
#' * `init`: the checks of the arguments, with the calls of `logpost`, `grad`
#'   and `hess` at `init`.
#' * `laplace`: the posterior mode and the Hessian there.
#' * `vb`: the VB mean correction, close to zero with `vb = FALSE`.
#' * `marginals`: the scans and the skew-normal fits, close to zero with
#'   `marginals = "gaussian"`.
#' * `diagnostics`: the diagnostics of the fit and its warning.
#' * `total`: the whole fit.
#'
#' The segments are disjoint, and `total` is their sum. Each stage includes
#' its progress output, so with `verbose = TRUE` they can be a little longer
#' than the times printed during the fit, which time the work alone.
#' [sampling()] and [loo()] run after the fit, so their time is not included.
#'
#' INLAvaan also exports a `timing()` generic, an S4 one. When both packages
#' are attached, the one attached last masks the other. With inlavanilla
#' attached last, its `timing()` passes an INLAvaan fit on to INLAvaan's.
#' With INLAvaan attached last, call `inlavanilla::timing(fit)` for an
#' inlavanilla fit.
#'
#' @param object An object of class `inla_vanilla`.
#' @param what The segments to return, by name, or `"all"` for every segment.
#'   Defaults to `"total"`.
#' @param ... Unused.
#'
#' @returns A named numeric vector of class `inla_vanilla_timing`, in
#'   seconds. It prints each time as [summary.inla_vanilla()] does: in
#'   milliseconds below a second, then in seconds, then in minutes and
#'   seconds.
#'
#' @seealso [summary.inla_vanilla()], which gives the time of each stage with
#'   the number of calls.
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
#' timing(fit)
#'
#' # Every stage
#' timing(fit, what = "all")
#'
#' # Some stages
#' timing(fit, what = c("vb", "marginals"))
#' @export
timing <- function(object, ...) {
  UseMethod("timing")
}

#' @rdname timing
#' @export
timing.default <- function(object, ...) {
  # INLAvaan's timing() is an S4 generic: pass its fits on to it when it is
  # loaded, so that attaching inlavanilla after INLAvaan does not break it
  inlavaan <- if (isS4(object)) inlavaan_timing()
  if (!is.null(inlavaan)) {
    return(inlavaan(object, ...))
  }
  abort(
    "There is no {.fn timing} method for {.obj_type_friendly {object}}."
  )
}

# INLAvaan's timing() generic when INLAvaan is loaded, or NULL.
inlavaan_timing <- function() {
  if (isNamespaceLoaded("INLAvaan")) {
    getExportedValue("INLAvaan", "timing")
  }
}

#' @rdname timing
#' @export
timing.inla_vanilla <- function(object, what = "total", ...) {
  tm <- object$timing
  available <- names(tm)
  if (!is.character(what) || length(what) == 0L || anyNA(what)) {
    abort(
      "{.arg what} must name timing segments, or be {.val all}."
    )
  }
  if (identical(what, "all")) {
    what <- available
  }
  unknown <- setdiff(what, available)
  if (length(unknown) > 0L) {
    abort(c(
      "Unknown timing segment{?s}: {.val {unknown}}.",
      "i" = "Use {.val all}, or any of {.val {available}}."
    ))
  }
  out <- tm[what]
  class(out) <- c("inla_vanilla_timing", "numeric")
  out
}

#' @export
print.inla_vanilla_timing <- function(x, ...) {
  print(
    vapply(unclass(x), format_time, character(1)),
    quote = FALSE,
    right = TRUE
  )
  invisible(x)
}
