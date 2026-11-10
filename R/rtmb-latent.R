#' Draw latent variables from an RTMB object
#'
#' Builds the `fun` of [sampling()] for a model whose latent variables are the
#' random effects of an RTMB or TMB object. For each draw of the parameters
#' from the fit, the function it returns draws the random effects from their
#' Gaussian approximation given those parameters, so that together the draws
#' follow the posterior of the latent variables.
#'
#' An object built with `random` integrates its random effects out by a
#' Laplace approximation. Each call of `obj$fn()` finds their mode given the
#' parameters, and the Hessian there is the precision of a Gaussian
#' approximation of their distribution given the parameters. A draw is the
#' mode plus Gaussian noise with that precision. The approximation is exact
#' when the random effects are Gaussian given the parameters, as in a linear
#' mixed model, and close when each has plenty of data. It is weakest for a
#' random effect with few discrete observations: see the article "Latent
#' variables" on the package website.
#'
#' The object must return a negative log-density, as RTMB and TMB expect, and
#' its parameters other than the random effects must be those of the fit, in
#' the same order and on the same scale. Both hold when the fit is made from
#' the object itself, with `logpost = function(theta) -obj$fn(theta)` and
#' `init = obj$par`. A fit made from another function of the same model, such
#' as a closed form with the random effects integrated out by hand, can use an
#' object built only for this purpose.
#'
#' The sparse Cholesky factor comes from the Matrix package, which RTMB and
#' TMB install. An object's tape does not survive [saveRDS()], so build the
#' object again in a new session.
#'
#' @param obj An object from `RTMB::MakeADFun()` or `TMB::MakeADFun()`, built
#'   with `random` naming the latent variables.
#'
#' @returns A function of the parameters, to pass as `fun` to [sampling()]. It
#'   returns one draw of every random effect of `obj`, in the object's order,
#'   named after the parameter it belongs to: `u[1]`, `u[2]` and so on for a
#'   vector `u`.
#'
#' @examplesIf requireNamespace("RTMB", quietly = TRUE)
#' # Random intercepts u, integrated out by RTMB's Laplace approximation
#' set.seed(1)
#' group <- rep(1:8, each = 5)
#' y <- 2 + rnorm(8)[group] + rnorm(40, sd = 0.5)
#' nll <- function(p) {
#'   -sum(RTMB::dnorm(y, p$mu + p$u[group], exp(p$log_sigma), log = TRUE)) -
#'     sum(RTMB::dnorm(p$u, 0, exp(p$log_tau), log = TRUE)) -
#'     RTMB::dnorm(p$mu, 0, 10, log = TRUE) -
#'     RTMB::dnorm(p$log_tau, 0, 2, log = TRUE) -
#'     RTMB::dnorm(p$log_sigma, 0, 2, log = TRUE)
#' }
#' obj <- RTMB::MakeADFun(
#'   nll,
#'   list(mu = 0, log_tau = 0, log_sigma = 0, u = numeric(8)),
#'   random = "u",
#'   silent = TRUE
#' )
#' fit <- inla_vanilla(
#'   function(theta) -obj$fn(theta),
#'   obj$par,
#'   grad = function(theta) -as.vector(obj$gr(theta)),
#'   silent = TRUE
#' )
#'
#' # One row per draw of the parameters, one column per random intercept
#' u <- sampling(fit, nsamp = 500, fun = rtmb_latent(obj))
#' round(colMeans(u), 2)
#' @export
rtmb_latent <- function(obj) {
  check_tmb_object(obj)
  env <- obj$env
  random <- env$random
  n_par <- length(obj$par)
  latent <- latent_names(names(env$last.par)[random], length(random))
  function(theta) {
    if (length(theta) != n_par) {
      abort(
        "{.arg obj} has {n_par} parameter{?s}, but the draw has
         {length(theta)}: its parameters must be those of the fit, in the
         same order.",
        call = NULL
      )
    }
    value <- obj$fn(unname(theta))
    par <- env$last.par
    step <- if (all(is.finite(value))) {
      tryCatch(
        precision_root(
          env$spHess(par, random = TRUE),
          stats::rnorm(length(random))
        ),
        error = function(e) NULL
      )
    }
    if (is.null(step) || !all(is.finite(step))) {
      abort(
        "The random effects of {.arg obj} have no Gaussian approximation at a
         draw of the parameters: its inner optimisation failed, or the
         Hessian there is not positive definite.",
        call = NULL
      )
    }
    stats::setNames(par[random] + step, latent)
  }
}

check_tmb_object <- function(obj) {
  ok <- is.list(obj) &&
    is.function(obj$fn) &&
    is.environment(obj$env) &&
    is.function(obj$env$spHess)
  if (!ok) {
    abort(
      "{.arg obj} must be an object from {.fn RTMB::MakeADFun} or
       {.fn TMB::MakeADFun}, not {.obj_type_friendly {obj}}."
    )
  }
  if (length(obj$env$random) == 0L) {
    abort(
      "{.arg obj} has no random effects. Build it with {.arg random}, naming
       the latent variables."
    )
  }
  invisible(obj)
}

# x = precision_root(H, e) has covariance H^-1 when e is standard normal. A
# sparse H = P' L L' P gives x = P' L^-T e, and a base matrix H = R' R gives
# x = R^-1 e.
precision_root <- function(H, e) {
  if (!is.matrix(H) && requireNamespace("Matrix", quietly = TRUE)) {
    # RTMB hands back the same matrix at every call, with new values, and
    # Cholesky() keeps its factor in the matrix it is given. So a matrix that
    # has been factorised once would give that first factor again: 1 * H is a
    # new matrix with no factor kept.
    ch <- Matrix::Cholesky(1 * H, perm = TRUE, LDL = FALSE)
    x <- Matrix::solve(ch, e, system = "Lt")
    return(as.vector(Matrix::solve(ch, x, system = "Pt")))
  }
  backsolve(chol(as.matrix(H)), e)
}

# TMB repeats a parameter's name for each of its elements: u, u, u becomes
# u[1], u[2], u[3], and a name that appears once is kept.
latent_names <- function(x, n) {
  if (length(x) != n) {
    x <- rep("latent", n)
  }
  count <- as.vector(table(x)[x])
  index <- stats::ave(seq_along(x), x, FUN = seq_along)
  ifelse(count > 1L, sprintf("%s[%d]", x, index), x)
}
