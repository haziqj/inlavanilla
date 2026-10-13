#' Draw from an approximate posterior
#'
#' Draws from the approximate joint posterior by a Gaussian copula: the
#' skew-normal marginals of the fit are coupled through the correlation matrix
#' \eqn{R} of the Laplace covariance. Each draw takes
#' \eqn{z \sim N(0, R)}, maps it to uniforms \eqn{u_j = \Phi(z_j)} and then to
#' each parameter's skew-normal quantile, \eqn{\vartheta_j = Q_j(u_j)}.
#'
#' For a fit with `marginals = "gaussian"` this is a draw from the joint
#' Gaussian \eqn{N(\tilde\vartheta, \Omega)}, centred at the VB mean (or at
#' the mode when `vb = FALSE`).
#'
#' Draws never call the log-posterior, so they are cheap: large `nsamp` costs
#' little. Use `fun` for natural-scale or derived quantities, such as
#' `exp(log_sigma)` or a covariance matrix built from several parameters.
#'
#' @param object An object to draw from.
#' @param nsamp Number of draws.
#' @param fun Optional function applied to each draw. It receives one draw as a
#'   named numeric vector, like `logpost`, and returns a numeric vector of a
#'   fixed length (named, ideally). When `NULL`, the parameters themselves are
#'   returned.
#' @param ... For the `inla_vanilla` method, further arguments passed to
#'   `fun`. As in [stats::optim()], a name that abbreviates `nsamp` or `fun`
#'   would be taken as that argument: a data argument named `n` would set the
#'   number of draws. Give such data a longer name, or wrap it in `fun`.
#'
#' @returns A numeric matrix with one row per draw: the parameters on the
#'   unconstrained scale, or the values of `fun`.
#'
#' @examples
#' set.seed(1)
#' y <- rnorm(30, mean = 5, sd = 2)
#' logpost <- function(theta, y) {
#'   sigma <- exp(theta[["log_sigma"]])
#'   sum(dnorm(y, theta[["mu"]], sigma, log = TRUE)) +
#'     dnorm(theta[["mu"]], 0, 10, log = TRUE) +
#'     dexp(sigma, log = TRUE) + theta[["log_sigma"]]
#' }
#' fit <- inla_vanilla(logpost, c(mu = 0, log_sigma = 0), y = y, silent = TRUE)
#'
#' draws <- sampling(fit, nsamp = 2000)
#' head(draws)
#'
#' # A derived quantity: the coefficient of variation sigma / mu
#' cv <- sampling(fit, nsamp = 2000, fun = function(theta) {
#'   c(cv = exp(theta[["log_sigma"]]) / theta[["mu"]])
#' })
#' quantile(cv, c(0.025, 0.5, 0.975))
#' @export
sampling <- function(object, ...) {
  UseMethod("sampling")
}

#' @rdname sampling
#' @export
sampling.inla_vanilla <- function(object, nsamp = 1000L, fun = NULL, ...) {
  if (
    !is.numeric(nsamp) || length(nsamp) != 1L || !is.finite(nsamp) || nsamp < 1
  ) {
    abort("{.arg nsamp} must be a single positive number.")
  }
  nsamp <- as.integer(nsamp)
  if (!is.null(fun)) {
    check_function(fun, "fun")
  }

  sn <- object$sn
  m <- nrow(sn)
  U <- chol(stats::cov2cor(object$cov))
  z <- matrix(stats::rnorm(nsamp * m), nsamp, m) %*% U
  theta <- matrix(0, nsamp, m, dimnames = list(NULL, rownames(sn)))
  for (j in seq_len(m)) {
    if (sn[j, "alpha"] == 0) {
      # a Gaussian marginal: the quantile map is linear
      theta[, j] <- sn[j, "xi"] + sn[j, "omega"] * z[, j]
    } else {
      # z -> Q_j(Phi(z)), the skew-normal quantile of each normal score
      map <- sn_score_map(sn[j, "xi"], sn[j, "omega"], sn[j, "alpha"])
      theta[, j] <- map(z[, j])
    }
  }
  if (is.null(fun)) {
    return(theta)
  }
  # Bind `...` here rather than forwarding it, so an argument meant for `fun`
  # cannot partially match a formal of apply_draws()
  fun_x <- function(t) fun(t, ...)
  apply_draws(theta, fun_x)
}

# Apply `fun` to each row of `theta` (passed as a named vector) and bind the
# results into a matrix.
apply_draws <- function(theta, fun) {
  first <- fun(theta[1, ])
  if (!is.numeric(first) && !is.logical(first)) {
    abort(
      "{.arg fun} must return a numeric vector, not {.obj_type_friendly {first}}."
    )
  }
  k <- length(first)
  if (k == 0L) {
    abort("{.arg fun} returned an empty vector.")
  }
  n <- nrow(theta)
  rest <- vapply(
    seq_len(n)[-1],
    function(i) {
      v <- fun(theta[i, ])
      if (length(v) != k) {
        abort(
          "{.arg fun} returned {length(v)} value{?s} for draw {i}, but {k}
           for the first draw; it must return the same length every time.",
          call = NULL
        )
      }
      as.numeric(v)
    },
    numeric(k)
  )
  # vapply() gives k x (n - 1), or a vector when k = 1; rows are draws
  out <- matrix(c(as.numeric(first), rest), nrow = n, ncol = k, byrow = TRUE)
  colnames(out) <- names(first)
  out
}
