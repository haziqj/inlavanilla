#' Leave-one-out cross-validation and WAIC for an INLAvanilla fit
#'
#' @description
#' `loo()` estimates how well the model predicts new data by leave-one-out
#' cross-validation (LOO), and `waic()` by the widely applicable information
#' criterion. Both work from the fit alone, with no refitting and no
#' posterior draws, by a second-order Taylor expansion of each unit's
#' log-likelihood (Alhyari et al., 2026), as in INLAvaan.
#'
#' Both need the log-likelihood split into one term per observation, or per
#' independent unit, which [inla_vanilla()] never sees. You supply it as
#' `loglik`: write the log-likelihood terms and the log-prior as separate
#' functions, and `logpost` as their sum.
#'
#' @details
#' Write \eqn{\ell_u(\vartheta) = \log p(y_u \mid \vartheta)} for the
#' log-likelihood of unit \eqn{u}. The leave-one-out predictive density of
#' unit \eqn{u}, its conditional predictive ordinate, is
#' \deqn{\mathrm{CPO}_u = p(y_u \mid y_{-u}) =
#'   1 / E[1 / p(y_u \mid \vartheta) \mid y],}
#' an expectation over the full posterior. The posterior is taken as the
#' Gaussian \eqn{N(\tilde\vartheta, \Omega)} of the fit, with the VB-corrected
#' mean (the mode when `vb = FALSE`) and the Laplace covariance, and
#' \eqn{\ell_u} is expanded to second order about \eqn{\tilde\vartheta}, with
#' score \eqn{s_u} and Hessian \eqn{H_u}. The expectation is then a Gaussian
#' integral:
#' \deqn{\log \mathrm{CPO}_u = \ell_u
#'   - \tfrac12 s_u^\top (\Omega^{-1} + H_u)^{-1} s_u
#'   + \tfrac12 \log |I + \Omega H_u|.}
#' The same expansion gives the log pointwise predictive density
#' \deqn{\mathrm{lpd}_u = \log E[p(y_u \mid \vartheta) \mid y] = \ell_u
#'   + \tfrac12 s_u^\top (\Omega^{-1} - H_u)^{-1} s_u
#'   - \tfrac12 \log |I - \Omega H_u|}
#' and the WAIC penalty
#' \eqn{p_{\mathrm{waic},u} = \mathrm{Var}[\ell_u(\vartheta)] =
#' s_u^\top \Omega s_u + \tfrac12 \mathrm{tr}[(\Omega H_u)^2]}. Then
#' \eqn{\mathrm{elpd}_{\mathrm{loo}} = \sum_u \log \mathrm{CPO}_u},
#' \eqn{p_{\mathrm{loo}} = \sum_u (\mathrm{lpd}_u - \log \mathrm{CPO}_u)},
#' \eqn{\mathrm{elpd}_{\mathrm{waic}} = \sum_u (\mathrm{lpd}_u -
#' p_{\mathrm{waic},u})}, and `looic` and `waic` are \eqn{-2} times the elpd.
#' Each standard error is \eqn{\sqrt{n \, \mathrm{Var}}} of the pointwise
#' terms over the \eqn{n} units, as in the loo package.
#'
#' The scores and Hessians come from finite differences of `loglik` along the
#' principal axes of \eqn{\Omega}, a hundredth of a posterior SD apart:
#' \eqn{m^2 + m + 1} calls of `loglik` for \eqn{m} parameters, each returning
#' every unit.
#'
#' **Existence and the leverage k.** For each unit, the eigenvalues \eqn{k}
#' of \eqn{-\Omega H_u} are the shares of the posterior precision the unit
#' carries along its directions. The second-order \eqn{\log \mathrm{CPO}_u}
#' exists only when every \eqn{k < 1}: at \eqn{k \ge 1}, the rest of the data
#' and the prior cannot identify some combination of the parameters without
#' the unit, and \eqn{E[1 / p(y_u \mid \vartheta)]} diverges. If any unit
#' fails, every LOO estimate falls back to the first-order
#' \eqn{\log \mathrm{CPO}_u = \ell_u - \tfrac12 s_u^\top \Omega s_u}, with a
#' warning naming the units, rather than mixing orders or dropping units. The
#' largest \eqn{k} of each unit is reported as `k`, a measure of its
#' influence. Likewise the second-order \eqn{\mathrm{lpd}_u} needs every
#' \eqn{k > -1}. Where that fails, WAIC falls back to first order, at which it
#' equals the first-order LOO, and `p_loo` uses the first-order term for that
#' unit.
#'
#' The expansion rests on the Gaussian approximation of the posterior, as the
#' VB correction does. A cross-check that uses the skew-normal marginals is
#' Pareto-smoothed importance sampling on copula draws:
#' `loo::loo(sampling(fit, 4000, fun = loglik, ...))`.
#'
#' @section Comparing models:
#' The results have the layout of the loo package's objects, so
#' `loo::loo_compare()` compares fits to the same units. By hand, with `d`
#' the difference of two fits' `pointwise[, "elpd_loo"]`, the difference in
#' elpd is `sum(d)` with standard error `sqrt(length(d)) * sd(d)`.
#'
#' @param object An object of class `inla_vanilla`.
#' @param loglik A function `loglik(theta, ...)` returning the log-likelihood
#'   terms, one per observation or independent unit, at the named parameter
#'   vector `theta` (on the scale of `logpost`). Their sum plus the log-prior
#'   is `logpost`. Names, if any, label the units.
#' @param ... Further arguments passed to `loglik`, such as the data.
#'
#' @returns `loo()` returns an object of class `inla_vanilla_loo`, which also
#'   has class `loo`: a list with
#'   * `estimates`: a matrix with rows `elpd_loo`, `p_loo` and `looic` and
#'     columns `Estimate` and `SE`;
#'   * `pointwise`: a matrix with one row per unit and columns `elpd_loo` (the
#'     log CPO), `p_loo`, `looic` and `k` (the unit's largest leverage);
#'   * `order`: the order of the expansion used, 2 or 1;
#'   * `first_order_units`: the units whose second-order term does not exist,
#'     if any;
#'   * `theta`: the point of the expansion, and `calls`: the number of calls
#'     of `loglik`;
#'   * `warnings`: the messages of any warnings raised.
#'
#'   `waic()` returns an object of class `inla_vanilla_waic` (also `waic` and
#'   `loo`) with the same layout, with rows and columns `elpd_waic`, `p_waic`
#'   and `waic`.
#'
#' @references
#' Alhyari, M., Jamil, H., Montcho, H. and Rue, H. (2026). Deterministic
#' leave-one-cluster-out cross-validation for multilevel Bayesian structural
#' equation models. \doi{10.48550/arXiv.2609.00670}
#'
#' Vehtari, A., Gelman, A. and Gabry, J. (2017). Practical Bayesian model
#' evaluation using leave-one-out cross-validation and WAIC. *Statistics and
#' Computing*, 27, 1413--1432. \doi{10.1007/s11222-016-9696-4}
#'
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and
#' widely applicable information criterion in singular learning theory.
#' *Journal of Machine Learning Research*, 11, 3571--3594.
#'
#' @examples
#' # Poisson counts: is a log-linear trend better than a constant rate?
#' y <- c(2, 3, 6, 7, 8, 9, 10, 12, 15)
#' x <- seq_along(y) - 5
#'
#' loglik <- function(theta, x, y) {
#'   dpois(y, exp(theta[["a"]] + theta[["b"]] * x), log = TRUE)
#' }
#' logprior <- function(theta) sum(dnorm(theta, 0, 5, log = TRUE))
#' logpost <- function(theta, x, y) sum(loglik(theta, x, y)) + logprior(theta)
#'
#' fit <- inla_vanilla(logpost, c(a = 0, b = 0), x = x, y = y)
#' loo(fit, loglik, x = x, y = y)
#' waic(fit, loglik, x = x, y = y)
#' @export
loo <- function(object, ...) {
  UseMethod("loo")
}

#' @rdname loo
#' @export
waic <- function(object, ...) {
  UseMethod("waic")
}

#' @rdname loo
#' @export
loo.inla_vanilla <- function(object, loglik, ...) {
  if (missing(loglik)) {
    abort(
      "{.arg loglik} is missing: a function returning one log-likelihood
       term per observation."
    )
  }
  check_function(loglik, "loglik")
  # Bind `...` here, so a data argument cannot partially match a formal of
  # an internal helper
  loglik_x <- function(theta) loglik(theta, ...)
  tp <- taylor_pointwise(object, loglik_x)
  loo_from_taylor(tp, silent = isTRUE(object$silent))
}

#' @rdname loo
#' @export
waic.inla_vanilla <- function(object, loglik, ...) {
  if (missing(loglik)) {
    abort(
      "{.arg loglik} is missing: a function returning one log-likelihood
       term per observation."
    )
  }
  check_function(loglik, "loglik")
  loglik_x <- function(theta) loglik(theta, ...)
  tp <- taylor_pointwise(object, loglik_x)
  waic_from_taylor(tp, silent = isTRUE(object$silent))
}

# The loo package's generics loo() and waic() mask these when it is attached
# after inlavanilla. Their first argument is `x`, so the same methods are
# registered for them through wrappers with that signature. (Through these
# generics a data argument named `x` would bind to the fit, which is why
# inlavanilla's own generics call it `object`.)
#' @rawNamespace S3method(loo::loo, inla_vanilla, loo_x_inla_vanilla)
loo_x_inla_vanilla <- function(x, loglik, ...) {
  loo.inla_vanilla(x, loglik, ...)
}

#' @rawNamespace S3method(loo::waic, inla_vanilla, waic_x_inla_vanilla)
waic_x_inla_vanilla <- function(x, loglik, ...) {
  waic.inla_vanilla(x, loglik, ...)
}

# The per-unit Taylor terms at the fit's Gaussian N(theta, Omega). Finite
# differences run in whitened coordinates u, with theta + L u and
# L L' = Omega, so they give the whitened score L' s_u and Hessian
# L' H_u L directly; every formula below is invariant to the choice of L.
# Steps of h = 0.01 posterior SDs: in these units the terms of an expansion
# shrink with the amount of data, so truncation error is negligible, and
# rounding error is about 1e-12 |l_u|.
taylor_pointwise <- function(fit, loglik, h = 0.01) {
  theta <- fit$mean_vb
  parnames <- names(theta)
  m <- length(theta)
  L <- t(chol(fit$cov))
  n_calls <- m * m + m + 1L
  calls <- 0L
  last_x <- theta
  n <- NA_integer_

  pb <- NULL
  if (!isTRUE(fit$silent)) {
    pb <- cli::cli_progress_bar(
      "Differentiating loglik",
      total = n_calls,
      clear = TRUE
    )
  }
  at <- function(u) {
    x <- theta + drop(L %*% u)
    names(x) <- parnames
    last_x <<- x
    calls <<- calls + 1L
    if (!is.null(pb)) {
      cli::cli_progress_update(id = pb)
    }
    v <- tryCatch(loglik(x), error = function(e) {
      where <- paste(parnames, "=", signif(last_x, 4), collapse = ", ")
      abort(
        c(
          "{.arg loglik} failed.",
          "i" = "The parameter value was {where}."
        ),
        parent = e,
        call = NULL
      )
    })
    if (is.matrix(v) && ncol(v) == 1L) {
      v <- v[, 1L]
    }
    if (!is.numeric(v) || is.matrix(v)) {
      abort(
        "{.arg loglik} must return a numeric vector, not
         {.obj_type_friendly {v}}."
      )
    }
    if (is.na(n)) {
      n <<- length(v)
    } else if (length(v) != n) {
      abort(
        "{.arg loglik} returned {length(v)} value{?s} here but {n} at the
         posterior mean; it must return one term per unit every time."
      )
    }
    if (!all(is.finite(v))) {
      bad <- which(!is.finite(v))[1L]
      abort(
        c(
          "{.arg loglik} is not finite for unit {bad} near the posterior
           mean.",
          "i" = "The parameter value was
                 {paste(parnames, '=', signif(x, 4), collapse = ', ')}."
        )
      )
    }
    v
  }

  l0 <- at(numeric(m))
  units <- names(l0)
  if (n < 2L) {
    abort(
      "{.arg loglik} returned a single value: it must return one
       log-likelihood term per observation (or unit), not their sum."
    )
  }
  fp <- fm <- matrix(0, n, m)
  for (i in seq_len(m)) {
    e <- numeric(m)
    e[i] <- h
    fp[, i] <- at(e)
    fm[, i] <- at(-e)
  }
  S <- (fp - fm) / (2 * h)
  H <- array(0, c(m, m, n))
  for (i in seq_len(m)) {
    H[i, i, ] <- (fp[, i] - 2 * l0 + fm[, i]) / h^2
  }
  if (m > 1L) {
    for (i in seq_len(m - 1L)) {
      for (j in (i + 1L):m) {
        e <- numeric(m)
        e[c(i, j)] <- h
        hij <- (at(e) +
          at(-e) -
          fp[, i] -
          fm[, i] -
          fp[, j] -
          fm[, j] +
          2 * l0) /
          (2 * h^2)
        H[i, j, ] <- hij
        H[j, i, ] <- hij
      }
    }
  }
  if (!is.null(pb)) {
    cli::cli_progress_done(id = pb)
  }

  # s' Omega s in whitened coordinates
  quad1 <- rowSums(S^2)
  log_cpo_2 <- lpd_2 <- k_max <- k_min <- ssq <- rep(NA_real_, n)
  for (u in seq_len(n)) {
    eg <- eigen(H[,, u], symmetric = TRUE)
    lam <- eg$values # eigenvalues of L' H_u L, so k = -lam
    proj2 <- drop(crossprod(eg$vectors, S[u, ]))^2
    k_max[u] <- max(-lam, 0)
    k_min[u] <- min(-lam)
    ssq[u] <- sum(lam^2)
    if (all(1 + lam > 0)) {
      log_cpo_2[u] <- l0[u] -
        0.5 * sum(proj2 / (1 + lam)) +
        0.5 * sum(log1p(lam))
    }
    if (all(1 - lam > 0)) {
      lpd_2[u] <- l0[u] + 0.5 * sum(proj2 / (1 - lam)) - 0.5 * sum(log1p(-lam))
    }
  }

  list(
    n = n,
    units = units,
    l0 = unname(l0),
    quad1 = quad1,
    ssq = ssq,
    log_cpo_1 = unname(l0) - 0.5 * quad1,
    lpd_1 = unname(l0) + 0.5 * quad1,
    log_cpo_2 = log_cpo_2,
    lpd_2 = lpd_2,
    k_max = k_max,
    k_min = k_min,
    theta = theta,
    calls = calls
  )
}

# Name units by label when loglik gave names, by position otherwise.
unit_labels <- function(tp, which) {
  if (is.null(tp$units)) as.character(which) else tp$units[which]
}

loo_se <- function(v) {
  sqrt(length(v) * stats::var(v))
}

# Assemble the LOO from the Taylor terms. A missing second-order log CPO
# means the unit's case-deletion integral really diverges, so a first-order
# stand-in for that unit alone would be wrong by an unbounded amount, and
# dropping it would score the model on fewer units: every estimate falls to
# first order instead. A missing second-order lpd is an artefact of the
# quadratic, so p_loo uses that unit's first-order term (as INLAvaan does).
loo_from_taylor <- function(tp, silent = FALSE) {
  cpo_ok <- !is.na(tp$log_cpo_2)
  lpd_ok <- !is.na(tp$lpd_2)
  second <- all(cpo_ok)
  warnings <- character(0)
  bad <- unit_labels(tp, which(!cpo_ok))
  if (!second) {
    n_bad <- length(bad)
    msg <- c(
      "Using the first-order LOO for all units.",
      "i" = "{n_bad} unit{?s} {?has/have} a leverage k of 1 or more:
             {.val {bad}}.",
      "i" = "{cli::qty(n_bad)}Without {?it/them}, the rest of the data and the
             prior cannot identify some combination of the parameters, and
             the second-order term does not exist."
    )
    warnings <- c(warnings, cli::ansi_strip(cli::format_inline(msg[[1]])))
    if (!silent) {
      warn(msg)
    }
  }
  log_cpo <- if (second) tp$log_cpo_2 else tp$log_cpo_1
  # lpd_1 - log_cpo_1 = s' Omega s
  p_u <- if (second) {
    ifelse(lpd_ok, tp$lpd_2 - tp$log_cpo_2, tp$quad1)
  } else {
    tp$quad1
  }
  est <- rbind(
    elpd_loo = c(sum(log_cpo), loo_se(log_cpo)),
    p_loo = c(sum(p_u), loo_se(p_u)),
    looic = c(-2 * sum(log_cpo), 2 * loo_se(log_cpo))
  )
  colnames(est) <- c("Estimate", "SE")
  pointwise <- cbind(
    elpd_loo = log_cpo,
    p_loo = p_u,
    looic = -2 * log_cpo,
    k = tp$k_max
  )
  rownames(pointwise) <- tp$units
  structure(
    list(
      estimates = est,
      pointwise = pointwise,
      order = if (second) 2L else 1L,
      first_order_units = bad,
      theta = tp$theta,
      calls = tp$calls,
      warnings = warnings
    ),
    class = c("inla_vanilla_loo", "loo")
  )
}

# Assemble the WAIC from the same terms. It reads no case-deletion term, so
# only the lpd must exist; where it does not, every estimate falls to first
# order, where elpd_waic equals the first-order LOO exactly.
waic_from_taylor <- function(tp, silent = FALSE) {
  lpd_ok <- !is.na(tp$lpd_2)
  second <- all(lpd_ok)
  warnings <- character(0)
  bad <- unit_labels(tp, which(!lpd_ok))
  if (!second) {
    msg <- c(
      "Using the first-order WAIC for all units.",
      "i" = "{length(bad)} unit{?s} {?has/have} a leverage k of -1 or less:
             {.val {bad}}, and the second-order lpd does not exist."
    )
    warnings <- c(warnings, cli::ansi_strip(cli::format_inline(msg[[1]])))
    if (!silent) {
      warn(msg)
    }
  }
  lpd <- if (second) tp$lpd_2 else tp$lpd_1
  p_w <- if (second) tp$quad1 + 0.5 * tp$ssq else tp$quad1
  elpd_u <- lpd - p_w
  est <- rbind(
    elpd_waic = c(sum(elpd_u), loo_se(elpd_u)),
    p_waic = c(sum(p_w), loo_se(p_w)),
    waic = c(-2 * sum(elpd_u), 2 * loo_se(elpd_u))
  )
  colnames(est) <- c("Estimate", "SE")
  pointwise <- cbind(elpd_waic = elpd_u, p_waic = p_w, waic = -2 * elpd_u)
  rownames(pointwise) <- tp$units
  structure(
    list(
      estimates = est,
      pointwise = pointwise,
      order = if (second) 2L else 1L,
      first_order_units = bad,
      theta = tp$theta,
      calls = tp$calls,
      warnings = warnings
    ),
    class = c("inla_vanilla_waic", "waic", "loo")
  )
}

#' @export
print.inla_vanilla_loo <- function(x, digits = 1, ...) {
  print_taylor_header(x, "Leave-one-out cross-validation", digits)
  pw <- x$pointwise
  worst <- which.max(pw[, "k"])
  unit <- if (is.null(rownames(pw))) worst else rownames(pw)[worst]
  k_worst <- pw[worst, "k"]
  if (x$order == 2L) {
    cat_status(
      "ok",
      "Second-order expansion for every unit. Largest leverage k = {fmt_num(k_worst)}, for unit {unit}."
    )
  } else {
    bad <- x$first_order_units
    cat_status(
      "warn",
      "First-order expansion: {length(bad)} unit{?s} {?has/have} a leverage k of 1 or more ({.val {bad}})."
    )
  }
  invisible(x)
}

#' @export
print.inla_vanilla_waic <- function(x, digits = 1, ...) {
  print_taylor_header(x, "WAIC", digits)
  if (x$order == 2L) {
    cat_status("ok", "Second-order expansion for every unit.")
  } else {
    bad <- x$first_order_units
    cat_status(
      "warn",
      "First-order expansion: {length(bad)} unit{?s} {?has/have} a leverage k of -1 or less ({.val {bad}})."
    )
  }
  invisible(x)
}

print_taylor_header <- function(x, title, digits) {
  cli::cat_rule(
    left = cli::style_bold(title),
    right = cli::col_grey(nrow(x$pointwise), " units"),
    line_col = "cyan"
  )
  cli::cat_line()
  cat_table(x$estimates, digits = digits + 1, scale = rep(1, 3))
  cli::cat_line()
}
