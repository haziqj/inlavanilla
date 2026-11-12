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
#' mean (the mode when `vb = "off"`) and the Laplace covariance, and
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
#' **The diagonal shortcut.** With `hessian = "diagonal"` only the
#' \eqn{2m + 1} calls along the axes are made. They give each score and the
#' diagonal of each Hessian in those axes, and so the trace
#' \eqn{\mathrm{tr}(\Omega H_u)}. To first order in the leverages,
#' \eqn{\log |I + \Omega H_u| = \mathrm{tr}(\Omega H_u)}, and then
#' \deqn{\log \mathrm{CPO}_u = \ell_u - \tfrac12 s_u^\top \Omega s_u
#'   + \tfrac12 \mathrm{tr}(\Omega H_u),}
#' with \eqn{p_u = s_u^\top \Omega s_u} for both criteria, which then agree
#' unit by unit. Summed over the units, the trace is minus the effective number
#' of parameters of the likelihood,
#' \eqn{p_D = \mathrm{tr}(-\Omega \sum_u H_u)}, so the total is the
#' first-order one less \eqn{p_D / 2}, which is what the first order alone
#' leaves out (Alhyari et al., 2026). The terms dropped are of second order in
#' the leverages, about \eqn{m^2 / n} in all for \eqn{n} units of similar
#' weight. Use it when every unit is a small part of the data and
#' \eqn{m^2} calls are too many. The `k` reported is then the sum of a unit's
#' leverages, \eqn{-\mathrm{tr}(\Omega H_u)}, and a unit at 1 or more gives
#' a warning.
#'
#' The calls use the `cores` of the fit, or the plan set with `future::plan()`
#' when the fit has none, as the scans of [inla_vanilla()] do (see its section
#' on parallel computation), with the same results as a serial run. Calls
#' expected to take under a second in all stay serial.
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
#' @section Silent fits:
#' On a fit made with `silent = TRUE`, `loo()` and `waic()` print nothing and
#' raise no warnings, their own or those of `loglik`; the warnings are kept
#' in the result's `warnings`, and a fallback to first order shows in `order`
#' and `first_order_units`.
#'
#' @section Other objects:
#' For any other object, `loo()` and `waic()` call the loo package's
#' functions of the same name, if it is installed, so attaching inlavanilla
#' after loo (or a package that uses it, such as rstan or brms) does not get
#' in their way.
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
#' @param ... Further arguments passed to `loglik`, such as the data. As in
#'   [stats::optim()], a name that abbreviates `object` or `loglik` (such as
#'   `l`) would be taken as that argument, so give such data a longer name.
#'   If the loo package is attached after inlavanilla, its generics
#'   `loo(x, ...)` and `waic(x, ...)` are the ones called, and they take an
#'   argument named `x` as the fit: then do not name data `x`, or call
#'   `inlavanilla::loo()`.
#' @param hessian `"full"` (the default) takes every second derivative of each
#'   unit's log-likelihood, with \eqn{m^2 + m + 1} calls of `loglik`.
#'   `"diagonal"` takes the scores and the trace only, with \eqn{2m + 1}
#'   calls: see the Details. It comes after the dots, so it must be named in
#'   full.
#'
#' @returns `loo()` returns an object of class `inla_vanilla_loo`, which also
#'   has class `loo`: a list with
#'   * `estimates`: a matrix with rows `elpd_loo`, `p_loo` and `looic` and
#'     columns `Estimate` and `SE`;
#'   * `pointwise`: a matrix with one row per unit and columns `elpd_loo` (the
#'     log CPO), `p_loo`, `looic` and `k` (the unit's largest leverage, or the
#'     sum of its leverages with `hessian = "diagonal"`);
#'   * `order`: the order of the expansion used, 2 or 1 (1 with
#'     `hessian = "diagonal"`, where the trace is added to the first order);
#'   * `hessian`: `"full"` or `"diagonal"`;
#'   * `first_order_units`: the units whose second-order term does not exist,
#'     if any, and with `hessian = "diagonal"`, `trace_units`: the units whose
#'     leverages add up to 1 or more;
#'   * `theta`: the point of the expansion, and `calls`: the number of calls
#'     of `loglik`;
#'   * `warnings`: the messages of any warnings raised.
#'
#'   `waic()` returns an object of class `inla_vanilla_waic` (also `waic` and
#'   `loo`) with the same layout, with rows and columns `elpd_waic`, `p_waic`
#'   and `waic`.
#'
#' @references
#' Alhyari, M., Jamil, H., Montcho, H., & Rue, H. (2026). *Deterministic
#' leave-one-cluster-out cross-validation for multilevel Bayesian models*
#' (arXiv:2609.00670 \[stat.ME\]). arXiv.
#' <https://doi.org/10.48550/arXiv.2609.00670>
#'
#' Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model
#' evaluation using leave-one-out cross-validation and WAIC. *Statistics and
#' Computing*, *27*(5), 1413--1432. <https://doi.org/10.1007/s11222-016-9696-4>
#'
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and
#' widely applicable information criterion in singular learning theory.
#' *Journal of Machine Learning Research*, *11*, 3571--3594.
#' <https://jmlr.org/papers/v11/watanabe10a.html>
#'
#' @examples
#' # Poisson counts: is a log-linear trend better than a constant rate?
#' y <- c(2, 3, 6, 7, 8, 9, 10, 12, 15)
#' time <- seq_along(y) - 5
#'
#' loglik <- function(theta, time, y) {
#'   dpois(y, exp(theta[["a"]] + theta[["b"]] * time), log = TRUE)
#' }
#' logprior <- function(theta) sum(dnorm(theta, 0, 5, log = TRUE))
#' logpost <- function(theta, time, y) {
#'   sum(loglik(theta, time, y)) + logprior(theta)
#' }
#'
#' fit <- inla_vanilla(logpost, c(a = 0, b = 0), time = time, y = y)
#' loo(fit, loglik, time = time, y = y)
#' waic(fit, loglik, time = time, y = y)
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
loo.default <- function(object, ...) {
  if (!requireNamespace("loo", quietly = TRUE)) {
    abort(
      c(
        "There is no {.fn loo} method for {.obj_type_friendly {object}}.",
        "i" = "Install the loo package for its {.fn loo} methods."
      )
    )
  }
  loo::loo(object, ...)
}

#' @rdname loo
#' @export
waic.default <- function(object, ...) {
  if (!requireNamespace("loo", quietly = TRUE)) {
    abort(
      c(
        "There is no {.fn waic} method for {.obj_type_friendly {object}}.",
        "i" = "Install the loo package for its {.fn waic} methods."
      )
    )
  }
  loo::waic(object, ...)
}

#' @rdname loo
#' @export
loo.inla_vanilla <- function(
  object,
  loglik,
  ...,
  hessian = c("full", "diagonal")
) {
  if (missing(loglik)) {
    abort(
      "{.arg loglik} is missing: a function returning one log-likelihood
       term per observation."
    )
  }
  check_function(loglik, "loglik")
  hessian <- check_choice(hessian, c("full", "diagonal"), "hessian")
  # Bind `...` here, so a data argument cannot partially match a formal of
  # an internal helper
  loglik_x <- function(theta) loglik(theta, ...)
  tp <- taylor_pointwise(object, loglik_x, hessian = hessian)
  loo_from_taylor(tp, silent = isTRUE(object$silent))
}

#' @rdname loo
#' @export
waic.inla_vanilla <- function(
  object,
  loglik,
  ...,
  hessian = c("full", "diagonal")
) {
  if (missing(loglik)) {
    abort(
      "{.arg loglik} is missing: a function returning one log-likelihood
       term per observation."
    )
  }
  check_function(loglik, "loglik")
  hessian <- check_choice(hessian, c("full", "diagonal"), "hessian")
  loglik_x <- function(theta) loglik(theta, ...)
  tp <- taylor_pointwise(object, loglik_x, hessian = hessian)
  waic_from_taylor(tp, silent = isTRUE(object$silent))
}

# The loo package's generics loo() and waic() mask these when it is attached
# after inlavanilla. Their first argument is `x`, so the same methods are
# registered for them through wrappers with that signature. (Through these
# generics a data argument named `x` would bind to the fit, which is why
# inlavanilla's own generics call it `object`.)
#' @rawNamespace S3method(loo::loo, inla_vanilla, loo_x_inla_vanilla)
loo_x_inla_vanilla <- function(
  x,
  loglik,
  ...,
  hessian = c("full", "diagonal")
) {
  loo.inla_vanilla(x, loglik, ..., hessian = hessian)
}

#' @rawNamespace S3method(loo::waic, inla_vanilla, waic_x_inla_vanilla)
waic_x_inla_vanilla <- function(
  x,
  loglik,
  ...,
  hessian = c("full", "diagonal")
) {
  waic.inla_vanilla(x, loglik, ..., hessian = hessian)
}

# The per-unit Taylor terms at the fit's Gaussian N(theta, Omega). Finite
# differences run in whitened coordinates u, with theta + L u and
# L L' = Omega, so they give the whitened score L' s_u and Hessian
# L' H_u L directly; every formula below is invariant to the choice of L.
# Steps of h = 0.01 posterior SDs: in these units the terms of an expansion
# shrink with the amount of data, so truncation error is negligible, and
# rounding error is about 1e-12 |l_u|. With hessian = "diagonal" only the
# centre and the points along the axes are evaluated, 2m + 1 in all, which
# give the scores and the diagonal of each whitened Hessian, and with it the
# trace of Omega H_u.
taylor_pointwise <- function(fit, loglik, h = 0.01, hessian = "full") {
  theta <- fit$mean_vb
  parnames <- names(theta)
  m <- length(theta)
  L <- t(chol(fit$cov))
  full <- hessian == "full"
  n_calls <- if (full) m * m + m + 1L else 2L * m + 1L
  calls <- 0L
  last_x <- theta
  n <- NA_integer_
  # Warnings from loglik are recorded once each and, for a silent fit,
  # muffled
  silent <- isTRUE(fit$silent)
  seen <- character()
  on_warning <- function(w) {
    msg <- warning_text(w)
    if (!msg %in% seen) {
      seen <<- c(seen, msg)
    }
    if (silent) {
      tryInvokeRestart("muffleWarning")
    }
  }

  prog <- progress_new(!silent)
  at <- function(u) {
    x <- theta + drop(L %*% u)
    names(x) <- parnames
    last_x <<- x
    calls <<- calls + 1L
    v <- withCallingHandlers(
      tryCatch(loglik(x), error = function(e) {
        where <- paste(parnames, "=", signif(last_x, 4), collapse = ", ")
        abort(
          c(
            "{.arg loglik} failed.",
            "i" = "The parameter value was {where}."
          ),
          parent = e,
          call = NULL
        )
      }),
      warning = on_warning
    )
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

  # The first calls are needed anyway, and say how long the rest will take
  U <- taylor_points(m, h, pairs = full)
  secs <- numeric(3L)
  t0 <- elapsed()
  l0 <- at(U[, 1L])
  secs[[1L]] <- elapsed() - t0
  units <- names(l0)
  if (n < 2L) {
    abort(
      "{.arg loglik} returned a single value: it must return one
       log-likelihood term per observation (or unit), not their sum."
    )
  }
  vals <- matrix(NA_real_, n, n_calls)
  vals[, 1L] <- l0
  for (k in 2:3) {
    t0 <- elapsed()
    vals[, k] <- at(U[, k])
    secs[[k]] <- elapsed() - t0
  }

  # On workers, if the fit has `cores` or the user a future plan, and the calls
  # are slow enough. Before the progress line starts, because it can print a
  # message
  data <- list(
    loglik = loglik,
    theta = theta,
    parnames = parnames,
    L = L,
    U = U,
    n = n
  )
  setup <- par_setup(
    list(loglik),
    taylor_probe,
    data,
    est_s = taylor_est_s(secs, n_calls),
    n = n_calls - 3L,
    note = function(text) progress_note(prog, text),
    label = "Differentiating loglik",
    what = "loglik",
    cores = fit$cores,
    how = "fit with the `cores` argument or call future::plan() first"
  )
  on.exit(setup$restore(), add = TRUE)
  progress_start(
    prog,
    "Differentiating loglik, {cli::pb_current}/{cli::pb_total} calls",
    total = n_calls
  )
  # A point ticks once, however many times its job reports it
  ticked <- logical(n_calls)
  tick_once <- function(k) {
    if (!ticked[[k]]) {
      ticked[[k]] <<- TRUE
      progress_tick(prog)
    }
  }
  for (k in 1:3) {
    tick_once(k)
  }
  done <- seq_len(n_calls) <= 3L
  if (setup$jobs > 1L) {
    items <- which(!done)
    groups <- lapply(
      par_groups(length(items), setup$jobs),
      function(i) items[i]
    )
    out <- withCallingHandlers(
      par_run(groups, taylor_job, data, setup, tick_once),
      warning = on_warning
    )
    for (g in seq_along(groups)) {
      if (!is.null(out[[g]])) {
        vals[, groups[[g]]] <- out[[g]]
        done[groups[[g]]] <- TRUE
        calls <- calls + length(groups[[g]])
      }
    }
  }
  # Every point a worker did not return, and with no workers every one, is
  # evaluated here, so an error reads as it does in a serial run
  for (k in which(!done)) {
    vals[, k] <- at(U[, k])
    tick_once(k)
  }

  fp <- vals[, 2L * seq_len(m), drop = FALSE]
  fm <- vals[, 2L * seq_len(m) + 1L, drop = FALSE]
  S <- (fp - fm) / (2 * h)
  # The diagonal of each unit's whitened Hessian, a row per unit. Its sum is
  # the trace of Omega H_u.
  D <- (fp - 2 * l0 + fm) / h^2
  trace <- rowSums(D)
  quad1 <- rowSums(S^2)
  if (!full) {
    progress_close(prog)
    return(list(
      n = n,
      units = units,
      l0 = unname(l0),
      quad1 = quad1,
      trace = trace,
      log_cpo_1 = unname(l0) - 0.5 * quad1,
      lpd_1 = unname(l0) + 0.5 * quad1,
      hessian = hessian,
      theta = theta,
      calls = calls,
      workers = setup$jobs,
      loglik_warnings = seen
    ))
  }
  H <- array(0, c(m, m, n))
  for (i in seq_len(m)) {
    H[i, i, ] <- D[, i]
  }
  if (m > 1L) {
    k <- 2L * m + 1L
    for (i in seq_len(m - 1L)) {
      for (j in (i + 1L):m) {
        hij <- (vals[, k + 1L] +
          vals[, k + 2L] -
          fp[, i] -
          fm[, i] -
          fp[, j] -
          fm[, j] +
          2 * l0) /
          (2 * h^2)
        H[i, j, ] <- hij
        H[j, i, ] <- hij
        k <- k + 2L
      }
    }
  }
  progress_close(prog)

  # quad1 is s' Omega s in whitened coordinates
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
    trace = trace,
    log_cpo_1 = unname(l0) - 0.5 * quad1,
    lpd_1 = unname(l0) + 0.5 * quad1,
    log_cpo_2 = log_cpo_2,
    lpd_2 = lpd_2,
    k_max = k_max,
    k_min = k_min,
    hessian = hessian,
    theta = theta,
    calls = calls,
    workers = setup$jobs,
    loglik_warnings = seen
  )
}

# The points of taylor_pointwise() in whitened coordinates, a column each: the
# centre, then +h and -h along each axis, then +h and -h along each pair of
# axes, m^2 + m + 1 in all. Without the pairs, the first 2m + 1 of them.
taylor_points <- function(m, h, pairs = TRUE) {
  U <- matrix(0, m, if (pairs) m * m + m + 1L else 2L * m + 1L)
  k <- 1L
  add <- function(e) {
    U[, k + 1L] <<- e
    U[, k + 2L] <<- -e
    k <<- k + 2L
  }
  for (i in seq_len(m)) {
    e <- numeric(m)
    e[i] <- h
    add(e)
  }
  if (pairs && m > 1L) {
    for (i in seq_len(m - 1L)) {
      for (j in (i + 1L):m) {
        e <- numeric(m)
        e[c(i, j)] <- h
        add(e)
      }
    }
  }
  U
}

# How long the calls of loglik would take in all, from the times of the first
# few. The shortest of them, because the byte compiler slows the first calls of
# a session.
taylor_est_s <- function(secs, n_calls) {
  min(secs) * n_calls
}

# loglik at point k of taylor_points()
taylor_value <- function(data, k) {
  x <- data$theta + drop(data$L %*% data$U[, k])
  names(x) <- data$parnames
  v <- data$loglik(x)
  if (is.matrix(v) && ncol(v) == 1L) v[, 1L] else v
}

# loglik at the centre, to compare a worker with the parent
taylor_probe <- function(ids, data) {
  unname(taylor_value(data, 1L))
}

# The points `ids` on a worker. A value that taylor_pointwise() would reject
# fails the job, so that the parent evaluates the point and words the error.
taylor_job <- function(ids, data) {
  vals <- matrix(NA_real_, data$n, length(ids))
  for (i in seq_along(ids)) {
    v <- taylor_value(data, ids[[i]])
    if (
      !is.numeric(v) ||
        is.matrix(v) ||
        length(v) != data$n ||
        !all(is.finite(v))
    ) {
      stop("loglik gave an unusable value")
    }
    vals[, i] <- v
    par_tick(ids[[i]])
  }
  vals
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
  if (identical(tp$hessian, "diagonal")) {
    return(trace_from_taylor(tp, "loo", silent))
  }
  cpo_ok <- !is.na(tp$log_cpo_2)
  lpd_ok <- !is.na(tp$lpd_2)
  second <- all(cpo_ok)
  warnings <- tp$loglik_warnings
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
      hessian = "full",
      first_order_units = bad,
      theta = tp$theta,
      calls = tp$calls,
      warnings = warnings
    ),
    class = c("inla_vanilla_loo", "loo")
  )
}

# The LOO and the WAIC from the scores and the trace alone. To first order
# in the leverages, log |I + Omega H_u| is tr(Omega H_u) and the quadratic
# form is s' Omega s, so
#   log CPO_u = l_u - s' Omega s / 2 + tr(Omega H_u) / 2,
#   lpd_u     = l_u + s' Omega s / 2 + tr(Omega H_u) / 2,
# and p_u = s' Omega s for both criteria, which then agree unit by unit.
# Summed over the units the trace is minus the effective number of parameters
# of the likelihood, tr(-Omega sum_u H_u), so this is the first-order total
# less half of it. What is left out is of second order in the leverages, so
# the sum of a unit's leverages, -tr(Omega H_u), is reported as its k, and a
# unit at 1 or more is warned about.
trace_from_taylor <- function(tp, what, silent = FALSE) {
  k_sum <- -tp$trace
  elpd_u <- tp$log_cpo_1 + 0.5 * tp$trace
  p_u <- tp$quad1
  warnings <- tp$loglik_warnings
  bad <- unit_labels(tp, which(k_sum >= 1))
  if (length(bad) > 0L) {
    n_bad <- length(bad)
    msg <- c(
      "The trace may not stand in for the second-order term.",
      "i" = "{n_bad} unit{?s} {?has/have} leverages that add up to 1 or
             more: {.val {bad}}.",
      "i" = "Use {.code hessian = \"full\"} to check {cli::qty(n_bad)}
             {?it/them}."
    )
    warnings <- c(warnings, cli::ansi_strip(cli::format_inline(msg[[1]])))
    if (!silent) {
      warn(msg)
    }
  }
  names_est <- if (what == "loo") {
    c("elpd_loo", "p_loo", "looic")
  } else {
    c("elpd_waic", "p_waic", "waic")
  }
  est <- rbind(
    c(sum(elpd_u), loo_se(elpd_u)),
    c(sum(p_u), loo_se(p_u)),
    c(-2 * sum(elpd_u), 2 * loo_se(elpd_u))
  )
  dimnames(est) <- list(names_est, c("Estimate", "SE"))
  pointwise <- cbind(elpd_u, p_u, -2 * elpd_u)
  colnames(pointwise) <- names_est
  if (what == "loo") {
    pointwise <- cbind(pointwise, k = k_sum)
  }
  rownames(pointwise) <- tp$units
  structure(
    list(
      estimates = est,
      pointwise = pointwise,
      order = 1L,
      hessian = "diagonal",
      first_order_units = character(),
      trace_units = bad,
      theta = tp$theta,
      calls = tp$calls,
      warnings = warnings
    ),
    class = if (what == "loo") {
      c("inla_vanilla_loo", "loo")
    } else {
      c("inla_vanilla_waic", "waic", "loo")
    }
  )
}

# Assemble the WAIC from the same terms. It reads no case-deletion term, so
# only the lpd must exist; where it does not, every estimate falls to first
# order, where elpd_waic equals the first-order LOO exactly.
waic_from_taylor <- function(tp, silent = FALSE) {
  if (identical(tp$hessian, "diagonal")) {
    return(trace_from_taylor(tp, "waic", silent))
  }
  lpd_ok <- !is.na(tp$lpd_2)
  second <- all(lpd_ok)
  warnings <- tp$loglik_warnings
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
      hessian = "full",
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
  if (identical(x$hessian, "diagonal")) {
    cat_status(
      if (length(x$trace_units) > 0L) "warn" else "ok",
      "First-order expansion with the trace of the second-order term. Largest summed leverage k = {fmt_num(k_worst)}, for unit {unit}."
    )
  } else if (x$order == 2L) {
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
  if (identical(x$hessian, "diagonal")) {
    cat_status(
      if (length(x$trace_units) > 0L) "warn" else "ok",
      "First-order expansion with the trace of the second-order term."
    )
  } else if (x$order == 2L) {
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
    width = out_width(),
    left = cli::style_bold(title),
    right = cli::col_grey(nrow(x$pointwise), " units"),
    line_col = vanilla_hex
  )
  cli::cat_line()
  cat_table(x$estimates, digits = digits + 1, scale = rep(1, 3))
  cli::cat_line()
}
