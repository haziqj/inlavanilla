################################################################################
#
# Benchmark of the simplex degree-5 rule for the VB correction (October 2026).
# vb_nodes() is now a degree-5 rule built on a regular simplex, with
# (m + 1)(m + 2) + 1 nodes, at most. With a gradient vb_correct() takes the
# expected gradient over its nodes; without one, Stein's identity on the same
# nodes gives it from log-posterior values alone. Before, both routes used
# the fully symmetric degree-5 rule on the principal axes, with 2m^2 + 1
# nodes: on the user's gradient, or on bias-corrected forward differences of
# the log-posterior (m + 1 calls each). That rule and those differences are
# rebuilt here and swapped into vb_correct() for the old routes.
#
# For each model and route the script reports the error of the shift against
# a high-accuracy reference, in posterior SDs, the calls of the VB stage, and
# its time against the mode-and-Hessian stage (laplace_fit()) of the same
# setting, warm. The reference solves the VB fixed point with the analytic
# gradient (central differences for the Weibull model) under tensor
# Gauss-Hermite quadrature for m <= 3, and 2^17 antithetic pairs of Monte
# Carlo draws otherwise. Takes about ten minutes, most of it in the old route
# without a gradient, and the reference, at m = 80. Run from the package
# root:
#   Rscript lab/vb-simplex-benchmark.R
#
################################################################################

## ----- Configuration ----------------------------------------------------------

suppressMessages(devtools::load_all(".", quiet = TRUE))
ns <- asNamespace("inlavanilla")
laplace_fit <- get("laplace_fit", ns)
vb_correct <- get("vb_correct", ns)
simplex_nodes <- get("vb_nodes", ns)
fd_grad <- get("fd_grad", ns)
ctrl <- get("vanilla_control", ns)(list())

## ----- Models ---------------------------------------------------------------------

models <- list()

set.seed(1)
yA <- rnorm(10, 1, 2)
models$normal_n10 <- list(
  init = c(mu = 0, ls = 0),
  lp = function(t) {
    s <- exp(t[[2]])
    sum(dnorm(yA, t[[1]], s, log = TRUE)) +
      dnorm(t[[1]], 0, 10, log = TRUE) +
      dexp(s, 0.5, log = TRUE) +
      t[[2]]
  },
  gr = function(t) {
    s <- exp(t[[2]])
    r <- yA - t[[1]]
    c(sum(r) / s^2 - t[[1]] / 100, -10 + sum(r^2) / s^2 - 0.5 * s + 1)
  }
)
models$betabinom <- list(
  init = c(eta = 0),
  lp = function(t) {
    p <- plogis(t[[1]])
    dbinom(2, 12, p, log = TRUE) + log(p * (1 - p))
  },
  gr = function(t) 3 - 14 * plogis(t[[1]])
)
yC <- c(0, 1, 0, 2, 1)
models$poisgamma <- list(
  init = c(ll = 0),
  lp = function(t) {
    l <- exp(t[[1]])
    sum(dpois(yC, l, log = TRUE)) + dgamma(l, 1, 1, log = TRUE) + t[[1]]
  },
  gr = function(t) sum(yC) - 6 * exp(t[[1]]) + 1
)
set.seed(5)
yF <- rweibull(15, shape = 1.5, scale = 2)
models$weibull <- list(
  init = c(lk = 0, ll = 0),
  lp = function(t) {
    sum(dweibull(yF, exp(t[[1]]), exp(t[[2]]), log = TRUE)) +
      sum(dnorm(t, 0, 2, log = TRUE))
  },
  gr = NULL
)

# Logistic regressions: the main vignette's mtcars model, those of
# lab/vb-benchmark.R (near separation from m = 20), and two larger ones with
# the linear predictor's SD held near 1.2
logit_spec <- function(X, yy) {
  force(X)
  force(yy)
  list(
    init = stats::setNames(rep(0, ncol(X)), paste0("b", seq_len(ncol(X)))),
    lp = function(t) {
      eta <- as.numeric(X %*% t)
      sum(yy * eta - log1p(exp(eta))) + sum(dnorm(t, 0, 2.5, log = TRUE))
    },
    gr = function(t) {
      as.numeric(crossprod(X, yy - plogis(as.numeric(X %*% t)))) - t / 6.25
    }
  )
}
models$mtcars <- logit_spec(
  cbind(1, scale(mtcars$wt), scale(mtcars$hp)),
  mtcars$am
)
for (spec in list(c(5, 40, 3), c(10, 80, 4), c(20, 150, 6), c(30, 250, 7))) {
  set.seed(spec[3])
  X <- cbind(1, matrix(rnorm(spec[2] * (spec[1] - 1)), spec[2]))
  yy <- rbinom(spec[2], 1, plogis(X %*% c(-0.5, rnorm(spec[1] - 1, 0, 0.8))))
  models[[paste0("logit_m", spec[1])]] <- logit_spec(X, yy)
}
for (m in c(40, 80)) {
  set.seed(100 + m)
  X <- cbind(1, matrix(rnorm(10 * m * (m - 1)), 10 * m))
  beta <- c(-0.5, rnorm(m - 1, 0, 1.2 / sqrt(m - 1)))
  models[[paste0("logit_m", m, "_mild")]] <- logit_spec(
    X,
    rbinom(10 * m, 1, plogis(X %*% beta))
  )
}

## ----- The old rule and routes ------------------------------------------------------

# The fully symmetric degree-5 rule on the principal axes, 2m^2 + 1 nodes
symmetric_nodes <- function(Sigma) {
  m <- ncol(Sigma)
  e <- eigen(Sigma, symmetric = TRUE)
  At <- t(sqrt(3) * sweep(e$vectors, 2, sqrt(pmax(e$values, 0)), "*"))
  nodes <- rbind(0, At, -At)
  weights <- c(1 + (m^2 - 7 * m) / 18, rep((4 - m) / 18, 2 * m))
  if (m > 1) {
    pairs <- utils::combn(m, 2)
    Ai <- At[pairs[1, ], , drop = FALSE]
    Aj <- At[pairs[2, ], , drop = FALSE]
    nodes <- rbind(nodes, Ai + Aj, Ai - Aj, -Ai + Aj, -Ai - Aj)
    weights <- c(weights, rep(1 / 36, 4 * ncol(pairs)))
  }
  list(nodes = nodes, weights = weights)
}

# vb_correct() on the old rule, with the user's gradient or, without one,
# forward differences a thousandth of a posterior SD long, their bias removed
# with the Laplace curvature
old_route <- function(lp, gr, lap) {
  if (is.null(gr)) {
    h <- 1e-3 * sqrt(diag(lap$Sigma))
    bias <- 0.5 * h * diag(lap$H)
    gr <- function(x) {
      f0 <- lp(x)
      g <- vapply(
        seq_along(x),
        function(j) {
          xp <- x
          xp[j] <- x[j] + h[j]
          (lp(xp) - f0) / h[j]
        },
        0
      )
      g + bias
    }
  }
  assignInNamespace("vb_nodes", symmetric_nodes, ns)
  on.exit(assignInNamespace("vb_nodes", simplex_nodes, ns))
  vb_correct(lp, gr, lap$mode, lap$Sigma, lap$R, lap$L, ctrl)
}

new_route <- function(lp, gr, lap) {
  vb_correct(lp, gr, lap$mode, lap$Sigma, lap$R, lap$L, ctrl)
}

## ----- Helpers ----------------------------------------------------------------------

# Seconds per run of f(), warm: repeated until half a second has passed. A
# run already timed at over a second (`first`) is not repeated.
time_it <- function(f, first = 0) {
  if (first > 1) {
    return(first)
  }
  f()
  n <- 0
  t0 <- proc.time()[[3]]
  repeat {
    f()
    n <- n + 1
    dt <- proc.time()[[3]] - t0
    if (dt > 0.5) break
  }
  dt / n
}

# The log-posterior and gradient with the parameters named and the calls
# counted
counting <- function(model) {
  n_lp <- 0L
  n_gr <- 0L
  nm <- names(model$init)
  list(
    lp = function(x) {
      names(x) <- nm
      n_lp <<- n_lp + 1L
      v <- model$lp(x)[[1]]
      if (is.na(v)) -Inf else v
    },
    gr = if (!is.null(model$gr)) {
      function(x) {
        names(x) <- nm
        n_gr <<- n_gr + 1L
        as.numeric(model$gr(x))
      }
    },
    calls = function() c(n_lp, n_gr),
    reset = function() {
      n_lp <<- 0L
      n_gr <<- 0L
    }
  )
}

# The VB shift under a high-accuracy rule, solved to 1e-10 SDs
reference_shift <- function(model, lap) {
  m <- length(lap$mode)
  sds <- sqrt(diag(lap$Sigma))
  L <- t(chol(lap$Sigma))
  gr <- if (is.null(model$gr)) {
    function(x) fd_grad(model$lp, x, 1e-5 * sds)
  } else {
    model$gr
  }
  if (m <= 3) {
    k <- c(80, 40, 20)[m]
    i <- seq_len(k - 1)
    J <- matrix(0, k, k)
    J[cbind(i, i + 1)] <- J[cbind(i + 1, i)] <- sqrt(i / 2)
    e <- eigen(J, symmetric = TRUE)
    z <- as.matrix(expand.grid(rep(list(e$values * sqrt(2)), m)))
    w <- apply(as.matrix(expand.grid(rep(list(e$vectors[1, ]^2), m))), 1, prod)
  } else {
    set.seed(42)
    z <- matrix(rnorm(2^17 * m), ncol = m)
    z <- rbind(z, -z)
    w <- rep(1 / nrow(z), nrow(z))
  }
  X <- z %*% t(L)
  shift <- numeric(m)
  for (it in 1:500) {
    mu <- lap$mode + shift
    g <- vapply(seq_len(nrow(X)), function(b) gr(mu + X[b, ]), numeric(m))
    step <- as.numeric(lap$Sigma %*% (matrix(g, nrow = m) %*% w))
    shift <- shift + step
    if (max(abs(step) / sds) < 1e-10) {
      return(shift)
    }
  }
  stop("the reference did not converge")
}

## ----- Benchmark ----------------------------------------------------------------------

rows <- list()
for (name in names(models)) {
  model <- models[[name]]
  w <- counting(model)
  # One Laplace fit per model, with the gradient where there is one, so that
  # both settings share the reference
  lap <- laplace_fit(w$lp, w$gr, model$init, ctrl)
  sds <- sqrt(diag(lap$Sigma))
  ref <- reference_shift(model, lap)
  settings <- if (is.null(model$gr)) "nograd" else c("grad", "nograd")
  for (setting in settings) {
    gr <- if (setting == "grad") w$gr else NULL
    t_lap <- time_it(function() laplace_fit(w$lp, gr, model$init, ctrl))
    for (route in c("old", "new")) {
      run <- if (route == "old") old_route else new_route
      w$reset()
      t0 <- proc.time()[[3]]
      vb <- run(w$lp, gr, lap)
      first <- proc.time()[[3]] - t0
      calls <- w$calls()
      t_vb <- time_it(function() run(w$lp, gr, lap), first)
      rows[[length(rows) + 1]] <- data.frame(
        model = name,
        m = length(lap$mode),
        setting = setting,
        route = route,
        shift = max(abs(ref) / sds),
        error = max(abs(vb$shift - ref) / sds),
        sweeps = vb$iterations,
        fallback = vb$fallback,
        calls = sum(calls),
        ratio = t_vb / t_lap
      )
    }
  }
  message(name, " done")
}
out <- do.call(rbind, rows)
saveRDS(out, file.path(tempdir(), "vb-simplex-benchmark.rds"))
for (setting in c("grad", "nograd")) {
  cat("\n", if (setting == "grad") "With" else "Without", "a gradient\n")
  wide <- reshape(
    out[out$setting == setting, names(out) != "setting"],
    idvar = c("model", "m", "shift"),
    timevar = "route",
    direction = "wide"
  )
  print(format(wide, digits = 2), row.names = FALSE)
}
