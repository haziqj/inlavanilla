################################################################################
#
# Shared set-up of the one-step VB lab idea: the models, the ground truth
# (quadrature or HMC), the candidate routes to the VB location, and the
# generators of the simulation. Sourced by the other scripts, which run from
# the package root. Reads the package's internals and changes nothing in it.
#
################################################################################

## ----- Configuration ---------------------------------------------------------

smoke <- "--smoke" %in% commandArgs(trailingOnly = TRUE)
idea_dir <- "lab/ideas/2026-11-one-step-vb"
cache_dir <- file.path(idea_dir, "cache") # gitignored
dir.create(cache_dir, showWarnings = FALSE)
suppressMessages(pkgload::load_all(".", quiet = TRUE))
ctrl <- vanilla_control(list())
n_chains <- 4
n_warm <- if (smoke) 200 else 1500
n_keep <- if (smoke) 500 else 12500
l_max <- 30

## ----- Models ----------------------------------------------------------------

models <- list()

set.seed(1)
yA <- rnorm(10, 1, 2)
norm_spec <- function(y) {
  force(y)
  n <- length(y)
  list(
    init = c(mu = 0, ls = 0),
    lp = function(t) {
      s <- exp(t[[2]])
      sum(dnorm(y, t[[1]], s, log = TRUE)) +
        dnorm(t[[1]], 0, 10, log = TRUE) +
        dexp(s, 0.5, log = TRUE) +
        t[[2]]
    },
    gr = function(t) {
      s <- exp(t[[2]])
      r <- y - t[[1]]
      c(sum(r) / s^2 - t[[1]] / 100, -n + sum(r^2) / s^2 - 0.5 * s + 1)
    }
  )
}
models$normal_n10 <- norm_spec(yA)
models$normal_n5 <- norm_spec(yA[1:5])

logit1_spec <- function(a, b) {
  force(a)
  force(b)
  list(
    init = c(eta = qlogis(a / (a + b))),
    lp = function(t) a * t[[1]] - (a + b) * log1p(exp(t[[1]])),
    gr = function(t) a - (a + b) * plogis(t[[1]])
  )
}
models$beta_3_11 <- logit1_spec(3, 11)
models$beta_2_30 <- logit1_spec(2, 30)
loggamma_spec <- function(a) {
  force(a)
  list(
    init = c(ll = log(a) + 0.1),
    lp = function(t) a * t[[1]] - exp(t[[1]]),
    gr = function(t) a - exp(t[[1]])
  )
}
models$loggamma_1 <- loggamma_spec(1)
models$loggamma_2 <- loggamma_spec(2)
models$loggamma_5 <- loggamma_spec(5)

sn_spec <- function(alpha) {
  force(alpha)
  list(
    init = c(x = 0.3),
    lp = function(t) {
      dnorm(t[[1]], log = TRUE) + pnorm(alpha * t[[1]], log.p = TRUE)
    },
    gr = function(t) {
      x <- t[[1]]
      -x +
        alpha *
          exp(dnorm(alpha * x, log = TRUE) - pnorm(alpha * x, log.p = TRUE))
    }
  )
}
models$skewnorm_a2 <- sn_spec(2)
models$skewnorm_a5 <- sn_spec(5)
models$skewnorm_a10 <- sn_spec(10)

models$boundary <- list(
  init = c(a = 1.5),
  lp = function(t) {
    x <- t[[1]] + 1
    if (x <= 0) -Inf else 3 * log(x) - x
  },
  gr = function(t) {
    x <- t[[1]] + 1
    if (x <= 0) NaN else 3 / x - 1
  }
)

banana_spec <- function(b) {
  force(b)
  list(
    init = c(x1 = 0.1, x2 = 0.1),
    lp = function(t) -0.5 * t[[1]]^2 - 0.5 * (t[[2]] - b * t[[1]]^2)^2,
    gr = function(t) {
      r <- t[[2]] - b * t[[1]]^2
      c(-t[[1]] + 2 * b * t[[1]] * r, -r)
    }
  )
}
models$banana_b03 <- banana_spec(0.3)
models$banana_b1 <- banana_spec(1)

set.seed(5)
yF <- rweibull(15, shape = 1.5, scale = 2)
models$weibull <- list(
  init = c(lk = 0, ll = 0),
  lp = function(t) {
    sum(dweibull(yF, exp(t[[1]]), exp(t[[2]]), log = TRUE)) +
      sum(dnorm(t, 0, 2, log = TRUE))
  },
  gr = function(t) {
    k <- exp(t[[1]])
    r <- yF / exp(t[[2]])
    c(
      sum(1 + k * log(r) - r^k * k * log(r)) - t[[1]] / 4,
      sum(-k + k * r^k) - t[[2]] / 4
    )
  }
)

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
set.seed(140)
X <- cbind(1, matrix(rnorm(400 * 39), 400))
beta <- c(-0.5, rnorm(39, 0, 1.2 / sqrt(39)))
models$logit_m40 <- logit_spec(X, rbinom(400, 1, plogis(X %*% beta)))

pois_spec <- function(X, yy) {
  force(X)
  force(yy)
  list(
    init = stats::setNames(rep(0, ncol(X)), paste0("b", seq_len(ncol(X)))),
    lp = function(t) {
      eta <- as.numeric(X %*% t)
      sum(yy * eta - exp(eta)) + sum(dnorm(t, 0, 2.5, log = TRUE))
    },
    gr = function(t) {
      as.numeric(crossprod(X, yy - exp(as.numeric(X %*% t)))) - t / 6.25
    }
  )
}
set.seed(11)
X <- cbind(1, matrix(rnorm(25 * 4), 25))
models$pois_m5 <- pois_spec(X, rpois(25, exp(-0.7 + 0.5 * X[, 2])))
set.seed(12)
X <- cbind(1, matrix(rnorm(60 * 14), 60))
models$pois_m15 <- pois_spec(X, rpois(60, exp(-0.7 + 0.5 * X[, 2])))

y8 <- c(28, 8, -3, 7, -1, 1, 18, 12)
s8 <- c(15, 10, 16, 11, 9, 11, 10, 18)
models$eight_schools <- list(
  init = c(mu = 5, lt = 1, stats::setNames(rep(0, 8), paste0("e", 1:8))),
  lp = function(t) {
    tau <- exp(t[[2]])
    e <- t[3:10]
    sum(dnorm(y8, t[[1]] + tau * e, s8, log = TRUE)) +
      dnorm(t[[1]], 0, 5, log = TRUE) +
      dcauchy(tau, 0, 5, log = TRUE) +
      t[[2]] +
      sum(dnorm(e, log = TRUE))
  },
  gr = function(t) {
    tau <- exp(t[[2]])
    e <- t[3:10]
    r <- (y8 - t[[1]] - tau * e) / s8^2
    c(
      sum(r) - t[[1]] / 25,
      sum(r * tau * e) + 1 - 2 * tau^2 / (25 + tau^2),
      r * tau - e
    )
  }
)

## ----- Political Democracy, with its cached HMC ------------------------------

hmc_file <- "lab/tk-width-poldem_cache/hmc.rds" # from lab/tk-width-poldem.R
if (file.exists(hmc_file)) {
  data(PoliticalDemocracy, package = "lavaan")
  ov <- c("x1", "x2", "x3", paste0("y", 1:8))
  Y <- as.matrix(PoliticalDemocracy[, ov])
  qmd <- readLines("vignettes/articles/rcpp.qmd")
  starts <- grep("^```\\{cat", qmd)
  cpp <- unlist(lapply(starts, function(s) {
    e <- s + which(qmd[(s + 1):length(qmd)] == "```")[1]
    qmd[(s + 1):(e - 1)]
  }))
  Rcpp::sourceCpp(code = paste(cpp, collapse = "\n"))
  hmc_pd <- readRDS(hmc_file)
  pd_names <- colnames(hmc_pd$chains[[1]])
  init_pd <- c(
    colMeans(Y),
    rep(1, 5),
    rep(0, 3),
    rep(0, 3),
    log(apply(Y, 2, var) / 2),
    rep(0, 6)
  )
  names(init_pd) <- pd_names
  models$poldem <- list(
    init = init_pd,
    lp = function(t) logpost_cpp(t, Y),
    gr = function(t) as.numeric(grad_cpp(t, Y)),
    draws = hmc_pd$chains
  )
}

## ----- Ground truth ----------------------------------------------------------

grid_truth <- function(lp, lap) {
  m <- length(lap$mode)
  L <- lap$L
  z1 <- if (m == 1) seq(-40, 40, by = 0.005) else seq(-30, 30, by = 0.1)
  z <- as.matrix(expand.grid(rep(list(z1), m)))
  theta <- sweep(z %*% t(L), 2, lap$mode, "+")
  v <- vapply(seq_len(nrow(theta)), function(i) lp(theta[i, ]), 0)
  v[is.na(v)] <- -Inf
  w <- exp(v - max(v[is.finite(v)]))
  w <- w / sum(w)
  list(mean = colSums(theta * w), se = 0)
}

hmc_truth <- function(lp, gr, lap) {
  m <- length(lap$mode)
  L <- lap$L
  centre <- lap$mode
  lp_u <- function(u) {
    v <- lp(centre + drop(L %*% u))
    if (is.finite(v)) v else -Inf
  }
  gr_u <- function(u) drop(crossprod(L, gr(centre + drop(L %*% u))))
  chain <- function(seed) {
    set.seed(seed)
    u <- rnorm(m) * 0.5
    lp0 <- lp_u(u)
    g <- gr_u(u)
    log_eps <- log(0.2)
    out <- matrix(NA_real_, n_keep, m)
    for (it in seq_len(n_warm + n_keep)) {
      e <- exp(log_eps) * runif(1, 0.8, 1.2)
      p <- rnorm(m)
      h0 <- lp0 - 0.5 * sum(p^2)
      u1 <- u
      g1 <- g
      p1 <- p + 0.5 * e * g1
      ok <- TRUE
      n_steps <- sample.int(l_max, 1)
      for (s in seq_len(n_steps)) {
        u1 <- u1 + e * p1
        g1 <- gr_u(u1)
        if (!all(is.finite(g1))) {
          ok <- FALSE
          break
        }
        if (s < n_steps) p1 <- p1 + e * g1
      }
      a <- 0
      if (ok) {
        p1 <- p1 + 0.5 * e * g1
        lp1 <- lp_u(u1)
        h1 <- lp1 - 0.5 * sum(p1^2)
        if (is.finite(h1)) a <- min(1, exp(h1 - h0))
      }
      if (runif(1) < a) {
        u <- u1
        lp0 <- lp1
        g <- g1
      }
      if (it <= n_warm) {
        log_eps <- log_eps + (a - 0.8) / sqrt(it + 10)
      } else {
        out[it - n_warm, ] <- centre + drop(L %*% u)
      }
    }
    out
  }
  draws_truth(lapply(seq_len(n_chains), chain), lap)
}

## ----- Routes ----------------------------------------------------------------

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
    gr = function(x) {
      names(x) <- nm
      n_gr <<- n_gr + 1L
      as.numeric(model$gr(x))
    },
    calls = function() c(n_lp, n_gr),
    reset = function() {
      n_lp <<- 0L
      n_gr <<- 0L
    }
  )
}

run_routes <- function(model, lap, w) {
  mode <- lap$mode
  Sigma <- lap$Sigma
  R <- lap$R
  L <- lap$L
  m <- length(mode)
  sds <- sqrt(diag(Sigma))
  gr <- w$gr
  lp <- w$lp
  rule <- vb_nodes(Sigma)
  Z <- rule$nodes
  wt <- rule$weights
  nn <- nrow(Z)
  maxstep <- 1.5

  cap <- function(step) {
    move <- max(abs(step) / sds)
    if (is.finite(move) && move > maxstep) step * maxstep / move else step
  }
  simplex_step <- function(shift) {
    mu <- mode + shift
    g0 <- gr(mu)
    s <- g0
    for (b in 2:nn) {
      s <- s + wt[b] * (gr(mu + Z[b, ]) - g0)
    }
    as.numeric(Sigma %*% s)
  }
  e <- eigen(Sigma, symmetric = TRUE)
  A <- sweep(e$vectors, 2, sqrt(pmax(e$values, 0)), "*")
  r_ax <- sqrt(3)
  axis_step <- function(shift) {
    mu <- mode + shift
    g0 <- gr(mu)
    s <- g0
    for (k in seq_len(m)) {
      s <- s +
        (gr(mu + r_ax * A[, k]) + gr(mu - r_ax * A[, k]) - 2 * g0) /
          (2 * r_ax^2)
    }
    as.numeric(Sigma %*% s)
  }
  stein_step <- function(shift) {
    mu <- mode + shift
    n_pairs <- (nn - 1) %/% 2
    plus <- seq_len(n_pairs)
    v <- vapply(2:nn, function(b) lp(mu + Z[b, ]), 0)
    as.numeric(crossprod(
      Z[1 + plus, , drop = FALSE],
      wt[1 + plus] * (v[plus] - v[n_pairs + plus])
    ))
  }
  anderson <- function(step1, step2, dx) {
    df <- step2 - step1
    denom <- sum(df^2)
    if (denom <= 0) {
      return(step2)
    }
    gam <- sum(step2 * df) / denom
    cand <- step2 - gam * (dx + df)
    move <- max(abs(step2) / sds)
    if (
      all(is.finite(cand)) && max(abs(cand) / sds) <= min(maxstep, 2 * move)
    ) {
      cand
    } else {
      step2
    }
  }
  quad_finish <- function(d1, d2) {
    u <- as.numeric(R %*% d1)
    v <- as.numeric(R %*% (d2 - d1))
    kappa <- sum(u * v) / sum(u^2)
    if (!is.finite(kappa) || kappa >= 0.25) {
      return(d2)
    }
    t <- if (abs(kappa) < 1e-12) 1 else (1 - sqrt(1 - 4 * kappa)) / (2 * kappa)
    d1 + t^2 * (d2 - d1)
  }

  out <- list()
  calls <- list()

  # The current solver, and the same with a tight tolerance
  w$reset()
  cur <- vb_correct(lp, gr, mode, Sigma, R, L, ctrl)
  out$current <- cur$shift
  calls$current <- w$calls()[2] + w$calls()[1]
  sweeps <- cur$iterations
  fallback <- cur$fallback

  # Simplex rule, two passes and finishes
  s1 <- cap(simplex_step(numeric(m)))
  out$one <- s1
  calls$one <- nn
  s2 <- cap(simplex_step(s1))
  out$two <- s1 + s2
  calls$two <- 2 * nn
  out$two_aa <- s1 + anderson(s1, s2, s1)
  calls$two_aa <- 2 * nn
  out$two_q <- quad_finish(s1, s1 + s2)
  calls$two_q <- 2 * nn

  # Axis rule, degree 3 and 2m + 1 gradients a pass
  a1 <- cap(axis_step(numeric(m)))
  out$axis1 <- a1
  calls$axis1 <- 2 * m + 1
  a2 <- cap(axis_step(a1))
  out$axis2 <- a1 + a2
  calls$axis2 <- 2 * (2 * m + 1)
  out$axis2_aa <- a1 + anderson(a1, a2, a1)
  calls$axis2_aa <- 2 * (2 * m + 1)
  ax <- a1 + a2
  n_ax <- 2
  for (it in 1:10) {
    st <- cap(axis_step(ax))
    ax <- ax + st
    n_ax <- n_ax + 1
    if (!is.finite(max(abs(st) / sds)) || max(abs(st) / sds) < 1e-3) break
  }
  out$axis_conv <- ax
  calls$axis_conv <- n_ax * (2 * m + 1)

  # The free trace term of gamma1_grad()
  Vscan <- sweep(Sigma, 2, sds, "/")
  trace_r <- function(j, d_outer = 0.01, h = 1e-5) {
    vj <- Vscan[, j]
    th <- mode + vj * d_outer
    tr <- 0
    for (k in seq_len(m)) {
      Lk <- L[, k]
      tr <- tr + sum(Lk * (gr(th - Lk * h) - gr(th + Lk * h))) / (2 * h)
    }
    -0.5 * (tr - m) / d_outer
  }
  tr <- sds * vapply(seq_len(m), trace_r, 0)
  tr <- cap(tr)
  out$trace <- tr
  calls$trace <- 0
  out$trace_one <- tr + cap(simplex_step(tr))
  calls$trace_one <- nn
  out$trace_axis <- tr + cap(axis_step(tr))
  calls$trace_axis <- 2 * m + 1

  # Skew-normal means minus the mode
  fit0 <- suppressWarnings(
    inla_vanilla(lp, model$init, gr, vb = FALSE, silent = TRUE)
  )
  mom <- sn_moments(fit0$sn[, "xi"], fit0$sn[, "omega"], fit0$sn[, "alpha"])
  sn <- as.numeric(mom$mean) - mode
  out$sn <- sn
  calls$sn <- 0
  out$sn_one <- sn + cap(simplex_step(sn))
  calls$sn_one <- nn
  out$sn_axis <- sn + cap(axis_step(sn))
  calls$sn_axis <- 2 * m + 1

  # No gradient: the current solver and one or two Stein passes
  w$reset()
  cur_ng <- vb_correct(lp, NULL, mode, Sigma, R, L, ctrl)
  out$ng_current <- cur_ng$shift
  calls$ng_current <- sum(w$calls())
  ng1 <- cap(stein_step(numeric(m)))
  out$ng_one <- ng1
  calls$ng_one <- nn - 1
  ng2 <- cap(stein_step(ng1))
  out$ng_two <- ng1 + ng2
  calls$ng_two <- 2 * (nn - 1)

  # The Laplace mode itself
  out$mode <- numeric(m)
  calls$mode <- 0

  list(
    shifts = out,
    calls = unlist(calls),
    sweeps = sweeps,
    fallback = fallback
  )
}

## ----- Truth from draws, and log Z -------------------------------------------

draws_truth <- function(chains, lap) {
  sds <- sqrt(diag(lap$Sigma))
  means <- vapply(chains, colMeans, numeric(length(sds)))
  all_draws <- do.call(rbind, chains)
  list(
    mean = rowMeans(means),
    se = max(apply(means, 1, sd) / sqrt(ncol(means)) / sds),
    cov = stats::cov(all_draws)
  )
}

# log Z by importance sampling, with a t proposal fitted to the draws
is_logz <- function(lp, mu, S, n = 30000, df = 8) {
  m <- length(mu)
  C <- t(chol(S * 1.1))
  z <- matrix(rnorm(n * m), n)
  u <- rchisq(n, df) / df
  x <- sweep(z %*% t(C) / sqrt(u), 2, mu, "+")
  q2 <- rowSums((z / sqrt(u))^2)
  log_q <- lgamma((df + m) / 2) -
    lgamma(df / 2) -
    0.5 * m * log(df * pi) -
    sum(log(diag(C))) -
    0.5 * (df + m) * log1p(q2 / df)
  v <- vapply(seq_len(n), function(i) lp(x[i, ]), 0) - log_q
  v[!is.finite(v)] <- -Inf
  top <- max(v)
  w <- exp(v - top)
  list(logz = top + log(mean(w)), ess = sum(w)^2 / sum(w^2))
}

## ----- Families --------------------------------------------------------------

gen_logit <- function(m, n) {
  function(seed) {
    set.seed(seed)
    X <- cbind(1, matrix(rnorm(n * (m - 1)), n))
    y <- rbinom(n, 1, plogis(X %*% c(-0.5, rnorm(m - 1, 0, 0.8))))
    logit_spec(X, y)
  }
}
gen_pois <- function(m, n) {
  function(seed) {
    set.seed(seed)
    X <- cbind(1, matrix(rnorm(n * (m - 1)), n))
    pois_spec(X, rpois(n, exp(-0.7 + 0.5 * X[, 2])))
  }
}
gen_glmm <- function(J, nj) {
  function(seed) {
    set.seed(seed)
    n <- J * nj
    g <- rep(seq_len(J), each = nj)
    x <- rnorm(n)
    u0 <- rnorm(J, 0, 0.8)
    y <- rbinom(n, 1, plogis(-0.5 + 0.7 * x + u0[g]))
    nm <- c("b0", "b1", "lt", paste0("u", seq_len(J)))
    list(
      init = stats::setNames(c(0, 0, -0.5, rep(0, J)), nm),
      lp = function(t) {
        tau <- exp(t[[3]])
        u <- t[-(1:3)]
        eta <- t[[1]] + t[[2]] * x + tau * u[g]
        sum(y * eta - log1p(exp(eta))) +
          sum(dnorm(t[1:2], 0, 2.5, log = TRUE)) -
          0.5 * tau^2 +
          t[[3]] -
          0.5 * sum(u^2)
      },
      gr = function(t) {
        tau <- exp(t[[3]])
        u <- t[-(1:3)]
        eta <- t[[1]] + t[[2]] * x + tau * u[g]
        r <- y - plogis(eta)
        c(
          sum(r) - t[[1]] / 6.25,
          sum(r * x) - t[[2]] / 6.25,
          sum(r * tau * u[g]) - tau^2 + 1,
          tau * as.numeric(tapply(r, factor(g, seq_len(J)), sum)) - u
        )
      }
    )
  }
}
gen_treg <- function(m, n, nu = 4) {
  function(seed) {
    set.seed(seed)
    X <- cbind(1, matrix(rnorm(n * (m - 2)), n))
    y <- as.numeric(X %*% c(1, rnorm(m - 2, 0, 0.5))) + rt(n, nu)
    nm <- c(paste0("b", seq_len(m - 1)), "ls")
    list(
      init = stats::setNames(c(median(y), rep(0, m - 2), 0), nm),
      lp = function(t) {
        b <- t[-m]
        z <- (y - as.numeric(X %*% b)) / exp(t[[m]])
        -n *
          t[[m]] -
          0.5 * (nu + 1) * sum(log1p(z^2 / nu)) +
          sum(dnorm(b, 0, 10, log = TRUE)) +
          dnorm(t[[m]], 0, 2, log = TRUE)
      },
      gr = function(t) {
        b <- t[-m]
        s <- exp(t[[m]])
        z <- (y - as.numeric(X %*% b)) / s
        k <- (nu + 1) * z / (nu + z^2)
        c(
          as.numeric(crossprod(X, k)) / s - b / 100,
          sum(k * z) - n - t[[m]] / 4
        )
      }
    )
  }
}
gen_cfa <- function(n, p = 6) {
  function(seed) {
    set.seed(seed)
    lam0 <- c(1, runif(p - 1, 0.6, 1.2))
    Y <- rnorm(n) %o% lam0 + matrix(rnorm(n * p, 0, 0.8), n)
    S <- crossprod(Y) / n
    nm <- c(paste0("l", 2:p), "lpsi", paste0("lth", 1:p))
    m <- 2 * p
    idx_l <- seq_len(p - 1)
    idx_psi <- p
    idx_th <- (p + 1):m
    build <- function(t) {
      lam <- c(1, t[idx_l])
      psi <- exp(t[[idx_psi]])
      th <- exp(t[idx_th])
      list(
        lam = lam,
        psi = psi,
        th = th,
        Sig = psi * tcrossprod(lam) + diag(th)
      )
    }
    list(
      init = stats::setNames(c(rep(1, p - 1), 0, rep(-0.5, p)), nm),
      lp = function(t) {
        b <- build(t)
        R <- tryCatch(chol(b$Sig), error = function(e) NULL)
        if (is.null(R)) {
          return(-Inf)
        }
        Ri <- backsolve(R, diag(p))
        -0.5 *
          n *
          (2 * sum(log(diag(R))) + sum(tcrossprod(Ri) * S)) +
          sum(dnorm(t[idx_l], 1, 1, log = TRUE)) +
          sum(dnorm(t[c(idx_psi, idx_th)], 0, 2, log = TRUE))
      },
      gr = function(t) {
        b <- build(t)
        Si <- solve(b$Sig)
        W <- Si %*% (b$Sig - S) %*% Si
        Wl <- as.numeric(W %*% b$lam)
        c(
          -n * b$psi * Wl[-1] - (t[idx_l] - 1),
          -0.5 * n * b$psi * sum(b$lam * Wl) - t[[idx_psi]] / 4,
          -0.5 * n * b$th * diag(W) - t[idx_th] / 4
        )
      }
    )
  }
}
