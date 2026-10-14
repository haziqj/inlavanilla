# Benchmark of integration rules for the VB mean correction, against a
# high-accuracy reference: tensor Gauss-Hermite quadrature for m <= 3, and
# 2^18 antithetic Monte Carlo draws otherwise. Reports, per model,
#   err = max_j |shift_rule_j - shift_ref_j| / sd_j
# and the number of gradient calls. Rules compared: the three-point
# Gauss-Hermite rule on each principal axis (2m + 1 nodes, INLAvaan's
# experimental option), INLAvaan's 64 and 128 Owen-scrambled Sobol points,
# and the fully symmetric degree-5 rule (2m^2 + 1 nodes) adopted by
# inlavanilla. In October 2026 the package moved to the degree-5 rule built on
# a regular simplex ((m + 1)(m + 2) + 1 nodes), which vb_nodes() now returns;
# it is the last column (see lab/vb-simplex-benchmark.R).
#
# Run from the package root with INLAVAAN set to a checkout of INLAvaan (for
# its stored Sobol table):
#   INLAVAAN=../INLAvaan Rscript lab/vb-benchmark.R
# Set VBTOL=1e-3 to count gradient calls at the package's default tolerance.
suppressMessages(devtools::load_all(".", quiet = TRUE))
load(file.path(Sys.getenv("INLAVAAN", "../INLAvaan"), "R", "sysdata.rda"))

ns <- asNamespace("inlavanilla")
pkg_nodes <- get("vb_nodes", ns)

# The fully symmetric degree-5 rule on the principal axes, 2m^2 + 1 nodes,
# which vb_nodes() returned before October 2026
deg5_nodes <- function(Sigma) {
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

# Three-point Gauss-Hermite rule on each principal axis (INLAvaan's
# vb_nodes_gauss_hermite())
gh_nodes <- function(Sigma) {
  m <- ncol(Sigma)
  e <- eigen(Sigma, symmetric = TRUE)
  axes <- sweep(e$vectors, 2, sqrt(pmax(e$values, 0)), "*")
  list(
    nodes = rbind(0, sqrt(3) * t(axes), -sqrt(3) * t(axes)),
    weights = c(1 - m / 3, rep(1 / 6, 2 * m))
  )
}

sobol_nodes <- function(n = 64L) {
  function(Sigma) {
    m <- ncol(Sigma)
    L <- t(chol(Sigma))
    z <- qnorm(SobolOwen[seq_len(n), seq_len(m), drop = FALSE])
    z <- sweep(z, 2, colMeans(z))
    list(nodes = rbind(0, z %*% t(L)), weights = rep(1 / (n + 1), n + 1))
  }
}

ref_nodes <- function(Sigma) {
  m <- ncol(Sigma)
  L <- t(chol(Sigma))
  if (m <= 3) {
    k <- if (m == 1) {
      80
    } else if (m == 2) {
      60
    } else {
      30
    }
    i <- seq_len(k - 1)
    b <- sqrt(i / 2)
    J <- matrix(0, k, k)
    J[cbind(i, i + 1)] <- b
    J[cbind(i + 1, i)] <- b
    e <- eigen(J, symmetric = TRUE)
    x <- e$values * sqrt(2)
    w <- e$vectors[1, ]^2
    grid <- as.matrix(expand.grid(rep(list(x), m)))
    wg <- apply(as.matrix(expand.grid(rep(list(w), m))), 1, prod)
    list(nodes = grid %*% t(L), weights = wg)
  } else {
    set.seed(42)
    z <- matrix(rnorm(2^17 * m), ncol = m)
    z <- rbind(z, -z)
    list(nodes = z %*% t(L), weights = rep(1 / nrow(z), nrow(z)))
  }
}

# vb_correct() uses the centre-difference form, which needs the first node
# to be the centre; prepend one with zero weight for the reference rule.
with_centre <- function(rule) {
  function(Sigma) {
    r <- rule(Sigma)
    list(nodes = rbind(0, r$nodes), weights = c(0, r$weights))
  }
}

run_rule <- function(rule, lp, gr, lap) {
  assignInNamespace("vb_nodes", rule, ns)
  on.exit(assignInNamespace("vb_nodes", pkg_nodes, ns))
  ctrl <- list(
    vb_maxit = 200L,
    vb_tol = as.numeric(Sys.getenv('VBTOL', '1e-9'))
  )
  n <- 0L
  grc <- function(x) {
    n <<- n + 1L
    gr(x)
  }
  out <- get("vb_correct", ns)(lp, grc, lap$mode, lap$Sigma, lap$R, lap$L, ctrl)
  list(shift = out$shift, ngrad = n, fallback = out$fallback)
}

models <- list()

# A. Normal, unknown mean and log sd, n = 10
set.seed(1)
yA <- rnorm(10, 1, 2)
models$normal_n10 <- list(
  init = c(mu = 0, ls = 0),
  lp = function(t) {
    s <- exp(t[2])
    sum(dnorm(yA, t[1], s, log = TRUE)) +
      dnorm(t[1], 0, 10, log = TRUE) +
      dexp(s, 0.5, log = TRUE) +
      t[2]
  },
  gr = function(t) {
    s <- exp(t[2])
    r <- yA - t[1]
    c(sum(r) / s^2 - t[1] / 100, -length(yA) + sum(r^2) / s^2 - 0.5 * s + 1)
  }
)

# B. Beta-binomial on the logit scale, y = 2 of n = 12, Beta(1, 1) prior
models$betabinom <- list(
  init = c(eta = 0),
  lp = function(t) {
    p <- plogis(t)
    dbinom(2, 12, p, log = TRUE) + dbeta(p, 1, 1, log = TRUE) + log(p * (1 - p))
  },
  gr = function(t) {
    p <- plogis(t)
    2 - 12 * p + (1 - 2 * p)
  }
)

# C. Poisson-gamma on the log scale, small counts
yC <- c(0, 1, 0, 2, 1)
models$poisgamma <- list(
  init = c(ll = 0),
  lp = function(t) {
    l <- exp(t)
    sum(dpois(yC, l, log = TRUE)) + dgamma(l, 1, 1, log = TRUE) + t
  },
  gr = function(t) {
    l <- exp(t)
    sum(yC) - length(yC) * l - l + 1
  }
)

# D/E. Logistic regression, m = 5 (n = 40) and m = 10 (n = 80)
logit_model <- function(m, n, seed) {
  set.seed(seed)
  X <- cbind(1, matrix(rnorm(n * (m - 1)), n))
  beta <- c(-0.5, rnorm(m - 1, 0, 0.8))
  yy <- rbinom(n, 1, plogis(X %*% beta))
  list(
    init = setNames(rep(0, m), paste0("b", seq_len(m))),
    lp = function(t) {
      eta <- as.numeric(X %*% t)
      sum(yy * eta - log1p(exp(eta))) + sum(dnorm(t, 0, 2.5, log = TRUE))
    },
    gr = function(t) {
      eta <- as.numeric(X %*% t)
      as.numeric(crossprod(X, yy - plogis(eta))) - t / 6.25
    }
  )
}
models$logit_m5 <- logit_model(5, 40, 3)
models$logit_m10 <- logit_model(10, 80, 4)
models$logit_m20 <- logit_model(20, 150, 6)
models$logit_m30 <- logit_model(30, 250, 7)
# Poisson regression with a heavier right tail, m = 8
set.seed(8)
Xp <- cbind(1, matrix(rnorm(30 * 7), 30))
yp <- rpois(30, exp(Xp %*% c(0.2, rnorm(7, 0, 0.3))))
models$pois_m8 <- list(
  init = setNames(rep(0, 8), paste0("b", 1:8)),
  lp = function(t) {
    eta <- as.numeric(Xp %*% t)
    sum(yp * eta - exp(eta)) + sum(dnorm(t, 0, 1, log = TRUE))
  },
  gr = function(t) {
    eta <- as.numeric(Xp %*% t)
    as.numeric(crossprod(Xp, yp - exp(eta))) - t
  }
)

# F. Weibull survival, log shape and log scale, n = 15
set.seed(5)
yF <- rweibull(15, shape = 1.5, scale = 2)
models$weibull <- list(
  init = c(lk = 0, ll = 0),
  lp = function(t) {
    k <- exp(t[1])
    l <- exp(t[2])
    sum(dweibull(yF, k, l, log = TRUE)) + sum(dnorm(t, 0, 2, log = TRUE))
  },
  gr = NULL
)

res <- list()
for (nm in names(models)) {
  mod <- models[[nm]]
  lp <- mod$lp
  gr <- mod$gr
  if (is.null(gr)) {
    gr <- function(x) get("fd_grad", ns)(lp, x, h = 1e-5)
  }
  ctrl <- list(nlminb = list(iter.max = 1000, eval.max = 2000))
  lap <- get("laplace_fit", ns)(lp, gr, mod$init, ctrl)
  sds <- sqrt(diag(lap$Sigma))
  ref <- run_rule(with_centre(ref_nodes), lp, gr, lap)
  gh <- run_rule(gh_nodes, lp, gr, lap)
  sob <- run_rule(sobol_nodes(64L), lp, gr, lap)
  sob16 <- run_rule(sobol_nodes(128L), lp, gr, lap)
  d5 <- run_rule(deg5_nodes, lp, gr, lap)
  sx <- run_rule(pkg_nodes, lp, gr, lap)
  err <- function(r) max(abs(r$shift - ref$shift) / sds)
  res[[nm]] <- data.frame(
    model = nm,
    m = length(mod$init),
    shift_ref = max(abs(ref$shift) / sds),
    err_gh = err(gh),
    err_sobol64 = err(sob),
    err_sobol128 = err(sob16),
    err_deg5 = err(d5),
    ngrad_deg5 = d5$ngrad,
    ngrad_gh = gh$ngrad,
    ngrad_sobol64 = sob$ngrad,
    err_simplex = err(sx),
    ngrad_simplex = sx$ngrad
  )
}
out <- do.call(rbind, res)
rownames(out) <- NULL
print(format(out, digits = 3))
