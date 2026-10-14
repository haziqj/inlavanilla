# Stage 1b: the variational Bayes mean correction. The Laplace Gaussian is
# recentred at mode + shift, keeping its covariance, where the shift maximises
# the Gaussian expectation of the log-posterior. At the optimum the expected
# gradient vanishes, E[grad log pi(mode + shift + z)] = 0 for z ~ N(0, Sigma).
#
# Ported from INLAvaan: the same fixed-point iteration with Anderson(1)
# acceleration and an nlminb() fallback. The Gaussian expectation uses a
# deterministic degree-5 cubature rule built on a regular simplex (see
# vb_nodes()) rather than INLAvaan's Sobol points: with a gradient, the
# expected gradient over its nodes; without one, Stein's identity on the same
# nodes, from log-posterior values alone.

# Node set for the Gaussian expectation: a degree-5 rule built on a regular
# simplex, on the principal axes of Sigma. In whitened coordinates its nodes
# are the centre and
#   +/- r a_k    at the m + 1 vertices a_k of a regular simplex (unit
#                vectors with a_k' a_l = -1/m),
#   +/- s b_kl   along the m (m + 1) / 2 edge directions
#                b_kl = (a_k + a_l) / |a_k + a_l|,
# (m + 1)(m + 2) nodes besides the centre, about half the 2m^2 of the fully
# symmetric degree-5 rule on the principal axes, which the package used
# before. The node set is invariant under the symmetries of the simplex and
# under -I, so by Sobolev's theorem the rule integrates every polynomial of
# degree five or less exactly once it integrates the invariant ones: 1,
# |z|^2, |z|^4 and p4(z) = sum_k (a_k' z)^4, whose Gaussian moments are 1, m,
# m(m + 2) and 3(m + 1). Given the radii, these fix the total weight of each
# shell. While the vertex shell's weight is not negative (m <= 7) both shells
# sit on the sphere of radius sqrt(m + 2), and every weight is positive.
# Beyond that the vertex weights are negative whatever the radii; r^2 = 3
# then keeps them small, and s^2 falls from 8.4 (m = 8) towards 6. At m = 7
# the vertex weight is exactly zero and those nodes are dropped. No node is
# more than 3 posterior SDs from the centre. Some weights are negative, so
# sums over the nodes are taken as differences, which keeps large terms from
# cancelling.
#
# For m <= 3 some directions coincide (a_2 = -a_1 at m = 1, a_k + a_l = -a_j
# at m = 2, a_1 + a_2 = -(a_3 + a_4) at m = 3), and each node is kept once.
# This gives the classical rules: the three-point Gauss-Hermite rule (m = 1),
# Radon's seven-point rule on a hexagon (m = 2), and the corners of a cube
# with the vertices of an octahedron (m = 3).
#
# The nodes come in antipodal pairs: the centre, then one node of each pair,
# then their negatives in the same order. The principal axes take a fixed
# sign (each with its largest entry positive), because the simplex is not
# symmetric under the sign of a single axis, so the rule must not depend on
# the signs that eigen() happens to return.
vb_nodes <- function(Sigma) {
  m <- ncol(Sigma)
  e <- eigen(Sigma, symmetric = TRUE)
  # Each axis with its largest entry positive
  top <- max.col(abs(t(e$vectors)), ties.method = "first")
  axes <- sweep(e$vectors, 2, sign(e$vectors[cbind(top, seq_len(m))]), "*")
  A <- sweep(axes, 2, sqrt(pmax(e$values, 0)), "*")
  rule <- simplex_rule(m)
  half <- rule$nodes[1 + seq_len(rule$n_pairs), , drop = FALSE] %*% t(A)
  list(nodes = rbind(0, half, -half), weights = rule$weights)
}

# The rule of vb_nodes() in whitened coordinates
simplex_rule <- function(m) {
  # Unit vertices: the rows of a Helmert contrast matrix, rescaled
  H <- matrix(0, m + 1, m)
  for (j in seq_len(m)) {
    H[seq_len(j), j] <- 1 / sqrt(j * (j + 1))
    H[j + 1, j] <- -j / sqrt(j * (j + 1))
  }
  vertices <- H * sqrt((m + 1) / m)
  p4 <- function(x) sum(as.numeric(vertices %*% x)^4)
  # One vertex and one edge direction per antipodal pair
  V <- if (m == 1) vertices[1, , drop = FALSE] else vertices
  edges <- if (m >= 4) {
    utils::combn(m + 1, 2)
  } else if (m == 3) {
    rbind(1, 2:4)
  } else {
    matrix(0L, 2, 0)
  }
  B <- vertices[edges[1, ], , drop = FALSE] +
    vertices[edges[2, ], , drop = FALSE]
  B <- B / sqrt(rowSums(B^2))

  # The weights of the two shells, W1 = a4 / r^4 and W2 = b4 / s^4
  if (nrow(B) == 0) {
    r2 <- s2 <- m + 2
    a4 <- m * (m + 2)
    b4 <- 0
  } else {
    alpha <- p4(V[1, ])
    beta <- p4(B[1, ])
    a4 <- (3 * (m + 1) - beta * m * (m + 2)) / (alpha - beta)
    b4 <- m * (m + 2) - a4
    if (a4 > -1e-8) {
      r2 <- s2 <- m + 2
    } else {
      r2 <- 3
      s2 <- b4 / (m - a4 / r2)
    }
  }
  w1 <- a4 / r2^2 / (2 * nrow(V))
  w2 <- if (nrow(B) > 0) b4 / s2^2 / (2 * nrow(B)) else numeric()
  half <- rbind(sqrt(r2) * V, sqrt(s2) * B)
  w <- c(rep(w1, nrow(V)), rep(w2, nrow(B)))
  keep <- abs(w) > 1e-12
  half <- half[keep, , drop = FALSE]
  w <- w[keep]
  list(
    nodes = rbind(0, half, -half),
    weights = c(1 - 2 * sum(w), w, w),
    n_pairs = nrow(half)
  )
}

vb_correct <- function(lp, gr, mode, Sigma, R, L, ctrl) {
  m <- length(mode)
  stein <- is.null(gr)
  rule <- vb_nodes(Sigma)
  Z <- rule$nodes
  w <- rule$weights
  nn <- nrow(Z)
  sds <- sqrt(diag(Sigma))

  # The log-posterior at every node but the centre. Without a gradient the
  # values of the last call are kept, so that the fallback's objective and
  # gradient at one point share them.
  at_shift <- NULL
  at_values <- NULL
  node_lp <- function(shift) {
    if (!stein || !identical(shift, at_shift)) {
      mu <- mode + shift
      v <- numeric(nn - 1)
      for (b in 2:nn) {
        v[b - 1] <- lp(mu + Z[b, ])
      }
      at_shift <<- shift
      at_values <<- v
    }
    at_values
  }

  # The fixed-point step Sigma E[grad log pi(mode + shift + z)] at a shift.
  # With a gradient, the expected gradient is the weighted mean over the
  # nodes, written as differences from the centre gradient: at most
  # (m + 1)(m + 2) + 1 gradient calls per sweep, and exact when the
  # log-posterior is a polynomial of degree six or less in whitened
  # coordinates. Without one it comes from Stein's identity,
  # E[z f(mu + z)] = Sigma E[grad f(mu + z)] for z ~ N(0, Sigma), on the same
  # nodes: the step is the weighted sum of the nodes times the log-posterior
  # there, taken over antipodal pairs as w_b z_b (f(mu + z_b) - f(mu - z_b)),
  # which cancels the size of the log-posterior. That takes at most
  # (m + 1)(m + 2) log-posterior calls per sweep, against m + 1 times as many
  # for the gradient form on finite differences, and z f(mu + z) has degree
  # five when the log-posterior has degree four, so the shift is then
  # exact.
  if (stein) {
    n_pairs <- (nn - 1) %/% 2
    plus <- seq_len(n_pairs)
    Zp <- Z[1 + plus, , drop = FALSE]
    wp <- w[1 + plus]
    step_at <- function(shift) {
      v <- node_lp(shift)
      as.numeric(crossprod(Zp, wp * (v[plus] - v[n_pairs + plus])))
    }
  } else {
    score <- function(shift) {
      mu <- mode + shift
      g0 <- gr(mu)
      s <- g0
      for (b in 2:nn) {
        s <- s + w[b] * (gr(mu + Z[b, ]) - g0)
      }
      s
    }
    step_at <- function(shift) as.numeric(Sigma %*% score(shift))
  }

  # Negative expected log-posterior, for the fallback optimiser. A failed node
  # fails the objective, since the negative centre weight would otherwise
  # reward a -Inf there.
  failed <- 1e40
  objective <- function(shift) {
    lps <- c(lp(mode + shift), node_lp(shift))
    if (!all(is.finite(lps))) {
      return(failed)
    }
    -(lps[1] + sum(w[-1] * (lps[-1] - lps[1])))
  }

  # The fallback can end with no mean at which logpost is finite at every
  # node, or at the edge of the region where it is, with the optimum beyond.
  # The rule does not define the correction there, and a mean on that edge
  # would only show where the nodes are, so the fit stops. The message names
  # a node around `mode + shift` where logpost is -Inf.
  stop_failed_nodes <- function(shift, pinned) {
    mu <- mode + shift
    lps <- vapply(seq_len(nn), function(b) lp(mu + Z[b, ]), 0)
    x <- mu + Z[which.max(!is.finite(lps)), ]
    at <- paste(names(mode), "=", signif(x, 4), collapse = ", ")
    reach <- with_sd(max(sweep(abs(Z), 2, sds, "/")))
    abort(
      c(
        "The VB mean correction cannot be computed.",
        "x" = if (pinned) {
          "It stopped at the edge of the region where {.arg logpost} is
           finite at all of its nodes, and its next step puts a node at
           {at}, where {.arg logpost} is {.code -Inf}."
        } else {
          "{.arg logpost} is {.code -Inf} at some of its nodes around every
           mean tried, for example at {at}."
        },
        "i" = "The nodes reach {reach} from the mean, and the correction
               needs {.arg logpost} finite at all of them.",
        "i" = "Reparameterise any bounded parameter so that every real value
               is allowed (for example its log or logit, with the
               log-Jacobian), or set {.code vb = FALSE}."
      ),
      call = quote(inla_vanilla())
    )
  }

  # Fixed-point iteration shift <- shift + step_at(shift): with the curvature
  # of the Laplace approximation, the Newton step for the stationarity
  # condition is one multiplication by Sigma. Convergence is judged on the
  # step in posterior SDs. A step over `maxstep` SDs, or a non-finite step,
  # means the quadratic premise has failed (a strongly skewed posterior), and
  # the solve falls back to nlminb().
  # The history (the step at each iteration, in posterior SDs, and whether
  # the Anderson secant replaced it) and the reason for a fallback are kept
  # for verbose output.
  maxstep <- 1
  shift <- numeric(m)
  iter <- 0L
  move <- Inf
  fallback <- FALSE
  rerun <- FALSE
  reason <- NULL
  steps <- numeric()
  accelerated <- logical()
  step_prev <- shift_prev <- NULL
  for (it in seq_len(ctrl$vb_maxit)) {
    step <- step_at(shift)
    if (!all(is.finite(step))) {
      fallback <- TRUE
      reason <- if (stein) {
        "the log-posterior was not finite at a node"
      } else {
        "the expected gradient was not finite"
      }
      break
    }
    move <- max(abs(step) / sds)
    steps[it] <- move
    if (move > maxstep) {
      fallback <- TRUE
      reason <- "a step was over 1 posterior SD"
      break
    }
    if (move < ctrl$vb_tol) {
      shift <- shift + step
      iter <- it
      break
    }
    # Anderson(1) acceleration: a secant estimate from consecutive steps
    # removes the dominant error mode. The plain step is kept whenever the
    # secant is degenerate or would move too far.
    update <- step
    accelerated[it] <- FALSE
    if (!is.null(step_prev)) {
      df <- step - step_prev
      dx <- shift - shift_prev
      denom <- sum(df^2)
      if (denom > 0) {
        gam <- sum(step * df) / denom
        cand <- step - gam * (dx + df)
        if (
          all(is.finite(cand)) &&
            max(abs(cand) / sds) <= min(maxstep, 2 * move)
        ) {
          update <- cand
          accelerated[it] <- TRUE
        }
      }
    }
    step_prev <- step
    shift_prev <- shift
    shift <- shift + update
    iter <- it
  }
  if (move >= ctrl$vb_tol) {
    fallback <- TRUE
    if (is.null(reason)) {
      reason <- "the iteration did not converge in `vb_maxit` steps"
    }
  }

  if (fallback) {
    # Optimise in whitened coordinates, where the problem is well conditioned.
    # The gradient there is -L' E[grad], which without a gradient is
    # -L' Sigma^-1 step = -R step. That is the expected gradient by Stein's
    # identity, which agrees with the gradient of the objective to the error of
    # the rule, and costs no calls beyond the objective's at the same point.
    #
    # The gradient is not finite where logpost is -Inf at a node, as the
    # objective is failed there, or where grad is not finite at one.
    # nlminb() asks for it only at its start and at the points it accepts,
    # which have a finite objective, but the start, the mode, has failed
    # nodes when logpost is -Inf within their reach. There the gradient is
    # the mean of the failed nodes in whitened coordinates, so that nlminb()
    # moves the mean away from them, back towards the region where every node
    # is finite. Finding the failed nodes takes no calls without a gradient,
    # and one more sweep with one.
    U <- Z %*% t(R) # whitened nodes, since L = R^-1
    away <- function(shift) {
      bad <- if (stein) {
        c(FALSE, !is.finite(node_lp(shift)))
      } else {
        mu <- mode + shift
        vapply(seq_len(nn), function(b) !all(is.finite(gr(mu + Z[b, ]))), NA)
      }
      if (any(bad)) colMeans(U[bad, , drop = FALSE]) else numeric(m)
    }
    fall_back <- function(start) {
      stats::nlminb(
        start = start,
        objective = function(d) objective(as.numeric(L %*% d)),
        gradient = function(d) {
          shift <- as.numeric(L %*% d)
          g <- if (stein) {
            -as.numeric(R %*% step_at(shift))
          } else {
            as.numeric(crossprod(L, -score(shift)))
          }
          if (all(is.finite(g))) g else away(shift)
        },
        control = list(rel.tol = 1e-8)
      )
    }
    nl <- fall_back(numeric(m))
    shift <- as.numeric(L %*% nl$par)
    iter <- nl$iterations
    if (nl$objective >= failed) {
      stop_failed_nodes(shift, pinned = FALSE)
    }
    # Where the fallback landed: the fixed-point step there, the Newton step
    # for the stationarity condition. One more sweep.
    step <- step_at(shift)
    # With a gradient the objective and its gradient agree, so a landing that
    # did not settle means nlminb() stopped short. It can report relative
    # convergence well short of the optimum after starting where a node
    # fails (the objective's 1e40 there). Rerun it once from where it landed,
    # as the mode search reruns nlminb() (R/laplace.R). Without a gradient a
    # rerun cannot help: its gradient, Stein's estimate, disagrees with the
    # objective by the rule's error (lab/DESIGN.md). Only an unsettled
    # landing with a gradient pays for it.
    rerun <- !stein &&
      all(is.finite(step)) &&
      max(abs(step) / sds) >= warn_levels$vb_step
    if (rerun) {
      nl <- fall_back(nl$par)
      shift <- as.numeric(L %*% nl$par)
      iter <- iter + nl$iterations
      step <- step_at(shift)
    }
    # A landing that did not settle may be pinned to the edge of the region
    # where every node is finite, with the optimum beyond it. nlminb() then
    # stops as close to the edge as its tolerance allows, so a thousandth of
    # the step still to go leads out of the region. A landing short of an
    # interior optimum is not that close to the edge. One more sweep, only
    # when the landing did not settle.
    if (
      all(is.finite(step)) &&
        max(abs(step) / sds) >= warn_levels$vb_step &&
        objective(shift + step / 1000) >= failed
    ) {
      stop_failed_nodes(shift + step / 1000, pinned = TRUE)
    }
  } else if (max(abs(step) / sds) >= warn_levels$vb_step) {
    # The iteration met a vb_tol above the settling level, so the step just
    # taken says little about the one still to go: take it. One more sweep,
    # which the default vb_tol never needs.
    step <- step_at(shift)
  }
  # The step still to go at the VB mean, per parameter in posterior SDs, as
  # the Newton step at the mode is for the mode: the step where the fallback
  # landed, or on the fixed-point path the last step taken, which was below
  # vb_tol and bounds the next (taken afresh when it was above the settling
  # level). A step that is not finite is Inf.
  step_sigma <- abs(step) / sds
  step_sigma[!is.finite(step_sigma)] <- Inf

  names(shift) <- names(mode)
  list(
    shift = shift,
    iterations = iter,
    fallback = fallback,
    step_sigma = stats::setNames(step_sigma, names(mode)),
    history = list(
      steps = steps,
      accelerated = accelerated,
      reason = reason,
      nlminb = if (fallback) nl$message,
      rerun = fallback && rerun
    ),
    # Joint KL divergence between the Laplace Gaussian and its shifted copy,
    # 0.5 shift' H shift. It is the location term the Laplace estimate of the
    # marginal likelihood misses, and is added to it.
    kld_global = 0.5 * sum((R %*% shift)^2)
  )
}
