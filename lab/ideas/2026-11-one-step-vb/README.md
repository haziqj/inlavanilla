# One-step VB correction

November 2026. The VB mean correction was the stage that made fits slow, at 1.2
to 21 times the mode and Hessian stage with a gradient. The owner's budget is
that VB should cost no more than that stage, and be as accurate as possible
within it. This idea asked whether the solve is needed at all.

## The idea

Minimising the KL divergence over the mean of a Gaussian with the Laplace
variance reduces to one equation, E[slope of the log-posterior] = 0 under that
Gaussian. In standardised units the Gaussian part of the slope averages to
exactly minus the shift, so the equation is d = phi(d), with phi the average of
the non-Gaussian part. The function phi is nearly flat (its slope is O(1/n)), so
the first step, d1 = phi(0), holds most of the movement. Its leading term is
half the skewness, and the error of one step against exact VB is O(n^-3/2) in
posterior SDs.

Exact VB is not the posterior mean: the two differ by O(n^-3/2) SDs as well. So
iterating to convergence exchanges one small remainder for another.

The cheap version of that step uses the three-point Gauss-Hermite rule along
each principal axis of Sigma: 2m + 1 gradients, about one numerical Hessian.

## The coauthor's remarks

- "The skew is quadratic in the mean": phi is quadratic in d up to O(n^-5/2).
- "Integrate twice": two passes reach exact VB to O(n^-5/2). True, and not
  needed for the posterior mean.
- "Lookup table", "once and for all": Dutta's (2024) h(s), the KL divergence
  between a skew-normal and a Gaussian. It is the value of the minimum, s^2/12
  to leading order, not its location, so it does not enter a mean-only solve.

## Decision

`vb = "onestep"` is the default, with `vb = "full"` keeping the iterative
solver and `vb = "off"` (or `FALSE`) leaving the correction out. Without a
gradient, one-step VB is one Stein pass on the degree-5 rule, since Stein's
identity needs that rule's degree. One step is accepted unless a node fails or
the step exceeds 2.5 SDs. In that case the same rule is iterated with capped
steps, within the budget.

The 2.5 SD level is a heuristic. At 2 SDs the guard fired on three large
logistic fits where the plain step was mostly better. At 3 SDs it missed the
strongly skewed skew-normal (a first step of 2.9 SDs).

## Files

- `common.R`: models, ground truth, candidate routes, generators.
- `routes-benchmark.R`: 24 posteriors, one dataset each.
- `routes-simulation.R`: nine families, 15 datasets each.
- `guard-and-budget.R`: guard levels, and cost against the budget.
- `RESULTS.md`: the numbers.

The scripts compare routes outside the package and call `vb_correct()` with its
default `method = "full"`, so they keep measuring the iterative solver.
