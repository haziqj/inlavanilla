# Results

Run on 8 November 2026, before the method was added to the package. The
criterion is the worst-parameter error in the posterior mean, in posterior SDs,
against quadrature (one or two parameters) or HMC. Monte Carlo standard errors
of the HMC means are 0.01 to 0.03 SDs.

## Accuracy and calls

Median over 15 datasets, except Political Democracy (real data, cached HMC). The
axis pass is guarded at 2 SDs here. "Now" is the iterative solver.

| Model | m | Now | Two passes | Skew-normal mean | Axis pass | Calls now | Calls, axis |
|---|---|---|---|---|---|---|---|
| CFA, n = 200 | 12 | 0.030 | 0.030 | 0.019 | 0.027 | 549 | 25 |
| CFA, n = 60 | 12 | 0.189 | 0.186 | 0.118 | 0.166 | 732 | 25 |
| Student-t regression | 5 | 0.092 | 0.093 | 0.032 | 0.049 | 172 | 11 |
| Poisson, low counts | 15 | 0.026 | 0.023 | 0.139 | 0.017 | 1,092 | 31 |
| Logistic, n = 80 | 10 | 0.103 | 0.095 | 0.097 | 0.032 | 532 | 21 |
| Logistic, n = 150 | 20 | 0.185 | 0.185 | 0.268 | 0.046 | 2,315 | 41 |
| Logistic, n = 600 | 20 | 0.029 | 0.029 | 0.017 | 0.017 | 1,389 | 41 |
| Logistic, n = 250 | 30 | 0.252 | 0.266 | 0.585 | 0.066 | 4,965 | 61 |
| Random-intercept logistic | 13 | 2.38 | 2.48 | 13.7 | 2.38 | 1,899 | 27 |
| Political Democracy SEM | 39 | 0.865 | 0.865 | 0.896 | 0.859 | 6,564 | 79 |

Across the 133 simulated datasets, the guarded axis pass against the iterative
solver: better by more than 0.01 SDs on 82, tied on 48, worse on 3 (largest
loss 0.14). The 2 SD guard fired on 13 (10 random-intercept, 3 logistic).

On Political Democracy the worst parameter is off for every route. The RMS
error over parameters is 0.156 for the axis pass and 0.169 for the solver.

Without a gradient, one Stein pass against the converged Stein solver: better on
98, tied on 27, worse on 8, at about a quarter of the calls.

Two passes stay within 0.02 SDs of the solver on every regular model, so
"integrate twice" reaches VB, but no closer to the posterior mean.

## One dimension

Log-gamma posterior, positions in Laplace SDs from the mode.

| a | Step 1 | Step 2 | Exact VB | True mean |
|---|---|---|---|---|
| 1 | -0.638 | -0.504 | -0.500 | -0.577 |
| 2 | -0.400 | -0.353 | -0.354 | -0.382 |
| 5 | -0.235 | -0.223 | -0.224 | -0.231 |
| 20 | -0.113 | -0.112 | -0.112 | -0.113 |

## Where one step needs help

| Posterior | Largest first step | Unguarded error | Solver | Capped axis iteration |
|---|---|---|---|---|
| Boundary (a + 1 ~ Gamma(4, 1)) | 1.7e8 SDs | 1.7e8 | 0.072 | 0.072 |
| Skew-normal, alpha = 10 | 2.9 SDs | 1.435 | 0.159 | 0.159 |
| Eight schools (non-centred) | 12.3 SDs | 10.2 | 0.975 | 2.53 |

Skew-normal, alpha = 5: 0.056 for the axis pass against 0.023 for the solver,
a small loss on strong skew in one dimension.

## Cost

VB time as a multiple of the mode and Hessian stage, with a gradient.

| Model | m | Solver | One axis pass |
|---|---|---|---|
| mtcars logistic | 3 | 1.21 | 0.10 |
| Logistic | 10 | 3.85 | 0.15 |
| CFA | 12 | 3.71 | 0.12 |
| Random-intercept logistic | 13 | 21.5 | 0.36 |
| Poisson | 15 | 5.30 | 0.12 |
| Logistic | 20 | 10.7 | 0.18 |
| Logistic | 30 | 16.3 | 0.25 |
| Political Democracy SEM | 39 | 12.5 | 0.15 |
| Logistic | 40 | 21.5 | 0.34 |

Without a gradient the solver runs at 0.75 to 2.4 times that stage, and one
Stein pass at 0.13 to 0.58.

More spending did not help: at m = 20, the worst-parameter error is 0.046 after
one axis pass and 0.225 after two. One pass of the degree-5 rule is further off
as well, and over budget above m of about 10.

## Log marginal likelihood

It changes by a median of 0.04 nats (90th percentile 0.41) between the capped
axis pass and the solver. The axis pass was closer to the importance-sampling
truth in 7 of 9 families, and 0.135 nats worse in the worst one.

## After implementation

`inla_vanilla()` itself, with the guard at 2.5 SDs and Gaussian marginals. VB
calls are the gradient calls over those of `vb = "off"`. The budget is all the
calls of the `vb = "off"` fit.

| Model | m | One-step | Full | Budget | Stopped at step 1 | Largest gap (SDs) |
|---|---|---|---|---|---|---|
| mtcars logistic | 3 | 7 | 60 | 92 | yes | 0.041 |
| Logistic | 20 | 41 | 2,315 | 478 | yes | 0.231 |
| Logistic | 30 | 61 | 3,972 | 622 | yes | 0.128 |
| Logistic | 40 | 81 | 5,169 | 768 | yes | 0.041 |
| Poisson | 15 | 31 | 1,365 | 352 | yes | 0.045 |
| CFA, n = 60 | 12 | 25 | 732 | 306 | yes | 0.026 |
| Random-intercept logistic | 13 | 27 | 1,688 | 345 | yes | 1.111 |
| Eight schools | 10 | 294 | 2,926 | 296 | no | 0.899 |
| Skew-normal, alpha = 10 | 1 | 18 | 18 | 43 | no | 0.000 |
| Political Democracy SEM | 39 | 79 | 6,564 | 921 | yes | 0.154 |

On the two hard models with large gaps (random-intercept logistic, eight
schools), neither method is close to the posterior mean: the errors are 1 SD or
more.

Larger logistic regressions (n = 10m), total time of the fit with Gaussian
marginals:

| m | `"off"` | `"onestep"` | `"full"` | VB gradients, one-step | VB gradients, full | Largest gap (SDs) |
|---|---|---|---|---|---|---|
| 80 | 0.10 s | 0.10 s | 1.90 s | 161 | 19,929 | 0.049 |
| 120 | 0.22 s | 0.27 s | 8.99 s | 241 | 44,289 | 0.072 |

`loo()` on logistic regressions, elpd under each method.

| Model | m | One-step | Full | Difference | SE of elpd |
|---|---|---|---|---|---|
| mtcars | 3 | -7.571 | -7.621 | 0.051 | 2.36 |
| Logistic | 5 | -26.649 | -26.656 | 0.006 | 3.65 |
| Logistic | 10 | -43.218 | -43.247 | 0.029 | 6.68 |
| Logistic, near separation | 20 | -46.198 | -46.768 | 0.570 | 9.57 |
| Logistic | 30 | -108.642 | -108.785 | 0.142 | 12.10 |

The two are compared with each other here, not with brute-force leave-one-out.

## ECPE against Stan

The DINA model of `vignettes/articles/cdm.qmd`, 63 parameters, no gradient, so
one-step is one Stein pass. Largest gap in the posterior mean from Stan, in
Stan's SDs, with Gaussian marginals centred on each mean.

| Block | Laplace | One-step | Full |
|---|---|---|---|
| Guessing | 0.256 | 0.102 | 0.121 |
| Slipping | 0.189 | 0.048 | 0.071 |
| Class proportions | 0.622 | 0.238 | 0.191 |

The VB stage takes 2.1 s against 14.6 s. The default fit takes 22.9 s and 45,148
calls against 35.7 s and 70,108. The skew-normal marginals do not depend on VB.

## Where one step is worse: diagnostic tests

Three latent class models of `vignettes/articles/diagnostic-tests.qmd`, no
gradient, against long Metropolis runs (means) and importance sampling (log
Z). Worst-parameter error in the mean, in SDs.

| Model | m | One-step | Full | Axis rule, numerical gradient |
|---|---|---|---|---|
| Conditionally independent | 11 | 0.16 | 0.18 | 0.14 |
| Random effects (LCRE) | 13 | 0.92 | 0.83 | 0.82 |
| Finite mixture | 13 | 0.36 | 0.31 | 0.30 |

On the random effects model the true shifts of the log SDs are 1.25 and 1.43
SDs, and every method is off by about 0.8 or more somewhere. One-step gets those
two right to 0.27, and is worse on the sensitivities and specificities (RMS over
parameters 0.50, against 0.36 for the full solve).

Log marginal likelihood, error against importance sampling (SE 0.004 to 0.022).

| Model | Laplace | One-step | Full |
|---|---|---|---|
| Conditionally independent | -0.35 | -0.08 | -0.14 |
| Random effects (LCRE) | -0.83 | +1.20 | -0.31 |
| Finite mixture | -0.70 | -0.28 | -0.32 |

The VB term is a leading-order term, and a first step of 1.5 SDs overstates it.
The article compares these models by their log marginal likelihoods, so it uses
`vb = "full"`. Differences from the random effects model, with its mirror image
counted: truth -0.49 and -2.66, full -0.33 and -2.68, one-step -1.77 and -4.13.

## Why one step fails, and the gain ratio

`why-one-step-fails.R`: 160 fits in 36 model families, each against the true
mean and log Z. One step assumes the expected log-posterior is a parabola with
the Laplace curvature. The gain ratio (`vb_gain_ratio()`) is the rise delivered
over the rise that parabola promises, on `axis_nodes()`, 4m + 2 calls of
`logpost`. Its median is 0.85, and it is flagged below 0.5.

Three ways to fail:

1. The step overshoots, because the parabola is wrong: boundary, eight schools,
   strongly skewed posteriors, LCRE (predicted +1.81, delivered -1.42). The full
   solve does better.
2. The posterior is not Gaussian-shaped: the random-intercept logistic, whose
   true mean is about 4 SDs from the mode. The ratio is about 0.1, and the full
   solve is worse.
3. A far tail the mode cannot see: the log variances of the SEMs. The parabola
   holds near the peak, the ratio is healthy, and only the mean gap fires.

As shipped, with the guard. "Caught" counts the fits with that failure that the
check flagged.

| Failure | Gradient: fits | Gain ratio caught | Mean gap > 0.25 caught | No gradient: fits | Gain ratio caught | Mean gap > 0.25 caught |
|---|---|---|---|---|---|---|
| Mean error above 0.2 SD | 18 | 11 | 14 | 19 | 10 | 14 |
| Mean error above 0.5 SD | 7 | 6 | 7 | 3 | 2 | 3 |
| Worse than the full solve by 0.05 SD | 4 | 2 | 3 | 4 | 3 | 4 |
| Log marginal likelihood off by over 1 nat | 11 | 1 | 10 | 19 | 4 | 15 |

The counts are of the fits that stopped at one step (151 with a gradient, 147
without). The rest went on to iterate, as the guard asked, and have no gain
ratio.

The gain ratio raised no false alarms: all 11 flagged fits with a gradient, and
all 10 without one, had a mean error above 0.2 SD. It fired on none of the 133
and 128 fits below that. The mean gap above 0.25 SD fires on 24 and 17 of them.
So the ratio is precise and the mean gap is sensitive, and neither replaces the
other. Most flags come from one family (random-intercept logistic: 9 of 11 with
a gradient, 6 of 10 without).

Rival checks, before the guard, on all 160 fits: the curvature along the step
above 1 caught 4 of the 16 mean errors above 0.5 SD. A step above 1 SD caught 15
of them but fired on 48 others.

Acting on the ratio does not pay. Using the full solve on the flagged fits helped
on 5 and hurt on 15 (all random-intercept), so the check warns and leaves the
fit alone.

Thresholds were chosen on these same fits, and 15 of the 36 families have 15
datasets each while the others have one.

## Not tested

`loo()` against brute-force leave-one-out.
