# Models R-INLA cannot fit: showcase candidates

Collected in October 2026, to be picked up later. The question: which models
do people want to fit in R-INLA but cannot, and can `inla_vanilla()` fit them?

## Positioning

inlavanilla will not beat R-INLA on a latent Gaussian model: its sparse C code,
its integration over the hyperparameters and its latent-field marginals are
years of engineering. Its case is the models outside that class, which R-INLA
cannot express:

- a nonlinear predictor (a parameter times a latent variable, a growth or
  concentration curve);
- a non-Gaussian latent variable (exponential inefficiency, discrete classes);
- a non-standard likelihood (sample selection, misclassification, mixtures);
- random effects on a dispersion or variance parameter;
- too many hyperparameters for R-INLA's integration.

Against general-purpose tools the competitor is Stan, and the offer is seconds
instead of minutes, deterministic results, and a log marginal likelihood and
LOO from one fit. With TMB or RTMB integrating out many random effects by
Laplace, the use case is prototyping: write the model in an afternoon, see
whether it is worth building properly.

How the latent variables are handled decides the route:

1. **Closed form or a sum** inside `logpost` (the stochastic frontier below,
   occupancy models, the DCM article).
2. **Quadrature** inside `logpost`, for one or two dimensions per unit (the
   binary factor article).
3. **Parameters**, when there are few of them (but see the Theoph finding:
   integrating them out is more accurate).
4. **RTMB's Laplace** over many random effects, leaving a modest number of
   outer parameters (up to about 30, since the calls grow with the square of
   the count and each is an inner optimisation).

## Prototyped

### Stochastic production frontier (`lab/frontier-exponential.R`)

- **Thread:** R-INLA discussion group, 7-8 October 2020, "Stochastic
  production frontiers",
  <https://groups.google.com/g/r-inla-discussion-group/c/7rd0KAncAlo> (read
  in full). The poster had the model in Stan. Håvard Rue: "We cannot do an
  exponential random effect like this." Finn Lindgren offered an inlabru
  copula trick, u = qexp(pnorm(w)), without the full INLA approximation.
- **Model:** $\log y_i = x_i^\top\beta - u_i + v_i$, $v_i \sim N(0,
  \sigma^2)$, $u_i \sim \text{Exp}(\lambda)$, efficiency $e^{-u_i}$. Priors as
  in the post: $\beta \sim N(0, 100^2)$, $\sigma \sim
  \text{InvGamma}(0.001, 0.001)$, $\lambda \sim \text{Exp}(-\log 0.875)$ (van
  den Broeck et al., 1994).
- **Route:** closed form. $p(\epsilon_i) = \lambda \exp(\lambda
  \epsilon_i + \lambda^2\sigma^2/2)\, \Phi(-\epsilon_i/\sigma -
  \lambda\sigma)$ (Aigner, Lovell & Schmidt, 1977), 6 parameters. Given
  $\epsilon_i$, $u_i$ is $N(-\epsilon_i - \lambda\sigma^2, \sigma^2)$
  truncated at 0, so the efficiencies come from `sampling(fun = ...)`.
- **Data:** `frontier::riceProdPhil`, 344 rice harvests (43 Philippine farms,
  8 years), Cobb-Douglas in area, labour and fertiliser. The thread had none.
- **Result:** matches the poster's Stan model (344 $u_i$ sampled): parameters
  within 0.03 SDs in mean and 0.07 in the 95% limits, SDs within 2.5%;
  efficiencies within 0.006 in mean and 0.02 in the 95% limits. 0.6 s against
  21 s of sampling and a minute of compiling. Every check passes.
- **Next:** a website article (Examples). Possible extensions: the panel
  version with one time-invariant $u_i$ per farm (also closed form), and
  exponential against half-normal inefficiency by log marginal likelihood (the
  frontier package fits the half-normal by ML).

### Random effects on the precision (`lab/dhglm-sleepstudy.R`)

- **Source:** Morales-Otero, Gómez-Rubio and Núñez-Antón (2022), "Fitting
  double hierarchical models with the integrated nested Laplace
  approximation", *Statistics and Computing*, 32(4), arXiv 2201.07493. Code
  and data: <https://github.com/becarioprecario/DHGLM-INLA> (folders
  `sleepstudy`, `sleepstudy_rcoef`, `Gaussian5`, `Poisson`, `NegBin`, each
  with a JAGS script). The abstract: DHGLMs do not fit within the models INLA
  can fit. Their workaround, INLA inside adaptive multiple importance sampling
  over the 18 log-precisions, took 21 rounds of 1,000 INLA fits. A first,
  vague proposal gave results they discarded. The second was centred at each
  subject's $\log(1/S_i^2)$, permuted 500 times to pick the start with the
  highest conditional marginal likelihood. The final effective sample size was
  2.0.
- **Model** (their equation 3, reaction time in seconds): $y_{ij} \sim
  N(\beta_0 + \beta_i\,\text{day}_{ij}, 1/\tau_i)$, $\log\tau_i = \gamma +
  u_i$, $\beta_i \sim N(0, 1/\tau_\beta)$, $u_i \sim N(0, 1/\tau_u)$, with
  $\beta_0, \gamma \sim N(0, \text{precision } 0.001)$, $\tau_\beta \sim
  \text{Gamma}(1, 0.00005)$ and $\tau_u \sim \text{Gamma}(0.01, 0.01)$. The
  paper's text gives $\text{Gamma}(1, 0.00005)$ for $\tau_u$, but its JAGS and
  INLA code use $\text{Gamma}(0.01, 0.01)$. The random slopes have mean 0, so
  they carry the average effect of the days too.
- **Route:** both random effects inside `logpost`. Given $\tau_i$, subject
  $i$'s responses are normal with covariance $I/\tau_i + d d^\top/\tau_\beta$,
  a rank-one update, so $\beta_i$ integrates out in closed form. Then
  $v_i = \log\tau_i$ integrates out by adaptive Gauss-Hermite (15 nodes,
  centred at the integrand's mode by Newton steps, scaled by its curvature).
  A fixed rule fails when $\tau_u$ is small. 4 outer parameters. Subject
  effects from `sampling(fun = ...)`: $v_i$ from its conditional density on a
  grid, then $\beta_i$ from its normal.
- **Data:** `lme4::sleepstudy`, 18 subjects, 10 days each.
- **Result:** matches Stan (non-centred, all 40 parameters). The 4 outer
  parameters are within 0.05 SDs in mean and 0.17 in the 95% limits, with
  SDs within 5%. The 36 subject effects are within 0.03 SDs in mean and 4% in
  SD. $\tau_u$ on its natural scale: mean 2.59 against Stan's 2.72 (the
  paper: 2.12 by AMIS, 2.66 by JAGS), SD 1.9 against 3.1, driven by its long
  right tail (the paper: 0.93 and 2.3). 0.27 s (389 calls) and 5 s for the
  subject draws. Stan samples this small model in 4.5 s too, so the contrast
  is with the paper's 21,000 INLA fits. Every check passes, and 41 nodes move
  `logpost` at the mode by 2e-8.
- **Next:** the strongest candidate for an article, since the authors are
  from the INLA community. Add the natural model with a fixed slope,
  $\beta_0 + \beta_1\,\text{day} + b_i\,\text{day}$, and compare the two by
  log marginal likelihood. The repo's other examples (`Gaussian5`, `Poisson`,
  `NegBin`) have JAGS references too.

### Population pharmacokinetics (`lab/theoph-rtmb.R`)

A modern application: pharmacometrics runs on NONMEM, Monolix and nlmixr2,
whose main method (FOCE) is a Laplace approximation over the subject effects.
Bayesian fits today mean Stan (Torsten, or brms's nonlinear formulas), which
takes minutes to hours. R-INLA cannot fit these models, since the random
effects enter the curve nonlinearly (inlabru's linearisation is meant for
mildly nonlinear predictors). Related threads, neither with data: a nonlinear
mixed model $y = a(1 - e^{-b^x})$ (2 November 2023,
<https://groups.google.com/g/r-inla-discussion-group/c/n6cLcss2jNI>, "cannot
work in pure inla", pointed to inlabru) and a logistic curve with random
effects on its location and scale (8 April 2021,
<https://groups.google.com/g/r-inla-discussion-group/c/YD_kJKfhc0o>).

- **Model:** one compartment, first-order absorption,
  $\mu_{ij} = \frac{D_i k_e k_{a,i}}{Cl_i (k_{a,i} - k_e)}(e^{-k_e t_{ij}} -
  e^{-k_{a,i} t_{ij}})$, $y_{ij} \sim N(\mu_{ij}, \sigma^2)$, with
  $\log k_{a,i} = \beta_a + u_{a,i}$ and $\log Cl_i = \beta_{Cl} + u_{Cl,i}$
  (Pinheiro & Bates, 2000, ch. 8). Priors $N(0, 10^2)$ on the fixed effects
  and Exp(1) on the three SDs.
- **Route:** RTMB's Laplace over the 24 subject effects (`random = c("u_ka",
  "u_cl")`), `inla_vanilla()` on the 6 population parameters with `-obj$fn`
  and `-obj$gr`. Subject effects from `sampling(fun = ...)`: at each draw,
  rerun `obj$fn`, then draw from $N(\hat u, H^{-1})$ with the inner mode
  (`obj$env$last.par[obj$env$random]`) and Hessian
  (`obj$env$spHess(par, random = TRUE)`).
- **Data:** `Theoph` (base R), 12 subjects, 132 concentrations.
- **Result:** 0.4 s against 6 s of Stan sampling. Population means within
  0.05 SDs, SDs of the log SDs within 4%, but the fixed effects log ka and log
  Cl 13% too narrow (`summary()` flags log ka's ridge bend, not log Cl's). 22
  of the 24 subject effects within 3% in SD. Log ka of subject 9 (Stan
  skewness 1.8) is 27% too narrow: a Gaussian at the inner mode cannot follow
  it, which is what R-INLA's simplified Laplace is for. No noise from the
  inner optimisation's warm start.
- **Finding:** integrating the effects out beat treating them as parameters
  (all 30 in `inla_vanilla()`): there, every subject effect came out 12% too
  narrow and the fit warned of a 1.4 SD VB shift, because the joint mode of a
  hierarchical model has its SDs too small. The nested structure is better,
  not only cheaper.
- **Next:** a dataset where Stan slows and the Laplace route should not.
  Candidates: `nlme::Phenobarb` (59 newborns, multiple doses, sparse
  sampling, a classic NONMEM example, needs dose superposition), the warfarin
  PK/PD data in nlmixr2data (to check), and a simulated population of several
  hundred subjects from the Theoph fit, for scaling. Compare with nlmixr2's
  FOCEi and with Stan. The skewed subject effects are the known weak point.

## Candidates not yet tried

From searches of the R-INLA discussion group, GitHub (r-inla and inlabru
issues and discussions) and papers whose repositories ship data and a JAGS or
Stan reference. Stack Overflow, Cross Validated and Reddit could not be
reached by the search tools. Except for the frontier thread, the details were
relayed by search agents: read each source before citing it.

### From GitHub and papers

- **Nonlinear mixed model with a Box-Cox effort correction.** inlabru
  Discussion #252, <https://github.com/inlabru-org/inlabru/discussions/252>
  (November 2024). Negative binomial counts, $\log\mu = \text{site}_i +
  \text{strat}_s + \text{trend}_s\,\text{year} + B_s(\text{eff}^{p_s} -
  1)/p_s$, with $B_s$ and $p_s$ random per stratum. inlabru works after a
  reformulation but crashes on about 30% of simulated data sets, and the
  linearisation breaks near $p = 0$ (Lindgren). Simulation code in the thread
  (10 strata, 15 sites each, 50 years). Route: RTMB's Laplace over about 190
  effects, 6-8 outer parameters. Same pattern as the Theoph prototype.
- **Hierarchical stacking.** <https://github.com/hrue/r-inla/discussions/112>
  and <https://github.com/inlabru-org/inlabru/discussions/281> (June-August
  2025). The likelihood of Yao, Pirš, Vehtari and Gelman (2022),
  $\sum_i \log\sum_k w_k(x_i)\, p_{k,-i}$, with softmax weights depending on
  covariates. Rue: not supported. Even INLA's new experimental `cloglike`
  interface could not express it (the poster's conclusion, August 2025).
  Public replica: <https://github.com/yao-yl/hierarchical-stacking-code>
  (`well_switch`, 3,020 Bangladesh wells, with a Stan reference), after
  computing the leave-one-out densities of K logistic models. 15-40
  parameters, direct `logpost` with a log-sum-exp per row.
- **INLA within Monte Carlo examples.** Berild et al. (2021), arXiv
  2103.02721, <https://github.com/berild/inla-mc>, data and JAGS code in the
  repo. Three small cases R-INLA can fit only conditional on some parameters:
  a zero-inflated Poisson with covariates on the zero probability (`fish.csv`,
  250 rows, 5 parameters, `pscl::zeroinfl` as an ML reference), a
  two-component Poisson mixture (`MixtureInf::earthquake`, 107 yearly counts,
  3 parameters, sum the classes), and a Gaussian RW2 regression with log
  variance linear in $x$ (`SemiPar::lidar`, 221 rows, 4 outer parameters plus
  the RW2 field, by RTMB or a closed-form Gaussian marginal).
- **N-mixture model with observation-level detection covariates.**
  <https://github.com/Peprah94/INLA_within_nimble> (Adjei and O'Hara, 2024):
  R-INLA fits N-mixture models except when detection covariates vary by
  visit. $y_{ij} \sim \text{Bin}(N_i, p_{ij})$, $N_i \sim
  \text{Poisson}(\lambda_i)$, on `unmarked::mallard` (239 sites). 8
  parameters, sum $N_i$ from 0 to a bound inside `logpost`.
- **Published "INLA cannot" statements, no code found:** Wilson and Wakefield,
  arXiv 2011.07129 (a nonlinear mean function), and Osgood-Zimmerman and
  Wakefield, arXiv 2103.09929, section 7 (joint cancer incidence and mortality
  with $q = pr$, which "prohibits" INLA; fitted with TMB).

### From the R-INLA discussion group

- **Heckman sample selection.** <https://groups.google.com/g/r-inla-discussion-group/c/ewl8Wtwoeg0>
  (14 May 2020, follow-ups 28 May 2022) and
  <https://groups.google.com/g/r-inla-discussion-group/c/7yc4Ff_9KuI> and
  <https://groups.google.com/g/r-inla-discussion-group/c/pjYYwsKVB_s> (May
  2022). Probit selection plus a Gaussian outcome with correlated errors. Rue:
  only a two-step procedure is possible. Attempts with `iidkd` and `copy` gave
  a correlation near 0.01 against 0.14-0.2 from ML. Closed form: not selected
  $\log\Phi(-z^\top\gamma)$, selected $\log\phi(r)/\sigma +
  \log\Phi((z^\top\gamma + \rho r)/\sqrt{1 - \rho^2})$ with $r = (y -
  x^\top\beta)/\sigma$. About 13 parameters. Public data,
  `sampleSelection::Mroz87` (753 women, 428 working), with an ML reference in
  `sampleSelection::selection()`. The skewed posterior of $\rho$ is the
  interesting output.
- **Beta regression with random dispersion.**
  <https://groups.google.com/g/r-inla-discussion-group/c/nVBanSVf2Go> (20 May
  2021, Murray Logan). Coral cover, $\text{Beta}(\mu_i\phi_t, (1 -
  \mu_i)\phi_t)$ with year-specific $\log\phi_t$ drawn from a common prior,
  plus reef effects. Lindgren: a hierarchical layer on the variances does not
  fit R-INLA's structure. His brms fit took about 11 hours. Route: RTMB's
  Laplace over about 165 random effects (reefs, years, $\log\phi_t$), about 5
  outer parameters. No data, simulate. The best showcase of the "random
  effects on a dispersion" blocker.
- **Shared scaled linear predictor across two likelihoods.**
  <https://groups.google.com/g/r-inla-discussion-group/c/y2gVEisEPb4> (May-June
  2025, Scott Foster). $y_1 \sim \text{Poisson}(e^{\alpha_1 + X\beta})$, $y_2
  \sim N(\alpha_2 + \gamma X\beta, \sigma^2)$, with a shared spatial effect in
  the real application. R-INLA can do it only with a fake-data and `copy`
  hack, and the poster reported odd estimates. Simulation code is in the
  thread. About 7 parameters (plus a field for RTMB).
- **Shared latent occupancy state.**
  <https://groups.google.com/g/r-inla-discussion-group/c/HSWmsCWK-iI> (March
  2019). $z_j \sim \text{Bernoulli}(\psi_j)$, $Y_j \sim \text{Bin}(5, z_j p)$,
  $W_j \sim \text{Poisson}(E_j z_j e^{\alpha})$. Rue: cannot be done as
  written. Sum $z_j$ out inside `logpost`. Simulate.
- **Misclassified binary outcome.**
  <https://groups.google.com/g/r-inla-discussion-group/c/OMeZa08jrRg> (5 June
  2023). $P(y = 1) = s\,\text{expit}(\eta) + (1 - c)(1 - \text{expit}(\eta))$
  with sensitivity $s$ and specificity $c$ as parameters. Lindgren: cannot be
  done in plain INLA. Closed form, informative priors on $s$ and $c$ needed.
- **Surveillance random walk scaled by sensitivities.**
  <https://groups.google.com/g/r-inla-discussion-group/c/u12aN7SyYos> (July
  2021). $i_t = i_{t-1} + o_t$, $Y_t \sim \text{Bin}(n_t,
  \text{expit}(\beta_Y i_t))$, $Z_t \sim \text{Bin}(m_t, \text{expit}(\beta_Z
  i_t))$. Unresolved `copy` attempt. RTMB's Laplace over the states, 3 outer
  parameters (fix one scale for identifiability).
- **Thinned Bernoulli.**
  <https://groups.google.com/g/r-inla-discussion-group/c/5Xl7ipYXZck> (31
  August 2021). $y \sim \text{Bernoulli}(q\,\text{expit}(\eta))$ with known
  $q$, plus a spatial field. Rue: the scaled link is not implemented. A
  non-spatial version is one line of `logpost`.

Too big or too vague for now: censored data needing a three-parameter Weibull
or a mixture (July 2026,
<https://groups.google.com/g/r-inla-discussion-group/c/AaJtm1gTcF0>), ODE
compartmental models (2025,
<https://groups.google.com/g/r-inla-discussion-group/c/BUJAz7f_G0Q> and
<https://groups.google.com/g/r-inla-discussion-group/c/6iEoDuYhhsM>, too many
parameters and need an AD-capable ODE solver), and a quadratic in a latent
spatial field (April-May 2022,
<https://groups.google.com/g/r-inla-discussion-group/c/Zi3S-mbORPU>).

## Notes on the RTMB route

- With random effects the sign matters: TMB minimises over them, so the
  objective is the negative log joint, and `inla_vanilla()` gets
  `function(theta) -obj$fn(theta)` and `-as.vector(obj$gr(theta))`. Without
  random effects the function can return the log-posterior and `obj$fn`,
  `obj$gr` go in as they are.
- Put the priors on the outer parameters, with their Jacobians, in the
  objective.
- Call RTMB's functions with `RTMB::` rather than attaching it: `library(RTMB)`
  masks `dnorm()`, `solve()`, `diag()` and others, and slows plain R code
  several times.
- An AD gradient costs about two objective calls, so it buys accuracy, not
  speed (`vignettes/articles/rtmb.qmd`).
