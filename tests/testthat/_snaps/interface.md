# validation errors are informative

    Code
      inla_vanilla(function(theta) 0, c(a = 0, b = 1, a = 2))
    Condition
      Error in `check_init()`:
      ! `init` has duplicated names: "a".

---

    Code
      inla_vanilla("lp", c(a = 0))
    Condition
      Error in `check_function()`:
      ! `logpost` must be a function, not a string.

---

    Code
      inla_vanilla(function(theta) theta, c(a = 0.5, b = 0.5))
    Condition
      Error in `check_logpost_at_init()`:
      ! `logpost` must return a single number.
      x `logpost(init)` returned a double vector of length 2.

---

    Code
      inla_vanilla(function(theta) -Inf, c(a = 0.5))
    Condition
      Error in `check_logpost_at_init()`:
      ! `logpost` is not finite at `init`.
      x `logpost(init)` returned -Inf.
      i Start inside the support of the posterior, and check that every parameter is on an unconstrained scale.

---

    Code
      inla_vanilla(lp_gauss, c(a = 0.5, b = 0.5), grad = function(theta) 1)
    Condition
      Error in `check_grad_at_init()`:
      ! `grad` must return a numeric vector of length 2.
      x `grad(init)` returned a number of length 1.

---

    Code
      inla_vanilla(lp_gauss, c(a = 0.5), control = list(foo = 1, ngrid = 21))
    Condition
      Error in `vanilla_control()`:
      ! Unknown `control` entries: "foo" and "ngrid".
      i Valid entries are "n_grid", "log_threshold", "vb_maxit", "vb_tol", and "nlminb".

# a wrong analytic gradient triggers the mismatch warning

    Code
      invisible(inla_vanilla(lp_gauss, c(a = 0.3, b = 0.2), grad = gr_wrong, verbose = FALSE))
    Condition
      Warning:
      `grad` does not match finite differences of `logpost` at `init`.
      x For b: `grad` gives 1.6, finite differences give 0.8.
      i A wrong gradient gives a wrong mode and wrong marginals.

# an improper posterior with no mode fails with the Hessian error

    Code
      inla_vanilla(function(theta) 0, c(a = 0.3), verbose = FALSE)
    Condition
      Error in `laplace_fit()`:
      ! The Hessian at the posterior mode is not positive definite.
      i The optimiser may have stopped short of the mode, or the posterior may be improper or not identified.
      i Optimiser message: <msg>

# print() shows the summary table and returns the fit invisibly

    Code
      print(fit_two)
    Output
      -- INLAvanilla -- <time> --
      2 parameters, skew-normal marginals, VB-corrected
      
                  Mean     SD    2.5%     50%  97.5%    Mode
      log_rate   0.595  0.263   0.047   0.607  1.078   0.631
      logit_p   -0.995  0.622  -2.294  -0.967  0.149  -0.910
      
      Log marginal likelihood: -14.93

# print(summary()) shows marginals, diagnostics and the evidence

    Code
      print(summary(fit_two))
    Output
      -- INLAvanilla summary ---------------------------------------------------------
      
      Posterior marginals
                  Mean     SD    2.5%     25%     50%     75%  97.5%    Mode
      log_rate   0.595  0.263   0.047   0.426   0.607   0.776  1.078   0.631
      logit_p   -0.995  0.622  -2.294  -1.396  -0.967  -0.566  0.149  -0.910
      
      Diagnostics
                VB shift   SN - VB  Shape    NMAD
      log_rate     -0.13  -0.00077   -1.4  0.0066
      logit_p      -0.12    -0.015   -1.4  0.0094
      
      v Mode found in 8 iterations (gradient <g> SDs).
      v VB correction in 3 iterations: max shift 0.13 SDs.
      v Skew-normal fits: max NMAD 0.0094.
      v Skew-normal and VB means agree to within 0.015 SDs.
      
      Log marginal likelihood: -14.935 (Laplace -14.950, VB term 0.015)
      <n> logpost calls in <time> (mode <time>, VB <time>, marginals <time>)

---

    Code
      print(summary(fit_novb))
    Output
      -- INLAvanilla summary ---------------------------------------------------------
      
      Posterior marginals
                  Mean     SD    2.5%     25%     50%     75%  97.5%    Mode
      log_rate   0.595  0.263   0.047   0.426   0.607   0.776  1.078   0.631
      logit_p   -0.995  0.622  -2.294  -1.396  -0.967  -0.566  0.149  -0.910
      
      Diagnostics
                VB shift  SN - VB  Shape    NMAD
      log_rate         0    -0.13   -1.4  0.0066
      logit_p          0    -0.13   -1.4  0.0094
      
      v Mode found in 8 iterations (gradient <g> SDs).
      i No VB correction (vb = FALSE).
      v Skew-normal fits: max NMAD 0.0094.
      
      Log marginal likelihood: -14.950
      <n> logpost calls in <time> (mode <time>, VB <time>, marginals <time>)

# summary() rejects probabilities outside (0, 1)

    Code
      summary(fit_two, probs = c(0.5, 1))
    Condition
      Error in `summary()`:
      ! `probs` must be probabilities strictly between 0 and 1.

# plot() rejects parameters that are not in the fit

    Code
      plot(fit_two, pars = c("sigma", "tau"))
    Condition
      Error in `plot()`:
      ! Unknown parameters: "sigma" and "tau".

