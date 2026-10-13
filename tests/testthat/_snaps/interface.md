# verbose and silent cannot both be TRUE

    Code
      inla_vanilla(lp_gauss, c(a = 0.5), verbose = TRUE, silent = TRUE)
    Condition
      Error:
      ! `verbose` and `silent` cannot both be `TRUE`.
      i `verbose = TRUE` prints a full diagnostic of the fit, and `silent = TRUE` prints nothing.

# validation errors are informative

    Code
      inla_vanilla("lp", c(a = 0))
    Condition
      Error:
      ! `logpost` must be a function, not a string.

---

    Code
      inla_vanilla(function(theta) theta, c(a = 0.5, b = 0.5))
    Condition
      Error:
      ! `logpost` must return a single number.
      x `logpost(init)` returned a double vector of length 2.

---

    Code
      inla_vanilla(function(theta) -Inf, c(a = 0.5))
    Condition
      Error:
      ! `logpost` is not finite at `init`.
      x `logpost(init)` returned -Inf.
      i Start inside the support of the posterior, and check that every parameter is on an unconstrained scale.

---

    Code
      inla_vanilla(lp_gauss, c(a = 0.5), control = list(foo = 1, ngrid = 21))
    Condition
      Error:
      ! Unknown `control` entries: "foo" and "ngrid".
      i Valid entries are "n_grid", "log_threshold", "vb_maxit", "vb_tol", and "nlminb".

# a wrong analytic gradient triggers the mismatch warning

    Code
      invisible(suppressMessages(inla_vanilla(lp_gauss, c(a = 0.3, b = 0.2), grad = gr_wrong)))
    Condition
      Warning:
      `grad` does not match finite differences of `logpost` at `init`.
      x For b: `grad` gives 1.6, finite differences give 0.8.
      i A wrong gradient gives a wrong mode and wrong marginals.

# an improper posterior with no mode fails with the Hessian error

    Code
      inla_vanilla(function(theta) 0, c(a = 0.3), silent = TRUE)
    Condition
      Error:
      ! The Hessian at the posterior mode is not positive definite.
      i The optimiser may have stopped short of the mode, or the posterior may be improper or not identified.
      i Optimiser message: <msg>

# verbose output keeps its layout

    Code
      invisible(inla_vanilla(lp_reg, init_reg, verbose = TRUE))
    Message
      i 2 parameters: intercept and slope
        Derivatives by finite differences of logpost
        logpost at init #
      v Found the posterior mode in <n> iteration(s) [<time>]
        nlminb: <message>, after <n> iteration(s) (<n> function and <n> gradient evaluations)
        logpost at the mode #, whitened gradient # SD (warning above 0.01 SD)
        Hessian by second differences of logpost: 6 calls, and <n> more to set the steps
        Posterior correlation matrix: largest correlation # (intercept, slope), smallest eigenvalue #
        Laplace log marginal likelihood #
        <n> logpost calls
      v VB mean correction, max shift # SD [<time>]
        Fixed-point steps, largest over the parameters (SD): <steps>
        Converged in <n> iteration(s), <anderson>; no fallback
        Largest shift # SD, for intercept
        KL term #, giving a log marginal likelihood of #
        <n> logpost calls
      v Fitted 2 skew-normal marginals [<time>]
                   Shape NMAD Bend Tails gamma1
        intercept # # # ok #
        slope # # #! ok #
        ! marks what summary() flags: NMAD above 0.05, bend above 0.1, heavy tails
        <n> logpost calls
      i Done: <n> logpost calls in <time>

---

    Code
      invisible(inla_vanilla(lp_reg, init_reg, grad = gr_reg, verbose = TRUE))
    Message
      i 2 parameters: intercept and slope
        Gradient from grad, checked at init: relative discrepancy # (warning above 0.001), <n> logpost calls
        logpost at init #
      v Found the posterior mode in <n> iteration(s) [<time>]
        nlminb: <message>, after <n> iteration(s) (<n> function and <n> gradient evaluations)
        logpost at the mode #, whitened gradient # SD (warning above 0.01 SD)
        Hessian by central differences of grad: 4 calls, and <n> logpost calls to set the steps
        Jacobian asymmetry # (warning above 0.001)
        Posterior correlation matrix: largest correlation # (intercept, slope), smallest eigenvalue #
        Laplace log marginal likelihood #
        <n> logpost and <n> grad calls
      v VB mean correction, max shift # SD [<time>]
        Fixed-point steps, largest over the parameters (SD): <steps>
        Converged in <n> iteration(s), <anderson>; no fallback
        Largest shift # SD, for intercept
        KL term #, giving a log marginal likelihood of #
        <n> grad calls
      v Fitted 2 skew-normal marginals [<time>]
                   Shape NMAD Bend Tails gamma1
        intercept # # # ok #
        slope # # #! ok #
        ! marks what summary() flags: NMAD above 0.05, bend above 0.1, heavy tails
        <n> logpost and <n> grad calls
      i Done: <n> logpost and <n> grad calls in <time>

# print() shows the summary table and returns the fit invisibly

    Code
      print(fit_reg)
    Output
      -- INLAvanilla --------------------------------------------- Skew-normal + VB --
      
                  Mean     SD    2.5%    50%  97.5%   Mode
      intercept  0.582  0.286  -0.025  0.598  1.098  0.632
      slope      0.594  0.222   0.174  0.588  1.047  0.576
      
      Log marginal likelihood: -16.25
      ! 1 check flagged by `summary()`.

# print(summary()) shows marginals, diagnostics and the evidence

    Code
      print(summary(fit_reg))
    Output
      -- INLAvanilla summary ------------------------------------- Skew-normal + VB --
      
      Posterior marginals
                  Mean     SD    2.5%    50%  97.5%   Mode
      intercept  0.582  0.286  -0.025  0.598  1.098  0.632
      slope      0.594  0.222   0.174  0.588  1.047  0.576
      
      Diagnostics
                 VB shift  SN - VB  Shape   NMAD   Bend
      intercept     -0.24  -0.0027   -1.6   0.01  0.037
      slope         0.087  -0.0069    1.1  0.015   0.14
      
      v Mode found in <n> iteration(s) (gradient <g> SDs).
      v VB correction in <n> iteration(s): max shift 0.24 SDs.
      ! The straight scan misses up to 0.14 nats at +-2 SDs for slope: the marginal
        may be too narrow.
      v Skew-normal fits: max NMAD 0.015.
      v Skew-normal and VB means agree to within 0.0069 SDs.
      
      Log marginal likelihood: -16.248 (Laplace -16.277, VB term 0.029)
      <n> logpost calls in <time> (mode <time>, VB <time>, marginals <time>)

---

    Code
      print(summary(fit_reg_novb))
    Output
      -- INLAvanilla summary ------------------------------------------ Skew-normal --
      
      Posterior marginals
                  Mean     SD    2.5%    50%  97.5%   Mode
      intercept  0.582  0.286  -0.025  0.598  1.098  0.632
      slope      0.594  0.222   0.174  0.588  1.047  0.576
      
      Diagnostics
                 Shape   NMAD   Bend
      intercept   -1.6   0.01  0.037
      slope        1.1  0.015   0.14
      
      v Mode found in <n> iteration(s) (gradient <g> SDs).
      i No VB correction (vb = FALSE).
      ! The straight scan misses up to 0.14 nats at +-2 SDs for slope: the marginal
        may be too narrow.
      v Skew-normal fits: max NMAD 0.015.
      
      Log marginal likelihood: -16.277
      <n> logpost calls in <time> (mode <time>, VB <time>, marginals <time>)

# summary() rejects probabilities outside (0, 1)

    Code
      summary(fit_reg, probs = c(0.5, 1))
    Condition
      Error:
      ! `probs` must be probabilities strictly between 0 and 1.

# plot() rejects parameters that are not in the fit

    Code
      plot(fit_reg, pars = c("sigma", "tau"))
    Condition
      Error:
      ! Unknown parameters: "sigma" and "tau".

