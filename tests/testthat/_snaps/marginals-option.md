# marginals takes one of two values, in full

    Code
      inla_vanilla(lp_reg, init_reg, marginals = "skew", verbose = FALSE)
    Condition
      Error:
      ! `marginals` must be one of "skewnormal" or "gaussian".

# print() and summary() describe Gaussian marginals

    Code
      print(fit_g_vb)
    Output
      -- INLAvanilla -- <time> --
      2 parameters, Gaussian marginals, VB-corrected
      
                  Mean     SD   2.5%    50%  97.5%   Mode
      intercept  0.582  0.276  0.041  0.582  1.124  0.582
      slope      0.596  0.226  0.153  0.596  1.038  0.596
      
      Log marginal likelihood: -16.25

---

    Code
      print(summary(fit_g_vb))
    Output
      -- INLAvanilla summary ---------------------------------------------------------
      
      Posterior marginals
                  Mean     SD   2.5%    25%    50%    75%  97.5%   Mode
      intercept  0.582  0.276  0.041  0.396  0.582  0.769  1.124  0.582
      slope      0.596  0.226  0.153  0.443  0.596  0.748  1.038  0.596
      
      Diagnostics
                 VB shift
      intercept     -0.24
      slope         0.087
      
      v Mode found in <n> iteration(s) (gradient <g> SDs).
      v VB correction in <n> iteration(s): max shift 0.24 SDs.
      i Gaussian marginals: no scans, so no skew-normal checks.
      
      Log marginal likelihood: -16.248 (Laplace -16.277, VB term 0.029)
      <n> logpost calls in <time> (mode <time>, VB <time>, marginals <time>)

---

    Code
      print(summary(fit_g))
    Output
      -- INLAvanilla summary ---------------------------------------------------------
      
      Posterior marginals
                  Mean     SD   2.5%    25%    50%    75%  97.5%   Mode
      intercept  0.648  0.276  0.106  0.461  0.648  0.834  1.190  0.648
      slope      0.576  0.226  0.133  0.424  0.576  0.728  1.019  0.576
      
      Diagnostics
      v Mode found in <n> iteration(s) (gradient <g> SDs).
      i No VB correction (vb = FALSE).
      i Gaussian marginals: no scans, so no skew-normal checks.
      
      Log marginal likelihood: -16.277
      <n> logpost calls in <time> (mode <time>, VB <time>, marginals <time>)

