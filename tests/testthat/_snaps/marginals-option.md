# marginals takes one of two values, or an abbreviation

    Code
      inla_vanilla(lp_reg, init_reg, marginals = "laplace", silent = TRUE)
    Condition
      Error:
      ! `marginals` must be one of "skewnormal" or "gaussian".

# print() and summary() describe Gaussian marginals

    Code
      print(fit_g_vb)
    Output
      -- INLAvanilla ------------------------------------------------- Laplace + VB --
      
                  Mean     SD   2.5%    50%  97.5%   Mode
      intercept  0.584  0.276  0.042  0.584  1.125  0.584
      slope      0.595  0.226  0.152  0.595  1.038  0.595
      
      Log marginal likelihood: -16.25

---

    Code
      print(summary(fit_g_vb))
    Output
      -- INLAvanilla summary ----------------------------------------- Laplace + VB --
      
      Log marginal likelihood: -16.249 (Laplace -16.277, VB term 0.029)
      <n> logpost calls in <time> (mode <time>, VB <time>, marginals <time>)
      
      -- Posterior marginals --
      
                  Mean     SD   2.5%    50%  97.5%   Mode
      intercept  0.584  0.276  0.042  0.584  1.125  0.584
      slope      0.595  0.226  0.152  0.595  1.038  0.595
      
      -- Diagnostics --
      
                 VB shift
      intercept    -0.23 
      slope         0.084
      
      v Mode found in <n> iteration(s) (Newton step <g> SD).
      v Hessian condition number 3.6, in SD units.
      v VB correction in <n> iteration(s): max shift 0.23 SD.
      i Gaussian marginals: no scans, so no skew-normal checks.

---

    Code
      print(summary(fit_g))
    Output
      -- INLAvanilla summary ---------------------------------------------- Laplace --
      
      Log marginal likelihood: -16.277
      <n> logpost calls in <time> (mode <time>, VB <time>, marginals <time>)
      
      -- Posterior marginals --
      
                  Mean     SD   2.5%    50%  97.5%   Mode
      intercept  0.648  0.276  0.106  0.648  1.190  0.648
      slope      0.576  0.226  0.133  0.576  1.019  0.576
      
      -- Diagnostics --
      
      v Mode found in <n> iteration(s) (Newton step <g> SD).
      v Hessian condition number 3.6, in SD units.
      i No VB correction (vb = FALSE).
      i Gaussian marginals: no scans, so no skew-normal checks.

