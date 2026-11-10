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
      -- INLAvanilla <version> -- Laplace + one-step VB --
      
      Log marginal likelihood: -16.25 | Converged in <n> iteration(s) (<time>)
      
                    Mean     SD   2.5%    50%  97.5%   Mode
      1 intercept  0.581  0.276  0.039  0.581  1.123  0.581
      2 slope      0.597  0.226  0.154  0.597  1.039  0.597

---

    Code
      print(summary(fit_g_vb))
    Output
      -- INLAvanilla <version> summary -- Laplace + one-step VB --
      
      Log marginal likelihood: -16.246 | Converged in <n> iteration(s) (<time>)
      Laplace -16.277 + VB term 0.031
      <n> logpost calls: mode <time>, VB <time>, marginals <time>
      
      -- Posterior marginals --
      
                    Mean     SD   2.5%    50%  97.5%   Mode
      1 intercept  0.581  0.276  0.039  0.581  1.123  0.581
      2 slope      0.597  0.226  0.154  0.597  1.039  0.597
      
      -- Diagnostics --
      
                   VB shift
      1 intercept    -0.24 
      2 slope         0.092
      
      v Mode found in <n> iteration(s) (Newton step <g> SD).
      v Hessian condition number 3.6, in SD units.
      v VB correction in one step: max shift 0.24 SD, gain ratio 0.977.
      i Gaussian marginals: no scans, so no skew-normal checks.

---

    Code
      print(summary(fit_g))
    Output
      -- INLAvanilla <version> summary -- Laplace --
      
      Log marginal likelihood: -16.277 | Converged in <n> iteration(s) (<time>)
      <n> logpost calls: mode <time>, VB <time>, marginals <time>
      
      -- Posterior marginals --
      
                    Mean     SD   2.5%    50%  97.5%   Mode
      1 intercept  0.648  0.276  0.106  0.648  1.190  0.648
      2 slope      0.576  0.226  0.133  0.576  1.019  0.576
      
      -- Diagnostics --
      
      v Mode found in <n> iteration(s) (Newton step <g> SD).
      v Hessian condition number 3.6, in SD units.
      i No VB correction (vb = "off").
      i Gaussian marginals: no scans, so no skew-normal checks.

