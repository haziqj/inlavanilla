# print() shows the estimates and the leverage

    Code
      print(res)
    Output
      -- Leave-one-out cross-validation ---------------------------------- 25 units --
      
                Estimate   SE
      elpd_loo     -35.6  4.6
      p_loo          4.4  1.4
      looic         71.1  9.3
      
      v Second-order expansion for every unit. Largest leverage k = 0.47, for unit 7.

---

    Code
      print(waic(fit_lin, ll_lin, X = lin_X, y = lin_y))
    Output
      -- WAIC ------------------------------------------------------------ 25 units --
      
                 Estimate   SE
      elpd_waic     -35.4  4.6
      p_waic          4.3  1.3
      waic           70.8  9.2
      
      v Second-order expansion for every unit.

# print() names the diagonal shortcut

    Code
      print(loo(fit_big, ll_lin, X = big_X, y = big_y, hessian = "diagonal"))
    Output
      -- Leave-one-out cross-validation --------------------------------- 400 units --
      
                Estimate    SE
      elpd_loo    -430.6  13.5
      p_loo          3.1   0.3
      looic        861.3  27.1
      
      v First-order expansion with the trace of the second-order term. Largest summed
        leverage k = 0.029, for unit 369.

---

    Code
      print(waic(fit_big, ll_lin, X = big_X, y = big_y, hessian = "diagonal"))
    Output
      -- WAIC ----------------------------------------------------------- 400 units --
      
                 Estimate    SE
      elpd_waic    -430.6  13.5
      p_waic          3.1   0.3
      waic          861.3  27.1
      
      v First-order expansion with the trace of the second-order term.

