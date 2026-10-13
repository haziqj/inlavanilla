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

