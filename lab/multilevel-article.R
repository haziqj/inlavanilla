################################################################################
#
# Precomputes the brms reference for the article "Multilevel models"
# (vignettes/articles/multilevel.qmd), so that the site build needs neither
# brms nor Stan. Fits lme4's sleepstudy with correlated random intercepts and
# slopes, on all 18 subjects and on the first 6, with long chains, and saves
# each parameter's posterior summary and a kernel density estimate to
# vignettes/articles/multilevel-brms.rds. The model and priors must match the
# article's logpost. Takes about two minutes, most of it compiling Stan. Run
# from the package root:
#   Rscript lab/multilevel-article.R
#
################################################################################

## ----- Configuration ----------------------------------------------------------

library(brms)
out_file <- "vignettes/articles/multilevel-brms.rds"
probs <- c(0.025, 0.5, 0.975)
brms_pars <- c(
  b0 = "b_Intercept",
  b1 = "b_Days",
  tau0 = "sd_Subject__Intercept",
  tau1 = "sd_Subject__Days",
  rho = "cor_Subject__Intercept__Days",
  sigma = "sigma"
)

## ----- Fit and summarise ------------------------------------------------------

fit_case <- function(n_subjects) {
  dat <- lme4::sleepstudy
  keep <- levels(dat$Subject)[seq_len(n_subjects)]
  dat <- droplevels(dat[dat$Subject %in% keep, ])
  fit <- brm(
    Reaction ~ 0 + Intercept + Days + (Days | Subject),
    data = dat,
    prior = c(
      prior(normal(0, 500), class = b),
      prior(normal(0, 50), class = sd),
      prior(normal(0, 50), class = sigma),
      prior(lkj(1), class = cor)
    ),
    chains = 4,
    cores = 4,
    iter = 7500,
    warmup = 2500,
    seed = 1,
    refresh = 0,
    silent = 2,
    control = list(adapt_delta = 0.999, max_treedepth = 12)
  )
  draws <- as.matrix(fit, variable = unname(brms_pars))
  colnames(draws) <- names(brms_pars)
  np <- nuts_params(fit)
  list(
    summary = data.frame(
      mean = colMeans(draws),
      sd = apply(draws, 2, sd),
      t(apply(draws, 2, quantile, probs)),
      check.names = FALSE
    ),
    density = lapply(
      as.data.frame(draws),
      function(d) density(d, n = 256)[c("x", "y")]
    ),
    draws = nrow(draws),
    rhat = max(rhat(fit, pars = unname(brms_pars))),
    divergences = sum(np$Value[np$Parameter == "divergent__"])
  )
}

ref <- list(
  subjects_18 = fit_case(18),
  subjects_6 = fit_case(6),
  brms_version = as.character(packageVersion("brms"))
)
saveRDS(ref, out_file, compress = "xz")
str(ref, max.level = 2)
