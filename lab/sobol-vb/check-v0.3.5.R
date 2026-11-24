################################################################################
#
# Check that sobol_vb() reproduces the VB shift of INLAvaan 0.3.5
#
# Fits three models with INLAvaan 0.3.5, takes from each fit its log-posterior,
# gradient, mode, covariance, factor L and saturated means, and compares the
# shift of sobol_vb() on those inputs with the one the fit found. Then it
# reports how far the shift moves with inlavanilla's factor cov_root() in
# place of INLAvaan's Cholesky factor. Run from the root of inlavanilla, with
# the source of INLAvaan 0.3.5 in the directory INLAVAAN_035, for example
#   git -C ~/GitHub/INLAvaan archive v0.3.5 | tar -x -C /tmp/inlavaan-0.3.5
#
################################################################################

## ----- Configuration ---------------------------------------------------------

src_035 <- Sys.getenv("INLAVAAN_035", "/tmp/inlavaan-0.3.5")
sobol_dir <- "lab/sobol-vb"
suppressMessages(pkgload::load_all(src_035, quiet = TRUE))
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = FALSE))
source(file.path(sobol_dir, "sobol-vb.R"))

## ----- What each fit holds at its end ----------------------------------------

# The closures and the VB inputs and result, from the frame of inlavaan()
keep <- c(
  "joint_lp",
  "joint_lp_grad",
  "theta_star",
  "Sigma_theta",
  "L",
  "fp_idx",
  "vb_shift",
  "vb_mcse"
)
invisible(suppressMessages(trace(
  "inlavaan",
  exit = substitute(
    assign(".fit_035", mget(KEEP), envir = globalenv()),
    list(KEEP = keep)
  ),
  where = asNamespace("INLAvaan"),
  print = FALSE
)))

hs <- "visual =~ x1 + x2 + x3
       textual =~ x4 + x5 + x6
       speed =~ x7 + x8 + x9"
poldem <- "ind60 =~ x1 + x2 + x3
           dem60 =~ y1 + y2 + y3 + y4
           dem65 =~ y5 + y6 + y7 + y8
           dem60 ~ ind60
           dem65 ~ ind60 + dem60
           y1 ~~ y5
           y2 ~~ y4 + y6
           y3 ~~ y7
           y4 ~~ y8
           y6 ~~ y8"
fits <- list(
  hs_cfa = function() {
    acfa(hs, lavaan::HolzingerSwineford1939, verbose = FALSE)
  },
  hs_cfa_means = function() {
    acfa(
      hs,
      lavaan::HolzingerSwineford1939,
      meanstructure = TRUE,
      verbose = FALSE
    )
  },
  poldem = function() {
    asem(poldem, lavaan::PoliticalDemocracy, verbose = FALSE)
  }
)

## ----- Compare ---------------------------------------------------------------

rows <- lapply(names(fits), function(nm) {
  set.seed(1)
  invisible(suppressWarnings(fits[[nm]]()))
  f <- .fit_035
  sds <- sqrt(diag(f$Sigma_theta))
  same_l <- sobol_vb(
    f$joint_lp,
    f$joint_lp_grad,
    f$theta_star,
    f$Sigma_theta,
    L = f$L,
    fixed = f$fp_idx
  )
  root <- sobol_vb(
    f$joint_lp,
    f$joint_lp_grad,
    f$theta_star,
    f$Sigma_theta,
    L = inlavanilla:::cov_root(f$Sigma_theta),
    fixed = f$fp_idx
  )
  data.frame(
    model = nm,
    m = length(f$theta_star),
    saturated_means = length(f$fp_idx),
    max_shift_sd = max(abs(f$vb_shift) / sds),
    same_L_diff_sd = max(abs(same_l$shift - f$vb_shift) / sds),
    same_L_mcse_diff_sd = max(abs(same_l$mcse - f$vb_mcse) / sds),
    cov_root_diff_sd = max(abs(root$shift - f$vb_shift) / sds)
  )
})
print(do.call(rbind, rows), digits = 3)

# Result, 2026-11-24 (R 4.6.1, lavaan 0.7.3, macOS): with INLAvaan's own
# factor the shift and its standard error are reproduced exactly, with a
# difference of 0 on all three models. With cov_root() the Sobol points are
# rotated differently, and the shift moves by the quadrature error of 64
# points: 0.068, 0.179 and 0.078 posterior SDs at most.
