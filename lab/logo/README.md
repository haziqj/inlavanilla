# Package logo

The approved artwork is `inlavanilla-photoreal-v4.png`. Earlier versions and
their generation prompts are retained here for reference.

From the package root, run `Rscript lab/logo/install-logo.R` to copy the
approved artwork to `man/figures/logo.png` and export the browser icons to
`pkgdown/favicon/`. This script requires magick and openssl; neither is a
runtime dependency of inlavanilla.

The heading in `README.Rmd` (and its generated `README.md`) displays the
package logo. pkgdown discovers the same image and copies the saved favicons
when building the site. No extra `_pkgdown.yml` configuration is needed.
After changing the artwork, regenerate the icons and run
`pkgdown::init_site()` before rebuilding an existing local site.

The original vector experiment lives in `logo.R` and `logo.svg`.
`logo.R` now writes `logo-vector.png` here, so it cannot overwrite the
approved package logo.
