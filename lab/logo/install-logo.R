# Run from the package root to install the approved logo and export pkgdown
# browser icons locally. Requires magick and openssl, only for this script.

source_logo <- "lab/logo/inlavanilla-photoreal-v4.png"
package_logo <- "man/figures/logo.png"
favicon_dir <- "pkgdown/favicon"

stopifnot(file.exists(source_logo))
dir.create(dirname(package_logo), recursive = TRUE, showWarnings = FALSE)
stopifnot(file.copy(source_logo, package_logo, overwrite = TRUE))
dir.create(favicon_dir, recursive = TRUE, showWarnings = FALSE)

logo <- magick::image_read(package_logo)
sizes <- c(
  "favicon-96x96.png" = 96,
  "apple-touch-icon.png" = 180,
  "web-app-manifest-192x192.png" = 192,
  "web-app-manifest-512x512.png" = 512
)
for (filename in names(sizes)) {
  icon <- magick::image_resize(
    logo,
    paste0(sizes[[filename]], "x", sizes[[filename]])
  )
  magick::image_write(icon, file.path(favicon_dir, filename), format = "png")
}

ico <- do.call(
  c,
  lapply(c(16, 32, 48), function(size) {
    magick::image_resize(logo, paste0(size, "x", size))
  })
)
magick::image_write(ico, file.path(favicon_dir, "favicon.ico"), format = "ico")

# pkgdown also links favicon.svg. Embed the PNG to preserve the photographic
# artwork; this is an SVG container, not a vector tracing.
png <- magick::image_write(
  magick::image_resize(logo, "512x512"),
  format = "png"
)
writeLines(
  sprintf(
    '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 512 512"><image width="512" height="512" href="data:image/png;base64,%s"/></svg>',
    openssl::base64_encode(png)
  ),
  file.path(favicon_dir, "favicon.svg")
)

writeLines(
  c(
    '{',
    '  "name": "inlavanilla",',
    '  "short_name": "inlavanilla",',
    '  "icons": [',
    '    {"src": "web-app-manifest-192x192.png", "sizes": "192x192", "type": "image/png"},',
    '    {"src": "web-app-manifest-512x512.png", "sizes": "512x512", "type": "image/png"}',
    '  ],',
    '  "theme_color": "#3A2415",',
    '  "background_color": "#FFFCF3",',
    '  "display": "standalone"',
    '}'
  ),
  file.path(favicon_dir, "site.webmanifest")
)
