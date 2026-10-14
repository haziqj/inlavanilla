# Run from the package root to install the approved logo and export pkgdown
# browser icons locally. Requires magick and openssl, only for this script.

source_logo <- "lab/logo/inlavanilla-photoreal-v4.png"
package_logo <- "man/figures/logo.png"
favicon_dir <- "pkgdown/favicon"

stopifnot(file.exists(source_logo))
dir.create(dirname(package_logo), recursive = TRUE, showWarnings = FALSE)
dir.create(favicon_dir, recursive = TRUE, showWarnings = FALSE)

# Ignore faint stray pixels outside the hex when trimming the generated canvas.
# Keep the high-resolution artwork for favicon exports; only the package copy
# uses the standard usethis 240 x 278 canvas (displayed at 120 x 139).
logo <- magick::image_trim(magick::image_read(source_logo), fuzz = 10)
package_image <- magick::image_extent(
  magick::image_resize(logo, "240x278"),
  "240x278",
  gravity = "center",
  color = "none"
)
magick::image_write(package_image, package_logo, format = "png")

# Browser icons need square canvases, with the hex fitted without stretching.
square_icon <- function(size) {
  geometry <- paste0(size, "x", size)
  magick::image_extent(
    magick::image_resize(logo, geometry),
    geometry,
    gravity = "center",
    color = "none"
  )
}

sizes <- c(
  "favicon-96x96.png" = 96,
  "apple-touch-icon.png" = 180,
  "web-app-manifest-192x192.png" = 192,
  "web-app-manifest-512x512.png" = 512
)
for (filename in names(sizes)) {
  icon <- square_icon(sizes[[filename]])
  magick::image_write(icon, file.path(favicon_dir, filename), format = "png")
}

ico <- do.call(
  c,
  lapply(c(16, 32, 48), square_icon)
)
magick::image_write(ico, file.path(favicon_dir, "favicon.ico"), format = "ico")

# pkgdown also links favicon.svg. Embed the PNG to preserve the photographic
# artwork; this is an SVG container, not a vector tracing.
png <- magick::image_write(
  square_icon(512),
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
