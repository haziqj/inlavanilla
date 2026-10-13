################################################################################
#
# Build the original vector logo concept, retained for reference.
# The approved package logo is inlavanilla-photoreal-v4.png, copied to
# man/figures/logo.png for the README and pkgdown.
#
# A skew-normal posterior sits on a vanilla pod, which serves as its x-axis,
# with the Laplace Gaussian dashed behind it and a vanilla orchid in bloom at
# the mode. The wordmark follows the INLAvaan sticker: the R-INLA letterforms
# for INLA, then Fira Code Bold. Every shape, the text included, is a path, so
# the SVG renders the same without the font installed.
#
################################################################################

## ----- Configuration ---------------------------------------------------------

library(magick)

# Run from the package root
svg_file <- "lab/logo/logo.svg"
png_file <- "lab/logo/logo-vector.png"

font_file <- systemfonts::match_fonts("Fira Code", weight = "bold")$path
if (!grepl("FiraCode", font_file)) {
  stop(
    "Install Fira Code (https://github.com/tonsky/FiraCode) to build the logo"
  )
}

hex_r <- 300
hex_w <- 2 * hex_r * cos(pi / 6)
border <- 8

# fmt: skip
pal <- list(
  bg_in     = "#FFFCF3",
  bg_out    = "#F6ECD2",
  border    = "#3A2415",
  teal      = "#479E98",  # INLAvaan's INLA teal
  lime      = "#BFCE52",  # INLAvaan's INLA lime
  laplace   = "#0F766E",  # plot.inla_vanilla() Laplace colour
  sn        = "#C2410C",  # plot.inla_vanilla() skew-normal colour
  bean_hi   = "#7A4B2A",
  bean      = "#4E2E19",
  bean_lo   = "#2A170C",
  leaf_hi   = "#A3C463",
  leaf_lo   = "#5E8A32",
  leaf_line = "#46702A",
  tepal_in  = "#C3C964",
  tepal_out = "#F3ECB0",
  tepal_ln  = "#8F8B3A",
  lip_out   = "#FCF4D2",
  lip_mid   = "#F8E39A",
  lip_in    = "#F2C23F",
  throat    = "#D98A1F",
  word      = "#4E2E19"
)
fan_cols <- c("#E3C063", "#E9A26C", "#86BDB7", "#C3CF77")

## ----- Helpers ---------------------------------------------------------------

f1 <- function(x) sub("\\.?0+$", "", sprintf("%.2f", x))

pts <- function(x, y) paste(f1(x), f1(y), sep = ",", collapse = " ")

path_line <- function(x, y, close = FALSE) {
  paste0("M", pts(x[1], y[1]), " L", pts(x[-1], y[-1]), if (close) " Z")
}

hex_points <- function(r) {
  a <- pi / 2 + (0:5) * pi / 3
  pts(hex_w / 2 + r * cos(a), hex_r - r * sin(a))
}

dsn <- function(t, xi, omega, alpha) {
  z <- (t - xi) / omega
  2 / omega * dnorm(z) * pnorm(alpha * z)
}

## =============================================================================
## Geometry
## =============================================================================

## ----- Vanilla pod as the x-axis ---------------------------------------------

pod_x0 <- 48
pod_x1 <- 468
pod_yc <- 391
pod_droop <- 7
pod_thick <- 13

pod_s <- function(x) (x - pod_x0) / (pod_x1 - pod_x0)
pod_mid <- function(x) pod_yc + pod_droop * (2 * pod_s(x) - 1)^2
pod_half <- function(x) {
  s <- pod_s(x)
  pod_thick * sin(pi * s)^0.55 * (1 - 0.25 * s)
}
pod_top <- function(x) pod_mid(x) - pod_half(x)

## ----- Posterior densities ---------------------------------------------------

sn_par <- c(xi = -1, omega = 1.55, alpha = 3.6)
sn_d <- function(t) dsn(t, sn_par[["xi"]], sn_par[["omega"]], sn_par[["alpha"]])
t_mode <- optimize(sn_d, c(-2, 2), maximum = TRUE)$maximum
d_max <- sn_d(t_mode)

mode_x <- hex_w / 2 - 8
px_per_unit <- 54
peak_px <- 200
to_x <- function(t) mode_x + (t - t_mode) * px_per_unit
tt <- t_mode +
  seq(pod_x0 + 12 - mode_x, pod_x1 - 16 - mode_x, length.out = 301) /
    px_per_unit
xx <- to_x(tt)
base_y <- pod_top(xx) + 2
to_y <- function(d) base_y - d / d_max * peak_px

# Laplace Gaussian: same mode and curvature as the skew-normal
h <- 1e-3
curv <- -diff(log(sn_d(t_mode + c(-h, 0, h))), differences = 2) / h^2
y_sn <- to_y(sn_d(tt))
y_lap <- to_y(d_max * exp(-0.5 * curv * (tt - t_mode)^2))
mode_y <- min(y_sn)

## ----- Faint fan of posteriors behind (as on the INLAvaan sticker) -----------

# A bloom of skew-normals whose modes sit near the main one, at assorted
# heights below the main peak. Keep only those whose tails come down onto the
# pod.
set.seed(11)
n_fan <- 16
fan_d <- list()
while (length(fan_d) < n_fan) {
  alpha <- runif(1, -1, 5)
  omega <- runif(1, 1, 2.2)
  m <- optimize(
    function(t) dsn(t, 0, omega, alpha),
    c(-3, 3),
    maximum = TRUE
  )$maximum
  xi <- t_mode - m + rnorm(1, 0, 0.45)
  d <- dsn(tt, xi, omega, alpha)
  if (max(d[1], d[length(d)]) < 0.01 * max(d)) {
    fan_d[[length(fan_d) + 1]] <- d / max(d) * d_max * runif(1, 0.4, 0.92)
  }
}
fan <- vapply(
  seq_len(n_fan),
  function(i) {
    sprintf(
      '<path d="%s" fill="none" stroke="%s" stroke-width="1.5" stroke-opacity="0.55"/>',
      path_line(xx, to_y(fan_d[[i]])),
      fan_cols[(i - 1) %% length(fan_cols) + 1]
    )
  },
  character(1)
)

## =============================================================================
## SVG pieces
## =============================================================================

## ----- Pod -------------------------------------------------------------------

px <- seq(pod_x0, pod_x1, length.out = 160)
pod_d <- paste(
  path_line(
    c(px, rev(px)),
    c(pod_mid(px) - pod_half(px), rev(pod_mid(px) + pod_half(px)))
  ),
  "Z"
)
wrinkle <- function(off, from, to, op) {
  w <- seq(from, to, length.out = 60)
  sprintf(
    '<path d="%s" fill="none" stroke="%s" stroke-width="1.2" stroke-linecap="round" stroke-opacity="%s"/>',
    path_line(w, pod_mid(w) + off * pod_half(w) + 0.6 * sin(w / 9)),
    pal$bean_hi,
    op
  )
}
pod_svg <- c(
  sprintf(
    '<path d="M%s Q%s %s" fill="none" stroke="%s" stroke-width="3.4" stroke-linecap="round"/>',
    pts(pod_x1 - 4, pod_mid(pod_x1)),
    pts(pod_x1 + 10, pod_mid(pod_x1) - 1),
    pts(pod_x1 + 13, pod_mid(pod_x1) - 10),
    pal$bean
  ),
  sprintf(
    '<path d="%s" fill="url(#pod)" stroke="%s" stroke-width="1.2"/>',
    pod_d,
    pal$bean_lo
  ),
  wrinkle(-0.5, pod_x0 + 28, pod_x1 - 36, 0.95),
  wrinkle(0.02, pod_x0 + 70, pod_x1 - 90, 0.5),
  wrinkle(0.48, pod_x0 + 120, pod_x1 - 56, 0.4)
)

## ----- Leaves and bud --------------------------------------------------------

leaf <- function(x, y, angle, len, wid) {
  body <- sprintf(
    "M0,0 C%s %s %s C%s %s 0,0 Z",
    pts(len * 0.22, -wid),
    pts(len * 0.7, -wid * 1.02),
    pts(len, -wid * 0.05),
    pts(len * 0.72, wid * 0.85),
    pts(len * 0.25, wid * 0.95)
  )
  rib <- sprintf(
    "M%s Q%s %s",
    pts(len * 0.04, 0),
    pts(len * 0.5, -wid * 0.12),
    pts(len * 0.93, -wid * 0.04)
  )
  sprintf(
    paste0(
      '<g transform="translate(%s) rotate(%s)">',
      '<path d="%s" fill="url(#leaf)" stroke="%s" stroke-width="1.3"/>',
      '<path d="%s" fill="none" stroke="%s" stroke-width="1.2" stroke-opacity="0.7"/>',
      "</g>"
    ),
    pts(x, y),
    f1(angle),
    body,
    pal$leaf_line,
    rib,
    pal$leaf_line
  )
}

bud <- function(x, y, angle, len) {
  sprintf(
    paste0(
      '<g transform="translate(%s) rotate(%s)">',
      '<path d="M0,0 Q%s %s" fill="none" stroke="%s" stroke-width="3" stroke-linecap="round"/>',
      '<path d="M%s C%s %s %s C%s %s %s Z" fill="url(#bud)" stroke="%s" stroke-width="1.2"/>',
      "</g>"
    ),
    pts(x, y),
    f1(angle),
    pts(len * 0.3, -6),
    pts(len * 0.5, -4),
    pal$leaf_lo,
    pts(len * 0.46, -4),
    pts(len * 0.62, -15),
    pts(len * 0.95, -12),
    pts(len, -4),
    pts(len * 0.95, 4),
    pts(len * 0.62, 6),
    pts(len * 0.46, -4),
    pal$leaf_line
  )
}

## ----- Vanilla orchid --------------------------------------------------------

tepal <- function(angle, len, wid, bend) {
  body <- sprintf(
    "M0,0 C%s %s %s C%s %s 0,0 Z",
    pts(len * 0.3, -wid),
    pts(len * 0.78, -wid * 0.75 + bend),
    pts(len, bend),
    pts(len * 0.78, wid * 0.6 + bend),
    pts(len * 0.3, wid * 0.95)
  )
  rib <- sprintf(
    "M%s Q%s %s",
    pts(len * 0.12, 0),
    pts(len * 0.55, bend * 0.4),
    pts(len * 0.9, bend * 0.95)
  )
  sprintf(
    paste0(
      '<g transform="rotate(%s)">',
      '<path d="%s" fill="url(#tepal)" stroke="%s" stroke-width="1.4" stroke-linejoin="round"/>',
      '<path d="%s" fill="none" stroke="%s" stroke-width="1" stroke-opacity="0.5"/>',
      "</g>"
    ),
    f1(angle),
    body,
    pal$tepal_ln,
    rib,
    pal$tepal_ln
  )
}

# The labellum is a trumpet rolled round the column, with a ruffled rim and a
# yellow throat.
lip_d <- paste(
  "M-15,-9 C-20,4 -31,20 -33,36",
  "Q-29,47 -20,43 Q-14,52 -6,46 Q0,54 7,46 Q15,52 21,43 Q30,47 33,36",
  "C31,20 20,4 15,-9 Q0,-16 -15,-9 Z"
)
veins <- vapply(
  c(-15, -8, 0, 8, 15),
  function(v) {
    sprintf(
      '<path d="M%s Q%s %s" fill="none" stroke="%s" stroke-width="1.1" stroke-opacity="0.55" stroke-linecap="round"/>',
      pts(v * 0.3, 16),
      pts(v * 0.85, 27),
      pts(v * 1.45, 40),
      pal$throat
    )
  },
  character(1)
)

flower_svg <- function(x, y, scale = 1) {
  c(
    sprintf(
      '<g transform="translate(%s) scale(%s)">',
      pts(x, y),
      f1(scale)
    ),
    tepal(-90, 86, 15, 0),
    tepal(-153, 78, 14, -7),
    tepal(-27, 78, 14, 7),
    tepal(146, 84, 15, 6),
    tepal(34, 84, 15, -6),
    sprintf(
      '<path d="%s" fill="url(#lip)" stroke="%s" stroke-width="1.4" stroke-linejoin="round"/>',
      lip_d,
      pal$tepal_ln
    ),
    '<ellipse cx="0" cy="9" rx="10" ry="13" fill="url(#throat)"/>',
    veins,
    sprintf(
      '<path d="M-6,-10 Q0,-15 6,-10 Q5,-1 0,1 Q-5,-1 -6,-10 Z" fill="%s" stroke="%s" stroke-width="0.9" stroke-opacity="0.6"/>',
      pal$lip_mid,
      pal$tepal_ln
    ),
    "</g>"
  )
}

## ----- Wordmark --------------------------------------------------------------

# Font units at size 1000, y up. INLA is drawn in the R-INLA style, with its
# own spacing, then "vanilla" follows in Fira Code Bold.
cap <- 710.77
adv <- 615.39
sw <- 142 # stem width
bw <- 118 # bar height
dw <- 172 # horizontal width of a diagonal stroke

rect <- function(x0, x1, y0, y1) {
  list(x = c(x0, x1, x1, x0), y = c(y0, y0, y1, y1))
}

# fmt: skip
inla_glyphs <- list(
  I = list(
    rect(30, 420, cap - bw, cap),
    rect(225 - sw / 2, 225 + sw / 2, 0, cap),
    rect(30, 420, 0, bw)
  ),
  N = list(
    rect(540, 540 + sw, 0, cap),
    list(x = c(540, 540 + dw, 1000, 1000 - dw), y = c(cap, cap, 0, 0)),
    rect(1000 - sw, 1000, cap - bw, cap)
  ),
  L = list(
    rect(1130, 1130 + sw, 0, cap),
    rect(1130, 1530, 0, bw)
  ),
  A = local({
    l <- 1620
    r <- 2270
    a0 <- 1915
    hc <- 0.45 * cap
    x_left <- function(y) l + (a0 - l) * y / cap
    list(
      list(x = c(a0, a0 + dw, r, r - dw), y = c(cap, cap, 0, 0)),
      list(x = c(l, l + dw, x_left(hc) + dw, x_left(hc)), y = c(0, 0, hc, hc)),
      rect(x_left(hc), 1975, hc - bw, hc)
    )
  })
)
inla_cols <- c(I = pal$teal, N = pal$lime, L = pal$teal, A = pal$teal)
vanilla_x <- 2330

shaped <- systemfonts::shape_string("vanilla", path = font_file, size = 1000)
outl <- systemfonts::glyph_outline(
  shaped$shape$index,
  path = font_file,
  size = 1000,
  tolerance = 0.5
)

word_size <- 52
word_scale <- word_size / 1000
word_w <- vanilla_x + 7 * adv - 30
word_x0 <- (hex_w - word_w * word_scale) / 2
word_y0 <- 473

to_page <- function(x, y) {
  list(x = word_x0 + x * word_scale, y = word_y0 - y * word_scale)
}

word_svg <- c(
  unlist(lapply(names(inla_glyphs), function(g) {
    vapply(
      inla_glyphs[[g]],
      function(p) {
        q <- to_page(p$x, p$y)
        sprintf(
          '<path d="%s" fill="%s"/>',
          path_line(q$x, q$y, TRUE),
          inla_cols[[g]]
        )
      },
      character(1)
    )
  })),
  sprintf(
    '<path d="%s" fill="%s"/>',
    paste(
      unlist(lapply(seq_len(7), function(k) {
        g <- outl[outl$glyph == k, ]
        vapply(
          split(g, g$contour),
          function(cc) {
            q <- to_page(vanilla_x + (k - 1) * adv + cc$x, cc$y)
            path_line(q$x, q$y, TRUE)
          },
          character(1)
        )
      })),
      collapse = " "
    ),
    pal$word
  )
)

## =============================================================================
## Assemble
## =============================================================================

grad <- function(id, stops, x2 = 0, y2 = 1) {
  paste0(
    sprintf(
      '<linearGradient id="%s" x1="0" y1="0" x2="%s" y2="%s">',
      id,
      x2,
      y2
    ),
    paste(stops, collapse = ""),
    "</linearGradient>"
  )
}
grad_stop <- function(offset, col, op = 1) {
  sprintf(
    '<stop offset="%s" stop-color="%s" stop-opacity="%s"/>',
    offset,
    col,
    op
  )
}

defs <- c(
  "<defs>",
  sprintf(
    '<clipPath id="hex"><polygon points="%s"/></clipPath>',
    hex_points(hex_r - border)
  ),
  sprintf(
    '<radialGradient id="bg" cx="50%%" cy="40%%" r="62%%">%s%s</radialGradient>',
    grad_stop(0, pal$bg_in),
    grad_stop(1, pal$bg_out)
  ),
  grad("snfill", c(grad_stop(0, pal$sn, 0.34), grad_stop(1, pal$sn, 0.05))),
  grad(
    "pod",
    c(
      grad_stop(0, pal$bean_hi),
      grad_stop(0.45, pal$bean),
      grad_stop(1, pal$bean_lo)
    )
  ),
  grad(
    "leaf",
    c(grad_stop(0, pal$leaf_lo), grad_stop(1, pal$leaf_hi)),
    x2 = 1,
    y2 = 0
  ),
  grad(
    "bud",
    c(grad_stop(0, pal$leaf_hi), grad_stop(1, "#E4E6A0")),
    x2 = 1,
    y2 = 0
  ),
  grad(
    "tepal",
    c(grad_stop(0, pal$tepal_in), grad_stop(0.6, pal$tepal_out)),
    x2 = 1,
    y2 = 0
  ),
  sprintf(
    '<radialGradient id="lip" cx="50%%" cy="20%%" r="80%%">%s%s%s</radialGradient>',
    grad_stop(0, pal$lip_in),
    grad_stop(0.5, pal$lip_mid),
    grad_stop(1, pal$lip_out)
  ),
  sprintf(
    '<radialGradient id="throat" cx="50%%" cy="35%%" r="60%%">%s%s</radialGradient>',
    grad_stop(0, pal$throat, 0.85),
    grad_stop(1, pal$lip_in, 0)
  ),
  "</defs>"
)

flower_x <- mode_x
flower_y <- mode_y - 40

svg <- c(
  sprintf(
    '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %s 600" width="%s" height="600">',
    f1(hex_w),
    f1(hex_w)
  ),
  "<title>inlavanilla</title>",
  defs,
  sprintf(
    '<polygon points="%s" fill="url(#bg)"/>',
    hex_points(hex_r - border / 2)
  ),
  '<g clip-path="url(#hex)">',
  fan,
  sprintf(
    '<path d="%s" fill="url(#snfill)"/>',
    paste(path_line(c(xx, rev(xx)), c(y_sn, rev(base_y) + 8)), "Z")
  ),
  sprintf(
    '<path d="%s" fill="none" stroke="%s" stroke-width="2.6" stroke-dasharray="8 6" stroke-linecap="round"/>',
    path_line(xx, y_lap),
    pal$laplace
  ),
  sprintf(
    '<path d="%s" fill="none" stroke="%s" stroke-width="4.2" stroke-linejoin="round" stroke-linecap="round"/>',
    path_line(xx, y_sn),
    pal$sn
  ),
  pod_svg,
  leaf(flower_x - 4, flower_y + 22, 164, 118, 22),
  leaf(flower_x + 8, flower_y + 18, -12, 104, 20),
  bud(flower_x + 14, flower_y - 8, -58, 92),
  flower_svg(flower_x, flower_y, 1.12),
  word_svg,
  "</g>",
  sprintf(
    '<polygon points="%s" fill="none" stroke="%s" stroke-width="%s" stroke-linejoin="round"/>',
    hex_points(hex_r - border / 2),
    pal$border,
    f1(border)
  ),
  "</svg>"
)

writeLines(svg, svg_file)

## ----- Render the PNG --------------------------------------------------------

# Rendered large and scaled down, at the 240 px width usethis::use_logo() uses
logo <- image_read(svg_file, density = 166)
logo <- image_resize(image_convert(logo, "png"), "240x278")
image_write(logo, png_file, depth = 8)
