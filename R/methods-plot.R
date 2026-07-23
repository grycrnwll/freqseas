# methods-plot.R
# Base-graphics plot() methods for the three freqseas classes, plus
# conditionally-registered ggplot2 autoplot() methods (ggplot2 is a Suggests
# dependency; every autoplot method starts with
# `rlang::check_installed("ggplot2")` and is registered via the
# `@exportS3Method ggplot2::autoplot` roxygen tag rather than an `@importFrom
# ggplot2` -- see the CRAN-safe Suggests-registration note on each autoplot
# method below).
#
# Palette and panel chrome follow the session's plotting conventions
# (seasonality_vintages/scripts/phase_demo.R): seasonal-bin ordinates /
# adjusted series in blue (#2a78d6), the removed seasonal component in
# orange (#eb6834), donor ordinates / the original series in a neutral
# muted ink (#898781), a recessive grid (#e1e0d9) and axis ink (#c3c2b7) on
# an off-white surface (#fcfcfb).

# ---------------------------------------------------------------------------
# Shared palette + data-prep helpers
# ---------------------------------------------------------------------------

.fs_pal <- list(
  blue    = "#2a78d6",  # seasonal-bin ordinates / adjusted series
  orange  = "#eb6834",  # removed seasonal component / EVT critical level
  muted   = "#898781",  # donor ordinates / original series
  grid    = "#e1e0d9",
  axis    = "#c3c2b7",
  surface = "#fcfcfb",
  ink     = "#0b0b0b",
  ink2    = "#52514e"
)

# `seas_sa$adjusted`/`$seasonal` are input-class-preserving (see
# new_seas_sa()): plain numeric, a `ts`, or -- for a tbl_ts-input fit
# produced through the collection layer -- a full `tbl_ts` re-clothed by
# `fs_reclothe_sa()` (R/collection.R), i.e. a multi-column tsibble, not a
# bare vector. `as.numeric()` works directly on the first two but errors on
# the third ("'list' object cannot be coerced to type 'double'"); this
# extracts the plotting-relevant numeric vector from all three uniformly.
.fs_as_numeric_series <- function(v) {
  if (inherits(v, "tbl_ts")) {
    mv <- tsibble::measured_vars(v)
    return(as.numeric(v[[mv[1L]]]))
  }
  as.numeric(v)
}

# Styled empty base-graphics panel: off-white surface, recessive grid, muted
# axis ink. Caller adds points/lines/legend afterwards. Mirrors the
# panel_base() convention in seasonality_vintages/scripts/phase_demo.R.
.fs_panel <- function(xlim, ylim, xlab, ylab, main, log = "") {
  graphics::plot(NA, xlim = xlim, ylim = ylim, log = log, xlab = xlab,
                 ylab = ylab, main = main, axes = FALSE)
  usr <- graphics::par("usr")
  x0 <- if (grepl("x", log, fixed = TRUE)) 10^usr[1] else usr[1]
  x1 <- if (grepl("x", log, fixed = TRUE)) 10^usr[2] else usr[2]
  y0 <- if (grepl("y", log, fixed = TRUE)) 10^usr[3] else usr[3]
  y1 <- if (grepl("y", log, fixed = TRUE)) 10^usr[4] else usr[4]
  graphics::rect(x0, y0, x1, y1, col = .fs_pal$surface, border = NA)
  graphics::grid(col = .fs_pal$grid, lty = 1, lwd = 0.6)
  graphics::axis(1, col = .fs_pal$axis, col.ticks = .fs_pal$axis)
  graphics::axis(2, col = .fs_pal$axis, col.ticks = .fs_pal$axis)
  graphics::box(col = .fs_pal$axis)
}

# Extracts the plotting-relevant data from a `seas_test` for the standardized
# periodogram views: omega/pi and the standardized periodogram on the
# positive-frequency grid, a per-ordinate group label ("seasonal" / "donor" /
# "other"), the bin-boundary positions (omega/pi), and the EVT rejection
# threshold on the periodogram's own scale (`max(I[J0]) + evt$critical`,
# the exact level a seasonal-bin ordinate would need to exceed to trigger
# detection at `alpha`, since the test statistic is
# `max(I[J1]) - max(I[J0])`).
#
# Coordinates: `donor_pool` (from `seas_test()`) holds positions into the
# length-Jmax positive-frequency-ARRAY (`Jmax = (n_e - 1) %/% 2`), the same
# array `partition$grid$omega_pos`'s first `Jmax` entries index (see
# R/seas_test.R's `pgram_pos <- pg$pgram_raw[grid$pos_idx[seq_len(Jmax)]]`),
# so `donor_pool` values are valid direct indices into the vectors built
# here. `partition$J1`/`partition$J0` are already given in these same
# grid-array coordinates.
.fs_periodogram_data <- function(tst) {
  grid     <- tst$partition$grid
  omega_pi <- grid$omega_pos / pi
  I_pos    <- tst$pgram$pgram_std[grid$pos_idx]
  J1       <- tst$partition$J1
  J0       <- tst$partition$J0
  donor    <- tst$donor_pool
  breaks_pi <- tst$partition$breaks / pi

  thresh <- if (length(J0) > 0L && is.finite(tst$evt$critical)) {
    max(I_pos[J0]) + tst$evt$critical
  } else {
    NA_real_
  }

  group <- rep("other", length(omega_pi))
  if (length(donor) > 0L) group[donor] <- "donor"
  if (length(J1) > 0L)    group[J1]    <- "seasonal"

  list(omega_pi = omega_pi, I = I_pos, group = group,
       breaks_pi = breaks_pi, thresh = thresh)
}

# Base-graphics rendering of the standardized periodogram view shared by
# plot.seas_test(), the seas_sa "periodogram" type, and the seas_collection
# small multiples. `mini = TRUE` renders a compact, axis-free panel (no bin
# boundaries, no legend) for the small-multiples grid.
.fs_draw_periodogram <- function(tst, main = NULL, mini = FALSE) {
  d <- .fs_periodogram_data(tst)
  ylim <- range(c(d$I, d$thresh), finite = TRUE)
  ylim[2] <- ylim[2] * 1.08

  if (mini) {
    graphics::plot(NA, xlim = c(0, 1), ylim = ylim, xlab = "", ylab = "",
                   main = main, axes = FALSE, cex.main = 0.8)
    graphics::box(col = .fs_pal$axis)
  } else {
    .fs_panel(xlim = c(0, 1), ylim = ylim,
              xlab = expression(omega / pi), ylab = "standardized periodogram",
              main = if (is.null(main)) {
                "Standardized periodogram (whitened residuals)"
              } else main)
    graphics::abline(v = d$breaks_pi, col = .fs_pal$axis, lty = 3, lwd = 1)
  }

  is_other <- d$group == "other"
  is_donor <- d$group == "donor"
  is_seas  <- d$group == "seasonal"
  pt_cex <- if (mini) 0.5 else 0.7
  graphics::points(d$omega_pi[is_other], d$I[is_other], col = .fs_pal$muted,
                   pch = 16, cex = pt_cex)
  graphics::points(d$omega_pi[is_donor], d$I[is_donor], col = .fs_pal$muted,
                   pch = 16, cex = pt_cex + 0.2)
  graphics::points(d$omega_pi[is_seas], d$I[is_seas], col = .fs_pal$blue,
                   pch = 16, cex = pt_cex + 0.2)

  if (is.finite(d$thresh)) {
    graphics::abline(h = d$thresh, col = .fs_pal$orange,
                     lty = 2, lwd = if (mini) 1 else 1.4)
  }

  if (!mini) {
    graphics::legend("topright", bty = "n", text.col = .fs_pal$ink2,
                     pch = c(16, 16, NA), lty = c(NA, NA, 2),
                     col = c(.fs_pal$blue, .fs_pal$muted, .fs_pal$orange),
                     legend = c("seasonal-bin ordinate", "donor ordinate",
                                "EVT critical level"))
  }
  invisible(NULL)
}

# Base-graphics rendering of the series view (original vs adjusted, plus the
# removed seasonal component) shared by plot.seas_sa() and the
# seas_collection small multiples. `two_panel = FALSE` (used for the
# small-multiples grid) draws only the overlay panel, without the seasonal
# component beneath it, and without axes/legend.
.fs_draw_series <- function(sa, main = NULL, mini = FALSE, two_panel = TRUE) {
  tst  <- sa$test
  orig <- .fs_as_numeric_series(tst$x)
  adj  <- .fs_as_numeric_series(sa$adjusted)
  seas <- .fs_as_numeric_series(sa$seasonal)
  tt <- if (!is.null(tst$tsp)) {
    seq(tst$tsp[1L], by = 1 / tst$tsp[3L], length.out = length(orig))
  } else {
    seq_along(orig)
  }

  if (mini) {
    ylim <- range(c(orig, adj), finite = TRUE)
    graphics::plot(tt, orig, type = "l", col = .fs_pal$muted, xlab = "",
                   ylab = "", main = main, axes = FALSE, cex.main = 0.8,
                   ylim = ylim)
    graphics::box(col = .fs_pal$axis)
    graphics::lines(tt, adj, col = .fs_pal$blue)
    return(invisible(NULL))
  }

  .fs_panel(xlim = range(tt), ylim = range(c(orig, adj), finite = TRUE),
            xlab = "time", ylab = "value",
            main = if (is.null(main)) "Original vs adjusted" else main)
  graphics::lines(tt, orig, col = .fs_pal$muted, lwd = 1.6)
  graphics::lines(tt, adj, col = .fs_pal$blue, lwd = 1.6)
  graphics::legend("topleft", bty = "n", text.col = .fs_pal$ink2, lwd = 1.6,
                   col = c(.fs_pal$muted, .fs_pal$blue),
                   legend = c("original", "adjusted"))

  if (two_panel) {
    .fs_panel(xlim = range(tt), ylim = range(seas, finite = TRUE),
              xlab = "time", ylab = "seasonal",
              main = "Removed seasonal component")
    graphics::lines(tt, seas, col = .fs_pal$orange, lwd = 1.6)
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# plot.seas_test
# ---------------------------------------------------------------------------

#' Plot a `seas_test` object
#'
#' The standardized periodogram of the whitened series against
#' \eqn{\omega/\pi}: seasonal-bin ordinates highlighted in blue, donor
#' ordinates in gray, the EVT rejection threshold as a dashed reference line
#' (`max(I[J0]) + evt$critical`, the level a seasonal-bin ordinate would
#' need to exceed to trigger detection at `alpha`), and dotted vertical
#' lines at the equal-width bin boundaries.
#'
#' @param x A `seas_test` object.
#' @param type Plot type. Only `"periodogram"` is currently implemented.
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#'
#' @examples
#' \donttest{
#' set.seed(2026)
#' n <- 120L; Phi <- 0.9
#' e <- stats::rnorm(n + 400)
#' z <- numeric(n + 400)
#' for (t in 6:(n + 400)) z[t] <- 0.5 * z[t - 1] + Phi * z[t - 4] -
#'   0.5 * Phi * z[t - 5] + e[t]
#' tst <- seas_test(z[401:(n + 400)], frequency = 4)
#' plot(tst)
#' }
#'
#' @export
plot.seas_test <- function(x, type = c("periodogram"), ...) {
  type <- match.arg(type)
  op <- graphics::par(bg = .fs_pal$surface, mar = c(4.0, 4.2, 2.6, 1.2),
                      col.axis = .fs_pal$ink2, col.lab = .fs_pal$ink2,
                      col.main = .fs_pal$ink, cex.main = 1.05,
                      mgp = c(2.4, 0.6, 0), tcl = -0.3)
  on.exit(graphics::par(op))
  .fs_draw_periodogram(x)
  invisible(x)
}

# ---------------------------------------------------------------------------
# plot.seas_sa
# ---------------------------------------------------------------------------

.plot_seas_sa_series <- function(x, ...) {
  op <- graphics::par(mfrow = c(2, 1), bg = .fs_pal$surface,
                      mar = c(3.0, 3.8, 2.2, 1.0), col.axis = .fs_pal$ink2,
                      col.lab = .fs_pal$ink2, col.main = .fs_pal$ink,
                      cex.main = 1.05, mgp = c(2.0, 0.6, 0), tcl = -0.3)
  on.exit(graphics::par(op))
  .fs_draw_series(x, two_panel = TRUE)
  invisible(x)
}

.plot_seas_sa_periodogram <- function(x, ...) {
  tst <- x$test
  d_forced <- if (isTRUE(tst$whitener$d == 1L)) "first" else "none"
  post <- tryCatch(
    seas_test(.fs_as_numeric_series(x$adjusted), frequency = tst$N,
              P = tst$P, M = tst$M, alpha = tst$alpha, d = d_forced,
              ar_max = tst$ar_max),
    error = function(e) NULL
  )
  op <- graphics::par(mfrow = c(1, 2), bg = .fs_pal$surface,
                      mar = c(4.0, 4.2, 2.6, 1.0), col.axis = .fs_pal$ink2,
                      col.lab = .fs_pal$ink2, col.main = .fs_pal$ink,
                      cex.main = 1.0, mgp = c(2.4, 0.6, 0), tcl = -0.3)
  on.exit(graphics::par(op))
  .fs_draw_periodogram(tst, main = "Before adjustment")
  if (!is.null(post)) {
    .fs_draw_periodogram(post, main = "After adjustment")
  } else {
    graphics::plot.new()
    graphics::title(main = "After adjustment (periodogram unavailable)")
  }
  invisible(x)
}

#' Plot a `seas_sa` object
#'
#' `type = "series"` (default): a 2-panel layout with the original series
#' (gray) and the adjusted series (blue) overlaid on top, and the removed
#' seasonal component (orange) beneath. `type = "periodogram"`: the
#' before-adjustment and after-adjustment standardized periodograms
#' side-by-side, drawn exactly as by [plot.seas_test()] (the "after" panel
#' re-runs [seas_test()] on the adjusted series so it can be rendered the
#' same way; if that re-test errors, e.g. on a degenerate adjusted series,
#' the panel is left blank with a note).
#'
#' @param x A `seas_sa` object.
#' @param type Plot type: `"series"` (default) or `"periodogram"`.
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#'
#' @examples
#' \donttest{
#' set.seed(2026)
#' n <- 120L; Phi <- 0.9
#' e <- stats::rnorm(n + 400)
#' z <- numeric(n + 400)
#' for (t in 6:(n + 400)) z[t] <- 0.5 * z[t - 1] + Phi * z[t - 4] -
#'   0.5 * Phi * z[t - 5] + e[t]
#' tst <- seas_test(z[401:(n + 400)], frequency = 4)
#' adj <- seas_ssi(tst, phase_rule = "minimum")
#' plot(adj)
#' plot(adj, type = "periodogram")
#' }
#'
#' @export
plot.seas_sa <- function(x, type = c("series", "periodogram"), ...) {
  type <- match.arg(type)
  switch(type,
    series      = .plot_seas_sa_series(x, ...),
    periodogram = .plot_seas_sa_periodogram(x, ...)
  )
  invisible(x)
}

# ---------------------------------------------------------------------------
# plot.seas_collection
# ---------------------------------------------------------------------------

#' Plot a `seas_collection` object
#'
#' A grid of small multiples, one panel per key: for a `type = "sa"`
#' collection, the original (gray) vs adjusted (blue) series overlay; for a
#' `type = "test"` collection, the mini standardized periodogram (seasonal
#' ordinates in blue). At most 12 keys are plotted; if `x` has more, a
#' warning is issued and only the first 12 (in `x`'s row order) are shown.
#'
#' @param x A `seas_collection`.
#' @param ... Ignored.
#' @param max_keys Maximum number of keys to plot. Default `12L`.
#'
#' @return `x`, invisibly.
#'
#' @examples
#' \dontrun{
#' plot(seas_test(some_keyed_tsibble))
#' }
#'
#' @export
plot.seas_collection <- function(x, ..., max_keys = 12L) {
  # `type` and `key_cols` are read BEFORE any row-subsetting below: `x[i, ]`
  # (tibble row-indexing) drops the custom `type` attribute and the
  # `seas_collection` class, so reading `attr(x, "type")` afterwards would
  # silently see `NULL` and mis-dispatch a >12-key "sa" collection to the
  # periodogram branch instead of the series-overlay branch.
  type     <- attr(x, "type")
  key_cols <- setdiff(names(x), "fit")

  n <- nrow(x)
  if (n > max_keys) {
    warning(sprintf(
      "`x` has %d keys; plotting only the first %d (max_keys = %d).",
      n, max_keys, max_keys))
    x <- x[seq_len(max_keys), ]
    n <- max_keys
  }

  key_label <- function(i) {
    if (length(key_cols) == 0L) return(paste0("series ", i))
    paste(unlist(x[i, key_cols]), collapse = "/")
  }

  mfrow <- grDevices::n2mfrow(n)
  op <- graphics::par(mfrow = mfrow, bg = .fs_pal$surface,
                      mar = c(2.4, 2.6, 1.8, 0.6), col.axis = .fs_pal$ink2,
                      col.main = .fs_pal$ink, cex.main = 0.85,
                      mgp = c(1.4, 0.4, 0), tcl = -0.25)
  on.exit(graphics::par(op))

  for (i in seq_len(n)) {
    fit <- x$fit[[i]]
    if (identical(type, "sa")) {
      .fs_draw_series(fit, main = key_label(i), mini = TRUE)
    } else {
      .fs_draw_periodogram(fit, main = key_label(i), mini = TRUE)
    }
  }
  invisible(x)
}

# ---------------------------------------------------------------------------
# autoplot() methods (ggplot2, Suggests -- conditional registration)
# ---------------------------------------------------------------------------
#
# ggplot2 is a Suggests dependency, not an Imports dependency: freqseas must
# load and work without it. Each method below therefore starts with
# `rlang::check_installed("ggplot2")` (an informative error if it is
# missing) and is registered via the `@exportS3Method ggplot2::autoplot`
# roxygen tag rather than `@importFrom ggplot2 autoplot` + `@export`: the
# latter would force ggplot2 into Imports. `@exportS3Method pkg::generic`
# (roxygen2 >= 7) emits `S3method(ggplot2::autoplot, <class>)` directly into
# NAMESPACE, which R registers for delayed dispatch -- activated only if/when
# a caller has ggplot2 loaded -- without freqseas ever attaching or importing
# it itself.

#' Plot a `seas_test` object with ggplot2
#'
#' The ggplot2 equivalent of [plot.seas_test()]: same content (standardized
#' periodogram, seasonal-bin ordinates in blue, donor ordinates in gray, the
#' EVT rejection threshold, dotted bin boundaries), `ggplot2::theme_minimal()`
#' chrome. Requires the ggplot2 package (Suggests).
#'
#' @param object A `seas_test` object.
#' @param type Plot type. Only `"periodogram"` is currently implemented.
#' @param ... Ignored.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(2026)
#'   n <- 120L; Phi <- 0.9
#'   e <- stats::rnorm(n + 400)
#'   z <- numeric(n + 400)
#'   for (t in 6:(n + 400)) z[t] <- 0.5 * z[t - 1] + Phi * z[t - 4] -
#'     0.5 * Phi * z[t - 5] + e[t]
#'   tst <- seas_test(z[401:(n + 400)], frequency = 4)
#'   ggplot2::autoplot(tst)
#' }
#' }
#'
#' @exportS3Method ggplot2::autoplot
autoplot.seas_test <- function(object, type = c("periodogram"), ...) {
  rlang::check_installed("ggplot2", reason = "to use `autoplot()` on a `seas_test`.")
  type <- match.arg(type)
  d <- .fs_periodogram_data(object)
  df <- tibble::tibble(
    omega_pi = d$omega_pi, I = d$I,
    group = factor(d$group, levels = c("other", "donor", "seasonal"))
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$omega_pi, y = .data$I,
                                        color = .data$group)) +
    ggplot2::geom_vline(xintercept = d$breaks_pi, linetype = "dotted",
                        color = .fs_pal$axis) +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(
      values = c(other = .fs_pal$muted, donor = .fs_pal$muted,
                 seasonal = .fs_pal$blue),
      breaks = c("seasonal", "donor"),
      labels = c("seasonal-bin ordinate", "donor ordinate")
    ) +
    ggplot2::labs(x = expression(omega / pi), y = "standardized periodogram",
                  title = "Standardized periodogram (whitened residuals)",
                  color = NULL) +
    ggplot2::theme_minimal()

  if (is.finite(d$thresh)) {
    p <- p + ggplot2::geom_hline(yintercept = d$thresh, linetype = "dashed",
                                 color = .fs_pal$orange)
  }
  p
}

.autoplot_seas_sa_series <- function(object) {
  tst  <- object$test
  orig <- .fs_as_numeric_series(tst$x)
  adj  <- .fs_as_numeric_series(object$adjusted)
  seas <- .fs_as_numeric_series(object$seasonal)
  tt <- if (!is.null(tst$tsp)) {
    seq(tst$tsp[1L], by = 1 / tst$tsp[3L], length.out = length(orig))
  } else {
    seq_along(orig)
  }

  df_top <- tibble::tibble(
    time  = rep(tt, 2L), value = c(orig, adj),
    series = factor(rep(c("original", "adjusted"), each = length(orig)),
                    levels = c("original", "adjusted")),
    panel = "series"
  )
  df_bot <- tibble::tibble(time = tt, value = seas, panel = "seasonal component")

  panel_levels <- c("series", "seasonal component")
  df_top$panel <- factor(df_top$panel, levels = panel_levels)
  df_bot$panel <- factor(df_bot$panel, levels = panel_levels)

  ggplot2::ggplot() +
    ggplot2::geom_line(data = df_top,
                       ggplot2::aes(x = .data$time, y = .data$value,
                                    color = .data$series)) +
    ggplot2::geom_line(data = df_bot,
                       ggplot2::aes(x = .data$time, y = .data$value),
                       color = .fs_pal$orange) +
    ggplot2::facet_wrap(~ .data$panel, ncol = 1, scales = "free_y") +
    ggplot2::scale_color_manual(values = c(original = .fs_pal$muted,
                                           adjusted = .fs_pal$blue)) +
    ggplot2::labs(x = "time", y = NULL, color = NULL) +
    ggplot2::theme_minimal()
}

.autoplot_seas_sa_periodogram <- function(object) {
  tst <- object$test
  d_forced <- if (isTRUE(tst$whitener$d == 1L)) "first" else "none"
  post <- tryCatch(
    seas_test(.fs_as_numeric_series(object$adjusted), frequency = tst$N,
              P = tst$P, M = tst$M, alpha = tst$alpha, d = d_forced,
              ar_max = tst$ar_max),
    error = function(e) NULL
  )

  d_before  <- .fs_periodogram_data(tst)
  df_before <- tibble::tibble(omega_pi = d_before$omega_pi, I = d_before$I,
                              stage = "before")
  if (!is.null(post)) {
    d_after  <- .fs_periodogram_data(post)
    df_after <- tibble::tibble(omega_pi = d_after$omega_pi, I = d_after$I,
                               stage = "after")
    df <- dplyr::bind_rows(df_before, df_after)
  } else {
    df <- df_before
  }
  df$stage <- factor(df$stage, levels = c("before", "after"))

  ggplot2::ggplot(df, ggplot2::aes(x = .data$omega_pi, y = .data$I,
                                   color = .data$stage)) +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::scale_color_manual(values = c(before = .fs_pal$muted,
                                           after = .fs_pal$blue)) +
    ggplot2::labs(x = expression(omega / pi), y = "standardized periodogram",
                  color = NULL,
                  title = "Standardized periodogram: before vs after") +
    ggplot2::theme_minimal()
}

#' Plot a `seas_sa` object with ggplot2
#'
#' The ggplot2 equivalent of [plot.seas_sa()]: same content (`type =
#' "series"`: original vs adjusted overlaid, seasonal component beneath, as
#' facets; `type = "periodogram"`: before/after standardized periodograms),
#' `ggplot2::theme_minimal()` chrome. Requires the ggplot2 package (Suggests).
#'
#' @param object A `seas_sa` object.
#' @param type Plot type: `"series"` (default) or `"periodogram"`.
#' @param ... Ignored.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(2026)
#'   n <- 120L; Phi <- 0.9
#'   e <- stats::rnorm(n + 400)
#'   z <- numeric(n + 400)
#'   for (t in 6:(n + 400)) z[t] <- 0.5 * z[t - 1] + Phi * z[t - 4] -
#'     0.5 * Phi * z[t - 5] + e[t]
#'   tst <- seas_test(z[401:(n + 400)], frequency = 4)
#'   adj <- seas_ssi(tst, phase_rule = "minimum")
#'   ggplot2::autoplot(adj)
#' }
#' }
#'
#' @exportS3Method ggplot2::autoplot
autoplot.seas_sa <- function(object, type = c("series", "periodogram"), ...) {
  rlang::check_installed("ggplot2", reason = "to use `autoplot()` on a `seas_sa`.")
  type <- match.arg(type)
  switch(type,
    series      = .autoplot_seas_sa_series(object),
    periodogram = .autoplot_seas_sa_periodogram(object)
  )
}

#' Plot a `seas_collection` object with ggplot2
#'
#' The ggplot2 equivalent of [plot.seas_collection()]: one facet per key
#' (`ggplot2::facet_wrap()`), the same per-key content as [plot.seas_sa()]'s
#' `"series"` type (for a `type = "sa"` collection) or [plot.seas_test()]'s
#' periodogram (for a `type = "test"` collection). At most `max_keys` keys
#' are plotted; if `object` has more, a warning is issued and only the first
#' `max_keys` (in `object`'s row order) are shown. Requires the ggplot2
#' package (Suggests).
#'
#' @param object A `seas_collection`.
#' @param ... Ignored.
#' @param max_keys Maximum number of keys to plot. Default `12L`.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   ggplot2::autoplot(seas_test(some_keyed_tsibble))
#' }
#' }
#'
#' @exportS3Method ggplot2::autoplot
autoplot.seas_collection <- function(object, ..., max_keys = 12L) {
  rlang::check_installed("ggplot2", reason = "to use `autoplot()` on a `seas_collection`.")

  # Read BEFORE any row-subsetting below -- see the matching note in
  # plot.seas_collection(): tibble row-indexing drops the `type` attribute
  # and the `seas_collection` class.
  type     <- attr(object, "type")
  key_cols <- setdiff(names(object), "fit")

  n <- nrow(object)
  if (n > max_keys) {
    warning(sprintf(
      "`object` has %d keys; plotting only the first %d (max_keys = %d).",
      n, max_keys, max_keys))
    object <- object[seq_len(max_keys), ]
    n <- max_keys
  }

  labels <- vapply(seq_len(n), function(i) {
    if (length(key_cols) == 0L) return(paste0("series ", i))
    paste(unlist(object[i, key_cols]), collapse = "/")
  }, character(1))

  if (identical(type, "sa")) {
    rows <- lapply(seq_len(n), function(i) {
      sa   <- object$fit[[i]]
      tst  <- sa$test
      orig <- .fs_as_numeric_series(tst$x)
      adj  <- .fs_as_numeric_series(sa$adjusted)
      tt <- if (!is.null(tst$tsp)) {
        seq(tst$tsp[1L], by = 1 / tst$tsp[3L], length.out = length(orig))
      } else {
        seq_along(orig)
      }
      tibble::tibble(
        key = labels[i], time = rep(tt, 2L), value = c(orig, adj),
        series = factor(rep(c("original", "adjusted"), each = length(orig)),
                        levels = c("original", "adjusted"))
      )
    })
    df <- dplyr::bind_rows(rows)
    ggplot2::ggplot(df, ggplot2::aes(x = .data$time, y = .data$value,
                                     color = .data$series)) +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(~ .data$key, scales = "free") +
      ggplot2::scale_color_manual(values = c(original = .fs_pal$muted,
                                             adjusted = .fs_pal$blue)) +
      ggplot2::labs(x = "time", y = NULL, color = NULL) +
      ggplot2::theme_minimal()
  } else {
    rows <- lapply(seq_len(n), function(i) {
      tst <- object$fit[[i]]
      d <- .fs_periodogram_data(tst)
      tibble::tibble(
        key = labels[i], omega_pi = d$omega_pi, I = d$I,
        group = factor(d$group, levels = c("other", "donor", "seasonal"))
      )
    })
    df <- dplyr::bind_rows(rows)
    ggplot2::ggplot(df, ggplot2::aes(x = .data$omega_pi, y = .data$I,
                                     color = .data$group)) +
      ggplot2::geom_point() +
      ggplot2::facet_wrap(~ .data$key, scales = "free_y") +
      ggplot2::scale_color_manual(
        values = c(other = .fs_pal$muted, donor = .fs_pal$muted,
                   seasonal = .fs_pal$blue),
        breaks = c("seasonal", "donor"),
        labels = c("seasonal-bin ordinate", "donor ordinate")
      ) +
      ggplot2::labs(x = expression(omega / pi), y = "standardized periodogram",
                    color = NULL) +
      ggplot2::theme_minimal()
  }
}
