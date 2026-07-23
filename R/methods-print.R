# methods-print.R
# Compact, informative print methods for the two single-series objects,
# following the print contract in notes/2026-07-22_freqseas_architecture.md
# sec.4: decision + p, spec with shoulder p and phase R, the phase rule tagged
# "[declared identification]", M and how it was chosen, the whitener line, and
# (for seas_sa) the post-adjustment detection p. summary()/plot() come in a
# later wave.

.fmt_p <- function(p) {
  if (is.null(p) || length(p) != 1L || is.na(p)) return("NA")
  sprintf("%.3g", p)
}

.fmt_R <- function(R) {
  if (is.null(R) || length(R) != 1L || is.na(R)) return("NA")
  sprintf("%.2f", R)
}

# The whitener line, with the seasonal-alignment trim appended only when the
# residual grid actually had to be shortened (n_trim > 0). The trim is a real,
# bounded edge effect: the leading n_trim + d + p observations of an adjusted
# series come back unadjusted, so it is reported rather than left implicit.
.whitener_line <- function(wh) {
  base <- sprintf("whitener: d = %d, AR(%d) on M0 ordinates", wh$d, wh$p)
  n_trim <- if (is.null(wh$n_trim)) 0L else wh$n_trim
  if (n_trim > 0L) {
    base <- sprintf("%s | seasonal alignment trimmed %d leading residual%s (first %d obs. unadjusted)",
                    base, n_trim, if (n_trim == 1L) "" else "s",
                    n_trim + wh$d + wh$p)
  }
  base
}

.M_line <- function(obj) {
  sprintf("M = %d (%s)", obj$M, obj$M_selection$descriptor)
}

#' Print a `seas_test` object
#'
#' Compact one-screen summary: the detection decision and p-value, the
#' line/band specification with its shoulder p-value and block-phase resultant,
#' the bin count and how it was chosen, and the whitener line.
#'
#' @param x A `seas_test` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.seas_test <- function(x, ...) {
  cat("freqseas seasonality test\n")
  decision_str <- if (isTRUE(x$decision)) "seasonal" else "not seasonal"
  cat(sprintf(" decision: %s (detection p = %s, alpha = %s)\n",
              decision_str, .fmt_p(x$evt$p), format(x$alpha)))
  cat(sprintf(" spec: %s (shoulder p = %s, phase R = %s)\n",
              x$spec$label, .fmt_p(x$spec$shoulder_p), .fmt_R(x$spec$phase_R)))
  cat(sprintf(" %s | %s\n", .M_line(x), .whitener_line(x$whitener)))
  if (!is.na(x$spec$note)) cat(sprintf(" note: %s\n", x$spec$note))
  # wave-3 addition (comb diagnostic): flag excess concentrated at a strict
  # subset of the seasonal harmonics (see seas_test.R's .harmonic_excess_table).
  if (!is.null(x$spec$comb_note) && !is.na(x$spec$comb_note)) {
    cat(sprintf(" comb: %s\n", x$spec$comb_note))
  }
  invisible(x)
}

#' Print a `seas_sa` object
#'
#' Compact one-screen summary of the adjustment: the embedded detection p and
#' specification, the declared phase rule (tagged `[declared identification]`),
#' the bin count and whitener line, and the post-adjustment detection p-value
#' (the self-check).
#'
#' @param x A `seas_sa` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.seas_sa <- function(x, ...) {
  tst <- x$test
  if (identical(x$spec, "none")) {
    cat("SSI seasonal adjustment (none -- not seasonal, input unchanged)\n")
  } else {
    cat(sprintf("SSI seasonal adjustment (%s)\n", x$spec))
  }
  cat(sprintf(" detection p = %s | spec: %s (shoulder p = %s, phase R = %s)\n",
              .fmt_p(tst$evt$p), tst$spec$label,
              .fmt_p(tst$spec$shoulder_p), .fmt_R(tst$spec$phase_R)))
  cat(sprintf(" phase rule: %s  [declared identification]\n", x$phase_rule))
  cat(sprintf(" %s | %s\n", .M_line(tst), .whitener_line(tst$whitener)))
  cat(sprintf(" post-adjustment detection p = %s\n", .fmt_p(x$post_evt)))
  invisible(x)
}
