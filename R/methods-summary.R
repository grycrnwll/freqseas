# methods-summary.R
# Richer, structured summaries of `seas_test`/`seas_sa`/`seas_collection`
# objects, on top of the compact print() contract in R/methods-print.R.
# `summary()` returns an invisible, structured list with its own S3 class
# (`summary.seas_test` / `summary.seas_sa` / `summary.seas_collection`) and a
# print method, the standard R pattern (cf. `summary.lm`). Reuses the
# `.fmt_p()` / `.fmt_R()` formatting helpers defined in R/methods-print.R
# rather than redefining them.
#
# Per the build brief: `spec$harmonic_table` is a field a sibling module may
# add concurrently (a per-harmonic diagnostic table). It does not exist yet
# on any `seas_test` produced by the current R/seas_test.R. Every place that
# reads it below uses a plain `$` lookup (which returns `NULL` for a missing
# list element) and an explicit `is.null()` guard, so this file works
# unchanged whether or not the field is present.

# ---------------------------------------------------------------------------
# summary.seas_test
# ---------------------------------------------------------------------------

#' Summarize a `seas_test` object
#'
#' A richer summary than [print.seas_test()]: the full detection record
#' (statistic, critical value, seasonal/nonseasonal ordinate counts), the
#' specification record (shoulder p, block-phase resultant), the whitener
#' record (differencing decision + reason, AR order/coefficients, the full
#' M0-Whittle BIC path), the `M` selection record, the donor pool size, and
#' -- if present -- the per-harmonic table (`x$spec$harmonic_table`, added by
#' a concurrently-developed module; silently omitted when absent).
#'
#' @param object A `seas_test` object.
#' @param ... Ignored.
#'
#' @return An object of class `"summary.seas_test"`, invisibly: a list with
#'   fields `decision`, `p`, `alpha`, `evt` (`statistic`, `p`, `critical`,
#'   `N1`, `N0`), `spec` (`label`, `shoulder_p`, `phase_R`, `note`), `N`,
#'   `P`, `M`, `M_selection`, `whitener` (`d`, `d_reason`, `p`, `ar`,
#'   `bic_path`, `sigma2`, `n_trim` = leading residuals dropped to align the
#'   residual grid to a whole number of seasonal cycles),
#'   `harmonic_table` (`NULL` if not present on
#'   `object`), `donor_pool_size`, `n_e`, and `call`.
#'
#' @examples
#' set.seed(2026)
#' n <- 120L; Phi <- 0.9
#' e <- stats::rnorm(n + 400)
#' z <- numeric(n + 400)
#' for (t in 6:(n + 400)) z[t] <- 0.5 * z[t - 1] + Phi * z[t - 4] -
#'   0.5 * Phi * z[t - 5] + e[t]
#' tst <- seas_test(z[401:(n + 400)], frequency = 4)
#' s <- summary(tst)
#' s
#'
#' @export
summary.seas_test <- function(object, ...) {
  wh <- object$whitener
  out <- list(
    decision        = isTRUE(object$decision),
    p               = object$evt$p,
    alpha           = object$alpha,
    evt             = object$evt,
    spec            = object$spec,
    N               = object$N,
    P               = object$P,
    M               = object$M,
    M_selection     = object$M_selection,
    whitener        = list(
      d        = wh$d,
      d_reason = wh$d_reason,
      p        = wh$p,
      ar       = wh$ar,
      bic_path = wh$bic_path,
      sigma2   = wh$sigma2,
      n_trim   = if (is.null(wh$n_trim)) 0L else wh$n_trim
    ),
    harmonic_table  = object$spec$harmonic_table,
    donor_pool_size = length(object$donor_pool),
    n_e             = length(object$e),
    call            = object$call
  )
  class(out) <- "summary.seas_test"
  out
}

#' Print a `summary.seas_test` object
#'
#' @param x A `summary.seas_test` object, as returned by [summary.seas_test()].
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.summary.seas_test <- function(x, ...) {
  cat("freqseas seasonality test -- summary\n\n")

  cat("Detection\n")
  cat(sprintf("  decision    : %s\n",
              if (x$decision) "seasonal" else "not seasonal"))
  cat(sprintf("  p-value     : %s  (alpha = %s)\n",
              .fmt_p(x$p), format(x$alpha)))
  cat(sprintf("  statistic   : %s  (critical = %s)\n",
              format(x$evt$statistic, digits = 4),
              format(x$evt$critical, digits = 4)))
  cat(sprintf("  ordinates   : N1 = %d seasonal, N0 = %d nonseasonal\n",
              x$evt$N1, x$evt$N0))

  cat("\nSpecification\n")
  cat(sprintf("  spec        : %s\n", x$spec$label))
  cat(sprintf("  shoulder p  : %s\n", .fmt_p(x$spec$shoulder_p)))
  cat(sprintf("  phase R     : %s\n", .fmt_R(x$spec$phase_R)))
  if (!is.null(x$spec$note) && !is.na(x$spec$note)) {
    cat(sprintf("  note        : %s\n", x$spec$note))
  }

  cat("\nBins\n")
  cat(sprintf("  M = %d (%s)\n", x$M, x$M_selection$descriptor))
  if (identical(x$M_selection$source, "suggest_M")) {
    cat(sprintf("  rho_M = %s, offset_u = %s\n",
                format(x$M_selection$rho_M, digits = 3),
                format(x$M_selection$offset_u, digits = 3)))
  }

  cat("\nWhitener\n")
  cat(sprintf("  d = %d (%s)\n", x$whitener$d, x$whitener$d_reason))
  cat(sprintf("  AR order p = %d\n", x$whitener$p))
  if (x$whitener$p > 0L) {
    cat(sprintf("  AR coefficients: %s\n",
                paste(format(x$whitener$ar, digits = 3), collapse = ", ")))
  }
  if (isTRUE(x$whitener$n_trim > 0L)) {
    cat(sprintf("  seasonal alignment: trimmed %d leading residual%s (first %d observations unadjusted)\n",
                x$whitener$n_trim, if (x$whitener$n_trim == 1L) "" else "s",
                x$whitener$n_trim + x$whitener$d + x$whitener$p))
  }
  bp <- x$whitener$bic_path
  if (!is.null(bp)) {
    cat("  BIC path (M0-Whittle):\n")
    cat(sprintf("    %s\n",
                paste(names(bp), format(bp, digits = 4), sep = " = ",
                      collapse = "   ")))
  }

  cat(sprintf("\nDonor pool size: %d (of %d whitened observations)\n",
              x$donor_pool_size, x$n_e))

  if (!is.null(x$harmonic_table)) {
    cat("\nPer-harmonic table\n")
    print(x$harmonic_table)
  }

  invisible(x)
}

# ---------------------------------------------------------------------------
# summary.seas_sa
# ---------------------------------------------------------------------------

#' Gain-surgery summary (internal)
#'
#' Reduces a `seas_sa`'s `Ghat` field (a full-circle numeric vector for the
#' band branch, or scalar `1` for the line/none branch -- see
#' [new_seas_sa()]) to the two headline numbers the summary reports: the
#' maximum gain applied, and the count of ordinates the surgery actually
#' touched (`Ghat > 1`, out of the full ordinate count). For the line/none
#' branch there is no continuous gain surface (the line branch nulls the
#' exact-harmonic ordinates directly; the none branch removes nothing), so
#' `ordinates_touched`/`n_ordinates` are `NA_integer_`.
#'
#' @param Ghat A `seas_sa$Ghat` field.
#' @return A list with `max_Ghat`, `ordinates_touched`, `n_ordinates`.
#' @keywords internal
#' @noRd
.gain_summary <- function(Ghat) {
  if (length(Ghat) > 1L) {
    list(
      max_Ghat          = max(Ghat),
      ordinates_touched = sum(Ghat > 1 + 1e-9),
      n_ordinates       = length(Ghat)
    )
  } else {
    list(
      max_Ghat          = as.numeric(Ghat),
      ordinates_touched = NA_integer_,
      n_ordinates        = NA_integer_
    )
  }
}

#' Summarize a `seas_sa` object
#'
#' A richer summary than [print.seas_sa()]: the embedded [summary.seas_test()]
#' record, the declared phase rule, a gain-surgery summary (maximum gain
#' applied and the count of ordinates it touched -- `NA` for the line/none
#' branch, which has no continuous gain surface), and the post-adjustment
#' detection p-value (the self-check).
#'
#' @param object A `seas_sa` object.
#' @param ... Ignored.
#'
#' @return An object of class `"summary.seas_sa"`, invisibly: a list with
#'   fields `test_summary` (a `summary.seas_test`), `spec` (`"band"`,
#'   `"line"`, or `"none"`), `phase_rule`, `gain` (`max_Ghat`,
#'   `ordinates_touched`, `n_ordinates`), `post_evt`, and `call`.
#'
#' @examples
#' set.seed(2026)
#' n <- 120L; Phi <- 0.9
#' e <- stats::rnorm(n + 400)
#' z <- numeric(n + 400)
#' for (t in 6:(n + 400)) z[t] <- 0.5 * z[t - 1] + Phi * z[t - 4] -
#'   0.5 * Phi * z[t - 5] + e[t]
#' tst <- seas_test(z[401:(n + 400)], frequency = 4)
#' adj <- seas_ssi(tst, phase_rule = "minimum")
#' summary(adj)
#'
#' @export
summary.seas_sa <- function(object, ...) {
  out <- list(
    test_summary = summary.seas_test(object$test),
    spec         = object$spec,
    phase_rule   = object$phase_rule,
    gain         = .gain_summary(object$Ghat),
    post_evt     = object$post_evt,
    call         = object$call
  )
  class(out) <- "summary.seas_sa"
  out
}

#' Print a `summary.seas_sa` object
#'
#' @param x A `summary.seas_sa` object, as returned by [summary.seas_sa()].
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.summary.seas_sa <- function(x, ...) {
  cat("SSI seasonal adjustment -- summary\n\n")
  cat(sprintf("Specification applied : %s\n", x$spec))
  cat(sprintf("Phase rule             : %s  [declared identification]\n\n",
              x$phase_rule))

  cat("Gain surgery\n")
  if (is.na(x$gain$ordinates_touched)) {
    cat("  (line/none spec: no continuous gain surface -- the exact-harmonic\n")
    cat("   ordinate(s) are nulled directly rather than gain-corrected)\n")
    cat(sprintf("  Ghat = %s\n", format(x$gain$max_Ghat, digits = 4)))
  } else {
    cat(sprintf("  max Ghat          : %s\n",
                format(x$gain$max_Ghat, digits = 4)))
    cat(sprintf("  ordinates touched : %d / %d\n",
                x$gain$ordinates_touched, x$gain$n_ordinates))
  }

  cat(sprintf("\nPost-adjustment detection p = %s\n", .fmt_p(x$post_evt)))

  cat("\n--- embedded seas_test summary ---\n")
  print(x$test_summary)

  invisible(x)
}

# ---------------------------------------------------------------------------
# summary.seas_collection
# ---------------------------------------------------------------------------

#' Summarize a `seas_collection` object
#'
#' Maps [summary.seas_test()] or [summary.seas_sa()] (as appropriate to the
#' collection's `type` attribute -- see `new_seas_collection()`) over every
#' row's fit and bundles the results with the key columns.
#'
#' @param object A `seas_collection`.
#' @param ... Ignored.
#'
#' @return An object of class `"summary.seas_collection"`, invisibly: a list
#'   with `type` (`"test"`/`"sa"`), `keys` (a tibble of the key columns), and
#'   `summaries` (a list, one `summary.seas_test`/`summary.seas_sa` per row,
#'   in the same order as `keys`).
#'
#' @examples
#' \dontrun{
#' summary(seas_test(some_keyed_tsibble))
#' }
#'
#' @export
summary.seas_collection <- function(object, ...) {
  type     <- attr(object, "type")
  key_cols <- setdiff(names(object), "fit")
  out <- list(
    type      = type,
    keys      = tibble::as_tibble(object[key_cols]),
    summaries = lapply(object$fit, summary)
  )
  class(out) <- "summary.seas_collection"
  out
}

#' Print a `summary.seas_collection` object
#'
#' @param x A `summary.seas_collection` object, as returned by
#'   [summary.seas_collection()].
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.summary.seas_collection <- function(x, ...) {
  n <- nrow(x$keys)
  cat(sprintf("A seas_collection: %d series (%s) -- summary\n\n", n, x$type))
  for (i in seq_len(n)) {
    key_label <- if (ncol(x$keys) > 0L) {
      paste(unlist(x$keys[i, ]), collapse = "/")
    } else {
      paste0("series ", i)
    }
    cat(sprintf("== %s ==\n", key_label))
    print(x$summaries[[i]])
    cat("\n")
  }
  invisible(x)
}
