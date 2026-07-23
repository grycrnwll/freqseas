# classes.R
# Internal S3 constructors and validators for the two user-facing objects the
# freqseas single-series workflow produces: `seas_test` (Stage 0 whiten + Stage
# A detect + Stage B specify) and `seas_sa` (Stage C adjust). The field lists
# below are the object-model CONTRACTS from notes/2026-07-22_freqseas_
# architecture.md sec.4; downstream accessors, print/summary/plot methods, and
# the seas_collection mapping all read these fields by name.

# ---------------------------------------------------------------------------
# seas_test
# ---------------------------------------------------------------------------

#' Construct a `seas_test` object (internal)
#'
#' Low-level constructor: assembles the list and stamps the class. No
#' computation, no coercion beyond what the caller supplies. Use
#' [seas_test()] for the user-facing entry point.
#'
#' @param x The original input series (class preserved: numeric or `ts`).
#' @param input_class Character scalar, the class of `x` on input
#'   (e.g. `"numeric"`, `"ts"`).
#' @param tsp_attr Numeric length-3 `tsp` vector when `x` is a `ts`, else
#'   `NULL`; used to re-clothe adjusted output.
#' @param N,P,M The seasonal period, fundamental multiplier, and bin count.
#' @param M_selection List recording how `M` was chosen (`source` is
#'   `"suggest_M"` or `"user"`, plus diagnostics).
#' @param partition List bundling the grid and bin/ordinate index sets
#'   (`grid`, `breaks`, `J1`, `J0`, `H`, `S`, `J1_full`, `J0_full`, `N1`,
#'   `N0`, `M1`, `M0`, `omega_G`, `band`).
#' @param whitener The `fs_whiten()` record (`d`, `d_reason`, `mu`, `p`,
#'   `ar`, `bic_path`, ...).
#' @param e Numeric vector of whitened residuals (length `n_e`).
#' @param E Complex vector, the raw DFT of `e`.
#' @param pgram The `fs_periodogram()` record for `e`.
#' @param donor_pool Integer positions (positive-frequency-array coordinates)
#'   of the donor ordinates.
#' @param donor_level Numeric donor reference power level.
#' @param tau_hat Numeric null scale used by the EVT test.
#' @param evt List with the detection `statistic`, `p`, and `critical`.
#' @param spec List with the specification `label` (`"band"`/`"line"`),
#'   `shoulder_p`, `phase_R`, and an optional `note`.
#' @param decision Logical; `TRUE` when detection rejects the white null.
#' @param alpha Numeric test level.
#' @param ar_max Integer AR-order cap used by the whitener (carried for the
#'   post-adjustment self-check).
#' @param call The matched call.
#'
#' @return An object of class `"seas_test"`.
#' @keywords internal
#' @noRd
new_seas_test <- function(x, input_class, tsp_attr, N, P, M, M_selection,
                          partition, whitener, e, E, pgram, donor_pool,
                          donor_level, tau_hat, evt, spec, decision, alpha,
                          ar_max, call) {
  structure(
    list(
      x           = x,
      input_class = input_class,
      tsp         = tsp_attr,
      N           = N,
      P           = P,
      M           = M,
      M_selection = M_selection,
      partition   = partition,
      whitener    = whitener,
      e           = e,
      E           = E,
      pgram       = pgram,
      donor_pool  = donor_pool,
      donor_level = donor_level,
      tau_hat     = tau_hat,
      evt         = evt,
      spec        = spec,
      decision    = decision,
      alpha       = alpha,
      ar_max      = ar_max,
      call        = call
    ),
    class = "seas_test"
  )
}

#' Validate a `seas_test` object (internal)
#'
#' Checks the presence and basic shape of the contract fields. Raises an
#' informative error on the first violation; returns the object invisibly on
#' success.
#'
#' @param x A candidate `seas_test` object.
#' @return `x`, invisibly, if valid.
#' @keywords internal
#' @noRd
validate_seas_test <- function(x) {
  required <- c("x", "input_class", "N", "P", "M", "M_selection", "partition",
                "whitener", "e", "E", "pgram", "donor_pool", "tau_hat", "evt",
                "spec", "decision", "alpha", "call")
  missing <- setdiff(required, names(x))
  if (length(missing) > 0L) {
    stop("`seas_test` is missing required field(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  if (!is.complex(x$E)) stop("`seas_test$E` must be complex.", call. = FALSE)
  if (length(x$E) != length(x$e)) {
    stop("`seas_test$E` and `$e` must have the same length.", call. = FALSE)
  }
  if (!x$spec$label %in% c("band", "line")) {
    stop("`seas_test$spec$label` must be \"band\" or \"line\".", call. = FALSE)
  }
  if (!is.logical(x$decision) || length(x$decision) != 1L) {
    stop("`seas_test$decision` must be a single logical.", call. = FALSE)
  }
  for (fld in c("statistic", "p", "critical")) {
    if (is.null(x$evt[[fld]])) {
      stop("`seas_test$evt$", fld, "` is missing.", call. = FALSE)
    }
  }
  invisible(x)
}

# ---------------------------------------------------------------------------
# seas_sa
# ---------------------------------------------------------------------------

#' Construct a `seas_sa` object (internal)
#'
#' Low-level constructor for the seasonal-adjustment result. Use [seas_ssi()]
#' or [seas_adjust()] for the user-facing entry points.
#'
#' @param test The embedded `seas_test` object the adjustment was built from.
#' @param adjusted The adjusted series (input class preserved).
#' @param seasonal The removed seasonal component (`x - adjusted`, input class
#'   preserved).
#' @param Ghat The gain applied (a full-circle numeric vector for the band
#'   branch; scalar `1` when nothing was removed or for the line branch).
#' @param theta The phase applied (a full-circle numeric vector for the band
#'   branch under the minimum-phase rule; scalar `0` otherwise).
#' @param phase_rule Character, `"minimum"` or `"zero"`, the declared
#'   identification.
#' @param post_evt Numeric, the detection p-value of the re-test on the
#'   adjusted series (self-check).
#' @param spec Character, the specification branch used (`"band"`, `"line"`,
#'   or `"none"` for a not-seasonal input).
#' @param call The matched call.
#'
#' @return An object of class `"seas_sa"`.
#' @keywords internal
#' @noRd
new_seas_sa <- function(test, adjusted, seasonal, Ghat, theta, phase_rule,
                        post_evt, spec, call) {
  structure(
    list(
      test       = test,
      adjusted   = adjusted,
      seasonal   = seasonal,
      Ghat       = Ghat,
      theta      = theta,
      phase_rule = phase_rule,
      post_evt   = post_evt,
      spec       = spec,
      call       = call
    ),
    class = "seas_sa"
  )
}

#' Validate a `seas_sa` object (internal)
#'
#' @param x A candidate `seas_sa` object.
#' @return `x`, invisibly, if valid.
#' @keywords internal
#' @noRd
validate_seas_sa <- function(x) {
  required <- c("test", "adjusted", "seasonal", "Ghat", "theta", "phase_rule",
                "post_evt", "spec", "call")
  missing <- setdiff(required, names(x))
  if (length(missing) > 0L) {
    stop("`seas_sa` is missing required field(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  if (!inherits(x$test, "seas_test")) {
    stop("`seas_sa$test` must be a `seas_test` object.", call. = FALSE)
  }
  if (!x$phase_rule %in% c("minimum", "zero")) {
    stop("`seas_sa$phase_rule` must be \"minimum\" or \"zero\".", call. = FALSE)
  }
  if (length(x$adjusted) != length(x$seasonal)) {
    stop("`seas_sa$adjusted` and `$seasonal` must have equal length.",
         call. = FALSE)
  }
  invisible(x)
}
