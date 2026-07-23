# benchmark.R -- yearly-total benchmarking: the NIPA-style constraint that
# complete-period totals are invariant to seasonal adjustment. Port of
# helper_functions/yearly_total_matching.R under the freqseas API.
#
# Convention preserved from the port (documented in @details):
#   * PROPORTIONAL (multiplicative) rescaling: each complete window of the
#     adjusted series is multiplied by a single factor so its total equals the
#     corresponding window's total in the original. (The port multiplies by
#     `sums_orig / sums_adj`; freqseas keeps that convention -- not additive.)
#   * Windows are POSITIONAL, taken from index 1 onward, length = one period
#     (`frequency` observations for `period = "year"`); there is no calendar
#     awareness. For a series that starts on the first sub-period of a year
#     (as all `ssi_examples` do), positional windows coincide with calendar
#     years.
#   * A trailing INCOMPLETE window (when the length is not a whole number of
#     periods) is LEFT UNADJUSTED, with a warning.

#' Benchmark a series to complete-period totals of another
#'
#' @description
#' Proportionally rescales `adjusted` so that each complete period's total
#' matches the corresponding period's total in `original`. This is the
#' NIPA-style benchmarking constraint that **annual totals are invariant to
#' seasonal adjustment**: a seasonally adjusted series should redistribute
#' activity *within* a year without changing the year's total.
#'
#' The rescaling is proportional (multiplicative), a port of the
#' `yearly_total_matching()` helper: within each complete period, every value
#' of `adjusted` is multiplied by the single factor
#' `sum(original_window) / sum(adjusted_window)`, so the rescaled window's total
#' equals the original window's total exactly.
#'
#' @details
#' Windowing is **positional**: consecutive, non-overlapping windows of one
#' period's length are taken from index 1 onward, with no calendar awareness
#' (matching the port). For `period = "year"` the window length is the annual
#' frequency (`4` for quarterly, `12` for monthly). A series that begins on the
#' first sub-period of a year -- e.g. a quarterly `ts` starting in Q1 -- has
#' positional windows that coincide with calendar years.
#'
#' If the length of the series is not a whole number of periods, the trailing
#' **incomplete** window is left unadjusted and a warning is issued (there is
#' no complete total to match it to). There is no leading partial window,
#' because windows start at index 1.
#'
#' @param adjusted The series to rescale: a [stats::ts()] (frequency inferred)
#'   or a bare numeric vector (supply `frequency`).
#' @param original The series providing the target complete-period totals; same
#'   length as `adjusted`. A `ts` or numeric vector.
#' @param period The period whose totals are constrained. Currently only
#'   `"year"` (annual totals), the default.
#' @param frequency Observations per year. Required when `adjusted` is a bare
#'   numeric vector; ignored when `adjusted` is a `ts` (its own frequency is
#'   used). If both `adjusted` and `original` are `ts`, their frequencies must
#'   agree.
#'
#' @return The benchmarked series, in the same class as `adjusted` (a `ts` with
#'   `adjusted`'s `tsp` if `adjusted` was a `ts`, else a plain numeric vector).
#'   It carries a `"benchmark"` attribute: a list with `period`, `frequency`,
#'   `scale_factors` (one multiplicative factor per complete window), and
#'   `n_tail_unadjusted` (the number of trailing observations left unchanged).
#'
#' @section Errors:
#' `benchmark_totals()` errors informatively when `adjusted` and `original`
#' differ in length, when either contains non-finite values, when a bare
#' numeric `adjusted` is given with no `frequency`, when two `ts` inputs have
#' mismatched frequencies, when the series is shorter than one complete period,
#' or when a complete window of `adjusted` sums to zero (proportional
#' benchmarking is then undefined -- there is no factor that makes a zero total
#' match a nonzero one).
#'
#' @examples
#' # Constrain a seasonal adjustment so annual totals match the original's:
#' orig    <- ssi_examples$x[[1]]                    # filter-world series (ts)
#' adj     <- seas_adjust(orig, phase_rule = "minimum")
#' matched <- benchmark_totals(adjusted(adj), orig)  # ts in, ts out
#'
#' # Annual totals of `matched` now equal those of `orig` exactly:
#' N  <- stats::frequency(orig)
#' yo <- colSums(matrix(as.numeric(orig),    nrow = N))
#' ym <- colSums(matrix(as.numeric(matched), nrow = N))
#' max(abs(ym - yo))                                 # ~0
#'
#' @export
benchmark_totals <- function(adjusted, original, period = c("year"),
                             frequency = NULL) {
  period <- match.arg(period)

  # --- resolve the annual frequency (window length for period = "year") ----
  adj_is_ts <- stats::is.ts(adjusted)
  org_is_ts <- stats::is.ts(original)
  adj_freq  <- if (adj_is_ts) stats::frequency(adjusted) else NA_real_
  org_freq  <- if (org_is_ts) stats::frequency(original) else NA_real_
  if (adj_is_ts && org_is_ts && !isTRUE(all.equal(adj_freq, org_freq))) {
    stop("`adjusted` and `original` are both `ts` but have different ",
         "frequencies (", adj_freq, " vs ", org_freq, ").", call. = FALSE)
  }
  freq <- if (adj_is_ts) {
    adj_freq
  } else if (!is.null(frequency)) {
    frequency
  } else if (org_is_ts) {
    org_freq
  } else {
    stop("`frequency` must be supplied when `adjusted` is a bare numeric ",
         "vector (e.g. `frequency = 4` for quarterly). `ts` inputs infer it.",
         call. = FALSE)
  }
  N <- as.integer(round(freq))
  if (length(N) != 1L || is.na(N) || N < 1L) {
    stop("resolved period length must be a single positive integer; got ",
         freq, ".", call. = FALSE)
  }

  a <- as.numeric(adjusted)
  o <- as.numeric(original)
  n <- length(a)

  # --- degenerate-input guards (informative) -------------------------------
  if (n != length(o)) {
    stop("`adjusted` and `original` must have the same length (got ", n,
         " and ", length(o), ").", call. = FALSE)
  }
  if (any(!is.finite(a)) || any(!is.finite(o))) {
    stop("`adjusted`/`original` contain non-finite values (NA/NaN/Inf); ",
         "clean the series first.", call. = FALSE)
  }
  n_full <- n %/% N
  if (n_full < 1L) {
    stop("the series is shorter than one complete ", period, " (length ", n,
         " < frequency ", N, "); nothing to benchmark.", call. = FALSE)
  }

  tail_len <- n - n_full * N
  if (tail_len > 0L) {
    warning(sprintf(
      paste0("length %d is not a whole number of %ss (frequency = %d); the ",
             "last %d observation(s) are left unadjusted."),
      n, period, N, tail_len), call. = FALSE)
  }

  # --- proportional matching on the complete windows -----------------------
  full_idx  <- seq_len(n_full * N)
  mat_adj   <- matrix(a[full_idx], nrow = N)
  mat_orig  <- matrix(o[full_idx], nrow = N)
  sums_adj  <- colSums(mat_adj)
  sums_orig <- colSums(mat_orig)

  sf <- sums_orig / sums_adj
  if (any(!is.finite(sf))) {
    bad <- which(!is.finite(sf))
    stop("cannot proportionally benchmark: complete-", period, " total(s) of ",
         "`adjusted` at window(s) ", paste(bad, collapse = ", "),
         " are zero, so no rescaling factor exists.", call. = FALSE)
  }

  matched <- a
  mat_matched <- sweep(mat_adj, 2, sf, `*`)
  matched[full_idx] <- as.numeric(mat_matched)

  # --- re-clothe into the input class (ts keeps its tsp) -------------------
  out <- if (adj_is_ts) {
    stats::ts(matched, start = stats::start(adjusted),
              frequency = stats::frequency(adjusted))
  } else {
    matched
  }
  attr(out, "benchmark") <- list(
    period = period, frequency = N,
    scale_factors = sf, n_tail_unadjusted = tail_len
  )
  out
}
