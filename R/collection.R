# collection.R -- the keyed layer: mapping seas_test()/seas_ssi()/
# seas_adjust() over the keys of a tsibble, and the seas_collection class
# that holds the per-key results.
#
# Ownership note: this file defines METHODS for the generics `seas_test()`,
# `seas_ssi()`, and `seas_adjust()` (owned by R/seas_test.R, R/seas_ssi.R,
# and R/seas_adjust.R respectively) plus `seas_collection` methods for the
# four accessor generics defined in R/accessors.R. It never defines those
# three generics itself; dispatch on `seas_test()`/`seas_ssi()`/
# `seas_adjust()` only works once the sibling modules that define them are
# loaded too.
#
# Design note (`template`): every seas_test produced here -- whether
# returned bare (a key-less tsibble) or embedded in a seas_collection row --
# carries an extra `template` field (the per-key sliced `tbl_ts`, as
# attached by `fs_extract()`) beyond the fields the class contract lists.
# This is the bridge that lets a `seas_test` with `input_class == "tbl_ts"`
# be turned back into a tsibble downstream: `fs_rebuild()` needs a template
# object, and nothing else in the seas_test contract carries one. `x` itself
# is left as the plain numeric vector `seas_test.default()` produced -- not
# eagerly rebuilt -- so the object stays safe to use numerically downstream
# (see whiten.R's documented expectation that a whitener's `x` is "already
# extracted from any ts/tsibble").
#
# This bridge is load-bearing, not speculative: R/seas_ssi.R's `.reclothe()`
# (the helper `seas_ssi.seas_test()` uses to re-clothe `adjusted`/`seasonal`)
# only special-cases `input_class == "ts"`; anything else, including our
# `"tbl_ts"` tag, falls through to a plain numeric vector. `fs_reclothe_sa()`
# below compensates for that gap at the two points this file controls
# (`seas_ssi.seas_collection()` and the bare-series branch of
# `seas_adjust.tbl_ts()`), using the `template` carried on each fit. The one
# path neither point covers is a *manual* two-stage call on a bare, unkeyed
# tsibble -- `seas_ssi(seas_test(unkeyed_tsibble, ...), phase_rule = ...)`
# called directly rather than through `seas_adjust()` -- which dispatches
# straight to `seas_ssi.seas_test()` with no method of ours in between; only
# a `tbl_ts` branch in `.reclothe()` itself (R/seas_ssi.R, not our file) can
# close that last gap. Flagged for the integration wave.

#' Re-clothe a `seas_sa`'s adjusted/seasonal fields into a template's class
#'
#' Compensates for `seas_ssi.seas_test()`'s `.reclothe()` helper (R/seas_ssi.R)
#' not knowing about `"tbl_ts"` input: it re-wraps `sa$adjusted` and
#' `sa$seasonal` (which `.reclothe()` leaves as plain numeric for any
#' `input_class` other than `"ts"`) via [fs_rebuild()], using `template` (as
#' attached to the originating `seas_test` by [seas_test.tbl_ts()]).
#'
#' Only rebuilds fields that are still plain numeric: if `.reclothe()` is
#' ever extended with its own `"tbl_ts"` branch, `sa$adjusted`/`sa$seasonal`
#' would already be `tbl_ts` objects, and re-running [fs_rebuild()] on those
#' (which coerces its first argument via `as.numeric()`) would error. This
#' keeps the compensation a safe no-op once/if that upstream gap closes.
#'
#' @param sa A `seas_sa` object.
#' @param template `NULL`, or the `template` object [fs_rebuild()] needs
#'   (typically `sa$test$template`). `NULL` is a no-op (returns `sa`
#'   unchanged), covering non-tsibble inputs.
#'
#' @return `sa`, with `adjusted`/`seasonal` re-clothed when `template` is not
#'   `NULL` and they are still plain numeric.
#'
#' @keywords internal
#' @noRd
fs_reclothe_sa <- function(sa, template) {
  if (is.null(template)) {
    return(sa)
  }
  if (is.numeric(sa$adjusted)) {
    sa$adjusted <- fs_rebuild(sa$adjusted, template)
  }
  if (is.numeric(sa$seasonal)) {
    sa$seasonal <- fs_rebuild(sa$seasonal, template)
  }
  sa
}

#' Construct a `seas_collection`
#'
#' Internal constructor for the `seas_collection` class: a tibble of key
#' columns plus a list-column `fit` (one fitted object -- a `seas_test` or a
#' `seas_sa` -- per key), tagged with a `type` attribute so methods (in
#' particular [print.seas_collection()]) know which.
#'
#' @param rows A tibble/data frame containing the key columns and a
#'   list-column named `fit`.
#' @param type `"test"` (rows hold `seas_test` fits) or `"sa"` (rows hold
#'   `seas_sa` fits).
#'
#' @return `rows`, with class `c("seas_collection", class(rows))` and a
#'   `type` attribute set to the matched `type`.
#'
#' @keywords internal
#' @noRd
new_seas_collection <- function(rows, type = c("test", "sa")) {
  type <- match.arg(type)
  if (!is.data.frame(rows) || !"fit" %in% names(rows)) {
    stop("`rows` must be a tibble/data frame with a list-column named `fit`.")
  }
  structure(
    rows,
    class = c("seas_collection", class(rows)),
    type  = type
  )
}

#' Test a tsibble for seasonality, mapping over keys
#'
#' The `tbl_ts` method for `seas_test()`: extracts one numeric series per key
#' (via `fs_extract()`), runs the single-series test
#' (`seas_test.default()`) on each, and returns either a bare `seas_test`
#' (an unkeyed tsibble holding a single series) or a `seas_collection` (a
#' keyed tsibble, one fit per key).
#'
#' @param x A `tbl_ts` (tsibble) with a regular index and exactly one
#'   measured variable.
#' @param frequency Ignored for `tbl_ts` input; the annual frequency is
#'   always inferred from the index (see `fs_extract()`). Present only so
#'   the method's signature matches the generic's.
#' @param ... Passed on to `seas_test.default()` for each series (e.g. `P`,
#'   `M`, `alpha`, `d`, `ar_max`, `band`).
#'
#' @return A `seas_test` (if `x` is unkeyed, i.e. has a single series) or a
#'   `seas_collection` of type `"test"` (if `x` has one or more key columns).
#'
#' @examples
#' \dontrun{
#' library(tsibble)
#' seas_test(some_keyed_tsibble)
#' }
#'
#' @export
seas_test.tbl_ts <- function(x, frequency = NULL, ...) {
  ex <- fs_extract(x, frequency = frequency)

  fits <- lapply(ex$series, function(s) {
    fit <- seas_test.default(s$x, frequency = s$N, ...)
    fit$input_class <- "tbl_ts"
    fit$template <- s$template
    fit
  })

  if (length(ex$series) == 1L && is.null(ex$series[[1]]$key)) {
    return(fits[[1]])
  }

  key_tbl <- dplyr::bind_rows(lapply(ex$series, function(s) s$key))
  rows <- dplyr::bind_cols(key_tbl, tibble::tibble(fit = fits))
  new_seas_collection(rows, type = "test")
}

#' Adjust a `seas_collection`, mapping over keys
#'
#' The `seas_collection` method for `seas_ssi()`: maps `seas_ssi()` over
#' each row's `fit` (a `seas_test`), producing a `seas_sa` per key. As with
#' the single-series `seas_ssi()`, `phase_rule` has no default and must be
#' supplied -- this method does not validate it itself; a missing or invalid
#' `phase_rule` surfaces as whatever error the per-series `seas_ssi()`
#' (`seas_ssi.seas_test()`) raises. Each row's `adjusted`/`seasonal` is then
#' re-clothed into the per-key tsibble via `fs_reclothe_sa()` (see the design
#' note at the top of this file).
#'
#' @param object A `seas_collection` of type `"test"`.
#' @param phase_rule `"minimum"` or `"zero"`; required, no default (see
#'   `seas_ssi()`).
#' @param ... Passed on to `seas_ssi()` for each row's fit.
#'
#' @return A `seas_collection` of type `"sa"`, one `seas_sa` per key.
#'
#' @examples
#' \dontrun{
#' tst <- seas_test(some_keyed_tsibble)
#' seas_ssi(tst, phase_rule = "minimum")
#' }
#'
#' @export
seas_ssi.seas_collection <- function(object, phase_rule, ...) {
  fits <- lapply(object$fit, function(fit) {
    sa <- seas_ssi(fit, phase_rule = phase_rule, ...)
    fs_reclothe_sa(sa, fit$template)
  })
  rows <- object
  rows$fit <- fits
  new_seas_collection(rows, type = "sa")
}

#' Adjust a tsibble in one call, mapping over keys
#'
#' The `tbl_ts` method for `seas_adjust()`: the one-call wrapper, `seas_ssi(
#' seas_test(x, ...), phase_rule = phase_rule)`, dispatched through
#' [seas_test.tbl_ts()] (and, for keyed input, [seas_ssi.seas_collection()],
#' which re-clothes each key's `adjusted`/`seasonal` itself). For an unkeyed
#' `x`, `seas_ssi()` dispatches straight to `seas_ssi.seas_test()` with no
#' method of ours in between, so this method re-clothes the bare result
#' itself via `fs_reclothe_sa()` (see the design note at the top of this
#' file).
#'
#' @param x A `tbl_ts` (tsibble) with a regular index and exactly one
#'   measured variable.
#' @param ... Passed on to `seas_test()` (e.g. `P`, `M`, `alpha`, `d`,
#'   `ar_max`, `band`).
#' @param phase_rule `"minimum"` (default) or `"zero"`; the declared
#'   identification passed to `seas_ssi()`.
#'
#' @return A `seas_sa` (if `x` is unkeyed) or a `seas_collection` of type
#'   `"sa"` (if `x` has one or more key columns).
#'
#' @examples
#' \dontrun{
#' library(tsibble)
#' seas_adjust(some_keyed_tsibble, phase_rule = "minimum")
#' }
#'
#' @export
seas_adjust.tbl_ts <- function(x, ..., phase_rule = "minimum") {
  result <- seas_ssi(seas_test(x, ...), phase_rule = phase_rule)
  if (inherits(result, "seas_sa")) {
    result <- fs_reclothe_sa(result, result$test$template)
  }
  result
}

#' Print a `seas_collection`
#'
#' A compact, one-row-per-key summary: the key columns, the seasonality
#' decision, the specification label, the bin count `M`, and (for a
#' `type = "sa"` collection) the declared `phase_rule`.
#'
#' @param x A `seas_collection`.
#' @param ... Passed on to the underlying tibble `print()` method.
#'
#' @return `x`, invisibly.
#'
#' @examples
#' \dontrun{
#' print(seas_test(some_keyed_tsibble))
#' }
#'
#' @export
print.seas_collection <- function(x, ...) {
  type <- attr(x, "type")
  cat(sprintf("A seas_collection: %d series (%s)\n", nrow(x), type))

  key_cols <- setdiff(names(x), "fit")
  dec <- lapply(x$fit, decision)

  m_val <- vapply(x$fit, function(f) {
    m <- if (identical(type, "sa")) f$test$M else f$M
    if (is.null(m)) NA_real_ else as.numeric(m)
  }, numeric(1))

  tbl <- tibble::as_tibble(x[key_cols])
  tbl$seasonal <- vapply(dec, function(d) {
    if (isTRUE(d$seasonal)) "seasonal" else "not seasonal"
  }, character(1))
  tbl$spec <- vapply(dec, function(d) as.character(d$spec), character(1))
  tbl$M    <- m_val

  if (identical(type, "sa")) {
    tbl$phase_rule <- vapply(
      x$fit, function(f) as.character(f$phase_rule), character(1)
    )
  }

  print(tbl, ...)
  invisible(x)
}

#' @rdname adjusted
#' @export
adjusted.seas_collection <- function(object, ...) {
  key_cols <- setdiff(names(object), "fit")
  out <- tibble::as_tibble(object[key_cols])
  out$adjusted <- lapply(object$fit, adjusted)
  out
}

#' @rdname seasonal
#' @export
seasonal.seas_collection <- function(object, ...) {
  key_cols <- setdiff(names(object), "fit")
  out <- tibble::as_tibble(object[key_cols])
  out$seasonal <- lapply(object$fit, seasonal)
  out
}

#' @rdname whitener
#' @export
whitener.seas_collection <- function(object, ...) {
  key_cols <- setdiff(names(object), "fit")
  out <- tibble::as_tibble(object[key_cols])
  out$whitener <- lapply(object$fit, whitener)
  out
}

#' @rdname decision
#' @export
decision.seas_collection <- function(object, ...) {
  key_cols <- setdiff(names(object), "fit")
  d <- dplyr::bind_rows(lapply(object$fit, decision))
  dplyr::bind_cols(tibble::as_tibble(object[key_cols]), d)
}
