# coerce.R -- input extraction and reconstruction for the three supported
# freqseas input classes (bare numeric, `ts`, and `tbl_ts`/tsibble).
#
# `fs_extract()` is the single entry point the rest of the package (in
# particular `seas_test.tbl_ts()` in R/collection.R) uses to turn any
# supported input into one or more plain numeric series plus enough context
# (`N`, `key`, `template`) to test/adjust each series and, later, rebuild
# output back into the input's original class via `fs_rebuild()`. Both
# functions are internal: they are the shared plumbing behind the exported
# `seas_test()`/`seas_ssi()`/`seas_adjust()` generics (owned elsewhere), not
# user-facing API.

#' Extract one or more numeric series from a supported freqseas input
#'
#' Normalizes a bare numeric vector, a [stats::ts()] object, or a keyed/
#' unkeyed tsibble (`tbl_ts`) into a common internal shape: a list of
#' per-series records, each carrying the numeric data, its inferred annual
#' frequency, an optional key (for keyed tsibbles), and a `template` object
#' that [fs_rebuild()] can later use to wrap adjusted/seasonal output back
#' into the input's original class.
#'
#' For a `tbl_ts`, the annual frequency `N` is inferred from the class of the
#' tsibble's index column: `tsibble::yearquarter()` indices give `N = 4`,
#' `tsibble::yearmonth()` indices give `N = 12`. `tsibble::yearweek()`
#' indices are explicitly unsupported (a calendar week is not an integer
#' fraction of a year), and any other index class (e.g. `Date`) errors with
#' guidance rather than guessing. The tsibble must have a regular index and
#' exactly one measured variable (select one first with `dplyr::select()` if
#' there are several); it is then split by key, one record per key
#' combination, in `tsibble::key_data()` order (tsibble's own canonical,
#' sorted key order -- this is what "key ordering" means here, since
#' tsibble does not preserve first-appearance order of key values).
#'
#' @param x A numeric vector, a [stats::ts()] object, or a `tbl_ts` (tsibble).
#' @param frequency Single positive number, observations per year. Required
#'   (and used) only when `x` is a bare numeric vector; ignored (frequency is
#'   always inferred from the object itself) for `ts` and `tbl_ts` inputs.
#'
#' @return A named list:
#'   \describe{
#'     \item{`series`}{A list of per-series records, each a list with
#'       `key` (`NULL` for an unkeyed input, or a one-row tibble of
#'       key-column values), `x` (numeric vector, the series data), `N`
#'       (single positive number, annual frequency), and `template` (`NULL`
#'       for a bare numeric input, else the object [fs_rebuild()] needs to
#'       wrap numeric output back into `x`'s original class: the original
#'       `ts` object, or the per-key sliced `tbl_ts`).}
#'     \item{`input_class`}{`"numeric"`, `"ts"`, or `"tbl_ts"`.}
#'   }
#'
#' @keywords internal
#' @noRd
fs_extract <- function(x, frequency = NULL) {
  if (inherits(x, "tbl_ts")) {
    return(fs_extract_tsibble(x))
  }
  if (stats::is.ts(x)) {
    return(fs_extract_ts(x))
  }
  if (is.numeric(x)) {
    return(fs_extract_numeric(x, frequency))
  }
  stop(
    "`x` must be a numeric vector, a `ts` object, or a `tbl_ts` (tsibble); ",
    "got an object of class `", paste(class(x), collapse = "/"), "`."
  )
}

#' Numeric branch of `fs_extract()`
#' @keywords internal
#' @noRd
fs_extract_numeric <- function(x, frequency) {
  if (is.null(frequency)) {
    stop(
      "`frequency` must be supplied for a bare numeric input (e.g. ",
      "`frequency = 12` for monthly data). `ts` and `tbl_ts` inputs infer ",
      "it automatically."
    )
  }
  if (length(frequency) != 1L || !is.finite(frequency) || frequency <= 0) {
    stop("`frequency` must be a single positive number.")
  }
  list(
    series = list(list(
      key      = NULL,
      x        = as.numeric(x),
      N        = frequency,
      template = NULL
    )),
    input_class = "numeric"
  )
}

#' `ts` branch of `fs_extract()`
#' @keywords internal
#' @noRd
fs_extract_ts <- function(x) {
  list(
    series = list(list(
      key      = NULL,
      x        = as.numeric(x),
      N        = stats::frequency(x),
      template = x
    )),
    input_class = "ts"
  )
}

#' Infer an annual frequency from a tsibble's index column
#' @keywords internal
#' @noRd
fs_index_frequency <- function(x) {
  idx <- x[[tsibble::index_var(x)]]
  if (inherits(idx, "yearquarter")) {
    return(4L)
  }
  if (inherits(idx, "yearmonth")) {
    return(12L)
  }
  if (inherits(idx, "yearweek")) {
    stop(
      "Weekly tsibble indices (`tsibble::yearweek()`) are unsupported: a ",
      "calendar week is not an integer fraction of a year, so an annual ",
      "seasonal frequency cannot be inferred. Aggregate to monthly or ",
      "quarterly first, or supply a `ts`/numeric input with an explicit ",
      "`frequency`."
    )
  }
  stop(
    "Cannot infer an annual frequency from an index of class `",
    paste(class(idx), collapse = "/"), "`. freqseas infers frequency from ",
    "`tsibble::yearquarter()` (N = 4) or `tsibble::yearmonth()` (N = 12) ",
    "indices. Convert the index first (e.g. `tsibble::yearmonth(date_col)`), ",
    "or supply a `ts`/numeric input with an explicit `frequency`."
  )
}

#' `tbl_ts` branch of `fs_extract()`
#' @keywords internal
#' @noRd
fs_extract_tsibble <- function(x) {
  if (!tsibble::is_regular(x)) {
    stop(
      "`x` has an irregular time index; freqseas requires a regular ",
      "tsibble (build one with `tsibble::fill_gaps()`, or supply a `ts`/",
      "numeric input instead)."
    )
  }

  N <- fs_index_frequency(x)

  mv <- tsibble::measured_vars(x)
  if (length(mv) != 1L) {
    if (length(mv) == 0L) {
      stop("`x` has no measured variables (only index/key columns).")
    }
    stop(
      "`x` has ", length(mv), " measured variables (",
      paste(mv, collapse = ", "), "); freqseas needs exactly one. Select ",
      "one first, e.g. `dplyr::select(x, ", mv[1], ")`."
    )
  }

  kd       <- tsibble::key_data(x)
  key_cols <- setdiff(names(kd), ".rows")
  has_key  <- length(key_cols) > 0L

  series <- vector("list", nrow(kd))
  for (i in seq_len(nrow(kd))) {
    rows <- kd[[".rows"]][[i]]
    sub  <- x[rows, ]
    key_val <- if (has_key) {
      tibble::as_tibble(kd[i, key_cols, drop = FALSE])
    } else {
      NULL
    }
    series[[i]] <- list(
      key      = key_val,
      x        = as.numeric(sub[[mv]]),
      N        = N,
      template = sub
    )
  }

  list(series = series, input_class = "tbl_ts")
}

#' Rebuild a numeric vector into a template's class
#'
#' The inverse half of `fs_extract()`'s work: given a plain numeric vector
#' and a `template` object (as attached to an `fs_extract()` series record),
#' wraps `values` back into the template's class -- a `ts` with the
#' template's `tsp` (start, end, frequency), or a `tbl_ts` with the
#' template's single measured column replaced. With `template = NULL` (the
#' bare-numeric case), `values` is returned as-is (coerced to `numeric`).
#'
#' @param values Numeric vector to wrap; its length must match `template`
#'   (`length(template)` for `ts`, `nrow(template)` for `tbl_ts`).
#' @param template `NULL`, a `ts` object, or a `tbl_ts` with exactly one
#'   measured variable, as produced by `fs_extract()`.
#'
#' @return `values`, either unchanged (plain numeric, when `template` is
#'   `NULL`), as a `ts` sharing `template`'s `tsp`, or as a `tbl_ts` sharing
#'   `template`'s index/keys with its measured column replaced by `values`.
#'
#' @keywords internal
#' @noRd
fs_rebuild <- function(values, template) {
  values <- as.numeric(values)

  if (is.null(template)) {
    return(values)
  }

  if (stats::is.ts(template)) {
    if (length(values) != length(template)) {
      stop(
        "`values` length (", length(values), ") does not match `template` ",
        "length (", length(template), ")."
      )
    }
    out <- template
    out[] <- values
    return(out)
  }

  if (inherits(template, "tbl_ts")) {
    mv <- tsibble::measured_vars(template)
    if (length(mv) != 1L) {
      stop("`template` must have exactly one measured variable.")
    }
    if (length(values) != nrow(template)) {
      stop(
        "`values` length (", length(values), ") does not match `template` ",
        "row count (", nrow(template), ")."
      )
    }
    template[[mv[1]]] <- values
    return(template)
  }

  values
}
