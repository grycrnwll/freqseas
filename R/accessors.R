# accessors.R -- exported generics for pulling pieces out of `seas_test`,
# `seas_sa`, and (via R/collection.R) `seas_collection` objects.
#
# This file owns the four accessor GENERICS (`adjusted()`, `seasonal()`,
# `whitener()`, `decision()`) plus their `seas_test`/`seas_sa` methods. The
# `seas_collection` methods live in R/collection.R, next to the class they
# dispatch on; both files' methods share these generics. This file does NOT
# define `seas_test()`/`seas_ssi()`/`seas_adjust()` -- those generics belong
# to the sibling modules (R/seas_test.R, R/seas_ssi.R, R/seas_adjust.R).

#' Extract the adjusted series
#'
#' The seasonally-adjusted series from a fitted object: for a `seas_sa`, its
#' `adjusted` field (the same class as the original input); for a
#' `seas_collection` of adjustments, a tibble of key columns plus a list
#' column of per-key adjusted series (see the `seas_collection` method in
#' `collection.R`).
#'
#' @param object A fitted freqseas object (`seas_sa`, or a `seas_collection`
#'   produced by [seas_adjust()]/`seas_ssi()`).
#' @param ... Passed on to methods; unused by the methods defined here.
#'
#' @return The adjusted series, in the same class as the original input for
#'   a single-series object, or a tibble (key columns + list-column) for a
#'   collection.
#'
#' @examples
#' fit <- list(
#'   test = list(x = 1:12, input_class = "numeric", alpha = 0.05,
#'               evt = list(p = 0.01), spec = list(label = "line")),
#'   adjusted = 1:12 - 0.1, seasonal = rep(0.1, 12),
#'   phase_rule = "minimum", post_evt = 0.6
#' )
#' class(fit) <- "seas_sa"
#' adjusted(fit)
#'
#' @export
adjusted <- function(object, ...) UseMethod("adjusted")

#' @rdname adjusted
#' @export
adjusted.seas_sa <- function(object, ...) object$adjusted

#' @rdname adjusted
#' @export
adjusted.seas_test <- function(object, ...) {
  stop(
    "A `seas_test` has no adjusted series (it only records the detection/",
    "specification stage). Run `seas_ssi()` or `seas_adjust()` first."
  )
}

#' @rdname adjusted
#' @export
adjusted.default <- function(object, ...) {
  stop(
    "`adjusted()` is not defined for an object of class `",
    paste(class(object), collapse = "/"),
    "`. Use it on a `seas_sa` or a `seas_collection` of adjustments."
  )
}

#' Extract the removed seasonal component
#'
#' The component removed by the adjustment: for a `seas_sa`, its `seasonal`
#' field (`x - adjusted`, same class as the original input); for a
#' `seas_collection` of adjustments, a tibble of key columns plus a list
#' column of per-key seasonal components.
#'
#' @param object A fitted freqseas object (`seas_sa`, or a `seas_collection`
#'   produced by [seas_adjust()]/`seas_ssi()`).
#' @param ... Passed on to methods; unused by the methods defined here.
#'
#' @return The seasonal component, in the same class as the original input
#'   for a single-series object, or a tibble (key columns + list-column) for
#'   a collection.
#'
#' @examples
#' fit <- list(
#'   test = list(x = 1:12, input_class = "numeric", alpha = 0.05,
#'               evt = list(p = 0.01), spec = list(label = "line")),
#'   adjusted = 1:12 - 0.1, seasonal = rep(0.1, 12),
#'   phase_rule = "minimum", post_evt = 0.6
#' )
#' class(fit) <- "seas_sa"
#' seasonal(fit)
#'
#' @export
seasonal <- function(object, ...) UseMethod("seasonal")

#' @rdname seasonal
#' @export
seasonal.seas_sa <- function(object, ...) object$seasonal

#' @rdname seasonal
#' @export
seasonal.seas_test <- function(object, ...) {
  stop(
    "A `seas_test` has no seasonal component (it only records the ",
    "detection/specification stage). Run `seas_ssi()` or `seas_adjust()` ",
    "first."
  )
}

#' @rdname seasonal
#' @export
seasonal.default <- function(object, ...) {
  stop(
    "`seasonal()` is not defined for an object of class `",
    paste(class(object), collapse = "/"),
    "`. Use it on a `seas_sa` or a `seas_collection` of adjustments."
  )
}

#' Extract the whitener record
#'
#' The M0-Whittle whitener record (differencing decision, AR order/
#' coefficients, BIC path) attached to a fitted object: for a `seas_test`,
#' its `whitener` field directly; for a `seas_sa`, its embedded test's
#' `whitener` field (`object$test$whitener`); for a `seas_collection`, a
#' tibble of key columns plus a list column of per-key whitener records.
#'
#' @param object A fitted freqseas object (`seas_test`, `seas_sa`, or a
#'   `seas_collection`).
#' @param ... Passed on to methods; unused by the methods defined here.
#'
#' @return The whitener record: a list with (at least) `d`, `d_reason`, `p`,
#'   `ar`, and `bic_path`, for a single-series object, or a tibble (key
#'   columns + list-column) for a collection.
#'
#' @examples
#' tst <- list(
#'   x = 1:12, input_class = "numeric", N = 12, P = 1, M = 3, alpha = 0.05,
#'   evt = list(statistic = 2.1, p = 0.01, critical = 1.8),
#'   spec = list(label = "line", shoulder_p = 0.2, phase_R = 0.9),
#'   whitener = list(d = 0, d_reason = "not near-integrated", p = 1,
#'                    ar = 0.4, bic_path = c(-10, -12, -11))
#' )
#' class(tst) <- "seas_test"
#' whitener(tst)
#'
#' @export
whitener <- function(object, ...) UseMethod("whitener")

#' @rdname whitener
#' @export
whitener.seas_test <- function(object, ...) object$whitener

#' @rdname whitener
#' @export
whitener.seas_sa <- function(object, ...) object$test$whitener

#' @rdname whitener
#' @export
whitener.default <- function(object, ...) {
  stop(
    "`whitener()` is not defined for an object of class `",
    paste(class(object), collapse = "/"),
    "`. Use it on a `seas_test`, a `seas_sa`, or a `seas_collection`."
  )
}

#' Extract the seasonality decision
#'
#' A one-row (single object) or one-row-per-key (`seas_collection`) tibble
#' summarizing the detection decision: whether the series is judged
#' seasonal (`evt$p < alpha`), the detection p-value, the specification
#' label (`"line"`/`"band"`), and the test level `alpha`.
#'
#' @param object A fitted freqseas object (`seas_test`, `seas_sa`, or a
#'   `seas_collection`).
#' @param ... Passed on to methods; unused by the methods defined here.
#'
#' @return A [tibble::tibble()] with columns `seasonal` (logical), `p`
#'   (numeric), `spec` (character, `"line"`/`"band"`), and `alpha`
#'   (numeric); one row for a single object, or one row per key (plus the
#'   key columns) for a `seas_collection`.
#'
#' @examples
#' tst <- list(
#'   x = 1:12, input_class = "numeric", N = 12, P = 1, M = 3, alpha = 0.05,
#'   evt = list(statistic = 2.1, p = 0.01, critical = 1.8),
#'   spec = list(label = "line", shoulder_p = 0.2, phase_R = 0.9),
#'   whitener = list(d = 0, d_reason = "not near-integrated", p = 1,
#'                    ar = 0.4, bic_path = c(-10, -12, -11))
#' )
#' class(tst) <- "seas_test"
#' decision(tst)
#'
#' @export
decision <- function(object, ...) UseMethod("decision")

#' @rdname decision
#' @export
decision.seas_test <- function(object, ...) {
  # Prefer a stored `decision` field when the object carries one (the actual
  # `seas_test` constructor does: `decision = evt$p < alpha` at fit time, and
  # future logic -- e.g. persisted panel decisions re-evaluated only at
  # annual reviews -- may make that field diverge from a naive recompute).
  # Fall back to the documented rule for minimal hand-built fixtures that
  # don't carry a `decision` field at all.
  seasonal <- if (!is.null(object$decision)) {
    isTRUE(object$decision)
  } else {
    object$evt$p < object$alpha
  }
  tibble::tibble(
    seasonal = seasonal,
    p        = object$evt$p,
    spec     = object$spec$label,
    alpha    = object$alpha
  )
}

#' @rdname decision
#' @export
decision.seas_sa <- function(object, ...) decision(object$test)

#' @rdname decision
#' @export
decision.default <- function(object, ...) {
  stop(
    "`decision()` is not defined for an object of class `",
    paste(class(object), collapse = "/"),
    "`. Use it on a `seas_test`, a `seas_sa`, or a `seas_collection`."
  )
}
