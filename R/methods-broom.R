# methods-broom.R
# broom-style tidiers for `seas_test`/`seas_sa`/`seas_collection`, built on
# the `generics` package (an Imports dependency, unlike ggplot2 -- see
# R/methods-plot.R for the Suggests-only autoplot() registration). The
# `generics` package supplies the bare `tidy()`/`glance()`/`augment()`
# generics that `broom` and its extension packages share; this file both
# defines the S3 methods and re-exports the three bare generic names (the
# standard "broom-extension" pattern: `#' @importFrom generics tidy` +
# `#' @export` on an assignment of the generic itself, so `library(freqseas)`
# alone makes bare `tidy()`/`glance()`/`augment()` calls dispatch correctly).
#
# Each individual method (`tidy.seas_test()`, `glance.seas_sa()`, ...) is
# tagged `@exportS3Method generics::tidy` (etc.), NOT a plain `@export`,
# mirroring the `@exportS3Method ggplot2::autoplot` pattern already used in
# R/methods-plot.R. This is load-bearing, not stylistic: `tidy <- generics::tidy`
# makes the LOCAL name resolve to the generic, but the generic's own closure
# environment is `namespace:generics`, and S3 dispatch consults the method
# table attached to THAT environment. A bare `@export` on `tidy.seas_test`
# lets roxygen2 infer an ambiguous `S3method(tidy, seas_test)` NAMESPACE
# directive that (verified empirically, both under `pkgload::load_all()` and
# a real `R CMD INSTALL`) registers the method into freqseas's OWN S3 table
# instead of `generics`'s -- so bare `tidy()`/`glance()`/`augment()` calls on
# a `seas_sa`/`seas_test`/`seas_collection` raise "no applicable method",
# even though `tidy.seas_test()` etc. work fine when called directly. The
# qualified `@exportS3Method generics::tidy` tag emits the correctly-qualified
# `S3method(generics::tidy, seas_test)` directive, which registers into the
# right table. Fixed 2026-07-22 while verifying the documentation vignettes
# (getting-started.Rmd exercises glance()/augment()); tags-only, no method
# body changed.
#
# `tidy.seas_test()` renders `object$spec$harmonic_table` (one row per row of
# that table) when present -- a per-harmonic diagnostic field a
# concurrently-developed sibling module may add to `spec`; it does not exist
# yet on any `seas_test` produced by the current R/seas_test.R, so the
# `is.null()` branch below is the one currently exercised.

#' Tidiers re-exported from the generics package
#'
#' freqseas re-exports the bare [generics::tidy()], [generics::glance()], and
#' [generics::augment()] generics, so that attaching the package
#' (`library(freqseas)`) is enough for bare `tidy()` / `glance()` / `augment()`
#' calls to dispatch to the freqseas methods. Each method is documented on its
#' own page (for example [tidy.seas_test()], [glance.seas_sa()], and
#' [augment.seas_sa()]).
#'
#' @name freqseas-reexports
#' @keywords internal
#' @usage NULL
#' @importFrom generics tidy
#' @export
tidy <- generics::tidy

#' @rdname freqseas-reexports
#' @usage NULL
#' @importFrom generics glance
#' @export
glance <- generics::glance

#' @rdname freqseas-reexports
#' @usage NULL
#' @importFrom generics augment
#' @export
augment <- generics::augment

# ---------------------------------------------------------------------------
# tidy()
# ---------------------------------------------------------------------------

#' Tidy a `seas_test` object
#'
#' One row per harmonic when `object$spec$harmonic_table` is present
#' (rendered as-is via [tibble::as_tibble()]); otherwise a single summary row
#' with the detection statistic, p-value, specification label, and bin count.
#'
#' @param x A `seas_test` object.
#' @param ... Ignored.
#'
#' @return A [tibble::tibble()]. Without a harmonic table: one row with
#'   columns `statistic`, `p.value`, `spec`, `M`. With a harmonic table: one
#'   row per harmonic, columns as provided by `spec$harmonic_table`.
#'
#' @examples
#' set.seed(2026)
#' n <- 120L; Phi <- 0.9
#' e <- stats::rnorm(n + 400)
#' z <- numeric(n + 400)
#' for (t in 6:(n + 400)) z[t] <- 0.5 * z[t - 1] + Phi * z[t - 4] -
#'   0.5 * Phi * z[t - 5] + e[t]
#' tst <- seas_test(z[401:(n + 400)], frequency = 4)
#' tidy(tst)
#'
#' @exportS3Method generics::tidy
tidy.seas_test <- function(x, ...) {
  if (!is.null(x$spec$harmonic_table)) {
    return(tibble::as_tibble(x$spec$harmonic_table))
  }
  tibble::tibble(
    statistic = x$evt$statistic,
    p.value   = x$evt$p,
    spec      = x$spec$label,
    M         = x$M
  )
}

#' Tidy a `seas_sa` object
#'
#' Delegates to [tidy.seas_test()] on the embedded test (`x$test`): the
#' adjustment stage does not add per-harmonic structure of its own.
#'
#' @param x A `seas_sa` object.
#' @param ... Passed on to [tidy.seas_test()].
#'
#' @return See [tidy.seas_test()].
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
#' tidy(adj)
#'
#' @exportS3Method generics::tidy
tidy.seas_sa <- function(x, ...) {
  tidy.seas_test(x$test, ...)
}

#' Tidy a `seas_collection` object
#'
#' Row-binds [tidy.seas_test()]/[tidy.seas_sa()] (as appropriate to each
#' row's fit) across every key, prefixing the key columns onto each row
#' (replicated across multiple rows per key when `tidy()` on a single fit
#' returns more than one row, e.g. a harmonic table).
#'
#' @param x A `seas_collection`.
#' @param ... Passed on to the per-key tidy method.
#'
#' @return A [tibble::tibble()]: the key columns followed by the columns
#'   from [tidy.seas_test()]/[tidy.seas_sa()], one or more rows per key.
#'
#' @examples
#' \dontrun{
#' tidy(seas_test(some_keyed_tsibble))
#' }
#'
#' @exportS3Method generics::tidy
tidy.seas_collection <- function(x, ...) {
  key_cols <- setdiff(names(x), "fit")
  rows <- lapply(seq_len(nrow(x)), function(i) {
    fit <- x$fit[[i]]
    td  <- if (inherits(fit, "seas_sa")) tidy.seas_sa(fit, ...) else tidy.seas_test(fit, ...)
    key_row <- tibble::as_tibble(x[i, key_cols, drop = FALSE])
    dplyr::bind_cols(key_row[rep(1L, nrow(td)), , drop = FALSE], td)
  })
  dplyr::bind_rows(rows)
}

# ---------------------------------------------------------------------------
# glance()
# ---------------------------------------------------------------------------

#' Glance at a `seas_test` object
#'
#' A one-row model-level summary: the detection decision, p-value,
#' specification label, bin count, differencing order, and AR order.
#'
#' @param x A `seas_test` object.
#' @param ... Ignored.
#'
#' @return A one-row [tibble::tibble()] with columns `decision` (logical),
#'   `p.value`, `spec`, `M`, `d`, `ar_order`.
#'
#' @examples
#' set.seed(2026)
#' n <- 120L; Phi <- 0.9
#' e <- stats::rnorm(n + 400)
#' z <- numeric(n + 400)
#' for (t in 6:(n + 400)) z[t] <- 0.5 * z[t - 1] + Phi * z[t - 4] -
#'   0.5 * Phi * z[t - 5] + e[t]
#' tst <- seas_test(z[401:(n + 400)], frequency = 4)
#' glance(tst)
#'
#' @exportS3Method generics::glance
glance.seas_test <- function(x, ...) {
  tibble::tibble(
    decision = isTRUE(x$decision),
    p.value  = x$evt$p,
    spec     = x$spec$label,
    M        = x$M,
    d        = x$whitener$d,
    ar_order = x$whitener$p
  )
}

#' Glance at a `seas_sa` object
#'
#' A one-row model-level summary: everything [glance.seas_test()] reports for
#' the embedded test, plus the declared phase rule, the post-adjustment
#' detection p-value (the self-check), and the maximum gain applied.
#'
#' @param x A `seas_sa` object.
#' @param ... Ignored.
#'
#' @return A one-row [tibble::tibble()] with columns `decision`, `p.value`,
#'   `spec`, `M`, `d`, `ar_order` (as [glance.seas_test()]), plus
#'   `phase_rule`, `post_p`, `max_gain`.
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
#' glance(adj)
#'
#' @exportS3Method generics::glance
glance.seas_sa <- function(x, ...) {
  base <- glance.seas_test(x$test)
  Ghat <- x$Ghat
  max_gain <- if (length(Ghat) > 1L) max(Ghat) else as.numeric(Ghat)
  dplyr::bind_cols(base, tibble::tibble(
    phase_rule = x$phase_rule,
    post_p     = x$post_evt,
    max_gain   = max_gain
  ))
}

#' Glance at a `seas_collection` object
#'
#' Row-binds [glance.seas_test()]/[glance.seas_sa()] (as appropriate to each
#' row's fit) across every key, prefixing the key columns onto each row.
#'
#' @param x A `seas_collection`.
#' @param ... Ignored.
#'
#' @return A one-row-per-key [tibble::tibble()]: the key columns followed by
#'   the columns from [glance.seas_test()]/[glance.seas_sa()].
#'
#' @examples
#' \dontrun{
#' glance(seas_test(some_keyed_tsibble))
#' }
#'
#' @exportS3Method generics::glance
glance.seas_collection <- function(x, ...) {
  key_cols <- setdiff(names(x), "fit")
  rows <- lapply(x$fit, function(fit) {
    if (inherits(fit, "seas_sa")) glance.seas_sa(fit) else glance.seas_test(fit)
  })
  dplyr::bind_cols(tibble::as_tibble(x[key_cols]), dplyr::bind_rows(rows))
}

# ---------------------------------------------------------------------------
# augment()
# ---------------------------------------------------------------------------

#' Augment a `seas_sa` object
#'
#' Observation-level output: the original series, the adjusted series, and
#' the removed seasonal component, one row per observation. For a `ts`- or
#' numeric-input adjustment, the index column is the (fractional-year, for
#' `ts`) time value; for a `tbl_ts`-input adjustment, the result is itself a
#' tsibble indexed on the original tsibble's index column.
#'
#' @param x A `seas_sa` object.
#' @param ... Ignored.
#'
#' @return A [tibble::tibble()] (or [tsibble::tsibble()] for `tbl_ts` input)
#'   with columns `index` (or the original tsibble's index column name), `x`
#'   (the original series), `adjusted`, `seasonal`.
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
#' augment(adj)
#'
#' @exportS3Method generics::augment
augment.seas_sa <- function(x, ...) {
  tst  <- x$test
  # `adjusted`/`seasonal` are input-class-preserving: plain numeric, `ts`, or
  # -- for a tbl_ts-input fit produced through the collection layer -- a full
  # `tbl_ts` re-clothed by `fs_reclothe_sa()` (R/collection.R). `as.numeric()`
  # errors on that last case (a multi-column tsibble is a list, not a plain
  # vector), so pull the measured column out first when it applies.
  as_num <- function(v) {
    if (inherits(v, "tbl_ts")) {
      mv <- tsibble::measured_vars(v)
      return(as.numeric(v[[mv[1L]]]))
    }
    as.numeric(v)
  }
  orig <- as_num(tst$x)
  adj  <- as_num(x$adjusted)
  seas <- as_num(x$seasonal)

  if (identical(tst$input_class, "tbl_ts") && !is.null(tst$template)) {
    idx_col  <- tsibble::index_var(tst$template)
    idx_vals <- tst$template[[idx_col]]
    out <- tibble::tibble(!!idx_col := idx_vals, x = orig, adjusted = adj,
                          seasonal = seas)
    return(tsibble::as_tsibble(out, index = !!rlang::sym(idx_col)))
  }

  idx <- if (identical(tst$input_class, "ts") && !is.null(tst$tsp)) {
    as.numeric(stats::time(stats::ts(orig, start = tst$tsp[1L],
                                     frequency = tst$tsp[3L])))
  } else {
    seq_along(orig)
  }
  tibble::tibble(index = idx, x = orig, adjusted = adj, seasonal = seas)
}

#' Augment a `seas_collection` object
#'
#' Row-binds [augment.seas_sa()] across every key, prefixing the key columns
#' onto each row. Requires a `type = "sa"` collection (the result of
#' `seas_ssi()`/`seas_adjust()`); a `type = "test"` collection has no
#' adjusted/seasonal series to augment.
#'
#' @param x A `seas_collection` of type `"sa"`.
#' @param ... Ignored.
#'
#' @return A [tibble::tibble()]: the key columns followed by the columns
#'   from [augment.seas_sa()], one row per observation per key.
#'
#' @examples
#' \dontrun{
#' augment(seas_ssi(seas_test(some_keyed_tsibble), phase_rule = "minimum"))
#' }
#'
#' @exportS3Method generics::augment
augment.seas_collection <- function(x, ...) {
  type <- attr(x, "type")
  if (!identical(type, "sa")) {
    stop(
      "`augment()` requires a `seas_collection` of adjustments ",
      "(type = \"sa\"); run `seas_ssi()`/`seas_adjust()` first.",
      call. = FALSE
    )
  }
  key_cols <- setdiff(names(x), "fit")
  rows <- lapply(seq_len(nrow(x)), function(i) {
    fit <- x$fit[[i]]
    # Coerced to a plain tibble (rather than left as a tsibble on tbl_ts
    # input): row-binding several keys' per-key tsibbles together would
    # violate tsibble's own index/key-uniqueness assumptions (the same time
    # index recurs once per key with no key column of its own). The
    # documented, and consistent-with-collection.R-precedent, return shape
    # for every seas_collection tidier is a flat keyed tibble.
    au <- tibble::as_tibble(augment.seas_sa(fit))
    key_row <- tibble::as_tibble(x[i, key_cols, drop = FALSE])
    dplyr::bind_cols(key_row[rep(1L, nrow(au)), , drop = FALSE], au)
  })
  dplyr::bind_rows(rows)
}
