# Tests for the methods/tidier layer: summary()/print() (R/methods-summary.R),
# plot()/autoplot() (R/methods-plot.R), and the broom generics
# tidy()/glance()/augment() (R/methods-broom.R), across seas_test, seas_sa,
# and seas_collection objects.

# --- assembly-time integration guard (see test-seas-ssi.R for the rationale):
# restore the real single-series + collection layers in case the
# concurrently-built test-collection.R clobbered the globals with `<<-`
# stubs when the full suite runs alphabetically before this file. No-op
# under a real package build.
local({
  r_dir <- file.path("..", "..", "R")
  if (file.exists(file.path(r_dir, "seas_test.R"))) {
    for (f in c("classes.R", "seas_test.R", "seas_ssi.R", "seas_adjust.R",
                "methods-print.R", "collection.R")) {
      sys.source(file.path(r_dir, f), envir = globalenv())
    }
  }
})

# Filter-world quarterly DGP: (1 - phi B)(1 - Phi B^s) x = eps -- a band
# spectrum with seasonal peaks at pi/2 and pi. Self-contained (see
# test-seas-test.R).
sim_filter_q <- function(n, Phi = 0.8, phi = 0.5, s = 4L, spin = 400L) {
  ntot <- spin + n
  eps  <- stats::rnorm(ntot)
  z    <- numeric(ntot)
  for (t in (s + 2L):ntot) {
    z[t] <- phi * z[t - 1L] + Phi * z[t - s] - phi * Phi * z[t - s - 1L] + eps[t]
  }
  out <- z[(spin + 1L):ntot]
  out - mean(out)
}

make_keyed_q <- function(n = 80L, seeds = c(101L, 102L, 103L)) {
  grps <- c("a", "b", "c")
  vals <- unlist(lapply(seeds, function(sd) {
    set.seed(sd)
    sim_filter_q(n, Phi = 0.85)
  }))
  tsibble::tsibble(
    qtr   = rep(tsibble::yearquarter("2000 Q1") + 0:(n - 1L), length(grps)),
    grp   = rep(grps, each = n),
    value = vals,
    key   = grp,
    index = qtr
  )
}

make_many_keys_q <- function(k = 13L, n = 24L) {
  keys <- sprintf("s%02d", seq_len(k))
  vals <- unlist(lapply(seq_len(k), function(i) {
    set.seed(2000L + i)
    stats::rnorm(n)
  }))
  tsibble::tsibble(
    qtr   = rep(tsibble::yearquarter("2000 Q1") + 0:(n - 1L), k),
    grp   = rep(keys, each = n),
    value = vals,
    key   = grp,
    index = qtr
  )
}

# ---- shared fixtures -------------------------------------------------------

set.seed(2026)
tst <- seas_test(sim_filter_q(120, Phi = 0.9), frequency = 4)
sa  <- seas_ssi(tst, phase_rule = "minimum")

coll_test <- seas_test(make_keyed_q())
coll_sa   <- seas_ssi(coll_test, phase_rule = "minimum")

# =============================================================================
# summary()
# =============================================================================

test_that("summary.seas_test has the expected class and fields", {
  s <- summary(tst)
  expect_s3_class(s, "summary.seas_test")
  expect_true(is.logical(s$decision))
  expect_true(is.numeric(s$p))
  expect_true(all(c("statistic", "p", "critical", "N1", "N0") %in% names(s$evt)))
  expect_true(all(c("label", "shoulder_p", "phase_R") %in% names(s$spec)))
  expect_true(all(c("d", "d_reason", "p", "ar", "bic_path") %in% names(s$whitener)))
  expect_equal(s$M, tst$M)
  # wave-3: seas_test() now produces spec$harmonic_table, and summary carries it.
  expect_s3_class(s$harmonic_table, "tbl_df")
  expect_true(all(c("omega", "label", "excess", "elevated") %in%
                    names(s$harmonic_table)))
  expect_type(s$donor_pool_size, "integer")
})

test_that("print.summary.seas_test runs without error and is more detailed than print()", {
  out <- capture.output(print(summary(tst)))
  expect_true(any(grepl("Detection", out)))
  expect_true(any(grepl("Specification", out)))
  expect_true(any(grepl("Bins", out)))
  expect_true(any(grepl("Whitener", out)))
  expect_true(any(grepl("Donor pool size", out)))
  expect_gt(length(out), length(capture.output(print(tst))))
})

test_that("summary.seas_sa has the expected class and fields", {
  s <- summary(sa)
  expect_s3_class(s, "summary.seas_sa")
  expect_s3_class(s$test_summary, "summary.seas_test")
  expect_equal(s$phase_rule, "minimum")
  expect_true(all(c("max_Ghat", "ordinates_touched", "n_ordinates") %in% names(s$gain)))
  expect_equal(s$post_evt, sa$post_evt)
})

test_that("print.summary.seas_sa runs without error and reports gain + post-adjustment p", {
  out <- capture.output(print(summary(sa)))
  expect_true(any(grepl("Gain surgery", out)))
  expect_true(any(grepl("Post-adjustment detection p", out)))
  expect_true(any(grepl("declared identification", out)))
  expect_true(any(grepl("embedded seas_test summary", out)))
})

test_that("summary.seas_collection maps over keys and prints without error", {
  s <- summary(coll_test)
  expect_s3_class(s, "summary.seas_collection")
  expect_equal(s$type, "test")
  expect_equal(nrow(s$keys), 3L)
  expect_length(s$summaries, 3L)
  expect_true(all(vapply(s$summaries, inherits, logical(1), what = "summary.seas_test")))

  out <- capture.output(print(s))
  expect_true(grepl("3 series", paste(out, collapse = "\n")))

  s_sa <- summary(coll_sa)
  expect_equal(s_sa$type, "sa")
  expect_true(all(vapply(s_sa$summaries, inherits, logical(1), what = "summary.seas_sa")))
  out_sa <- capture.output(print(s_sa))
  expect_true(any(grepl("Gain surgery", out_sa)))
})

# =============================================================================
# plot()
# =============================================================================

test_that("plot.seas_test runs into a png device without error", {
  f <- tempfile(fileext = ".png")
  grDevices::png(f)
  expect_no_error(plot(tst))
  grDevices::dev.off()
  expect_true(file.exists(f))
})

test_that("plot.seas_sa runs both types into a png device without error", {
  f <- tempfile(fileext = ".png")
  grDevices::png(f)
  expect_no_error(plot(sa))
  expect_no_error(plot(sa, type = "periodogram"))
  grDevices::dev.off()
  expect_true(file.exists(f))
})

test_that("plot.seas_collection runs for a test-type and a sa-type 3-key collection", {
  f <- tempfile(fileext = ".png")
  grDevices::png(f)
  expect_no_error(plot(coll_test))
  expect_no_error(plot(coll_sa))
  grDevices::dev.off()
  expect_true(file.exists(f))
})

test_that("plot.seas_collection warns and subsets beyond max_keys", {
  many <- make_many_keys_q(13L, 24L)
  coll <- seas_test(many)
  expect_equal(nrow(coll), 13L)

  f <- tempfile(fileext = ".png")
  grDevices::png(f)
  expect_warning(plot(coll), "12")
  grDevices::dev.off()
  expect_true(file.exists(f))
})

test_that("plot.seas_collection keeps the \"sa\" series branch past the max_keys subset", {
  # Regression test: `type` (and `key_cols`) must be read from the FULL `x`
  # BEFORE the >max_keys subsetting, because tibble row-indexing (`x[i, ]`)
  # drops the `seas_collection` class and its `type` attribute. Reading
  # `type` after subsetting would silently fall through to the periodogram
  # branch for a "sa" collection, which reaches into fields (`$partition`,
  # `$pgram`) a `seas_sa` fit doesn't have and errors out on a non-finite
  # `ylim` inside `plot.default()`.
  many    <- make_many_keys_q(13L, 24L)
  many_sa <- suppressMessages(seas_ssi(seas_test(many), phase_rule = "minimum"))
  expect_equal(attr(many_sa, "type"), "sa")

  f <- tempfile(fileext = ".png")
  grDevices::png(f)
  expect_warning(plot(many_sa), "12")
  grDevices::dev.off()
  expect_true(file.exists(f))
})

# =============================================================================
# autoplot() (ggplot2, Suggests)
# =============================================================================

test_that("autoplot methods return ggplot objects that actually build", {
  testthat::skip_if_not_installed("ggplot2")

  # expect_s3_class() alone only checks that ggplot() was constructed; the
  # aes()/facet_wrap() expressions (in particular the `.data$panel`/
  # `.data$key` facet formulas) are lazily evaluated and are only forced by
  # ggplot_build(). Build every one so a broken facet/aes/scale surfaces here.
  p1 <- ggplot2::autoplot(tst)
  p2 <- ggplot2::autoplot(sa)
  p3 <- ggplot2::autoplot(sa, type = "periodogram")
  p4 <- ggplot2::autoplot(coll_test)
  p5 <- ggplot2::autoplot(coll_sa)

  for (p in list(p1, p2, p3, p4, p5)) {
    expect_s3_class(p, "ggplot")
    expect_no_error(ggplot2::ggplot_build(p))
  }
})

test_that("plot.seas_collection/autoplot.seas_collection keep the \"sa\" series branch past max_keys (ggplot2)", {
  testthat::skip_if_not_installed("ggplot2")
  many    <- make_many_keys_q(13L, 24L)
  many_sa <- suppressMessages(seas_ssi(seas_test(many), phase_rule = "minimum"))
  expect_warning(p <- ggplot2::autoplot(many_sa), "12")
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("autoplot.seas_collection warns and subsets beyond max_keys", {
  testthat::skip_if_not_installed("ggplot2")
  many <- make_many_keys_q(13L, 24L)
  coll <- seas_test(many)
  expect_warning(p <- ggplot2::autoplot(coll), "12")
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})

# =============================================================================
# tidy() / glance() / augment()  (generics package)
# =============================================================================

test_that("tidy.seas_test renders the per-harmonic table (wave-3)", {
  td <- generics::tidy(tst)
  expect_s3_class(td, "tbl_df")
  # seas_test() now attaches spec$harmonic_table, so tidy() switches to one row
  # per harmonic (quarterly P=1 -> pi/2 and pi) instead of the one-row summary.
  expect_equal(nrow(td), nrow(tst$spec$harmonic_table))
  expect_true(all(c("omega", "label", "excess", "elevated") %in% names(td)))
})

test_that("tidy.seas_sa delegates to the embedded test", {
  td <- generics::tidy(sa)
  expect_equal(td, generics::tidy(tst))
})

test_that("glance.seas_test returns the documented one-row shape", {
  gl <- generics::glance(tst)
  expect_s3_class(gl, "tbl_df")
  expect_equal(nrow(gl), 1L)
  expect_true(all(c("decision", "p.value", "spec", "M", "d", "ar_order") %in% names(gl)))
})

test_that("glance.seas_sa adds phase_rule/post_p/max_gain to the test glance", {
  gl <- generics::glance(sa)
  expect_true(all(c("decision", "p.value", "spec", "M", "d", "ar_order",
                     "phase_rule", "post_p", "max_gain") %in% names(gl)))
  expect_equal(gl$phase_rule, "minimum")
  expect_equal(gl$post_p, sa$post_evt)
})

test_that("augment.seas_sa returns index/x/adjusted/seasonal with the right length", {
  au <- generics::augment(sa)
  expect_true(all(c("index", "x", "adjusted", "seasonal") %in% names(au)))
  expect_equal(nrow(au), length(sa$test$x))
  expect_equal(au$x, as.numeric(sa$test$x))
  expect_equal(au$adjusted, as.numeric(sa$adjusted))
  expect_equal(au$seasonal, as.numeric(sa$seasonal))
})

test_that("collection tidiers row-bind per key and carry key columns", {
  td_c <- generics::tidy(coll_test)
  expect_true("grp" %in% names(td_c))
  # wave-3: each key's tidy() is now the per-harmonic table, so the collection
  # tidy has (n harmonics) rows per key.
  per_key <- nrow(coll_test$fit[[1]]$spec$harmonic_table)
  expect_equal(nrow(td_c), 3L * per_key)
  expect_setequal(td_c$grp, c("a", "b", "c"))

  gl_c <- generics::glance(coll_test)
  expect_true(all(c("grp", "decision", "p.value", "spec", "M") %in% names(gl_c)))
  expect_equal(nrow(gl_c), 3L)

  gl_c_sa <- generics::glance(coll_sa)
  expect_true(all(c("grp", "phase_rule", "post_p", "max_gain") %in% names(gl_c_sa)))
  expect_equal(nrow(gl_c_sa), 3L)

  au_c <- generics::augment(coll_sa)
  # tbl_ts-input augment() keeps the ORIGINAL index column name ("qtr", from
  # make_keyed_q()'s fixture) rather than a generic "index" -- see
  # augment.seas_sa()'s tbl_ts branch.
  expect_true(all(c("grp", "qtr", "x", "adjusted", "seasonal") %in% names(au_c)))
  expect_equal(nrow(au_c), 3L * 80L)
})

test_that("augment.seas_collection errors informatively for a test-type collection", {
  expect_error(generics::augment(coll_test), "seas_ssi|seas_adjust")
})
