# Tests for R/coerce.R: fs_extract() / fs_rebuild(), the input-normalization
# and reconstruction pair behind the tbl_ts collection layer.

test_that("fs_extract: bare numeric requires and uses `frequency`", {
  ex <- fs_extract(1:24, frequency = 12)
  expect_equal(ex$input_class, "numeric")
  expect_length(ex$series, 1L)
  expect_null(ex$series[[1]]$key)
  expect_equal(ex$series[[1]]$x, as.numeric(1:24))
  expect_equal(ex$series[[1]]$N, 12)
  expect_null(ex$series[[1]]$template)
})

test_that("fs_extract: bare numeric without `frequency` errors", {
  expect_error(fs_extract(1:24), "frequency")
  expect_error(fs_extract(1:24, frequency = NULL), "frequency")
})

test_that("fs_extract: bare numeric rejects a bad `frequency`", {
  expect_error(fs_extract(1:24, frequency = -1), "positive")
  expect_error(fs_extract(1:24, frequency = c(4, 12)), "single")
})

test_that("fs_extract/fs_rebuild: ts roundtrip preserves tsp exactly", {
  x <- stats::ts(stats::rnorm(37), start = c(2001, 3), frequency = 12)
  ex <- fs_extract(x)

  expect_equal(ex$input_class, "ts")
  expect_length(ex$series, 1L)
  expect_null(ex$series[[1]]$key)
  expect_equal(ex$series[[1]]$N, 12)
  expect_equal(ex$series[[1]]$x, as.numeric(x))
  expect_identical(ex$series[[1]]$template, x)

  new_vals <- as.numeric(x) + 1
  rebuilt <- fs_rebuild(new_vals, ex$series[[1]]$template)

  expect_true(stats::is.ts(rebuilt))
  expect_identical(stats::tsp(rebuilt), stats::tsp(x))
  expect_equal(as.numeric(rebuilt), new_vals)
})

test_that("fs_rebuild: NULL template returns a plain numeric vector", {
  out <- fs_rebuild(c(1, 2, 3), NULL)
  expect_type(out, "double")
  expect_false(inherits(out, "ts"))
  expect_equal(out, c(1, 2, 3))
})

test_that("fs_rebuild: ts length mismatch errors", {
  x <- stats::ts(1:10, frequency = 4)
  expect_error(fs_rebuild(1:5, x), "length")
})

test_that("fs_extract: quarterly tsibble infers N = 4", {
  q <- tsibble::tsibble(
    qtr = tsibble::yearquarter("2020 Q1") + 0:11,
    value = as.numeric(1:12),
    index = qtr
  )
  ex <- fs_extract(q)

  expect_equal(ex$input_class, "tbl_ts")
  expect_length(ex$series, 1L)
  expect_null(ex$series[[1]]$key)
  expect_equal(ex$series[[1]]$N, 4L)
  expect_equal(ex$series[[1]]$x, as.numeric(1:12))
  expect_s3_class(ex$series[[1]]$template, "tbl_ts")
})

test_that("fs_extract: monthly tsibble infers N = 12", {
  m <- tsibble::tsibble(
    mth = tsibble::yearmonth("2020 Jan") + 0:23,
    value = as.numeric(1:24),
    index = mth
  )
  ex <- fs_extract(m)

  expect_equal(ex$input_class, "tbl_ts")
  expect_equal(ex$series[[1]]$N, 12L)
  expect_equal(ex$series[[1]]$x, as.numeric(1:24))
})

test_that("fs_extract/fs_rebuild: tsibble roundtrip replaces the measured column", {
  q <- tsibble::tsibble(
    qtr = tsibble::yearquarter("2020 Q1") + 0:7,
    value = as.numeric(1:8),
    index = qtr
  )
  ex <- fs_extract(q)
  tmpl <- ex$series[[1]]$template

  new_vals <- as.numeric(1:8) * 10
  rebuilt <- fs_rebuild(new_vals, tmpl)

  expect_s3_class(rebuilt, "tbl_ts")
  expect_equal(rebuilt$value, new_vals)
  expect_equal(rebuilt$qtr, q$qtr)
})

test_that("fs_extract: weekly tsibble index errors as unsupported", {
  w <- tsibble::tsibble(
    wk = tsibble::yearweek("2020 W01") + 0:9,
    value = as.numeric(1:10),
    index = wk
  )
  expect_error(fs_extract(w), "unsupported")
})

test_that("fs_extract: Date (or other unrecognized) index errors with guidance", {
  d <- tsibble::tsibble(
    dt = as.Date("2020-01-01") + 0:9,
    value = as.numeric(1:10),
    index = dt
  )
  expect_error(fs_extract(d), "infer")
})

test_that("fs_extract: multi-measure tsibble errors, naming the candidates", {
  q <- tsibble::tsibble(
    qtr = tsibble::yearquarter("2020 Q1") + 0:7,
    a = as.numeric(1:8),
    b = as.numeric(8:1),
    index = qtr
  )
  expect_error(fs_extract(q), "a")
  expect_error(fs_extract(q), "b")
})

test_that("fs_extract: keyed tsibble with 3 keys yields 3 entries in canonical key order", {
  kx <- tsibble::tsibble(
    qtr   = rep(tsibble::yearquarter("2020 Q1") + 0:3, 3),
    grp   = rep(c("c", "a", "b"), each = 4),
    value = as.numeric(seq_len(12)),
    key   = grp,
    index = qtr
  )
  ex <- fs_extract(kx)

  expect_equal(ex$input_class, "tbl_ts")
  expect_length(ex$series, 3L)

  # tsibble's own canonical key order is alphabetical, not insertion order --
  # that IS the "preserved" order fs_extract must reproduce.
  grp_order <- vapply(ex$series, function(s) s$key$grp, character(1))
  expect_equal(grp_order, c("a", "b", "c"))

  for (s in ex$series) {
    expect_s3_class(s$key, "tbl_df")
    expect_equal(names(s$key), "grp")
    expect_equal(s$N, 4L)
    expect_length(s$x, 4L)
    expect_s3_class(s$template, "tbl_ts")
  }

  # values line up with the right key
  a_series <- ex$series[[which(grp_order == "a")]]
  expect_equal(a_series$x, kx$value[kx$grp == "a"])
})

test_that("fs_extract: irregular tsibble errors", {
  irr <- tsibble::tsibble(
    dt = as.Date(c("2020-01-01", "2020-01-03", "2020-01-10")),
    value = c(1, 2, 3),
    index = dt,
    regular = FALSE
  )
  expect_error(fs_extract(irr), "regular")
})

test_that("fs_extract: rejects unsupported input classes", {
  expect_error(fs_extract(list(1, 2, 3)), "numeric|ts|tbl_ts")
})
