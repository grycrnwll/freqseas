# Tests for R/benchmark.R: proportional yearly-total benchmarking (the
# NIPA-style constraint that complete-period totals are invariant to seasonal
# adjustment). Covers both input conventions (ts and numeric + frequency), the
# trailing-incomplete-window edge convention, informative errors on degenerate
# inputs, and the end-to-end invariant on the bundled `ssi_examples` filter
# world (benchmarking a seas_adjust result restores x's annual totals).

# --- assembly-time integration guard (see test-seas-ssi.R for the rationale):
# the concurrently-built test-collection.R installs `<<-` stubs of the real
# seas_test()/seas_ssi()/seas_adjust() layer into the global environment. Under
# the source-into-globalenv verification harness those stubs would clobber the
# real functions the invariant test needs. Re-source the whole R/ tree to
# restore them. No-op under a real package build (R CMD check): the source tree
# is absent and testthat isolates each file's environment there.
local({
  r_dir <- file.path("..", "..", "R")
  if (file.exists(file.path(r_dir, "seas_test.R"))) {
    for (f in list.files(r_dir, pattern = "\\.R$", full.names = TRUE)) {
      sys.source(f, envir = globalenv())
    }
  }
})

# Per-period totals via the same positional windowing benchmark_totals uses
# (matrix colSums, not calendar-aware aggregation), so the assertion cannot
# drift on start quarter.
year_totals <- function(z, N) colSums(matrix(as.numeric(z), nrow = N))

test_that("complete-year totals match the original exactly (ts input, ts output)", {
  set.seed(101)
  orig <- stats::ts(100 + stats::rnorm(24), frequency = 4)   # 6 complete years
  adj  <- orig + stats::rnorm(24)                            # a perturbed series
  m <- benchmark_totals(adj, orig)

  expect_true(stats::is.ts(m))
  expect_equal(stats::tsp(m), stats::tsp(orig))
  expect_equal(year_totals(m, 4), year_totals(orig, 4), tolerance = 1e-10)

  bm <- attr(m, "benchmark")
  expect_equal(bm$period, "year")
  expect_equal(bm$frequency, 4L)
  expect_equal(bm$n_tail_unadjusted, 0L)
  expect_length(bm$scale_factors, 6L)
})

test_that("numeric + frequency input gives the same result as the ts input", {
  set.seed(101)
  orig <- 100 + stats::rnorm(24)
  adj  <- orig + stats::rnorm(24)
  m_num <- benchmark_totals(adj, orig, frequency = 4)
  m_ts  <- benchmark_totals(stats::ts(adj, frequency = 4),
                            stats::ts(orig, frequency = 4))

  expect_true(is.numeric(m_num))
  expect_false(stats::is.ts(m_num))
  expect_equal(year_totals(m_num, 4), year_totals(orig, 4), tolerance = 1e-10)
  expect_equal(as.numeric(m_num), as.numeric(m_ts), tolerance = 1e-12)
})

test_that("a trailing incomplete period is left unadjusted, with a warning", {
  orig <- 100 + as.numeric(1:22)          # 22 obs: 5 complete years + 2 tail
  adj  <- orig + 0.5
  expect_warning(m <- benchmark_totals(adj, orig, frequency = 4),
                 "unadjusted|whole number")

  # the 2 trailing observations are untouched
  expect_equal(tail(as.numeric(m), 2), tail(adj, 2))
  expect_equal(attr(m, "benchmark")$n_tail_unadjusted, 2L)
  # the 20 complete-window observations were benchmarked
  expect_equal(year_totals(as.numeric(m)[1:20], 4),
               year_totals(orig[1:20], 4), tolerance = 1e-10)
})

test_that("proportional (not additive) rescaling: totals scale multiplicatively", {
  # additive matching would ADD (orig_sum - adj_sum)/N to each obs, leaving the
  # within-year *shape* unchanged; proportional matching multiplies, so the
  # ratio of any two values within a window is preserved.
  orig <- c(10, 20, 30, 40)
  adj  <- c(1, 2, 3, 4)                    # same shape, tenth the level
  m <- suppressWarnings(benchmark_totals(adj, orig, frequency = 4))
  expect_equal(as.numeric(m), orig, tolerance = 1e-10)   # x10 factor
  expect_equal(attr(m, "benchmark")$scale_factors, 10, tolerance = 1e-10)
})

test_that("degenerate inputs error informatively", {
  expect_error(benchmark_totals(1:8, 1:10, frequency = 4), "same length")
  expect_error(benchmark_totals(as.numeric(1:8), as.numeric(1:8)), "frequency")
  expect_error(
    benchmark_totals(c(as.numeric(1:7), NA), as.numeric(1:8), frequency = 4),
    "finite")
  # a complete window of `adjusted` sums to zero -> no proportional factor
  expect_error(
    benchmark_totals(c(1, -1, 2, -2, 5, 5, 5, 5), c(1, 1, 1, 1, 2, 2, 2, 2),
                     frequency = 4),
    "zero")
  # shorter than one complete period
  expect_error(benchmark_totals(as.numeric(1:3), as.numeric(1:3), frequency = 4),
               "shorter")
  # mismatched ts frequencies
  expect_error(
    benchmark_totals(stats::ts(1:24, frequency = 4),
                     stats::ts(1:24, frequency = 12)),
    "different")
})

test_that("NIPA invariant: benchmarking a seas_adjust result restores x's annual totals", {
  # load the bundled known-truth data (not lazy-loaded under the source-into-
  # global verification harness; available as package data under R CMD check).
  if (!exists("ssi_examples")) {
    rda <- file.path("..", "..", "data", "ssi_examples.rda")
    if (file.exists(rda)) load(rda) else skip("ssi_examples not available")
  }

  fw <- ssi_examples[ssi_examples$world == "filter", ]
  x  <- fw$x[[1]]
  N  <- stats::frequency(x)

  adj <- suppressMessages(seas_adjust(x, phase_rule = "minimum"))
  a   <- adjusted(adj)
  # before benchmarking, adjustment redistributes within years -> annual totals
  # generally move away from x's:
  expect_gt(max(abs(year_totals(a, N) - year_totals(x, N))), 1e-6)

  m <- benchmark_totals(a, x)
  # after benchmarking, every complete-year total equals x's exactly:
  expect_equal(year_totals(m, N), year_totals(x, N), tolerance = 1e-8)
  expect_true(stats::is.ts(m))
})
