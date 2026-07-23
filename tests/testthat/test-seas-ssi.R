# Tests for the adjustment orchestration (R/seas_ssi.R): the required-phase-rule
# teaching error, the not-seasonal passthrough, band/line surgery correctness,
# the post-adjustment self-check, and input-class preservation.

# --- assembly-time integration guard -------------------------------------
# The concurrently-built collection-layer test file (test-collection.R) installs
# STUB seas_test()/seas_ssi()/seas_adjust() into the GLOBAL environment via `<<-`
# so it can exercise R/collection.R's key-mapping in isolation. Under the
# source()-into-globalenv verification harness those stubs overwrite the real
# single-series implementations for every later test file. Re-source the real
# layer here so this file exercises the real functions. Strict no-op under a
# real package build (R CMD check): the source tree is absent (file.exists is
# FALSE) and testthat already isolates each test file's environment there.
local({
  r_dir <- file.path("..", "..", "R")
  if (file.exists(file.path(r_dir, "seas_test.R"))) {
    for (f in c("classes.R", "seas_test.R", "seas_ssi.R", "seas_adjust.R",
                "methods-print.R")) {
      sys.source(file.path(r_dir, f), envir = globalenv())
    }
  }
})

# Filter-world quarterly band DGP (self-contained; see test-seas-test.R).
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

test_that("missing phase_rule raises the two-worlds teaching error", {
  set.seed(2026)
  tst <- seas_test(sim_filter_q(120, Phi = 0.9), frequency = 4)
  err <- tryCatch(seas_ssi(tst), error = function(e) conditionMessage(e))
  expect_match(err, "minimum")
  expect_match(err, "zero")
  expect_match(err, "two-worlds")
})

test_that("a not-seasonal input is returned unchanged, loudly", {
  set.seed(4)
  x   <- stats::rnorm(120)
  tst <- seas_test(x, frequency = 4)
  expect_false(tst$decision)

  expect_message(res <- seas_ssi(tst, phase_rule = "zero"), "UNCHANGED")
  expect_identical(as.numeric(res$adjusted), x)
  expect_true(all(as.numeric(res$seasonal) == 0))
  expect_equal(res$Ghat, 1)
  expect_equal(res$theta, 0)
})

test_that("band path: adjusted + seasonal reconstructs x and seasonal is nonzero", {
  set.seed(5001)
  x   <- sim_filter_q(120, Phi = 0.6)
  adj <- seas_ssi(seas_test(x, frequency = 4), phase_rule = "minimum")

  expect_equal(adj$spec, "band")
  recon <- as.numeric(adj$adjusted) + as.numeric(adj$seasonal)
  expect_lt(max(abs(recon - x)), 1e-10)
  expect_gt(stats::sd(as.numeric(adj$seasonal)), 0)
  # adjusted actually differs from the input (something was removed)
  expect_gt(max(abs(as.numeric(adj$adjusted) - x)), 1e-6)
})

test_that("line path: the removed seasonal component matches the sinusoid", {
  set.seed(7)
  tt   <- 0:119
  sinu <- 5 * cos(pi / 2 * tt + 0.4)
  x    <- sinu + stats::rnorm(120, sd = 1)
  adj  <- seas_ssi(seas_test(x, frequency = 4), phase_rule = "minimum")

  expect_equal(adj$spec, "line")
  expect_gt(cor(as.numeric(adj$seasonal), sinu), 0.9)
})

test_that("post-adjustment detection p exceeds alpha for most band-path seeds (self-check)", {
  ok <- 0L
  for (r in 1:10) {
    set.seed(5000L + r)
    x   <- sim_filter_q(120, Phi = 0.6)
    adj <- suppressMessages(
      seas_adjust(x, frequency = 4, phase_rule = "minimum"))
    if (isTRUE(adj$spec == "band") && !is.na(adj$post_evt) && adj$post_evt > 0.05) {
      ok <- ok + 1L
    }
  }
  expect_gte(ok, 8L)
})

test_that("class is preserved: numeric in -> numeric out, ts in -> ts out (tsp kept)", {
  set.seed(5001)
  xv      <- sim_filter_q(120, Phi = 0.6)
  adj_num <- seas_ssi(seas_test(xv, frequency = 4), phase_rule = "minimum")
  expect_true(is.numeric(adj_num$adjusted))
  expect_false(stats::is.ts(adj_num$adjusted))
  expect_true(is.numeric(adj_num$seasonal))

  xts    <- stats::ts(xv, frequency = 4, start = c(1990, 1))
  adj_ts <- seas_ssi(seas_test(xts), phase_rule = "minimum")
  expect_true(stats::is.ts(adj_ts$adjusted))
  expect_true(stats::is.ts(adj_ts$seasonal))
  expect_equal(stats::tsp(adj_ts$adjusted), stats::tsp(xts))
})
