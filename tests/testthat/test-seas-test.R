# Tests for the detect + specify orchestration (R/seas_test.R): end-to-end
# behaviour on seeded filter-world / white-noise / line series, input-class
# handling, informative errors, and the print methods.

# --- assembly-time integration guard (see test-seas-ssi.R for the rationale):
# restore the real single-series layer in case the concurrently-built
# test-collection.R clobbered the globals with `<<-` stubs. No-op under a real
# package build.
local({
  r_dir <- file.path("..", "..", "R")
  if (file.exists(file.path(r_dir, "seas_test.R"))) {
    for (f in c("classes.R", "seas_test.R", "seas_ssi.R", "seas_adjust.R",
                "methods-print.R")) {
      sys.source(file.path(r_dir, f), envir = globalenv())
    }
  }
})

# Filter-world quarterly DGP: (1 - phi B)(1 - Phi B^s) x = eps -- a band
# spectrum with seasonal peaks at pi/2 and pi. Deliberately self-contained
# (no module internals) so a pipeline bug cannot hide in the fixture too.
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

test_that("filter-world quarterly series is detected, specified band, M = suggest_M", {
  set.seed(2026)
  x   <- sim_filter_q(120, Phi = 0.9)
  tst <- seas_test(x, frequency = 4)

  expect_s3_class(tst, "seas_test")
  expect_true(tst$decision)
  expect_lt(tst$evt$p, 0.05)
  expect_equal(tst$spec$label, "band")
  expect_equal(tst$M, 5)
  expect_equal(tst$M_selection$source, "suggest_M")
  expect_true(is.finite(tst$spec$phase_R))          # block-phase resultant present
})

test_that("white noise is not seasonal at the seeded alpha", {
  set.seed(4)
  x   <- stats::rnorm(120)
  tst <- seas_test(x, frequency = 4)

  expect_false(tst$decision)
  expect_gt(tst$evt$p, 0.05)
})

test_that("a deterministic sinusoid plus noise is specified line", {
  set.seed(7)
  tt  <- 0:119
  x   <- 5 * cos(pi / 2 * tt + 0.4) + stats::rnorm(120, sd = 1)
  tst <- seas_test(x, frequency = 4)

  expect_true(tst$decision)
  expect_equal(tst$spec$label, "line")
  expect_gt(length(tst$partition$H), 0)             # exact harmonic on-grid
})

test_that("ts input preserves class and infers frequency", {
  set.seed(2026)
  x   <- stats::ts(sim_filter_q(120, Phi = 0.9), frequency = 4, start = c(1990, 1))
  tst <- seas_test(x)

  expect_equal(tst$input_class, "ts")
  expect_equal(tst$N, 4L)
  expect_true(tst$decision)
})

test_that("input errors are informative", {
  # numeric with no frequency
  expect_error(seas_test(stats::rnorm(120)), "frequency")
  # too short: length < 3 * frequency
  expect_error(seas_test(stats::rnorm(10), frequency = 4), "at least|required|length")
  # non-finite values
  expect_error(seas_test(c(stats::rnorm(119), NA_real_), frequency = 4), "finite")
})

test_that("print methods run and contain the contract strings", {
  set.seed(2026)
  x   <- sim_filter_q(120, Phi = 0.9)
  tst <- seas_test(x, frequency = 4)

  out_t <- capture.output(print(tst))
  expect_true(any(grepl("seasonality test", out_t)))
  expect_true(any(grepl("spec:", out_t)))
  expect_true(any(grepl("whitener", out_t)))

  adj   <- seas_ssi(tst, phase_rule = "minimum")
  out_a <- capture.output(print(adj))
  expect_true(any(grepl("declared identification", out_a)))
  expect_true(any(grepl("phase rule: minimum", out_a)))
  expect_true(any(grepl("post-adjustment detection p", out_a)))
})

# ===========================================================================
# WAVE-3 additions: comb diagnostic (harmonic_table + comb_note) and the
# whiten_exclusion / whiten_guard threading. (Extend, not replace.)
# ===========================================================================

# helper: pure nonseasonal AR(2) cycle at pi/2 (roots at angle pi/2, radius rho)
sim_ar2_cycle <- function(n, rho, seed, spin = 400L) {
  set.seed(seed)
  a2 <- -rho^2                                  # x_t = -rho^2 x_{t-2} + e
  ntot <- spin + n; e <- stats::rnorm(ntot); x <- numeric(ntot)
  for (t in 3:ntot) x[t] <- a2 * x[t - 2L] + e[t]
  x[(spin + 1L):ntot]
}

test_that("harmonic_table is present with the documented columns", {
  set.seed(2026)
  tst <- seas_test(sim_filter_q(120, Phi = 0.9), frequency = 4)
  ht <- tst$spec$harmonic_table
  expect_s3_class(ht, "tbl_df")
  expect_true(all(c("omega", "label", "excess", "elevated") %in% names(ht)))
  expect_type(ht$elevated, "logical")
  # quarterly P=1 -> harmonics pi/2 and pi (Nyquist kept)
  expect_setequal(ht$label, c("pi/2", "pi"))
})

test_that("a true seasonal band elevates all harmonics and sets no comb_note", {
  set.seed(2026)
  tst <- seas_test(sim_filter_q(120, Phi = 0.9), frequency = 4)
  expect_true(tst$decision)
  expect_true(all(tst$spec$harmonic_table$elevated))
  expect_true(is.na(tst$spec$comb_note))
})

test_that("a pure AR(2) cycle at pi/2 that rejects gets a pi/2-only comb_note (default config)", {
  # THE two-worlds payoff, under the SHIPPED default (guard(3)): even when the
  # whitener cannot fully whiten a narrow rho=0.95 cycle and it false-alarms, the
  # comb diagnostic catches that the excess is concentrated at pi/2 only (pi
  # stays flat) -- flagging a nonseasonal cycle rather than seasonality. Seed 1
  # rejects under the default; skip guard keeps the test robust to a numeric nudge.
  tst <- seas_test(sim_ar2_cycle(120, rho = 0.95, seed = 1), frequency = 4, M = 5)
  expect_equal(tst$whitener$exclusion, "guard")   # confirm the default path
  skip_if_not(isTRUE(tst$decision), "cycle did not reject under the default on this seed")
  elev <- tst$spec$harmonic_table$elevated
  labs <- tst$spec$harmonic_table$label
  expect_true(elev[labs == "pi/2"])          # pi/2 elevated
  expect_false(elev[labs == "pi"])           # pi (Nyquist) not elevated
  expect_false(is.na(tst$spec$comb_note))
  expect_match(tst$spec$comb_note, "pi/2")
  expect_match(tst$spec$comb_note, "nonseasonal cycle")
  # and it renders in print()
  expect_true(any(grepl("comb:", capture.output(print(tst)))))
})

test_that("the comb diagnostic mechanism also holds under the 'bins' geometry", {
  # same pure cycle under the historical geometry (where it false-alarms more
  # readily): the diagnostic still isolates the pi/2-only excess.
  tst <- seas_test(sim_ar2_cycle(120, rho = 0.95, seed = 101), frequency = 4,
                   M = 5, whiten_exclusion = "bins")
  skip_if_not(isTRUE(tst$decision), "cycle did not reject under bins on this seed")
  labs <- tst$spec$harmonic_table$label
  elev <- tst$spec$harmonic_table$elevated
  expect_true(elev[labs == "pi/2"])
  expect_false(elev[labs == "pi"])
  expect_match(tst$spec$comb_note, "pi/2")
})

test_that("whiten_exclusion / whiten_guard thread into the whitener record", {
  set.seed(2026)
  x <- sim_filter_q(120, Phi = 0.9)
  # defaults: guard(3)
  tst_def <- seas_test(x, frequency = 4)
  expect_equal(tst_def$whitener$exclusion, "guard")
  expect_equal(tst_def$whitener$guard, 3L)
  # explicit overrides are honored
  tst_bins <- seas_test(x, frequency = 4, whiten_exclusion = "bins")
  expect_equal(tst_bins$whitener$exclusion, "bins")
  tst_g2 <- seas_test(x, frequency = 4, whiten_guard = 2L)
  expect_equal(tst_g2$whitener$guard, 2L)
})
