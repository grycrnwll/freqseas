# Regression tests for seasonal alignment of the surgery grid.
#
# Everything downstream of fs_whiten() -- fourier_grid(), mbin_partition(),
# index_sets(), the periodogram and fs_surgery() -- lives on the length-n_e
# Fourier grid of the AR residual. Nothing used to constrain n_e = n - d - p to
# be a multiple of the seasonal period N, and after differencing plus an AR fit
# it essentially never was. Off-grid, the seasonal harmonics are not Fourier
# ordinates: the exact-harmonic set H collapses to at most {Nyquist}, a
# deterministic line leaks through the Dirichlet kernel and elevates the
# shoulder ordinates so the Stage B shoulder EVT calls a line a band, and the
# band surgery then cannot remove the leaked line.
#
# The fix trims the residual vector to a whole number of seasonal cycles from
# the FRONT, holding the n_trim < N leading residuals out of the surgery and
# splicing them back in fs_recolor(). These tests lock the behavioural contract
# (public API where possible) rather than the trimming mechanics.
#
# See notes/2026-07-23_airpassengers_residual_seasonality.md.

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

# ---------------------------------------------------------------------------
# Test-side helpers (deliberately independent of module internals).
# ---------------------------------------------------------------------------

# Deterministic-line monthly DGP: one pure sinusoid at the seasonal fundamental
# plus light noise. The SAME seed for every n, so the first min(n) noise draws
# are shared and length is the only thing that varies across the calls below.
sim_line_monthly <- function(n, amp = 3, sd = 0.5, seed = 4242L) {
  set.seed(seed)
  stats::ts(amp * cos(2 * pi * (1:n) / 12) + stats::rnorm(n, sd = sd),
            frequency = 12)
}

# LOCAL excess of periodogram mass at a seasonal harmonic: the periodogram at
# omega = 2*pi*h/N over the median periodogram at the near-neighbour ordinates
# (offsets +-2..+-5). The GLOBAL median understates the local floor wherever
# trend mass sits nearby, and `post_evt` -- computed on the re-whitened residual
# grid, where any leakage is spread -- is not a proxy for raw-scale seasonal
# mass at all: the note records post_evt = 0.9994 alongside a 2.3x local peak.
# Requires length(z) %% N == 0 so that the harmonic is an exact ordinate.
local_excess <- function(z, h = 1L, N = 12L) {
  z <- as.numeric(z)
  n <- length(z)
  stopifnot(n %% N == 0L)
  I <- Mod(stats::fft(z - mean(z)))^2 / n        # 0-indexed ordinates
  k <- n * h / N
  nb <- k + c(-5:-2, 2:5)
  I[k + 1L] / stats::median(I[nb + 1L])
}

# The lengths used for the alignment sweep: n - d - p lands on several
# different residues mod 12, so the unaligned grid would have been misaligned
# in several different ways.
align_lengths <- c(144L, 145L, 146L, 147L, 150L)

# ---------------------------------------------------------------------------
# Test 1 -- the surgery grid is a whole number of seasonal cycles
# ---------------------------------------------------------------------------

test_that("the surgery grid holds a whole number of seasonal cycles", {
  residues <- integer(0)
  for (n in align_lengths) {
    tst <- seas_test(sim_line_monthly(n))
    wh  <- whitener(tst)

    # the grid everything downstream lives on is aligned ...
    expect_equal(length(tst$e) %% 12L, 0L,
                 info = sprintf("n = %d", n))
    expect_equal(length(tst$e), wh$n_e, info = sprintf("n = %d", n))
    # the downstream objects really do live on that same aligned grid
    expect_equal(tst$pgram$n, length(tst$e), info = sprintf("n = %d", n))
    expect_equal(tst$partition$grid$n, length(tst$e), info = sprintf("n = %d", n))

    # ... and the UNALIGNED length it would otherwise have had is recorded
    # from d and p (not from n_e, which is post-trim), so this test keeps
    # exercising misalignment even if the d/p selection shifts.
    residues <- c(residues, (n - wh$d - wh$p) %% 12L)
  }
  # the sweep really does cover several distinct misalignments
  expect_gt(length(unique(residues)), 1L)
})

# ---------------------------------------------------------------------------
# Test 2 -- the exact-harmonic set is complete
# ---------------------------------------------------------------------------

test_that("all six monthly seasonal harmonics are exact Fourier ordinates", {
  for (n in align_lengths) {
    tst <- seas_test(sim_line_monthly(n))
    # h = 1..6 for N = 12, Nyquist (h = 6) included. Off-grid this collapsed to
    # {Nyquist} for even n_e and to the empty set for odd n_e.
    expect_equal(length(tst$partition$H), 6L, info = sprintf("n = %d", n))
    # H must index the exact harmonic ordinates of the aligned grid
    n_e <- length(tst$e)
    expect_setequal(tst$partition$H, as.integer(n_e * (1:6) / 12))
  }
})

# ---------------------------------------------------------------------------
# Test 3 -- specification is length-invariant for a deterministic line
#   The headline regression. Same DGP, same seed, only n varies. Before the fix
#   n = 144/145/146 gave line / band / band, and post_evt 0.853 / 6.0e-15 /
#   3.5e-05 -- correctness was a lottery on n_e mod 12.
# ---------------------------------------------------------------------------

test_that("a deterministic line is specified 'line' at every series length", {
  for (n in c(144L, 145L, 146L)) {
    tst <- seas_test(sim_line_monthly(n))
    expect_true(tst$decision, info = sprintf("n = %d", n))
    expect_equal(tst$spec$label, "line", info = sprintf("n = %d", n))
  }
})

test_that("line adjustment leaves no detected seasonality at any length", {
  for (n in c(144L, 145L, 146L)) {
    adj <- suppressMessages(
      seas_ssi(seas_test(sim_line_monthly(n)), phase_rule = "minimum"))
    expect_equal(adj$spec, "line", info = sprintf("n = %d", n))
    expect_gt(adj$post_evt, 0.05)
  }
})

# ---------------------------------------------------------------------------
# Test 4 -- log(AirPassengers), the reported case
# ---------------------------------------------------------------------------

test_that("log(AirPassengers) is specified 'line' with a complete harmonic set", {
  tst <- seas_test(log(datasets::AirPassengers))
  expect_true(tst$decision)
  # the phase resultant pointed at a frozen pattern all along; off-grid leakage
  # was what drove the shoulder EVT to reject and label it "band".
  expect_equal(tst$spec$label, "line")
  expect_equal(length(tst$partition$H), 6L)
  expect_equal(length(tst$e) %% 12L, 0L)
})

test_that("log(AirPassengers) passes its own post-adjustment self-check", {
  adj <- suppressMessages(seas_adjust(log(datasets::AirPassengers)))
  expect_equal(adj$spec, "line")
  expect_gt(adj$post_evt, 0.05)             # was 0.0101 (< alpha) before the fix
})

# ---------------------------------------------------------------------------
# Test 5 -- raw-scale seasonal mass is actually gone
#   This is the assertion that answers the original bug report. post_evt is not
#   a sufficient proxy (see local_excess() above), so measure the local excess
#   on the adjusted series directly.
# ---------------------------------------------------------------------------

test_that("adjustment removes the raw-scale peak at the seasonal fundamental", {
  x   <- log(datasets::AirPassengers)
  adj <- suppressMessages(seas_adjust(x))

  # calibration: the helper must see the peak that is there to begin with,
  # otherwise "small excess after adjustment" would pass for the wrong reason.
  expect_gt(local_excess(x, h = 1L), 10)          # measured ~16.7 on the input

  # the contract: after adjustment the fundamental is not a local peak.
  # Raw 16.7; the pre-fix band path left 4.3; alignment gives ~1.2.
  expect_lt(local_excess(as.numeric(adj$adjusted), h = 1L), 2)
})

# ---------------------------------------------------------------------------
# Test 6 -- the round-trip invariant survives the head splice
#   fs_recolor() must still invert fs_whiten() exactly now that it is handed a
#   trimmed innovation vector and has to splice wh$e_head back on. This is the
#   invariant most at risk from the trim.
# ---------------------------------------------------------------------------

test_that("fs_recolor inverts fs_whiten exactly on the trimmed grid (d = 1, p >= 1)", {
  set.seed(61)
  # I(1) with AR(1) increments: d = "auto" fires and p > 0 is selected.
  x  <- cumsum(as.numeric(stats::arima.sim(list(ar = 0.6), n = 149)))
  wh <- fs_whiten(x, N = 12L, P = 1L, M = 19L, d = "auto", ar_max = 3L)
  expect_equal(wh$d, 1L)
  expect_gte(wh$p, 1L)
  expect_gt(wh$n_trim, 0L)                       # the splice is exercised
  expect_length(wh$e_head, wh$n_trim)

  xr <- fs_recolor(wh$e, wh, x)
  expect_length(xr, length(x))
  expect_lt(max(abs(xr - x)), 1e-10)
  # anchor convention is untouched by the trim
  k <- wh$d + wh$p
  expect_equal(xr[seq_len(k)], x[seq_len(k)])
})

test_that("fs_recolor inverts fs_whiten exactly on the trimmed grid (d = 0, p = 0)", {
  set.seed(62)
  x  <- stats::rnorm(130)                        # white noise -> p = 0
  wh <- fs_whiten(x, N = 12L, P = 1L, M = 19L, d = "none", ar_max = 3L)
  expect_equal(wh$d, 0L)
  expect_equal(wh$p, 0L)
  expect_gt(wh$n_trim, 0L)                       # 130 is not a multiple of 12
  expect_length(wh$e_head, wh$n_trim)

  xr <- fs_recolor(wh$e, wh, x)
  expect_lt(max(abs(xr - x)), 1e-10)
})

# ---------------------------------------------------------------------------
# Test 7 -- the trim is bounded and reported
#   Up to N - 1 leading residuals go unadjusted. That cost must be bounded,
#   recorded on the whitener record, and visible in the output: the head is
#   held out of the surgery, so it comes back bit-identical to the input.
# ---------------------------------------------------------------------------

test_that("the trim is bounded by N and recorded on the whitener record", {
  for (n in align_lengths) {
    wh <- whitener(seas_test(sim_line_monthly(n)))
    expect_true(is.numeric(wh$n_trim) || is.integer(wh$n_trim))
    expect_gte(wh$n_trim, 0L)
    expect_lt(wh$n_trim, 12L)                    # strictly less than N
    expect_length(wh$e_head, wh$n_trim)
    expect_equal(wh$n_e + wh$n_trim + wh$d + wh$p, n)
  }
})

test_that("the unadjusted head of the adjusted series equals the input exactly", {
  x   <- log(datasets::AirPassengers)
  fit <- suppressMessages(seas_adjust(x))
  wh  <- whitener(fit)

  skip_if_not(isTRUE(wh$n_trim > 0L),
              "grid was already aligned; no head is held out")

  k <- wh$n_trim + wh$d + wh$p
  expect_lt(wh$n_trim, 12L)
  expect_equal(as.numeric(adjusted(fit))[seq_len(k)], as.numeric(x)[seq_len(k)])
  # and the surgery did change the rest of the series
  expect_gt(max(abs(as.numeric(adjusted(fit)) - as.numeric(x))), 1e-6)
})

# ---------------------------------------------------------------------------
# Test 8 -- the Nyquist comb tooth on the band path
#   Alignment has a second-order consequence. Trimming n_e to a whole number of
#   seasonal cycles forces n_e EVEN whenever N is even, so omega = pi now always
#   lands exactly on the Nyquist ordinate and is always an exact seasonal
#   harmonic. fs_gain() has no periodogram-based estimate to offer there -- its
#   `pgram_pos` has length Jmax = (n_e - 1) %/% 2, which for even n_e has no
#   Nyquist entry -- so it returns the neutral gain 1. Untreated, the band
#   branch therefore leaves the ENTIRE comb tooth at pi passing through
#   unadjusted; seas_ssi() fills the gain in from the raw periodogram.
#
#   Before the alignment trim this was harmless: omega = pi was generically OFF
#   the grid and its mass sat on ordinates that J1 already covered. The defect
#   was live for a while with no test asserting it directly -- it surfaced only
#   through band-path post-adjustment self-checks in other files, which is a
#   slow and indirect way to learn about it. These two blocks assert it at the
#   point of failure.
# ---------------------------------------------------------------------------

# Filter-world quarterly band DGP, (1 - phi B)(1 - Phi B^s) with s = 4: a band
# spectrum whose seasonal comb carries a tooth at omega = pi. Self-contained
# and seeded, like the sibling files' fixtures.
sim_band_q <- function(n, seed, Phi = 0.8, phi = 0.5, s = 4L, spin = 400L) {
  set.seed(seed)
  ntot <- spin + n
  eps  <- stats::rnorm(ntot)
  z    <- numeric(ntot)
  for (t in (s + 2L):ntot) {
    z[t] <- phi * z[t - 1L] + Phi * z[t - s] - phi * Phi * z[t - s - 1L] + eps[t]
  }
  out <- z[(spin + 1L):ntot]
  stats::ts(out - mean(out), frequency = 4)
}

# Seeds whose raw Nyquist tooth is unambiguous (local excess 8.9 - 98), so the
# "it was removed" assertion cannot pass merely because there was nothing there.
nyq_seeds <- c(5003L, 5005L, 5006L, 5009L)

test_that("band path applies a real gain at the Nyquist comb tooth", {
  for (sd in nyq_seeds) {
    x   <- sim_band_q(120L, seed = sd)
    tst <- seas_test(x)
    expect_equal(tst$spec$label, "band", info = sprintf("seed = %d", sd))

    n_e <- length(tst$e)
    expect_equal(n_e %% 4L, 0L, info = sprintf("seed = %d", sd))
    # even N and an aligned grid => omega = pi IS an exact seasonal harmonic
    nyq_j <- n_e %/% 2L
    expect_true(nyq_j %in% tst$partition$H, info = sprintf("seed = %d", sd))

    adj <- suppressMessages(seas_ssi(tst, phase_rule = "minimum"))
    expect_length(adj$Ghat, n_e)
    # The contract: a gain was actually estimated at pi. Left at fs_gain()'s
    # neutral value this is EXACTLY 1 and the tooth passes through untouched,
    # so the separation here needs no tolerance at all.
    expect_gt(adj$Ghat[nyq_j + 1L], 1)
  }
})

test_that("band adjustment removes the raw-scale peak at omega = pi", {
  for (sd in nyq_seeds) {
    x   <- sim_band_q(120L, seed = sd)
    adj <- suppressMessages(seas_adjust(x, phase_rule = "minimum"))

    # calibration: the tooth is genuinely present first (h = N/2 = 2 -> pi)
    expect_gt(local_excess(x, h = 2L, N = 4L), 4)
    # contract: it is gone afterwards. Note the failure mode is not a mild
    # under-correction -- with the tooth untreated the band surgery strips the
    # surrounding shoulders and leaves the peak MORE isolated, so the excess
    # rises rather than falls (measured 27x - 271x against 0.07 - 2.3 here).
    expect_lt(local_excess(as.numeric(adj$adjusted), h = 2L, N = 4L), 4)
  }
})
