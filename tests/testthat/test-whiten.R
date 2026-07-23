# Tests for the M0-Whittle whitener / recolorer (R/whiten.R).
#
# Tier 1: round-trip exactness, d-auto behavior, whitening quality, stationarity,
#         edge cases.
# Tier 2 (the package's signature guarantee): the M0-Whittle whitener cannot
#         absorb a planted seasonal peak, contrasted against an unrestricted
#         full-sample Whittle fit that can.

# ---------------------------------------------------------------------------
# Independent test-side helpers (deliberately NOT reusing module internals so a
# partition bug cannot hide in both the code under test and the verification).
# ---------------------------------------------------------------------------

# Seasonal-bin / donor-bin periodogram mass ratio of a residual series.
# Equal-width bins on (0, pi]; seasonal bins are those containing BOTH pi/2 and
# pi (quarterly, exclude_nyquist = FALSE) -- forgetting the Nyquist peak at pi
# would drop bin 5 and undercount seasonal mass.
seasonal_donor_ratio <- function(e, N = 4L, P = 1L, M = 5L) {
  n    <- length(e)
  Jmax <- (n - 1L) %/% 2L
  j    <- seq_len(Jmax)
  w    <- 2 * pi * j / n
  I    <- (Mod(stats::fft(e - mean(e)))^2 / n)[j + 1L]
  breaks <- seq(0, pi, length.out = M + 1L)
  bin_pos <- cut(w, breaks, right = TRUE, include.lowest = FALSE, labels = FALSE)
  # seasonal frequencies for quarterly, P = 1: pi/2 and pi (Nyquist kept)
  seas_w   <- 2 * pi * seq_len((P * N) %/% 2L) / (P * N)
  seas_bin <- unique(cut(seas_w, breaks, right = TRUE, include.lowest = FALSE,
                         labels = FALSE))
  J1 <- which(bin_pos %in% seas_bin)
  J0 <- which(!(bin_pos %in% seas_bin))
  mean(I[J1]) / mean(I[J0])
}

# Unrestricted full-sample Whittle AR(p) fit: the SAME Whittle machinery as the
# module but evaluated on ALL positive Fourier ordinates instead of the M0 set.
# The only difference from fs_whiten is the ordinate set -- isolating the
# objective as the mechanism of seasonality blindness.
unrestricted_whittle_ar <- function(y0, p) {
  n    <- length(y0)
  Jmax <- (n - 1L) %/% 2L
  j    <- seq_len(Jmax)
  w    <- 2 * pi * j / n
  I    <- (Mod(stats::fft(y0))^2 / n)[j + 1L]
  obj  <- function(theta) {
    kappa <- pmax(pmin(tanh(theta), 1 - 1e-6), -(1 - 1e-6))
    phi   <- .fs_pacf_to_ar(kappa)
    g     <- .fs_ar_gain(phi, w)
    s2    <- mean(I / g)
    length(I) * log(s2) + sum(log(g))
  }
  best_Q <- Inf; best_par <- rep(0, p)
  starts <- c(list(rep(0, p)), list(rep(0.5, p)), list(rep(-0.5, p)))
  for (s in starts) {
    fit <- tryCatch(suppressWarnings(
      stats::optim(s, obj, method = "BFGS",
                   control = list(reltol = 1e-10, maxit = 500L))),
      error = function(e) NULL)
    if (!is.null(fit) && is.finite(fit$value) && fit$value < best_Q) {
      best_Q <- fit$value; best_par <- fit$par
    }
  }
  .fs_pacf_to_ar(pmax(pmin(tanh(best_par), 1 - 1e-6), -(1 - 1e-6)))
}

roots_outside <- function(ar) {
  if (length(ar) == 0L) return(TRUE)
  all(Mod(polyroot(c(1, -ar))) > 1)
}

# ---------------------------------------------------------------------------
# Test 1 -- round-trip exactness (fs_recolor inverts fs_whiten to 1e-10)
# ---------------------------------------------------------------------------

test_that("round trip is exact for the d = 0 path (forced 'none')", {
  set.seed(11)
  x  <- as.numeric(stats::arima.sim(list(ar = c(0.5, -0.3)), n = 120))
  wh <- fs_whiten(x, N = 4L, P = 1L, M = 5L, d = "none", ar_max = 3L)
  xr <- fs_recolor(wh$e, wh, x)

  expect_equal(wh$d, 0L)
  expect_length(xr, length(x))
  expect_lt(max(abs(xr - x)), 1e-10)
  # first (d + p) values equal x exactly (anchor convention)
  k <- wh$d + wh$p
  expect_equal(xr[seq_len(k)], x[seq_len(k)])
  expect_equal(wh$anchors, x[seq_len(k)])
})

test_that("round trip is exact for the d = 1 path ('auto' fires and forced 'first')", {
  set.seed(12)
  # I(1) with AR(1) increments: strongly near-integrated (lag-1 ACF ~ 1, trips
  # auto) and diff(x) has AR structure, exercising the d = 1, p > 0 path.
  x_rw <- cumsum(as.numeric(stats::arima.sim(list(ar = 0.5), n = 120)))
  expect_gt(stats::acf(x_rw, lag.max = 1L, plot = FALSE)$acf[2L], 0.95)

  wh_auto <- fs_whiten(x_rw, N = 4L, P = 1L, M = 5L, d = "auto", ar_max = 3L)
  expect_equal(wh_auto$d, 1L)                       # near-integrated -> differenced
  xr_auto <- fs_recolor(wh_auto$e, wh_auto, x_rw)
  expect_length(xr_auto, length(x_rw))
  expect_lt(max(abs(xr_auto - x_rw)), 1e-10)
  k <- wh_auto$d + wh_auto$p
  expect_equal(xr_auto[seq_len(k)], x_rw[seq_len(k)])

  wh_first <- fs_whiten(x_rw, N = 4L, P = 1L, M = 5L, d = "first", ar_max = 3L)
  expect_equal(wh_first$d, 1L)
  xr_first <- fs_recolor(wh_first$e, wh_first, x_rw)
  expect_lt(max(abs(xr_first - x_rw)), 1e-10)
})

test_that("round trip holds when p = 0 is selected", {
  set.seed(13)
  x  <- stats::rnorm(120)                            # white noise -> p should be 0
  wh <- fs_whiten(x, N = 4L, P = 1L, M = 5L, d = "none", ar_max = 3L)
  xr <- fs_recolor(wh$e, wh, x)
  expect_lt(max(abs(xr - x)), 1e-10)
})

# ---------------------------------------------------------------------------
# Test 2 -- d-auto decision and d_reason strings
# ---------------------------------------------------------------------------

test_that("d = 'auto' picks 0 for white noise and 1 for a random walk", {
  set.seed(21)
  x_wn <- stats::rnorm(120)
  wh_wn <- fs_whiten(x_wn, N = 4L, P = 1L, M = 5L, d = "auto")
  expect_equal(wh_wn$d, 0L)
  expect_true(nzchar(wh_wn$d_reason))
  expect_match(wh_wn$d_reason, "ACF")
  expect_match(wh_wn$d_reason, "no differencing")

  set.seed(22)
  x_rw <- cumsum(stats::rnorm(120))
  wh_rw <- fs_whiten(x_rw, N = 4L, P = 1L, M = 5L, d = "auto")
  expect_equal(wh_rw$d, 1L)
  expect_true(nzchar(wh_rw$d_reason))
  expect_match(wh_rw$d_reason, "differenced once")

  # forced modes report their own reason
  wh_none  <- fs_whiten(x_rw, N = 4L, P = 1L, M = 5L, d = "none")
  wh_first <- fs_whiten(x_wn, N = 4L, P = 1L, M = 5L, d = "first")
  expect_equal(wh_none$d, 0L);  expect_match(wh_none$d_reason, "none")
  expect_equal(wh_first$d, 1L); expect_match(wh_first$d_reason, "first")
})

# ---------------------------------------------------------------------------
# Test 3 -- whitening quality on a seasonal-free AR(1) series
# ---------------------------------------------------------------------------

test_that("whitened AR(1) residuals are spectrally flat and p is usually 1", {
  reps   <- 20L
  ratios <- numeric(reps)
  p_sel  <- integer(reps)
  for (r in seq_len(reps)) {
    set.seed(2000L + r)
    x  <- as.numeric(stats::arima.sim(list(ar = 0.6), n = 120))
    wh <- fs_whiten(x, N = 4L, P = 1L, M = 5L, d = "none", ar_max = 3L)
    ratios[r] <- seasonal_donor_ratio(wh$e)
    p_sel[r]  <- wh$p
    expect_length(wh$bic_path, 4L)
    expect_true(all(is.finite(wh$bic_path)))
    expect_equal(unname(which.min(wh$bic_path)) - 1L, wh$p)   # selection = argmin BIC
    expect_true(roots_outside(wh$ar))
  }
  # flat spectrum: seasonal/donor mass ratio near 1 for every rep, and on average
  expect_true(all(ratios >= 0.4 & ratios <= 2.5))
  expect_gt(mean(ratios), 0.4)
  expect_lt(mean(ratios), 2.5)
  # "usually 1": a clear majority of reps select AR(1)
  expect_gte(sum(p_sel == 1L), 12L)
})

# ---------------------------------------------------------------------------
# Test 4 -- THE PLANTED-PEAK BLINDNESS GUARANTEE (tier 2)
#   DGP: (1 - 0.5 B)(1 - 0.9 B^4) x = eps  (quarterly), i.e. AR(5) with a strong
#   seasonal factor. The M0 whitener must NOT absorb the seasonality; an
#   unrestricted AR(5) does.
# ---------------------------------------------------------------------------

test_that("M0-Whittle whitener does not absorb a planted seasonal peak (guarantee)", {
  reps     <- 20L
  ar_dgp   <- c(0.5, 0, 0, 0.9, -0.45)         # (1-0.5B)(1-0.9B^4)
  ratio_m0 <- numeric(reps)
  for (r in seq_len(reps)) {
    set.seed(1000L + r)
    x  <- as.numeric(stats::arima.sim(list(ar = ar_dgp), n = 120,
                                      n.start = 400, sd = 1))
    # exercise the DEFAULT path: d = "auto". The seasonal (1-0.9B^4) factor puts
    # its mass at lags 4, 8, ..., contributing ~0 to lag-1 correlation, so auto
    # (driven by lag-1 ACF) correctly chooses d = 0 here -- no differencing.
    wh <- fs_whiten(x, N = 4L, P = 1L, M = 5L, d = "auto", ar_max = 3L)
    expect_equal(wh$d, 0L)                        # auto leaves this series in levels
    ratio_m0[r] <- seasonal_donor_ratio(wh$e)
    expect_true(roots_outside(wh$ar))            # test 5: every fit stationary
  }
  # residual seasonal mass stays large: the whitener left the peak intact
  expect_gte(sum(ratio_m0 > 5), 18L)
})

test_that("contrast: an unrestricted full-sample AR(5) DOES absorb the peak", {
  reps      <- 20L
  ar_dgp    <- c(0.5, 0, 0, 0.9, -0.45)
  ratio_unr <- numeric(reps)
  for (r in seq_len(reps)) {
    set.seed(1000L + r)                          # SAME series as the M0 test
    x   <- as.numeric(stats::arima.sim(list(ar = ar_dgp), n = 120,
                                       n.start = 400, sd = 1))
    y0  <- x - mean(x)
    phi <- unrestricted_whittle_ar(y0, p = 5L)   # unrestricted, all ordinates
    e   <- .fs_ar_residuals(y0, phi)
    ratio_unr[r] <- seasonal_donor_ratio(e)
  }
  # the failure mode this design prevents: the peak is largely gone
  expect_lt(stats::median(ratio_unr), 2.5)
})

# ---------------------------------------------------------------------------
# Test 5 -- stationarity across heterogeneous fits
# ---------------------------------------------------------------------------

test_that("fitted AR roots lie outside the unit circle across series types", {
  make <- list(
    ar1  = function() as.numeric(stats::arima.sim(list(ar = 0.8), n = 120)),
    ar2n = function() as.numeric(stats::arima.sim(list(ar = c(0.2, -0.5)), n = 120)),
    seas = function() as.numeric(stats::arima.sim(
             list(ar = c(0.5, 0, 0, 0.9, -0.45)), n = 120, n.start = 400)),
    wn   = function() stats::rnorm(120),
    rw   = function() cumsum(stats::rnorm(120))
  )
  for (nm in names(make)) {
    for (r in 1:4) {
      set.seed(3000L + r)
      x  <- make[[nm]]()
      wh <- fs_whiten(x, N = 4L, P = 1L, M = 5L, d = "auto", ar_max = 3L)
      expect_true(roots_outside(wh$ar),
                  info = sprintf("series = %s, rep = %d, p = %d", nm, r, wh$p))
    }
  }
})

# ---------------------------------------------------------------------------
# Test 6 -- edge cases error informatively
# ---------------------------------------------------------------------------

test_that("too-short and constant series raise informative errors", {
  expect_error(fs_whiten(stats::rnorm(10), N = 4L, P = 1L, M = 5L),
               regexp = "too short", ignore.case = TRUE)
  expect_error(fs_whiten(rep(5, 120), N = 4L, P = 1L, M = 5L),
               regexp = "constant", ignore.case = TRUE)
  # a linear trend is constant after differencing -> caught on the forced path
  expect_error(fs_whiten(seq_len(120), N = 4L, P = 1L, M = 5L, d = "first"),
               regexp = "no variation", ignore.case = TRUE)
})

# ===========================================================================
# WAVE-3 additions: de-biased Whittle expected periodogram, and the
# guard/bins exclusion geometries. (Extend, not replace.)
# ===========================================================================

# Test 7 -- the de-biased expected periodogram helper (.fs_expected_pgram)
#   It must equal E[I(w_j)] = sum_{|k|<n} (1-|k|/n) rho(k) e^{-ikw_j} exactly at
#   the Fourier frequencies (validated against a brute-force reference), obey the
#   Parseval invariant mean(Ebar) = rho(0) = 1, and -- the point of debiasing --
#   sit BELOW the infinite-sample AR gain at a finite-n spectral peak.
test_that(".fs_expected_pgram matches the brute-force expected periodogram", {
  ref <- function(phi, n) {
    rho <- stats::ARMAacf(ar = phi, lag.max = n - 1L)
    w   <- 2 * pi * (0:(n - 1)) / n
    vapply(w, function(wj) {
      k <- -(n - 1):(n - 1)
      Re(sum((1 - abs(k) / n) * rho[abs(k) + 1L] * exp(-1i * k * wj)))
    }, numeric(1))
  }
  for (phi in list(0.5, c(0.2, -0.5), c(0, -0.7225), c(0.5, 0, 0, 0.9, -0.45))) {
    n <- 40L
    Ebar <- .fs_expected_pgram(phi, n)
    expect_lt(max(abs(Ebar - ref(phi, n))), 1e-9)
    expect_equal(mean(Ebar), 1, tolerance = 1e-10)   # Parseval invariant
  }
  # p = 0 (white): flat expected periodogram
  expect_equal(.fs_expected_pgram(numeric(0), 32L), rep(1, 32L))
})

test_that("de-biased expected periodogram is below the AR gain at a finite-n peak", {
  phi <- c(0, -0.7225); n <- 120L            # AR(2) cycle peaking at pi/2 (j = 30)
  w    <- 2 * pi * (0:(n - 1)) / n
  Ebar <- .fs_expected_pgram(phi, n)
  gain <- .fs_ar_gain(phi, w)
  jpk  <- 30L + 1L                            # peak ordinate in 1-indexed DFT coords
  # both normalized to their own mean; the debiased peak is blurred (lower)
  expect_lt(max(Ebar) / mean(Ebar), max(gain) / mean(gain))
})

# Test 8 -- exclusion geometries (guard default vs bins)
test_that("guard vs bins exclusion: both round-trip and record their config", {
  set.seed(31)
  x  <- as.numeric(stats::arima.sim(list(ar = c(0.5, -0.3)), n = 120))
  wh_g <- fs_whiten(x, N = 4L, P = 1L, M = 5L, d = "none")               # default guard(3)
  wh_b <- fs_whiten(x, N = 4L, P = 1L, M = 5L, d = "none", exclusion = "bins")

  # default config is recorded on the record
  expect_equal(wh_g$exclusion, "guard")
  expect_equal(wh_g$guard, 3L)
  expect_equal(wh_b$exclusion, "bins")

  # both invert exactly
  expect_lt(max(abs(fs_recolor(wh_g$e, wh_g, x) - x)), 1e-10)
  expect_lt(max(abs(fs_recolor(wh_b$e, wh_b, x) - x)), 1e-10)
})

test_that("the guard set keeps strictly more ordinates than the bins (M0) set", {
  g   <- fourier_grid(118L)
  part <- mbin_partition(5L, 1L, 4L, g$omega_pos, band = "full", exclude_nyquist = FALSE)
  isx <- index_sets(part, g)
  jfit_guard <- .fs_guard_include(g, P = 1L, N = 4L, guard = 3L)
  # guard removes only ~ (2*guard+1) ordinates near each harmonic; bins removes
  # whole seasonal bins, so the M0 (bins) set is a strict subset of what guard keeps.
  expect_true(length(jfit_guard) > length(isx$J0))
  expect_true(all(isx$J0 %in% jfit_guard))
  # guard(0) excludes only the exact on-grid harmonic ordinates (pi/2 at j=29+..,
  # Nyquist at j=59): a wider band than bins.
  jfit_g0 <- .fs_guard_include(g, P = 1L, N = 4L, guard = 0L)
  expect_true(length(jfit_g0) > length(jfit_guard))
})

test_that("the blindness guarantee holds under BOTH exclusion geometries (guarantee)", {
  # Test 4's planted comb DGP, asserted for guard AND bins: neither absorbs the
  # concentrated seasonal comb (guard protects the harmonic tips; bins excludes
  # the whole seasonal bins). Extends the guarantee across the new default.
  reps   <- 20L
  ar_dgp <- c(0.5, 0, 0, 0.9, -0.45)
  ratio_guard <- numeric(reps); ratio_bins <- numeric(reps)
  for (r in seq_len(reps)) {
    set.seed(1000L + r)
    x <- as.numeric(stats::arima.sim(list(ar = ar_dgp), n = 120,
                                     n.start = 400, sd = 1))
    ratio_guard[r] <- seasonal_donor_ratio(
      fs_whiten(x, N = 4L, P = 1L, M = 5L, d = "auto", exclusion = "guard")$e)
    ratio_bins[r] <- seasonal_donor_ratio(
      fs_whiten(x, N = 4L, P = 1L, M = 5L, d = "auto", exclusion = "bins")$e)
  }
  expect_gte(sum(ratio_guard > 5), 18L)
  expect_gte(sum(ratio_bins  > 5), 18L)
})
