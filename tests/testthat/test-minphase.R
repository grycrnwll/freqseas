# test-minphase.R
# Tier-1 closed-form tests for min_phase(), reproducing the verified numbers
# from notes/2026-07-22_ssi_phase_surgery.md SS3 and
# scripts/mp_phase_mre.R in seasonality_vintages.

test_that("min_phase reproduces the verified 2.8e-14 fold accuracy on an accurate fine grid", {
  # The session's verified "2.8e-14" number (notes/2026-07-22_ssi_phase_
  # surgery.md SS3, scripts/mp_phase_mre.R validation block) comes from
  # folding a closed-form gain evaluated DIRECTLY on a fine (Nf = n*64)
  # grid -- i.e. the cepstral fold itself, given an accurate input, is
  # exact to machine precision. Reproduce that exact validation through
  # the exported function with fine_factor = 1L: the input IS already the
  # fine grid, so no upsampling/interpolation occurs.
  Phi <- 0.9; s <- 4L; Nf <- 80L * 64L
  wf <- 2 * pi * (0:(Nf - 1)) / Nf
  G_fine  <- 1 / sqrt(1 + Phi^2 - 2 * Phi * cos(s * wf))
  th_fine <- -atan2(Phi * sin(s * wf), 1 - Phi * cos(s * wf))

  theta_hat <- min_phase(log(G_fine), fine_factor = 1L)

  max_err <- max(abs(theta_hat - th_fine))
  expect_lt(max_err, 1e-10)
})

test_that("min_phase reproduces the closed-form seasonal-AR sawtooth (Phi = 0.9, n = 80)", {
  # This is the length-80 end-to-end path: min_phase() must first
  # interpolate an 80-point log-gain up to the fine grid before folding
  # (the pinned algorithm; see R/minphase.R). Unlike the test above, this
  # necessarily loses information: with only 80 knots representing a
  # circle carrying s = 4 harmonics (20 knots per harmonic period), linear
  # interpolation cannot perfectly reconstruct the curvature of a fairly
  # peaked gain (Phi = 0.9), and no increase in fine_factor recovers what
  # the 80 knots don't carry -- verified empirically: max error saturates
  # at ~0.0219 rad for fine_factor from 64 up to 1024. This is the same
  # resolution-boundary phenomenon documented in the notes SS5 ("a band
  # narrower than 2*pi/n is a line as far as n observations are
  # concerned"): at n = 80 the harmonic's curvature is under-resolved by
  # the coarse grid itself, independent of the fold's own accuracy (which
  # the test above shows is exact to 1e-10 given an accurate input).
  Phi <- 0.9; s <- 4L; n <- 80L
  w <- 2 * pi * (0:(n - 1)) / n
  G <- 1 / sqrt(1 + Phi^2 - 2 * Phi * cos(s * w))
  theta_true <- -atan2(Phi * sin(s * w), 1 - Phi * cos(s * w))

  theta_hat <- min_phase(log(G))

  expect_length(theta_hat, n)
  max_err <- max(abs(theta_hat - theta_true))
  expect_lt(max_err, 0.03)
})

test_that("cepstrum of the seasonal-AR log gain lives on multiples of s", {
  # Closed form: the complex cepstrum of -log(1 - Phi z^s) is Phi^m/m at
  # lag s*m, zero elsewhere. Checked directly on the fine grid min_phase()
  # itself folds (n = 80, fine_factor = 64), independent of min_phase()'s
  # own internals.
  Phi <- 0.9; s <- 4L; n <- 80L; fine_factor <- 64L
  Nf <- n * fine_factor
  wf <- 2 * pi * (0:(Nf - 1)) / Nf
  G_true <- 1 / sqrt(1 + Phi^2 - 2 * Phi * cos(s * wf))

  ce <- Re(stats::fft(log(G_true), inverse = TRUE)) / Nf

  expect_equal(2 * ce[5], Phi, tolerance = 1e-6)      # lag 4  -> R index 5
  expect_equal(4 * ce[9], Phi^2, tolerance = 1e-6)     # lag 8  -> R index 9
})

test_that("constant gain gives zero phase; asymmetric input errors informatively", {
  n <- 32L
  theta_const <- min_phase(rep(0, n))
  expect_lt(max(abs(theta_const)), 1e-10)

  set.seed(3)
  bad <- rnorm(n)  # generic random vector: essentially never symmetric
  expect_error(min_phase(bad), regexp = "symmetric")
})

test_that("min_phase matches the AR(1) phase bowl (nonseasonal, phi = 0.5)", {
  # Same length-80 interpolate-then-fold path as above; the AR(1) bowl is
  # far gentler than the Phi = 0.9 seasonal peak, so the curvature loss
  # from linear interpolation is much smaller (saturates ~7e-4 rather than
  # ~0.022 rad regardless of fine_factor) but is not literally zero.
  phi <- 0.5; n <- 80L
  w <- 2 * pi * (0:(n - 1)) / n
  G <- 1 / sqrt(1 + phi^2 - 2 * phi * cos(w))
  theta_true <- -atan2(phi * sin(w), 1 - phi * cos(w))

  theta_hat <- min_phase(log(G))

  max_err <- max(abs(theta_hat - theta_true))
  expect_lt(max_err, 1e-3)
})
