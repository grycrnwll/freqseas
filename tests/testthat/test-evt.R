# Tests for the periodogram / EVT-null machinery: closed-form conventions plus
# a fixed-seed statistical acceptance check on the detection test's size.

test_that("fs_periodogram preserves its statistical conventions", {
  set.seed(1)
  x  <- stats::rnorm(64)
  pg <- fs_periodogram(x, demean = TRUE, standardize = TRUE)

  expect_equal(pg$n, 64L)
  expect_equal(pg$tau_hat, 1)                       # standardized -> unit scale
  expect_equal(pg$pgram_raw, (Mod(pg$dft_raw)^2) / pg$n)
  expect_equal(pg$pgram_std, (Mod(pg$dft_std)^2) / pg$n)
  # Parseval: sum of standardized periodogram ordinates equals the sum of
  # squared standardized deviations, which is exactly (n - 1) = 63.
  expect_equal(sum(pg$pgram_std), 63, tolerance = 1e-8)

  # Without standardizing, tau_hat is the demeaned variance.
  pg2 <- fs_periodogram(x, demean = TRUE, standardize = FALSE)
  expect_equal(pg2$tau_hat, stats::var(x - mean(x)))
  expect_true(is.na(pg2$sd))
})

test_that("evt_critical and evt_pvalue agree at the rejection boundary", {
  N1 <- 8L
  N0 <- 32L
  alpha <- 0.05
  tau <- 1
  cval <- evt_critical(N1, N0, alpha = alpha, tau_hat = tau,
                       alternative = "greater")
  # The p-value evaluated exactly at the critical value equals alpha.
  expect_equal(evt_pvalue(cval, N1, N0, tau_hat = tau, alternative = "greater"),
               alpha, tolerance = 1e-8)
  # Larger statistic -> smaller p-value.
  expect_lt(evt_pvalue(cval + 1, N1, N0, tau_hat = tau), alpha)
  # p-values are clipped to [0, 1].
  expect_gte(evt_pvalue(-50, N1, N0), 0)
  expect_lte(evt_pvalue(50, N1, N0), 1)
})

test_that("EVT detection test holds its size under a white-noise null", {
  # Fixed-seed statistical acceptance: assemble the statistic inline (as the
  # reference seasEVT does) and check the empirical rejection rate at
  # alpha = 0.05 lands in [0.02, 0.09].
  reps  <- 200L
  n     <- 80L
  alpha <- 0.05

  g     <- fourier_grid(n)
  part  <- mbin_partition(5, 1, 4, g$omega_pos, band = "full",
                          exclude_nyquist = TRUE)
  isx   <- index_sets(part, g)

  set.seed(123)
  rej <- 0L
  for (r in seq_len(reps)) {
    x     <- stats::rnorm(n)
    pg    <- fs_periodogram(x, demean = TRUE, standardize = TRUE)
    I_pos <- pg$pgram_std[g$pos_idx]
    Delta <- max(I_pos[isx$J1]) - max(I_pos[isx$J0])
    pv    <- evt_pvalue(Delta, isx$N1, isx$N0, tau_hat = pg$tau_hat,
                        alternative = "greater")
    if (pv < alpha) rej <- rej + 1L
  }
  rate <- rej / reps

  expect_gte(rate, 0.02)
  expect_lte(rate, 0.09)
})
