# Tier-1 closed-form tests for the frequency-grid and bin machinery.

test_that("fourier_grid includes Nyquist for even n and drops it for odd n", {
  g <- fourier_grid(80)
  expect_equal(g$n, 80L)
  expect_equal(g$n_spec, 40L)
  expect_length(g$omega_pos, 40L)
  # Nyquist present: last positive ordinate is exactly pi.
  expect_equal(g$omega_pos[g$n_spec], pi)
  expect_equal(g$idx_nyq, 41L)               # full-DFT position n/2 + 1
  expect_equal(g$pos_idx, 2:41)
  expect_equal(g$omega_pos, 2 * pi * (1:40) / 80)

  go <- fourier_grid(81)
  expect_equal(go$n_spec, 40L)               # floor(81/2)
  expect_true(is.na(go$idx_nyq))             # no Nyquist for odd n
  expect_true(max(go$omega_pos) < pi)
  expect_true(all(go$omega_pos > 0 & go$omega_pos < pi))

  # Small even n: ordinates land on {pi/4, pi/2, 3pi/4, pi}.
  g8 <- fourier_grid(8)
  expect_equal(g8$omega_pos, c(pi/4, pi/2, 3*pi/4, pi))
  expect_equal(g8$idx_nyq, 5L)

  expect_error(fourier_grid(7), "single integer")
})

test_that("omega_seasonal(P=1, N=4) yields {pi/2, pi}, and {pi/2} without Nyquist", {
  expect_equal(omega_seasonal(1, 4), c(pi/2, pi))            # default includes Nyquist
  expect_equal(omega_seasonal(1, 4, exclude_nyquist = TRUE), pi/2)
})

test_that("mbin_partition M=5 gives seasonal bins {3,5}; M=4 warns on boundary", {
  g <- fourier_grid(80)

  # M = 5: pi/2 -> bin 3, pi -> bin 5, neither on an interior break -> no warning.
  part5 <- expect_no_warning(mbin_partition(5, 1, 4, g$omega_pos))
  expect_equal(part5$M1, c(3L, 5L))
  expect_equal(part5$P, 1L)
  expect_equal(part5$N, 4L)

  # M = 4: pi/2 lands exactly on the interior break pi/2 -> boundary warning.
  expect_warning(mbin_partition(4, 1, 4, g$omega_pos), "interior bin boundary")
})

test_that("index_sets isolates exact harmonics, shoulders, and mirrored full-DFT set", {
  g    <- fourier_grid(80)
  part <- mbin_partition(5, 1, 4, g$omega_pos)   # seasonal bins {3,5}
  isx  <- index_sets(part, g)

  # On-grid harmonics for n = 80: omega = pi/2 at j = 20, omega = pi at j = 40.
  expect_equal(isx$H, c(20L, 40L))

  # Shoulders are the seasonal band minus the exact harmonics.
  expect_equal(isx$S, setdiff(isx$J1, isx$H))
  expect_false(any(c(20L, 40L) %in% isx$S))

  # Cardinalities are consistent.
  expect_equal(isx$N1, length(isx$J1))
  expect_equal(isx$N0, length(isx$J0))

  # J1_full carries each seasonal ordinate AND its conjugate mirror (k <-> n-k+2).
  n <- g$n
  for (j in isx$J1) {
    k <- g$pos_idx[j]
    expect_true(k %in% isx$J1_full)
    expect_true((n - k + 2L) %in% isx$J1_full)
  }
  expect_equal(isx$J1_full, sort(unique(isx$J1_full)))
  # The Nyquist ordinate (full-DFT position 41) is self-conjugate: appears once.
  expect_equal(sum(isx$J1_full == 41L), 1L)

  # J0_full is positive-side only (no negative mirrors): all <= n_spec + 1.
  expect_true(all(isx$J0_full <= g$n_spec + 1L))
  expect_equal(isx$J0_full, g$pos_idx[isx$J0])
})

test_that("suggest_M reproduces the verified offset-free and offset selections", {
  # Offset-free default: M = 5, seasonal-bin share 0.2, no offset.
  sm0 <- suggest_M(P = 1, N = 4, n_obs = 60)
  expect_equal(sm0$M, 5)
  expect_equal(sm0$rho_M, 0.2)
  expect_equal(sm0$offset_u, 0)
  expect_equal(sm0$alpha, 0)
  expect_equal(sm0$H, 1L)

  # Shipped reference offset grid can realize M = 4 with a half-bin offset
  # (documented only; not honored by the package's offset-free binning).
  sm1 <- suggest_M(P = 1, N = 4, n_obs = 60, offset_grid = c(0, 0.5))
  expect_equal(sm1$M, 4)
  expect_equal(sm1$offset_u, 0.5)
  expect_equal(sm1$rho_M, 0.25)
})
