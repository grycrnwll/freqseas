# test-gain.R
# Tests for fs_donor_pool() and fs_gain() (internal, R/gain.R) against a
# fixture periodogram with known elevated J1 ordinates.

test_that("fs_donor_pool selects J0 ordinates at/below the donor quantile", {
  Jmax <- 20L
  set.seed(1)
  pgram_pos <- rep(1, Jmax)
  J1 <- c(9L, 10L, 11L)
  J0 <- setdiff(1:Jmax, J1)

  pool <- fs_donor_pool(pgram_pos, J0, donor_quantile = 0.9)
  expect_true(all(pool %in% J0))
  expect_gte(length(pool), 5L)

  # quantile filter leaves < 5 donors -> informative error
  J0b <- 1:6
  pg2 <- c(1, 2, 3, 4, 5, 6)
  expect_error(fs_donor_pool(pg2, J0b, donor_quantile = 0.5), "Too few donors")

  # too few nonseasonal ordinates to begin with -> informative error
  expect_error(fs_donor_pool(pgram_pos, J0 = 1:3), "Too few nonseasonal ordinates")
})

test_that("fs_gain: floor at 1, smoothing confined to J1 expanded by one, symmetric full circle", {
  n <- 41L; Jmax <- 20L
  J1 <- c(9L, 10L, 11L)
  J0 <- setdiff(1:Jmax, J1)

  pgram_pos <- rep(1, Jmax)          # flat donor level everywhere in J0
  pgram_pos[9]  <- 9
  pgram_pos[10] <- 25
  pgram_pos[11] <- 9

  pool <- fs_donor_pool(pgram_pos, J0, donor_quantile = 0.9)
  g <- fs_gain(pgram_pos, sets = list(J1 = J1), donor_pool = pool, n = n)

  expect_equal(g$donor_level, 1)     # every J0 value is exactly 1

  # never amplify: gain >= 1 everywhere
  expect_true(all(g$Ghat_pos >= 1 - 1e-12))
  expect_true(all(g$Ghat_full >= 1 - 1e-12))

  # == 1 (log-gain 0) outside the expanded bin {8, ..., 12}
  expanded <- 8:12
  outside <- setdiff(1:Jmax, expanded)
  expect_equal(g$Ghat_pos[outside], rep(1, length(outside)))

  # independent (loop-based) reference implementation of the same spec:
  # floor(pgram/donor, 1) -> log(gain) -> 3-pt moving average over J1
  # expanded by 1 each side -> identity elsewhere.
  G2_ref <- rep(1, Jmax); G2_ref[J1] <- pgram_pos[J1] / g$donor_level
  G2_ref <- pmax(G2_ref, 1)
  lg_ref <- log(G2_ref) / 2
  sm_ref <- numeric(Jmax)
  for (j in expanded) {
    lo <- max(1L, j - 1L); hi <- min(Jmax, j + 1L)
    sm_ref[j] <- mean(lg_ref[lo:hi])
  }
  expect_equal(g$Ghat_pos[expanded], exp(sm_ref[expanded]), tolerance = 1e-10)

  # full-circle vector: symmetric, DC = 1
  pos_idx <- (1:Jmax) + 1L
  neg_idx <- n + 1L - (1:Jmax)
  expect_length(g$Ghat_full, n)
  expect_equal(g$Ghat_full[pos_idx], g$Ghat_full[neg_idx])
  expect_equal(g$Ghat_full[1], 1)
  expect_equal(g$Ghat_full[pos_idx], g$Ghat_pos)
})

test_that("fs_gain errors informatively on structurally bad input", {
  Jmax <- 20L
  pgram_pos <- rep(1, Jmax)
  expect_error(fs_gain(pgram_pos, sets = list(J1 = integer(0)), donor_pool = 1:10, n = 41L),
               "J1.*empty")
  expect_error(fs_gain(pgram_pos, sets = list(J1 = c(5, 6)), donor_pool = integer(0), n = 41L),
               "donor_pool.*empty")
})
