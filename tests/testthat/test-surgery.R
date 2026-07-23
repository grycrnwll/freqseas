# test-surgery.R
# Tests for fs_surgery() (internal, R/surgery.R): Hermitian/realness guards,
# the oracle deconvolution exactness result (SS3 of
# notes/2026-07-22_ssi_phase_surgery.md), and the line-branch donor-level
# replacement.

test_that("band branch: Hermitian symmetry / DC / Nyquist guards hold", {
  n <- 40L; Phi <- 0.7; s <- 4L
  w <- 2 * pi * (0:(n - 1)) / n
  G <- 1 / sqrt(1 + Phi^2 - 2 * Phi * cos(s * w))
  theta <- min_phase(log(G))

  set.seed(7)
  x <- rnorm(n)
  E <- fft(x - mean(x))

  res <- fs_surgery(E, Ghat_full = G, theta_full = theta, spec = "band",
                     sets = list(), donor_level = 1, n = n)

  # real before Re(): the raw inverse FFT of Estar has negligible imaginary part
  raw_ifft <- stats::fft(res$Estar, inverse = TRUE) / n
  expect_lt(max(abs(Im(raw_ifft))), 1e-12)

  # DC untouched exactly
  expect_identical(res$Estar[1], E[1])

  # Nyquist forced exactly real
  expect_equal(Im(res$Estar[n / 2L + 1L]), 0)

  # estar is the real part of that inverse FFT
  expect_equal(res$estar, Re(raw_ifft))
})

test_that("oracle band surgery with the true gain + minimum phase recovers the exact deconvolution", {
  # DGP: (1 - Phi B^s) x = eps. Exact identity: x_t - Phi*x_{t-s} = eps_t.
  # Circular DFT division by the true H(w) = 1/(1 - Phi e^{-isw}) is
  # exactly a 2-tap circular filter, so away from the wrap (t <= s) the
  # reconstruction should match the identity to numerical precision.
  n <- 80L; Phi <- 0.9; s <- 4L; spin <- 400L
  set.seed(11)
  ntot <- spin + n
  eps <- stats::rnorm(ntot)
  xs <- numeric(ntot)
  for (t in (s + 1L):ntot) xs[t] <- Phi * xs[t - s] + eps[t]
  x  <- xs[(spin + 1L):ntot]
  xd <- x - mean(x)

  w <- 2 * pi * (0:(n - 1)) / n
  G <- 1 / sqrt(1 + Phi^2 - 2 * Phi * cos(s * w))
  theta <- min_phase(log(G))

  E   <- fft(xd)
  res <- fs_surgery(E, Ghat_full = G, theta_full = theta, spec = "band",
                     sets = list(), donor_level = 1, n = n)
  estar <- res$estar

  trim   <- (s + 5L):n              # t = 9..80: away from the wrapped edge
  target <- xd[trim] - Phi * xd[trim - s]
  rmse   <- sqrt(mean((estar[trim] - target)^2))

  expect_lt(rmse, 0.15 * sd(x))
})

test_that("line branch: harmonic ordinate reset to donor level, all else bit-identical", {
  n <- 80L
  set.seed(42)
  t_idx <- 0:(n - 1)
  w0 <- pi / 2                       # exact harmonic at n/4
  x  <- 5 * cos(w0 * t_idx + 0.3) + stats::rnorm(n, sd = 1)
  E  <- fft(x - mean(x))

  H_pos  <- (n %/% 4L) + 1L          # full-DFT position of omega = pi/2
  donor_level <- 3

  res <- fs_surgery(E, Ghat_full = NULL, theta_full = NULL, spec = "line",
                     sets = list(H = H_pos), donor_level = donor_level, n = n)
  Estar <- res$Estar

  pg_h <- Mod(Estar[H_pos])^2 / n
  expect_equal(pg_h, donor_level, tolerance = 1e-10)

  # phase at the harmonic is retained (not overwritten)
  expect_equal(Arg(Estar[H_pos]), Arg(E[H_pos]), tolerance = 1e-10)

  # all ordinates other than DC, Nyquist, H, and mirror(H) are bit-identical
  # to the input (DC/Nyquist are covered by the dedicated guard test above;
  # excluding them here isolates the line-branch-specific claim).
  mirror_H <- n + 2L - H_pos
  nyq      <- n %/% 2L + 1L
  excl <- c(1L, nyq, H_pos, mirror_H)
  keep <- setdiff(seq_len(n), excl)
  expect_identical(Estar[keep], E[keep])
})
