# whiten.R -- de-biased-Whittle whitening and exact recoloring.
#
# The design guarantee of this module: the fitted "nonseasonal" AR cannot
# absorb a concentrated seasonal comb, no matter its order. Protection comes
# from the OBJECTIVE (a Whittle likelihood evaluated only on a seasonality-
# blind ordinate set), not from any order restriction. At quarterly frequency
# an AR(2) with negative lag-2 coefficient peaks at pi/2 and an AR(1) with
# negative coefficient peaks at pi, so capping the order is not protection --
# restricting the ordinate set the likelihood sees is.
#
# Two refinements over the original M0-bin design (both measured in
# tests/testthat/test-attribution.R):
#   1. The model in the Whittle objective is the DE-BIASED expected periodogram
#      (Sykulski et al. 2019) -- the finite-n expectation E[I(w)], not the
#      infinite-sample AR spectral density -- so the fit no longer under-
#      extrapolates a peak's height.
#   2. The default exclusion is "guard" (only the harmonic tips), not "bins"
#      (whole seasonal bins). Guarding the tips keeps a genuine comb out of the
#      fit while leaving a broad nonseasonal peak's shoulders visible, which
#      resolves the size/power tension between the two pure designs.
# See notes/2026-07-22_ssi_phase_surgery.md sec.1 and the architecture note sec.5/8.
#
# Base R + stats only. The M0 ordinate geometry is delegated to R/partition.R
# (fourier_grid / mbin_partition / index_sets), so partition.R must be on the
# search path before this file is used (guaranteed inside the package
# namespace, and by alphabetical source() order for the test runner).

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Durbin-Levinson: partial autocorrelations (reflection coefficients) in (-1, 1)
# to AR coefficients in the (1 - sum phi_k B^k) convention. |kappa_k| < 1 for all
# k guarantees a stationary AR (roots strictly outside the unit circle).
.fs_pacf_to_ar <- function(kappa) {
  p <- length(kappa)
  if (p == 0L) return(numeric(0))
  phi <- kappa[1L]
  if (p >= 2L) {
    for (m in 2:p) {
      prev   <- phi
      phi    <- numeric(m)
      phi[m] <- kappa[m]
      for (i in seq_len(m - 1L)) phi[i] <- prev[i] - kappa[m] * prev[m - i]
    }
  }
  phi
}

# AR gain g(w; phi) = 1 / |1 - sum_k phi_k e^{-i k w}|^2 at frequencies w.
# (The infinite-sample AR spectral shape; retained as a general utility and used
# by the test-side unrestricted-Whittle contrast. The whitening objective itself
# now uses the de-biased expected periodogram, not this gain -- see below.)
.fs_ar_gain <- function(phi, w) {
  val <- complex(length.out = length(w), real = 1, imaginary = 0)
  for (k in seq_along(phi)) val <- val - phi[k] * exp(-1i * k * w)
  1 / (Mod(val)^2)
}

# De-biased Whittle expected periodogram (Sykulski, Olhede, Guillaumin, Lilly &
# Early 2019, Biometrika). For an AR(p) with coefficients `phi` and UNIT process
# variance (autocorrelations rho(k) from ARMAacf), the expected periodogram of a
# length-n_y sample is
#     Ebar(w_j) = sum_{|k| < n_y} (1 - |k|/n_y) rho(k) e^{-i k w_j},
# the DFT of the triangularly (Bartlett-) weighted acvf. At the Fourier
# frequencies w_j = 2*pi*j/n_y this is EXACT and n_y-periodic in k, so the whole
# two-sided weighted acvf folds into a length-n_y real vector whose FFT gives
# Ebar at every Fourier ordinate j = 0..n_y-1 in one pass. Ebar is the finite-n
# expectation E[I(w_j)] (unlike the AR gain, which is its n -> inf limit); using
# it as the model in the Whittle objective removes the finite-sample blur/bias
# that made a bin-restricted fit under-extrapolate a peak's height. mean(Ebar)
# over all j equals rho(0) = 1 by Parseval (a useful invariant).
.fs_expected_pgram <- function(phi, n_y) {
  if (length(phi) == 0L) return(rep.int(1, n_y))     # white: flat expected pgram
  rho   <- stats::ARMAacf(ar = phi, lag.max = n_y - 1L)   # lags 0..n_y-1, rho[1]=1
  m     <- seq_len(n_y - 1L)
  gfold <- numeric(n_y)
  gfold[1L]     <- rho[1L]                            # k = 0
  gfold[m + 1L] <- (1 - m / n_y) * rho[m + 1L] +      # positive lags
                   (m / n_y)     * rho[n_y - m + 1L]  # folded negative lags
  Re(stats::fft(gfold))
}

# Concentrated de-biased Whittle objective on the included ordinate set. `I_fit`
# are the periodogram ordinates the fit sees; `jfit` are their harmonic indices
# (into 1..n_spec, i.e. w = 2*pi*jfit/n_y), so the expected periodogram at those
# ordinates is `.fs_expected_pgram(phi, n_y)[jfit + 1L]`. The single scale (the
# process variance) is concentrated out exactly as before: s2 = mean(I / Ebar),
# profiled objective = length * log(s2) + sum(log(Ebar)).
.fs_whittle_Q <- function(phi, I_fit, jfit, n_y) {
  Ebar <- .fs_expected_pgram(phi, n_y)[jfit + 1L]
  s2   <- mean(I_fit / Ebar)
  length(I_fit) * log(s2) + sum(log(Ebar))
}

# Fit a stationary AR(p) by minimizing the concentrated de-biased Whittle
# objective over the tanh-PACF parameterization. The objective is multimodal, so
# optim is run from a deterministic set of starts -- the origin, an optional
# caller-supplied warm start (typically the order p-1 optimum, which yields a
# parsimony-preferring, monotone-Q path), and the sign-pattern grid
# {-1.5, +1.5}^p -- and the global best is kept. A single cold start lands in
# poor local minima and produces erratic, non-monotone BIC paths; the multistart
# recovers the true optimum reliably. Deterministic (no RNG): reproducible
# independent of the caller's random stream. p == 0 is the empty (white) model.
# `sigma2` is returned as the INNOVATION variance: the concentrated scale s2
# estimates the process variance (Ebar is normalized to rho(0) = 1), and
# innovation variance = process variance * prod(1 - kappa_k^2) by Levinson-Durbin.
.fs_fit_ar_p <- function(p, I_fit, jfit, n_y, warm = NULL) {
  if (p == 0L) {
    s2 <- mean(I_fit)
    return(list(phi = numeric(0), Q = length(I_fit) * log(s2), sigma2 = s2,
                par = numeric(0)))
  }
  obj <- function(theta) {
    # tanh maps R -> (-1, 1) so every candidate is stationary; the clamp is a
    # numerical margin keeping roots strictly outside the unit circle even if
    # optim wanders toward |kappa| = 1.
    kappa <- pmax(pmin(tanh(theta), 1 - 1e-6), -(1 - 1e-6))
    .fs_whittle_Q(.fs_pacf_to_ar(kappa), I_fit, jfit, n_y)
  }
  grid   <- as.matrix(expand.grid(rep(list(c(-1.5, 1.5)), p)))
  starts <- c(list(rep(0, p)),
              lapply(seq_len(nrow(grid)), function(i) as.numeric(grid[i, ])))
  if (!is.null(warm)) starts <- c(list(warm), starts)

  best_Q <- Inf; best_par <- rep(0, p)
  for (s in starts) {
    fit <- tryCatch(
      stats::optim(s, obj, method = "BFGS",
                   control = list(reltol = 1e-11, maxit = 500L)),
      error = function(e) NULL
    )
    if (!is.null(fit) && is.finite(fit$value) && fit$value < best_Q) {
      best_Q   <- fit$value
      best_par <- fit$par
    }
  }
  kappa  <- pmax(pmin(tanh(best_par), 1 - 1e-6), -(1 - 1e-6))
  phi    <- .fs_pacf_to_ar(kappa)
  Ebar   <- .fs_expected_pgram(phi, n_y)[jfit + 1L]
  s2     <- mean(I_fit / Ebar)               # concentrated process variance
  sigma2 <- s2 * prod(1 - kappa^2)           # innovation variance (Levinson-Durbin)
  list(phi = phi, Q = best_Q, sigma2 = sigma2, par = best_par)
}

# Guard-mode included-ordinate set: every positive Fourier ordinate EXCEPT those
# within `guard` Fourier steps of an exact seasonal harmonic. This excludes only
# the comb's concentrated mass (the harmonic tips) while leaving the broad
# shoulders of a genuine nonseasonal peak at a seasonal frequency visible to the
# fit -- so the whitener can explain such a peak (curing false alarms) yet stays
# blind to a concentrated comb (which lives inside the guarded neighborhoods).
.fs_guard_include <- function(grid, P, N, guard) {
  omega_h <- omega_seasonal(P, N, exclude_nyquist = FALSE)
  step    <- 2 * pi / grid$n
  tol     <- step * 1e-6
  excl    <- logical(grid$n_spec)
  for (wh in omega_h) {
    excl <- excl | (abs(grid$omega_pos - wh) <= guard * step + tol)
  }
  which(!excl)
}

# Conditional (time-domain) AR residuals e_t = z_t - sum_k phi_k z_{t-k},
# for t = (p+1)..length(z). Length length(z) - p. No circular whitening.
.fs_ar_residuals <- function(z, phi) {
  p  <- length(phi)
  nz <- length(z)
  if (nz <= p) return(numeric(0))
  e <- numeric(nz - p)
  for (t in (p + 1L):nz) {
    acc <- z[t]
    for (k in seq_len(p)) acc <- acc - phi[k] * z[t - k]
    e[t - p] <- acc
  }
  e
}

# ---------------------------------------------------------------------------
# fs_whiten
# ---------------------------------------------------------------------------

#' De-biased-Whittle whitening of a possibly-seasonal series
#'
#' Fits a small, provably-nonseasonal AR filter and returns the whitened
#' residuals together with everything \code{fs_recolor()} needs to invert the
#' transform exactly. The AR order and coefficients are chosen by a Whittle
#' likelihood evaluated \emph{only} on nonseasonal-bin Fourier ordinates (the
#' M0 set), so the filter cannot absorb energy at seasonal frequencies
#' regardless of the order selected. This is a deliberate design guarantee: at
#' small seasonal periods an order restriction alone is not seasonality
#' blindness (an AR(2) with negative lag-2 coefficient peaks at pi/2), whereas
#' restricting the likelihood's ordinate set is.
#'
#' The returned residual vector is trimmed from the front to a whole number of
#' seasonal cycles, so \code{n_e} is a multiple of \code{N} rather than
#' \code{length(x) - d - p}. Every downstream stage works on the Fourier grid of
#' that vector, and the exact seasonal harmonics are ordinates of that grid only
#' under this alignment. The trimmed head is returned in \code{e_head} and
#' spliced back by \code{fs_recolor()}, so the round trip stays exact; the cost
#' is that the leading \code{n_trim + d + p} observations of an adjusted series
#' are left unadjusted.
#'
#' @param x Numeric vector, already extracted from any \code{ts}/\code{tsibble}
#'   upstream. Must be finite, non-constant, and of length at least \code{3 * N}.
#' @param N Positive integer seasonal period (observations per year, e.g. 4 for
#'   quarterly, 12 for monthly).
#' @param P Positive integer fundamental-per-year multiplier (default \code{1L}).
#' @param M Integer (>= 2) bin count for the equal-width partition of \code{(0, pi]}
#'   used to build the nonseasonal (M0) ordinate set for the Whittle fit. The
#'   partition is built on the \emph{differenced, demeaned} series' Fourier grid.
#' @param d Differencing rule. One of \code{"auto"} (difference once iff the
#'   lag-1 sample autocorrelation of \code{x} exceeds 0.95), \code{"none"}
#'   (force \code{d = 0}), or \code{"first"} (force \code{d = 1}).
#' @param ar_max Non-negative integer, maximum AR order considered by the
#'   de-biased-Whittle BIC search over \code{p = 0, 1, ..., ar_max} (default
#'   \code{3L}).
#' @param exclusion Which ordinates the whitening fit is blinded to. One of
#'   \code{"guard"} (default) or \code{"bins"}. \code{"guard"} excludes only the
#'   ordinates within \code{guard} Fourier steps of each exact seasonal harmonic
#'   (the comb's concentrated mass), leaving the broad shoulders of a nonseasonal
#'   peak at a seasonal frequency visible to the fit. \code{"bins"} excludes all
#'   seasonal-\emph{bin} ordinates (the historical behavior). See Details.
#' @param guard Non-negative integer; for \code{exclusion = "guard"}, the
#'   half-width (in Fourier steps) of the excluded neighborhood around each exact
#'   seasonal harmonic. Default \code{3L} (the value selected by the experiment
#'   below). Ignored when \code{exclusion = "bins"}.
#'
#' @details
#' The whitener's blindness to seasonality comes from the ordinate set its
#' de-biased Whittle likelihood evaluates on, not from any order restriction.
#' The two exclusion geometries trade off in opposite directions. The default
#' (\code{exclusion = "guard"}, \code{guard = 3}) was chosen by a fixed-seed
#' experiment (quarterly, n = 120, 100 reps; see
#' \code{tests/testthat/test-attribution.R}) minimizing the case-A false-alarm
#' rate subject to power > 0.95 on a true seasonal band and comb-absorption
#' bin-ratio > 10. The measured matrix, at \eqn{\alpha = 0.05}:
#' \tabular{lllll}{
#'   \strong{config} \tab \strong{A1 rej} \tab \strong{A2 rej} \tab \strong{B power} \tab \strong{C bin-ratio} \cr
#'   \code{bins}     \tab 0.56 \tab 0.61 \tab 1.00 \tab 21.0 \cr
#'   \code{guard(2)} \tab 0.07 \tab 0.30 \tab 0.93 \tab 16.1 \cr
#'   \code{guard(3)} \tab 0.11 \tab 0.44 \tab 0.98 \tab 18.8 \cr
#' }
#' where A1/A2 are pure nonseasonal AR(2) cycles at \eqn{\pi/2}
#' (\eqn{\rho = 0.85 / 0.95}; correct answer ~0.05), B is a true AR(2)\eqn{\times}SAR(1)
#' band (correct answer ~1.0), and C is comb absorption (a genuine comb must
#' stay in the residual, ratio high). \code{guard(2)} gives the best size, but its
#' band power is borderline at 100 reps (0.93) and re-estimating at 400 reps puts
#' it at 0.91 -- robustly below the > 0.95 constraint -- so it is disqualified.
#' \code{guard(3)} is selected: it cuts the pure-cycle false alarm sharply versus
#' \code{bins} (A1 0.56 -> 0.11) while keeping power (0.965 at 400 reps) and comb
#' protection. The narrow \eqn{\rho = 0.95} cycle (A2) remains the hard case
#' (0.44): its peak is barely wider than the guard, so some residual mass
#' survives. Use \code{exclusion = "bins"} to recover the historical geometry.
#'
#' @return A list with fields:
#'   \describe{
#'     \item{d}{Integer differencing order applied (0 or 1).}
#'     \item{d_reason}{Human-readable string recording why \code{d} was chosen
#'       (e.g. \code{"lag-1 ACF 0.97 > 0.95 (auto: differenced once)"}).}
#'     \item{mu}{Numeric mean removed from the (differenced) series.}
#'     \item{p}{Selected AR order.}
#'     \item{ar}{Numeric length-\code{p} vector of AR coefficients in the
#'       \code{(1 - sum phi_k B^k)} convention (\code{numeric(0)} if \code{p = 0}).}
#'     \item{sigma2}{Innovation variance at the selected order (the concentrated
#'       de-biased-Whittle process-variance scale times
#'       \eqn{\prod_k (1 - \kappa_k^2)}, Levinson-Durbin).}
#'     \item{bic_path}{Named numeric vector of de-biased-Whittle BIC values, one per candidate
#'       order \code{0..ar_max} (names \code{"0".."ar_max"}).}
#'     \item{e}{Numeric vector of conditional AR residuals (the whitened series)
#'       on the seasonally aligned grid, length \code{n_e}.}
#'     \item{n_e}{Integer length of \code{e}: the largest multiple of \code{N}
#'       not exceeding \code{length(y0) - p}. Aligning the residual grid to a
#'       whole number of seasonal cycles is what makes the exact seasonal
#'       harmonics \code{2 * pi * h / N} Fourier ordinates of the grid every
#'       downstream stage works on.}
#'     \item{e_head}{The \code{n_trim} leading conditional residuals held out of
#'       the aligned grid (\code{numeric(0)} when the grid was already aligned).
#'       \code{fs_recolor()} splices these back, so the round trip stays exact.}
#'     \item{n_trim}{Integer number of leading residuals trimmed, in
#'       \code{0:(N - 1)}. The first \code{n_trim + d + p} observations of an
#'       adjusted series are therefore left unadjusted.}
#'     \item{anchors}{The first \code{d + p} values of \code{x} (level anchors
#'       that \code{fs_recolor()} reproduces exactly).}
#'     \item{N, P, M}{The seasonal period, multiplier, and bin count, echoed for
#'       downstream use.}
#'     \item{exclusion, guard}{The exclusion mode and guard half-width used for
#'       the fit's ordinate set, echoed for the post-adjustment self-check.}
#'   }
#'
#' @keywords internal
fs_whiten <- function(x, N, P = 1L, M, d = c("auto", "none", "first"),
                      ar_max = 3L, exclusion = c("guard", "bins"), guard = 3L) {
  d_arg     <- match.arg(d)
  exclusion <- match.arg(exclusion)

  # --- validation --------------------------------------------------------
  x <- as.numeric(x)
  if (length(x) < 1L || any(!is.finite(x))) {
    stop("`x` must be a finite numeric vector.", call. = FALSE)
  }
  if (length(N) != 1L || !is.finite(N) || N < 2L) {
    stop("`N` must be a single integer >= 2.", call. = FALSE)
  }
  N <- as.integer(N)
  if (length(P) != 1L || !is.finite(P) || P < 1L) {
    stop("`P` must be a single positive integer.", call. = FALSE)
  }
  P <- as.integer(P)
  if (length(M) != 1L || !is.finite(M) || M < 2L) {
    stop("`M` must be a single integer >= 2.", call. = FALSE)
  }
  M <- as.integer(M)
  if (length(ar_max) != 1L || !is.finite(ar_max) || ar_max < 0L) {
    stop("`ar_max` must be a single non-negative integer.", call. = FALSE)
  }
  ar_max <- as.integer(ar_max)
  if (length(guard) != 1L || !is.finite(guard) || guard < 0L) {
    stop("`guard` must be a single non-negative integer.", call. = FALSE)
  }
  guard <- as.integer(guard)

  if (length(x) < 3L * N) {
    stop(sprintf(
      "Series too short: length(x) = %d < 3 * N = %d. Need at least three full years to fit the de-biased-Whittle whitener.",
      length(x), 3L * N), call. = FALSE)
  }
  if (stats::sd(x) == 0) {
    stop("`x` is constant (zero variance); nothing to whiten.", call. = FALSE)
  }

  # --- differencing decision --------------------------------------------
  r1 <- stats::acf(x, lag.max = 1L, plot = FALSE, demean = TRUE)$acf[2L]
  if (d_arg == "auto") {
    if (is.finite(r1) && r1 > 0.95) {
      d <- 1L
      d_reason <- sprintf("lag-1 ACF %.2f > 0.95 (auto: differenced once)", r1)
    } else {
      d <- 0L
      d_reason <- sprintf("lag-1 ACF %.2f <= 0.95 (auto: no differencing)", r1)
    }
  } else if (d_arg == "first") {
    d <- 1L
    d_reason <- "d = 'first' (one difference forced)"
  } else {
    d <- 0L
    d_reason <- "d = 'none' (differencing suppressed)"
  }

  y  <- if (d == 1L) diff(x) else x
  mu <- mean(y)
  y0 <- y - mu
  n_y <- length(y0)

  if (stats::sd(y0) < 1e-12 * (1 + abs(mu))) {
    stop("(Differenced) series has no variation; nothing to whiten.",
         call. = FALSE)
  }

  # --- M0 ordinate set and periodogram of y0 ----------------------------
  # Partition machinery unified with R/partition.R (wave 2). The M0 (donor)
  # ordinate set is the nonseasonal-bin positive Fourier ordinates. Nyquist is
  # treated as seasonal (exclude_nyquist = FALSE) so its quarterly peak at
  # omega = pi is kept out of the fit. fourier_grid()'s positive grid includes
  # the Nyquist ordinate for even n_y, but that ordinate lands in a seasonal
  # bin and is therefore excluded from J0 -- so the resulting (I0, w0) are
  # identical to the pre-unification duplicated-logic set.
  grid  <- fourier_grid(n_y)
  part  <- mbin_partition(M, P, N, grid$omega_pos, band = "full",
                          exclude_nyquist = FALSE)
  isx   <- index_sets(part, grid)
  pgram <- Mod(stats::fft(y0))^2 / n_y            # positions 1..n_y (DC at 1)
  I_pos <- pgram[grid$pos_idx]                    # positive-frequency ordinates

  # --- included ordinate set for the de-biased Whittle fit --------------
  # "bins":  exclude all seasonal-BIN ordinates (the M0/donor set). Wide bins
  #          hide a genuine nonseasonal peak's shoulders, so the fit under-
  #          extrapolates its height.
  # "guard": exclude only ordinates within `guard` Fourier steps of an exact
  #          seasonal harmonic (the comb's concentrated mass), leaving the broad
  #          peak shoulders visible to the fit.
  jfit <- if (exclusion == "bins") isx$J0 else .fs_guard_include(grid, P, N, guard)
  if (length(jfit) < 1L) {
    stop("The whitening fit's ordinate set is empty (all positive ordinates ",
         "were excluded). Reduce `guard`, widen `M`, or use a longer series.",
         call. = FALSE)
  }
  I_fit <- I_pos[jfit]

  # --- de-biased-Whittle BIC search over p = 0..ar_max ------------------
  fits     <- vector("list", ar_max + 1L)
  bic_path <- numeric(ar_max + 1L)
  warm     <- NULL
  for (pp in 0:ar_max) {
    f <- .fs_fit_ar_p(pp, I_fit, jfit, n_y, warm = warm)
    fits[[pp + 1L]]    <- f
    bic_path[pp + 1L]  <- f$Q + pp * log(length(I_fit))
    warm <- c(f$par, 0)                 # warm-start the next order at this optimum
  }
  names(bic_path) <- as.character(0:ar_max)

  p_sel <- unname(which.min(bic_path)) - 1L    # bare integer (bic_path is named)
  best  <- fits[[p_sel + 1L]]
  phi   <- best$phi

  # --- conditional residuals (time domain, no circular whitening) -------
  e <- .fs_ar_residuals(y0, phi)

  # --- seasonal alignment of the residual (surgery) grid ----------------
  # Everything downstream -- fourier_grid(), mbin_partition(), index_sets(),
  # the periodogram and fs_surgery() -- lives on the length-n_e Fourier grid,
  # and the exact seasonal harmonics 2*pi*h/N are Fourier ordinates of that
  # grid only when N divides n_e. Differencing and the AR fit consume d + p
  # observations, so an unaligned grid was the generic case: the exact-harmonic
  # set H collapsed to at most {Nyquist}, a deterministic seasonal line leaked
  # across the spectrum (spuriously rejecting the shoulder test and so being
  # misclassified as a band), and bin-local gain surgery could not remove the
  # leaked line. Trimming to a whole number of seasonal cycles restores the
  # geometry the rest of the pipeline assumes. The trim is taken from the FRONT:
  # the recent end of the series is the end that matters for seasonal
  # adjustment, so the held-out head is the oldest (< N) residuals. They are
  # returned in `e_head` and spliced back by fs_recolor(), which keeps the round
  # trip exact -- the cost is that the first n_trim + d + p observations of an
  # adjusted series are left unadjusted.
  n_full <- length(e)
  n_use  <- N * (n_full %/% N)      # largest whole number of seasonal cycles
  n_trim <- n_full - n_use          # leading residuals held out of the grid (< N)
  if (n_use < 2L * N) {
    stop(sprintf(
      paste0("Series too short once the residual grid is seasonally aligned: ",
             "%d residuals (length(x) = %d, d = %d, AR order p = %d) trim to ",
             "%d, below the required 2 * N = %d. Supply a longer series."),
      n_full, length(x), d, p_sel, n_use, 2L * N), call. = FALSE)
  }
  e_head <- e[seq_len(n_trim)]      # untouched head, spliced back by fs_recolor()
  e      <- e[(n_trim + 1L):n_full]

  list(
    d        = d,
    d_reason = d_reason,
    mu       = mu,
    p        = p_sel,
    ar       = phi,
    sigma2   = best$sigma2,
    bic_path = bic_path,
    e         = e,
    n_e       = n_use,
    e_head    = e_head,
    n_trim    = n_trim,
    anchors   = x[seq_len(d + p_sel)],
    N         = N,
    P         = P,
    M         = M,
    exclusion = exclusion,
    guard     = guard
  )
}

# ---------------------------------------------------------------------------
# fs_recolor
# ---------------------------------------------------------------------------

#' Invert de-biased-Whittle whitening exactly
#'
#' Reconstructs a full-length series from (possibly adjusted) innovations by
#' running the AR recursion forward and undoing the differencing, using the
#' whitening record from \code{fs_whiten()}. The first \code{d + p} values are
#' held at the original level anchors, so the reconstruction reproduces the
#' input exactly at the endpoints regardless of the supplied innovations. When
#' \code{estar} is the stored \code{wh$e}, the round trip reproduces \code{x}
#' to numerical precision.
#'
#' In practice the leading run of unchanged values is longer than the \code{d + p}
#' anchors: the \code{wh$n_trim} residuals that \code{fs_whiten()} held out to
#' align the grid to a whole number of seasonal cycles are spliced back
#' unmodified from \code{wh$e_head}, so the first \code{wh$n_trim + d + p} values
#' of the reconstruction equal those of \code{x} whatever \code{estar} contains.
#' The guaranteed anchor convention is still \code{d + p}; the extra
#' \code{wh$n_trim} values are unchanged because that head is never surgered.
#'
#' @param estar Numeric vector of (adjusted) innovations on the seasonally
#'   aligned residual grid, length \code{wh$n_e}, aligned to residual times
#'   \code{(p + wh$n_trim + 1)..length(y0)}. The \code{wh$n_trim} leading
#'   residuals that \code{fs_whiten()} held out of that grid are taken
#'   unchanged from \code{wh$e_head} and prepended here, so the AR recursion
#'   still runs over all \code{length(y0) - p} innovations.
#' @param wh The list returned by \code{fs_whiten()} (supplies \code{d},
#'   \code{mu}, \code{p}, \code{ar}).
#' @param x The original numeric input to \code{fs_whiten()}; supplies the
#'   demeaned initial values and the level anchor \code{x[1]} for
#'   de-differencing.
#'
#' @return A numeric vector the same length as \code{x}. Its first
#'   \code{d + p} values equal those of \code{x} — the guaranteed anchor
#'   convention, which holds for any \code{estar}. Because the
#'   \code{wh$n_trim} leading residuals held out of the seasonally aligned
#'   grid are spliced back unmodified, the first \code{wh$n_trim + d + p}
#'   values in fact equal those of \code{x} as well; that longer run is a
#'   consequence of the trim (the head is never surgered), not part of the
#'   anchor guarantee.
#'
#' @keywords internal
fs_recolor <- function(estar, wh, x) {
  x     <- as.numeric(x)
  estar <- as.numeric(estar)
  d  <- wh$d
  p  <- wh$p
  mu <- wh$mu
  ar <- wh$ar

  y   <- if (d == 1L) diff(x) else x
  y0  <- y - mu
  n_y <- length(y0)

  # `estar` lives on the seasonally aligned residual grid; the leading residuals
  # fs_whiten() trimmed to reach that alignment are spliced back unchanged so
  # the recursion below sees all n_y - p innovations in their original order.
  # `e_head` is NULL on hand-built/legacy whitener records: treat as length 0.
  e_head <- if (is.null(wh$e_head)) numeric(0) else as.numeric(wh$e_head)
  n_e    <- if (is.null(wh$n_e)) n_y - p - length(e_head) else wh$n_e

  if (length(estar) != n_e) {
    stop(sprintf(
      paste0("`estar` has length %d but %d innovations are expected on the ",
             "seasonally aligned residual grid (%d leading residual(s) were ",
             "trimmed and are supplied by `wh$e_head`)."),
      length(estar), n_e, length(e_head)), call. = FALSE)
  }
  estar <- c(e_head, estar)

  ystar <- numeric(n_y)
  if (p > 0L) ystar[seq_len(p)] <- y0[seq_len(p)]
  if (n_y > p) {
    for (t in (p + 1L):n_y) {
      acc <- estar[t - p]
      for (k in seq_len(p)) acc <- acc + ar[k] * ystar[t - k]
      ystar[t] <- acc
    }
  }

  y_rec <- ystar + mu
  if (d == 1L) {
    stats::diffinv(y_rec, differences = 1L, xi = x[1L])
  } else {
    y_rec
  }
}
