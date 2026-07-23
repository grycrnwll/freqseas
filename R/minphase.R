# minphase.R
# Cepstral minimum-phase reconstruction of a filter's phase from its gain
# alone. See notes/2026-07-22_ssi_phase_surgery.md (SS3-4) for the theory
# and the verified closed-form numbers this port reproduces, and
# scripts/mp_phase_mre.R in seasonality_vintages for the original working
# code (mp_phase()) this function ports exactly.

#' Minimum-phase reconstruction from a gain function
#'
#' @description
#' A gain function alone does not determine a filter. For any linear filter
#' \eqn{H(\omega) = G(\omega) e^{i\theta(\omega)}}, multiplying by an
#' all-pass factor (magnitude 1 everywhere, arbitrary phase) leaves the gain
#' \eqn{G(\omega) = |H(\omega)|} unchanged but changes the filter, hence
#' changes which time-domain sequence it produces. Estimating \eqn{G} from a
#' periodogram (as \code{fs_gain()} does) and then only replacing magnitudes
#' is therefore an underdetermined inverse problem: infinitely many filters
#' share that gain, and "keep the original phase" is merely one arbitrary
#' choice among them, not a principled one.
#'
#' \strong{Minimum phase} resolves the degeneracy by adding a filter-class
#' assumption. Among all filters with a given gain, exactly one is both
#' causal (its output at time \eqn{t} depends only on inputs at or before
#' \eqn{t}) and causally invertible (its inverse is also causal): the
#' minimum-phase filter. For a real, causal, minimum-phase filter the phase
#' is not a free parameter -- it is pinned to the gain by the discrete
#' Hilbert-transform relation \eqn{\theta(\omega) = -\mathcal{H}\{\log
#' G\}(\omega)}, i.e. phase is the harmonic conjugate of log gain. This is
#' the same Hardy-space relation that makes the real and imaginary parts of
#' a one-sided (causal) sequence's transform Hilbert-transform pairs of one
#' another; here it is applied to the cepstrum (the log-spectrum's own
#' Fourier transform) rather than to the sequence itself.
#'
#' \code{min_phase()} computes \eqn{\theta} by cepstral folding (Oppenheim &
#' Schafer's minimum-phase reconstruction): the real cepstrum of
#' \code{log_gain} is computed via an inverse FFT, and then \strong{folded}
#' -- the anticausal (negative-lag) half of the cepstrum is added onto the
#' causal (positive-lag) half, and the anticausal half is discarded (set to
#' zero). This fold is exactly where causality enters the computation. A
#' symmetric (two-sided) cepstrum is the cepstrum of a symmetric, zero-phase
#' filter: no phase distortion, but noncausal (it needs future values --
#' the classical X-11-style centered moving average, which revises as new
#' data extend the available future window). A one-sided cepstrum is the
#' cepstrum of a causal, minimum-phase filter: no future values are needed
#' (no such revisions), at the cost of the phase rotation this function
#' computes. The choice between them is exactly the revision-vs-timing
#' tradeoff familiar from seasonal adjustment practice.
#'
#' Because the fold only uses the evenness of \code{log_gain} (true of any
#' real filter's gain) and standard Fourier identities, this reconstruction
#' is not specific to any parametric model: it consumes any log-integrable
#' gain (the Wold/Szego condition), including gains that do not correspond
#' to a finite-order ARMA filter.
#'
#' @param log_gain Numeric vector of length \code{n}: the natural-log gain
#'   \eqn{\log G(\omega_j)} on the FULL DFT circle, ordinate
#'   \eqn{j = 0, \dots, n-1} in \code{stats::fft()}'s 1-indexed convention
#'   (\code{log_gain[1]} is \eqn{\omega = 0}, \code{log_gain[2]} is
#'   \eqn{\omega = 2\pi/n}, ..., \code{log_gain[n]} is
#'   \eqn{\omega = 2\pi(n-1)/n}). Gain is an even function of frequency for
#'   any real filter, so \code{log_gain} is required to be symmetric:
#'   \code{log_gain[j + 1] == log_gain[n - j + 1]} for every
#'   \code{j = 1, ..., n - 1} (equivalently, \code{log_gain[j] ==
#'   log_gain[n - j]} in the 0-indexed ordinate labeling used above). This
#'   is validated (to a numerical tolerance) and an informative error is
#'   raised on violation -- most commonly caused by passing a
#'   positive-frequencies-only vector instead of the full circle.
#' @param fine_factor Single positive integer, default \code{64L}.
#'   \code{log_gain} is first linearly interpolated (periodic: the circle is
#'   closed by appending the \eqn{j = 0} value at \eqn{\omega = 2\pi}) onto a
#'   finer circle of length \code{n * fine_factor}, the cepstral fold is
#'   performed on that fine grid, and the resulting phase is subsampled back
#'   to the original \code{n} ordinates (every \code{fine_factor}-th fine
#'   point). Interpolating before folding reduces cepstral aliasing for
#'   gains with sharp, narrow peaks (e.g. near-unit-root seasonal AR gains);
#'   it refines the numerical accuracy of the reconstruction, it does not
#'   change which filter is being reconstructed.
#'
#' @return Numeric vector of length \code{n}: the phase
#'   \eqn{\theta(\omega_j)} (radians) of the unique minimum-phase filter
#'   with gain \code{exp(log_gain)}, evaluated on the same ordinate grid as
#'   the input.
#'
#' @examples
#' # Reconstruct the phase of a seasonal AR(4)-type filter
#' # H(w) = 1 / (1 - Phi * exp(-4i * w)) from its gain alone -- the sawtooth
#' # phase this filter rotates through at each seasonal harmonic. The
#' # closed form is theta(w) = -atan2(Phi * sin(4w), 1 - Phi * cos(4w)).
#' Phi <- 0.9; n <- 80L
#' w  <- 2 * pi * (0:(n - 1)) / n
#' G  <- 1 / sqrt(1 + Phi^2 - 2 * Phi * cos(4 * w))
#' theta_true <- -atan2(Phi * sin(4 * w), 1 - Phi * cos(4 * w))
#' theta_hat  <- min_phase(log(G))
#' # At only n = 80 knots for a fairly peaked gain (Phi = 0.9), linear
#' # interpolation to the fine grid (see `fine_factor`) under-resolves the
#' # peak's curvature and the error saturates around 0.02 rad no matter how
#' # large `fine_factor` is set (verified in tests/testthat/test-minphase.R):
#' # the fold itself is exact to < 1e-10 given an accurate input (also
#' # tested), so this ~0.02 rad is a resolution-boundary artifact of the
#' # n = 80 grid, not the cepstral fold. It shrinks as n grows.
#' max(abs(theta_hat - theta_true))  # ~0.02 rad at n = 80
#'
#' @export
min_phase <- function(log_gain, fine_factor = 64L) {
  n <- length(log_gain)
  if (n < 2L) stop("`log_gain` must have length >= 2.")
  if (!is.numeric(log_gain) || anyNA(log_gain) || any(!is.finite(log_gain))) {
    stop("`log_gain` must be a finite numeric vector.")
  }
  if (length(fine_factor) != 1L || !is.finite(fine_factor) || fine_factor < 1L) {
    stop("`fine_factor` must be a single positive integer.")
  }
  fine_factor <- as.integer(round(fine_factor))

  # Symmetry check: log_gain[j] == log_gain[n - j] for j = 0..n-1 (0-indexed
  # ordinate labels); in R's 1-indexing the partner of position p is
  # position 1 for p = 1 (DC, self-symmetric), and n - p + 2 for p = 2..n.
  partner <- c(1L, rev(seq(2L, n)))
  if (!isTRUE(all.equal(log_gain, log_gain[partner], tolerance = 1e-8))) {
    stop(
      "`log_gain` is not symmetric: log_gain[j] must equal log_gain[n - j] ",
      "for j = 0..n-1 (a real filter's gain is even in frequency). This ",
      "usually means a positive-frequencies-only gain was passed instead ",
      "of the gain on the FULL DFT circle (length n, DC at position 1)."
    )
  }

  # ---- interpolate onto a fine circle grid (linear, periodic) ----
  Nf <- n * fine_factor
  wn <- 2 * pi * (0:(n - 1)) / n
  wf <- 2 * pi * (0:(Nf - 1)) / Nf
  lg_fine <- stats::approx(
    x    = c(wn, 2 * pi),
    y    = c(log_gain, log_gain[1]),
    xout = wf,
    rule = 2
  )$y

  # ---- cepstral fold: real cepstrum -> one-sided (causal) fold -> phase ----
  ce  <- Re(stats::fft(lg_fine, inverse = TRUE)) / Nf
  wgt <- numeric(Nf)
  wgt[1] <- 1
  if (Nf %% 2L == 0L) {
    wgt[2:(Nf / 2L)]  <- 2
    wgt[Nf / 2L + 1L] <- 1
  } else {
    wgt[2:((Nf + 1L) / 2L)] <- 2
  }
  logH <- stats::fft(ce * wgt)
  phase_fine <- Im(logH)

  # ---- sample back to the original n ordinates ----
  sub <- 1L + seq(0L, Nf - 1L, by = fine_factor)
  phase_fine[sub]
}
