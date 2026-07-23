# gain.R
# Donor pool selection and gain estimation for the SSI surgery step
# (internal). Ports the conventions of
# seasonality_vintages/helper_functions/donor_pool.R and
# replace_ordinates.R, adapted to positive-frequency-only arrays and to
# gain estimation (rather than ordinate replacement) as the surgery target.
# See notes/2026-07-22_ssi_phase_surgery.md SS3 for the theory.
#
# Coordinate convention used throughout this file (and by fs_surgery()'s
# `sets$J1`/`sets$H` where relevant): "positive-frequency-array coordinates"
# means an index into a length-Jmax vector holding one value per positive
# ordinate j = 1..Jmax (Jmax = (n - 1) %/% 2), the same ordering as
# `fourier_grid(n)$omega_pos` / `index_sets()$J1`/`$J0`. This is NOT the
# same coordinate system as positions in the full length-n DFT vector
# (which is what `sets$H` uses in surgery.R, and what `index_sets()` calls
# `J1_full`/`J0_full`).

#' Donor pool of nonseasonal Fourier ordinates for gain estimation (internal)
#'
#' @description
#' Selects the low-power subset of the nonseasonal (M0) positive-frequency
#' ordinates to use as the reference ("donor") level against which
#' seasonal-bin (J1) ordinates are compared when estimating gain. An
#' ordinate qualifies as a donor if its periodogram value is at or below
#' the \code{donor_quantile} quantile of the M0 periodogram distribution --
#' this excludes any M0 ordinate that happens to carry unusually high power
#' (leakage from a nearby seasonal peak, or any other spectral feature)
#' from inflating the donor (reference) level.
#'
#' @param pgram_pos Numeric vector, the periodogram on the positive-frequency
#'   ordinates only (length \code{Jmax}; positive-frequency-array
#'   coordinates -- see file header).
#' @param J0 Integer vector, positions (into \code{pgram_pos}, i.e.
#'   positive-frequency-array coordinates) of the nonseasonal (M0)
#'   ordinates.
#' @param donor_quantile Single numeric in (0, 1], default 0.9. Upper
#'   cutoff on the M0 periodogram distribution.
#'
#' @return Integer vector, a subset of \code{J0}: positions (in
#'   \code{pgram_pos} coordinates) of the accepted donor ordinates. Errors
#'   informatively if fewer than 5 donors are available, either because
#'   \code{J0} itself is too small or because the quantile filter leaves
#'   too few.
fs_donor_pool <- function(pgram_pos, J0, donor_quantile = 0.9) {
  min_donors <- 5L
  if (length(donor_quantile) != 1L || !is.finite(donor_quantile) ||
      donor_quantile <= 0 || donor_quantile > 1) {
    stop("`donor_quantile` must be a single value in (0, 1].")
  }
  if (length(J0) < min_donors) {
    stop(sprintf(
      "Too few nonseasonal ordinates (%d < %d) to build a donor pool.",
      length(J0), min_donors
    ))
  }
  if (any(J0 < 1L | J0 > length(pgram_pos))) {
    stop("`J0` contains positions outside `pgram_pos`.")
  }

  pg0 <- pgram_pos[J0]
  if (any(!is.finite(pg0))) {
    stop("Non-finite periodogram values among nonseasonal (J0) ordinates.")
  }

  thr  <- stats::quantile(pg0, probs = donor_quantile, names = FALSE)
  pool <- J0[pg0 <= thr]

  if (length(pool) < min_donors) {
    stop(sprintf(
      "Too few donors after quantile filter (%d < %d). Raise `donor_quantile`.",
      length(pool), min_donors
    ))
  }
  pool
}

#' Gain of the filter being removed, estimated from donor-normalized power (internal)
#'
#' @description
#' Estimates the (linear-scale) gain \eqn{|H(\omega)| \ge 1} of the
#' seasonal filter that \code{fs_surgery()} removes, from the ratio of
#' seasonal-bin (J1) periodogram power to the donor-pool reference level
#' (\code{donor_level}). \code{fs_gain} always estimates gain for the
#' filter being REMOVED (hence \eqn{\ge 1}); \code{fs_surgery} is the one
#' that divides by it. Gain is floored at 1 -- the surgery never amplifies
#' an ordinate relative to the donor level -- then lightly smoothed in log
#' space (so the correction is not driven by single-ordinate periodogram
#' noise) and expanded to a full-circle, Hermitian-symmetric vector ready
#' for \code{fs_surgery()}.
#'
#' The smoothing is a length-3 symmetric moving average applied to
#' \eqn{\log(G) = \log(G^2)/2}, restricted to the seasonal bins (J1)
#' expanded by one ordinate on each side (so the correction tapers smoothly
#' into the donor level rather than presenting a hard log-gain cliff at the
#' bin edge); ordinates more than one step outside J1 are left at log-gain
#' 0 (gain 1, i.e. untouched).
#'
#' @param pgram_pos Numeric vector, periodogram on positive-frequency
#'   ordinates only (length \code{Jmax}; see \code{fs_donor_pool()}).
#' @param sets List with (at least) element \code{J1}: integer vector,
#'   positions (into \code{pgram_pos}, positive-frequency-array
#'   coordinates) of the seasonal-bin ordinates.
#' @param donor_pool Integer vector, positions (into \code{pgram_pos}) of
#'   the donor ordinates, as returned by \code{fs_donor_pool()}.
#' @param n Single positive integer, the full series/DFT length (\code{n}
#'   in \code{stats::fft()} terms; \code{pgram_pos} must have
#'   \code{Jmax = (n - 1) \%/\% 2} entries).
#'
#' @return Named list:
#'   \item{Ghat_full}{Numeric vector, length \code{n}: linear-scale gain on
#'     the full DFT circle, Hermitian-symmetric
#'     (\code{Ghat_full[j] == Ghat_full[n - j]}), DC = 1, Nyquist (if
#'     \code{n} is even) = 1 (the Nyquist ordinate is outside \code{pgram_pos}
#'     so it gets the neutral value here; when it is a seasonal harmonic the
#'     caller must fill it in -- see Details), ready to divide directly into a
#'     DFT.}
#'   \item{Ghat_pos}{Numeric vector, length \code{Jmax}: linear-scale gain
#'     on the positive-frequency axis only (the source for the positive
#'     half of \code{Ghat_full}).}
#'   \item{donor_level}{Single numeric: mean periodogram over the donor
#'     pool -- the reference power level that defines "gain 1".}
#'
#' @details
#' \code{pgram_pos} has length \code{Jmax = (n - 1) \%/\% 2}, which for even
#' \code{n} excludes the Nyquist ordinate, so no periodogram-based gain
#' estimate is available for it here; \code{Ghat_full} assigns it the neutral
#' (donor-level, gain 1) value, consistent with the "1 for all other j" rule
#' applied to every ordinate outside J1.
#'
#' Note that \code{omega = pi} \emph{is} a seasonal harmonic in the freqseas
#' pipeline whenever the seasonal period \code{N} is even: the partition
#' machinery is called with \code{exclude_nyquist = FALSE} throughout, and
#' \code{fs_whiten()} aligns the residual length to a multiple of \code{N}, so
#' for even \code{N} the Nyquist ordinate always exists and always carries a
#' comb tooth. The neutral value returned here is therefore \emph{not} the
#' right gain for that ordinate. \code{seas_ssi()}'s band branch fills it in
#' from the raw periodogram after calling this function; a caller using
#' \code{fs_gain()} directly must do the same.
fs_gain <- function(pgram_pos, sets, donor_pool, n) {
  if (length(n) != 1L || !is.finite(n) || n < 2L) {
    stop("`n` must be a single integer >= 2.")
  }
  n    <- as.integer(round(n))
  Jmax <- (n - 1L) %/% 2L
  if (length(pgram_pos) != Jmax) {
    stop(sprintf(
      "`pgram_pos` has length %d but n = %d implies Jmax = %d positive ordinates.",
      length(pgram_pos), n, Jmax
    ))
  }

  J1 <- sort(unique(sets$J1))
  if (length(J1) == 0L) {
    stop("`sets$J1` is empty; no seasonal ordinates to estimate gain for.")
  }
  if (any(J1 < 1L | J1 > Jmax)) {
    stop("`sets$J1` contains positions outside `pgram_pos`.")
  }
  if (length(donor_pool) == 0L) stop("`donor_pool` is empty.")
  if (any(donor_pool < 1L | donor_pool > Jmax)) {
    stop("`donor_pool` contains positions outside `pgram_pos`.")
  }

  donor_level <- mean(pgram_pos[donor_pool])
  if (!is.finite(donor_level) || donor_level <= 0) {
    stop("`donor_level` (mean periodogram over the donor pool) is non-finite or non-positive.")
  }

  # ---- raw squared gain: elevated in J1, 1 (donor level) elsewhere ----
  G2 <- rep(1, Jmax)
  G2[J1] <- pgram_pos[J1] / donor_level
  G2 <- pmax(G2, 1)              # floor: never amplify

  lg <- log(G2) / 2               # log gain (G2 is squared gain)

  # ---- smooth within J1 expanded by 1 ordinate each side; identity elsewhere ----
  J1_exp <- sort(unique(pmin(pmax(c(J1 - 1L, J1, J1 + 1L), 1L), Jmax)))
  lg_left  <- c(NA_real_, lg[-Jmax])
  lg_right <- c(lg[-1L], NA_real_)
  sm <- rowMeans(cbind(lg_left, lg, lg_right), na.rm = TRUE)

  lg_smoothed <- numeric(Jmax)
  lg_smoothed[J1_exp] <- sm[J1_exp]

  Ghat_pos <- exp(lg_smoothed)

  # ---- build the full-circle Hermitian-symmetric gain vector ----
  Ghat_full <- numeric(n)
  Ghat_full[1] <- 1                          # DC
  pos_idx <- (1:Jmax) + 1L                    # positions of j = 1..Jmax in full-DFT coords
  neg_idx <- n + 1L - (1:Jmax)                 # their Hermitian mirrors
  Ghat_full[pos_idx] <- Ghat_pos
  Ghat_full[neg_idx] <- Ghat_pos
  if (n %% 2L == 0L) {
    Ghat_full[n %/% 2L + 1L] <- 1              # Nyquist: not covered by J1 -> donor level
  }

  list(Ghat_full = Ghat_full, Ghat_pos = Ghat_pos, donor_level = donor_level)
}
