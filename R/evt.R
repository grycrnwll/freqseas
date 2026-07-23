# Periodogram and extreme-value null machinery for the freqseas detection test.
#
# The detection statistic is a max-periodogram contrast between seasonal and
# nonseasonal Fourier ordinates; under a white null its distribution is
# approximately logistic. This file supplies the pieces (standardized
# periodogram, the logistic critical value, the logistic p-value); the wave-2
# `seas_test()` orchestrator assembles the statistic from them.

#' Series normalization for the standardized periodogram
#'
#' Subtract the sample mean (optional), then divide by the sample standard
#' deviation (optional). Returns the normalized series plus the scalars the
#' downstream test and null distribution need. After demean + divide-by-sd the
#' periodogram is on a unit-variance scale and the logistic null uses
#' `tau_hat = 1`; with `standardize = FALSE`, `tau_hat = var(x_dm)` so the null
#' can be rescaled.
#'
#' @param x Numeric vector or `ts`.
#' @param demean Logical; subtract `mean(x)` before any scaling. Default `TRUE`.
#' @param standardize Logical; divide by the sd of the (possibly demeaned)
#'   series. Default `TRUE`.
#'
#' @return A named list: `x` (normalized series), `sd` (scale used, `NA` when
#'   `standardize = FALSE`), `tau_hat` (`1` if standardized, else `var(x_dm)`),
#'   and the echoed `demean` / `standardize` flags.
#'
#' @keywords internal
#' @noRd
fs_normalize <- function(x, demean = TRUE, standardize = TRUE) {
  x_num <- as.numeric(x)
  if (length(x_num) < 2L) stop("`x` must have length >= 2.")
  if (any(!is.finite(x_num))) stop("`x` contains non-finite values.")

  if (isTRUE(demean)) {
    x_dm <- x_num - mean(x_num)
  } else {
    x_dm <- x_num
  }

  if (isTRUE(standardize)) {
    sx <- stats::sd(x_dm)
    if (!is.finite(sx) || sx <= 0) {
      stop("Series has zero (or undefined) standard deviation after demeaning.")
    }
    x_out   <- x_dm / sx
    tau_hat <- 1
  } else {
    sx      <- NA_real_
    x_out   <- x_dm
    tau_hat <- stats::var(x_dm)
  }

  list(
    x           = x_out,
    sd          = sx,
    tau_hat     = tau_hat,
    demean      = isTRUE(demean),
    standardize = isTRUE(standardize)
  )
}

#' DFT and periodogram (raw and standardized)
#'
#' Compute the DFT and periodogram of a real series, returning both raw and
#' standardized variants in one pass. The detection test uses the standardized
#' periodogram (so the statistic lives on a comparable scale across series); the
#' adjustment uses the raw DFT (so phase and donor magnitudes round-trip). The
#' pre-FFT scaling is delegated to [fs_normalize()] as a single source of truth.
#'
#' @param x Numeric vector or `ts` (length `>= 8`).
#' @param demean Logical; subtract the mean before the FFT. Default `TRUE`.
#' @param standardize Logical; divide by sd after demeaning. Default `TRUE`.
#'
#' @return A named list:
#'   \describe{
#'     \item{`n`}{Length of `x`.}
#'     \item{`dft_raw`}{Complex FFT of raw `x`.}
#'     \item{`pgram_raw`}{`Mod(dft_raw)^2 / n` (length `n`).}
#'     \item{`dft_std`}{Complex FFT of the normalized series.}
#'     \item{`pgram_std`}{`Mod(dft_std)^2 / n` (length `n`).}
#'     \item{`sd`}{Scale from [fs_normalize()] (`NA` if not standardized).}
#'     \item{`tau_hat`}{`1` if standardized, else `var(x_dm)`.}
#'     \item{`demean`, `standardize`}{Echoed flags.}
#'   }
#'
#' @keywords internal
#' @noRd
fs_periodogram <- function(x, demean = TRUE, standardize = TRUE) {
  x_num <- as.numeric(x)
  n     <- length(x_num)
  if (n < 8L) stop("Series too short (need at least 8 observations).")

  # raw
  dft_raw   <- stats::fft(x_num)
  pgram_raw <- (Mod(dft_raw)^2) / n

  # standardized -- delegate to fs_normalize()
  norm      <- fs_normalize(x_num, demean = demean, standardize = standardize)
  dft_std   <- stats::fft(norm$x)
  pgram_std <- (Mod(dft_std)^2) / n

  list(
    n           = n,
    dft_raw     = dft_raw,
    pgram_raw   = pgram_raw,
    dft_std     = dft_std,
    pgram_std   = pgram_std,
    sd          = norm$sd,
    tau_hat     = norm$tau_hat,
    demean      = norm$demean,
    standardize = norm$standardize
  )
}

#' Logistic critical value for the EVT seasonality statistic
#'
#' Under the white-noise null the max-contrast statistic
#' `Delta = max_{J1} I - max_{J0} I` is approximately
#' `Logistic(location = tau_hat * log(N1 / N0), scale = tau_hat)`. This returns
#' the critical value at level `alpha`.
#'
#' @param N1,N0 Seasonal / nonseasonal Fourier-ordinate counts (positive).
#' @param alpha Test level in `(0, 1)`. Default `0.05`.
#' @param tau_hat Null scale; `1` for a standardized periodogram, else
#'   `var(x_dm)` (see [fs_periodogram()]). Default `1`.
#' @param alternative `"greater"` (default), `"less"`, or `"two.sided"`.
#'
#' @return Single numeric critical value `c_alpha`: reject when the observed
#'   `Delta` exceeds it (`"greater"`), is below it (`"less"`), or as the upper
#'   bound for `"two.sided"`.
#'
#' @keywords internal
#' @noRd
evt_critical <- function(N1, N0,
                         alpha       = 0.05,
                         tau_hat     = 1,
                         alternative = c("greater", "two.sided", "less")) {
  alternative <- match.arg(alternative)
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single value in (0, 1).")
  }
  if (N1 < 1L || N0 < 1L) stop("`N1` and `N0` must be positive.")
  if (!is.finite(tau_hat) || tau_hat <= 0) stop("`tau_hat` must be positive.")

  location <- tau_hat * log(N1 / N0)
  scale    <- tau_hat

  switch(
    alternative,
    "greater"   = stats::qlogis(1 - alpha,     location = location, scale = scale),
    "less"      = stats::qlogis(alpha,         location = location, scale = scale),
    "two.sided" = stats::qlogis(1 - alpha / 2, location = location, scale = scale)
  )
}

#' Logistic p-value for the EVT seasonality statistic
#'
#' The p-value companion to [evt_critical()], using the same logistic null.
#'
#' @param Delta Observed test statistic.
#' @param N1,N0 Seasonal / nonseasonal Fourier-ordinate counts.
#' @param tau_hat Null scale (see [evt_critical()]). Default `1`.
#' @param alternative `"greater"` (default), `"less"`, or `"two.sided"`.
#'
#' @return Numeric p-value clipped to `[0, 1]`.
#'
#' @keywords internal
#' @noRd
evt_pvalue <- function(Delta, N1, N0,
                       tau_hat     = 1,
                       alternative = c("greater", "two.sided", "less")) {
  alternative <- match.arg(alternative)
  if (!is.finite(tau_hat) || tau_hat <= 0) stop("`tau_hat` must be positive.")
  location <- tau_hat * log(N1 / N0)
  scale    <- tau_hat

  cdf <- stats::plogis(Delta, location = location, scale = scale)
  p <- switch(
    alternative,
    "greater"   = 1 - cdf,
    "less"      = cdf,
    "two.sided" = 2 * min(cdf, 1 - cdf)
  )
  max(min(p, 1), 0)
}
