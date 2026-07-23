# surgery.R
# DFT-side removal of the estimated seasonal filter from a whitened
# residual spectrum (internal). See notes/2026-07-22_ssi_phase_surgery.md
# SS4-6 for the band/line, minimum-phase/zero-phase theory this
# implements, and helper_functions/replace_ordinates.R in
# seasonality_vintages for the ordinate-replacement conventions (donor
# magnitude at retained phase, conjugate-symmetry enforcement) this ports.
#
# Coordinate convention: `sets$H` is given in FULL-DFT-vector coordinates
# (1-indexed, DC at position 1) -- the same coordinate system as `E`,
# `Ghat_full`, and `theta_full` -- restricted to the positive-frequency
# side only. This differs from `sets$J1` in gain.R, which is given in
# positive-frequency-ARRAY coordinates (indices into a length-Jmax vector).
# The two files use different coordinate systems because gain.R only ever
# indexes into positive-frequency-only arrays (pgram_pos), while
# surgery.R operates directly on the full-length complex DFT.

#' DFT surgery: remove the estimated filter from a whitened residual spectrum (internal)
#'
#' @description
#' Applies seasonal-filter removal to the full DFT of a whitened residual
#' series, in one of two specifications (see \code{spec}):
#' \itemize{
#'   \item \strong{"band"}: the seasonal feature is a stochastic band
#'     around each harmonic. Every ordinate is divided by the estimated
#'     gain and de-rotated by \code{theta_full}. Passing
#'     \code{theta_full = min_phase(log(Ghat_full))} implements the
#'     filter-world (minimum-phase) removal rule; passing a vector of
#'     zeros implements the components-world (zero-phase, Wiener-style)
#'     rule. \code{fs_surgery()} itself is agnostic to which rule produced
#'     \code{theta_full} -- the phase-rule declaration happens upstream,
#'     where \code{theta_full} is built.
#'   \item \strong{"line"}: the seasonal feature is a deterministic
#'     sinusoid frozen at the exact harmonic ordinates. Only those
#'     ordinates (\code{sets$H}) and their conjugate mirrors are touched;
#'     their magnitude is reset to the donor (nonseasonal reference) level
#'     and their phase is left alone -- a deterministic line has no filter
#'     phase to de-rotate. \code{Ghat_full} and \code{theta_full} are
#'     ignored in this branch.
#' }
#'
#' Guards are enforced identically in both branches, exactly as specified
#' and tested: the DC ordinate is left byte-for-byte identical to the
#' input, the Nyquist ordinate (if \code{n} is even) is forced to be
#' exactly real, and Hermitian symmetry is rebuilt -- for every
#' positive-side ordinate the branch touched, its negative-side mirror is
#' overwritten with its exact conjugate, rather than trusting elementwise
#' arithmetic on \code{Ghat_full}/\code{theta_full} (which need not be
#' symmetric to full floating-point precision) or the untouched input's
#' own symmetry (\code{stats::fft()} of a real vector is Hermitian only up
#' to floating-point rounding, not bit-for-bit). Ordinates the branch did
#' not touch are left as literal copies of \code{E} -- this is what makes
#' the \code{"line"} branch's untouched ordinates bit-identical to the
#' input.
#'
#' @param E Complex vector of length \code{n}: the full DFT (as from
#'   \code{stats::fft()}) of the whitened residual series.
#' @param Ghat_full Numeric vector of length \code{n}: linear-scale gain on
#'   the full DFT circle (e.g. \code{fs_gain()$Ghat_full}). Used only by
#'   the \code{"band"} branch.
#' @param theta_full Numeric vector of length \code{n}: phase (radians) on
#'   the full DFT circle to remove (e.g.
#'   \code{min_phase(log(Ghat_full))} for the minimum-phase rule, or a
#'   zero vector for the zero-phase rule). Used only by the \code{"band"}
#'   branch.
#' @param spec Character, one of \code{"band"} or \code{"line"}: which
#'   surgery specification to apply.
#' @param sets List with (at least) element \code{H} for the
#'   \code{"line"} branch: integer vector of positions, in full-DFT
#'   coordinates (1-indexed, DC at 1), of the exact-harmonic ordinates on
#'   the positive side of the circle only (negative-side mirrors are
#'   derived internally -- do not include them in \code{sets$H}). Ignored
#'   by the \code{"band"} branch.
#' @param donor_level Single positive numeric, the nonseasonal reference
#'   periodogram level (e.g. \code{fs_gain()$donor_level}). Used only by
#'   the \code{"line"} branch: the target DFT magnitude is
#'   \code{sqrt(donor_level * n)}, since the periodogram convention here is
#'   \eqn{I_k = |X_k|^2 / n}, so a magnitude of
#'   \eqn{\sqrt{n \cdot \code{donor\_level}}} is exactly the one whose
#'   periodogram value equals \code{donor_level}.
#' @param n Single positive integer, the DFT length (\code{length(E)}).
#'
#' @return Named list:
#'   \item{Estar}{Complex vector, length \code{n}: the surgered DFT --
#'     exactly Hermitian-symmetric where rebuilt, DC identical to
#'     \code{E[1]}, Nyquist (if present) exactly real.}
#'   \item{estar}{Numeric vector, length \code{n}: the real-valued
#'     reconstructed series, \code{Re(stats::fft(Estar, inverse = TRUE)) /
#'     n}.}
fs_surgery <- function(E, Ghat_full, theta_full, spec = c("band", "line"),
                        sets, donor_level, n) {
  spec <- match.arg(spec)
  if (length(n) != 1L || !is.finite(n) || n < 2L) {
    stop("`n` must be a single integer >= 2.")
  }
  n <- as.integer(round(n))
  if (!is.complex(E) || length(E) != n) {
    stop("`E` must be a complex vector of length `n` (the DFT of a real series).")
  }

  half    <- n %/% 2L
  is_even <- (n %% 2L == 0L)
  pos_range <- 2:(if (is_even) half else half + 1L)  # positive-frequency positions, full-DFT coords

  if (spec == "band") {
    if (length(Ghat_full) != n) stop("`Ghat_full` must have length `n`.")
    if (length(theta_full) != n) stop("`theta_full` must have length `n`.")
    if (any(!is.finite(Ghat_full)) || any(Ghat_full <= 0)) {
      stop("`Ghat_full` must be strictly positive and finite.")
    }
    if (any(!is.finite(theta_full))) stop("`theta_full` must be finite.")

    Estar   <- E * (1 / Ghat_full) * exp(-1i * theta_full)
    touched <- pos_range
  } else {
    if (is.null(sets) || is.null(sets$H)) {
      stop("`sets$H` (exact-harmonic ordinate positions, full-DFT coordinates) ",
           "is required for spec = \"line\".")
    }
    H <- as.integer(sets$H)
    if (length(H) == 0L) stop("`sets$H` is empty; nothing to surgery for spec = \"line\".")
    if (any(H < 2L | H > n)) stop("`sets$H` contains positions outside [2, n].")
    if (length(donor_level) != 1L || !is.finite(donor_level) || donor_level <= 0) {
      stop("`donor_level` must be a single positive finite value.")
    }

    Estar <- E
    target_mag <- sqrt(donor_level * n)
    ph <- Arg(Estar[H])
    ph[!is.finite(ph)] <- 0
    Estar[H] <- complex(modulus = target_mag, argument = ph)
    touched <- H
  }

  # ---- guards (both branches) ----
  # Hermitian symmetry: for every touched positive-side position, rebuild
  # its mirror as the exact conjugate. DC and Nyquist are excluded here
  # and handled by the dedicated guards below.
  touched <- setdiff(unique(touched), 1L)
  if (is_even) touched <- setdiff(touched, half + 1L)
  if (length(touched) > 0L) {
    mirror <- n + 2L - touched
    Estar[mirror] <- Conj(Estar[touched])
  }

  Estar[1] <- E[1]                                     # DC: untouched, exactly

  if (is_even) {
    nyq <- half + 1L
    Estar[nyq] <- complex(real = Re(Estar[nyq]), imaginary = 0)  # Nyquist: forced real
  }

  estar <- Re(stats::fft(Estar, inverse = TRUE)) / n
  list(Estar = Estar, estar = estar)
}
