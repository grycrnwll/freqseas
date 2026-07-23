# seas_ssi.R
# Stage C orchestration: turn a `seas_test` into a `seas_sa` by removing the
# estimated seasonal filter on the whitened-residual spectrum, recoloring, and
# re-testing the adjusted output (self-check). The phase rule is a declared
# IDENTIFICATION choice, not a tuning knob -- it has no default here. See
# notes/2026-07-22_ssi_phase_surgery.md sec.4/6 and
# notes/2026-07-22_freqseas_architecture.md sec.2.

# Re-clothe a numeric result into the input's class (numeric stays numeric; ts
# regains its tsp; a keyed/unkeyed tsibble is rebuilt from its per-series
# template via fs_rebuild()). The tbl_ts branch closes the "last gap" flagged in
# collection.R: a manual two-stage call on a bare tsibble --
# seas_ssi(seas_test(tsibble, ...), ...) -- dispatches straight here with no
# collection method in between, so without this branch adjusted/seasonal would
# fall through to plain numeric. fs_reclothe_sa() (collection.R) remains a safe
# no-op afterwards: it only rebuilds fields still `is.numeric()`.
.reclothe <- function(vec, object) {
  if (identical(object$input_class, "ts") && !is.null(object$tsp)) {
    return(stats::ts(as.numeric(vec), start = object$tsp[1L],
                     frequency = object$tsp[3L]))
  }
  if (identical(object$input_class, "tbl_ts") && !is.null(object$template)) {
    return(fs_rebuild(as.numeric(vec), object$template))
  }
  as.numeric(vec)
}

# The three-line teaching error for a missing phase rule.
.phase_rule_error <- function() {
  stop(
    "`phase_rule` must be declared explicitly -- it is an identification ",
    "choice, not a default:\n",
    "  \"minimum\" = filter world: seasonality distorts the series' own ",
    "innovations; removal is deconvolution (gain + minimum-phase).\n",
    "  \"zero\"    = components world: seasonality is an independent added ",
    "component; removal is Wiener-style shrinkage (gain only).\n",
    "  See vignette(\"two-worlds\") for the identification story.",
    call. = FALSE
  )
}

#' Seasonal adjustment from a `seas_test` (declared-phase SSI surgery)
#'
#' @description
#' `seas_ssi()` performs Stage C of the workflow: given a completed
#' [seas_test()], it removes the estimated seasonal filter and returns a
#' `seas_sa`. The removal happens on the DFT of the whitened residuals and is
#' inverted through the whitener's recoloring, so seasonality can only re-enter
#' through the cleaned innovations (never via fitted values).
#'
#' `phase_rule` is a required IDENTIFICATION declaration, not a tuning knob:
#' \itemize{
#'   \item `"minimum"` -- filter world (shared innovations): removal is
#'     deconvolution, dividing by the gain AND de-rotating by the implied
#'     minimum phase.
#'   \item `"zero"` -- components world (independent seasonal): removal is
#'     Wiener-style shrinkage, dividing by the gain only (zero phase).
#' }
#' It is always echoed in print/summary tagged `[declared identification]`.
#'
#' If the test decided the series is not seasonal, the input is returned
#' unchanged (with a loud message) wrapped in a `seas_sa` with unit gain and
#' zero phase. Otherwise the band branch estimates the gain from the donor pool
#' and applies the phase rule; the line (deterministic-sinusoid) branch resets
#' the exact-harmonic ordinates to the donor level and leaves phase alone.
#'
#' @param object A `seas_test` object (or, via the collection layer, a
#'   `seas_collection`).
#' @param phase_rule Required. One of `"minimum"` or `"zero"` (see Description).
#'   There is deliberately no default; omitting it raises a teaching error that
#'   points to `vignette("two-worlds")`.
#' @param donor_quantile Upper cutoff on the M0 periodogram distribution when
#'   selecting the donor (reference) ordinates. Default `0.9`.
#' @param ... Passed to methods.
#'
#' @return An object of class `"seas_sa"`: a list with `test` (the embedded
#'   `seas_test`), `adjusted` (input class preserved), `seasonal`
#'   (`x - adjusted`), `Ghat` (gain applied), `theta` (phase applied),
#'   `phase_rule`, `post_evt` (detection p-value of the re-test on the adjusted
#'   series -- a self-check), `spec` (`"band"`, `"line"`, or `"none"`), and
#'   `call`.
#'
#' @examples
#' set.seed(1)
#' n <- 120L; Phi <- 0.8
#' e <- stats::rnorm(n + 400)
#' z <- numeric(n + 400)
#' for (t in 6:(n + 400)) z[t] <- 0.5 * z[t - 1] + Phi * z[t - 4] -
#'   0.5 * Phi * z[t - 5] + e[t]
#' x   <- stats::ts(z[401:(n + 400)], frequency = 4)
#' tst <- seas_test(x)
#' adj <- seas_ssi(tst, phase_rule = "minimum")
#' adj$post_evt          # detection p-value after adjustment (self-check)
#'
#' @export
seas_ssi <- function(object, ...) {
  UseMethod("seas_ssi")
}

#' @rdname seas_ssi
#' @export
seas_ssi.seas_test <- function(object, phase_rule, donor_quantile = 0.9, ...) {
  if (missing(phase_rule) || is.null(phase_rule)) .phase_rule_error()
  phase_rule <- match.arg(phase_rule, c("minimum", "zero"))
  call <- match.call()

  x_num <- as.numeric(object$x)

  # --- not seasonal: return the input unchanged, loudly -------------------
  if (!isTRUE(object$decision)) {
    message("seas_ssi: detection did not reject the white-noise null ",
            sprintf("(p = %.3g >= alpha = %.3g); ", object$evt$p, object$alpha),
            "returning the input series UNCHANGED (no adjustment performed).")
    adjusted <- .reclothe(x_num, object)
    seasonal <- .reclothe(rep(0, length(x_num)), object)
    obj <- new_seas_sa(
      test = object, adjusted = adjusted, seasonal = seasonal,
      Ghat = 1, theta = 0, phase_rule = phase_rule,
      post_evt = object$evt$p, spec = "none", call = call
    )
    return(validate_seas_sa(obj))
  }

  # --- seasonal: set up the donor pool on the raw n_e-grid periodogram ----
  grid <- object$partition$grid
  n_e  <- length(object$e)
  Jmax <- (n_e - 1L) %/% 2L
  pgram_pos <- object$pgram$pgram_raw[grid$pos_idx[seq_len(Jmax)]]
  J1_gain   <- object$partition$J1[object$partition$J1 <= Jmax]
  J0_gain   <- object$partition$J0[object$partition$J0 <= Jmax]
  donor_pool  <- fs_donor_pool(pgram_pos, J0_gain, donor_quantile = donor_quantile)
  donor_level <- mean(pgram_pos[donor_pool])

  spec_label <- object$spec$label

  if (spec_label == "band") {
    gn         <- fs_gain(pgram_pos, sets = list(J1 = J1_gain),
                          donor_pool = donor_pool, n = n_e)
    Ghat_full  <- gn$Ghat_full
    donor_level <- gn$donor_level
    theta_full <- if (phase_rule == "minimum") {
      min_phase(log(Ghat_full))
    } else {
      numeric(n_e)                       # zero phase = components-world rule
    }
    surg     <- fs_surgery(object$E, Ghat_full = Ghat_full,
                           theta_full = theta_full, spec = "band",
                           sets = list(), donor_level = donor_level, n = n_e)
    Ghat_out <- Ghat_full
    theta_out <- theta_full
  } else {
    # line branch: null the exact-harmonic ordinates to the donor level.
    H_full <- grid$pos_idx[object$partition$H]
    if (length(H_full) == 0L) {
      # decision = seasonal but no on-grid exact harmonic to null: nothing to
      # remove under the line specification. Return unchanged, loudly.
      message("seas_ssi: line specification but no on-grid exact-harmonic ",
              "ordinate was identified; returning the input UNCHANGED.")
      adjusted <- .reclothe(x_num, object)
      seasonal <- .reclothe(rep(0, length(x_num)), object)
      obj <- new_seas_sa(
        test = object, adjusted = adjusted, seasonal = seasonal,
        Ghat = 1, theta = 0, phase_rule = phase_rule,
        post_evt = object$evt$p, spec = "line", call = call
      )
      return(validate_seas_sa(obj))
    }
    surg      <- fs_surgery(object$E, Ghat_full = NULL, theta_full = NULL,
                            spec = "line", sets = list(H = H_full),
                            donor_level = donor_level, n = n_e)
    Ghat_out  <- 1
    theta_out <- 0
  }

  # --- recolor to the original scale, class preserved ---------------------
  adjusted_num <- fs_recolor(surg$estar, object$whitener, x_num)
  seasonal_num <- x_num - adjusted_num
  adjusted     <- .reclothe(adjusted_num, object)
  seasonal     <- .reclothe(seasonal_num, object)

  # --- post-adjustment EVT self-check (same M/d/ar_max/exclusion) ---------
  d_forced   <- if (isTRUE(object$whitener$d == 1L)) "first" else "none"
  wh_excl    <- if (is.null(object$whitener$exclusion)) "guard" else object$whitener$exclusion
  wh_guard   <- if (is.null(object$whitener$guard))     3L      else object$whitener$guard
  post <- seas_test(adjusted_num, frequency = object$N, P = object$P,
                    M = object$M, alpha = object$alpha, d = d_forced,
                    ar_max = object$ar_max,
                    whiten_exclusion = wh_excl, whiten_guard = wh_guard)
  post_evt <- post$evt$p

  obj <- new_seas_sa(
    test = object, adjusted = adjusted, seasonal = seasonal,
    Ghat = Ghat_out, theta = theta_out, phase_rule = phase_rule,
    post_evt = post_evt, spec = spec_label, call = call
  )
  validate_seas_sa(obj)
}
