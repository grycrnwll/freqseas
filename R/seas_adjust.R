# seas_adjust.R
# The one-call wrapper: seas_adjust(x) == seas_ssi(seas_test(x), phase_rule).
# A single code path -- the wrapper IS the two calls. The ONLY asymmetry with
# the two-stage path is the default: at THIS wrapper layer `phase_rule`
# defaults to "minimum" (and is printed loudly as a declared identification),
# whereas seas_ssi() itself REQUIRES an explicit `phase_rule`. The default is a
# convenience of the wrapper, never a claim that "minimum" is generically
# correct -- the filter-versus-components identification is not testable at
# second order (see notes/2026-07-22_freqseas_architecture.md sec.2/6 and
# notes/2026-07-22_ssi_phase_surgery.md sec.4).

#' One-call seasonal adjustment
#'
#' @description
#' Convenience wrapper equal to `seas_ssi(seas_test(x, ...), phase_rule =
#' phase_rule, donor_quantile = donor_quantile)`. It runs the full detect ->
#' specify -> identify -> operate workflow in a single call.
#'
#' \strong{Default asymmetry (read this).} At this wrapper layer `phase_rule`
#' defaults to `"minimum"` (the filter-world / deconvolution rule) so that the
#' one-call path is usable out of the box, and the choice is echoed in every
#' print/summary tagged `[declared identification]`. The lower-level
#' [seas_ssi()] deliberately has \emph{no} default and forces you to declare
#' the rule, because filter-world versus components-world is an identification
#' commitment that is not testable at second order -- choosing it is choosing
#' which counterfactual you claim to recover. Pass `phase_rule = "zero"` for the
#' components-world (Wiener-shrinkage, classic-SSI) rule.
#'
#' @param x A numeric vector (with `frequency`) or a `stats::ts` object.
#' @param phase_rule Identification rule, `"minimum"` (default) or `"zero"`.
#'   See [seas_ssi()].
#' @param donor_quantile Upper cutoff on the M0 periodogram distribution for the
#'   donor pool. Default `0.9`.
#' @param ... Passed to [seas_test()] (e.g. `frequency`, `P`, `M`, `alpha`,
#'   `d`, `ar_max`, `band`, `whiten_exclusion`, `whiten_guard`).
#'
#' @return A `seas_sa` object; see [seas_ssi()].
#'
#' @examples
#' set.seed(1)
#' n <- 120L; Phi <- 0.8
#' e <- stats::rnorm(n + 400)
#' z <- numeric(n + 400)
#' for (t in 6:(n + 400)) z[t] <- 0.5 * z[t - 1] + Phi * z[t - 4] -
#'   0.5 * Phi * z[t - 5] + e[t]
#' x <- stats::ts(z[401:(n + 400)], frequency = 4)
#'
#' adj_min  <- seas_adjust(x)                      # default: minimum phase
#' adj_zero <- seas_adjust(x, phase_rule = "zero") # components-world rule
#'
#' @export
seas_adjust <- function(x, phase_rule = "minimum", donor_quantile = 0.9, ...) {
  UseMethod("seas_adjust")
}

#' @rdname seas_adjust
#' @export
seas_adjust.default <- function(x, phase_rule = "minimum", donor_quantile = 0.9,
                                ...) {
  phase_rule <- match.arg(phase_rule, c("minimum", "zero"))
  tst <- seas_test.default(x, ...)
  seas_ssi(tst, phase_rule = phase_rule, donor_quantile = donor_quantile)
}

#' @rdname seas_adjust
#' @export
seas_adjust.ts <- function(x, phase_rule = "minimum", donor_quantile = 0.9,
                           ...) {
  phase_rule <- match.arg(phase_rule, c("minimum", "zero"))
  tst <- seas_test.ts(x, ...)
  seas_ssi(tst, phase_rule = phase_rule, donor_quantile = donor_quantile)
}
