# seas_test.R
# Stage 0 (whiten) + Stage A (detect) + Stage B (specify) orchestration, and
# the user-facing `seas_test()` generic with `.default` (numeric + frequency)
# and `.ts` methods. The `.tbl_ts` (keyed tsibble) method is added by the
# collection layer. See notes/2026-07-22_freqseas_architecture.md sec.2/4/5 and
# notes/2026-07-22_ssi_phase_surgery.md sec.5/6 for the tests and pipeline this
# implements.

# ---------------------------------------------------------------------------
# Internal: block-phase resultant (descriptive specification corroborator)
# ---------------------------------------------------------------------------

# Per-block DFT phase at the fundamental seasonal harmonic, summarized as the
# circular resultant length R = |mean(exp(i * phases))| over K blocks. A frozen
# (deterministic-line) pattern gives R ~ 1; a stochastic band re-randomizes the
# per-block phase, giving R well below 1. Purely descriptive: stored and
# printed, never a decision input in v1. Block length is rounded DOWN to a
# multiple of N so the fundamental lands on a block-DFT ordinate.
.block_phase_R <- function(e, N, P, n_e) {
  K <- max(2L, min(8L, as.integer(floor(n_e / (5L * N)))))
  L <- as.integer(floor(n_e / K))
  L <- L - (L %% N)                       # round DOWN to a multiple of N
  if (L < 2L * N) return(NA_real_)        # too short to see the pattern twice
  omega_fund <- 2 * pi / (P * N)          # smallest seasonal frequency
  k_b <- as.integer(round(omega_fund * L / (2 * pi)))
  if (k_b < 1L || k_b >= L) return(NA_real_)
  phases <- numeric(K)
  for (b in seq_len(K)) {
    idx <- ((b - 1L) * L + 1L):(b * L)
    blk <- e[idx]
    blk <- blk - mean(blk)
    D   <- stats::fft(blk)
    phases[b] <- Arg(D[k_b + 1L])
  }
  Mod(mean(exp(1i * phases)))
}

# ---------------------------------------------------------------------------
# Internal: per-harmonic comb diagnostic
# ---------------------------------------------------------------------------

# Format an angular frequency on (0, pi] as a small rational multiple of pi
# ("pi/2", "pi", "2*pi/3", ...) for the harmonic table / comb note.
.omega_label <- function(w) {
  r <- w / pi
  for (q in 1:24) {
    num <- r * q
    if (abs(num - round(num)) < 1e-6) {
      num <- round(num)
      if (q == 1L)   return(if (num == 1L) "pi" else sprintf("%d*pi", num))
      if (num == 1L) return(sprintf("pi/%d", q))
      return(sprintf("%d*pi/%d", num, q))
    }
  }
  sprintf("%.3f*pi", r)
}

# Per-harmonic excess table, consistent with the EVT machinery: for each exact
# seasonal harmonic (Nyquist included), the standardized-periodogram max in a
# +-1-ordinate neighborhood of the harmonic minus the donor-set (J0) max -- the
# same "seasonal max - nonseasonal max" contrast the detection test uses, but
# localized to one harmonic. A harmonic is "elevated" when that localized excess
# is itself significant at `alpha` under the logistic EVT null with the
# neighborhood's ordinate count. Returns a tibble(omega, label, excess, elevated).
.harmonic_excess_table <- function(I_pos, omega_pos, J0, N0, tau_hat, alpha,
                                   P, N) {
  omega_h <- omega_seasonal(P, N, exclude_nyquist = FALSE)
  n_spec  <- length(omega_pos)
  donor_max <- max(I_pos[J0])
  rows <- lapply(omega_h, function(wh) {
    center <- which.min(abs(omega_pos - wh))
    nbhd   <- unique(pmin(pmax(center + (-1:1), 1L), n_spec))
    excess <- max(I_pos[nbhd]) - donor_max
    p_h    <- evt_pvalue(excess, length(nbhd), N0, tau_hat = tau_hat,
                         alternative = "greater")
    tibble::tibble(omega = wh, label = .omega_label(wh),
                   excess = excess, elevated = p_h < alpha)
  })
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# Internal: the Stage 0 + A + B pipeline on a bare numeric series
# ---------------------------------------------------------------------------

.seas_test_core <- function(x_orig, x_num, input_class, tsp_attr, N, P, M,
                            alpha, d, ar_max, band, call,
                            whiten_exclusion = "guard", whiten_guard = 3L) {
  # --- input validation (informative, before any heavy machinery) --------
  N <- as.integer(round(N))
  if (length(N) != 1L || is.na(N) || N < 2L) {
    stop("`frequency` (N) must be a single integer >= 2.", call. = FALSE)
  }
  P <- as.integer(round(P))
  if (length(P) != 1L || is.na(P) || P < 1L) {
    stop("`P` must be a single positive integer.", call. = FALSE)
  }
  if (any(!is.finite(x_num))) {
    stop("`x` contains non-finite values (NA/NaN/Inf); clean the series first.",
         call. = FALSE)
  }
  n <- length(x_num)
  if (n < 3L * N) {
    stop(sprintf(
      "`x` has length %d but at least 3 * frequency = %d observations are required.",
      n, 3L * N), call. = FALSE)
  }
  if (stats::sd(x_num) == 0) {
    stop("`x` is constant (zero variance); nothing to test.", call. = FALSE)
  }

  # --- M selection (offset-free suggest_M unless user-supplied) ----------
  user_M <- !is.null(M)
  if (user_M) {
    M <- as.integer(round(M))
    if (length(M) != 1L || is.na(M) || M < 2L) {
      stop("`M` must be a single integer >= 2.", call. = FALSE)
    }
    M_for_whiten <- M
  } else {
    # First-pass M from the raw length; suggest_M is offset-free and its M is
    # (in the seasonal regime) insensitive to n_obs, so this equals the
    # n_e-based value below in practice. Re-whitening guards the rare mismatch.
    M_for_whiten <- suggest_M(P, N, n_obs = n, offset_grid = 0)$M
  }

  # --- Stage 0: de-biased-Whittle whitening ------------------------------
  wh  <- fs_whiten(x_num, N = N, P = P, M = M_for_whiten, d = d, ar_max = ar_max,
                   exclusion = whiten_exclusion, guard = whiten_guard)
  n_e <- wh$n_e

  if (user_M) {
    M_selection <- list(source = "user", M = M, descriptor = "user")
  } else {
    sm <- suggest_M(P, N, n_obs = n_e, offset_grid = 0)
    M  <- sm$M
    M_selection <- list(source = "suggest_M", M = M, rho_M = sm$rho_M,
                        H = sm$H, offset_u = sm$offset_u,
                        descriptor = "suggest_M, offset-free")
    if (M != M_for_whiten) {                     # rare: re-whiten at final M
      wh  <- fs_whiten(x_num, N = N, P = P, M = M, d = d, ar_max = ar_max,
                       exclusion = whiten_exclusion, guard = whiten_guard)
      n_e <- wh$n_e
    }
  }
  e <- wh$e

  # --- grid / partition / index sets on the n_e grid ---------------------
  grid <- fourier_grid(n_e)
  part <- mbin_partition(M, P, N, grid$omega_pos, band = band,
                         exclude_nyquist = FALSE)
  isx  <- index_sets(part, grid)

  # --- periodogram of the whitened residuals -----------------------------
  pg      <- fs_periodogram(e, demean = TRUE, standardize = TRUE)
  I_pos   <- pg$pgram_std[grid$pos_idx]           # standardized, positive side
  tau_hat <- pg$tau_hat

  # --- Stage A: EVT detection (max seasonal vs max nonseasonal) ----------
  Delta    <- max(I_pos[isx$J1]) - max(I_pos[isx$J0])
  evt_p    <- evt_pvalue(Delta, isx$N1, isx$N0, tau_hat = tau_hat,
                         alternative = "greater")
  critical <- evt_critical(isx$N1, isx$N0, alpha = alpha, tau_hat = tau_hat,
                           alternative = "greater")
  decision <- evt_p < alpha

  # --- donor pool + level (raw scale; recorded, reused by seas_ssi) ------
  Jmax      <- (n_e - 1L) %/% 2L
  pgram_pos <- pg$pgram_raw[grid$pos_idx[seq_len(Jmax)]]
  J0_gain   <- isx$J0[isx$J0 <= Jmax]
  donor_pool <- tryCatch(
    fs_donor_pool(pgram_pos, J0_gain, donor_quantile = 0.9),
    error = function(e) integer(0)
  )
  donor_level <- if (length(donor_pool) > 0L) mean(pgram_pos[donor_pool]) else NA_real_

  # --- Stage B: shoulder EVT (line vs band) + block-phase resultant ------
  S <- isx$S
  spec_note <- NA_character_
  if (length(S) == 0L) {
    # every seasonal ordinate is an exact harmonic (e.g. tiny grid-aligned n):
    # there are no shoulders to elevate -> a line by construction.
    spec_label <- "line"
    shoulder_p <- NA_real_
    spec_note  <- "no shoulder ordinates (all seasonal ordinates are exact harmonics); spec forced to line"
  } else {
    Delta_sh   <- max(I_pos[S]) - max(I_pos[isx$J0])
    shoulder_p <- evt_pvalue(Delta_sh, length(S), isx$N0, tau_hat = tau_hat,
                             alternative = "greater")
    spec_label <- if (shoulder_p < alpha) "band" else "line"
  }
  phase_R <- .block_phase_R(e, N = N, P = P, n_e = n_e)

  # --- comb diagnostic: per-harmonic excess table + optional comb note ----
  harmonic_table <- .harmonic_excess_table(
    I_pos, grid$omega_pos, isx$J0, isx$N0, tau_hat, alpha, P = P, N = N)
  n_elev <- sum(harmonic_table$elevated)
  # If detection REJECTS but the excess is elevated at a strict, nonempty subset
  # of the harmonics, the "seasonality" may be a nonseasonal cycle sitting at one
  # seasonal frequency rather than a full comb -- flag it.
  comb_note <- NA_character_
  if (isTRUE(decision) && n_elev >= 1L && n_elev < nrow(harmonic_table)) {
    elev_labels <- harmonic_table$label[harmonic_table$elevated]
    comb_note <- sprintf(
      paste0("excess concentrated at %s only: consistent with a nonseasonal ",
             "cycle at that frequency rather than seasonality; see the ",
             "two-worlds vignette"),
      paste(elev_labels, collapse = ", "))
  }

  spec <- list(label = spec_label, shoulder_p = shoulder_p, phase_R = phase_R,
               note = spec_note, harmonic_table = harmonic_table,
               comb_note = comb_note)

  partition <- list(
    grid    = grid,
    breaks  = part$breaks,
    band    = part$band,
    J1      = isx$J1, J0 = isx$J0, H = isx$H, S = isx$S,
    J1_full = isx$J1_full, J0_full = isx$J0_full,
    N1      = isx$N1, N0 = isx$N0,
    M1      = part$M1, M0 = part$M0,
    omega_G = part$omega_G
  )

  evt <- list(statistic = Delta, p = evt_p, critical = critical,
              N1 = isx$N1, N0 = isx$N0)

  obj <- new_seas_test(
    x = x_orig, input_class = input_class, tsp_attr = tsp_attr,
    N = N, P = P, M = M, M_selection = M_selection,
    partition = partition, whitener = wh, e = e, E = pg$dft_raw, pgram = pg,
    donor_pool = donor_pool, donor_level = donor_level, tau_hat = tau_hat,
    evt = evt, spec = spec, decision = decision, alpha = alpha,
    ar_max = as.integer(round(ar_max)), call = call
  )
  validate_seas_test(obj)
}

# ---------------------------------------------------------------------------
# Generic + methods
# ---------------------------------------------------------------------------

#' Frequency-domain seasonality test (detect + specify)
#'
#' @description
#' Runs the freqseas Stage 0 + A + B pipeline on a single series: M0-Whittle
#' whitening (Stage 0), an extreme-value detection test for excess mass at the
#' seasonal frequencies against a white-noise null (Stage A), and a line-versus-
#' band specification test built from a shoulder EVT and a descriptive
#' block-phase resultant (Stage B). Nothing is adjusted here; a not-seasonal
#' series still returns a full object with the decision recorded.
#'
#' `seas_test()` is an S3 generic. The `.default` method takes a bare numeric
#' vector and a stated `frequency`; the `.ts` method infers `frequency` from
#' the series. A `.tbl_ts` (keyed tsibble) method is provided by the collection
#' layer.
#'
#' @param x A numeric vector (with `frequency`) or a `stats::ts` object.
#' @param frequency Observations per year (e.g. 4 quarterly, 12 monthly).
#'   Required for a bare numeric `x`; inferred for `ts`.
#' @param P Positive integer fundamental-per-year multiplier. Default `1`.
#' @param M Number of equal-width frequency bins on `(0, pi]`. Default `NULL`,
#'   in which case the offset-free [suggest_M()] chooses it (recorded in the
#'   result).
#' @param alpha Test level for the detection and shoulder tests. Default `0.05`.
#' @param d Differencing rule passed to the whitener: `"auto"` (difference once
#'   iff near-integrated), `"none"`, or `"first"`. Default `"auto"`.
#' @param ar_max Maximum AR order for the de-biased-Whittle BIC search. Default
#'   `3`.
#' @param band Bin-restriction rule: `"full"` (default) or
#'   `"between_seasonal_extremes"`.
#' @param whiten_exclusion Ordinate-exclusion geometry for the whitening fit,
#'   passed to [fs_whiten()]: `"guard"` (default) or `"bins"`. See
#'   [fs_whiten()]'s Details for the size/power tradeoff and the measured
#'   numbers behind the default.
#' @param whiten_guard Guard half-width (Fourier steps) for
#'   `whiten_exclusion = "guard"`. Default `2`.
#' @param ... Passed to methods.
#'
#' @return An object of class `"seas_test"`: a list with fields `x`,
#'   `input_class`, `N`, `P`, `M`, `M_selection`, `partition` (grid, breaks,
#'   `J1`/`J0`/`H`/`S` and their full-DFT counterparts, `M1`/`M0`, `omega_G`),
#'   `whitener` (the [fs_whiten()] record: `d`, `d_reason`, `mu`, `p`, `ar`,
#'   `bic_path`, ...), `e` (whitened residuals), `E` (their DFT), `pgram`,
#'   `donor_pool`, `donor_level`, `tau_hat`, `evt` (`statistic`, `p`,
#'   `critical`), `spec` (`label` = `"line"`/`"band"`, `shoulder_p`,
#'   `phase_R`, `harmonic_table` = a per-harmonic excess tibble with columns
#'   `omega`/`label`/`excess`/`elevated`, and `comb_note` = a character flag,
#'   `NA` unless the excess is elevated at only a strict subset of harmonics),
#'   `decision` (logical), `alpha`, `ar_max`, and `call`.
#'
#' @examples
#' # Quarterly series with a seasonal band -> detected, specified "band".
#' set.seed(1)
#' n <- 120L; Phi <- 0.8
#' e <- stats::rnorm(n + 400)
#' z <- numeric(n + 400)
#' for (t in 6:(n + 400)) z[t] <- 0.5 * z[t - 1] + Phi * z[t - 4] -
#'   0.5 * Phi * z[t - 5] + e[t]
#' x <- stats::ts(z[401:(n + 400)], frequency = 4)
#' tst <- seas_test(x)
#' tst$decision
#' tst$spec$label
#'
#' @export
seas_test <- function(x, ...) {
  UseMethod("seas_test")
}

#' @rdname seas_test
#' @export
seas_test.default <- function(x, frequency = NULL, P = 1, M = NULL,
                              alpha = 0.05, d = c("auto", "none", "first"),
                              ar_max = 3, band = c("full", "between_seasonal_extremes"),
                              whiten_exclusion = c("guard", "bins"),
                              whiten_guard = 3L,
                              ...) {
  d                <- match.arg(d)
  band             <- match.arg(band)
  whiten_exclusion <- match.arg(whiten_exclusion)
  if (!is.numeric(x)) {
    stop("`seas_test.default()` expects a numeric `x`; supply a `ts` (or a ",
         "keyed tsibble) for the class-aware methods.", call. = FALSE)
  }
  if (is.null(frequency)) {
    stop("`frequency` is required when `x` is a bare numeric vector ",
         "(observations per year, e.g. 4 for quarterly).", call. = FALSE)
  }
  .seas_test_core(
    x_orig = as.numeric(x), x_num = as.numeric(x), input_class = "numeric",
    tsp_attr = NULL, N = frequency, P = P, M = M, alpha = alpha, d = d,
    ar_max = ar_max, band = band, call = match.call(),
    whiten_exclusion = whiten_exclusion, whiten_guard = whiten_guard
  )
}

#' @rdname seas_test
#' @export
seas_test.ts <- function(x, P = 1, M = NULL, alpha = 0.05,
                         d = c("auto", "none", "first"), ar_max = 3,
                         band = c("full", "between_seasonal_extremes"),
                         whiten_exclusion = c("guard", "bins"),
                         whiten_guard = 3L, ...) {
  d                <- match.arg(d)
  band             <- match.arg(band)
  whiten_exclusion <- match.arg(whiten_exclusion)
  N                <- stats::frequency(x)
  .seas_test_core(
    x_orig = x, x_num = as.numeric(x), input_class = "ts",
    tsp_attr = stats::tsp(x), N = N, P = P, M = M, alpha = alpha, d = d,
    ar_max = ar_max, band = band, call = match.call(),
    whiten_exclusion = whiten_exclusion, whiten_guard = whiten_guard
  )
}
