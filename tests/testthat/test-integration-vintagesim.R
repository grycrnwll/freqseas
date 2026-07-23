# TIER 4 integration tests: the freqseas workflow against vintagesim, the
# controlled vintage-simulation environment whose latent nonseasonal
# counterfactual `xNS` is known by construction (see the vintagesim README and
# notes/2026-07-22_freqseas_architecture.md sec.8 tier 4). vintagesim is a
# Suggested dependency, so every test is gated on
# skip_if_not_installed("vintagesim") + skip_on_cran().
#
# The DGP vintagesim implements is the FILTER world (shared innovations:
# (1 - phi B)(1 - Phi B^s)(x - alpha) = eps, (1 - phi B)(xNS - alpha) = eps), so
# the minimum-phase rule is the correct identification here.

# --- assembly-time integration guard (see test-seas-ssi.R for the rationale).
# This file sorts AFTER test-collection.R, whose `<<-` stubs of the real
# seas_test()/seas_ssi()/seas_adjust() layer would, under the source-into-global
# verification harness, still be installed when we run. The keyed-tsibble path
# (c) additionally needs the REAL collection layer (seas_test.tbl_ts,
# seas_ssi.seas_collection, seas_adjust.tbl_ts, fs_extract/fs_rebuild), which
# the 5-file subset used elsewhere does not restore -- so re-source the whole
# R/ tree. No-op under a real package build (source tree absent).
local({
  r_dir <- file.path("..", "..", "R")
  if (file.exists(file.path(r_dir, "seas_test.R"))) {
    for (f in list.files(r_dir, pattern = "\\.R$", full.names = TRUE)) {
      sys.source(f, envir = globalenv())
    }
  }
})

# Small shared simulation spec (fixed seed). J = 4 components, 4 sample years,
# 6 history years -> 10 years = 40 quarterly reference periods per component.
.vs_spec <- function(seed = 101L) {
  vintagesim::new_simulation_spec(
    J = 4, sample_years = 4, history_years = 6, seed = seed
  )
}

# The final-vintage snapshot's seasonal component series, index-ordered.
.seasonal_component <- function(snap, sid) {
  df <- tibble::as_tibble(snap)
  df <- dplyr::filter(df, process == "seasonal",
                      series_type == "component", series_id == sid)
  dplyr::arrange(df, period)
}

# That component's TRUE nonseasonal counterfactual from the latent truth.
.nonseasonal_truth <- function(env, sid) {
  df <- tibble::as_tibble(env$truth)
  df <- dplyr::filter(df, process == "nonseasonal",
                      series_type == "component", series_id == sid)
  dplyr::arrange(df, period)
}

test_that("tier 4: minimum-phase adjustment beats no adjustment vs the nonseasonal truth", {
  skip_on_cran()
  skip_if_not_installed("vintagesim")

  spec <- .vs_spec(101L)
  env  <- vintagesim::simulate_vintage_environment(spec)
  v    <- vintagesim::vintage_range(spec)
  snap <- vintagesim::get_snapshot(env, v[length(v)])

  sid <- "C01"
  sx  <- .seasonal_component(snap, sid)
  tx  <- .nonseasonal_truth(env, sid)
  common <- dplyr::intersect(sx$period, tx$period)
  xv  <- sx$value[match(common, sx$period)]        # observed seasonal series
  txv <- tx$value[match(common, tx$period)]        # nonseasonal truth (xNS)

  adj <- suppressMessages(
    seas_adjust(stats::ts(xv, frequency = 4), phase_rule = "minimum"))
  av  <- as.numeric(adj$adjusted)
  expect_true(adj$test$decision)                   # detected seasonal

  # xNS is NOT recoverable in LEVELS: the shared-innovation seasonal filter has
  # the same gain at frequency zero as at the seasonal frequencies, so no
  # adjustment converges to xNS in levels even in population (vintagesim README,
  # "Key design facts"). The honest tier-4 comparison is therefore on demeaned
  # series -- exactly the metric the package's own tier-3 test-worlds.R uses.
  dm <- function(z) z - mean(z)
  rmse_adj_dm  <- sqrt(mean((dm(av) - dm(txv))^2))
  rmse_none_dm <- sqrt(mean((dm(xv) - dm(txv))^2))
  rmse_adj_lvl  <- sqrt(mean((av - txv)^2))        # reported for honesty
  rmse_none_lvl <- sqrt(mean((xv - txv)^2))
  message(sprintf(
    "[integration a] %s vs xNS: demeaned RMSE adj=%.3f none=%.3f | level RMSE adj=%.3f none=%.3f",
    sid, rmse_adj_dm, rmse_none_dm, rmse_adj_lvl, rmse_none_lvl))

  # tier-4 acceptance: adjustment reduces the (demeaned) distance to truth.
  expect_lt(rmse_adj_dm, rmse_none_dm)
})

test_that("post-adjustment seasonality test fails to reject in most components (self-check)", {
  skip_on_cran()
  skip_if_not_installed("vintagesim")

  # Pool the four components across a small fixed set of seeds (20 series). A
  # single seed's four components is a noisy estimate of "most" (2-4 of 4);
  # pooling gives a stable clean fraction. All seeds and the adjustment path are
  # deterministic (no RNG in adjustment), so the pooled count is reproducible.
  ids   <- c("C01", "C02", "C03", "C04")
  seeds <- c(101L, 202L, 303L, 404L, 505L)
  per_seed <- integer(0)
  pvals    <- numeric(0)
  for (sd in seeds) {
    spec <- .vs_spec(sd)
    env  <- vintagesim::simulate_vintage_environment(spec)
    v    <- vintagesim::vintage_range(spec)
    snap <- vintagesim::get_snapshot(env, v[length(v)])
    ps <- vapply(ids, function(sid) {
      sx  <- .seasonal_component(snap, sid)
      adj <- suppressMessages(
        seas_adjust(stats::ts(sx$value, frequency = 4), phase_rule = "minimum"))
      seas_test(as.numeric(adj$adjusted), frequency = 4)$evt$p
    }, numeric(1))
    pvals    <- c(pvals, ps)
    per_seed <- c(per_seed, sum(ps > 0.05))
  }

  n_clean <- sum(pvals > 0.05)
  n_tot   <- length(pvals)
  message(sprintf(
    "[integration b] post-adjustment p > 0.05 in %d/%d pooled series (%.0f%%); per-seed clean: %s",
    n_clean, n_tot, 100 * n_clean / n_tot, paste(per_seed, collapse = ", ")))

  # "p > alpha in most cases": assert a strict majority of the pooled series are
  # left non-seasonal by the adjustment (actual at these seeds: ~70%).
  expect_gt(n_clean, n_tot / 2)
})

test_that("keyed-tsibble path maps over 3 components and reconstructs each key", {
  skip_on_cran()
  skip_if_not_installed("vintagesim")

  spec <- .vs_spec(101L)
  env  <- vintagesim::simulate_vintage_environment(spec)
  v    <- vintagesim::vintage_range(spec)
  snap <- vintagesim::get_snapshot(env, v[length(v)])

  keys <- c("C01", "C02", "C03")
  tb <- tibble::as_tibble(snap)
  tb <- dplyr::filter(tb, process == "seasonal", series_type == "component",
                      series_id %in% keys)
  tb <- dplyr::select(tb, series_id, period, value)
  tb <- tsibble::as_tsibble(tb, key = series_id, index = period)

  coll <- suppressMessages(seas_adjust(tb, phase_rule = "minimum"))

  # collection shape: one seas_sa per key
  expect_s3_class(coll, "seas_collection")
  expect_equal(attr(coll, "type"), "sa")
  expect_equal(nrow(coll), length(keys))

  adj_tb <- adjusted(coll)
  sea_tb <- seasonal(coll)
  expect_true(all(c("series_id", "adjusted") %in% names(adj_tb)))

  for (i in seq_len(nrow(coll))) {
    key <- coll$series_id[i]
    a <- adj_tb$adjusted[[i]]
    s <- sea_tb$seasonal[[i]]
    expect_s3_class(a, "tbl_ts")                   # class preserved: tsibble out

    av <- as.numeric(a[[tsibble::measured_vars(a)]])
    sv <- as.numeric(s[[tsibble::measured_vars(s)]])
    xv <- dplyr::arrange(dplyr::filter(tibble::as_tibble(tb), series_id == key),
                         period)$value
    # adjusted + seasonal == x, per key (surgery identity)
    expect_lt(max(abs((av + sv) - xv)), 1e-8)
  }
})
