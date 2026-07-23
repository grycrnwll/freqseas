# test-worlds.R -- TIER 3: the package's central scientific regression.
#
# Ports the logic of scripts/phase_worlds.R (notes/2026-07-22_ssi_phase_
# surgery.md sec.4) but THROUGH THE PACKAGE API: estimated donor-level gain and
# the full whiten -> surgery -> recolor path, not the oracle-gain circular
# filter of the original script. The claim is the two-worlds RMSE ordering:
#   * filter world  (shared innovations): phase_rule = "minimum" beats "zero";
#   * components world (independent seasonal): "zero" beats "minimum";
# and both rules beat no adjustment against the known counterfactual xNS.
#
# Margins are smaller than the session's oracle table because the gain is
# ESTIMATED at n = 120; the thresholds below (win rates 0.70 / 0.60, strict
# mean orderings) were chosen after running and hold with margin. Fully
# deterministic: fixed seeds, and the adjustment path uses no RNG.

# --- assembly-time integration guard (see test-seas-ssi.R for the rationale):
# restore the real single-series layer in case the concurrently-built
# test-collection.R clobbered the globals with `<<-` stubs. No-op under a real
# package build.
local({
  r_dir <- file.path("..", "..", "R")
  if (file.exists(file.path(r_dir, "seas_test.R"))) {
    for (f in c("classes.R", "seas_test.R", "seas_ssi.R", "seas_adjust.R",
                "methods-print.R")) {
      sys.source(file.path(r_dir, f), envir = globalenv())
    }
  }
})

# Two-worlds DGP (self-contained port of sim_worlds() from phase_worlds.R):
#   ns : AR(1) nonseasonal counterfactual (shared eps stream)
#   f  : filter world  x = H(B) xNS, one innovation stream
#   c  : components world x = xNS + s, s an independent seasonal AR
sim_worlds <- function(n, phi = 0.5, Phi = 0.9, s = 4L, sig_eta = 1,
                       spin = 400L) {
  ntot <- spin + n
  keep <- (spin + 1L):ntot
  eps  <- stats::rnorm(ntot)
  eta  <- stats::rnorm(ntot, sd = sig_eta)
  d_ns <- numeric(ntot)
  for (t in 2:ntot) d_ns[t] <- phi * d_ns[t - 1L] + eps[t]
  d_f <- numeric(ntot)
  for (t in (s + 2L):ntot) {
    d_f[t] <- phi * d_f[t - 1L] + Phi * d_f[t - s] -
      phi * Phi * d_f[t - s - 1L] + eps[t]
  }
  sc <- numeric(ntot)
  for (t in (s + 1L):ntot) sc[t] <- Phi * sc[t - s] + eta[t]
  xNS <- d_ns[keep]; xF <- d_f[keep]; xC <- (d_ns + sc)[keep]
  list(ns = xNS - mean(xNS), f = xF - mean(xF), c = xC - mean(xC))
}

# RMSE vs the counterfactual on a trimmed window (both demeaned; the phase rule
# is level-invariant, and recolor anchors the first few points to x).
rmse_vs <- function(a, b, trim) {
  a <- a - mean(a); b <- b - mean(b)
  sqrt(mean((a - b)[trim]^2))
}

NREP <- 30L
TRIM <- 13:120

test_that("filter world: minimum-phase beats zero-phase (win rate and mean RMSE)", {
  Fmin <- Fzero <- Fnone <- numeric(NREP)
  for (r in seq_len(NREP)) {
    set.seed(3000L + r)
    p     <- sim_worlds(120)
    amin  <- suppressMessages(seas_adjust(p$f, frequency = 4, phase_rule = "minimum"))
    azero <- suppressMessages(seas_adjust(p$f, frequency = 4, phase_rule = "zero"))
    Fmin[r]  <- rmse_vs(as.numeric(amin$adjusted),  p$ns, TRIM)
    Fzero[r] <- rmse_vs(as.numeric(azero$adjusted), p$ns, TRIM)
    Fnone[r] <- rmse_vs(p$f, p$ns, TRIM)
  }

  win_min <- mean(Fmin < Fzero)
  expect_gte(win_min, 0.70)                     # minimum wins >= 70% of reps
  expect_lt(mean(Fmin), mean(Fzero))            # strict mean ordering
  expect_lt(mean(Fmin),  mean(Fnone))           # both rules beat no adjustment
  expect_lt(mean(Fzero), mean(Fnone))
})

test_that("components world: zero-phase beats minimum-phase (win rate and mean RMSE)", {
  Cmin <- Czero <- Cnone <- numeric(NREP)
  for (r in seq_len(NREP)) {
    set.seed(3000L + r)
    p     <- sim_worlds(120)
    amin  <- suppressMessages(seas_adjust(p$c, frequency = 4, phase_rule = "minimum"))
    azero <- suppressMessages(seas_adjust(p$c, frequency = 4, phase_rule = "zero"))
    Cmin[r]  <- rmse_vs(as.numeric(amin$adjusted),  p$ns, TRIM)
    Czero[r] <- rmse_vs(as.numeric(azero$adjusted), p$ns, TRIM)
    Cnone[r] <- rmse_vs(p$c, p$ns, TRIM)
  }

  win_zero <- mean(Czero < Cmin)
  expect_gte(win_zero, 0.60)                    # zero wins >= 60% of reps
  expect_lt(mean(Czero), mean(Cmin))            # strict mean ordering
  expect_lt(mean(Czero), mean(Cnone))           # both rules beat no adjustment
  expect_lt(mean(Cmin),  mean(Cnone))
})
