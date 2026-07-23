# data-raw/make_ssi_examples.R
# ---------------------------------------------------------------------------
# Regenerate the bundled `ssi_examples` known-truth dataset (data/ssi_examples.rda).
#
# `ssi_examples` is a 3-row tibble, one row per identification "world", each
# carrying a simulated quarterly series `x` (n = 120) together with its TRUE
# nonseasonal counterfactual `truth` -- a benchmark no real dataset can supply.
# The three worlds are the ones the package's identification story turns on
# (see notes/2026-07-22_ssi_phase_surgery.md sec.4 and
# notes/2026-07-22_freqseas_architecture.md sec.7):
#
#   "filter"     -- shared-innovation (deconvolution) world.
#   "components" -- independent-seasonal (Wiener-shrinkage) world.
#   "line"       -- deterministic (fixed-pattern) seasonality.
#
# This script is SEEDED and RERUNNABLE: run it with
#   Rscript data-raw/make_ssi_examples.R
# from the package root and it reproduces data/ssi_examples.rda bit-for-bit,
# printing the defining-identity cross-checks for each world (all ~1e-15).
#
# DGP calibration (the session's): nonseasonal AR(1) coefficient phi = 0.5;
# seasonal coefficient Phi = 0.9; quarterly period s = 4; deterministic-line
# amplitudes A12 = 1.5 (annual harmonic, cos + sin) and A3 = 2.0 (semiannual /
# Nyquist harmonic). All innovation streams are standard normal (sd = 1).
# Every series is generated over a spin-up of 400 quarters (discarded) followed
# by the 120 retained quarters, so the retained window is free of the zero
# initial conditions and the defining identities hold to machine precision.
# ---------------------------------------------------------------------------

# No package dependencies beyond `tibble` (used only to assemble the result and
# to save it). Base R for all simulation.
suppressWarnings(suppressMessages(library(tibble)))

## --- fixed generation constants --------------------------------------------
SEED  <- 20260722L    # master seed; all three worlds draw from this one stream
N     <- 120L         # retained quarterly observations
FREQ  <- 4L           # quarterly
SPIN  <- 400L         # spin-up quarters (discarded)
PHI   <- 0.5          # nonseasonal AR(1) coefficient (phi)
BIGPHI<- 0.9          # seasonal AR coefficient (Phi)
S     <- 4L           # seasonal period
A12   <- 1.5          # deterministic annual-harmonic amplitude
A3    <- 2.0          # deterministic semiannual (Nyquist) amplitude
START <- c(2000L, 1L) # ts start: 2000 Q1 .. 2029 Q4

ntot <- SPIN + N
keep <- (SPIN + 1L):ntot

as_q <- function(v) stats::ts(as.numeric(v), start = START, frequency = FREQ)

set.seed(SEED)

# ===========================================================================
# World 1 -- FILTER (shared innovations, deconvolution world)
#   (1 - phi B)(1 - Phi B^s) x   = eps          (the observed seasonal series)
#   (1 - phi B)               xNS = eps          (the nonseasonal counterfactual,
#                                                 SAME innovations eps)
#   Expanded recursion:
#     x_t   = phi x_{t-1} + Phi x_{t-4} - phi*Phi x_{t-5} + eps_t
#     xNS_t = phi xNS_{t-1} + eps_t
#   Defining identity (generation cross-check): xNS_t = x_t - Phi x_{t-4}
#   (apply (1 - phi B) to both definitions; the seasonal factor cancels).
# ===========================================================================
eps <- stats::rnorm(ntot)
xf  <- numeric(ntot)          # the seasonal series
xns <- numeric(ntot)          # the nonseasonal counterfactual (truth)
for (t in 2L:ntot) xns[t] <- PHI * xns[t - 1L] + eps[t]
for (t in (S + 2L):ntot) {
  xf[t] <- PHI * xf[t - 1L] + BIGPHI * xf[t - S] - PHI * BIGPHI * xf[t - S - 1L] +
    eps[t]
}
x_filter     <- xf[keep]
truth_filter <- xns[keep]
# identity check on the retained window (needs the lag-4 of the full array)
xf_lag4      <- c(rep(NA_real_, S), utils::head(xf, ntot - S))
resid_filter <- (xns - (xf - BIGPHI * xf_lag4))[keep]
params_filter <- list(
  world = "filter", phi = PHI, Phi = BIGPHI, s = S, sd_eps = 1,
  spin = SPIN, n = N, frequency = FREQ, seed = SEED, start = START,
  equation = "(1 - phi B)(1 - Phi B^s) x = eps ; (1 - phi B) xNS = eps",
  identity = "xNS_t = x_t - Phi * x_{t-4}"
)

# ===========================================================================
# World 2 -- COMPONENTS (independent seasonal, Wiener-shrinkage world)
#   truth a : (1 - phi B) a = u                  (nonseasonal AR(1), the truth)
#   seasonal s : (1 - Phi B^s) s = eta           (INDEPENDENT seasonal AR)
#   x = a + s
#   Cross-checks: x_t - truth_t = s_t exactly; (1 - Phi B^s) s recovers eta;
#   (1 - phi B) a recovers u.
# ===========================================================================
u   <- stats::rnorm(ntot)
eta <- stats::rnorm(ntot)
a   <- numeric(ntot)          # nonseasonal AR(1) part (truth)
sc  <- numeric(ntot)          # independent seasonal AR part
for (t in 2L:ntot)       a[t]  <- PHI    * a[t - 1L]  + u[t]
for (t in (S + 1L):ntot) sc[t] <- BIGPHI * sc[t - S]  + eta[t]
x_components     <- (a + sc)[keep]
truth_components <- a[keep]
sc_lag4  <- c(rep(NA_real_, S), utils::head(sc, ntot - S))
a_lag1   <- c(NA_real_, utils::head(a, ntot - 1L))
resid_components <- (x_components - truth_components) - sc[keep]   # == 0
eta_resid <- ((sc - BIGPHI * sc_lag4) - eta)[keep]                 # == 0
u_resid   <- ((a  - PHI    * a_lag1)  - u  )[keep]                 # == 0
params_components <- list(
  world = "components", phi = PHI, Phi = BIGPHI, s = S,
  sd_u = 1, sd_eta = 1, spin = SPIN, n = N, frequency = FREQ,
  seed = SEED, start = START,
  equation = "a: (1 - phi B) a = u ; s: (1 - Phi B^s) s = eta ; x = a + s",
  identity = "x_t - truth_t = s_t  (truth = a)",
  seasonal = sc[keep]   # the true (independent) seasonal component, s = x - truth
)

# ===========================================================================
# World 3 -- LINE (deterministic fixed quarterly pattern)
#   truth a : (1 - phi B) a = u                  (nonseasonal AR(1), the truth)
#   pattern p_t = A12 cos(pi t / 2) + A12 sin(pi t / 2) + A3 cos(pi t)
#              (a FIXED, non-random quarterly pattern: A12 at the annual
#               harmonic pi/2, A3 at the semiannual/Nyquist harmonic pi)
#   x = a + p
#   Cross-check: x_t - truth_t = p_t exactly (deterministic).
# ===========================================================================
u2 <- stats::rnorm(ntot)
a2 <- numeric(ntot)
for (t in 2L:ntot) a2[t] <- PHI * a2[t - 1L] + u2[t]
tt   <- seq_len(ntot)
patt <- A12 * cos(pi * tt / 2) + A12 * sin(pi * tt / 2) + A3 * cos(pi * tt)
x_line       <- (a2 + patt)[keep]
truth_line   <- a2[keep]
pattern_keep <- patt[keep]
resid_line   <- (x_line - truth_line) - pattern_keep              # == 0
params_line <- list(
  world = "line", phi = PHI, A12 = A12, A3 = A3, s = S, sd_u = 1,
  spin = SPIN, n = N, frequency = FREQ, seed = SEED, start = START,
  equation = "a: (1 - phi B) a = u ; p_t = A12 cos(pi t/2) + A12 sin(pi t/2) + A3 cos(pi t) ; x = a + p",
  identity = "x_t - truth_t = p_t  (truth = a)",
  pattern  = pattern_keep
)

# ===========================================================================
# Assemble and save
# ===========================================================================
ssi_examples <- tibble::tibble(
  world  = c("filter", "components", "line"),
  x      = list(as_q(x_filter),     as_q(x_components),     as_q(x_line)),
  truth  = list(as_q(truth_filter), as_q(truth_components), as_q(truth_line)),
  params = list(params_filter,      params_components,      params_line)
)

# --- print the defining-identity cross-checks (all should be ~1e-15) -------
cat("ssi_examples known-truth cross-checks (max abs residual on n = 120):\n")
cat(sprintf("  filter     : |xNS - (x - Phi x_{t-4})|      = %.3e\n",
            max(abs(resid_filter))))
cat(sprintf("  components  : |x - truth - s|               = %.3e\n",
            max(abs(resid_components))))
cat(sprintf("  components  : |(1 - Phi B^s) s - eta|        = %.3e\n",
            max(abs(eta_resid))))
cat(sprintf("  components  : |(1 - phi B) a - u|            = %.3e\n",
            max(abs(u_resid))))
cat(sprintf("  line        : |x - truth - pattern|          = %.3e\n",
            max(abs(resid_line))))

out_dir <- "data"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
save(ssi_examples, file = file.path(out_dir, "ssi_examples.rda"),
     compress = "xz", version = 2)
cat(sprintf("\nSaved %s (%d rows).\n",
            file.path(out_dir, "ssi_examples.rda"), nrow(ssi_examples)))
