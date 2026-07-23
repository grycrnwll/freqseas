# Residual seasonality after `seas_adjust()` — diagnosis

Reported: `adj_air <- seas_adjust(log(AirPassengers))` leaves visible mass at the
seasonal frequencies. Confirmed, and the package's own self-check agrees:
`adj_air$post_evt = 0.0101 < alpha`.

Three defects, one of which is the root cause.

---

## A. Root cause: the surgery grid is not seasonally aligned

`fs_whiten()` returns conditional AR residuals of length

```
n_e = n - d - p
```

(`whiten.R:387-389`) and every downstream object — `fourier_grid()`,
`mbin_partition()`, `index_sets()`, the periodogram, and `fs_surgery()` — lives
on that length-`n_e` Fourier grid. Nothing constrains `n_e` to be a multiple of
`N`, and after differencing plus an AR fit it essentially never is.

For `log(AirPassengers)`: `n = 144`, `d = 1`, `p = 1` → **`n_e = 142`**, and the
seasonal harmonic ordinates land at

```
n_e * h / 12,  h = 1..6  =  11.833  23.667  35.5  47.333  59.167  71
```

Only `h = 6` (Nyquist) is an integer. Three consequences:

**A1 — the exact-harmonic set `H` collapses.** `tst$partition$H` is `71`, i.e.
Nyquist alone (`|H| = 1`). Nyquist is on-grid whenever `n_e` is even; for odd
`n_e`, `|H| = 0`. So the `"line"` branch of `seas_ssi()` has essentially nothing
to null. Forcing the line branch on the AirPassengers fit changes the output not
at all — `post_evt = 6.6e-10`, and the adjusted periodogram is equal to the raw
one to four decimals. The `length(H_full) == 0` guard at `seas_ssi.R:154` does
not fire in the even case, so the branch silently "adjusts" by touching only the
Nyquist ordinate.

**A2 — a deterministic line is misclassified as a band.** Off-grid, a line leaks
across the whole spectrum via the Dirichlet kernel, elevating the shoulder
ordinates; the Stage B shoulder EVT then rejects. AirPassengers shows
`phase_R = 0.996` (frozen pattern) alongside `shoulder_p = 5.9e-10` ("band") —
the two diagnostics contradict each other, and the shoulder is the artifact.

**A3 — band surgery cannot remove a leaked line.** Bin-local gain division
flattens `J1`; the `1/Δ` leakage tails outside `J1` survive and reconstruct a
sizeable fraction of the sinusoid.

### Controlled demonstration

Same DGP throughout — `3*cos(2*pi*t/12) + N(0, 0.5^2)`, a pure deterministic
line — varying only the series length:

| `n` | `n_e` | `n_e %% 12 == 0` | spec | shoulder p | phase_R | \|H\| | `post_evt` |
|---|---|---|---|---|---|---|---|
| 144 | 144 | TRUE  | line | 0.280 | 1.000 | 6 | **0.853** |
| 120 | 120 | TRUE  | line | 0.333 | 1.000 | 6 | **0.976** |
| 146 | 146 | FALSE | band | 0 | 1.000 | 1 | **3.5e-05** |
| 145 | 145 | FALSE | band | 0 | 1.000 | 0 | **6.0e-15** |

Correctness is a lottery on `n_e mod N`. Aligned: correct specification, clean
adjustment. Misaligned: wrong specification and an adjustment that leaves
detectable seasonality.

---

## B. `fs_gain()`'s log-space smoothing under-corrects sharp peaks

`gain.R:164-175` smooths the log-gain with a length-3 moving average. At an
isolated peak the two neighbours sit near log-gain 0, so the smoothed value is
roughly `log(G)/3` — the applied gain is the **cube root** of the estimated one,
and the power reduction is `G^(2/3)` instead of `G^2`.

Measured on `log(AirPassengers)`, at the ordinates carrying the seasonal peaks:

| ordinate | omega | I/donor | G raw | G smoothed | power kept (raw) | power kept (smoothed) |
|---|---|---|---|---|---|---|
| 12 | 0.531 | 292.3 | 17.10 | 3.97 | 0.0034 | **0.0635** |
| 24 | 1.062 | 262.0 | 16.19 | 6.60 | 0.0038 | 0.0229 |
| 47 | 2.080 |  80.5 |  8.97 | 3.11 | 0.0124 | **0.1031** |
| 59 | 2.611 |  29.2 |  5.41 | 1.76 | 0.0342 | **0.3247** |

At the fundamental the smoothing keeps 6.3% of the peak power where the
unsmoothed gain would keep 0.34% — a 19x under-correction.

Counterfactual, identical pipeline with the smoothing removed:

```
post_evt, smoothed gain    = 0.0101      (residual seasonality detected)
post_evt, unsmoothed gain  = 0.9994      (clean)
```

and the raw-scale periodogram at the fundamental drops from 0.1663 to 0.0915
(raw series: 0.8969).

The docstring justifies the smoothing as noise control ("so the correction is
not driven by single-ordinate periodogram noise"). That is reasonable for a
genuinely broad band and actively harmful for a narrow one. Note that fixing A
partly defuses B: on an aligned grid a line is routed to the `"line"` branch,
which does not use `fs_gain()` at all.

---

## C. Minor: the empty-`H` guard is unreachable for even `n_e`

`seas_ssi.R:153-167` returns the input unchanged with a message when
`length(H_full) == 0`. Because Nyquist is always on-grid for even `n_e`, `H` is
never empty there — it is `{Nyquist}` — so the guard is skipped and the line
branch proceeds to do nothing useful and report `spec = "line"` as though it had
worked. Whatever is done about A, this branch should report honestly when the
harmonic set it needs is not representable on the grid.

---

## Proposed fix for A — prototyped and verified on the reported series

Align the surgery grid by trimming the residual vector to a whole number of
seasonal cycles, and splice the untouched leading residuals back before
recoloring:

```r
n_use <- N * (n_e %/% N)          # largest aligned length
off   <- n_e - n_use              # leading residuals left untouched (< N)
e_use <- e[(off + 1L):n_e]
# grid / partition / periodogram / surgery all on e_use (length n_use)
estar <- c(e[seq_len(off)], estar_use)   # fs_recolor() still sees n_y - p
```

`fs_recolor()` requires `length(estar) == n_y - p` exactly (`whiten.R:450`), so
the splice is what keeps the reconstruction valid.

### Result on `log(AirPassengers)`

`n_e = 142` → `n_use = 132`, `off = 10`. The specification flips to the branch
the phase diagnostic was pointing at all along, and all six harmonics become
representable:

```
spec = line   shoulder p = 0.23   |H| = 6   detection p = 1.1e-11
```

The right acceptance metric is **local** excess, not the global median: near
`omega = 0.5` the raw series carries 0.03–0.12 at neighbouring nonseasonal
ordinates (trend), so the global median understates the local floor. Ratio of
the periodogram at each harmonic to the median over ordinates ±2..±5:

| h | omega | raw | current | unsmoothed gain | **aligned** |
|---|---|---|---|---|---|
| 1 | 0.524 | 16.7 | **4.3** | 2.3 | **1.2** |
| 2 | 1.047 |  8.9 | 0.0 | 0.3 | 0.2 |
| 3 | 1.571 |  8.8 | 2.7 | 2.1 | 1.6 |
| 4 | 2.094 | 13.3 | 3.1 | 1.3 | 1.8 |
| 5 | 2.618 |  5.5 | 2.7 | 1.3 | 1.5 |

The current path leaves a 4.3x local peak at the fundamental — that is the mass
the report is about. Alignment removes it (1.2x, i.e. flat). Total power at the
six seasonal ordinates: raw 1.180 → current 0.211 → unsmoothed 0.126 →
**aligned 0.090**.

> **Correction (post-implementation).** This section originally reported
> `post_evt` 0.0101 → 0.421. The shipped value is **0.208**, and 0.421 was an
> artifact of the prototype, not a shortfall in the fix. The prototype pushed the
> adjusted series back through the *unfixed* pipeline for its self-check, so it
> re-tested on an unaligned n_e = 142 grid; the real implementation re-tests on
> the aligned n_e = 132 grid. Different test, identical adjustment — the
> raw-scale local excess is 1.16 against the prototype's 1.2. Note what this
> incidentally shows: `post_evt` is itself grid-dependent, which is a second-order
> instance of the same defect and another reason to treat it as a diagnostic
> rather than a test.

Note the two failure modes of the current band path are visible together at the
fundamental's neighbourhood: it **under-removes** at the peak (0.1663, vs 0.0572
aligned) while **over-removing** in the shoulders (ordinate 11: 0.0332, vs 0.0998
aligned and 0.1181 raw) — bin-wide gain division strips legitimate nonseasonal
trend mass. The line branch on the aligned grid touches only the six harmonics
and leaves the rest bit-identical.

`post_evt` alone would not have caught this: dropping the gain smoothing on the
misaligned grid drives `post_evt` to 0.9994 while still leaving a 2.3x local peak
at the fundamental. The self-check is run on the re-whitened residual grid, where
the leakage is spread; it is not a proxy for raw-scale seasonal mass.

### Defect D — the Nyquist comb tooth, created BY the alignment fix

Not anticipated anywhere above; found during implementation, and worth recording
because it is the kind of thing this fix class generates.

Aligning `n_e` to a multiple of `N` forces `n_e` **even** whenever `N` is even.
So `omega = pi` now lands exactly on the Nyquist ordinate and is *always* an
exact seasonal harmonic. But `fs_gain()` returns gain 1 there: its `pgram_pos`
has length `Jmax = (n_e - 1) %/% 2`, which for even `n_e` contains no Nyquist
entry, so there is no periodogram-based estimate to work from. The band branch
therefore left the **entire comb tooth at pi untouched**.

Before alignment this was harmless — `omega = pi` was generically off-grid and
its mass sat on ordinates `J1` did cover. The fix created the exposure.

Measured: `test-seas-ssi.R:93` (band-path post-adjustment self-check) fell from
9/10 to 4/10 seeds clean, and `test-integration-vintagesim.R:127` from 12/20 to
9/20. Fixed in `seas_ssi.R`'s band branch by filling the Nyquist gain from the
raw periodogram under the same donor-ratio rule every other seasonal ordinate
gets, guarded on membership in `H` so an odd `N` is left alone.

Two things worth remembering about how this was caught. It surfaced only through
**downstream self-checks in unrelated files** — no test targeted the Nyquist
tooth, and the dedicated alignment regression suite passed throughout. And the
fill-in applies the *unsmoothed* gain `sqrt(max(I / donor, 1))`, while every
other harmonic still receives the cube-rooted smoothed gain of defect B. The
band path now runs two different gain rules at different ordinates. Defensible
(a symmetric 3-point smoother is not well defined at the edge) but it should be
a stated choice, not a side effect.

### Cost

Up to `N - 1` leading residuals go unadjusted. For AirPassengers that is the
first 12 observations (10 trimmed + the 2 `d + p` anchors) — a full year, 8% of
the sample, verified identical to the input. Trimming from the front rather than
the back keeps the recent end adjusted, which is the end that matters for most
applications, but the unadjusted head must be recorded in the whitener record and
surfaced in `print`/`summary`. Bounded and known, unlike the current behaviour.

Time-domain sanity: `sd` 0.4415 → 0.4076, `max|adjusted - raw| = 0.342`, anchors
preserved.

## Reproduction

`scratchpad/repro.R`, `diag.R`, `diag2.R`, `diag3.R` in the session scratchpad.
