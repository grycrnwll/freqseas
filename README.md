# freqseas

<!-- badges: start -->
[![R-CMD-check](https://github.com/OWNER/freqseas/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/OWNER/freqseas/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

freqseas is a frequency-domain toolkit for **detect**ing, **specify**ing,
**identify**ing, and **operat**ing on seasonality in time series, organized
around the workflow's own epistemic hierarchy rather than treating "adjust
the series" as one undifferentiated step. An extreme-value test detects
excess mass at seasonal frequencies against a white-noise null (testable);
shoulder and block-phase tests distinguish a fixed deterministic pattern
from a stochastic band (testable); a *declared* phase rule fixes the
filter-versus-components identification, because that choice is not testable
at second order and freqseas refuses to pretend otherwise; and
seasonal-subseries-inspired gain surgery, with the phase rule's implied
minimum-phase or zero-phase correction and exact whitener recoloring,
performs the adjustment. Inputs may be `ts`, numeric with a stated
`frequency`, or keyed `tsibble` panels, mapped over independently per key.

## Install

freqseas is not yet on CRAN.

```r
# from GitHub
remotes::install_github("OWNER/freqseas")

# from a local clone
remotes::install_local("path/to/freqseas")
```

## Quick start

```r
library(freqseas)

fit <- seas_adjust(UKgas)
fit
#> SSI seasonal adjustment (band)
#>  detection p = 1.41e-11 | spec: band (shoulder p = 1.41e-11, phase R = 0.96)
#>  phase rule: minimum  [declared identification]
#>  M = 5 (suggest_M, offset-free) | whitener: d = 0, AR(3) on M0 ordinates
#>  post-adjustment detection p = 0.677

glance(fit)
#> # A tibble: 1 x 9
#>   decision  p.value spec      M     d ar_order phase_rule post_p max_gain
#>   <lgl>       <dbl> <chr> <int> <int>    <int> <chr>       <dbl>    <dbl>
#> 1 TRUE     1.41e-11 band      5     0        3 minimum     0.677     17.3
```

`seas_adjust(x)` is the one-call wrapper for the whole
detect-specify-identify-operate workflow; `seas_test()` + `seas_ssi()` is the
equivalent two-stage path when you want to inspect the test before
committing to an adjustment. `adjusted()`, `seasonal()`, `whitener()`, and
`decision()` pull the pieces back out; `tidy()`/`glance()`/`augment()`
(from the `generics` package) give tibble-shaped views; `plot()`/`autoplot()`
cover both base graphics and ggplot2.

### `phase_rule` is declared, not defaulted (at the low level)

A periodogram only ever identifies a filter's *gain*; it says nothing about
its *phase*, and infinitely many filters share a gain. freqseas resolves
that with two named, mutually exclusive identification choices:

* **`"minimum"`** -- filter world: seasonality distorts the series' own
  innovations (shared-innovation filter). Removal is deconvolution: divide
  the gain **and** de-rotate by the implied minimum phase.
* **`"zero"`** -- components world: seasonality is an independent, additive
  component. Removal is Wiener-style shrinkage: divide the gain, leave phase
  alone. This is the classic published-SSI rule.

`seas_adjust()` defaults `phase_rule` to `"minimum"` so the one-call path is
usable out of the box -- but the lower-level `seas_ssi()` has **no default**
and raises a teaching error if you omit it, because filter-versus-components
is a modelling commitment the data cannot decide on their own. Every
print/summary tags the choice `[declared identification]` for exactly this
reason. See `vignette("two-worlds")` for the full identification story,
including a verified demonstration of what it costs to get it backwards.

## Vignettes

* [`vignette("getting-started")`](vignettes/getting-started.Rmd) -- the
  one-call and two-stage paths, accessors, broom generics, and a keyed-panel
  teaser, on `UKgas`.
* [`vignette("two-worlds")`](vignettes/two-worlds.Rmd) -- the identification
  story: why gain alone underdetermines the filter, minimum versus zero
  phase, a verified miniature demonstration, the deterministic-line case,
  and the comb diagnostic.
* [`vignette("the-machinery")`](vignettes/the-machinery.Rmd) -- bins, donors,
  the EVT detection statistic, the M0/guard-restricted de-biased-Whittle
  whitener, `suggest_M()`, `min_phase()`, the surgery, and recoloring.
* [`vignette("panels-and-vintages")`](vignettes/panels-and-vintages.Rmd) --
  keyed tsibbles end to end, decision persistence under concurrent
  adjustment, a `vintagesim` round trip, and `benchmark_totals()`.

## Status

Pre-CRAN, under active development. The R-CMD-check badge above will go live
once the CI workflow (`.github/workflows/R-CMD-check.yaml`) and the GitHub
remote are in place.
