# globals.R -- register non-standard-evaluation symbols so R CMD check's
# "checking R code for possible problems" pass does not flag them as
# "no visible binding for global variable" / "no visible global function
# definition". These are the only two NSE names in the package:
#
#   .data  the rlang/ggplot2 pronoun used inside the `ggplot2::aes()` mappings
#          of the autoplot methods (R/methods-plot.R); resolved by ggplot2's
#          tidy-evaluation data mask at draw time, never as an ordinary binding.
#   :=     the rlang walrus used once in `augment.seas_sa()` (R/methods-broom.R)
#          to name a tsibble column from a run-time string (`!!idx_col := ...`).
#
# Declaring them here keeps the static-analysis pass clean without importing
# anything at package level for the NSE itself (ggplot2 stays a Suggests-only,
# check_installed()-guarded dependency). `globalVariables()` comes from utils.

#' @importFrom utils globalVariables
NULL

globalVariables(c(".data", ":="))
