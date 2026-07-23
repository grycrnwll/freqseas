# Tests for R/accessors.R: the four exported generics (adjusted(), seasonal(),
# whitener(), decision()) and their seas_test/seas_sa methods.
#
# Fixtures are hand-built against the PINNED CLASS CONTRACT (see the
# architecture note and the build brief) rather than constructed via the
# sibling's seas_test()/seas_ssi(), so these tests carry no dependency on
# R/classes.R, R/seas_test.R, or R/seas_ssi.R.

make_seas_test <- function(p = 0.01, alpha = 0.05, label = "line") {
  structure(
    list(
      x           = as.numeric(1:24),
      input_class = "numeric",
      N           = 12,
      P           = 1,
      M           = 5,
      alpha       = alpha,
      evt         = list(statistic = 3.2, p = p, critical = 2.1),
      spec        = list(label = label, shoulder_p = 0.2, phase_R = 0.8),
      whitener    = list(
        d        = 0,
        d_reason = "not near-integrated",
        p        = 1,
        ar       = 0.4,
        bic_path = c(-10, -12, -11)
      ),
      e    = as.numeric(1:24) - mean(1:24),
      call = quote(seas_test(x))
    ),
    class = "seas_test"
  )
}

make_seas_sa <- function(test = make_seas_test(), phase_rule = "minimum") {
  structure(
    list(
      test       = test,
      adjusted   = test$x - 0.1,
      seasonal   = rep(0.1, length(test$x)),
      Ghat       = rep(1.05, 12),
      theta      = rep(-0.2, 12),
      phase_rule = phase_rule,
      post_evt   = list(p = 0.6),
      call       = quote(seas_ssi(tst, phase_rule = "minimum"))
    ),
    class = "seas_sa"
  )
}

# ---- adjusted() ----------------------------------------------------------

test_that("adjusted() returns the adjusted field of a seas_sa", {
  sa <- make_seas_sa()
  expect_equal(adjusted(sa), sa$adjusted)
})

test_that("adjusted() on a bare seas_test errors informatively", {
  tst <- make_seas_test()
  expect_error(adjusted(tst), "seas_ssi|seas_adjust")
})

test_that("adjusted() on an unsupported class errors informatively", {
  other <- structure(list(), class = "something_else")
  expect_error(adjusted(other), "adjusted")
})

# ---- seasonal() -----------------------------------------------------------

test_that("seasonal() returns the seasonal field of a seas_sa", {
  sa <- make_seas_sa()
  expect_equal(seasonal(sa), sa$seasonal)
})

test_that("seasonal() on a bare seas_test errors informatively", {
  tst <- make_seas_test()
  expect_error(seasonal(tst), "seas_ssi|seas_adjust")
})

test_that("seasonal() on an unsupported class errors informatively", {
  other <- structure(list(), class = "something_else")
  expect_error(seasonal(other), "seasonal")
})

# ---- whitener() -----------------------------------------------------------

test_that("whitener() returns the whitener record of a seas_test", {
  tst <- make_seas_test()
  expect_equal(whitener(tst), tst$whitener)
})

test_that("whitener() on a seas_sa reaches through to test$whitener", {
  sa <- make_seas_sa()
  expect_equal(whitener(sa), sa$test$whitener)
})

test_that("whitener() on an unsupported class errors informatively", {
  other <- structure(list(), class = "something_else")
  expect_error(whitener(other), "whitener")
})

# ---- decision() -------------------------------------------------------

test_that("decision() on a seas_test returns the documented one-row tibble", {
  tst <- make_seas_test(p = 0.01, alpha = 0.05, label = "band")
  d <- decision(tst)

  expect_s3_class(d, "tbl_df")
  expect_equal(nrow(d), 1L)
  expect_named(d, c("seasonal", "p", "spec", "alpha"))
  expect_true(d$seasonal)  # p < alpha -> seasonal
  expect_equal(d$p, 0.01)
  expect_equal(d$spec, "band")
  expect_equal(d$alpha, 0.05)
})

test_that("decision() applies the p < alpha rule correctly when not seasonal", {
  tst <- make_seas_test(p = 0.5, alpha = 0.05, label = "line")
  d <- decision(tst)
  expect_false(d$seasonal)
})

test_that("decision() prefers a stored decision field over recomputing", {
  # The real seas_test constructor stores `decision` at fit time; a fixture
  # that carries a `decision` field diverging from a naive p < alpha
  # recompute (e.g. a persisted panel decision) must win over the fallback.
  tst <- make_seas_test(p = 0.5, alpha = 0.05, label = "line")
  tst$decision <- TRUE
  d <- decision(tst)
  expect_true(d$seasonal)
})

test_that("decision() on a seas_sa delegates to its embedded test", {
  tst <- make_seas_test(p = 0.02, alpha = 0.05, label = "line")
  sa  <- make_seas_sa(test = tst)
  expect_equal(decision(sa), decision(tst))
})

test_that("decision() on an unsupported class errors informatively", {
  expect_error(decision(structure(list(), class = "something_else")), "decision")
})
