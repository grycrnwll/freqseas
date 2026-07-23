# Tests for R/collection.R: the keyed-tsibble mapping layer (seas_test.tbl_ts,
# seas_ssi.seas_collection, seas_adjust.tbl_ts, print.seas_collection, and the
# seas_collection accessor methods).
#
# WAVE-3 CLEANUP: the obsolete `<<-` stubs that shadowed the real
# seas_test()/seas_ssi()/seas_adjust() generics have been REMOVED so these tests
# exercise the real single-series pipeline. The value-specific assertions that
# were coupled to the old stub (adjusted == x - 0.1, spec == "line",
# post_evt = list(p = ...)) are replaced by STRUCTURAL assertions on
# collection.R's actual job: per-key extraction, canonical key ordering,
# re-clothing into tsibbles, and error propagation. Fixtures are now real
# filter-world seasonal series long enough (60 obs/key) that the SA branch's
# donor pool is well populated.

# --- assembly-time integration guard (mirrors test-seas-test.R): source the
# real package layer so the file also runs in isolation (test_file), not only
# under the verification harness that pre-sources R/. No-op under a real build.
local({
  r_dir <- file.path("..", "..", "R")
  if (file.exists(file.path(r_dir, "collection.R"))) {
    for (f in list.files(r_dir, pattern = "\\.R$", full.names = TRUE)) {
      sys.source(f, envir = globalenv())
    }
  }
})

# ---- real seasonal fixture DGP (filter-world quarterly, band spectrum) -----
# Deterministic (seeded internally); the same generator test-seas-test.R uses,
# known to be detected seasonal and specified "band".
sim_seasonal_q <- function(n, seed, Phi = 0.9, phi = 0.5, s = 4L, spin = 300L) {
  set.seed(seed)
  ntot <- spin + n
  eps  <- stats::rnorm(ntot)
  z    <- numeric(ntot)
  for (t in (s + 2L):ntot) {
    z[t] <- phi * z[t - 1L] + Phi * z[t - s] - phi * Phi * z[t - s - 1L] + eps[t]
  }
  out <- z[(spin + 1L):ntot]
  out - mean(out)
}

# fixed per-key seeds so tests can regenerate the exact series they assert on
SEED_UNKEYED <- 11L
SEED_A <- 21L; SEED_B <- 22L; SEED_C <- 23L
NQ <- 60L

make_unkeyed_q <- function() {
  tsibble::tsibble(
    qtr   = tsibble::yearquarter("2000 Q1") + 0:(NQ - 1L),
    value = sim_seasonal_q(NQ, SEED_UNKEYED),
    index = qtr
  )
}

make_keyed_q <- function() {
  # insertion order c, a, b -> tsibble sorts keys to a, b, c on construction
  vc <- sim_seasonal_q(NQ, SEED_C)
  va <- sim_seasonal_q(NQ, SEED_A)
  vb <- sim_seasonal_q(NQ, SEED_B)
  tsibble::tsibble(
    qtr   = rep(tsibble::yearquarter("2000 Q1") + 0:(NQ - 1L), 3),
    grp   = rep(c("c", "a", "b"), each = NQ),
    value = c(vc, va, vb),
    key   = grp,
    index = qtr
  )
}

# ---- seas_test.tbl_ts: bare vs. collection routing -------------------

test_that("seas_test.tbl_ts on an unkeyed tsibble returns a bare seas_test", {
  fit <- seas_test(make_unkeyed_q())

  expect_s3_class(fit, "seas_test")
  expect_false(inherits(fit, "seas_collection"))
  expect_equal(fit$input_class, "tbl_ts")
  expect_s3_class(fit$template, "tbl_ts")
  # x carried through as the plain numeric series (not eagerly rebuilt)
  expect_equal(fit$x, sim_seasonal_q(NQ, SEED_UNKEYED))
})

test_that("seas_test.tbl_ts on a keyed tsibble returns a 3-row seas_collection", {
  coll <- seas_test(make_keyed_q())

  expect_s3_class(coll, "seas_collection")
  expect_equal(nrow(coll), 3L)
  expect_equal(attr(coll, "type"), "test")
  expect_true(all(vapply(coll$fit, inherits, logical(1), what = "seas_test")))

  # tsibble's canonical (sorted) key order, not insertion order
  expect_equal(coll$grp, c("a", "b", "c"))

  # each fit carries its per-key numeric series; grp "a" is the second block
  # written (seed A) before tsibble's construction-time key sort.
  expect_equal(coll$fit[[1]]$x, sim_seasonal_q(NQ, SEED_A))
})

test_that("seas_collection type attribute is correct for test vs sa", {
  coll_test <- seas_test(make_keyed_q())
  coll_sa   <- seas_ssi(coll_test, phase_rule = "minimum")

  expect_equal(attr(coll_test, "type"), "test")
  expect_equal(attr(coll_sa, "type"), "sa")
})

# ---- print.seas_collection (snapshot-free: grepl checks) ------------------

test_that("print.seas_collection renders header and per-key content", {
  coll <- seas_test(make_keyed_q())
  out <- paste(capture.output(print(coll)), collapse = "\n")

  expect_true(grepl("seas_collection", out))
  expect_true(grepl("3 series", out))
  expect_true(grepl("test", out))
  expect_true(grepl("seasonal", out))     # decision column ("seasonal"/"not seasonal")
  expect_true(grepl("spec", out))         # spec column present
  expect_true(grepl("a", out))
  expect_true(grepl("b", out))
  expect_true(grepl("c", out))
})

test_that("print.seas_collection includes phase_rule for a sa-type collection", {
  coll_sa <- seas_ssi(seas_test(make_keyed_q()), phase_rule = "zero")
  out <- paste(capture.output(print(coll_sa)), collapse = "\n")

  expect_true(grepl("sa", out))
  expect_true(grepl("zero", out))
})

# ---- accessors on collections ------------------------------------------

test_that("decision() on a test-type collection returns a keyed tibble", {
  coll <- seas_test(make_keyed_q())
  d <- decision(coll)

  expect_s3_class(d, "tbl_df")
  expect_equal(nrow(d), 3L)
  expect_true(all(c("grp", "seasonal", "p", "spec", "alpha") %in% names(d)))
  expect_equal(d$grp, c("a", "b", "c"))
  expect_type(d$seasonal, "logical")
  # strong filter-world fixtures -> all detected seasonal
  expect_true(all(d$seasonal))
})

test_that("whitener() on a collection returns a keyed list-column", {
  coll <- seas_test(make_keyed_q())
  w <- whitener(coll)

  expect_equal(nrow(w), 3L)
  expect_true(all(c("grp", "whitener") %in% names(w)))
  expect_type(w$whitener, "list")
  # real whitener record carries a human-readable d_reason string
  expect_true(is.character(w$whitener[[1]]$d_reason))
  expect_true(nzchar(w$whitener[[1]]$d_reason))
})

test_that("adjusted()/seasonal() on a test-type collection error per-key", {
  coll <- seas_test(make_keyed_q())
  expect_error(adjusted(coll), "seas_ssi|seas_adjust")
  expect_error(seasonal(coll), "seas_ssi|seas_adjust")
})

test_that("adjusted()/seasonal() on a sa-type collection return re-clothed tsibbles", {
  coll_sa <- seas_ssi(seas_test(make_keyed_q()), phase_rule = "minimum")
  adj <- adjusted(coll_sa)
  sea <- seasonal(coll_sa)

  expect_equal(nrow(adj), 3L)
  expect_true(all(c("grp", "adjusted") %in% names(adj)))
  # fs_reclothe_sa() (and/or .reclothe()'s tbl_ts branch) must rewrap each key's
  # plain-numeric output back into a tbl_ts using the template carried on each
  # fit -- the point of the "tbl_ts" compensation.
  expect_true(all(vapply(adj$adjusted, inherits, logical(1), what = "tbl_ts")))
  expect_true(all(vapply(sea$seasonal, inherits, logical(1), what = "tbl_ts")))

  # adjusted + seasonal must reconstruct the original per-key series exactly
  # (seasonal is defined as x - adjusted).
  ai      <- which(coll_sa$grp == "a")
  a_adj   <- adj[adj$grp == "a", ]$adjusted[[1]]$value
  a_sea   <- sea[sea$grp == "a", ]$seasonal[[1]]$value
  a_fit_x <- coll_sa$fit[[ai]]$test$x
  expect_equal(a_adj + a_sea, a_fit_x)
})

# ---- seas_ssi.seas_collection -------------------------------------------

test_that("seas_ssi.seas_collection requires phase_rule and propagates the error", {
  coll <- seas_test(make_keyed_q())
  expect_error(seas_ssi(coll), "phase_rule")
})

test_that("seas_ssi.seas_collection maps seas_ssi() over every row", {
  coll    <- seas_test(make_keyed_q())
  coll_sa <- seas_ssi(coll, phase_rule = "minimum")

  expect_s3_class(coll_sa, "seas_collection")
  expect_equal(nrow(coll_sa), 3L)
  expect_true(all(vapply(coll_sa$fit, inherits, logical(1), what = "seas_sa")))
  expect_true(all(vapply(coll_sa$fit, function(f) f$phase_rule == "minimum", logical(1))))
})

# ---- seas_adjust.tbl_ts (one-call wrapper) -------------------------------

test_that("seas_adjust.tbl_ts on an unkeyed tsibble returns a re-clothed seas_sa", {
  sa <- seas_adjust(make_unkeyed_q(), phase_rule = "zero")

  expect_s3_class(sa, "seas_sa")
  expect_false(inherits(sa, "seas_collection"))
  expect_equal(sa$phase_rule, "zero")
  # bare-series re-clothing: .reclothe()'s tbl_ts branch (and the belt-and-braces
  # fs_reclothe_sa() in seas_adjust.tbl_ts) return a tbl_ts, not plain numeric.
  expect_s3_class(sa$adjusted, "tbl_ts")
  expect_s3_class(sa$seasonal, "tbl_ts")
  expect_equal(nrow(sa$adjusted), NQ)
  # adjusted + seasonal reconstructs the input series
  expect_equal(sa$adjusted$value + sa$seasonal$value, sim_seasonal_q(NQ, SEED_UNKEYED))
})

test_that("seas_adjust.tbl_ts on a keyed tsibble returns a sa-type seas_collection", {
  coll_sa <- seas_adjust(make_keyed_q(), phase_rule = "minimum")

  expect_s3_class(coll_sa, "seas_collection")
  expect_equal(attr(coll_sa, "type"), "sa")
  expect_equal(nrow(coll_sa), 3L)
})

# ---- new_seas_collection validation --------------------------------------

test_that("new_seas_collection requires a `fit` list-column", {
  expect_error(
    new_seas_collection(tibble::tibble(grp = c("a", "b")), type = "test"),
    "fit"
  )
})

test_that("new_seas_collection rejects a non-data-frame `rows`", {
  expect_error(new_seas_collection(list(fit = list(1)), type = "test"), "fit")
})
