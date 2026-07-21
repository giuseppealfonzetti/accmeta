test_that("set_prior stores what it is given", {
  p <- set_prior()
  expect_s3_class(p, "accmeta_prior")
  expect_named(p, c("DEGREES", "SCALE"))
  # the default is penalised
  expect_identical(p$DEGREES, 5)
  expect_identical(formals(set_prior)$DEGREES, 5)
  expect_identical(p$SCALE, Inf)
  # four and infinity mean flat
  expect_identical(set_prior(4)$DEGREES, 4)
  expect_identical(set_prior(DEGREES = 6)$DEGREES, 6)
  expect_identical(set_prior(SCALE = 100)$SCALE, 100)
})

test_that("set_prior rejects impossible values", {
  expect_error(set_prior(DEGREES = 3), ">= 4")
  expect_error(set_prior(DEGREES = 4.5), "round")
  expect_error(set_prior(DEGREES = c(4, 5)), "length")
  expect_error(set_prior(SCALE = 0), "> 0")
  expect_error(set_prior(SCALE = -1), "> 0")
  expect_error(set_prior(SCALE = "big"), "is.numeric")
})

test_that("the fits want a prior object", {
  th <- c(2.94, -2.2, -0.4, 0.0953, 0.4, -0.5108, 0.3, 0.2, -0.6931)
  set.seed(1)
  d <- set_meta_data(sim_data(10, th, rep(100, 10)), CC = 0.5)
  expect_error(fit_tlmm(d, PRIOR = list(DEGREES = 4, SCALE = Inf)), "accmeta_prior")
  expect_error(fit_tglmm(d, PRIOR = 0), "accmeta_prior")
  # all three share one default
  expect_identical(formals(fit_tlmm)$PRIOR, quote(set_prior()))
  expect_identical(formals(fit_tglmm)$PRIOR, quote(set_prior()))
  expect_identical(formals(fit_ib)$PRIOR, quote(set_prior()))
  # and the penalty reaches the fit
  expect_equal(fit_tlmm(d)$NLL, fit_tlmm(d, PRIOR = set_prior(5))$NLL)
  expect_false(isTRUE(all.equal(fit_tlmm(d)$NLL, fit_tlmm(d, PRIOR = set_prior(4))$NLL)))
})
