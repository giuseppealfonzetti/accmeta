th <- c(2.94, -2.2, -0.4, 0.0953, 0.4, -0.5108, 0.3, 0.2, -0.6931)

test_that("the fits carry a class", {
  set.seed(1)
  d <- set_meta_data(sim_data(10, th, rep(100, 10)), CC = 0.5)
  expect_s3_class(fit_tlmm(d), c("accmeta_tlmm", "accmeta_fit"), exact = TRUE)
  expect_s3_class(fit_tglmm(d), c("accmeta_tglmm", "accmeta_fit"), exact = TRUE)
  expect_s3_class(
    fit_ib(d, H = 10, MAX_ITER = 2, SEEDS = 1:10),
    c("accmeta_ib", "accmeta_fit"),
    exact = TRUE
  )
})

test_that("coef unpacks the working vector", {
  set.seed(1)
  d <- set_meta_data(sim_data(10, th, rep(100, 10)), CC = 0.5)
  f <- fit_tlmm(d)
  expect_identical(coef(f), theta2list(f$THETA))
  expect_named(coef(f)$MU, c("eta", "xi", "gamma"))
  expect_identical(dim(coef(f)$SIGMA), c(3L, 3L))

  g <- fit_tglmm(d)
  expect_identical(coef(g), theta2list(g$THETA))

  # the corrected estimate, not the auxiliary
  b <- fit_ib(d, H = 10, MAX_ITER = 2, SEEDS = 1:10)
  expect_identical(coef(b), theta2list(b$THETA))
  expect_false(isTRUE(all.equal(coef(b)$MU, theta2list(b$PI_HAT)$MU)))
})

test_that("printing is short and invisible", {
  set.seed(1)
  d <- set_meta_data(sim_data(10, th, rep(100, 10)), CC = 0.5)
  f <- fit_tlmm(d)
  out <- capture.output(res <- print(f))
  expect_lt(length(out), 6)
  expect_identical(res, f)
  expect_match(out[1], "accmeta_tlmm")
  # the TMB object stays hidden
  expect_false(any(grepl("OBJ|function", out)))

  b <- fit_ib(d, H = 10, MAX_ITER = 2, SEEDS = 1:10)
  ib_out <- capture.output(print(b))
  expect_match(ib_out[1], "accmeta_ib")
  expect_match(ib_out[1], b$STOP)
})
