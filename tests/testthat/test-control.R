test_that("set_ib_control stores what it is given and fills the rest", {
  ctrl <- set_ib_control()
  expect_s3_class(ctrl, "accmeta_ib_control")
  expect_named(ctrl, c(
    "H", "MAX_ITER", "TOL", "STEP",
    "PLATEAU", "PLATEAU_PVALUE", "PLATEAU_WINDOW"
  ))
  # defaults
  expect_identical(ctrl$H, 100)
  expect_identical(ctrl$MAX_ITER, 25)
  expect_equal(ctrl$TOL, 0.2 / sqrt(100))
  expect_identical(ctrl$STEP, 0.1)
  expect_true(ctrl$PLATEAU)
  # TOL default tracks H unless overridden
  expect_equal(set_ib_control(H = 400)$TOL, 0.2 / sqrt(400))
  # a partial call fills the missing options
  expect_identical(set_ib_control(MAX_ITER = 5)$MAX_ITER, 5)
  expect_identical(set_ib_control(MAX_ITER = 5)$H, 100)
})

test_that("set_ib_control rejects impossible values", {
  expect_error(set_ib_control(H = 1), "H >= 2")
  expect_error(set_ib_control(MAX_ITER = 0), "MAX_ITER >= 1")
  expect_error(set_ib_control(TOL = 0), "TOL > 0")
  expect_error(set_ib_control(STEP = 2), "STEP <= 1")
  expect_error(set_ib_control(STEP = 0), "STEP > 0")
  expect_error(set_ib_control(PLATEAU = "yes"), "is.logical")
  expect_error(set_ib_control(PLATEAU_PVALUE = 1), "PLATEAU_PVALUE < 1")
  expect_error(set_ib_control(PLATEAU_WINDOW = 2), "PLATEAU_WINDOW >= 3")
})

test_that("fit_ib wants a control object and honours PLATEAU = FALSE", {
  th <- c(2.94, -2.2, -0.4, 0.0953, 0.4, -0.5108, 0.3, 0.2, -0.6931)
  set.seed(1)
  d <- set_meta_data(sim_data(15, th, rep(100, 15)), CC = 0.5)
  expect_identical(formals(fit_ib)$CONTROL, quote(set_ib_control()))
  expect_error(fit_ib(d, CONTROL = list(H = 10)), "accmeta_ib_control")
  # with the plateau rule off, only tol/maxit can end the recursion
  f <- fit_ib(
    d,
    CONTROL = set_ib_control(H = 10, MAX_ITER = 25, TOL = 1e-12, PLATEAU = FALSE),
    SEEDS = 1:10
  )
  expect_identical(f$STOP, "maxit")
})
