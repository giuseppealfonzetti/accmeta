test_that("set_ib_control stores what it is given and fills the rest", {
  ctrl <- set_ib_control()
  expect_s3_class(ctrl, "accmeta_ib_control")
  expect_named(ctrl, c(
    "H", "MAX_ITER", "TOL", "STEP", "PATIENCE", "UPDATE", "NCORES"
  ))
  # defaults
  expect_identical(ctrl$H, 100)
  expect_identical(ctrl$MAX_ITER, 25)
  expect_identical(ctrl$TOL, 0.5)
  expect_identical(ctrl$STEP, 0.1)
  expect_identical(ctrl$PATIENCE, 5L)
  expect_identical(ctrl$UPDATE, "fixedpoint")
  expect_identical(ctrl$NCORES, 1L)
  expect_identical(set_ib_control(UPDATE = "lm")$UPDATE, "lm")
  # a partial call fills the missing options
  expect_identical(set_ib_control(MAX_ITER = 5)$MAX_ITER, 5)
  expect_identical(set_ib_control(MAX_ITER = 5)$H, 100)
})

test_that("set_ib_control rejects impossible values", {
  expect_error(set_ib_control(H = 1), "H >= 2")
  expect_error(set_ib_control(MAX_ITER = 0), "MAX_ITER >= 1")
  expect_error(set_ib_control(TOL = 0), "TOL > 0")
  expect_error(set_ib_control(TOL = 1), "TOL < 1")
  expect_error(set_ib_control(TOL = 1.5), "TOL < 1")
  expect_error(set_ib_control(STEP = 2), "STEP <= 1")
  expect_error(set_ib_control(STEP = 0), "STEP > 0")
  expect_error(set_ib_control(PATIENCE = 0), "PATIENCE >= 1")
  expect_error(set_ib_control(UPDATE = "nope"), "should be one of")
  expect_error(set_ib_control(NCORES = "two"), "is.numeric")
  expect_error(set_ib_control(NCORES = 0), "NCORES >= 1")
})

test_that("fit_ib wants a control object and runs to budget on unreachable TOL", {
  th <- c(2.94, -2.2, -0.4, 0.0953, 0.4, -0.5108, 0.3, 0.2, -0.6931)
  set.seed(1)
  d <- set_meta_data(sim_data(15, th, rep(100, 15)), CC = 0.5)
  expect_identical(formals(fit_ib)$CONTROL, quote(set_ib_control()))
  expect_error(fit_ib(d, CONTROL = list(H = 10)), "accmeta_ib_control")
  # patience disabled (>= MAX_ITER) + a near-1 alpha (unreachable) -> runs to budget
  f <- fit_ib(
    d,
    CONTROL = set_ib_control(
      H = 10, MAX_ITER = 6, TOL = 0.999, PATIENCE = 6
    ),
    SEEDS = 1:10
  )
  expect_identical(f$STOP, "maxit")
  # patience stops early and returns the best-statistic iterate
  g <- fit_ib(
    d,
    CONTROL = set_ib_control(H = 10, MAX_ITER = 25, TOL = 0.999, PATIENCE = 2),
    SEEDS = 1:10
  )
  expect_identical(g$STOP, "stall")
  expect_false(g$CONVERGED)
  expect_equal(g$THETA, g$PATH[which.min(g$PROGRESS), ])
})
