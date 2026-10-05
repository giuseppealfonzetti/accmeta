test_that("set_ib_control stores what it is given and fills the rest", {
  ctrl <- set_ib_control()
  expect_s3_class(ctrl, "accmeta_ib_control")
  expect_named(ctrl, c(
    "H", "MAX_ITER", "TOL", "STEP", "PATIENCE", "TERMINATION",
    "PRECISION", "BOOST", "BOOST_FACTOR", "MAX_H", "NCORES"
  ))
  expect_identical(ctrl$H, 100)
  expect_identical(ctrl$MAX_ITER, 100)
  expect_identical(ctrl$TOL, 0.05)
  expect_identical(ctrl$STEP, 1)
  expect_identical(ctrl$PATIENCE, 5L)
  expect_identical(ctrl$NCORES, 1L)
  expect_identical(ctrl$TERMINATION, "confidence")
  expect_identical(ctrl$PRECISION, 1)
  expect_false(ctrl$BOOST)
  expect_identical(ctrl$BOOST_FACTOR, 2)
  expect_identical(ctrl$MAX_H, 500)
  expect_no_warning(set_ib_control(PRECISION = 0.1, BOOST = TRUE))
  hot <- set_ib_control(TERMINATION = "hotelling")
  expect_identical(hot$TOL, 0.5)
  expect_false(hot$BOOST)
  expect_identical(
    set_ib_control(H = 300, TERMINATION = "confidence")$TOL,
    0.05
  )
  expect_identical(
    set_ib_control(H = 300, TERMINATION = "confidence", TOL = 0.01)$TOL,
    0.01
  )
  expect_warning(set_ib_control(PRECISION = 0.1), "cannot pass")
  expect_no_warning(set_ib_control(H = 10, TERMINATION = "max_iter"))
  expect_identical(set_ib_control(MAX_ITER = 5)$MAX_ITER, 5)
  expect_identical(set_ib_control(MAX_ITER = 5)$H, 100)
})

test_that("set_ib_control rejects impossible values", {
  expect_error(set_ib_control(H = 9), "H > 9")
  expect_error(set_ib_control(MAX_ITER = 0), "MAX_ITER >= 1")
  expect_error(set_ib_control(TOL = 0), "TOL > 0")
  expect_error(set_ib_control(TOL = 1), "TOL < 1")
  expect_error(set_ib_control(TOL = 1.5), "TOL < 1")
  expect_error(set_ib_control(STEP = 2), "STEP <= 1")
  expect_error(set_ib_control(STEP = 0), "STEP > 0")
  expect_error(set_ib_control(PATIENCE = 0), "PATIENCE >= 1")
  expect_error(set_ib_control(TERMINATION = "nope"), "should be one of")
  expect_error(set_ib_control(PRECISION = 0), "PRECISION > 0")
  expect_error(
    set_ib_control(TERMINATION = "hotelling", BOOST = TRUE),
    "TERMINATION == \"confidence\""
  )
  expect_error(
    set_ib_control(TERMINATION = "max_iter", BOOST = TRUE),
    "TERMINATION == \"confidence\""
  )
  expect_error(set_ib_control(BOOST = NA), "isTRUE")
  expect_error(set_ib_control(H = 100, MAX_H = 50), "MAX_H >= H")
  expect_error(set_ib_control(BOOST_FACTOR = 1), "BOOST_FACTOR > 1")
  expect_error(set_ib_control(NCORES = "two"), "is.numeric")
  expect_error(set_ib_control(NCORES = 0), "NCORES >= 1")
})

test_that("fit_ib wants a control object and runs to budget on unreachable TOL", {
  th <- c(2.94, -2.2, -0.4, 0.0953, 0.4, -0.5108, 0.3, 0.2, -0.6931)
  set.seed(1)
  d <- set_meta_data(sim_data(15, th, rep(100, 15)), CC = 0.5)
  expect_identical(formals(fit_ib)$CONTROL, quote(set_ib_control()))
  expect_error(fit_ib(d, CONTROL = list(H = 10)), "accmeta_ib_control")
  f <- fit_ib(
    d,
    CONTROL = set_ib_control(
      H = 10, MAX_ITER = 6, TOL = 0.999, STEP = 0.5, PATIENCE = 6,
      TERMINATION = "hotelling"
    ),
    SEEDS = 1:10
  )
  expect_identical(f$STOP, "maxit")
  g <- fit_ib(
    d,
    CONTROL = set_ib_control(
      H = 10, MAX_ITER = 25, TOL = 0.999, STEP = 0.5, PATIENCE = 2,
      TERMINATION = "hotelling"
    ),
    SEEDS = 1:10
  )
  expect_identical(g$STOP, "stall")
  expect_false(g$CONVERGED)
  expect_equal(g$THETA, g$PATH[which.min(g$PROGRESS), ])
})
