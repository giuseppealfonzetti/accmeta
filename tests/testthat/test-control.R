test_that("set_ib_control stores what it is given and fills the rest", {
  ctrl <- set_ib_control()
  expect_s3_class(ctrl, "accmeta_ib_control")
  expect_named(ctrl, c(
    "H", "MAX_ITER", "STEP", "PATIENCE", "TERMINATION", "PRECISION", "NCORES",
    "VCOV_H", "VCOV_N_NODES", "VCOV_N_ITER", "BOX_SD", "BOX_COR", "BOX_MU"
  ))
  expect_identical(ctrl$H, 100)
  expect_null(ctrl$MAX_ITER)
  expect_identical(ctrl$BOX_SD, c(0.1, 10))
  expect_identical(ctrl$BOX_COR, 0.99)
  expect_identical(ctrl$STEP, 1)
  expect_identical(ctrl$PATIENCE, 30L)
  expect_identical(ctrl$NCORES, 1L)
  expect_identical(ctrl$TERMINATION, "mahalanobis")
  expect_identical(ctrl$PRECISION, 1)
  expect_no_warning(set_ib_control(H = 10, TERMINATION = "max_iter"))
  expect_no_warning(set_ib_control(H = 10))
  expect_identical(
    set_ib_control(TERMINATION = "mahalanobis", PRECISION = 0.5)$PRECISION,
    0.5
  )
  expect_identical(set_ib_control(MAX_ITER = 5)$MAX_ITER, 5)
  expect_identical(set_ib_control(MAX_ITER = 5)$H, 100)
})

test_that("set_ib_control rejects impossible values", {
  expect_error(set_ib_control(H = 9), "H > 9")
  expect_error(set_ib_control(MAX_ITER = 0), "MAX_ITER")
  expect_error(set_ib_control(STEP = 2), "STEP <= 1")
  expect_error(set_ib_control(STEP = 0), "STEP > 0")
  expect_error(set_ib_control(PATIENCE = 0), "PATIENCE >= 1")
  expect_error(set_ib_control(TERMINATION = "nope"), "should be one of")
  expect_error(set_ib_control(TERMINATION = "hotelling"), "should be one of")
  expect_error(set_ib_control(PRECISION = 0), "PRECISION > 0")
  expect_identical(set_ib_control(H = 1000)$H, 1000)
  expect_error(set_ib_control(NCORES = "two"), "is.numeric")
  expect_error(set_ib_control(NCORES = 0), "NCORES >= 1")
  expect_error(set_ib_control(BOX_SD = c(1, 0.5)), "BOX_SD\\[1\\] < BOX_SD\\[2\\]")
  expect_error(set_ib_control(BOX_SD = c(0, 1)), "BOX_SD\\[1\\] > 0")
  expect_error(set_ib_control(BOX_COR = 1), "BOX_COR < 1")
  expect_error(set_ib_control(BOX_MU = 0), "BOX_MU > 0")
  expect_error(set_ib_control(BOUNDARY = "project"), "unused argument")
})

test_that("fit_ib wants a control object and runs to budget on unreachable PRECISION", {
  th <- c(2.94, -2.2, -0.4, 0.0953, 0.4, -0.5108, 0.3, 0.2, -0.6931)
  set.seed(1)
  d <- set_meta_data(sim_data(15, th, rep(100, 15)), CC = 0.5)
  expect_identical(formals(fit_ib)$CONTROL, quote(set_ib_control()))
  expect_error(fit_ib(d, CONTROL = list(H = 10)), "accmeta_ib_control")
  f <- fit_ib(
    d,
    CONTROL = set_ib_control(
      H = 10, MAX_ITER = 6, PRECISION = 1e-12, STEP = 0.5, PATIENCE = 6
    ),
    SEEDS = 1:10
  )
  expect_identical(f$STOP, "maxit")
  g <- fit_ib(
    d,
    CONTROL = set_ib_control(
      H = 10, MAX_ITER = 25, PRECISION = 1e-12, STEP = 0.5, PATIENCE = 2
    ),
    SEEDS = 1:10
  )
  expect_identical(g$STOP, "stall")
  expect_false(g$CONVERGED)
  expect_equal(g$THETA, g$PATH[which.min(g$PROGRESS), ])
})
