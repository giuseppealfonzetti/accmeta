th <- c(2.94, -2.2, -0.4, 0.0953, 0.4, -0.5108, 0.3, 0.2, -0.6931)

test_that("plot draws every panel and returns the fit", {
  set.seed(1)
  d <- set_meta_data(sim_data(15, th, sample(40:200, 15, TRUE)), CC = 0.5)
  f <- fit_ib(d, CONTROL = set_ib_control(H = 20), SEEDS = 1:20)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  for (w in list("convergence", "contraction", "path", NULL)) {
    out <- withVisible(if (is.null(w)) plot(f) else plot(f, WHICH = w))
    expect_false(out$visible)
    expect_identical(out$value, f)
  }
  expect_error(plot(f, WHICH = "other"))
})

test_that("plot handles boundary stops and bivariate fits", {
  set.seed(3)
  d <- set_meta_data(sim_data(15, th, sample(40:200, 15, TRUE)), CC = 0.5)
  b <- fit_ib(d, CONTROL = set_ib_control(H = 20), SEEDS = 1:20)
  expect_identical(b$STOP, "boundary")
  v <- fit_ib(
    d,
    CONTROL = set_ib_control(H = 20, MAX_ITER = 3, TERMINATION = "max_iter"),
    SEEDS = 1:20,
    MODEL = "bivariate"
  )
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_no_error(plot(b))
  expect_no_error(plot(v))
})
