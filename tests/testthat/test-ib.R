th <- c(2.94, -2.2, -0.4, 0.0953, 0.4, -0.5108, 0.3, 0.2, -0.6931)

test_that("the same seeds give the same answer", {
  set.seed(1)
  d <- set_meta_data(sim_data(15, th, rep(100, 15)), CC = 0.5)
  a <- fit_ib(d, H = 10, MAX_ITER = 2, SEEDS = 1:10)
  b <- fit_ib(d, H = 10, MAX_ITER = 2, SEEDS = 1:10)
  expect_identical(a$THETA, b$THETA)
  expect_identical(a$PATH, b$PATH)
  # seeds drive the answer
  other <- fit_ib(d, H = 10, MAX_ITER = 2, SEEDS = 101:110)
  expect_false(isTRUE(all.equal(a$THETA, other$THETA)))
})

test_that("the generator is left untouched", {
  set.seed(1)
  d <- set_meta_data(sim_data(15, th, rep(100, 15)), CC = 0.5)
  set.seed(42)
  before <- get(".Random.seed", .GlobalEnv)
  fit_ib(d, H = 10, MAX_ITER = 2)
  expect_identical(get(".Random.seed", .GlobalEnv), before)
})

test_that("a corrected object is required", {
  set.seed(1)
  x <- sim_data(15, th, rep(100, 15))
  d <- set_meta_data(x, CC = 0.5)
  expect_error(fit_ib(suppressMessages(set_meta_data(x)), H = 10), "CC = 0.5")
  expect_error(fit_ib(x, H = 10), "accmeta_data")
  expect_error(fit_ib(d, H = 10, PRIOR = list()), "accmeta_prior")
  expect_error(fit_ib(d, H = 1), "H >= 2")
})

test_that("the result carries its path and diagnostics", {
  set.seed(1)
  d <- set_meta_data(sim_data(15, th, rep(100, 15)), CC = 0.5)
  f <- fit_ib(d, H = 10, MAX_ITER = 3, SEEDS = 1:10)
  expect_named(f, c(
    "THETA", "PI_HAT", "N_ITER", "CONVERGED", "STOP", "RESIDUAL",
    "PROGRESS", "PATH", "FAIL", "DEGEN", "SEEDS"
  ))
  expect_length(f$THETA, 9)
  expect_identical(nrow(f$PATH), f$N_ITER + 1L)
  expect_length(f$FAIL, f$N_ITER)
  expect_length(f$DEGEN, f$N_ITER)
  expect_length(f$PROGRESS, f$N_ITER)
  # the recursion starts uncorrected
  expect_equal(f$PATH[1, ], f$PI_HAT)
  expect_equal(f$PI_HAT, fit_tlmm(d, PRIOR = set_prior())$THETA)
  expect_equal(f$THETA, f$PATH[nrow(f$PATH), ])
})

test_that("STOP says which rule ended it", {
  set.seed(1)
  d <- set_meta_data(sim_data(15, th, rep(100, 15)), CC = 0.5)

  out_of_budget <- fit_ib(d, H = 10, MAX_ITER = 3, TOL = 1e-12, SEEDS = 1:10)
  expect_identical(out_of_budget$STOP, "maxit")
  expect_false(out_of_budget$CONVERGED)

  met <- fit_ib(d, H = 10, MAX_ITER = 3, TOL = 1e3, SEEDS = 1:10)
  expect_identical(met$STOP, "tol")
  expect_true(met$CONVERGED)
  expect_identical(met$N_ITER, 1L)

  # flat rule needs eleven iterations
  flat <- fit_ib(d, H = 10, MAX_ITER = 25, TOL = 1e-12, SEEDS = 1:10)
  expect_true(flat$STOP %in% c("plateau", "maxit"))
  if (identical(flat$STOP, "plateau")) expect_gt(flat$N_ITER, 10L)
})

test_that("a failed fit is redrawn, not dropped", {
  set.seed(1)
  d <- set_meta_data(sim_data(15, th, rep(100, 15)), CC = 0.5)
  tries <- 0
  real <- fit_tlmm
  local_mocked_bindings(fit_tlmm = function(...) {
    tries <<- tries + 1
    if (tries %% 3 == 0) stop("no fit")
    real(...)
  })
  f <- fit_ib(d, H = 10, MAX_ITER = 1, SEEDS = 1:10)
  expect_identical(sum(f$FAIL), 0L)
  expect_true(all(is.finite(f$THETA)))
  # more fits than replicates
  expect_gt(tries, 11)
})

test_that("the correction moves the estimate", {
  set.seed(7)
  d <- set_meta_data(sim_data(20, th, rep(100, 20)), CC = 0.5)
  f <- fit_ib(d, H = 30, MAX_ITER = 4, SEEDS = 1:30)
  expect_true(all(is.finite(f$THETA)))
  expect_false(isTRUE(all.equal(f$THETA, f$PI_HAT)))
  expect_identical(sum(f$FAIL), 0L)
  # penalty keeps fits interior
  expect_lt(max(f$DEGEN), 0.5)
})
