th <- c(2.94, -2.2, -0.4, 0.0953, 0.4, -0.5108, 0.3, 0.2, -0.6931)

# mode sits far from start
far <- matrix(
  c(54L, 12L, 48L, 286L), 1, 4,
  dimnames = list(NULL, c("n11", "n10", "n01", "n00"))
)

test_that("the mode is found for a far study", {
  d <- set_meta_data(far)
  obj <- fit_tglmm(
    d, THETA_START = th, N_NODES = 15L, PRIOR = set_prior(4),
    CONTROL = list(maxeval = 1)
  )$OBJ
  # value from importance sampling
  expect_equal(obj$fn(th), 355.1859, tolerance = 1e-5)
  expect_lt(max(abs(as.numeric(obj$gr(th)))), 100)
})

test_that("INNER_GRAD_MAX shows inner convergence", {
  d <- set_meta_data(far)
  slack <- fit_tglmm(
    d, THETA_START = th, N_NODES = 15L, N_ITER = 10L, CONTROL = list(maxeval = 1)
  )$OBJ
  tight <- fit_tglmm(
    d, THETA_START = th, N_NODES = 15L, N_ITER = 1L, CONTROL = list(maxeval = 1)
  )$OBJ
  expect_lt(slack$report(th)$INNER_GRAD_MAX, 1e-6)
  # too few steps show up
  expect_gt(tight$report(th)$INNER_GRAD_MAX, 1e-4)
})

test_that("the quadrature has enough nodes", {
  set.seed(1)
  d <- suppressMessages(set_meta_data(sim_data(5, th, rep(100, 5))))
  ten <- fit_tglmm(
    d, THETA_START = th, N_NODES = 10L, CONTROL = list(maxeval = 1)
  )$OBJ$fn(th)
  twenty <- fit_tglmm(
    d, THETA_START = th, N_NODES = 20L, CONTROL = list(maxeval = 1)
  )$OBJ$fn(th)
  expect_equal(ten, twenty, tolerance = 1e-6)
})

test_that("the prior adds exactly the Wishart terms", {
  set.seed(1)
  d <- suppressMessages(set_meta_data(sim_data(5, th, rep(100, 5))))
  S <- theta2list(th)$SIGMA
  # the baseline must be flat
  flat <- fit_tglmm(
    d, THETA_START = th, N_NODES = 15L, PRIOR = set_prior(4, Inf),
    CONTROL = list(maxeval = 1)
  )$OBJ$fn(th)
  scaled <- fit_tglmm(
    d, THETA_START = th, N_NODES = 15L, PRIOR = set_prior(4, 1),
    CONTROL = list(maxeval = 1)
  )$OBJ$fn(th)
  barrier <- fit_tglmm(
    d, THETA_START = th, N_NODES = 15L, PRIOR = set_prior(9, Inf),
    CONTROL = list(maxeval = 1)
  )$OBJ$fn(th)
  # scale leaves the trace term
  expect_equal(scaled - flat, sum(diag(S)) / 2)
  # degrees leaves the barrier
  expect_equal(barrier - flat, -(9 - 4) / 2 * log(det(S)))
})

test_that("fit_tlmm reproduces its recorded fit", {
  set.seed(1)
  d <- set_meta_data(sim_data(50, th, rep(100, 50)), CC = 0.5)
  f <- fit_tlmm(d, PRIOR = set_prior(4))
  expect_equal(f$NLL, 157.105, tolerance = 1e-5)
  expect_equal(f$THETA[1:3], c(2.48487, -1.91188, -0.37203), tolerance = 1e-4)
  expect_lt(max(abs(as.numeric(f$OBJ$gr(f$THETA)))), 1e-4)
  expect_equal(fit_tlmm(d, PRIOR = set_prior(10, 50))$NLL, 169.3108, tolerance = 1e-5)
})

test_that("fit_tglmm recovers what generated the data", {
  set.seed(3)
  d <- suppressMessages(set_meta_data(sim_data(150, th, rep(400, 150))))
  f <- fit_tglmm(d, N_NODES = 8L)
  expect_lt(max(abs(as.numeric(f$OBJ$gr(f$THETA)))), 1e-3)
  expect_lt(f$OBJ$report(f$THETA)$INNER_GRAD_MAX, 1e-5)
  expect_equal(f$THETA[1:3], th[1:3], tolerance = 0.1)
  # exact must beat approximate
  a <- fit_tlmm(suppressMessages(set_meta_data(d$tab, CC = 0.5)))
  expect_lt(sum(abs(f$THETA[1:3] - th[1:3])), sum(abs(a$THETA[1:3] - th[1:3])))
})

test_that("CONTROL reaches the optimiser", {
  set.seed(1)
  d <- set_meta_data(sim_data(20, th, rep(100, 20)), CC = 0.5)
  # three means budget exhausted
  expect_identical(fit_tlmm(d, CONTROL = list(maxeval = 1))$CONVERGENCE, 3L)
  expect_false(fit_tlmm(d)$CONVERGENCE == 3L)
})
