th <- c(2.94, -2.2, -0.4, 0.0953, 0.4, -0.5108, 0.3, 0.2, -0.6931)

test_that("the same seeds give the same answer", {
  set.seed(1)
  d <- set_meta_data(sim_data(15, th, rep(100, 15)), CC = 0.5)
  a <- fit_ib(d, CONTROL = set_ib_control(H = 10, MAX_ITER = 2, TERMINATION = "max_iter"), SEEDS = 1:10)
  b <- fit_ib(d, CONTROL = set_ib_control(H = 10, MAX_ITER = 2, TERMINATION = "max_iter"), SEEDS = 1:10)
  expect_identical(a$THETA, b$THETA)
  expect_identical(a$PATH, b$PATH)
  other <- fit_ib(d, CONTROL = set_ib_control(H = 10, MAX_ITER = 2, TERMINATION = "max_iter"), SEEDS = 101:110)
  expect_false(isTRUE(all.equal(a$THETA, other$THETA)))
})

test_that("a corrected object is required", {
  set.seed(1)
  x <- sim_data(15, th, rep(100, 15))
  d <- set_meta_data(x, CC = 0.5)
  expect_error(
    fit_ib(suppressMessages(set_meta_data(x))),
    "continuity correction"
  )
  expect_error(fit_ib(x), "accmeta_data")
  expect_error(fit_ib(d, PRIOR = list()), "accmeta_prior")
  expect_error(set_ib_control(H = 9), "H > 9")
})

test_that("the result carries its path and diagnostics", {
  set.seed(1)
  d <- set_meta_data(sim_data(15, th, rep(100, 15)), CC = 0.5)
  f <- fit_ib(d, CONTROL = set_ib_control(H = 10, MAX_ITER = 3, TERMINATION = "max_iter"), SEEDS = 1:10)
  expect_named(f, c(
    "THETA", "PI_HAT", "N_ITER", "CONVERGED", "STOP", "RESIDUAL",
    "PROGRESS", "THRESHOLD", "PATH", "FAIL", "DEGEN", "HALVED", "STEPSIZE",
    "GAP", "SE", "H", "SEEDS", "PRIOR", "CONTROL"
  ))
  expect_length(f$THETA, 9)
  expect_identical(nrow(f$PATH), f$N_ITER + 1L)
  expect_length(f$FAIL, f$N_ITER)
  expect_length(f$DEGEN, f$N_ITER)
  expect_length(f$PROGRESS, f$N_ITER)
  expect_length(f$HALVED, f$N_ITER)
  expect_length(f$STEPSIZE, f$N_ITER)
  expect_length(f$THRESHOLD, f$N_ITER)
  expect_identical(dim(f$GAP), c(f$N_ITER, 9L))
  expect_identical(dim(f$SE), c(f$N_ITER, 9L))
  expect_equal(f$PATH[1, ], f$PI_HAT)
  expect_identical(sum(f$HALVED), 0L)
  expect_equal(f$PI_HAT, fit_tlmm(d, PRIOR = set_prior())$THETA)
  expect_equal(f$THETA, f$PATH[which.min(f$PROGRESS), ])
})

test_that("NCORES > 1 refuses to run over a pre-existing pool", {
  skip_if_not_installed("mirai")
  set.seed(1)
  d <- set_meta_data(sim_data(15, th, rep(100, 15)), CC = 0.5)
  mirai::daemons(1)
  on.exit(mirai::daemons(0), add = TRUE)
  expect_error(
    fit_ib(d, CONTROL = set_ib_control(H = 20, NCORES = 2, TERMINATION = "max_iter")),
    "Tear down"
  )
})

test_that("NCORES > 1 over mirai matches the serial fit", {
  skip_if_not_installed("mirai")
  skip_if_not(
    "accmeta" %in% rownames(utils::installed.packages()),
    "accmeta must be installed for mirai daemons"
  )
  set.seed(1)
  d <- set_meta_data(sim_data(15, th, rep(100, 15)), CC = 0.5)
  ctrl <- function(ncores) {
    set_ib_control(H = 20, MAX_ITER = 3, NCORES = ncores, TERMINATION = "max_iter")
  }
  ser <- fit_ib(d, WORKPAR = "Joe", SEEDS = 1:20, CONTROL = ctrl(1))
  par <- fit_ib(d, WORKPAR = "Joe", SEEDS = 1:20, CONTROL = ctrl(2))
  expect_equal(par$THETA, ser$THETA)
  expect_equal(par$PATH, ser$PATH)
  boosted <- function(ncores) {
    set.seed(7)
    fit_ib(
      d, SEEDS = 1:30,
      CONTROL = set_ib_control(
        H = 30, MAX_ITER = 6, STEP = 1, TOL = 0.01, PRECISION = 1,
        BOOST = TRUE, MAX_H = 120, PATIENCE = 6, NCORES = ncores
      )
    )
  }
  bs <- boosted(1)
  bp <- boosted(2)
  expect_gt(max(bs$H), 30)
  expect_identical(bp$H, bs$H)
  expect_identical(bp$SEEDS, bs$SEEDS)
  expect_equal(bp$PATH, bs$PATH)
  expect_identical(mirai::status()$connections, 0L)
})

test_that("STOP says which rule ended it", {
  set.seed(1)
  d <- set_meta_data(sim_data(15, th, rep(100, 15)), CC = 0.5)

  out_of_budget <- fit_ib(
    d, CONTROL = set_ib_control(
      H = 30, MAX_ITER = 3, TOL = 0.999, TERMINATION = "hotelling"
    ), SEEDS = 1:30
  )
  expect_identical(out_of_budget$STOP, "maxit")
  expect_false(out_of_budget$CONVERGED)

  met <- fit_ib(
    d, CONTROL = set_ib_control(
      H = 30, MAX_ITER = 3, TOL = 1e-20, TERMINATION = "hotelling"
    ), SEEDS = 1:30
  )
  expect_identical(met$STOP, "tol")
  expect_true(met$CONVERGED)
  expect_identical(met$N_ITER, 1L)

  flat <- fit_ib(
    d, CONTROL = set_ib_control(
      H = 30, MAX_ITER = 25, TOL = 0.999, STEP = 0.1,
      TERMINATION = "hotelling"
    ), SEEDS = 1:30
  )
  expect_true(flat$STOP %in% c("stall", "maxit"))
  expect_false(flat$CONVERGED)
})

test_that("confidence termination is ergm's ellipsoid inclusion test", {
  set.seed(1)
  d <- set_meta_data(sim_data(15, th, rep(100, 15)), CC = 0.5)
  run <- function(prec) {
    fit_ib(
      d,
      CONTROL = set_ib_control(
        H = 30, MAX_ITER = 1, TOL = 0.01, TERMINATION = "confidence",
        PRECISION = prec
      ),
      SEEDS = 1:30
    )
  }
  f <- run(10)
  h_ok <- 30 - f$FAIL[1]
  crit <- 9 * (h_ok - 1) / (h_ok - 9) * qf(0.99, 9, h_ok - 9)
  expect_equal(
    f$THRESHOLD[1],
    (h_ok - 1) * (sqrt(10) - sqrt(crit / (h_ok - 1)))^2
  )
  r <- sqrt(f$PROGRESS[1] / (h_ok - 1))
  expect_identical(
    f$STOP == "tol",
    r < sqrt(10) && (h_ok - 1) * (sqrt(10) - r)^2 > crit
  )
  tight <- suppressWarnings(run(1))
  expect_identical(tight$THRESHOLD[1], 0)
  expect_false(tight$CONVERGED)
})

test_that("BOOST grows H up to MAX_H and extends the seeds", {
  set.seed(1)
  d <- set_meta_data(sim_data(15, th, rep(100, 15)), CC = 0.5)
  f <- fit_ib(
    d,
    CONTROL = set_ib_control(
      H = 30, MAX_ITER = 6, STEP = 1, TOL = 0.01, TERMINATION = "confidence",
      PRECISION = 1, BOOST = TRUE, MAX_H = 120, PATIENCE = 6
    ),
    SEEDS = 1:30
  )
  expect_length(f$H, f$N_ITER)
  expect_identical(f$H[1], 30)
  expect_false(is.unsorted(f$H))
  expect_gt(max(f$H), 30)
  expect_lte(max(f$H), 120)
  expect_identical(f$SEEDS[1:30], 1:30)
  expect_length(f$SEEDS, tail(f$H, 1))
  expect_true(f$CONVERGED)
  expect_equal(f$THETA, f$PATH[f$N_ITER, ])
  slow <- fit_ib(
    d,
    CONTROL = set_ib_control(
      H = 30, MAX_ITER = 6, STEP = 1, TOL = 0.01, TERMINATION = "confidence",
      PRECISION = 1, BOOST = TRUE, BOOST_FACTOR = 1.1, MAX_H = 120,
      PATIENCE = 6
    ),
    SEEDS = 1:30
  )
  expect_gt(max(slow$H), 30)
  expect_true(all(tail(slow$H, -1) <= ceiling(1.1 * head(slow$H, -1))))
  wide <- fit_ib(
    d,
    CONTROL = set_ib_control(
      H = 30, MAX_ITER = 6, STEP = 1, TOL = 0.01, PRECISION = 1,
      BOOST = TRUE, BOOST_FACTOR = 10, MAX_H = 1000, PATIENCE = 6
    ),
    SEEDS = 1:30
  )
  k <- which(diff(wide$H) > 0)
  expect_gt(length(k), 0)
  h_ok <- wide$H[k] - wide$FAIL[k]
  r2 <- wide$PROGRESS[k] / (h_ok - 1)
  need <- ib_crit(h_ok, 0.01) / ((h_ok - 1) * (1 - sqrt(r2))^2)
  expect_true(all(r2 < 1))
  expect_identical(wide$H[k + 1], ceiling(wide$H[k] * pmin(need, 10)))
  g <- suppressWarnings(fit_ib(
    d,
    CONTROL = set_ib_control(
      H = 30, MAX_ITER = 2, TERMINATION = "confidence", PRECISION = 1,
      BOOST = FALSE
    ),
    SEEDS = 1:30
  ))
  expect_identical(g$H, c(30, 30))
  expect_identical(g$SEEDS, 1:30)
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
  f <- fit_ib(d, CONTROL = set_ib_control(H = 10, MAX_ITER = 1, TERMINATION = "max_iter"), SEEDS = 1:10)
  expect_identical(sum(f$FAIL), 0L)
  expect_true(all(is.finite(f$THETA)))
  expect_gt(tries, 11)
})

test_that("the correction moves the estimate", {
  set.seed(7)
  d <- set_meta_data(sim_data(20, th, rep(100, 20)), CC = 0.5)
  f <- fit_ib(d, CONTROL = set_ib_control(H = 30, MAX_ITER = 4, TERMINATION = "max_iter"), SEEDS = 1:30)
  expect_true(all(is.finite(f$THETA)))
  expect_false(isTRUE(all.equal(f$THETA, f$PI_HAT)))
  expect_identical(sum(f$FAIL), 0L)
  expect_lt(max(f$DEGEN), 0.5)
})

test_that("a boundary start is projected inward", {
  sd_true <- sqrt(c(1.2, 0.5, 0.25))
  cor_true <- matrix(c(1, -0.6, 0.7, -0.6, 1, -0.7, 0.7, -0.7, 1), 3, 3)
  tv <- list2theta(list(
    MU = c(2.94, -2.20, -0.405),
    SIGMA = diag(sd_true) %*% cor_true %*% diag(sd_true)
  ))
  set.seed(123)
  ss <- sample(40:200, 15, TRUE)
  set.seed(15)
  d <- set_meta_data(sim_data(15, tv, ss), CC = 0.5)

  flat <- fit_tlmm(d, PRIOR = set_prior(4))$THETA
  ev <- eigen(theta2list(flat)$SIGMA, symmetric = TRUE, only.values = TRUE)$values
  expect_lt(min(ev), 1e-4)

  f <- suppressWarnings(
    fit_ib(
      d, CONTROL = set_ib_control(H = 10, MAX_ITER = 3, TERMINATION = "max_iter"),
      PRIOR = set_prior(4), SEEDS = 1:10
    )
  )
  expect_true(all(is.finite(f$THETA)))
  expect_equal(f$PI_HAT, flat)
  expect_false(isTRUE(all.equal(f$PATH[1, ], f$PI_HAT)))
  start_ev <- eigen(
    theta2list(f$PATH[1, ])$SIGMA,
    symmetric = TRUE, only.values = TRUE
  )$values
  expect_gte(min(start_ev), 1e-4 - 1e-9)
})

test_that("a runaway update is halved, and a boundary root stops early", {
  sd_true <- sqrt(c(1.2, 0.5, 0.25))
  cor_true <- matrix(c(1, -0.6, 0.7, -0.6, 1, -0.7, 0.7, -0.7, 1), 3, 3)
  tv <- list2theta(list(
    MU = c(2.94, -2.20, -0.405),
    SIGMA = diag(sd_true) %*% cor_true %*% diag(sd_true)
  ))
  set.seed(123)
  ss <- sample(40:200, 15, TRUE)
  set.seed(1)
  d <- set_meta_data(sim_data(15, tv, ss), CC = 0.5)
  f <- suppressWarnings(
    fit_ib(
      d, CONTROL = set_ib_control(H = 20, MAX_ITER = 6, STEP = 1),
      PRIOR = set_prior(4), SEEDS = 1:20
    )
  )
  expect_true(all(is.finite(f$THETA)))
  expect_gt(sum(f$HALVED), 0)
  expect_identical(f$STOP, "boundary")
  expect_true(all(tail(f$HALVED, 2) > 0))
  expect_false(f$CONVERGED)
  expect_equal(f$THETA, f$PATH[which.min(f$PROGRESS), ])
})

test_that("max_iter runs the whole budget past the early stops", {
  sd_true <- sqrt(c(1.2, 0.5, 0.25))
  cor_true <- matrix(c(1, -0.6, 0.7, -0.6, 1, -0.7, 0.7, -0.7, 1), 3, 3)
  tv <- list2theta(list(
    MU = c(2.94, -2.20, -0.405),
    SIGMA = diag(sd_true) %*% cor_true %*% diag(sd_true)
  ))
  set.seed(123)
  ss <- sample(40:200, 15, TRUE)
  set.seed(1)
  d <- set_meta_data(sim_data(15, tv, ss), CC = 0.5)
  f <- suppressWarnings(
    fit_ib(
      d,
      CONTROL = set_ib_control(
        H = 20, MAX_ITER = 6, STEP = 1, TERMINATION = "max_iter"
      ),
      PRIOR = set_prior(4), SEEDS = 1:20
    )
  )
  expect_identical(f$N_ITER, 6L)
  expect_identical(f$STOP, "maxit")
  expect_false(f$CONVERGED)
  expect_true(all(is.na(f$THRESHOLD)))
  expect_equal(f$THETA, f$PATH[which.min(f$PROGRESS), ])
})
