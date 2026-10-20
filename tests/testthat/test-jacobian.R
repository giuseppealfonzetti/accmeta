th <- c(2.94, -2.2, -0.4, 0.0953, 0.4, -0.5108, 0.3, 0.2, -0.6931)

ib_small <- function(VCOV = FALSE) {
  set.seed(1)
  d <- set_meta_data(sim_data(10, th, rep(100, 10)), CC = 0.5)
  fit_ib(
    d,
    CONTROL = set_ib_control(
      H = 20, MAX_ITER = 2, TERMINATION = "max_iter",
      VCOV_H = 30, VCOV_N_NODES = 5L
    ),
    SEEDS = 1:20,
    VCOV = VCOV
  )
}

test_that("ib_jacobian works in the IB working coordinates", {
  joe <- ib_small()
  pb <- joe
  pb$WORKPAR <- "PinheiroBates"
  J_joe <- ib_jacobian(joe, H = 30, N_NODES = 5L, SEEDS = 1:30)
  J_pb <- ib_jacobian(pb, H = 30, N_NODES = 5L, SEEDS = 1:30)
  expect_identical(J_joe$H_OK, 30L)
  expect_identical(dim(J_joe$J), c(9L, 9L))
  expect_true(all(is.finite(J_joe$J)) && all(is.finite(J_joe$V)))
  expect_equal(t(apply(J_pb$PI, 1, theta2joe)), J_joe$PI, tolerance = 1e-4)
})

test_that("ib_jacobian over mirai matches the serial run", {
  skip_if_not_installed("mirai")
  skip_if_not(
    "accmeta" %in% rownames(utils::installed.packages()),
    "accmeta must be installed for mirai daemons"
  )
  fit <- ib_small()
  ser <- ib_jacobian(fit, H = 30, N_NODES = 5L, NCORES = 1L, SEEDS = 1:30)
  par <- ib_jacobian(fit, H = 30, N_NODES = 5L, NCORES = 2L, SEEDS = 1:30)
  expect_equal(par$J, ser$J)
  expect_identical(mirai::status()$connections, 0L)
})

test_that("the Jacobian is added by fit_ib or update, and dropped by update", {
  plain <- ib_small()
  expect_null(plain$JACOBIAN)
  expect_identical(ib_small(VCOV = TRUE)$JACOBIAN$H_OK, 30L)
  fit <- update(plain, VCOV = TRUE)
  expect_identical(fit$JACOBIAN$H_OK, 30L)
  expect_identical(fit$THETA, plain$THETA)
  expect_null(update(fit, VCOV = FALSE)$JACOBIAN)
  expect_error(update(fit, H = 3), "only takes VCOV")
})

test_that("vcov is the sandwich on the working scale", {
  expect_error(vcov(ib_small()), "No Jacobian")
  fit <- ib_small(VCOV = TRUE)
  jac <- fit$JACOBIAN
  j_inv <- solve(jac$J)
  v <- vcov(fit)
  expect_equal(
    unname(v),
    (1 + 1 / fit$CONTROL$H) * j_inv %*% jac$V %*% t(j_inv)
  )
  expect_true(isSymmetric(v))
  expect_identical(rownames(v)[c(1, 9)], c("mu_eta", "z_rho23_1"))
})

test_that("confint gives Wald and basic intervals", {
  fit <- ib_small(VCOV = TRUE)
  wald <- confint(fit)
  basic <- confint(fit, METHOD = "basic")
  se <- sqrt(diag(vcov(fit)))
  expect_equal(
    unname(wald[, 1]),
    theta2joe(fit$THETA) - stats::qnorm(0.975) * unname(se)
  )
  for (ci in list(wald, basic)) {
    expect_identical(dim(ci), c(9L, 2L))
    expect_identical(colnames(ci), c("2.5 %", "97.5 %"))
    expect_true(all(ci[, 1] < ci[, 2]))
  }
  expect_identical(rownames(confint(fit, "mu_eta")), "mu_eta")
})

test_that("ib_jacobian leaves the caller's RNG state untouched", {
  fit <- ib_small()
  withr::local_seed(1, .rng_kind = "L'Ecuyer-CMRG")
  before <- .Random.seed
  ib_jacobian(fit, H = 30, N_NODES = 5L, SEEDS = 1:30)
  expect_identical(.Random.seed, before)
})
