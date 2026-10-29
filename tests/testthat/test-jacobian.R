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

test_that("confint back-transforms to the natural scale", {
  fit <- ib_small(VCOV = TRUE)
  rho <- cov2cor(theta2list(fit$THETA)$SIGMA)[2, 3]
  for (m in c("wald", "basic")) {
    w <- confint(fit, METHOD = m)
    n <- confint(fit, METHOD = m, SCALE = "natural")
    expect_identical(rownames(n)[c(4, 9)], c("sd_eta", "rho_xi_gamma"))
    expect_equal(unname(n[1:3, ]), unname(w[1:3, ]))
    expect_equal(unname(n[4:6, ]), unname(exp(w[4:6, ])))
    expect_equal(unname(n[7:8, ]), unname(tanh(w[7:8, ])))
    expect_true(n[9, 1] < rho && rho < n[9, 2])
    expect_true(all(abs(n[7:9, ]) < 1))
    p <- confint(fit, METHOD = m, SCALE = "reporting")
    expect_identical(rownames(p)[1:3], c("sensitivity", "specificity", "prevalence"))
    expect_equal(unname(p[4:9, ]), unname(n[4:9, ]))
    expect_equal(unname(p[c(1, 3), ]), unname(plogis(n[c(1, 3), ])))
    expect_equal(unname(p[2, ]), unname(1 - plogis(n[2, 2:1])))
    expect_true(all(p[1:3, 1] < p[1:3, 2] & p[1:3, ] > 0 & p[1:3, ] < 1))
  }
  zr <- function(x) {
    r <- tanh(x[7:9])
    atanh(r[1] * r[2] + r[3] * sqrt((1 - r[1]^2) * (1 - r[2]^2)))
  }
  joe <- theta2joe(fit$THETA)
  grad <- sapply(1:9, function(k) {
    e <- replace(numeric(9), k, 1e-6)
    (zr(joe + e) - zr(joe - e)) / 2e-6
  })
  se <- sqrt(drop(crossprod(grad, vcov(fit) %*% grad)))
  expect_equal(
    unname(confint(fit, SCALE = "natural")[9, ]),
    tanh(atanh(rho) + stats::qnorm(c(0.025, 0.975)) * se),
    tolerance = 1e-6
  )
  pb <- fit
  pb$WORKPAR <- "PinheiroBates"
  expect_error(confint(pb, SCALE = "natural"), "Joe")
})

test_that("ib_jacobian leaves the caller's RNG state untouched", {
  fit <- ib_small()
  withr::local_seed(1, .rng_kind = "L'Ecuyer-CMRG")
  before <- .Random.seed
  ib_jacobian(fit, H = 30, N_NODES = 5L, SEEDS = 1:30)
  expect_identical(.Random.seed, before)
})
