th <- c(2.94, -2.2, -0.4, 0.0953, 0.4, -0.5108, 0.3, 0.2, -0.6931)
s12 <- c(-0.8, -0.3, -0.3)
th12 <- list2theta(list(
  MU = theta2list(th)$MU,
  SIGMA = theta2list(th)$SIGMA,
  SKEW = s12
))

test_that("theta2list and list2theta round-trip on the natural scale", {
  li <- theta2list(th12)
  expect_named(li, c("MU", "SIGMA", "SKEW"))
  expect_equal(unname(li$SKEW), s12)
  expect_equal(li$SIGMA, theta2list(th)$SIGMA)
  expect_equal(list2theta(li), th12)
  expect_equal(joe2theta(theta2joe(th12)), th12)
  expect_error(
    list2theta(modifyList(li, list(SKEW = c(-0.8, -0.5, 0)))),
    "not admissible"
  )
})

test_that("theta2dp gives the direct parameters of sn", {
  li <- theta2list(th12)
  dp <- theta2dp(th12)
  ref <- sn::cp2dp(list(mean = li$MU, var.cov = li$SIGMA, gamma1 = s12), "SN")
  expect_equal(unname(dp$xi), unname(ref[[1]]))
  expect_equal(unname(dp$Omega), unname(ref$Omega))
  expect_equal(unname(dp$alpha), unname(ref$alpha))
  cp <- sn::dp2cp(dp, "SN")
  expect_equal(unname(cp$mean), unname(li$MU))
  expect_equal(unname(cp$var.cov), unname(li$SIGMA))
  expect_equal(unname(cp$gamma1), s12)
  expect_equal(theta2dp(th)$alpha, c(eta = 0, xi = 0, gamma = 0))
})

test_that("THETA is the original parameterisation of sn", {
  op <- sn::dp2op(theta2dp(th12), "SN")
  expect_equal(unname(op$xi), th12[1:3])
  expect_equal(unname(op$Psi), unname(theta2list(th12[1:9])$SIGMA))
  expect_equal(unname(op$lambda), th12[10:12])
})


test_that("sim_data draws skew-normal random effects", {
  set.seed(1)
  e <- suppressMessages(
    set_meta_data(sim_data(5e4, th12, rep(1e8, 5e4)), CC = 0.5)
  )$est
  li <- theta2list(th12)
  z <- sweep(e, 2, colMeans(e))
  expect_equal(unname(colMeans(e)), th[1:3], tolerance = 1e-2)
  expect_equal(unname(cov(e)), unname(li$SIGMA), tolerance = 2e-2)
  expect_equal(unname(colMeans(z^3) / colMeans(z^2)^1.5), s12, tolerance = 0.1)
})
