th5 <- c(2.94, -2.2, 0.0953, 0.4, -0.5108)

biv_data <- function(CC = 0.5, N = 20) {
  set.seed(4)
  ni <- sample(40:200, N, TRUE)
  n1 <- round(ni * stats::runif(N, 0.2, 0.6))
  set_meta_data(sim_data(N, th5, ni, N1 = n1), CC = CC, VERBOSE = FALSE)
}

test_that("length-5 vectors map like the leading block of length 9", {
  li <- theta2list(th5)
  expect_named(li$MU, c("eta", "xi"))
  expect_identical(dim(li$SIGMA), c(2L, 2L))
  expect_equal(list2theta(li), th5)
  expect_equal(joe2theta(theta2joe(th5)), th5)
  expect_equal(
    theta2joe(th5),
    theta2joe(embed(th5, "PinheiroBates"))[c(1, 2, 4, 5, 7)]
  )
  expect_error(sim_data(3, th5, rep(100, 3)), "N1")
})

test_that("the default prior has p + 2 degrees", {
  d <- biv_data()
  expect_identical(
    fit_mle(d, "bivariate", "lmm", NULL, set_prior(), list(), "PinheiroBates")$PRIOR$DEGREES,
    4
  )
  expect_identical(fit_tlmm(set_meta_data(d$tab, CC = 0.5))$PRIOR$DEGREES, 5)
  expect_error(fit_tlmm(set_meta_data(d$tab, CC = 0.5), PRIOR = set_prior(3)), "at least 4")
})

test_that("BLMM sums the bivariate normal log-densities", {
  d <- biv_data()
  o <- make_obj(d, th5, set_prior(3, Inf), "PinheiroBates", "tlmm")
  li <- theta2list(th5)
  ll <- sum(sapply(seq_len(d$n_studies), function(i) {
    s <- li$SIGMA + diag(d$wvar[i, 1:2])
    r <- d$est[i, 1:2] - li$MU
    -0.5 * (2 * log(2 * pi) + determinant(s)$modulus + sum(r * solve(s, r)))
  }))
  expect_equal(o$fn(o$par), -ll)
  b <- fit_mle(d, "bivariate", "lmm", NULL, set_prior(3, Inf), list(), "PinheiroBates")
  expect_s3_class(b, c("accmeta_blmm", "accmeta_fit"), exact = TRUE)
  expect_lt(max(abs(b$OBJ$gr(b$OBJ$env$last.par.best))), 1e-3)
})

test_that("BGLMM integrates one study like integrate()", {
  d <- set_meta_data(biv_data()$tab[1, , drop = FALSE])
  o <- make_obj(d, th5, set_prior(3, Inf), "PinheiroBates", "tglmm", 20L)
  L <- t(chol(theta2list(th5)$SIGMA))
  y <- d$counts[1, ]
  f2 <- function(u1, u2) {
    e <- th5[1] + L[1, 1] * u1
    x <- th5[2] + L[2, 1] * u1 + L[2, 2] * u2
    exp(
      stats::dbinom(y[1], y[4], stats::plogis(e), log = TRUE) +
        stats::dbinom(y[2], y[5], stats::plogis(x), log = TRUE)
    ) *
      stats::dnorm(u1) *
      stats::dnorm(u2)
  }
  inner <- function(u1) {
    sapply(u1, function(a) {
      stats::integrate(function(u2) f2(a, u2), -10, 10, rel.tol = 1e-12)$value
    })
  }
  ll <- log(stats::integrate(inner, -10, 10, rel.tol = 1e-12)$value) -
    lchoose(y[4], y[1]) -
    lchoose(y[5], y[2])
  expect_equal(-o$fn(o$par), unname(ll), tolerance = 1e-6)
})

test_that("bivariate IB fits and reports five parameters", {
  d <- biv_data()
  f <- fit_ib(
    d,
    CONTROL = set_ib_control(H = 20, MAX_ITER = 3, TERMINATION = "max_iter", VCOV_H = 30, VCOV_N_NODES = 5L),
    SEEDS = 1:20,
    VCOV = TRUE,
    MODEL = "bivariate"
  )
  expect_length(f$THETA, 5)
  expect_identical(f$MODEL, "bivariate")
  expect_identical(dim(vcov(f)), c(5L, 5L))
  nat <- confint(f, SCALE = "natural")
  expect_identical(
    rownames(nat),
    c("mu_eta", "mu_xi", "sd_eta", "sd_xi", "rho_eta_xi")
  )
  rep <- confint(f, SCALE = "reporting")
  expect_identical(rownames(rep)[1:2], c("sensitivity", "specificity"))
  expect_true(all(rep[1:2, ] > 0 & rep[1:2, ] < 1 & rep[1:2, 1] < rep[1:2, 2]))
  expect_equal(unname(rep[3:5, ]), unname(nat[3:5, ]))
  expect_output(print(f), "SP")
})

test_that("fit_meta_model reaches every model and estimator", {
  th <- c(2.94, -2.2, -0.4, 0.0953, 0.4, -0.5108, 0.3, 0.2, -0.6931)
  set.seed(1)
  design <- rep(c("cohort", "case-control"), 8)
  x <- sim_data(16, th, rep(100, 16), N1 = ifelse(design == "cohort", NA, 40))
  d <- set_meta_data(x, CC = 0.5, DESIGN = design)
  d0 <- set_meta_data(x, VERBOSE = FALSE, DESIGN = design)
  for (m in c("trivariate", "bivariate")) {
    p <- if (m == "trivariate") 9L else 5L
    l <- fit_meta_model(d, m, "lmm")
    g <- fit_meta_model(d0, m, "glmm", CONTROL = set_control("glmm", N_NODES = 5L))
    b <- fit_meta_model(
      d,
      m,
      "ib",
      CONTROL = set_control("ib", H = 10, MAX_ITER = 2, TERMINATION = "max_iter", SEEDS = 1:10)
    )
    for (f in list(l, g, b)) {
      expect_length(f$THETA, p)
      expect_identical(f$MODEL, m)
    }
    expect_identical(c(l$ESTIMATOR, g$ESTIMATOR, b$ESTIMATOR), c("lmm", "glmm", "ib"))
    expect_identical(dim(vcov(l)), c(p, p))
  }
  expect_identical(
    fit_meta_model(
      d0,
      "trivariate",
      "glmm",
      CONTROL = set_control("glmm", WORKPAR = "PinheiroBates", N_NODES = 5L)
    )$THETA,
    fit_tglmm(d0, N_NODES = 5L)$THETA
  )
  expect_error(set_control("lmm", N_NODES = 5L), "N_NODES not used")
  expect_error(set_control("glmm", H = 10, SEEDS = 1:10), "SEEDS, H not used")
  expect_error(set_control("ib", H = 5), "H > 9")
  expect_error(fit_meta_model(d, "trivariate", "glmm", CONTROL = set_control("lmm")), "set_control\\('glmm'\\)")
  expect_identical(set_control("lmm")$OPTIM, list(maxeval = 1000))
  expect_identical(set_control("ib")$IB, set_ib_control())
  cc_only <- set_meta_data(x, CC = 0.5, DESIGN = rep("case-control", 16))
  expect_error(fit_meta_model(cc_only, "trivariate", "lmm"), "one cohort")
  expect_length(fit_meta_model(cc_only, "bivariate", "lmm")$THETA, 5)
})
