th <- c(2.94, -2.2, -0.4, 0.0953, 0.4, -0.5108, 0.3, 0.2, -0.6931)

mixed <- function(CC = 0.5, N = 12) {
  set.seed(1)
  design <- rep(c("cohort", "case-control"), length.out = N)
  ni <- sample(40:200, N, TRUE)
  n1 <- ifelse(design == "case-control", round(ni / 2), NA)
  set_meta_data(
    sim_data(N, th, ni, N1 = n1),
    CC = CC,
    VERBOSE = FALSE,
    DESIGN = design
  )
}

test_that("case-control studies carry no prevalence", {
  d <- mixed()
  cc <- d$design == "case-control"
  expect_true(all(d$counts[cc, c("y_gamma", "den_gamma")] == 0))
  expect_true(all(is.na(d$est[cc, "gamma"])))
  expect_true(all(is.na(d$rates[cc, "PREV"])))
  expect_true(all(is.finite(d$est[!cc, ])))
  expect_identical(set_meta_data(d$tab, CC = 0.5)$design, rep("cohort", 12))
  expect_error(set_meta_data(d$tab, DESIGN = "cohort"))
  expect_error(set_meta_data(d$tab, DESIGN = rep("other", 12)))
})

test_that("sim_data keeps the diseased of case-control studies", {
  ni <- rep(100, 6)
  n1 <- c(NA, 30, NA, 50, NA, 70)
  set.seed(3)
  x <- sim_data(6, th, ni, N1 = n1)
  expect_identical(unname((x[, "n11"] + x[, "n01"])[c(2, 4, 6)]), c(30, 50, 70))
  set.seed(3)
  a <- sim_data(6, th, ni)
  set.seed(3)
  expect_identical(sim_data(6, th, ni, N1 = rep(NA, 6)), a)
})

test_that("TLMM uses the bivariate marginal for case-control studies", {
  d <- mixed()
  obj <- make_obj(d, th, set_prior(4, Inf), "PinheiroBates", "tlmm")
  li <- theta2list(th)
  ll <- sum(sapply(seq_len(d$n_studies), function(i) {
    k <- if (d$design[i] == "cohort") 1:3 else 1:2
    s <- li$SIGMA[k, k] + diag(d$wvar[i, k], length(k))
    r <- d$est[i, k] - li$MU[k]
    -0.5 * (length(k) * log(2 * pi) + determinant(s)$modulus + sum(r * solve(s, r)))
  }))
  expect_equal(obj$fn(th), -ll)
})

test_that("TGLMM integrates out the prevalence of case-control studies exactly", {
  d <- mixed(CC = 0, N = 4)
  cc <- make_obj(d, th, set_prior(4, Inf), "Joe", "tglmm", 7L)
  gh <- statmod::gauss.quad(7, "hermite")
  full <- TMB::MakeADFun(
    data = list(
      MODEL = "tglmm",
      y = as.integer(t(d$counts[, 1:3])),
      den = as.integer(t(d$counts[, 4:6])),
      f = as.factor(rep(1:4, each = 3)),
      niter = 10L,
      ws = exp(gh$nodes^2) * gh$weights,
      z = gh$nodes,
      DEGREES = 4,
      SCALE = Inf,
      WORKPAR = "Joe",
      COHORT = rep(1L, 4)
    ),
    parameters = list(MU = th[1:3], ALPHA = theta2joe(th)[4:9]),
    DLL = "accmeta",
    silent = TRUE
  )
  w <- theta2joe(th)
  expect_equal(cc$fn(w), full$fn(w), tolerance = 1e-10)
  expect_equal(cc$gr(w), full$gr(w), tolerance = 1e-8)
})

test_that("all estimators fit mixed designs", {
  d <- mixed()
  expect_true(all(is.finite(init_theta(d))))
  expect_true(all(is.finite(fit_tlmm(d)$THETA)))
  g <- fit_tglmm(mixed(CC = 0), N_NODES = 5L)
  expect_true(all(is.finite(g$THETA)))
  expect_true(all(is.finite(confint(g))))
  b <- fit_ib(
    d,
    CONTROL = set_ib_control(H = 20, MAX_ITER = 2, TERMINATION = "max_iter"),
    SEEDS = 1:20
  )
  expect_true(all(is.finite(b$THETA)))
  s <- withr::with_seed(1, sim_like(d, th))
  expect_identical(s$design, d$design)
  cc <- d$design == "case-control"
  expect_identical(s$margins[cc, "n1"], d$margins[cc, "n1"])
})
