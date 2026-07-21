th <- c(2.94, -2.2, -0.4, 0.0953, 0.4, -0.5108, 0.3, 0.2, -0.6931)

test_that("summaries match their definitions", {
  set.seed(1)
  x <- sim_data(20, th, rep(100, 20))
  d <- set_meta_data(x, CC = 0.5)
  n1 <- x[, "n11"] + x[, "n01"]
  n0 <- x[, "n10"] + x[, "n00"]

  expect_s3_class(d, "accmeta_data")
  expect_identical(d$n_studies, 20L)
  expect_identical(d$CC, 0.5)
  expect_equal(d$tab, x)
  expect_equal(unname(d$margins[, "n1"]), unname(n1))
  expect_equal(unname(d$margins[, "n0"]), unname(n0))
  expect_equal(unname(d$margins[, "n"]), unname(rowSums(x)))

  expect_equal(unname(d$est[, "eta"]), unname(qlogis((x[, "n11"] + .5) / (n1 + 1))))
  expect_equal(unname(d$est[, "xi"]), unname(qlogis((x[, "n10"] + .5) / (n0 + 1))))
  expect_equal(unname(d$est[, "gamma"]), unname(qlogis((n1 + 1) / (rowSums(x) + 2))))
  expect_equal(
    unname(d$wvar[, "var_eta"]),
    unname(1 / (x[, "n11"] + .5) + 1 / (x[, "n01"] + .5))
  )
  expect_equal(
    unname(d$est),
    unname(qlogis(cbind(d$rates[, "SE"], 1 - d$rates[, "SP"], d$rates[, "PREV"])))
  )
})

test_that("counts are the raw table", {
  set.seed(1)
  x <- sim_data(20, th, rep(100, 20))
  d <- set_meta_data(x, CC = 0.5)
  n1 <- x[, "n11"] + x[, "n01"]

  expect_equal(unname(d$counts[, "y_eta"]), unname(x[, "n11"]))
  expect_equal(unname(d$counts[, "y_gamma"]), unname(n1))
  expect_equal(unname(d$counts[, "den_gamma"]), unname(rowSums(x)))
  expect_equal(d$counts, round(d$counts))
  # correction leaves counts alone
  expect_equal(suppressMessages(set_meta_data(x))$counts, d$counts)
})

test_that("input shape is forgiving, names are not", {
  set.seed(1)
  x <- sim_data(20, th, rep(100, 20))
  d <- set_meta_data(x, CC = 0.5)

  expect_equal(set_meta_data(x[, c(3, 1, 4, 2)], CC = 0.5), d)
  expect_equal(set_meta_data(as.data.frame(x), CC = 0.5), d)
  expect_message(u <- set_meta_data(unname(x), CC = 0.5), "n11, n10, n01, n00")
  expect_equal(unname(u$est), unname(d$est))

  colnames(x)[1] <- "wrong"
  expect_error(set_meta_data(x), "n11")
})

test_that("row names carry through", {
  set.seed(1)
  x <- sim_data(3, th, rep(100, 3))
  rownames(x) <- c("Smith", "Jones", "Lee")
  d <- set_meta_data(x, CC = 0.5)
  expect_identical(rownames(d$est), rownames(x))
  expect_identical(rownames(d$rates), rownames(x))
  expect_identical(rownames(d$counts), rownames(x))
})

test_that("bad input is rejected", {
  set.seed(1)
  x <- sim_data(20, th, rep(100, 20))
  expect_error(set_meta_data(x[, 1:3]))
  expect_error(set_meta_data(unname(x[, 1:3])), "4")
  expect_error(set_meta_data(x + 0.5), "round")
  expect_error(set_meta_data(replace(x, 1, -1L)), ">= 0")
  expect_error(set_meta_data(x, CC = -1), ">= 0")
  expect_error(set_meta_data(x, CC = c(0.5, 0.5)))
  expect_error(set_meta_data(matrix(0L, 1, 4)), "at least one observation")
  expect_identical(formals(set_meta_data)$CC, 0)
})

test_that("a correction keeps an empty cell finite", {
  set.seed(1)
  x <- sim_data(20, th, rep(100, 20))
  x[1, ] <- c(0L, 0L, 10L, 90L)
  d <- set_meta_data(x, CC = 0.5)
  expect_true(all(is.finite(d$est)))
  expect_true(all(is.finite(d$wvar)))
  expect_silent(set_meta_data(x, CC = 0.5))
  expect_message(bare <- set_meta_data(x, CC = 0), "infinite")
  expect_false(all(is.finite(bare$est)))
})

test_that("the empty cell message names studies", {
  set.seed(4)
  x <- sim_data(5, th, rep(100, 5))
  expect_silent(set_meta_data(x))
  x[c(2, 5), ] <- c(0L, 0L, 0L, 0L, 10L, 10L, 90L, 90L)
  expect_message(set_meta_data(x), "studies 2, 5")
  rownames(x) <- paste0("s", 1:5)
  expect_message(set_meta_data(x), "s2, s5")
  expect_message(set_meta_data(unname(x[2, , drop = FALSE])), "study 1")
})

test_that("fitting functions want the object", {
  set.seed(1)
  x <- sim_data(20, th, rep(100, 20))
  d <- set_meta_data(x, CC = 0.5)

  expect_error(init_theta(x), "accmeta_data")
  expect_error(fit_tlmm(x), "accmeta_data")
  expect_error(fit_tglmm(x), "accmeta_data")
  expect_error(init_theta(list()), "accmeta_data")
  expect_error(init_theta(unclass(d)), "accmeta_data")

  # dropped slots must not pass
  expect_error(init_theta(`[[<-`(d, "tab", NULL)), "is.matrix")
  expect_error(fit_tlmm(`[[<-`(d, "est", NULL)), "is.matrix")
  expect_error(fit_tglmm(`[[<-`(d, "counts", NULL)), "is.matrix")
})

test_that("fit_tlmm refuses infinite estimates", {
  set.seed(1)
  x <- sim_data(20, th, rep(100, 20))
  x[1, ] <- c(0L, 0L, 10L, 90L)
  expect_error(fit_tlmm(suppressMessages(set_meta_data(x))), "is.finite")
})

test_that("init_theta reads raw counts only", {
  set.seed(1)
  x <- sim_data(20, th, rep(100, 20))
  x[1, ] <- c(0L, 0L, 10L, 90L)
  st <- init_theta(suppressMessages(set_meta_data(x)))

  expect_true(all(is.finite(st)))
  # object correction is ignored
  expect_identical(init_theta(set_meta_data(x, CC = 1)), st)
  expect_identical(init_theta(set_meta_data(x, CC = 0.5)), st)
  expect_equal(st[1:3], unname(colMeans(set_meta_data(x, CC = 0.5)$est)))
  expect_error(
    suppressMessages(init_theta(suppressMessages(set_meta_data(x)), CC = 0)),
    "is.finite"
  )
})

test_that("init_theta gives a usable start", {
  set.seed(1)
  d <- set_meta_data(sim_data(20, th, rep(100, 20)), CC = 0.5)
  st <- init_theta(d)
  expect_length(st, 9)
  expect_true(all(is.finite(st)))
  expect_true(all(eigen(theta2list(st)$SIGMA, symmetric = TRUE)$values > 0))
  expect_equal(st[1:3], unname(colMeans(d$est)))
})

test_that("print returns its argument", {
  set.seed(1)
  d <- set_meta_data(sim_data(20, th, rep(100, 20)), CC = 0.5)
  expect_output(out <- print(d), "accmeta_data")
  expect_identical(out, d)
})
