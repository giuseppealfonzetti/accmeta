#' Unpack the working parameter vector
#'
#' Maps the working parameter vector `THETA` onto its interpretable
#' components: the random-effects mean vector and covariance matrix.
#'
#' @param THETA Numeric vector of length 9. Entries 1 to 3 are the random-effects
#'   means \eqn{(\bar\eta, \bar\xi, \bar\gamma)}. Entries 4 to 9 are the
#'   log-Cholesky factor of \eqn{\Sigma_3} in row-major lower-triangular order
#'   \eqn{(log(L_{11}), L_{21}, log(L_{22}), L_{31}, L_{32}, log(L_{33})}).
#'
#' @return A list with two components: `MU`, the length-3 mean vector, and
#'   `SIGMA`, the 3x3 positive-definite covariance matrix. Both are named after
#'   the logit-scale components `eta`, `xi` and `gamma`.
#'
#' @seealso [list2theta()] for the inverse map.
#'
#' @examples
#' theta2list(c(2.94, -2.2, -0.4, 0.0953, 0.4, -0.5108, 0.3, 0.2, -0.6931))
#'
#' @export
theta2list <- function(THETA) {
  stopifnot(is.numeric(THETA), length(THETA) == 9)
  nm <- c("eta", "xi", "gamma")
  L <- matrix(0, 3, 3)
  L[lower.tri(L, diag = TRUE)] <- THETA[c(4, 5, 7, 6, 8, 9)]
  diag(L) <- exp(diag(L))
  mu <- THETA[1:3]
  names(mu) <- nm
  sigma <- tcrossprod(L)
  dimnames(sigma) <- list(nm, nm)
  list(MU = mu, SIGMA = sigma)
}


#' Pack the working parameter vector
#'
#' Inverse of [theta2list()]: packs the random-effects mean vector and covariance
#' matrix into the working parameter vector `THETA`.
#'
#' @param LIST A list with components `MU`, a numeric vector of length 3, and
#'   `SIGMA`, a 3x3 symmetric positive-definite matrix.
#'
#' @return A numeric vector of length 9, in the layout documented for
#'   [theta2list()].
#'
#' @seealso [theta2list()] for the inverse map.
#'
#' @examples
#' nm <- c("eta", "xi", "gamma")
#' li <- list(
#'   MU = setNames(c(2.94, -2.2, -0.4), nm),
#'   SIGMA = matrix(c(1.21, 0.44, 0.33,
#'                    0.44, 0.52, 0.24,
#'                    0.33, 0.24, 0.38), 3, 3, dimnames = list(nm, nm))
#' )
#' list2theta(li)
#' all.equal(theta2list(list2theta(li)), li)
#'
#' @export
list2theta <- function(LIST) {
  stopifnot(
    is.numeric(LIST$MU),
    length(LIST$MU) == 3,
    identical(dim(LIST$SIGMA), c(3L, 3L)),
    isSymmetric(LIST$SIGMA)
  )
  L <- t(chol(LIST$SIGMA))
  diag(L) <- log(diag(L))
  # THETA carries no names
  unname(c(LIST$MU, L[lower.tri(L, diag = TRUE)][c(1, 2, 4, 3, 5, 6)]))
}

#' Set the prior on the random-effects covariance
#'
#' Control the prior object for \eqn{\Sigma_3 \sim W(\nu, A I_3)} to be passed to [fit_tlmm()] or [fit_tglmm()] through the `PRIOR` argument. The prior is
#' \eqn{\Sigma_3 \sim W(\nu, A I_3)}. Setting \eqn{\nu=4} and \eqn{A=Inf} correspond to maximum likelihood estimation.
#'
#' @param DEGREES Degrees of freedom \eqn{\nu} of the Wishart prior. The default
#'   keeps the covariance off the boundary of the parameter space, which the
#'   unpenalised fit reaches on small or sparse data. Use `set_prior(4)` for the
#'   flat prior.
#' @param SCALE Scale \eqn{A}. Represents a soft ceiling on the variance scale.
#'
#' @return An object of class `accmeta_prior`: a list with `DEGREES` and
#'   `SCALE`, held as given and passed to the template unchanged.
#'
#' @examples
#' set_prior()
#' set_prior(DEGREES = 4)
#' set_prior(DEGREES = 5, SCALE = 100)
#'
#' @export
set_prior <- function(DEGREES = 5, SCALE = Inf) {
  stopifnot(
    is.numeric(DEGREES),
    length(DEGREES) == 1,
    DEGREES >= 4,
    DEGREES == round(DEGREES),
    is.numeric(SCALE),
    length(SCALE) == 1,
    SCALE > 0
  )
  out <- list(DEGREES = DEGREES, SCALE = SCALE)
  class(out) <- "accmeta_prior"
  return(out)
}

#' Map the working vector to Joe's unconstrained parameters
#'
#' @param THETA Numeric vector of length 9 (see [theta2list()]).
#'
#' @return A numeric vector of length 9:
#'   entries 1 to 3 the means (copied verbatim), 4 to 6 the log marginal SDs
#'   \eqn{\log\sqrt{\mathrm{diag}\,\Sigma_3}}, 7 to 9 the Fisher-z values
#'   \eqn{(\mathrm{atanh}\,\rho_{12}, \mathrm{atanh}\,\rho_{13},
#'   \mathrm{atanh}\,\rho_{23\mid1})}. Correlations are clamped to
#'   \eqn{\pm(1 - 10^{-10})} before `atanh`, so `joe2theta(theta2joe())` is
#'   exact for interior values but not bit-exact at a \eqn{\pm 1} boundary.
#'
#' @seealso [joe2theta()] for the inverse map.
#' @export
theta2joe <- function(THETA) {
  li <- theta2list(THETA)
  s <- sqrt(diag(li$SIGMA))
  cormat <- li$SIGMA / tcrossprod(s)
  r12 <- cormat[1, 2]
  r13 <- cormat[1, 3]
  r23 <- cormat[2, 3]
  r23g1 <- (r23 - r12 * r13) / sqrt((1 - r12^2) * (1 - r13^2))
  r <- pmin(pmax(c(r12, r13, r23g1), -1 + 1e-10), 1 - 1e-10)
  unname(c(li$MU, log(s), atanh(r)))
}

#' Map Joe's unconstrained parameters back to the working vector
#'
#' @param JOEPAR Numeric vector of length 9 in the layout returned by
#'   [theta2joe()].
#'
#' @return A numeric vector of length 9 (the working vector `THETA`).
#'
#' @seealso [theta2joe()] for the inverse map.
#' @export
joe2theta <- function(JOEPAR) {
  stopifnot(is.numeric(JOEPAR), length(JOEPAR) == 9)
  s <- exp(JOEPAR[4:6])
  z <- tanh(JOEPAR[7:9])
  r12 <- z[1]
  r13 <- z[2]
  r23g1 <- z[3]
  r23 <- r12 * r13 + r23g1 * sqrt((1 - r12^2) * (1 - r13^2))
  cormat <- matrix(c(1, r12, r13, r12, 1, r23, r13, r23, 1), 3, 3)
  list2theta(list(MU = JOEPAR[1:3], SIGMA = cormat * tcrossprod(s)))
}


#' Set the control parameters for the iterative bootstrap
#'
#' Build the `CONTROL` object passed to [fit_ib()].
#'
#' @param H Number of datasets simulated per iteration.
#' @param MAX_ITER Maximum number of iterations.
#' @param TOL Significance level \eqn{\alpha \in (0, 1)} of the convergence
#'   test (see `TERMINATION`). If `NULL`, `0.5` for `"hotelling"` and `0.05`
#'   for `"confidence"`. Ignored for `"max_iter"`.
#' @param STEP Damping factor \eqn{\gamma \in (0, 1]}.
#' @param PATIENCE Stop after this many iterations with no improvement in the
#'   best convergence statistic (returning the best iterate seen). Ignored for
#'   `"max_iter"`.
#' @param TERMINATION Convergence test on the Hotelling \eqn{T^2} statistc constructed from the IB gap:
#'   `"hotelling"` stops when \eqn{H_0: E[gap] = 0} is not rejected;
#'   `"confidence"` (default) stops when the
#'   \eqn{1 - \alpha} Hotelling confidence ellipsoid of the gap lies inside the
#'   tolerance region \eqn{\delta^\top \Sigma^{-1} \delta \le} `PRECISION`,
#'   with \eqn{\Sigma} the covariance of a single simulated estimate;
#'   `"max_iter"` runs all `MAX_ITER` iterations with no test and returns the
#'   best iterate seen.
#' @param PRECISION Squared equivalence margin of the `"confidence"` test, in
#'   squared standard deviations of the estimator.
#' @param BOOST If `TRUE`, increases `H` when the `"confidence"` test fails.
#'   Only for `TERMINATION = "confidence"`.
#' @param BOOST_FACTOR Maximum multiplicative growth of `H` when `BOOST = TRUE`.
#' @param MAX_H Upper bound on `H` when `BOOST = TRUE`.
#' @param NCORES Number of cores to be passed to [mirai::daemons()]. Default `1` runs serially.
#'   Values `> 1` for parallel computations across the H simulated datasets at each iteration.
#'
#' @return An object of class `accmeta_ib_control` to be passed to [fit_ib()] via `CONTROL` argumnet.
#'
#' @seealso [fit_ib()].
#'
#' @examples
#' set_ib_control()
#'
#' @export
set_ib_control <- function(
  H = 100,
  MAX_ITER = 100,
  TOL = NULL,
  STEP = 1,
  PATIENCE = 5L,
  TERMINATION = c("confidence", "hotelling", "max_iter"),
  PRECISION = 1,
  BOOST = FALSE,
  BOOST_FACTOR = 2,
  MAX_H = 500,
  NCORES = 1L
) {
  TERMINATION <- match.arg(TERMINATION)
  if (is.null(TOL)) {
    TOL <- if (TERMINATION == "hotelling") 0.5 else 0.05
  }
  stopifnot(
    is.numeric(H),
    length(H) == 1,
    H > 9,
    is.numeric(MAX_ITER),
    length(MAX_ITER) == 1,
    MAX_ITER >= 1,
    is.numeric(TOL),
    length(TOL) == 1,
    TOL > 0,
    TOL < 1,
    is.numeric(STEP),
    length(STEP) == 1,
    STEP > 0,
    STEP <= 1,
    is.numeric(PATIENCE),
    length(PATIENCE) == 1,
    PATIENCE >= 1,
    is.numeric(PRECISION),
    length(PRECISION) == 1,
    PRECISION > 0,
    isTRUE(BOOST) || isFALSE(BOOST),
    !BOOST || TERMINATION == "confidence",
    is.numeric(BOOST_FACTOR),
    length(BOOST_FACTOR) == 1,
    BOOST_FACTOR > 1,
    is.numeric(MAX_H),
    length(MAX_H) == 1,
    MAX_H >= H,
    is.numeric(NCORES),
    length(NCORES) == 1,
    NCORES >= 1
  )
  h_top <- if (BOOST) MAX_H else H
  if (
    TERMINATION == "confidence" &&
      (h_top - 1) * PRECISION <= ib_crit(h_top, TOL)
  ) {
    warning(
      "the confidence test cannot pass with H = ",
      h_top,
      " and PRECISION = ",
      PRECISION,
      " ((H - 1) * PRECISION must exceed ",
      signif(ib_crit(h_top, TOL), 3),
      ").",
      call. = FALSE
    )
  }
  out <- list(
    H = H,
    MAX_ITER = MAX_ITER,
    TOL = TOL,
    STEP = STEP,
    PATIENCE = PATIENCE,
    TERMINATION = TERMINATION,
    PRECISION = PRECISION,
    BOOST = BOOST,
    BOOST_FACTOR = BOOST_FACTOR,
    MAX_H = MAX_H,
    NCORES = NCORES
  )
  class(out) <- "accmeta_ib_control"
  return(out)
}
