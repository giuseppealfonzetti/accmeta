#' Unpack the parameter vector
#'
#' @param THETA Numeric vector of length 9 (means of `eta`, `xi`, `gamma`, then
#'   log-Cholesky entries `log L11`, `L21`, `log L22`, `L31`, `L32`, `log L33`),
#'   5 (means of `eta`, `xi`, then `log L11`, `L21`, `log L22`) or 12 (skew-normal:
#'   the length-9 layout plus three slants).
#'
#' @return List with `MU` (named numeric, length 3 or 2) and `SIGMA` (named
#'   covariance matrix); for length 12 also `SKEW` (named numeric, length 3).
#'
#' @seealso [list2theta()] for the inverse map.
#'
#' @examples
#' theta2list(c(2.94, -2.2, -0.4, 0.0953, 0.4, -0.5108, 0.3, 0.2, -0.6931))
#'
#' @export
theta2list <- function(THETA) {
  stopifnot(is.numeric(THETA), length(THETA) %in% c(5, 9, 12))
  if (length(THETA) == 5) {
    li <- theta2list(embed(THETA, "PinheiroBates"))
    return(list(MU = li$MU[1:2], SIGMA = li$SIGMA[1:2, 1:2]))
  }
  nm <- c("eta", "xi", "gamma")
  L <- matrix(0, 3, 3)
  L[lower.tri(L, diag = TRUE)] <- THETA[c(4, 5, 7, 6, 8, 9)]
  diag(L) <- exp(diag(L))
  mu <- THETA[1:3]
  names(mu) <- nm
  sigma <- tcrossprod(L)
  dimnames(sigma) <- list(nm, nm)
  if (length(THETA) == 9) {
    return(list(MU = mu, SIGMA = sigma))
  }
  delta <- THETA[10:12] * sqrt(diag(sigma))
  vc <- sigma + (1 - 2 / pi) * tcrossprod(delta)
  m <- sqrt(2 / pi) * delta / sqrt(diag(vc))
  list(
    MU = mu + sqrt(2 / pi) * delta,
    SIGMA = vc,
    SKEW = stats::setNames((4 - pi) / 2 * m^3, nm)
  )
}


#' Pack the parameter vector
#'
#' @param LIST List with `MU` (numeric, length 3 or 2), `SIGMA` (symmetric
#'   positive-definite matrix, 3 x 3 or 2 x 2) and optionally `SKEW` (numeric,
#'   length 3).
#'
#' @return Numeric vector of length 9, 5 or 12, in the layout of [theta2list()].
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
  if (length(LIST$MU) == 2 && is.null(LIST$SKEW)) {
    stopifnot(identical(dim(LIST$SIGMA), c(2L, 2L)), isSymmetric(LIST$SIGMA))
    s <- diag(3)
    s[1:2, 1:2] <- LIST$SIGMA
    out <- list2theta(list(MU = c(LIST$MU, 0), SIGMA = s))
    return(project(out, "PinheiroBates"))
  }
  stopifnot(
    is.numeric(LIST$MU),
    length(LIST$MU) == 3,
    identical(dim(LIST$SIGMA), c(3L, 3L)),
    isSymmetric(LIST$SIGMA),
    is.null(LIST$SKEW) || (is.numeric(LIST$SKEW) && length(LIST$SKEW) == 3)
  )
  if (is.null(LIST$SKEW)) {
    L <- t(chol(LIST$SIGMA))
    diag(L) <- log(diag(L))
    # THETA carries no names
    return(unname(c(
      LIST$MU,
      L[lower.tri(L, diag = TRUE)][c(1, 2, 4, 3, 5, 6)]
    )))
  }
  delta <- skew_delta(LIST$SIGMA, LIST$SKEW)
  t <- skew_t(LIST$SIGMA, LIST$SKEW)
  if (!is.finite(t) || t >= 1) {
    stop("skewness not admissible for this covariance", call. = FALSE)
  }
  psi <- LIST$SIGMA - (1 - 2 / pi) * tcrossprod(delta)
  c(
    list2theta(list(
      MU = LIST$MU - sqrt(2 / pi) * delta,
      SIGMA = (psi + t(psi)) / 2
    )),
    unname(delta / sqrt(diag(psi)))
  )
}

skew_delta <- function(SIGMA, SKEW) {
  sign(SKEW) *
    sqrt(diag(SIGMA)) *
    sqrt(pi / 2) *
    (2 * abs(SKEW) / (4 - pi))^(1 / 3)
}

skew_t <- function(SIGMA, SKEW) {
  delta <- skew_delta(SIGMA, SKEW)
  (1 - 2 / pi) * sum(delta * solve(SIGMA, delta))
}

#' Skew-normal direct parameters
#'
#' @param THETA Numeric vector of length 9 or 12, in the layout of
#'   [theta2list()].
#'
#' @return List with `xi`, `Omega` and `alpha`, as used by [sn::dmsn()]; `alpha`
#'   is zero for length 9.
#'
#' @export
theta2dp <- function(THETA) {
  li <- theta2list(THETA[1:9])
  if (length(THETA) == 9) {
    return(list(xi = li$MU, Omega = li$SIGMA, alpha = 0 * li$MU))
  }
  delta <- THETA[10:12] * sqrt(diag(li$SIGMA))
  om <- li$SIGMA + tcrossprod(delta)
  oid <- solve(om, delta)
  list(
    xi = li$MU,
    Omega = om,
    alpha = sqrt(diag(om)) *
      oid /
      sqrt(max(1 - sum(delta * oid), .Machine$double.eps))
  )
}

#' Set the prior on the random-effects covariance
#'
#' @param DEGREES Integer, at least 3, Wishart degrees of freedom, or `NULL` for
#'   the number of random effects plus 2.
#' @param SCALE Positive numeric, Wishart scale. `Inf` with `DEGREES` equal to the
#'   number of random effects plus 1 gives maximum likelihood.
#'
#' @return Object of class `accmeta_prior`, a list with `DEGREES` and `SCALE`.
#'
#' @examples
#' set_prior()
#' set_prior(DEGREES = 4)
#' set_prior(DEGREES = 5, SCALE = 100)
#'
#' @export
set_prior <- function(DEGREES = NULL, SCALE = 100) {
  if (!is.null(DEGREES)) {
    stopifnot(
      is.numeric(DEGREES),
      length(DEGREES) == 1,
      DEGREES >= 3,
      DEGREES == round(DEGREES)
    )
  }
  stopifnot(
    is.numeric(SCALE),
    length(SCALE) == 1,
    SCALE > 0
  )
  out <- list(DEGREES = DEGREES, SCALE = SCALE)
  class(out) <- "accmeta_prior"
  return(out)
}

#' Map the parameter vector to Joe coordinates
#'
#' @param THETA Numeric vector of length 9, 5 or 12, in the layout of
#'   [theta2list()].
#'
#' @return Numeric vector of the same length: means, log standard deviations,
#'   Fisher-z correlations, then any skew entries unchanged.
#'
#' @seealso [joe2theta()] for the inverse map.
#'
#' @export
theta2joe <- function(THETA) {
  if (length(THETA) == 5) {
    return(project(theta2joe(embed(THETA, "PinheiroBates")), "Joe"))
  }
  li <- theta2list(THETA[1:9])
  s <- sqrt(diag(li$SIGMA))
  cormat <- li$SIGMA / tcrossprod(s)
  r <- pmin(pmax(c(cormat[1, 2], cormat[1, 3]), -1 + 1e-10), 1 - 1e-10)
  r23g1 <- (cormat[2, 3] - r[1] * r[2]) / sqrt((1 - r[1]^2) * (1 - r[2]^2))
  r <- pmin(pmax(c(r, r23g1), -1 + 1e-10), 1 - 1e-10)
  unname(c(li$MU, log(s), atanh(r), THETA[-(1:9)]))
}

#' Map Joe coordinates to the parameter vector
#'
#' @param JOEPAR Numeric vector of length 9, 5 or 12, in the layout returned by
#'   [theta2joe()].
#'
#' @return Numeric vector of the same length, in the layout of [theta2list()].
#'
#' @seealso [theta2joe()] for the inverse map.
#'
#' @export
joe2theta <- function(JOEPAR) {
  stopifnot(is.numeric(JOEPAR), length(JOEPAR) %in% c(5, 9, 12))
  if (length(JOEPAR) == 5) {
    return(project(joe2theta(embed(JOEPAR, "Joe")), "PinheiroBates"))
  }
  s <- exp(JOEPAR[4:6])
  z <- tanh(JOEPAR[7:9])
  r12 <- z[1]
  r13 <- z[2]
  r23g1 <- z[3]
  r23 <- r12 * r13 + r23g1 * sqrt((1 - r12^2) * (1 - r13^2))
  cormat <- matrix(c(1, r12, r13, r12, 1, r23, r13, r23, 1), 3, 3)
  c(
    list2theta(list(MU = JOEPAR[1:3], SIGMA = cormat * tcrossprod(s))),
    JOEPAR[-(1:9)]
  )
}


#' Set the iterative bootstrap options
#'
#' @param H Integer greater than 9, simulated datasets per iteration.
#' @param MAX_ITER Positive integer, maximum iterations, or `NULL` for 100.
#' @param STEP Numeric in (0, 1], damping factor.
#' @param PATIENCE Positive integer, iterations without improvement before
#'   stopping.
#' @param TERMINATION `"mahalanobis"` (convergence test) or `"max_iter"` (run
#'   all iterations).
#' @param PRECISION Positive numeric, tolerance of the `"mahalanobis"` test.
#' @param NCORES Positive integer, parallel workers; `1` runs serially.
#' @param VCOV_H Integer greater than 18, simulated datasets for the Jacobian.
#' @param VCOV_N_NODES Integer, at least 2, quadrature nodes per dimension for
#'   the Jacobian.
#' @param VCOV_N_ITER Positive integer, Newton steps for the Jacobian.
#' @param BOX_SD Numeric vector of length 2, lower and upper bound of the
#'   standard deviations.
#' @param BOX_COR Numeric in (0, 1), bound on the absolute correlations.
#' @param BOX_MU Positive numeric, bound on the absolute means.
#'
#' @return Object of class `accmeta_ib_control`.
#'
#' @seealso [fit_ib()].
#'
#' @examples
#' set_ib_control()
#'
#' @export
set_ib_control <- function(
  H = 100,
  MAX_ITER = NULL,
  STEP = 1,
  PATIENCE = 30L,
  TERMINATION = c("mahalanobis", "max_iter"),
  PRECISION = 1,
  NCORES = 1L,
  VCOV_H = 1000,
  VCOV_N_NODES = 7L,
  VCOV_N_ITER = 10L,
  BOX_SD = c(0.1, 10),
  BOX_COR = 0.99,
  BOX_MU = 20
) {
  TERMINATION <- match.arg(TERMINATION)
  stopifnot(
    is.numeric(H),
    length(H) == 1,
    H > 9,
    is.null(MAX_ITER) ||
      (is.numeric(MAX_ITER) && length(MAX_ITER) == 1 && MAX_ITER >= 1),
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
    is.numeric(NCORES),
    length(NCORES) == 1,
    NCORES >= 1,
    is.numeric(VCOV_H),
    length(VCOV_H) == 1,
    VCOV_H > 18,
    is.numeric(VCOV_N_NODES),
    length(VCOV_N_NODES) == 1,
    VCOV_N_NODES >= 2,
    is.numeric(VCOV_N_ITER),
    length(VCOV_N_ITER) == 1,
    VCOV_N_ITER >= 1,
    is.numeric(BOX_SD),
    length(BOX_SD) == 2,
    BOX_SD[1] > 0,
    BOX_SD[1] < BOX_SD[2],
    is.numeric(BOX_COR),
    length(BOX_COR) == 1,
    BOX_COR > 0,
    BOX_COR < 1,
    is.numeric(BOX_MU),
    length(BOX_MU) == 1,
    BOX_MU > 0
  )
  out <- list(
    H = H,
    MAX_ITER = MAX_ITER,
    STEP = STEP,
    PATIENCE = PATIENCE,
    TERMINATION = TERMINATION,
    PRECISION = PRECISION,
    NCORES = NCORES,
    VCOV_H = VCOV_H,
    VCOV_N_NODES = VCOV_N_NODES,
    VCOV_N_ITER = VCOV_N_ITER,
    BOX_SD = BOX_SD,
    BOX_COR = BOX_COR,
    BOX_MU = BOX_MU
  )
  class(out) <- "accmeta_ib_control"
  return(out)
}

free_idx <- function(WORKPAR) {
  if (WORKPAR == "Joe") c(1, 2, 4, 5, 7) else c(1, 2, 4, 5, 6)
}

embed <- function(X, WORKPAR) {
  replace(numeric(9), free_idx(WORKPAR), X)
}

project <- function(X, WORKPAR) {
  X[free_idx(WORKPAR)]
}

theta_model <- function(THETA) {
  if (length(THETA) == 5) "bivariate" else "trivariate"
}

flat_prior <- function(THETA) {
  set_prior(if (length(THETA) == 5) 3 else 4, Inf)
}

resolve_prior <- function(PRIOR, P) {
  stopifnot(inherits(PRIOR, "accmeta_prior"))
  if (is.null(PRIOR$DEGREES)) {
    PRIOR$DEGREES <- P + 2
  }
  if (PRIOR$DEGREES < P + 1) {
    stop(
      "DEGREES must be at least ",
      P + 1,
      " for ",
      P,
      " random effects.",
      call. = FALSE
    )
  }
  PRIOR
}

make_obj <- function(
  DATA,
  THETA,
  PRIOR,
  WORKPAR,
  MODEL,
  N_NODES = 15L,
  N_ITER = 10L
) {
  biv <- length(THETA) == 5
  work <- ib_to_work(THETA, WORKPAR)
  map <- list()
  if (biv) {
    work <- embed(work, WORKPAR)
    free <- seq_len(9) %in% free_idx(WORKPAR)
    map <- list(
      MU = factor(c(1, 2, NA)),
      ALPHA = factor(ifelse(free[4:9], seq_len(6), NA))
    )
    DATA$counts[, c("y_gamma", "den_gamma")] <- 0
    DATA$est[, "gamma"] <- NA
    DATA$wvar[, "var_gamma"] <- NA
  }
  cohort <- if (biv || is.null(DATA$design)) {
    rep(as.integer(!biv), DATA$n_studies)
  } else {
    as.integer(DATA$design == "cohort")
  }
  prior <- list(
    DEGREES = as.numeric(PRIOR$DEGREES) + biv,
    SCALE = as.numeric(PRIOR$SCALE),
    WORKPAR = WORKPAR,
    COHORT = cohort
  )
  data <- if (MODEL == "tlmm") {
    c(
      list(
        MODEL = "tlmm",
        EST = replace(DATA$est, is.na(DATA$est), 0),
        WVAR = replace(DATA$wvar, is.na(DATA$wvar), 0)
      ),
      prior
    )
  } else {
    gh <- statmod::gauss.quad(N_NODES, "hermite")
    c(
      list(
        MODEL = "tglmm",
        y = as.integer(t(DATA$counts[, 1:3, drop = FALSE])),
        den = as.integer(t(DATA$counts[, 4:6, drop = FALSE])),
        f = as.factor(rep(seq_len(DATA$n_studies), each = 3)),
        niter = as.integer(N_ITER),
        ws = as.numeric(exp(gh$nodes^2) * gh$weights),
        z = as.numeric(gh$nodes)
      ),
      prior
    )
  }
  TMB::MakeADFun(
    data = data,
    parameters = list(MU = work[1:3], ALPHA = work[4:9]),
    map = map,
    DLL = "accmeta",
    silent = TRUE
  )
}
