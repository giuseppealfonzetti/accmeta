#' Moment-based starting value
#'
#' Builds a starting `THETA` from the data. The means are the column means of
#' the transformed estimates; the covariance is their sample covariance less the
#' average within-study variance, since the marginal covariance of the estimates
#' is \eqn{\Sigma_3 + \Gamma_{i,3}}. The result is symmetrised and its
#' eigenvalues floored to keep it positive-definite.
#'
#' @param DATA An `accmeta_data` object, as returned by [set_meta_data()]. Only
#'   its raw counts are read, so an object built without a continuity correction
#'   is fine here.
#' @param CC Continuity correction used to build the starting value, applied to
#'   the raw counts of `DATA` independently of the correction recorded on it.
#'   It exists to keep the logit transform finite when a study has an empty
#'   cell, and affects the starting value only, never the fit. The default keeps
#'   the starting value a fixed function of the counts, so that fits at
#'   different corrections all set out from the same point.
#' @param MIN_VAR Lower bound applied to the eigenvalues of the covariance.
#' @param SHRINK Proportion of the average within-study variance subtracted from
#'   the sample covariance. `0` leaves the sample covariance untouched.
#'
#' @return A numeric vector of length 9, in the layout documented for
#'   [theta2list()].
#'
#' @examples
#' th <- c(2.94, -2.2, -0.4, 0.0953, 0.4, -0.5108, 0.3, 0.2, -0.6931)
#' set.seed(1)
#' init_theta(set_meta_data(sim_data(15, th, rep(100, 15))))
#'
#' @export
init_theta <- function(DATA, CC = 0.5, MIN_VAR = 1e-4, SHRINK = 0.5) {
  stopifnot(
    inherits(DATA, "accmeta_data"),
    is.matrix(DATA$tab),
    is.numeric(MIN_VAR),
    MIN_VAR > 0,
    is.numeric(SHRINK),
    SHRINK >= 0
  )
  # rebuild with own correction
  dat <- set_meta_data(DATA$tab, CC = CC)
  stopifnot(is.finite(dat$est))
  est <- dat$est
  wvar <- dat$wvar

  # moment-based estimates
  S <- stats::cov(est) - SHRINK * diag(colMeans(wvar), 3, 3)
  S <- (S + t(S)) / 2
  eig <- eigen(S, symmetric = TRUE)
  S <- eig$vectors %*% diag(pmax(eig$values, MIN_VAR), 3, 3) %*% t(eig$vectors)
  out <- list2theta(list(MU = unname(colMeans(est)), SIGMA = (S + t(S)) / 2))
  return(out)
}
