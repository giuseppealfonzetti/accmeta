#' Compute a starting value
#'
#' @param DATA An `accmeta_data` object from [set_meta_data()].
#' @param CC Non-negative numeric, continuity correction used for the starting
#'   value only.
#' @param MIN_VAR Positive numeric, lower bound on the covariance eigenvalues.
#' @param SHRINK Non-negative numeric, share of the mean within-study variance
#'   subtracted from the sample covariance.
#' @param MODEL `"trivariate"` or `"bivariate"`.
#'
#' @return Numeric vector of length 9 (trivariate) or 5 (bivariate), in the
#'   layout of [theta2list()].
#'
#' @examples
#' th <- c(2.94, -2.2, -0.4, 0.0953, 0.4, -0.5108, 0.3, 0.2, -0.6931)
#' set.seed(1)
#' init_theta(set_meta_data(sim_data(15, th, rep(100, 15))))
#'
#' @export
init_theta <- function(
  DATA,
  CC = 0.5,
  MIN_VAR = 1e-4,
  SHRINK = 0.5,
  MODEL = c("trivariate", "bivariate")
) {
  MODEL <- match.arg(MODEL)
  k <- if (MODEL == "bivariate") 1:2 else 1:3
  stopifnot(
    inherits(DATA, "accmeta_data"),
    is.matrix(DATA$tab),
    is.numeric(MIN_VAR),
    MIN_VAR > 0,
    is.numeric(SHRINK),
    SHRINK >= 0
  )
  # rebuild with own correction
  dat <- set_meta_data(DATA$tab, CC = CC, DESIGN = DATA$design)
  stopifnot(
    all(is.finite(dat$est) | is.na(dat$est)),
    MODEL == "bivariate" || sum(!is.na(dat$est[, 3])) >= 2
  )
  est <- dat$est[, k, drop = FALSE]
  wvar <- dat$wvar[, k, drop = FALSE]

  # moment-based estimates
  S <- stats::cov(est, use = if (anyNA(est)) "pairwise.complete.obs" else "everything") -
    SHRINK * diag(colMeans(wvar, na.rm = TRUE), length(k), length(k))
  S <- (S + t(S)) / 2
  eig <- eigen(S, symmetric = TRUE)
  S <- eig$vectors %*% diag(pmax(eig$values, MIN_VAR), length(k), length(k)) %*% t(eig$vectors)
  out <- list2theta(list(MU = unname(colMeans(est, na.rm = TRUE)), SIGMA = (S + t(S)) / 2))
  return(out)
}
