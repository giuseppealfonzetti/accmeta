#' Fit the asymptotic normal approximation model
#'
#' Trivariate linear mixed model. A within-study normal approximation
#' makes the likelihood available in closed form.
#'
#' @param DATA An `accmeta_data` object, as returned by [set_meta_data()].
#' @param THETA_START Numeric vector of length 9 giving the starting value. If
#'   `NULL`, [init_theta()] is used.
#' @param PRIOR Prior on random effects covariance matrix, as returned by
#'   [set_prior()]. Under any prior but the flat one, `NLL` is a penalised
#'   objective rather than a log-likelihood.
#' @param CONTROL List of control parameters passed to [ucminf::ucminf()]; see
#'   its documentation for the accepted entries.
#'
#' @param WORKPAR Coordinates the optimiser works in for the covariance:
#'   `"PinheiroBates"` (log-Cholesky, default) or `"Joe"` (log standard
#'   deviations and Fisher-z partial correlations, see [theta2joe()]). The
#'   objective is a function of \eqn{\Sigma_3} alone, so the fit is the same
#'   under any prior with `DEGREES > 4`; under the flat prior the two can stop
#'   at different points of a boundary fit. `THETA` is returned in the layout of
#'   [theta2list()] either way, while `OBJ` takes `WORKPAR` coordinates: under
#'   `"Joe"`, evaluate it at `theta2joe(THETA)`.
#'
#' @return A list with components `THETA`, the fitted parameter vector in the
#'   layout documented for [theta2list()]; `CONVERGENCE`, the optimiser
#'   convergence code; `NLL`, the negative log-likelihood at the optimum;
#'   `OBJ`, the TMB object, retained for [TMB::sdreport()]; and `PRIOR`, the
#'   prior used.
#'
#' @examples
#' th <- c(2.94, -2.2, -0.4, 0.0953, 0.4, -0.5108, 0.3, 0.2, -0.6931)
#' set.seed(1)
#' fit_tlmm(set_meta_data(sim_data(50, th, rep(100, 50)), CC = 0.5))$THETA
#'
#' @export
fit_tlmm <- function(
  DATA,
  THETA_START = NULL,
  PRIOR = set_prior(),
  CONTROL = list(maxeval = 1000),
  WORKPAR = c("PinheiroBates", "Joe")
) {
  WORKPAR <- match.arg(WORKPAR)
  stopifnot(
    inherits(DATA, "accmeta_data"),
    is.matrix(DATA$est),
    is.finite(DATA$est),
    inherits(PRIOR, "accmeta_prior")
  )
  if (is.null(THETA_START)) {
    THETA_START <- init_theta(DATA)
  }
  stopifnot(is.numeric(THETA_START), length(THETA_START) == 9)
  obj <- TMB::MakeADFun(
    data = list(
      MODEL = "tlmm",
      EST = DATA$est,
      WVAR = DATA$wvar,
      DEGREES = as.numeric(PRIOR$DEGREES),
      SCALE = as.numeric(PRIOR$SCALE),
      WORKPAR = WORKPAR
    ),
    parameters = list(
      MU = THETA_START[1:3],
      ALPHA = if (WORKPAR == "Joe") {
        theta2joe(THETA_START)[4:9]
      } else {
        THETA_START[4:9]
      }
    ),
    DLL = "accmeta",
    silent = TRUE
  )
  est <- ucminf::ucminf(
    par = obj$par,
    fn = obj$fn,
    gr = obj$gr,
    control = CONTROL
  )
  out <- list(
    THETA = unname(c(est$par[1:3], obj$report(est$par)$LOGCHOL)),
    CONVERGENCE = est$convergence,
    NLL = est$value,
    OBJ = obj,
    PRIOR = PRIOR
  )
  class(out) <- c("accmeta_tlmm", "accmeta_fit")
  return(out)
}
