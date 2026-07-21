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
#' @return A list with components `THETA`, the fitted parameter vector in the
#'   layout documented for [theta2list()]; `CONVERGENCE`, the optimiser
#'   convergence code; `NLL`, the negative log-likelihood at the optimum; and
#'   `OBJ`, the TMB object, retained for [TMB::sdreport()].
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
  CONTROL = list(maxeval = 1000)
) {
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
      SCALE = as.numeric(PRIOR$SCALE)
    ),
    parameters = list(MU = THETA_START[1:3], ALPHA = THETA_START[4:9]),
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
    THETA = unname(est$par),
    CONVERGENCE = est$convergence,
    NLL = est$value,
    OBJ = obj
  )
  class(out) <- c("accmeta_tlmm", "accmeta_fit")
  return(out)
}
