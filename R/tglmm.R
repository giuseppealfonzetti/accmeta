#' Fit the exact model
#'
#' Trivariate generalized linear mixed model. The
#' likelihood has no closed form, so each study contributes a trivariate integral
#' evaluated by adaptive Gauss-Hermite quadrature.
#'
#' @param DATA An `accmeta_data` object, as returned by [set_meta_data()].
#' @param THETA_START Numeric vector of length 9 giving the starting value. If
#'   `NULL`, [init_theta()] is used.
#' @param N_NODES Number of quadrature nodes per dimension. The integrand is
#'   evaluated `N_NODES^3` times per study per likelihood evaluation.
#' @param N_ITER Number of Newton steps used to locate the mode of the integrand.
#' @param PRIOR Prior on \eqn{\Sigma_3}, as returned by [set_prior()]. Under any
#'   prior but the flat one, `NLL` is a penalised objective rather than a
#'   log-likelihood.
#' @param CONTROL List of control parameters passed to [ucminf::ucminf()]; see
#'   its documentation for the accepted entries. The default raises `maxeval`
#'   above the ucminf default, which this likelihood exhausts before converging.
#'
#' @param WORKPAR Coordinates the optimiser works in for the covariance:
#'   `"PinheiroBates"` (log-Cholesky, default) or `"Joe"` (log standard
#'   deviations and Fisher-z partial correlations, see [theta2joe()]).
#' @return A list with components `THETA`, the fitted parameter vector in the
#'   layout documented for [theta2list()]; `CONVERGENCE`, the optimiser
#'   convergence code; `NLL`, the negative log-likelihood at the optimum;
#'   `OBJ`, the TMB object, retained for [TMB::sdreport()]; `PRIOR`, the prior
#'   used; and `DATA`, `WORKPAR`, `N_NODES` and `N_ITER`, the inputs, retained
#'   for [vcov()][vcov.accmeta_fit].
#'
#' @examples
#' \dontrun{
#' # N_NODES^3 integrand evaluations per study per likelihood call
#' th <- c(2.94, -2.2, -0.4, 0.0953, 0.4, -0.5108, 0.3, 0.2, -0.6931)
#' set.seed(1)
#' fit_tglmm(set_meta_data(sim_data(15, th, rep(100, 15))))$THETA
#' }
#'
#' @export
fit_tglmm <- function(
  DATA,
  THETA_START = NULL,
  N_NODES = 15L,
  N_ITER = 10L,
  PRIOR = set_prior(),
  CONTROL = list(maxeval = 1000),
  WORKPAR = c("PinheiroBates", "Joe")
) {
  WORKPAR <- match.arg(WORKPAR)
  stopifnot(
    inherits(DATA, "accmeta_data"),
    is.matrix(DATA$counts),
    inherits(PRIOR, "accmeta_prior")
  )
  if (is.null(THETA_START)) {
    THETA_START <- init_theta(DATA)
  }
  stopifnot(
    is.numeric(THETA_START),
    length(THETA_START) == 9,
    is.numeric(N_NODES),
    length(N_NODES) == 1,
    N_NODES >= 2,
    is.numeric(N_ITER),
    length(N_ITER) == 1,
    N_ITER >= 1
  )
  obj <- make_obj(DATA, THETA_START, PRIOR, WORKPAR, "tglmm", N_NODES, N_ITER)
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
    PRIOR = PRIOR,
    DATA = DATA,
    WORKPAR = WORKPAR,
    N_NODES = N_NODES,
    N_ITER = N_ITER
  )
  class(out) <- c("accmeta_tglmm", "accmeta_fit")
  return(out)
}
