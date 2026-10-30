#' Fit the trivariate generalized linear mixed model
#'
#' @param DATA An `accmeta_data` object from [set_meta_data()].
#' @param THETA_START Numeric vector of length 9, or `NULL` for [init_theta()].
#' @param N_NODES Integer, at least 2, quadrature nodes per dimension.
#' @param N_ITER Positive integer, Newton steps for the integrand mode.
#' @param PRIOR An `accmeta_prior` object from [set_prior()].
#' @param CONTROL List passed to [ucminf::ucminf()].
#' @param WORKPAR `"PinheiroBates"` or `"Joe"`, covariance coordinates.
#'
#' @return Object of class `accmeta_tglmm`, a list with `THETA` (numeric, length
#'   9), `CONVERGENCE` (integer), `NLL` (numeric), `OBJ` (TMB object), `PRIOR`,
#'   `DATA`, `WORKPAR`, `MODEL` and `ESTIMATOR` (character), `N_NODES` and
#'   `N_ITER` (integer).
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
  fit_mle(
    DATA,
    "trivariate",
    "glmm",
    THETA_START,
    PRIOR,
    CONTROL,
    match.arg(WORKPAR),
    N_NODES,
    N_ITER
  )
}
