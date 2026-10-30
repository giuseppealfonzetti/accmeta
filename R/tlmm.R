#' Fit the trivariate linear mixed model
#'
#' @param DATA An `accmeta_data` object from [set_meta_data()] with finite `est`.
#' @param THETA_START Numeric vector of length 9, or `NULL` for [init_theta()].
#' @param PRIOR An `accmeta_prior` object from [set_prior()].
#' @param CONTROL List passed to [ucminf::ucminf()].
#' @param WORKPAR `"PinheiroBates"` or `"Joe"`, covariance coordinates.
#'
#' @return Object of class `accmeta_tlmm`, a list with `THETA` (numeric, length
#'   9), `CONVERGENCE` (integer), `NLL` (numeric), `OBJ` (TMB object), `PRIOR`,
#'   `DATA`, `WORKPAR`, `MODEL` and `ESTIMATOR` (character).
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
  fit_mle(
    DATA,
    "trivariate",
    "lmm",
    THETA_START,
    PRIOR,
    CONTROL,
    match.arg(WORKPAR)
  )
}

fit_mle <- function(
  DATA,
  MODEL,
  ESTIMATOR,
  THETA_START = NULL,
  PRIOR = set_prior(),
  CONTROL = list(maxeval = 1000),
  WORKPAR = "Joe",
  N_NODES = 15L,
  N_ITER = 10L
) {
  biv <- MODEL == "bivariate"
  p <- if (biv) 5 else 9
  k <- if (biv) 1:2 else 1:3
  stopifnot(
    inherits(DATA, "accmeta_data"),
    is.matrix(DATA$counts),
    is.matrix(DATA$est)
  )
  if (!biv && !any(DATA$counts[, "den_gamma"] > 0)) {
    stop(
      "the trivariate model needs at least one cohort study; ",
      "use MODEL = 'bivariate'.",
      call. = FALSE
    )
  }
  if (ESTIMATOR == "lmm") {
    stopifnot(all(is.finite(DATA$est[, k]) | is.na(DATA$est[, k])))
  }
  PRIOR <- resolve_prior(PRIOR, length(k))
  if (is.null(THETA_START)) {
    THETA_START <- init_theta(DATA, MODEL = MODEL)
  }
  stopifnot(
    is.numeric(THETA_START),
    length(THETA_START) == p,
    is.numeric(N_NODES),
    length(N_NODES) == 1,
    N_NODES >= 2,
    is.numeric(N_ITER),
    length(N_ITER) == 1,
    N_ITER >= 1
  )
  obj <- make_obj(
    DATA,
    THETA_START,
    PRIOR,
    WORKPAR,
    paste0("t", ESTIMATOR),
    N_NODES,
    N_ITER
  )
  est <- ucminf::ucminf(
    par = obj$par,
    fn = obj$fn,
    gr = obj$gr,
    control = CONTROL
  )
  out <- list(
    THETA = unname(c(
      est$par[seq_along(k)],
      obj$report(est$par)$LOGCHOL[seq_len(p - length(k))]
    )),
    CONVERGENCE = est$convergence,
    NLL = est$value,
    OBJ = obj,
    PRIOR = PRIOR,
    DATA = DATA,
    WORKPAR = WORKPAR,
    MODEL = MODEL,
    ESTIMATOR = ESTIMATOR
  )
  if (ESTIMATOR == "glmm") {
    out$N_NODES <- N_NODES
    out$N_ITER <- N_ITER
  }
  class(out) <- c(
    paste0("accmeta_", if (biv) "b" else "t", ESTIMATOR),
    "accmeta_fit"
  )
  return(out)
}
