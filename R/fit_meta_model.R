#' Fit a meta-analysis model of diagnostic accuracy
#'
#' @param DATA An `accmeta_data` object from [set_meta_data()], built with
#'   `CC > 0` for `"ib"` and `"lmm"` when cells are empty.
#' @param MODEL `"trivariate"` or `"bivariate"`.
#' @param ESTIMATOR `"ib"`, `"glmm"` or `"lmm"`.
#' @param PRIOR An `accmeta_prior` object from [set_prior()].
#' @param CONTROL An `accmeta_control` object from [set_control()], built for the
#'   same `ESTIMATOR`.
#'
#' @return Object of class `accmeta_ib` for `"ib"`, otherwise `accmeta_tlmm`,
#'   `accmeta_tglmm`, `accmeta_blmm` or `accmeta_bglmm`, with the components
#'   listed in [fit_ib()] and [fit_tglmm()].
#'
#' @examples
#' th <- c(2.94, -2.2, -0.4, 0.0953, 0.4, -0.5108, 0.3, 0.2, -0.6931)
#' set.seed(1)
#' design <- rep(c("cohort", "case-control"), 10)
#' x <- sim_data(20, th, rep(100, 20), N1 = ifelse(design == "cohort", NA, 40))
#' d <- set_meta_data(x, CC = 0.5, DESIGN = design)
#' fit <- fit_meta_model(d, MODEL = "trivariate", ESTIMATOR = "lmm")
#' confint(fit, SCALE = "reporting")
#' g <- fit_meta_model(
#'   set_meta_data(x, DESIGN = design),
#'   MODEL = "bivariate",
#'   ESTIMATOR = "glmm",
#'   CONTROL = set_control("glmm", N_NODES = 7L)
#' )
#' coef(g)
#'
#' @export
fit_meta_model <- function(
  DATA,
  MODEL = c("trivariate", "bivariate"),
  ESTIMATOR = c("ib", "glmm", "lmm"),
  PRIOR = set_prior(),
  CONTROL = set_control(ESTIMATOR)
) {
  MODEL <- match.arg(MODEL)
  ESTIMATOR <- match.arg(ESTIMATOR)
  stopifnot(inherits(CONTROL, "accmeta_control"))
  if (CONTROL$ESTIMATOR != ESTIMATOR) {
    stop(
      "CONTROL was built for ESTIMATOR = '",
      CONTROL$ESTIMATOR,
      "': use set_control('",
      ESTIMATOR,
      "').",
      call. = FALSE
    )
  }
  if (ESTIMATOR == "ib") {
    return(fit_ib(
      DATA,
      CONTROL = CONTROL$IB,
      PRIOR = PRIOR,
      WORKPAR = CONTROL$WORKPAR,
      SEEDS = CONTROL$SEEDS,
      VCOV = CONTROL$VCOV,
      MODEL = MODEL
    ))
  }
  fit_mle(
    DATA,
    MODEL,
    ESTIMATOR,
    THETA_START = NULL,
    PRIOR = PRIOR,
    CONTROL = CONTROL$OPTIM,
    WORKPAR = CONTROL$WORKPAR,
    N_NODES = CONTROL$N_NODES,
    N_ITER = CONTROL$N_ITER
  )
}


#' Set the estimator options
#'
#' @param ESTIMATOR `"ib"`, `"glmm"` or `"lmm"`. Setting an option of
#'   another estimator is an error.
#' @param WORKPAR `"Joe"` or `"PinheiroBates"`, covariance coordinates. All
#'   estimators.
#' @param MAXEVAL Positive integer, maximum optimiser evaluations. Default 1000.
#'   `"glmm"` and `"lmm"`.
#' @param N_NODES Integer, at least 2, quadrature nodes per dimension. Default 15.
#'   `"glmm"`.
#' @param N_ITER Positive integer, Newton steps for the integrand mode. Default
#'   10. `"glmm"`.
#' @param SEEDS Integer vector of length `H`, or `NULL` to draw them. `"ib"`.
#' @param VCOV Logical, compute the Jacobian needed by [vcov()] and [confint()].
#'   Default `FALSE`. `"ib"`.
#'
#' @examples
#' set_control("lmm")
#' set_control("glmm", N_NODES = 7L)
#' set_control("ib", H = 200, VCOV = TRUE)
#'
#' @inheritParams set_ib_control
#' @export
set_control <- function(
  ESTIMATOR = c("ib", "glmm", "lmm"),
  WORKPAR = c("Joe", "PinheiroBates"),
  MAXEVAL = NULL,
  N_NODES = NULL,
  N_ITER = NULL,
  SEEDS = NULL,
  VCOV = NULL,
  H = NULL,
  MAX_ITER = NULL,
  STEP = NULL,
  PATIENCE = NULL,
  TERMINATION = NULL,
  TOL = NULL,
  NCORES = NULL,
  VCOV_H = NULL,
  VCOV_N_NODES = NULL,
  VCOV_N_ITER = NULL,
  BOX_SD = NULL,
  BOX_COR = NULL,
  BOX_MU = NULL
) {
  ESTIMATOR <- match.arg(ESTIMATOR)
  WORKPAR <- match.arg(WORKPAR)
  ib_opts <- names(formals(set_ib_control))
  allowed <- switch(
    ESTIMATOR,
    lmm = "MAXEVAL",
    glmm = c("MAXEVAL", "N_NODES", "N_ITER"),
    ib = c("SEEDS", "VCOV", ib_opts)
  )
  opts <- c("MAXEVAL", "N_NODES", "N_ITER", "SEEDS", "VCOV", ib_opts)
  given <- Filter(Negate(is.null), mget(opts))
  bad <- setdiff(names(given), allowed)
  if (length(bad) > 0) {
    stop(
      paste(bad, collapse = ", "),
      " not used by ESTIMATOR = '",
      ESTIMATOR,
      "'.",
      call. = FALSE
    )
  }
  out <- list(ESTIMATOR = ESTIMATOR, WORKPAR = WORKPAR)
  if (ESTIMATOR == "ib") {
    stopifnot(is.null(VCOV) || isTRUE(VCOV) || isFALSE(VCOV))
    out$IB <- do.call(set_ib_control, given[intersect(names(given), ib_opts)])
    out$SEEDS <- SEEDS
    out$VCOV <- isTRUE(VCOV)
  } else {
    if (is.null(MAXEVAL)) {
      MAXEVAL <- 1000
    }
    stopifnot(is.numeric(MAXEVAL), length(MAXEVAL) == 1, MAXEVAL >= 1)
    out$OPTIM <- list(maxeval = MAXEVAL)
    out$N_NODES <- if (is.null(N_NODES)) 15L else N_NODES
    out$N_ITER <- if (is.null(N_ITER)) 10L else N_ITER
  }
  class(out) <- "accmeta_control"
  out
}
