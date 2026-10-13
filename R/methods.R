#' Print method for `accmeta_data` objects
#'
#' @param x An `accmeta_data` object.
#' @param ... Ignored.
#'
#' @rdname set_meta_data
#' @export
print.accmeta_data <- function(x, ...) {
  cat(
    "<accmeta_data>",
    x$n_studies,
    if (x$n_studies == 1) "study," else "studies,",
    "continuity correction",
    format(x$CC),
    "\n"
  )
  s <- cbind(x$rates, n = x$margins[, "n"])
  print(apply(s, 2, function(z) {
    c(min = min(z), median = stats::median(z), max = max(z))
  }))
  invisible(x)
}


#' Extract estimated parameters
#'
#' Extract random-effects mean vector and covariance matrix.
#'
#' @param object A fitted object, as returned by [fit_tlmm()], [fit_tglmm()] or
#'   [fit_ib()].
#' @param ... Ignored.
#'
#' @return A list with `MU` and `SIGMA`, named as documented for [theta2list()].
#'
#' @examples
#' th <- c(2.94, -2.2, -0.4, 0.0953, 0.4, -0.5108, 0.3, 0.2, -0.6931)
#' set.seed(1)
#' d <- set_meta_data(sim_data(15, th, rep(100, 15)), CC = 0.5)
#' fit <- fit_tlmm(d)
#' fit
#' coef(fit)
#'
#' @importFrom stats coef
#' @export
coef.accmeta_fit <- function(object, ...) {
  theta2list(object$THETA)
}


#' @param x A fitted object.
#'
#' @rdname coef.accmeta_fit
#' @export
print.accmeta_fit <- function(x, ...) {
  cat(
    paste0("<", class(x)[1], ">"),
    "convergence",
    x$CONVERGENCE,
    " NLL",
    format(x$NLL),
    "\n"
  )
  mu <- coef(x)$MU
  print(round(
    c(
      SE = stats::plogis(mu[[1]]),
      SP = 1 - stats::plogis(mu[[2]]),
      PREV = stats::plogis(mu[[3]])
    ),
    3
  ))
  invisible(x)
}


#' @rdname coef.accmeta_fit
#' @export
print.accmeta_ib <- function(x, ...) {
  cat(
    "<accmeta_ib>",
    x$N_ITER,
    if (x$N_ITER == 1) "iteration," else "iterations,",
    "stopped on",
    x$STOP,
    "\n"
  )
  mu <- coef(x)$MU
  print(round(
    c(
      SE = stats::plogis(mu[[1]]),
      SP = 1 - stats::plogis(mu[[2]]),
      PREV = stats::plogis(mu[[3]])
    ),
    3
  ))
  invisible(x)
}


#' Inference for IB fits
#'
#' Covariance matrix and confidence intervals for a [fit_ib()] estimate.
#'
#' @param object An `accmeta_ib` object, as returned by [fit_ib()].
#' @param VCOV If `TRUE`, computes and stores the Jacobian; if `FALSE`, drops
#'   it. The IB is not refitted.
#' @param METHOD `"wald"` or `"basic"`.
#' @param parm Parameters to report, by name or position. All by default.
#' @param level Confidence level.
#' @param ... Ignored by `vcov()` and `confint()`; not allowed by `update()`.
#'
#' @return `update()` returns the fit with its Jacobian added or removed.
#'   `vcov()` returns a named 9x9 matrix. `confint()` returns a matrix with one
#'   row per parameter and the lower and upper limits as columns.
#'
#' @examples
#' \dontrun{
#' th <- c(2.94, -2.2, -0.4, 0.0953, 0.4, -0.5108, 0.3, 0.2, -0.6931)
#' set.seed(1)
#' d <- set_meta_data(sim_data(15, th, sample(40:200, 15, TRUE)), CC = 0.5)
#' fit <- fit_ib(d, VCOV = TRUE, CONTROL = set_ib_control(NCORES = 4))
#' vcov(fit)
#' confint(fit)
#' confint(fit, METHOD = "basic")
#' }
#'
#' @importFrom stats vcov
#' @export
vcov.accmeta_ib <- function(object, ...) {
  jac <- object$JACOBIAN
  if (is.null(jac)) {
    stop(
      "No Jacobian in the fit: use fit_ib(VCOV = TRUE) or update(fit, VCOV = TRUE).",
      call. = FALSE
    )
  }
  j_inv <- solve(jac$J)
  out <- (1 + 1 / object$H[length(object$H)]) * j_inv %*% jac$V %*% t(j_inv)
  nm <- ib_work_names(object$WORKPAR)
  dimnames(out) <- list(nm, nm)
  out
}


#' @rdname vcov.accmeta_ib
#' @importFrom stats update
#' @export
update.accmeta_ib <- function(object, VCOV = TRUE, ...) {
  if (...length() > 0) {
    stop("update() on an accmeta_ib fit only takes VCOV.", call. = FALSE)
  }
  object$JACOBIAN <- if (VCOV) ib_jacobian(object)
  object
}


#' @rdname vcov.accmeta_ib
#' @importFrom stats confint
#' @export
confint.accmeta_ib <- function(
  object,
  parm,
  level = 0.95,
  METHOD = c("wald", "basic"),
  ...
) {
  METHOD <- match.arg(METHOD)
  v <- vcov(object)
  theta <- if (object$WORKPAR == "Joe") {
    theta2joe(object$THETA)
  } else {
    object$THETA
  }
  a <- (1 - level) / 2
  ci <- if (METHOD == "wald") {
    theta + outer(sqrt(diag(v)), stats::qnorm(c(a, 1 - a)))
  } else {
    jac <- object$JACOBIAN
    pi_dev <- sweep(jac$PI, 2, colMeans(jac$PI))
    theta_dev <- sqrt(1 + 1 / object$H[length(object$H)]) *
      t(solve(jac$J, t(pi_dev)))
    theta -
      t(apply(
        theta_dev,
        2,
        stats::quantile,
        probs = c(1 - a, a),
        names = FALSE
      ))
  }
  dimnames(ci) <- list(
    rownames(v),
    paste(format(100 * c(a, 1 - a), trim = TRUE), "%")
  )
  if (missing(parm)) ci else ci[parm, , drop = FALSE]
}
