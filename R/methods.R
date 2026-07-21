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
