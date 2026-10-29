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
#' @param SCALE `"working"` for intervals in the coordinates of `WORKPAR`, or
#'   `"natural"` for the means, standard deviations and correlations.
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
  out <- (1 + 1 / object$CONTROL$H) * j_inv %*% jac$V %*% t(j_inv)
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
  SCALE = c("working", "natural"),
  ...
) {
  METHOD <- match.arg(METHOD)
  SCALE <- match.arg(SCALE)
  if (SCALE == "natural" && object$WORKPAR != "Joe") {
    stop("SCALE = 'natural' needs a fit with WORKPAR = 'Joe'.", call. = FALSE)
  }
  v <- vcov(object)
  dev <- if (METHOD == "basic") {
    jac <- object$JACOBIAN
    pi_dev <- sweep(jac$PI, 2, colMeans(jac$PI))
    sqrt(1 + 1 / object$CONTROL$H) * t(solve(jac$J, t(pi_dev)))
  }
  ci <- ci_table(
    ib_to_work(object$THETA, object$WORKPAR),
    v,
    level,
    SCALE == "natural",
    dev
  )
  if (missing(parm)) ci else ci[parm, , drop = FALSE]
}


#' @param TYPE `"sandwich"` or `"hessian"` .
#' @param SCALE `"working"` for intervals in the coordinates of the fit's
#'   `WORKPAR`, or `"natural"` for the means, standard deviations and
#'   correlations.
#' @param parm Parameters to report, by name or position. All by default.
#' @param level Confidence level.
#'
#' @rdname coef.accmeta_fit
#' @importFrom stats vcov
#' @export
vcov.accmeta_fit <- function(object, TYPE = c("sandwich", "hessian"), ...) {
  fit_vcov(object, match.arg(TYPE), object$WORKPAR)
}


#' @rdname coef.accmeta_fit
#' @importFrom stats confint
#' @export
confint.accmeta_fit <- function(
  object,
  parm,
  level = 0.95,
  TYPE = c("sandwich", "hessian"),
  SCALE = c("working", "natural"),
  ...
) {
  TYPE <- match.arg(TYPE)
  SCALE <- match.arg(SCALE)
  workpar <- if (SCALE == "natural") "Joe" else object$WORKPAR
  ci <- ci_table(
    ib_to_work(object$THETA, workpar),
    fit_vcov(object, TYPE, workpar),
    level,
    SCALE == "natural"
  )
  if (missing(parm)) ci else ci[parm, , drop = FALSE]
}


fit_vcov <- function(FIT, TYPE, WORKPAR) {
  model <- if (inherits(FIT, "accmeta_tlmm")) "tlmm" else "tglmm"
  build <- function(DATA, PRIOR) {
    make_obj(
      DATA,
      FIT$THETA,
      PRIOR,
      WORKPAR,
      model,
      FIT$N_NODES,
      FIT$N_ITER
    )
  }
  par <- ib_to_work(FIT$THETA, WORKPAR)
  h <- build(FIT$DATA, FIT$PRIOR)$he(par)
  if (rcond(h) < .Machine$double.eps) {
    stop(
      "singular Hessian: the fit is on the boundary of the parameter space. ",
      "A prior with DEGREES > 4 (see set_prior()) keeps it inside.",
      call. = FALSE
    )
  }
  out <- solve(h)
  if (TYPE == "sandwich") {
    s <- sapply(seq_len(FIT$DATA$n_studies), function(i) {
      d <- set_meta_data(
        FIT$DATA$tab[i, , drop = FALSE],
        CC = FIT$DATA$CC,
        VERBOSE = FALSE
      )
      build(d, set_prior(4, Inf))$gr(par)
    })
    out <- out %*% tcrossprod(s) %*% out
  }
  nm <- ib_work_names(WORKPAR)
  dimnames(out) <- list(nm, nm)
  out
}


ci_table <- function(THETA, V, LEVEL, NATURAL, DEV = NULL) {
  g <- diag(9)
  nm <- rownames(V)
  if (NATURAL) {
    r <- tanh(THETA[7:9])
    q <- sqrt((1 - r[1]^2) * (1 - r[2]^2))
    rho <- r[1] * r[2] + r[3] * q
    g[9, ] <- c(
      rep(0, 6),
      r[2] * (1 - r[1]^2) - r[1] * r[3] * q,
      r[1] * (1 - r[2]^2) - r[2] * r[3] * q,
      (1 - r[3]^2) * q
    ) /
      (1 - rho^2)
    THETA[9] <- atanh(rho)
    nm <- c(
      "mu_eta",
      "mu_xi",
      "mu_gamma",
      "sd_eta",
      "sd_xi",
      "sd_gamma",
      "rho_eta_xi",
      "rho_eta_gamma",
      "rho_xi_gamma"
    )
  }
  a <- (1 - LEVEL) / 2
  ci <- if (is.null(DEV)) {
    THETA + outer(sqrt(diag(g %*% V %*% t(g))), stats::qnorm(c(a, 1 - a)))
  } else {
    THETA -
      t(apply(
        DEV %*% t(g),
        2,
        stats::quantile,
        probs = c(1 - a, a),
        names = FALSE
      ))
  }
  if (NATURAL) {
    ci[4:6, ] <- exp(ci[4:6, ])
    ci[7:9, ] <- tanh(ci[7:9, ])
  }
  dimnames(ci) <- list(
    nm,
    paste(format(100 * c(a, 1 - a), trim = TRUE), "%")
  )
  ci
}
