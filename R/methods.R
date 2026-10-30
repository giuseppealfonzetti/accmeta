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
    c(min = min(z, na.rm = TRUE), median = stats::median(z, na.rm = TRUE), max = max(z, na.rm = TRUE))
  }))
  invisible(x)
}


#' Methods for fitted objects
#'
#' @param object A fitted object from [fit_meta_model()], [fit_tlmm()],
#'   [fit_tglmm()] or [fit_ib()].
#' @param ... Ignored.
#'
#' @return `coef()`: list with `MU` (named numeric, length 3 or 2) and `SIGMA`
#'   (named covariance matrix). `vcov()`: named numeric matrix, one row and column
#'   per parameter. `confint()`: numeric matrix with one row per parameter and two
#'   columns. `print()`: `x`, invisibly.
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
  rates <- stats::plogis(unname(c(mu[1], -mu[2], mu[3])))
  names(rates) <- c("SE", "SP", "PREV")
  print(round(rates[seq_along(mu)], 3))
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
  rates <- stats::plogis(unname(c(mu[1], -mu[2], mu[3])))
  names(rates) <- c("SE", "SP", "PREV")
  print(round(rates[seq_along(mu)], 3))
  invisible(x)
}


#' Inference for IB fits
#'
#' @param object An `accmeta_ib` object; `vcov()` and `confint()` need its
#'   Jacobian.
#' @param VCOV Logical, add (`TRUE`) or drop (`FALSE`) the Jacobian.
#' @param METHOD `"wald"` or `"basic"`.
#' @param SCALE `"working"`, `"natural"` or `"reporting"`; the last two need
#'   `WORKPAR = "Joe"`.
#' @param parm Character or integer vector, parameters to report. All by default.
#' @param level Numeric in (0, 1), confidence level.
#' @param ... Ignored by `vcov()` and `confint()`; not allowed by `update()`.
#'
#' @return `vcov()`: named numeric matrix, one row and column per parameter.
#'   `confint()`: numeric matrix with one row per parameter and two columns.
#'   `update()`: the `accmeta_ib` object.
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
  nm <- ib_work_names(object$WORKPAR, length(object$THETA))
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
  SCALE = c("working", "natural", "reporting"),
  ...
) {
  METHOD <- match.arg(METHOD)
  SCALE <- match.arg(SCALE)
  if (SCALE != "working" && object$WORKPAR != "Joe") {
    stop("SCALE = '", SCALE, "' needs a fit with WORKPAR = 'Joe'.", call. = FALSE)
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
    SCALE,
    dev
  )
  if (missing(parm)) ci else ci[parm, , drop = FALSE]
}


#' @param TYPE `"sandwich"` or `"hessian"`.
#' @param SCALE `"working"`, `"natural"` or `"reporting"`.
#' @param parm Character or integer vector, parameters to report. All by default.
#' @param level Numeric in (0, 1), confidence level.
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
  SCALE = c("working", "natural", "reporting"),
  ...
) {
  TYPE <- match.arg(TYPE)
  SCALE <- match.arg(SCALE)
  workpar <- if (SCALE == "working") object$WORKPAR else "Joe"
  ci <- ci_table(
    ib_to_work(object$THETA, workpar),
    fit_vcov(object, TYPE, workpar),
    level,
    SCALE
  )
  if (missing(parm)) ci else ci[parm, , drop = FALSE]
}


fit_vcov <- function(FIT, TYPE, WORKPAR) {
  model <- if (identical(FIT$ESTIMATOR, "glmm")) "tglmm" else "tlmm"
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
      "A prior with larger DEGREES (see set_prior()) keeps it inside.",
      call. = FALSE
    )
  }
  out <- solve(h)
  if (TYPE == "sandwich") {
    s <- sapply(seq_len(FIT$DATA$n_studies), function(i) {
      d <- set_meta_data(
        FIT$DATA$tab[i, , drop = FALSE],
        CC = FIT$DATA$CC,
        VERBOSE = FALSE,
        DESIGN = FIT$DATA$design[i]
      )
      build(d, flat_prior(FIT$THETA))$gr(par)
    })
    out <- out %*% tcrossprod(s) %*% out
  }
  nm <- ib_work_names(WORKPAR, length(FIT$THETA))
  dimnames(out) <- list(nm, nm)
  out
}


ci_table <- function(THETA, V, LEVEL, SCALE, DEV = NULL) {
  p <- length(THETA)
  m <- if (p == 5) 2 else 3
  g <- diag(p)
  nm <- rownames(V)
  if (SCALE != "working" && p == 9) {
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
  }
  if (SCALE != "working") {
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
    if (p == 5) {
      nm <- project(nm, "Joe")
    }
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
  if (SCALE != "working") {
    ci[m + seq_len(m), ] <- exp(ci[m + seq_len(m), ])
    ci[(2 * m + 1):p, ] <- tanh(ci[(2 * m + 1):p, ])
  }
  if (SCALE == "reporting") {
    ci[seq_len(m), ] <- stats::plogis(ci[seq_len(m), ])
    ci[2, ] <- 1 - ci[2, 2:1]
    nm[seq_len(m)] <- c("sensitivity", "specificity", "prevalence")[seq_len(m)]
  }
  dimnames(ci) <- list(
    nm,
    paste(format(100 * c(a, 1 - a), trim = TRUE), "%")
  )
  ci
}
