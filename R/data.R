#' Simulate misclassification tables
#'
#' Generates synthetic diagnostic accuracy data from the exact TGLMM model. For each
#' study a trivariate normal random effect \eqn{(\eta_i, \xi_i, \gamma_i)} is
#' drawn, followed by nested binomial draws: the diseased total from
#' \eqn{\gamma_i}, the true positives from \eqn{\eta_i}, and the false positives
#' from \eqn{\xi_i}.
#'
#' @param N_STUDIES Number of studies to generate. A single positive integer.
#' @param THETA Numeric vector of length 9. See [theta2list()].
#' @param N_I Numeric vector of length `N_STUDIES` giving the total sample size
#'   \eqn{n_i} of each study. Entries must be positive integers.
#'
#' @return An integer matrix with `N_STUDIES` rows and columns `n11`, `n10`,
#'   `n01`, `n00`, the misclassification counts of each study. Each row sums to
#'   the corresponding entry of `N_I`.
#'
#' @examples
#' th <- c(2.94, -2.2, -0.4, 0.0953, 0.4, -0.5108, 0.3, 0.2, -0.6931)
#' set.seed(1)
#' sim_data(5, th, rep(100, 5))
#' sim_data(5, th, sample(40:200, 5, TRUE))
#'
#' @export
sim_data <- function(N_STUDIES, THETA, N_I) {
  stopifnot(
    is.numeric(N_STUDIES),
    length(N_STUDIES) == 1,
    N_STUDIES >= 1,
    is.numeric(N_I),
    length(N_I) == N_STUDIES,
    all(N_I > 0),
    all(N_I == round(N_I))
  )
  par <- theta2list(THETA)
  U <- matrix(stats::rnorm(3 * N_STUDIES), N_STUDIES, 3) %*%
    chol(par$SIGMA) +
    rep(par$MU, each = N_STUDIES)
  n1 <- stats::rbinom(N_STUDIES, N_I, stats::plogis(U[, 3]))
  n11 <- stats::rbinom(N_STUDIES, n1, stats::plogis(U[, 1]))
  n10 <- stats::rbinom(N_STUDIES, N_I - n1, stats::plogis(U[, 2]))
  out <- cbind(n11 = n11, n10 = n10, n01 = n1 - n11, n00 = N_I - n1 - n10)
  return(out)
}


#' Validate misclassification counts and pre-compute study summaries
#'
#' Checks the four misclassification counts of each study and pre-computes study
#' summaries used by [fit_tlmm()] or [fit_tglmm()].
#'
#' @param MATRIX An integer matrix or data frame with four columns holding the
#'   misclassification counts of each study. If the
#'   columns are named, they must be named `n11`, `n10`, `n01`, `n00`.
#'   If they are unnamed, that order is assumed.
#' @param CC Non-negative continuity correction, added to each of the four cells
#'   of every study. Needed by [fit_tlmm()] when zero cells are observed.
#' @param VERBOSE allow messages.
#'
#' @return An object of class `accmeta_data`, a list with components:
#'   \describe{
#'     \item{`tab`}{the raw counts, columns `n11`, `n10`, `n01`, `n00`.}
#'     \item{`margins`}{the raw margins, columns `n1`, `n0`, `n`.}
#'     \item{`counts`}{raw binomial responses and denominators for the exact
#'       model, columns `y_eta`, `y_xi`, `y_gamma`, `den_eta`, `den_xi`,
#'       `den_gamma`.}
#'     \item{`est`}{logit of `SE`, of `1 - SP` and of `PREV`, columns `eta`, `xi`, `gamma`.}
#'     \item{`wvar`}{within-study variances, columns `var_eta`, `var_xi`,
#'       `var_gamma`.}
#'     \item{`rates`}{per-study accuracy, columns `SE`, `SP`, `PREV`.}
#'     \item{`n_studies`}{the number of studies.}
#'     \item{`CC`}{the continuity correction applied.}
#'   }
#'
#' @seealso [sim_data()] for the input layout, [fit_tlmm()] and [fit_tglmm()]
#'   for the fitting functions that consume the result.
#'
#' @examples
#' th <- c(2.94, -2.2, -0.4, 0.0953, 0.4, -0.5108, 0.3, 0.2, -0.6931)
#' set.seed(4)
#' x <- sim_data(5, th, rep(100, 5))
#' d <- set_meta_data(x)
#' d
#' d$rates
#' # a correction shifts the logit-scale estimates towards the centre
#' cbind(d$est[, "eta"], set_meta_data(x, CC = 0.5)$est[, "eta"])
#'
#' @export
set_meta_data <- function(MATRIX, CC = 0, VERBOSE = TRUE) {
  stopifnot(
    is.matrix(MATRIX) || is.data.frame(MATRIX),
    is.numeric(CC),
    length(CC) == 1,
    is.finite(CC),
    CC >= 0
  )
  if (is.data.frame(MATRIX)) {
    MATRIX <- as.matrix(MATRIX)
  }
  cells <- c("n11", "n10", "n01", "n00")
  stopifnot(
    is.numeric(MATRIX),
    nrow(MATRIX) >= 1,
    ncol(MATRIX) == 4
  )
  if (is.null(colnames(MATRIX))) {
    message("unnamed columns: assuming the order n11, n10, n01, n00")
    colnames(MATRIX) <- cells
  } else if (!setequal(colnames(MATRIX), cells)) {
    stop("columns of MATRIX must be named n11, n10, n01, n00")
  }

  # subset raw data matrix
  tab <- MATRIX[, cells, drop = FALSE]
  stopifnot(
    all(is.finite(tab)),
    all(tab >= 0),
    all(tab == round(tab))
  )

  # margins
  n1 <- tab[, "n11"] + tab[, "n01"]
  n0 <- tab[, "n10"] + tab[, "n00"]
  margins <- cbind(n1 = n1, n0 = n0, n = n1 + n0)
  if (any(margins[, "n"] == 0)) {
    stop("every study must have at least one observation")
  }

  # data for tglmm
  counts <- cbind(
    y_eta = tab[, "n11"],
    y_xi = tab[, "n10"],
    y_gamma = n1,
    den_eta = n1,
    den_xi = n0,
    den_gamma = n1 + n0
  )

  # data after cc
  m <- tab + CC
  c1 <- m[, "n11"] + m[, "n01"]
  c0 <- m[, "n10"] + m[, "n00"]
  y <- cbind(m[, "n11"], m[, "n10"], c1)
  den <- cbind(c1, c0, c1 + c0)
  empty <- which(rowSums(y == 0 | y == den) > 0)
  if (length(empty) > 0 & VERBOSE) {
    lab <- if (is.null(rownames(tab))) empty else rownames(tab)[empty]
    message(
      "empty cell in stud",
      if (length(empty) == 1) "y " else "ies ",
      paste(lab, collapse = ", "),
      ": 'est' and 'wvar' are infinite. Rebuild with CC > 0 if needed."
    )
  }

  p <- y / den

  # data for tlmm
  est <- cbind(
    eta = stats::qlogis(p[, 1]),
    xi = stats::qlogis(p[, 2]),
    gamma = stats::qlogis(p[, 3])
  )
  wvar <- cbind(
    var_eta = 1 / y[, 1] + 1 / (den[, 1] - y[, 1]),
    var_xi = 1 / y[, 2] + 1 / (den[, 2] - y[, 2]),
    var_gamma = 1 / y[, 3] + 1 / (den[, 3] - y[, 3])
  )
  rates <- cbind(SE = p[, 1], SP = 1 - p[, 2], PREV = p[, 3])

  # return output
  out <- list(
    tab = tab,
    margins = margins,
    counts = counts,
    est = est,
    wvar = wvar,
    rates = rates,
    n_studies = nrow(tab),
    CC = CC
  )
  out[1:6] <- lapply(out[1:6], function(x) {
    rownames(x) <- rownames(MATRIX)
    x
  })
  class(out) <- "accmeta_data"
  return(out)
}
