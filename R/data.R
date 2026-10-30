#' Simulate misclassification tables
#'
#' @param N_STUDIES Positive integer, number of studies.
#' @param THETA Numeric vector of length 9 (trivariate), 5 (bivariate) or 12
#'   (skew-normal, needs 'sn'), in the layout of [theta2list()].
#' @param N_I Positive integer vector of length `N_STUDIES`, study sizes.
#' @param DF Numeric greater than 2, degrees of freedom of Student t random
#'   effects; `Inf` for normal. Length-9 `THETA` only.
#' @param N1 Optional integer vector of length `N_STUDIES`, number of diseased of
#'   case-control studies and `NA` for cohort studies. Required, without `NA`, for
#'   a length-5 `THETA`.
#'
#' @return Integer matrix with `N_STUDIES` rows and columns `n11`, `n10`, `n01`,
#'   `n00`.
#'
#' @examples
#' th <- c(2.94, -2.2, -0.4, 0.0953, 0.4, -0.5108, 0.3, 0.2, -0.6931)
#' set.seed(1)
#' sim_data(5, th, rep(100, 5))
#' sim_data(5, th, sample(40:200, 5, TRUE))
#'
#' @export
sim_data <- function(N_STUDIES, THETA, N_I, DF = Inf, N1 = NULL) {
  stopifnot(
    is.numeric(N_STUDIES),
    length(N_STUDIES) == 1,
    N_STUDIES >= 1,
    is.numeric(N_I),
    length(N_I) == N_STUDIES,
    all(N_I > 0),
    all(N_I == round(N_I)),
    is.numeric(DF),
    length(DF) == 1,
    DF > 2,
    is.null(N1) || (length(N1) == N_STUDIES && all(is.na(N1) | (N1 >= 0 & N1 <= N_I)))
  )
  if (length(THETA) == 5) {
    if (is.null(N1) || anyNA(N1)) {
      stop("a bivariate THETA needs N1 for every study.", call. = FALSE)
    }
    THETA <- embed(THETA, "PinheiroBates")
  }
  if (is.finite(DF) && length(THETA) == 12) {
    stop("t random effects need a length-9 THETA.", call. = FALSE)
  }
  par <- theta2list(THETA)
  if (length(THETA) == 9) {
    U <- matrix(stats::rnorm(3 * N_STUDIES), N_STUDIES, 3) %*% chol(par$SIGMA)
    if (is.finite(DF)) {
      U <- U * sqrt((DF - 2) / stats::rchisq(N_STUDIES, DF))
    }
    U <- U + rep(par$MU, each = N_STUDIES)
  } else {
    if (!requireNamespace("sn", quietly = TRUE)) {
      stop("skew-normal random effects need the 'sn' package.", call. = FALSE)
    }
    dp <- theta2dp(THETA)
    dp <- list(xi = c(0, 0, 0), Omega = dp$Omega, alpha = dp$alpha)
    U <- sn::rmsn(N_STUDIES, dp = dp) +
      rep(par$MU - sn::dp2cp(dp, "SN", upto = 1)$mean, each = N_STUDIES)
  }
  n1 <- stats::qbinom(stats::runif(N_STUDIES), N_I, stats::plogis(U[, 3]))
  if (!is.null(N1)) {
    n1[!is.na(N1)] <- N1[!is.na(N1)]
  }
  n11 <- stats::qbinom(stats::runif(N_STUDIES), n1, stats::plogis(U[, 1]))
  n10 <- stats::qbinom(stats::runif(N_STUDIES), N_I - n1, stats::plogis(U[, 2]))
  out <- cbind(n11 = n11, n10 = n10, n01 = n1 - n11, n00 = N_I - n1 - n10)
  return(out)
}


#' Build the data object
#'
#' @param MATRIX Matrix or data frame of non-negative integer counts with four
#'   columns named `n11`, `n10`, `n01`, `n00` (taken in this order if unnamed).
#' @param CC Non-negative numeric, continuity correction added to every cell.
#' @param VERBOSE Logical, print messages.
#' @param DESIGN Optional character vector with one entry per study, `"cohort"`
#'   or `"case-control"`. Default all `"cohort"`.
#'
#' @return Object of class `accmeta_data`, a list with:
#'   \describe{
#'     \item{`tab`}{numeric matrix of counts, columns `n11`, `n10`, `n01`, `n00`.}
#'     \item{`margins`}{numeric matrix, columns `n1`, `n0`, `n`.}
#'     \item{`counts`}{numeric matrix, columns `y_eta`, `y_xi`, `y_gamma`,
#'       `den_eta`, `den_xi`, `den_gamma`.}
#'     \item{`est`}{numeric matrix, columns `eta`, `xi`, `gamma`.}
#'     \item{`wvar`}{numeric matrix, columns `var_eta`, `var_xi`, `var_gamma`.}
#'     \item{`rates`}{numeric matrix, columns `SE`, `SP`, `PREV`.}
#'     \item{`n_studies`}{integer.}
#'     \item{`CC`}{numeric.}
#'     \item{`design`}{character vector.}
#'   }
#'   The `gamma`, `var_gamma` and `PREV` entries are `NA` for case-control studies.
#'
#' @seealso [sim_data()] for the input layout, [fit_tlmm()] and [fit_tglmm()]
#'   for the fitting functions that consume the result.
#'
#' @examples
#' th <- c(2.94, -2.2, -0.4, 0.0953, 0.4, -0.5108, 0.3, 0.2, -0.6931)
#' set.seed(2)
#' x <- sim_data(5, th, rep(100, 5))
#' d <- set_meta_data(x)
#' d
#' d$rates
#' # a correction shifts the logit-scale estimates towards the centre
#' cbind(d$est[, "eta"], set_meta_data(x, CC = 0.5)$est[, "eta"])
#'
#' @export
set_meta_data <- function(MATRIX, CC = 0, VERBOSE = TRUE, DESIGN = NULL) {
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
  if (is.null(DESIGN)) {
    DESIGN <- rep("cohort", nrow(tab))
  }
  stopifnot(
    length(DESIGN) == nrow(tab),
    all(DESIGN %in% c("cohort", "case-control"))
  )
  cc <- DESIGN == "case-control"

  # data for tglmm
  counts <- cbind(
    y_eta = tab[, "n11"],
    y_xi = tab[, "n10"],
    y_gamma = ifelse(cc, 0, n1),
    den_eta = n1,
    den_xi = n0,
    den_gamma = ifelse(cc, 0, n1 + n0)
  )

  # data after cc
  m <- tab + CC
  c1 <- m[, "n11"] + m[, "n01"]
  c0 <- m[, "n10"] + m[, "n00"]
  y <- cbind(m[, "n11"], m[, "n10"], c1)
  den <- cbind(c1, c0, c1 + c0)
  bad <- y == 0 | y == den
  bad[cc, 3] <- FALSE
  empty <- which(rowSums(bad) > 0)
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
  est[cc, "gamma"] <- NA
  wvar[cc, "var_gamma"] <- NA
  rates <- cbind(SE = p[, 1], SP = 1 - p[, 2], PREV = p[, 3])
  rates[cc, "PREV"] <- NA

  # return output
  out <- list(
    tab = tab,
    margins = margins,
    counts = counts,
    est = est,
    wvar = wvar,
    rates = rates,
    n_studies = nrow(tab),
    CC = CC,
    design = DESIGN
  )
  out[1:6] <- lapply(out[1:6], function(x) {
    rownames(x) <- rownames(MATRIX)
    x
  })
  class(out) <- "accmeta_data"
  return(out)
}
