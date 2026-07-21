#' Unpack the working parameter vector
#'
#' Maps the working parameter vector `THETA` onto its interpretable
#' components: the random-effects mean vector and covariance matrix.
#'
#' @param THETA Numeric vector of length 9. Entries 1 to 3 are the random-effects
#'   means \eqn{(\bar\eta, \bar\xi, \bar\gamma)}. Entries 4 to 9 are the
#'   log-Cholesky factor of \eqn{\Sigma_3} in row-major lower-triangular order
#'   \eqn{(log(L_{11}), L_{21}, log(L_{22}), L_{31}, L_{32}, log(L_{33})}).
#'
#' @return A list with two components: `MU`, the length-3 mean vector, and
#'   `SIGMA`, the 3x3 positive-definite covariance matrix.
#'
#' @seealso [list2theta()] for the inverse map.
#'
#' @examples
#' theta2list(c(2.94, -2.2, -0.4, 0.0953, 0.4, -0.5108, 0.3, 0.2, -0.6931))
#'
#' @export
theta2list <- function(THETA) {
  stopifnot(is.numeric(THETA), length(THETA) == 9)
  L <- matrix(0, 3, 3)
  L[lower.tri(L, diag = TRUE)] <- THETA[c(4, 5, 7, 6, 8, 9)]
  diag(L) <- exp(diag(L))
  list(MU = THETA[1:3], SIGMA = tcrossprod(L))
}


#' Pack the working parameter vector
#'
#' Inverse of [theta2list()]: packs the random-effects mean vector and covariance
#' matrix into the working parameter vector `THETA`.
#'
#' @param LIST A list with components `MU`, a numeric vector of length 3, and
#'   `SIGMA`, a 3x3 symmetric positive-definite matrix.
#'
#' @return A numeric vector of length 9, in the layout documented for
#'   [theta2list()].
#'
#' @seealso [theta2list()] for the inverse map.
#'
#' @examples
#' li <- list(
#'   MU = c(2.94, -2.2, -0.4),
#'   SIGMA = matrix(c(1.21, 0.44, 0.33,
#'                    0.44, 0.52, 0.24,
#'                    0.33, 0.24, 0.38), 3, 3)
#' )
#' list2theta(li)
#' all.equal(theta2list(list2theta(li)), li)
#'
#' @export
list2theta <- function(LIST) {
  stopifnot(
    is.numeric(LIST$MU),
    length(LIST$MU) == 3,
    identical(dim(LIST$SIGMA), c(3L, 3L)),
    isSymmetric(LIST$SIGMA)
  )
  L <- t(chol(LIST$SIGMA))
  diag(L) <- log(diag(L))
  c(LIST$MU, L[lower.tri(L, diag = TRUE)][c(1, 2, 4, 3, 5, 6)])
}

#' Set the prior on the random-effects covariance
#'
#' Control the prior object for \eqn{\Sigma_3 \sim W(\nu, A I_3)} to be passed to [fit_tlmm()] or [fit_tglmm()] through the `PRIOR` argument. The prior is
#' \eqn{\Sigma_3 \sim W(\nu, A I_3)}. Setting \eqn{\nu=4} and \eqn{A=Inf} correspond to maximum likelihood estimation.
#'
#' @param DEGREES Degrees of freedom \eqn{\nu} of the Wishart prior.
#' @param SCALE Scale \eqn{A}. Represents a soft ceiling on the variance scale.
#'
#' @return An object of class `accmeta_prior`: a list with `DEGREES` and
#'   `SCALE`, held as given and passed to the template unchanged.
#'
#' @examples
#' set_prior()
#' set_prior(DEGREES = 5, SCALE = 100)
#'
#' @export
set_prior <- function(DEGREES = 4, SCALE = Inf) {
  stopifnot(
    is.numeric(DEGREES),
    length(DEGREES) == 1,
    DEGREES >= 4,
    DEGREES == round(DEGREES),
    is.numeric(SCALE),
    length(SCALE) == 1,
    SCALE > 0
  )
  out <- list(DEGREES = DEGREES, SCALE = SCALE)
  class(out) <- "accmeta_prior"
  return(out)
}
