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
  U <- matrix(rnorm(3 * N_STUDIES), N_STUDIES, 3) %*%
    chol(par$SIGMA) +
    rep(par$MU, each = N_STUDIES)
  n1 <- rbinom(N_STUDIES, N_I, plogis(U[, 3]))
  n11 <- rbinom(N_STUDIES, n1, plogis(U[, 1]))
  n10 <- rbinom(N_STUDIES, N_I - n1, plogis(U[, 2]))
  cbind(n11 = n11, n10 = n10, n01 = n1 - n11, n00 = N_I - n1 - n10)
}
