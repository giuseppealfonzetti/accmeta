#' Iterative bootstrap bias correction for TLMM
#'
#' Estimates from TLMM estimates might contain bias from up to three distinct sources:
#' 1) within-study normal approximation; 2) continuity correction; 3) Wishart prior.
#' Iterative bootstrap removes them all.
#'
#' @param DATA An `accmeta_data` object, as returned by [set_meta_data()],
#'   with `CC > 0`.
#' @param H Number of datasets simulated per iteration.
#' @param MAX_ITER Maximum number of iterations.
#' @param TOL Convergence tolerance on
#'   \eqn{\max_j |\hat\pi_j - \overline\pi_{H,j}(\theta)|}, the departure from the
#'   matching equation that defines the estimator.
#' @param STEP Damping factor \eqn{\gamma \in (0, 1]}.
#' @param PRIOR Prior on \eqn{\Sigma_3}, as returned by [set_prior()].
#' @param SEEDS Integer vector of length `H` seeding the simulated datasets. If
#'   `NULL`, drawn once and then held fixed.
#'
#' @details
#' The `H` seeds are drawn once and reused for every \eqn{\theta} and at every
#' iteration. The state of the random number generator is restored on exit.
#'
#' The `PRIOR` and the continuity correction recorded on `DATA` reach \eqn{\hat\pi}
#' and every simulated fit.
#'
#' Differently from [fit_tlmm()] and [fit_tglmm()], a penalised `PRIOR` is the default here.
#'
#' @return A list with components `THETA`, the bias-corrected estimate;
#'   `PI_HAT`, the uncorrected TLMM estimate; `N_ITER`, the number of iterations
#'   used; `CONVERGED`, whether the tolerance was met; `STOP`, which rule ended
#'   the recursion, one of `"tol"`, `"plateau"` or
#'   `"maxit"`; `RESIDUAL`, the largest absolute departure from the matching
#'   equation at the last iteration; `PROGRESS`, that departure at every
#'   iteration; `PATH`, the iterates, one row per step; `FAIL`, the number of
#'   simulated fits that failed every retry, per iteration; `DEGEN`, the
#'   proportion of the `H` simulated fits whose \eqn{\Sigma_3} was degenerate,
#'   per iteration; and `SEEDS`.
#'
#' @seealso [fit_tlmm()] for the auxiliary estimator and [set_prior()] for the
#'   prior specification.
#'
#' @examples
#' th <- c(2.94, -2.2, -0.4, 0.0953, 0.4, -0.5108, 0.3, 0.2, -0.6931)
#' set.seed(1)
#' x <- sim_data(15, th, sample(40:200, 15, TRUE))
#' fit <- fit_ib(set_meta_data(x, CC = 0.5), H = 20, MAX_ITER = 3)
#' rbind(TLMM = fit$PI_HAT, IB = fit$THETA)[, 1:3]
#'
#' @export
fit_ib <- function(
  DATA,
  H = 100,
  MAX_ITER = 25,
  TOL = 0.2 / sqrt(H),
  STEP = 1,
  PRIOR = set_prior(DEGREES = 5),
  SEEDS = NULL
) {
  stopifnot(
    inherits(DATA, "accmeta_data"),
    is.matrix(DATA$tab),
    inherits(PRIOR, "accmeta_prior"),
    is.numeric(H),
    length(H) == 1,
    H >= 2,
    is.numeric(MAX_ITER),
    length(MAX_ITER) == 1,
    MAX_ITER >= 1,
    is.numeric(TOL),
    length(TOL) == 1,
    TOL > 0,
    is.numeric(STEP),
    length(STEP) == 1,
    STEP > 0,
    STEP <= 1,
    is.null(SEEDS) || (is.numeric(SEEDS) && length(SEEDS) == H)
  )
  if (DATA$CC <= 0) {
    stop(
      "DATA must carry a continuity correction, since every simulated dataset ",
      "is transformed the same way: set_meta_data(MATRIX, CC = 0.5)",
      call. = FALSE
    )
  }
  if (exists(".Random.seed", .GlobalEnv)) {
    old <- get(".Random.seed", .GlobalEnv)
    on.exit(assign(".Random.seed", old, .GlobalEnv), add = TRUE)
  }
  if (is.null(SEEDS)) {
    SEEDS <- sample.int(.Machine$integer.max, H)
  }

  n_studies <- DATA$n_studies
  n_i <- DATA$margins[, "n"]
  pi_hat <- fit_tlmm(DATA, PRIOR = PRIOR)$THETA

  theta <- pi_hat
  path <- matrix(NA, MAX_ITER + 1, 9)
  path[1, ] <- theta
  fail <- integer(MAX_ITER)
  degen <- numeric(MAX_ITER)
  progress <- rep(NA, MAX_ITER)
  converged <- FALSE
  residual <- NA
  stop_rule <- "maxit"
  k <- 0L

  for (i in seq_len(MAX_ITER)) {
    sim <- matrix(NA, H, 9)
    for (h in seq_len(H)) {
      set.seed(SEEDS[h])
      # retry on failure
      for (attempt in seq_len(10L)) {
        d <- set_meta_data(sim_data(n_studies, theta, n_i), CC = DATA$CC)
        f <- try(
          fit_tlmm(d, THETA_START = theta, PRIOR = PRIOR),
          silent = TRUE
        )
        if (!inherits(f, "try-error") && all(is.finite(f$THETA))) {
          sim[h, ] <- f$THETA
          break
        }
      }
    }
    ok <- stats::complete.cases(sim)
    fail[i] <- sum(!ok)
    if (!any(ok)) {
      stop("all simulated TLMM fits failed at iteration ", i, call. = FALSE)
    }
    degen[i] <- mean(apply(sim[ok, , drop = FALSE], 1, function(t) {
      min(
        eigen(
          theta2list(t)$SIGMA,
          symmetric = TRUE,
          only.values = TRUE
        )$values
      ) <
        1e-8
    }))

    # tolerance rule on update
    gap <- pi_hat - colMeans(sim[ok, , drop = FALSE])
    theta <- theta + STEP * gap
    path[i + 1, ] <- theta
    k <- i
    residual <- max(abs(gap))
    progress[i] <- residual
    if (residual < TOL) {
      converged <- TRUE
      stop_rule <- "tol"
      break
    }

    # flat progress curve, after ib
    if (k > 10L) {
      y <- progress[k:(k - 10L)]
      x <- k:(k - 10L)
      if (summary(stats::lm(y ~ x))$coefficients[2, 4] > 0.2) {
        stop_rule <- "plateau"
        break
      }
    }
  }

  if (sum(fail[seq_len(k)]) > 0) {
    warning(
      sum(fail[seq_len(k)]),
      " simulated TLMM fits failed and were excluded; ",
      "the average is then over a selected subset"
    )
  }

  list(
    THETA = path[k + 1, ],
    PI_HAT = pi_hat,
    N_ITER = k,
    CONVERGED = converged,
    STOP = stop_rule,
    RESIDUAL = residual,
    PROGRESS = progress[seq_len(k)],
    PATH = path[seq_len(k + 1), , drop = FALSE],
    FAIL = fail[seq_len(k)],
    DEGEN = degen[seq_len(k)],
    SEEDS = SEEDS
  )
}
