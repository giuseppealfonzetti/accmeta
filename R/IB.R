#' Iterative bootstrap bias correction for TLMM
#'
#' @param DATA An `accmeta_data` object, as returned by [set_meta_data()],
#'   with `CC > 0`.
#' @param CONTROL Iterative-bootstrap control object from [set_ib_control()].
#' @param PRIOR Prior on \eqn{\Sigma_3}, as returned by [set_prior()].
#' @param WORKPAR Working scale for the matching equation and the
#'   update: `"PinheiroBates"` based on the log-Cholesky decomposition,
#'   or `"Joe"` (default) based on unconstrained partial correlations. Also
#'   passed to the auxiliary [fit_tlmm()] fits.
#' @param SEEDS Integer vector of length `H` seeding the simulated datasets. If
#'   `NULL`, drawn once and then held fixed, so serial and parallel runs agree.
#'   When `NULL`, the seeds are drawn from the caller's RNG, which is the only
#'   change `fit_ib()` makes to the caller's RNG state.
#' @param VCOV If `TRUE`, also computes the Jacobian at the estimate, stored as `JACOBIAN`, which
#'   [vcov()][vcov.accmeta_ib] and [confint()][confint.accmeta_ib] need. Its settings are the
#'   `VCOV_*` entries of [set_ib_control()]. It can be added later with
#'   [update()][update.accmeta_ib].
#'
#' @return An `accmeta_ib` object. `STOP` says why the iteration ended: `"converged"`
#'   (the convergence test passed), `"stall"` (no improvement for `PATIENCE` iterations),
#'   `"boundary"` (the next update would leave the box of [set_ib_control()], in which case
#'   `THETA` is the last iterate inside it) or `"maxit"`, in which case `THETA` is the last
#'   iterate. The TLMM start is clamped into the box. `RATIO` is the ratio of successive steps,
#'   an empirical contraction rate.
#'
#' @seealso [fit_tlmm()] for the auxiliary estimator, [set_ib_control()] for the
#'   control settings, and [set_prior()] for the prior specification.
#'
#' @examples
#' th <- c(2.94, -2.2, -0.4, 0.0953, 0.4, -0.5108, 0.3, 0.2, -0.6931)
#' set.seed(1)
#' x <- sim_data(15, th, sample(40:200, 15, TRUE))
#' fit <- fit_ib(set_meta_data(x, CC = 0.5), CONTROL = set_ib_control(H = 20, MAX_ITER = 3))
#' rbind(TLMM = fit$PI_HAT, IB = fit$THETA)[, 1:3]
#'
#' @export
fit_ib <- function(
  DATA,
  CONTROL = set_ib_control(),
  PRIOR = set_prior(),
  WORKPAR = c("Joe", "PinheiroBates"),
  SEEDS = NULL,
  VCOV = FALSE
) {
  # Setup
  n_params <- 9
  # run the whole budget: no convergence or stall stops
  budget_only <- CONTROL$TERMINATION == "max_iter"

  WORKPAR <- match.arg(WORKPAR)
  stopifnot(
    inherits(DATA, "accmeta_data"),
    is.matrix(DATA$tab),
    inherits(PRIOR, "accmeta_prior"),
    inherits(CONTROL, "accmeta_ib_control"),
    is.null(SEEDS) || (is.numeric(SEEDS) && length(SEEDS) == CONTROL$H),
    isTRUE(VCOV) || isFALSE(VCOV)
  )

  # manage seeds for reproducibility between serial and parallel exec
  if (is.null(SEEDS)) {
    SEEDS <- sample.int(.Machine$integer.max, CONTROL$H)
  }

  # setup parallal cluster via mirai package
  use_parallel <- CONTROL$NCORES > 1L
  if (use_parallel) {
    ib_start_daemons(CONTROL$NCORES, CONTROL$H)
    on.exit(mirai::daemons(0), add = TRUE)
  }
  if (is.null(CONTROL$MAX_ITER)) {
    CONTROL$MAX_ITER <- if (DATA$n_studies < 15) 3 else 100
    if (DATA$n_studies < 15) {
      warning(
        "fewer than 15 studies: the IB is stopped after 3 steps, where further ",
        "iterations tend to overcorrect; set MAX_ITER to override.",
        call. = FALSE
      )
    }
  }
  if (DATA$CC <= 0) {
    stop(
      "DATA must carry a continuity correction. See `?set_meta_data()`",
      call. = FALSE
    )
  }

  n_i <- DATA$margins[, "n"]

  # initial estimate
  pi_hat <- fit_tlmm(DATA, PRIOR = PRIOR, WORKPAR = WORKPAR)$THETA
  pi_hat_work <- ib_to_work(pi_hat, WORKPAR)
  box <- ib_box(WORKPAR, CONTROL)
  theta <- ib_from_work(pmin(pmax(pi_hat_work, box$L), box$U), WORKPAR)

  # initialise path tracking
  path <- matrix(NA, CONTROL$MAX_ITER + 1, n_params)
  path[1, ] <- theta
  fail <- integer(CONTROL$MAX_ITER)
  degen <- numeric(CONTROL$MAX_ITER)
  ratio <- rep(NA, CONTROL$MAX_ITER)
  gap_path <- matrix(NA, CONTROL$MAX_ITER, n_params)
  se_path <- matrix(NA, CONTROL$MAX_ITER, n_params)
  progress <- rep(NA, CONTROL$MAX_ITER)
  threshold_path <- rep(NA, CONTROL$MAX_ITER)

  # auxiliary quantities
  best_theta <- theta
  best_val <- Inf
  best_t2 <- Inf
  best_iter <- 0
  converged <- FALSE
  stop_rule <- "maxit"
  filled <- 1

  step_last <- NULL

  # root finding loop
  for (iter in seq_len(CONTROL$MAX_ITER)) {
    # compute estimator correction
    gap_summary <- ib_gap(
      THETA = theta,
      SEEDS = SEEDS,
      H = CONTROL$H,
      N_STUDIES = DATA$n_studies,
      N_I = n_i,
      CC = DATA$CC,
      PRIOR = PRIOR,
      USE_PARALLEL = use_parallel,
      PI_HAT_WORK = pi_hat_work,
      WORKPAR = WORKPAR,
      MIN_EIG = 1e-3
    )
    if (is.null(gap_summary)) {
      stop("all simulated TLMM fits failed at iteration ", iter, call. = FALSE)
    }

    h_ok <- gap_summary$H_OK

    # track gap-related quantities
    gap_path[iter, ] <- gap_summary$GAP
    se_path[iter, ] <- gap_summary$SE
    fail[iter] <- gap_summary$FAIL
    degen[iter] <- gap_summary$DEGEN

    # check enough simulations are ok
    if (h_ok <= n_params) {
      stop(
        "only ",
        h_ok,
        " of ",
        CONTROL$H,
        " simulated fits ",
        "succeeded at iteration ",
        iter,
        "; the Mahalanobis test needs more ",
        "than ",
        n_params,
        ".",
        call. = FALSE
      )
    }

    # update and residual
    res <- gap_summary$GAP
    step_k <- CONTROL$STEP * res
    if (!is.null(step_last)) {
      ratio[iter] <- sqrt(sum(step_k^2)) / sqrt(sum(step_last^2))
    }
    step_last <- step_k

    # add small diagonal constant for stability
    live <- diag(gap_summary$V) > 1e-10 * max(diag(gap_summary$V))
    vcov_reg <- gap_summary$V[live, live, drop = FALSE] +
      diag(1e-8 * pmax(diag(gap_summary$V)[live], 1e-12), sum(live))

    # mahalanobis distance of the residual for the convergence test
    t2 <- drop(crossprod(res[live], solve(vcov_reg, res[live])))
    threshold <- if (budget_only) NA else (h_ok - 1) * CONTROL$PRECISION / h_ok

    r2 <- t2 / (h_ok - 1)

    # track t2-related quantities
    progress[iter] <- t2
    threshold_path[iter] <- threshold
    if (r2 < best_val) {
      best_val <- r2
      best_t2 <- t2
      best_theta <- theta
      best_iter <- iter
    }

    # stop by convergence test
    if (!budget_only && t2 < threshold) {
      best_theta <- theta
      best_t2 <- t2
      converged <- TRUE
      stop_rule <- "converged"
      break
    }

    # stop before an update leaving the box
    work <- ib_to_work(theta, WORKPAR) + step_k
    if (any(work < box$L | work > box$U)) {
      best_theta <- theta
      best_t2 <- t2
      stop_rule <- "boundary"
      break
    }

    # update estimates
    theta <- ib_from_work(work, WORKPAR)
    path[iter + 1, ] <- theta
    filled <- iter + 1L

    # stop by patience on stall updateds
    if (!budget_only && iter - best_iter >= CONTROL$PATIENCE) {
      stop_rule <- "stall"
      break
    }
  }

  if (stop_rule == "maxit" && !budget_only) {
    best_theta <- path[filled, ]
  }

  # report explicitely failed fits
  if (sum(fail[seq_len(iter)]) > 0) {
    warning(
      sum(fail[seq_len(iter)]),
      " simulated TLMM fits failed and were excluded."
    )
  }

  out <- list(
    THETA = best_theta,
    PI_HAT = pi_hat,
    N_ITER = iter,
    CONVERGED = converged,
    STOP = stop_rule,
    RESIDUAL = best_t2,
    PROGRESS = progress[seq_len(iter)],
    THRESHOLD = threshold_path[seq_len(iter)],
    PATH = path[seq_len(filled), , drop = FALSE],
    FAIL = fail[seq_len(iter)],
    DEGEN = degen[seq_len(iter)],
    RATIO = ratio[seq_len(iter)],
    GAP = gap_path[seq_len(iter), , drop = FALSE],
    SE = se_path[seq_len(iter), , drop = FALSE],
    SEEDS = SEEDS,
    PRIOR = PRIOR,
    CONTROL = CONTROL,
    WORKPAR = WORKPAR,
    DATA = DATA
  )
  class(out) <- c("accmeta_ib", "accmeta_fit")
  if (VCOV) {
    if (use_parallel) {
      mirai::daemons(0)
    }
    out$JACOBIAN <- ib_jacobian(out)
  }
  return(out)
}


ib_start_daemons <- function(NCORES, N) {
  if (!requireNamespace("mirai", quietly = TRUE)) {
    stop("NCORES > 1 needs the 'mirai' package.", call. = FALSE)
  }
  if (mirai::daemons_set()) {
    stop(
      "Close your existing pool with mirai::daemons(0) before proceeding",
      call. = FALSE
    )
  }
  mirai::daemons(min(NCORES, N))
}

ib_box <- function(WORKPAR, CONTROL) {
  mu <- rep(CONTROL$BOX_MU, 3)
  sd <- log(CONTROL$BOX_SD)
  if (WORKPAR == "Joe") {
    z <- rep(atanh(CONTROL$BOX_COR), 3)
    return(list(L = c(-mu, rep(sd[1], 3), -z), U = c(mu, rep(sd[2], 3), z)))
  }
  off <- CONTROL$BOX_SD[2]
  list(
    L = c(-mu, sd[1], -off, sd[1], -off, -off, sd[1]),
    U = c(mu, sd[2], off, sd[2], off, off, sd[2])
  )
}

ib_to_work <- function(THETA, WORKPAR) {
  if (WORKPAR == "Joe") theta2joe(THETA) else THETA
}

ib_from_work <- function(WORK, WORKPAR) {
  if (WORKPAR == "Joe") joe2theta(WORK) else WORK
}

# single fit helper function
ib_one_fit <- function(REP, SEEDS, THETA, N_STUDIES, N_I, CC, PRIOR, WORKPAR) {
  withr::local_seed(
    SEEDS[REP],
    .rng_kind = "Mersenne-Twister",
    .rng_normal_kind = "Inversion"
  )

  # attempts loop to defend from bad sims
  for (attempt in seq_len(10L)) {
    d <- accmeta::set_meta_data(
      accmeta::sim_data(N_STUDIES, THETA, N_I),
      CC = CC
    )
    f <- try(
      accmeta::fit_tlmm(
        d,
        THETA_START = THETA,
        PRIOR = PRIOR,
        WORKPAR = WORKPAR
      ),
      silent = TRUE
    )
    if (!inherits(f, "try-error") && all(is.finite(f$THETA))) {
      return(f$THETA)
    }
  }
  rep(NA, 9)
}

# helper function to compute the ib correction term at a given iteration
ib_gap <- function(
  THETA,
  SEEDS,
  H,
  N_STUDIES,
  N_I,
  CC,
  PRIOR,
  USE_PARALLEL,
  PI_HAT_WORK,
  WORKPAR,
  MIN_EIG
) {
  rows <- if (USE_PARALLEL) {
    mirai::mirai_map(
      seq_len(H),
      ib_one_fit,
      .args = list(
        SEEDS = SEEDS,
        THETA = THETA,
        N_STUDIES = N_STUDIES,
        N_I = N_I,
        CC = CC,
        PRIOR = PRIOR,
        WORKPAR = WORKPAR
      )
    )[]
  } else {
    lapply(seq_len(H), function(REP) {
      ib_one_fit(
        REP = REP,
        SEEDS = SEEDS,
        THETA = THETA,
        N_STUDIES = N_STUDIES,
        N_I = N_I,
        CC = CC,
        PRIOR = PRIOR,
        WORKPAR = WORKPAR
      )
    })
  }
  sim <- do.call(rbind, rows)
  ok <- stats::complete.cases(sim)
  if (!any(ok)) {
    return(NULL)
  }
  sim_ok <- sim[ok, , drop = FALSE]
  degen <- mean(apply(sim_ok, 1, function(t) {
    min(
      eigen(theta2list(t)$SIGMA, symmetric = TRUE, only.values = TRUE)$values
    ) <
      MIN_EIG
  }))
  sim_ok <- t(apply(sim_ok, 1, ib_to_work, WORKPAR = WORKPAR))
  h_ok <- nrow(sim_ok)
  list(
    GAP = PI_HAT_WORK - colMeans(sim_ok),
    SE = pmax(apply(sim_ok, 2, stats::sd) / sqrt(h_ok), 1e-5),
    V = stats::cov(sim_ok) / h_ok,
    H_OK = h_ok,
    FAIL = sum(!ok),
    DEGEN = degen
  )
}
