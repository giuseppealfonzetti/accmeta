#' Iterative bootstrap bias correction for TLMM
#'
#' @param DATA An `accmeta_data` object, as returned by [set_meta_data()],
#'   with `CC > 0`.
#' @param CONTROL Iterative-bootstrap control object from [set_ib_control()].
#' @param PRIOR Prior on \eqn{\Sigma_3}, as returned by [set_prior()].
#' @param WORKPAR Working scale for the matching equation and the
#'   update: `"PinheiroBates"` based on the log-Cholesky decomposition,
#'   or `"Joe"` (default) based on unconstrained partial correlations.
#' @param SEEDS Integer vector of length `H` seeding the simulated datasets. If
#'   `NULL`, drawn once and then held fixed. With `BOOST`, the seeds for up to
#'   `MAX_H` datasets are appended once, before any fit, so serial and parallel
#'   runs agree; the returned `SEEDS` are the ones used.
#'
#' @return An `accmeta_ib` object. `STOP` says why the iteration ended: `"tol"`
#'   (the convergence test passed), `"boundary"` (no valid generative model available for next iteration),
#'   `"stall"` (no improvement for `PATIENCE` iterations) or `"maxit"`.
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
  SEEDS = NULL
) {
  # Setup
  min_eig <- 1e-3
  n_params <- 9
  # run the whole budget: no convergence, stall or boundary stops
  budget_only <- CONTROL$TERMINATION == "max_iter"

  WORKPAR <- match.arg(WORKPAR)
  stopifnot(
    inherits(DATA, "accmeta_data"),
    is.matrix(DATA$tab),
    inherits(PRIOR, "accmeta_prior"),
    inherits(CONTROL, "accmeta_ib_control"),
    is.null(SEEDS) || (is.numeric(SEEDS) && length(SEEDS) == CONTROL$H)
  )

  # manage seeds for reproducibility between serial and parallel exec
  if (is.null(SEEDS)) {
    SEEDS <- sample.int(.Machine$integer.max, CONTROL$H)
  }

  if (CONTROL$BOOST) {
    SEEDS <- c(
      SEEDS,
      sample.int(.Machine$integer.max, CONTROL$MAX_H - CONTROL$H)
    )
  }

  # setup parallal cluster via mirai package
  use_parallel <- CONTROL$NCORES > 1L
  if (use_parallel) {
    if (!requireNamespace("mirai", quietly = TRUE)) {
      stop("NCORES > 1 needs the 'mirai' package.", call. = FALSE)
    }
    if (mirai::status()$connections > 0L) {
      stop(
        "Tear down your existing pool with mirai::daemons(0).",
        call. = FALSE
      )
    }
    mirai::daemons(min(CONTROL$NCORES, CONTROL$H))
    on.exit(mirai::daemons(0), add = TRUE)
  }
  if (DATA$CC <= 0) {
    stop(
      "DATA must carry a continuity correction. See `?set_meta_data()`",
      call. = FALSE
    )
  }

  n_i <- DATA$margins[, "n"]

  # initial estimate
  pi_hat <- fit_tlmm(DATA, PRIOR = PRIOR)$THETA
  pi_hat_work <- if (WORKPAR == "Joe") theta2joe(pi_hat) else pi_hat
  theta <- project_pd(pi_hat, min_eig)

  # initialise path tracking
  path <- matrix(NA, CONTROL$MAX_ITER + 1, n_params)
  path[1, ] <- theta
  fail <- integer(CONTROL$MAX_ITER)
  degen <- numeric(CONTROL$MAX_ITER)
  halved <- integer(CONTROL$MAX_ITER)
  stepsize <- rep(NA, CONTROL$MAX_ITER)
  gap_path <- matrix(NA, CONTROL$MAX_ITER, n_params)
  se_path <- matrix(NA, CONTROL$MAX_ITER, n_params)
  progress <- rep(NA, CONTROL$MAX_ITER)
  threshold_path <- rep(NA, CONTROL$MAX_ITER)
  h_path <- rep(NA, CONTROL$MAX_ITER)

  # auxiliary quantities
  best_theta <- theta
  best_val <- Inf
  best_t2 <- Inf
  best_iter <- 0
  converged <- FALSE
  stop_rule <- "maxit"
  filled <- 1

  gap_summary <- NULL

  # sample size boosting (confidence termination), as in ergm
  h <- CONTROL$H
  gap_last <- NULL
  not_improved <- rep(FALSE, 4)

  # root finding loop
  for (iter in seq_len(CONTROL$MAX_ITER)) {
    if (is.null(gap_summary)) {
      # compute estimator correction
      gap_summary <- ib_gap(
        THETA = theta,
        SEEDS = SEEDS,
        H = h,
        N_STUDIES = DATA$n_studies,
        N_I = n_i,
        CC = DATA$CC,
        PRIOR = PRIOR,
        USE_PARALLEL = use_parallel,
        PI_HAT_WORK = pi_hat_work,
        WORKPAR = WORKPAR,
        MIN_EIG = min_eig
      )
    }
    if (is.null(gap_summary)) {
      stop("all simulated TLMM fits failed at iteration ", iter, call. = FALSE)
    }

    h_ok <- gap_summary$H_OK

    # track gap-related quantities
    h_path[iter] <- h
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
        h,
        " simulated fits ",
        "succeeded at iteration ",
        iter,
        "; the Hotelling test needs more ",
        "than ",
        n_params,
        ".",
        call. = FALSE
      )
    }

    # add small diagonal constant for stability
    vcov_reg <- gap_summary$V +
      diag(1e-8 * pmax(diag(gap_summary$V), 1e-12), n_params)

    # compute hotelling t2 for convergence test
    t2 <- drop(crossprod(
      gap_summary$GAP,
      solve(vcov_reg, gap_summary$GAP)
    ))

    # threshold on t2. Inspired by {ergm} stopping criteria for MCMLE
    crit <- ib_crit(h_ok, CONTROL$TOL, n_params)
    threshold <- switch(
      CONTROL$TERMINATION,
      hotelling = crit,
      confidence = (h_ok - 1) *
        max(sqrt(CONTROL$PRECISION) - sqrt(crit / (h_ok - 1)), 0)^2,
      residual = (h_ok - 1) * CONTROL$PRECISION / h_ok,
      max_iter = NA
    )

    r2 <- t2 / (h_ok - 1)
    if (!is.null(gap_last)) {
      r2_last <- drop(crossprod(gap_last, solve(vcov_reg, gap_last))) /
        (h_ok - 1)
      not_improved <- c(not_improved[-1], r2 >= r2_last)
    }
    gap_last <- gap_summary$GAP

    # track t2-related quantities
    progress[iter] <- t2
    threshold_path[iter] <- threshold
    if (r2 < best_val) {
      best_val <- r2
      best_t2 <- t2
      best_theta <- theta
      best_iter <- iter
    }

    # update estimates
    theta_work <- if (WORKPAR == "Joe") theta2joe(theta) else theta
    step <- ib_step_fixedpoint(
      THETA_WORK = theta_work,
      GAP = gap_summary$GAP,
      STEP = CONTROL$STEP,
      WORKPAR = WORKPAR,
      MIN_EIG = min_eig
    )

    # store update-realted quantities
    halved[iter] <- step$HALVED
    stepsize[iter] <- step$STEPSIZE
    gap_summary <- NULL

    # no positive-definite step left
    if (is.null(step$THETA)) {
      stop_rule <- "boundary"
      break
    }

    # store thate
    theta <- step$THETA
    path[iter + 1, ] <- theta
    filled <- iter + 1L

    # stop by convergence test
    if (!budget_only && t2 < threshold) {
      best_theta <- path[iter, ]
      best_t2 <- t2
      converged <- TRUE
      stop_rule <- "tol"
      break
    }

    # root on the boundary
    if (
      !budget_only &&
        iter > 1 &&
        halved[iter] > 0 &&
        halved[iter - 1] > 0 &&
        best_iter < iter
    ) {
      stop_rule <- "boundary"
      break
    }

    # increase H by BOOST_FACTOR when:
    # 1) T2 inside tolerance region but its confidence region is not;
    # 2) when updates stall with T2 outside the tolerance region
    if (CONTROL$BOOST && h < CONTROL$MAX_H) {
      inside <- r2 < CONTROL$PRECISION
      boost <- if (inside) {
        min(
          crit /
            ((h_ok - 1) * (sqrt(CONTROL$PRECISION) - sqrt(r2))^2),
          CONTROL$BOOST_FACTOR
        )
      } else if (sum(not_improved) > 1) {
        not_improved[] <- FALSE
        CONTROL$BOOST_FACTOR
      } else {
        1
      }
      if (boost > 1) {
        h_new <- min(ceiling(h * boost), CONTROL$MAX_H)
        h <- h_new
        gap_summary <- NULL
        if (inside) {
          best_iter <- iter
        }
      }
    }

    # stop by patience on stall updateds
    if (!budget_only && iter - best_iter >= CONTROL$PATIENCE) {
      stop_rule <- "stall"
      break
    }
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
    HALVED = halved[seq_len(iter)],
    STEPSIZE = stepsize[seq_len(iter)],
    GAP = gap_path[seq_len(iter), , drop = FALSE],
    SE = se_path[seq_len(iter), , drop = FALSE],
    H = h_path[seq_len(iter)],
    SEEDS = SEEDS[seq_len(h)],
    PRIOR = PRIOR,
    CONTROL = CONTROL
  )
  class(out) <- c("accmeta_ib", "accmeta_fit")
  return(out)
}


# critical value of the Hotelling t2 at level TOL with H draws
ib_crit <- function(H, TOL, N_PARAMS = 9) {
  N_PARAMS *
    (H - 1) /
    (H - N_PARAMS) *
    stats::qf(1 - TOL, N_PARAMS, H - N_PARAMS)
}

# ensure pd reff sigma
project_pd <- function(THETA, MIN_EIG) {
  li <- theta2list(THETA)
  eig <- eigen(li$SIGMA, symmetric = TRUE)
  if (min(eig$values) >= MIN_EIG) {
    return(THETA)
  }
  sigma <- eig$vectors %*%
    diag(pmax(eig$values, MIN_EIG), 3, 3) %*%
    t(eig$vectors)
  list2theta(list(MU = li$MU, SIGMA = (sigma + t(sigma)) / 2))
}

# single fit helper function
ib_one_fit <- function(REP, SEEDS, THETA, N_STUDIES, N_I, CC, PRIOR) {
  set.seed(SEEDS[REP], kind = "Mersenne-Twister", normal.kind = "Inversion")

  # attempts loop to defend from bad sims
  for (attempt in seq_len(10L)) {
    d <- accmeta::set_meta_data(
      accmeta::sim_data(N_STUDIES, THETA, N_I),
      CC = CC
    )
    f <- try(
      accmeta::fit_tlmm(d, THETA_START = THETA, PRIOR = PRIOR),
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
        PRIOR = PRIOR
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
        PRIOR = PRIOR
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
  if (WORKPAR == "Joe") {
    sim_ok <- t(apply(sim_ok, 1, theta2joe))
  }
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

# update with eventual stepsize halvening
ib_step_fixedpoint <- function(THETA_WORK, GAP, STEP, WORKPAR, MIN_EIG) {
  step <- STEP
  halved <- 0L
  cand <- NULL
  repeat {
    proposal <- THETA_WORK + step * GAP
    trial <- try(
      if (WORKPAR == "Joe") joe2theta(proposal) else proposal,
      silent = TRUE
    )
    valid <- !inherits(trial, "try-error") && all(is.finite(trial))
    if (valid) {
      sigma <- theta2list(trial)$SIGMA
      valid <- all(is.finite(sigma)) &&
        min(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values) >=
          MIN_EIG &&
        !inherits(try(chol(sigma), silent = TRUE), "try-error")
    }
    if (valid) {
      cand <- trial
      break
    }
    step <- step / 2
    halved <- halved + 1L
    if (step < STEP * 2^-19) {
      break
    }
  }
  list(THETA = cand, HALVED = halved, STEPSIZE = step)
}
