#' Fit by iterative bootstrap
#'
#' @param DATA An `accmeta_data` object from [set_meta_data()] with `CC > 0`.
#' @param CONTROL An `accmeta_ib_control` object from [set_ib_control()].
#' @param PRIOR An `accmeta_prior` object from [set_prior()].
#' @param WORKPAR `"Joe"` or `"PinheiroBates"`, covariance coordinates.
#' @param SEEDS Integer vector of length `CONTROL$H`, or `NULL` to draw them from
#'   the caller's random number generator.
#' @param VCOV Logical, also compute the Jacobian, stored as `JACOBIAN`.
#' @param MODEL `"trivariate"` or `"bivariate"`.
#'
#' @return Object of class `accmeta_ib`, a list with `THETA` (numeric, length 9
#'   or 5), `PI_HAT` (numeric, auxiliary estimate), `N_ITER` (integer),
#'   `CONVERGED` (logical), `STOP` (`"converged"`, `"stall"`, `"boundary"`
#'   or `"maxit"`), `RESIDUAL` (numeric), `PROGRESS`, `THRESHOLD`, `FAIL`,
#'   `DEGEN`, `RATIO` (numeric vectors, one entry per iteration), `PATH`, `GAP`,
#'   `SE` (numeric matrices, one row per iteration), `SEEDS`, `PRIOR`, `CONTROL`,
#'   `WORKPAR`, `DATA`, `MODEL`, `ESTIMATOR`, and `JACOBIAN` (list) when
#'   `VCOV = TRUE`.
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
  VCOV = FALSE,
  MODEL = c("trivariate", "bivariate")
) {
  MODEL <- match.arg(MODEL)
  # run the whole budget: no convergence or stall stops
  budget_only <- CONTROL$TERMINATION == "max_iter"

  WORKPAR <- match.arg(WORKPAR)
  stopifnot(
    inherits(DATA, "accmeta_data"),
    is.matrix(DATA$tab),
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
    CONTROL$MAX_ITER <- 100
  }
  if (DATA$CC <= 0) {
    stop(
      "DATA must carry a continuity correction. See `?set_meta_data()`",
      call. = FALSE
    )
  }

  # initial estimate
  PRIOR <- resolve_prior(PRIOR, if (MODEL == "bivariate") 2 else 3)
  pi_hat <- fit_mle(
    DATA,
    MODEL,
    "lmm",
    NULL,
    PRIOR,
    list(maxeval = 1000),
    WORKPAR
  )$THETA
  n_params <- length(pi_hat)
  pi_hat_work <- ib_to_work(pi_hat, WORKPAR)
  box <- ib_box(WORKPAR, CONTROL, n_params)
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
      DATA = DATA,
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
    delta_m <- if (iter > 1) sqrt(t2 / progress[iter - 1]) else NA
    threshold <- if (budget_only) {
      NA
    } else {
      (h_ok - 1) * CONTROL$TOL / h_ok * (1 - delta_m)^2
    }

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
    if (!budget_only && isTRUE(delta_m < 1 && t2 < threshold)) {
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
    DATA = DATA,
    MODEL = MODEL,
    ESTIMATOR = "ib"
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

ib_box <- function(WORKPAR, CONTROL, P = 9) {
  if (P == 5) {
    return(lapply(ib_box(WORKPAR, CONTROL), project, WORKPAR = WORKPAR))
  }
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

sim_like <- function(DATA, THETA) {
  m <- DATA$margins
  design <- if (is.null(DATA$design)) rep("cohort", DATA$n_studies) else DATA$design
  n1 <- if (length(THETA) == 5) {
    m[, "n1"]
  } else {
    ifelse(design == "case-control", m[, "n1"], NA)
  }
  set_meta_data(
    sim_data(DATA$n_studies, THETA, m[, "n"], N1 = n1),
    CC = DATA$CC,
    DESIGN = design
  )
}

# single fit helper function
ib_one_fit <- function(REP, SEEDS, THETA, DATA, PRIOR, WORKPAR) {
  withr::local_seed(
    SEEDS[REP],
    .rng_kind = "Mersenne-Twister",
    .rng_normal_kind = "Inversion"
  )

  # attempts loop to defend from bad sims
  for (attempt in seq_len(10L)) {
    d <- sim_like(DATA, THETA)
    f <- try(
      fit_mle(
        d,
        theta_model(THETA),
        "lmm",
        THETA,
        PRIOR,
        list(maxeval = 1000),
        WORKPAR
      ),
      silent = TRUE
    )
    if (!inherits(f, "try-error") && all(is.finite(f$THETA))) {
      return(f$THETA)
    }
  }
  rep(NA, length(THETA))
}

# helper function to compute the ib correction term at a given iteration
ib_gap <- function(
  THETA,
  SEEDS,
  H,
  DATA,
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
        DATA = DATA,
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
        DATA = DATA,
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
