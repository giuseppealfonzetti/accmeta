ib_jacobian_draw <- function(
  REP,
  SEEDS,
  THETA,
  DATA,
  PRIOR,
  WORKPAR,
  N_NODES,
  N_ITER
) {
  withr::local_seed(
    SEEDS[REP],
    .rng_kind = "Mersenne-Twister",
    .rng_normal_kind = "Inversion"
  )
  work <- function(t) if (WORKPAR == "Joe") theta2joe(t) else t
  theta_work <- work(THETA)

  for (attempt in seq_len(10L)) {
    d <- set_meta_data(
      sim_data(DATA$n_studies, THETA, DATA$margins[, "n"]),
      CC = DATA$CC
    )
    f <- try(
      fit_tlmm(d, THETA_START = THETA, PRIOR = PRIOR, WORKPAR = WORKPAR),
      silent = TRUE
    )
    if (!inherits(f, "try-error") && all(is.finite(f$THETA))) {
      break
    }
    f <- NULL
  }
  if (is.null(f)) {
    return(rep(NA, 18))
  }

  obj <- make_obj(d, THETA, set_prior(4, Inf), WORKPAR, "tglmm", N_NODES, N_ITER)
  c(work(f$THETA), -as.numeric(obj$gr(theta_work)))
}

ib_jacobian <- function(
  FIT,
  H = FIT$CONTROL$VCOV_H,
  N_NODES = FIT$CONTROL$VCOV_N_NODES,
  N_ITER = FIT$CONTROL$VCOV_N_ITER,
  NCORES = FIT$CONTROL$NCORES,
  SEEDS = NULL
) {
  stopifnot(inherits(FIT, "accmeta_ib"), H > 18)
  if (is.null(SEEDS)) {
    SEEDS <- sample.int(.Machine$integer.max, H)
  }
  args <- list(
    SEEDS = SEEDS,
    THETA = FIT$THETA,
    DATA = FIT$DATA,
    PRIOR = FIT$PRIOR,
    WORKPAR = FIT$WORKPAR,
    N_NODES = N_NODES,
    N_ITER = N_ITER
  )
  rows <- if (NCORES > 1L) {
    ib_start_daemons(NCORES, H)
    on.exit(mirai::daemons(0), add = TRUE)
    mirai::mirai_map(seq_len(H), ib_jacobian_draw, .args = args)[]
  } else {
    lapply(seq_len(H), function(REP) {
      do.call(ib_jacobian_draw, c(REP = REP, args))
    })
  }
  sim <- do.call(rbind, rows)
  sim <- sim[rowSums(!is.finite(sim)) == 0, , drop = FALSE]
  pi_work <- sim[, 1:9]
  score <- sim[, 10:18]
  list(
    J = stats::cov(pi_work, score),
    V = stats::cov(pi_work),
    PI = pi_work,
    SCORE = score,
    H_OK = nrow(sim),
    SEEDS = SEEDS,
    WORKPAR = FIT$WORKPAR
  )
}

ib_work_names <- function(WORKPAR) {
  if (WORKPAR == "Joe") {
    c(
      "mu_eta",
      "mu_xi",
      "mu_gamma",
      "log_sd_eta",
      "log_sd_xi",
      "log_sd_gamma",
      "z_rho12",
      "z_rho13",
      "z_rho23_1"
    )
  } else {
    c(
      "mu_eta",
      "mu_xi",
      "mu_gamma",
      "log_L11",
      "L21",
      "log_L22",
      "L31",
      "L32",
      "log_L33"
    )
  }
}
