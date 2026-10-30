#' Diagnostic plots for IB fits
#'
#' @param x An `accmeta_ib` object.
#' @param WHICH Character vector, any of `"convergence"`, `"contraction"` and
#'   `"path"`. All by default.
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#'
#' @seealso [fit_ib()]
#'
#' @examples
#' th <- c(2.94, -2.2, -0.4, 0.0953, 0.4, -0.5108, 0.3, 0.2, -0.6931)
#' set.seed(1)
#' d <- set_meta_data(sim_data(15, th, sample(40:200, 15, TRUE)), CC = 0.5)
#' fit <- fit_ib(d, CONTROL = set_ib_control(H = 20), SEEDS = 1:20)
#' plot(fit)
#'
#' @importFrom graphics plot
#' @export
plot.accmeta_ib <- function(
  x,
  WHICH = c("convergence", "contraction", "path"),
  ...
) {
  tinyplot::tinytheme("clean2")
  WHICH <- match.arg(WHICH, several.ok = TRUE)
  k <- seq_along(x$PROGRESS) - 1
  if ("convergence" %in% WHICH) {
    conv <- data.frame(
      iterate = c(k, k),
      value = c(x$PROGRESS, x$THRESHOLD),
      series = rep(c("step statistic", "threshold"), each = length(k))
    )
    tinyplot::tinyplot(
      value ~ iterate | series,
      data = conv[is.finite(conv$value) & conv$value > 0, ],
      type = "o",
      lty = "by",
      log = "y",
      ylab = "Mahalanobis step",
      main = paste0("Stop: ", x$STOP, " after ", x$N_ITER, " iterations")
    )
  }
  if ("contraction" %in% WHICH && length(k) > 1) {
    contr <- data.frame(
      iterate = k[-1],
      rate = sqrt(x$PROGRESS[-1] / x$PROGRESS[-length(k)])
    )
    tinyplot::tinyplot(
      rate ~ iterate,
      data = contr,
      type = "o",
      ylim = range(c(contr$rate, 1), finite = TRUE),
      ylab = "Estimated contraction rate",
      main = "Contraction"
    )
    graphics::abline(h = 1, lty = 2)
  }
  if ("path" %in% WHICH) {
    p <- length(x$THETA)
    nm <- ib_work_names(x$WORKPAR, p)
    work <- t(apply(x$PATH, 1, ib_to_work, WORKPAR = x$WORKPAR))
    path <- data.frame(
      iterate = rep(seq_len(nrow(work)) - 1, times = p),
      value = c(work),
      param = factor(rep(nm, each = nrow(work)), nm),
      series = "iterate"
    )
    if (x$STOP == "boundary") {
      box <- ib_box(x$WORKPAR, x$CONTROL, p)
      last <- ib_to_work(x$THETA, x$WORKPAR)
      rej <- last + x$CONTROL$STEP * x$GAP[x$N_ITER, ]
      out <- which(rej < box$L | rej > box$U)
      bound <- ifelse(rej < box$L, box$L, box$U)[out]
      path <- rbind(
        path,
        data.frame(
          iterate = nrow(work),
          value = rej,
          param = factor(nm, nm),
          series = "rejected update"
        ),
        data.frame(
          iterate = rep(c(0, nrow(work)), each = length(out)),
          value = rep(bound, 2),
          param = factor(rep(nm[out], 2), nm),
          series = "box bound"
        )
      )
    }
    lv <- c("iterate", "rejected update", "box bound")
    keep <- lv %in% path$series
    path$series <- factor(path$series, lv[keep])
    tinyplot::tinyplot(
      value ~ iterate | series,
      facet = ~param,
      data = path,
      type = "o",
      pch = c(16, 4, NA)[keep],
      lty = c(1, 0, 2)[keep],
      facet.args = list(free = TRUE),
      ylab = "Working coordinate",
      main = "Paths"
    )
  }
  invisible(x)
}
