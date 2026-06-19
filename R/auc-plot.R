# Visualisation of the AUC cumulative-burden estimand.

#' Extract mean cumulative score and AUC-ratio curves from a fit
#'
#' Computes, on a common time grid, the expected cumulative severity score in
#' each arm (the quantity whose area is the AUC) and the AUC ratio as a function
#' of the restriction time. This is the data underlying [auc_plot()]; it is
#' exposed so users can build their own figures.
#'
#' @param x An object of class `"aucfit"` from [auc_fit()].
#' @param ratio_grid Optional numeric vector of restriction times at which to
#'   evaluate the AUC ratio. Defaults to 200 points spanning `[tau/10, tau]`
#'   (the ratio is numerically unstable as the restriction time approaches 0).
#'
#' @return A list with two data frames:
#'   \describe{
#'     \item{curves}{`arm`, `time`, `score` (= \eqn{m - \sum_k S_k(t)}), the
#'       expected cumulative score in each arm.}
#'     \item{ratio}{`time`, `ratio` (= \eqn{\mathrm{AUC}_{trt}(t)/\mathrm{AUC}_{ctrl}(t)}).}
#'   }
#'
#' @seealso [auc_plot()]
#' @export
auc_curve_data <- function(x, ratio_grid = NULL) {
  stopifnot(inherits(x, "aucfit"))
  m   <- x$n_states
  tau <- x$tau

  # ---- expected cumulative score curves, by arm ----------------------------
  score_curve <- function(side) {
    times <- sort(unique(unlist(lapply(x$endpoints, function(ep) {
      d <- ep[[side]]; d$time[d$ind == 1]
    }))))
    times <- times[times <= tau]
    if (length(times) == 0L) times <- tau
    surv_mat <- vapply(x$endpoints, function(ep) {
      d  <- ep[[side]]
      sf <- survival::survfit(survival::Surv(d$time, d$ind) ~ 1)
      summary(sf, times = times, extend = TRUE)$surv
    }, numeric(length(times)))
    if (is.null(dim(surv_mat))) surv_mat <- matrix(surv_mat, nrow = length(times))
    data.frame(time = times, score = m - rowSums(surv_mat))
  }

  c_trt  <- score_curve("trt")
  c_ctrl <- score_curve("ctrl")
  c_trt$arm  <- x$arm_levels["treatment"]
  c_ctrl$arm <- x$arm_levels["control"]
  curves <- rbind(c_ctrl[c("arm", "time", "score")],
                  c_trt [c("arm", "time", "score")])

  # ---- AUC ratio as a function of the restriction time ---------------------
  if (is.null(ratio_grid))
    ratio_grid <- seq(tau / 10, tau, length.out = 200)
  ratio_grid <- ratio_grid[ratio_grid > 0 & ratio_grid <= tau]

  rmst_sum <- function(side) {
    mat <- vapply(x$endpoints, function(ep) {
      d <- ep[[side]]
      .km_rmst(d$time, d$ind, ratio_grid)
    }, numeric(length(ratio_grid)))
    if (is.null(dim(mat))) mat <- matrix(mat, nrow = length(ratio_grid))
    rowSums(mat)
  }
  auc_trt  <- m * ratio_grid - rmst_sum("trt")
  auc_ctrl <- m * ratio_grid - rmst_sum("ctrl")
  ratio <- data.frame(time = ratio_grid, ratio = auc_trt / auc_ctrl)

  list(curves = curves, ratio = ratio)
}

#' Plot the mean cumulative disease-burden curves and AUC ratio
#'
#' Reproduces the key figure for the AUC estimand: the expected cumulative
#' severity score over time in each arm (primary axis), optionally overlaid with
#' the AUC ratio as a function of the restriction time (secondary axis).
#'
#' @param x An object of class `"aucfit"` from [auc_fit()].
#' @param show_ratio Logical; overlay the AUC-ratio curve on a secondary axis.
#'   Default `TRUE`.
#' @param ratio_grid Optional numeric vector of times for the ratio curve (see
#'   [auc_curve_data()]).
#' @param colors Length-2 vector of colours for control and treatment arms.
#' @param ratio_color Colour of the AUC-ratio curve.
#' @param xlab,ylab,ratio_lab Axis labels.
#'
#' @return A \pkg{ggplot2} object, which can be further customised or saved with
#'   [ggplot2::ggsave()].
#'
#' @seealso [auc_fit()], [auc_curve_data()]
#'
#' @examples
#' sim <- sim_cumburden(n = 400, tau = 6, seed = 1)
#' fit <- auc_fit(sim, tau = 6)
#' p <- auc_plot(fit)
#' # print(p)
#'
#' @export
auc_plot <- function(x, show_ratio = TRUE, ratio_grid = NULL,
                     colors = c("#2C7FB8", "#D95F0E"),
                     ratio_color = "grey50",
                     xlab = "Time from baseline",
                     ylab = "Expected cumulative score",
                     ratio_lab = "AUC ratio") {
  stopifnot(inherits(x, "aucfit"))
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for auc_plot(). Please install it.",
         call. = FALSE)

  dat <- auc_curve_data(x, ratio_grid = ratio_grid)
  curves <- dat$curves
  ratio  <- dat$ratio

  ctrl_lab <- x$arm_levels["control"]
  trt_lab  <- x$arm_levels["treatment"]
  curves$arm <- factor(curves$arm, levels = c(ctrl_lab, trt_lab))
  col_map <- stats::setNames(colors, c(ctrl_lab, trt_lab))

  y_max <- max(curves$score, na.rm = TRUE)
  r_max <- max(1, max(ratio$ratio, na.rm = TRUE))

  p <- ggplot2::ggplot() +
    ggplot2::geom_step(
      data = curves,
      ggplot2::aes(x = .data$time, y = .data$score, color = .data$arm),
      linewidth = 0.6
    ) +
    ggplot2::scale_color_manual(name = NULL, values = col_map) +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")

  if (show_ratio) {
    p <- p +
      ggplot2::geom_line(
        data = ratio,
        ggplot2::aes(x = .data$time, y = .data$ratio * y_max / r_max),
        color = ratio_color, linewidth = 0.6, linetype = "dashed"
      ) +
      ggplot2::scale_y_continuous(
        name = ylab,
        limits = c(0, y_max),
        sec.axis = ggplot2::sec_axis(
          ~ . * r_max / y_max, name = ratio_lab
        )
      )
  }

  p
}
