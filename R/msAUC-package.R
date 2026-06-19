#' msAUC: Area-Under-the-Curve Estimation of Cumulative Disease Burden
#'
#' Tools to quantify cumulative disease burden in trials with a progressive,
#' ordinal multistate outcome and a terminal event (death). The central
#' estimand is the area under the mean cumulative severity-score curve (AUC)
#' over a restriction window; treatment effects are reported as an AUC ratio and
#' an AUC difference, with a per-state decomposition. The user can also select
#' which ordinal states enter the composite.
#'
#' @section Main functions:
#' \describe{
#'   \item{[auc_fit()]}{Estimate the AUC ratio, difference, and decomposition.}
#'   \item{[auc_plot()]}{Plot the mean cumulative score curves and AUC ratio.}
#'   \item{[auc_curve_data()]}{Extract the underlying curve data.}
#'   \item{[sim_cumburden()]}{Simulate example multistate data.}
#' }
#'
#' @section Data format:
#' Functions expect long multistate data (one or more rows per subject) with a
#' status code that is `0` for censoring, `1, ..., K` for ordinal non-fatal
#' states in increasing severity, and `K + 1` for death. This matches the input
#' format of the \pkg{rmt} package.
#'
#' @seealso
#' For a step-by-step worked example, see the package tutorial: run
#' `vignette("msAUC-tutorial")` (or `browseVignettes("msAUC")` to list
#' it). The tutorial is built only when the package is installed with
#' `build_vignettes = TRUE`; its source file is
#' `vignettes/msAUC-tutorial.Rmd` in the GitHub repository.
#'
#' @importFrom ggplot2 .data
#' @keywords internal
"_PACKAGE"
