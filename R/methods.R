# S3 methods for objects of class "aucfit".

#' @export
print.aucfit <- function(x, digits = 4, ...) {
  cat("Cumulative disease-burden AUC analysis\n")
  cat(sprintf("  Restriction time tau = %s\n", format(x$tau, digits = digits)))
  cat(sprintf("  States included: %s  (scored 1..%d, death highest)\n",
              paste(x$states, collapse = ", "), x$n_states))
  cat(sprintf("  Arms: treatment = '%s' (n = %d), control = '%s' (n = %d)\n",
              x$arm_levels["treatment"], x$n["treatment"],
              x$arm_levels["control"],   x$n["control"]))
  cat("\n")
  cat(sprintf("  AUC (treatment) = %s\n", format(x$auc["treatment"], digits = digits)))
  cat(sprintf("  AUC (control)   = %s\n", format(x$auc["control"],   digits = digits)))
  cat("\n")
  ci_pct <- round(100 * x$conf.level)
  cat(sprintf("  AUC ratio       = %s   %d%% CI (%s, %s)   p = %s\n",
              format(x$ratio["estimate"], digits = digits), ci_pct,
              format(x$ratio["conf.low"],  digits = digits),
              format(x$ratio["conf.high"], digits = digits),
              format.pval(x$ratio["p.value"], digits = digits)))
  cat(sprintf("  AUC difference  = %s   %d%% CI (%s, %s)   p = %s\n",
              format(x$difference["estimate"], digits = digits), ci_pct,
              format(x$difference["conf.low"],  digits = digits),
              format(x$difference["conf.high"], digits = digits),
              format.pval(x$difference["p.value"], digits = digits)))
  cat("\n  A ratio < 1 / difference < 0 indicates lower burden under treatment.\n")
  invisible(x)
}

#' Summarise a cumulative-burden AUC fit
#'
#' @param object An object of class `"aucfit"`.
#' @param ... Ignored.
#' @return The `object`, invisibly, after printing the overall effect and the
#'   per-state decomposition.
#' @export
summary.aucfit <- function(object, ...) {
  print(object)
  cat("\nPer-state decomposition")
  cat("\n(event-free time gained under treatment; sums to -AUC difference)\n\n")
  dec <- object$decomposition
  dec_disp <- data.frame(
    state    = dec$state,
    score    = dec$score,
    estimate = round(dec$estimate, 4),
    se       = round(dec$se, 4),
    z        = round(dec$z, 3),
    p.value  = format.pval(dec$p.value, digits = 3)
  )
  print(dec_disp, row.names = FALSE)
  cat(sprintf("\n  Sum of contributions = %.4f  (= -AUC difference = %.4f)\n",
              sum(dec$estimate), -object$difference["estimate"]))
  invisible(object)
}

#' @rdname auc_plot
#' @param ... Additional arguments passed to [auc_plot()].
#' @export
plot.aucfit <- function(x, ...) {
  auc_plot(x, ...)
}
