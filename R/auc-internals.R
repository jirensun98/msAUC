# Internal (non-exported) helpers for the AUC cumulative-burden estimator.
# These are deliberately written in base R + the 'survival' package only, so
# that the package stays lightweight and passes R CMD check cleanly (no calls
# to unexported functions such as survival:::survmean, and no heavy
# dependencies such as 'pec').

#' Restricted mean survival time from a Kaplan-Meier curve
#'
#' Area under the Kaplan-Meier survival curve from 0 to each value in `at`.
#' This reproduces `survival:::survmean(km, rmean = tau)[["rmean"]]`: the curve
#' is held constant at its last observed value if `at` extends beyond the last
#' event time.
#'
#' @param time numeric event/censoring times.
#' @param ind  event indicator (1 = event, 0 = censored).
#' @param at   numeric vector of upper integration limits.
#' @return numeric vector of RMST values, one per element of `at`.
#' @keywords internal
#' @noRd
.km_rmst <- function(time, ind, at) {
  km <- survival::survfit(survival::Surv(time, ind) ~ 1)
  # Step function: S(0) = 1, dropping at each km$time, right-continuous.
  ktime <- c(0, km$time)
  ksurv <- c(1, km$surv)

  vapply(at, function(tau) {
    if (tau <= 0) return(0)
    # widths of each step interval, capped at tau
    upper <- pmin(ktime, tau)
    # length of interval [ktime[j], ktime[j+1]) intersected with [0, tau]
    n <- length(ktime)
    width <- pmin(c(ktime[-1], tau), tau) - upper
    width[width < 0] <- 0
    sum(ksurv * width)
  }, numeric(1))
}

#' Baseline survival from a null Cox model at a set of times
#'
#' Returns the survival curve from `coxph(Surv(time, ind) ~ 1)` evaluated at
#' `at`. This is the quantity used in the denominator of the influence
#' function. It is computed via the survival package directly (equivalent to
#' the value returned by `pec::predictSurvProb()` on a null Cox fit) so that
#' tie handling matches a standard Cox baseline estimate.
#'
#' @param time numeric event/censoring times.
#' @param ind  event indicator (1 = event, 0 = censored).
#' @param at   numeric vector of times at which to evaluate the survival curve.
#' @return numeric survival probabilities, one per element of `at`.
#' @keywords internal
#' @noRd
.cox_baseline_surv <- function(time, ind, at) {
  fit <- survival::coxph(survival::Surv(time, ind) ~ 1)
  summary(survival::survfit(fit), times = at, extend = TRUE)$surv
}

#' Influence function for the restricted mean event-free time of one transition
#'
#' For a single (composite) time-to-event endpoint within a single arm, returns
#' the per-subject influence-function contributions for the area under the
#' survival curve (restricted mean survival time, RMST) up to `tau`, together
#' with the point estimate `U` = RMST(tau).
#'
#' The estimator follows the martingale representation used in Sun et al. (2025)
#' for area-under-the-curve estimands: for subject i,
#' \deqn{\xi_i = \sum_{u \le \min(X_i, \tau)} w(u)\,[dN_i(u) - dN(u)/Y(u)] / S(u),}
#' where \eqn{w(u) = \int_u^\tau S(v)\,dv}, \eqn{S} is the survival function and
#' \eqn{Y(u)} the at-risk set.
#'
#' @param time numeric event/censoring times (one row per subject).
#' @param ind  event indicator (1 = event, 0 = censored).
#' @param tau  numeric restriction time.
#' @return list with elements `influence` (named numeric vector, one per
#'   subject, in input order) and `U` (scalar RMST estimate).
#' @keywords internal
#' @noRd
.auc_influence <- function(time, ind, tau) {
  n <- length(time)

  # Unique event times (ascending), restricted to <= tau.
  ev <- sort(unique(time[ind == 1]))
  ev <- ev[ev <= tau]

  U <- .km_rmst(time, ind, tau)

  if (length(ev) == 0L) {
    # No events at or before tau: influence is exactly zero for everyone.
    return(list(influence = rep(0, n), U = U))
  }

  # At-risk count and event count at each unique event time.
  Y          <- vapply(ev, function(u) sum(time >= u), numeric(1))
  d          <- vapply(ev, function(u) sum(time == u & ind == 1), numeric(1))
  event_rate <- d / Y

  # Survival in the denominator: baseline survival from a null Cox model,
  # matching pec::predictSurvProb() on coxph(~ 1).
  surv_prob <- .cox_baseline_surv(time, ind, ev)

  # Weight w(u) = RMST(tau) - RMST(u) = integral_u^tau S(v) dv.
  w <- U - .km_rmst(time, ind, ev)

  # Per-subject influence, summing martingale increments over event times the
  # subject was at risk for (u <= time_i, u <= tau).
  #   contrib_{i,u} = 1{time_i >= u} * w(u) * (1{time_i == u & event} - rate(u)) / S(u)
  scale_term <- w / surv_prob                       # length = #event times
  influence  <- numeric(n)
  for (j in seq_along(ev)) {
    at_risk    <- time >= ev[j]
    own_event  <- (time == ev[j]) & (ind == 1)
    increment  <- (own_event - event_rate[j])
    influence  <- influence + at_risk * scale_term[j] * increment
  }

  list(influence = influence, U = U)
}
