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

#' Influence function for the restricted mean event-free time of one transition
#'
#' For a single (composite) time-to-event endpoint within a single arm, returns
#' the per-subject influence-function contributions for the area under the
#' survival curve (restricted mean survival time, RMST) up to `tau`, together
#' with the point estimate `U` = RMST(tau).
#'
#' For subject i,
#' \deqn{\xi_i = \sum_{u \le \min(X_i, \tau)} \frac{n\,w(u)}{Y(u) - d(u)}
#'   \{dN_i(u) - d(u)/Y(u)\},}
#' where \eqn{w(u) = \int_u^\tau \hat S(v)\,dv}, \eqn{Y(u)} is the number at
#' risk and \eqn{d(u)} the number of events at \eqn{u}. The weight
#' \eqn{n/\{Y(u) - d(u)\} = 1/[\hat\pi(u)\{1 - \Delta\hat\Lambda(u)\}]}, with
#' \eqn{\hat\pi(u) = Y(u)/n} the at-risk proportion, is the first-order form
#' that remains valid with tied event times: with it, \eqn{\sum_i \xi_i^2/n^2} equals
#' the Greenwood-type variance of the Kaplan-Meier restricted mean. Terms with
#' \eqn{Y(u) = d(u)} have \eqn{w(u) = 0} and are set to zero. (Versions <= 0.1.0
#' divided by the estimated survival function instead, which understates the
#' standard error when subjects are censored before `tau`.)
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

  # Weight w(u) = RMST(tau) - RMST(u) = integral_u^tau S(v) dv.
  w <- U - .km_rmst(time, ind, ev)

  # n * w(u) / {Y(u) - d(u)}  =  w(u) / [pi_hat(u) * {1 - dLambda_hat(u)}].
  # If Y(u) == d(u) the KM curve is zero from u onward, so w(u) == 0.
  scale_term <- ifelse(Y > d, n * w / pmax(Y - d, 1), 0)

  # Per-subject influence, summing martingale increments over event times the
  # subject was at risk for (u <= time_i, u <= tau).
  #   contrib_{i,u} = 1{time_i >= u} * scale(u) * (1{time_i == u & event} - rate(u))
  influence  <- numeric(n)
  for (j in seq_along(ev)) {
    at_risk    <- time >= ev[j]
    own_event  <- (time == ev[j]) & (ind == 1)
    increment  <- (own_event - event_rate[j])
    influence  <- influence + at_risk * scale_term[j] * increment
  }

  list(influence = influence, U = U)
}
