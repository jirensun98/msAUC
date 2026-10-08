#' Area under the mean cumulative disease-burden curve (AUC estimand)
#'
#' Estimates the area-under-the-curve (AUC) summary of cumulative disease
#' burden for a progressive multistate process, comparing two treatment arms.
#' The AUC integrates the mean cumulative severity score over a restriction
#' window \eqn{[0, \tau]} and summarises the treatment effect as a ratio and a
#' difference. The user chooses which ordinal disease states to include in the
#' composite; the included states are scored \eqn{1, \dots, m} in increasing
#' severity, where \eqn{m} is the number of states retained.
#'
#' @details
#' The data must be in the long multistate format also used by the \pkg{rmt}
#' package: one or more rows per subject, where `status` is `0` for censoring,
#' `1, ..., K` for the ordinal non-fatal states in increasing severity, and
#' `K + 1` for death. Each subject's records end with either death or a
#' censoring row.
#'
#' For each retained state \eqn{s} the estimator forms the "time to reach state
#' \eqn{s} or worse" endpoint (a subject who dies, or who is first observed in a
#' more severe state, is treated as having reached \eqn{s} at that time), then
#' estimates its restricted mean survival time (RMST) by the area under the
#' Kaplan-Meier curve. With \eqn{m} retained states the arm-specific AUC is
#' \deqn{\mathrm{AUC}_a = m\tau - \sum_{k=1}^{m} \mathrm{RMST}_k^{(a)},}
#' interpreted as the cumulative time lost to disease progression.
#' Standard errors use the martingale influence-function representation of each
#' RMST, in the form that remains valid with tied event times (for a single state it
#' reproduces the Greenwood-type variance of the Kaplan-Meier restricted mean);
#' the two arms are treated as independent.
#'
#' Selecting fewer states reproduces the endpoint-definition sensitivity
#' analysis in the source manuscript: dropping the mildest states lowers
#' \eqn{m} and re-scores the remaining states, while death always carries the
#' top score.
#'
#' @param data A data frame in long multistate format.
#' @param id,time,status,arm Column names (character scalars) in `data` giving
#'   the subject identifier, the event/censoring time, the status code
#'   (`0` censored, `1..K` non-fatal states, `K+1` death), and the treatment
#'   arm. Defaults are `"id"`, `"time"`, `"status"`, `"arm"`.
#' @param tau Numeric restriction time. Required.
#' @param states Optional integer vector of status codes to include in the
#'   composite, in any order. Defaults to all positive status codes observed in
#'   the data. Must be a subset of the observed positive codes.
#' @param trt The value of `arm` denoting the active treatment (numerator of the
#'   ratio). Defaults to the larger of the two arm values (e.g. `1` when arm is
#'   coded `0`/`1`). The other value is the control/reference arm.
#' @param conf.level Confidence level for intervals. Default `0.95`.
#'
#' @return An object of class `"aucfit"`: a list with components
#'   \describe{
#'     \item{ratio}{Named numeric vector: AUC ratio, SE, CI bounds, and p-value
#'       (treatment vs control). A ratio below 1 indicates lower burden under
#'       treatment.}
#'     \item{difference}{Named numeric vector: AUC difference (treatment minus
#'       control), SE, CI bounds, and p-value, in score-years. A negative
#'       difference indicates reduced burden under treatment.}
#'     \item{auc}{Named numeric vector of the two arm-specific AUC values.}
#'     \item{decomposition}{Data frame of per-state contributions: event-free
#'       time gained under treatment (treatment minus control RMST) for each
#'       retained transition, with SE, z, and p. The contributions sum to the
#'       absolute burden reduction (minus the AUC difference).}
#'     \item{endpoints}{Internal list of the per-arm composite survival data for
#'       each retained state (used by [auc_plot()]).}
#'     \item{tau, states, n_states, conf.level, arm_levels, n}{Settings and
#'       sample sizes used in the fit.}
#'   }
#'
#' @seealso [auc_plot()] to visualise the mean cumulative score curves and the
#'   AUC ratio over time; [sim_cumburden()] to simulate example data.
#'
#' @references
#' Claggett B, Tian L, Fu H, et al. Quantifying the totality of treatment effect
#' with multiple event-time observations in the presence of a terminal event.
#' \emph{Statistics in Medicine} 2018; 37(25): 3589--3598.
#'
#' Sun J, et al. Improve the precision of area under the curve
#' estimation for recurrent events through covariate adjustment.
#' \emph{Statistics in Medicine} 2025; 44(15--17): e70187.
#'
#' @examples
#' sim <- sim_cumburden(n = 400, tau = 6, seed = 1)
#' fit <- auc_fit(sim, tau = 6)
#' fit
#'
#' # Endpoint-definition sensitivity: drop the two mildest states.
#' fit2 <- auc_fit(sim, tau = 6, states = c(3, 4, 5))
#' fit2$ratio
#'
#' @export
auc_fit <- function(data, id = "id", time = "time", status = "status",
                    arm = "arm", tau, states = NULL, trt = NULL,
                    conf.level = 0.95) {

  if (missing(tau) || length(tau) != 1L || !is.finite(tau) || tau <= 0)
    stop("`tau` must be a single positive number.", call. = FALSE)
  if (!is.data.frame(data))
    stop("`data` must be a data frame.", call. = FALSE)
  for (nm in c(id, time, status, arm))
    if (!nm %in% names(data))
      stop(sprintf("Column '%s' not found in `data`.", nm), call. = FALSE)
  if (conf.level <= 0 || conf.level >= 1)
    stop("`conf.level` must be between 0 and 1.", call. = FALSE)

  idv  <- data[[id]]
  tmv  <- as.numeric(data[[time]])
  stv  <- data[[status]]
  arv  <- data[[arm]]

  if (any(!is.finite(tmv)))
    stop("`time` contains missing or non-finite values.", call. = FALSE)
  if (any(is.na(stv)) || any(stv != round(stv)) || any(stv < 0))
    stop("`status` must be a non-negative integer (0 = censored).", call. = FALSE)
  stv <- as.integer(stv)

  # --- treatment arm coding --------------------------------------------------
  arm_vals <- sort(unique(arv))
  if (length(arm_vals) != 2L)
    stop("`arm` must take exactly two distinct values.", call. = FALSE)
  if (is.null(trt)) trt <- arm_vals[2L]
  if (!trt %in% arm_vals)
    stop("`trt` must be one of the two values of `arm`.", call. = FALSE)
  ctrl <- setdiff(arm_vals, trt)

  # arm must be constant within subject
  idf <- factor(idv)
  arm_by_id <- tapply(as.character(arv), idf, function(x) x[1L])
  chk <- tapply(as.character(arv), idf, function(x) length(unique(x)))
  if (any(chk > 1L))
    stop("`arm` must be constant within each `id`.", call. = FALSE)

  # --- states to include -----------------------------------------------------
  obs_states <- sort(unique(stv[stv > 0]))
  if (length(obs_states) == 0L)
    stop("No events found (all `status` values are 0).", call. = FALSE)
  if (is.null(states)) {
    states <- obs_states
  } else {
    states <- sort(unique(as.integer(states)))
    if (!all(states %in% obs_states))
      stop("`states` must be a subset of the observed positive status codes: ",
           paste(obs_states, collapse = ", "), ".", call. = FALSE)
  }
  m <- length(states)

  # Per-subject follow-up end (censoring/death anchor).
  fu <- tapply(tmv, idf, max)

  # Build, for each retained state, the per-subject "time to state >= s" data.
  ids       <- levels(idf)
  arm_id    <- arm_by_id[ids]
  is_trt    <- arm_id == as.character(trt)
  n_trt     <- sum(is_trt)
  n_ctrl    <- sum(!is_trt)

  make_endpoint <- function(s) {
    qt   <- ifelse(stv >= s, tmv, Inf)
    minq <- tapply(qt, idf, min)[ids]           # Inf if subject never reaches s
    reached <- is.finite(minq)
    Ts   <- ifelse(reached, minq, fu[ids])
    list(time = as.numeric(Ts), ind = as.integer(reached))
  }

  endpoints <- vector("list", m)
  U_trt  <- numeric(m); U_ctrl <- numeric(m)
  IF_trt  <- matrix(0, nrow = n_trt,  ncol = m)
  IF_ctrl <- matrix(0, nrow = n_ctrl, ncol = m)

  for (j in seq_len(m)) {
    ep <- make_endpoint(states[j])
    tt <- ep$time; ii <- ep$ind

    inf_t <- .auc_influence(tt[is_trt],  ii[is_trt],  tau)
    inf_c <- .auc_influence(tt[!is_trt], ii[!is_trt], tau)

    U_trt[j]  <- inf_t$U
    U_ctrl[j] <- inf_c$U
    IF_trt[, j]  <- inf_t$influence
    IF_ctrl[, j] <- inf_c$influence

    endpoints[[j]] <- list(
      state = states[j],
      score = j,
      trt   = data.frame(time = tt[is_trt],  ind = ii[is_trt]),
      ctrl  = data.frame(time = tt[!is_trt], ind = ii[!is_trt])
    )
  }

  # --- aggregate to AUC, ratio, difference -----------------------------------
  Ut <- sum(U_trt); Uc <- sum(U_ctrl)
  AUC_trt  <- m * tau - Ut
  AUC_ctrl <- m * tau - Uc
  if (AUC_trt <= 0 || AUC_ctrl <= 0)
    stop("Non-positive AUC encountered; check `tau` and the data.", call. = FALSE)

  IF_tot_t <- rowSums(IF_trt)
  IF_tot_c <- rowSums(IF_ctrl)
  # Influence values sum to zero within arm, so sum(IF^2)/n^2 is the variance
  # estimator of the accompanying manuscript (Supplementary Section S1).
  varU_t <- sum(IF_tot_t^2) / n_trt^2
  varU_c <- sum(IF_tot_c^2) / n_ctrl^2

  z <- stats::qnorm(1 - (1 - conf.level) / 2)

  # Ratio (treatment / control) on the log scale.
  R <- AUC_trt / AUC_ctrl
  var_logR <- varU_t / AUC_trt^2 + varU_c / AUC_ctrl^2
  se_logR  <- sqrt(var_logR)
  ci_R     <- exp(log(R) + c(-1, 1) * z * se_logR)
  se_R     <- R * se_logR
  p_R      <- 2 * stats::pnorm(abs(log(R) / se_logR), lower.tail = FALSE)

  # Difference (treatment - control).
  D    <- AUC_trt - AUC_ctrl
  varD <- varU_t + varU_c
  se_D <- sqrt(varD)
  ci_D <- D + c(-1, 1) * z * se_D
  p_D  <- 2 * stats::pnorm(abs(D / se_D), lower.tail = FALSE)

  # --- per-state decomposition (event-free time gained under treatment) ------
  comp_est <- U_trt - U_ctrl            # RMST_trt - RMST_ctrl per state
  comp_se  <- sqrt(colSums(IF_trt^2) / n_trt^2 +
                   colSums(IF_ctrl^2) / n_ctrl^2)
  comp_z   <- comp_est / comp_se
  comp_p   <- 2 * stats::pnorm(abs(comp_z), lower.tail = FALSE)

  decomposition <- data.frame(
    state    = states,
    score    = seq_len(m),
    estimate = comp_est,
    se       = comp_se,
    z        = comp_z,
    p.value  = comp_p,
    row.names = NULL
  )

  ratio <- c(estimate = R, se = se_R,
             conf.low = ci_R[1], conf.high = ci_R[2], p.value = p_R)
  difference <- c(estimate = D, se = se_D,
                  conf.low = ci_D[1], conf.high = ci_D[2], p.value = p_D)
  auc_vals <- c(treatment = AUC_trt, control = AUC_ctrl)

  structure(
    list(
      ratio         = ratio,
      difference    = difference,
      auc           = auc_vals,
      decomposition = decomposition,
      endpoints     = endpoints,
      tau           = tau,
      states        = states,
      n_states      = m,
      conf.level    = conf.level,
      arm_levels    = c(treatment = as.character(trt), control = as.character(ctrl)),
      n             = c(treatment = n_trt, control = n_ctrl)
    ),
    class = "aucfit"
  )
}
