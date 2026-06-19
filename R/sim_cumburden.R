#' Simulate progressive multistate disease-burden data
#'
#' Generates two-arm data in the long multistate format expected by [auc_fit()]
#' from a continuous-time progressive (Markov) multistate model: subjects
#' advance through `K` ordinal non-fatal states of increasing severity and may
#' die (state `K + 1`) from any state, subject to random dropout and
#' administrative censoring. A treatment effect is induced by multiplying the
#' baseline transition intensities, so that an effective treatment reduces
#' cumulative disease burden and yields an AUC ratio below 1.
#'
#' @details
#' Subjects are allocated to the active arm (`arm = 1`) and the control arm
#' (`arm = 0`) in fixed numbers determined by `ratio` (not by independent
#' coin-flips), so the realized allocation matches the requested ratio up to
#' rounding. With `ratio = 1` and even `n` the split is exactly balanced.
#'
#' Each subject begins event-free (state 0). While in state \eqn{j}, two
#' independent exponential clocks compete:
#' \itemize{
#'   \item progression to state \eqn{j + 1} with intensity
#'     \eqn{\lambda_{j+1} \times f_p}, for \eqn{j = 0, \dots, K-1};
#'   \item death (transition to state \eqn{K + 1}) with intensity
#'     \eqn{\mu_j \times f_d}.
#' }
#' The earlier event occurs, time advances, and the process repeats from the new
#' state until death, random dropout (an exponential time with rate
#' `dropout_rate`), or administrative censoring at `tau`. Sojourn times are
#' therefore exponential and the model is time-homogeneous.
#'
#' The treatment effect enters only through the multipliers \eqn{f_p} and
#' \eqn{f_d}, which equal `hr_progression` and `hr_death` in the active arm and
#' `1` in the control arm. **The hazard ratio is thus defined at the level of
#' each transition intensity, not at the level of any composite endpoint.** A
#' value below 1 slows that family of transitions. Two consequences are worth
#' noting:
#' \itemize{
#'   \item If `mortality` is constant across states, the all-cause mortality
#'     hazard is constant and the overall-survival hazard ratio equals
#'     `hr_death` exactly.
#'   \item For the progression composites ("time to reach state \eqn{k} or
#'     worse") there is no single constant marginal hazard ratio, because those
#'     endpoints combine several competing and sequential transitions; a common
#'     multiplier uniformly delays the whole disease trajectory rather than
#'     implying one hazard ratio per composite.
#' }
#'
#' Setting `hr` alone applies a common multiplier to every transition (the
#' simplest case). Supplying `hr_progression` and `hr_death` separately lets the
#' treatment act differently on progression and on mortality. Supplying vector
#' `progression` and/or `mortality` lets the baseline hazards vary by severity
#' (for example, faster progression and higher mortality in more severe states).
#'
#' This is an illustrative generator for the examples, tests, and vignette; it
#' is not intended to reproduce any particular trial.
#'
#' @param n Total number of subjects (both arms combined).
#' @param K Number of ordinal non-fatal states (death is state `K + 1`).
#'   Default `4`.
#' @param tau Administrative censoring time (maximum follow-up). Default `6`.
#' @param ratio Allocation ratio of active to control subjects. `1` (the
#'   default) gives 1:1; `2` gives 2:1 active:control; `0.5` gives 1:2. The
#'   number allocated to the active arm is `round(n * ratio / (ratio + 1))`.
#' @param hr Convenience hazard ratio applied to all transition intensities in
#'   the active arm; used as the default for both `hr_progression` and
#'   `hr_death`. Default `0.8`.
#' @param hr_progression Hazard ratio applied to the progression intensities
#'   (state \eqn{j} to \eqn{j+1}) in the active arm. Defaults to `hr`.
#' @param hr_death Hazard ratio applied to the death intensities in the active
#'   arm. Defaults to `hr`.
#' @param progression Control-arm progression rates -- how quickly patients move
#'   to the next, more severe state, before any treatment effect is applied.
#'   Either a single rate used for every step, or a length-`K` vector giving a
#'   separate rate for each transition (state 0 to 1, 1 to 2, ..., K-1 to K). Use
#'   a vector to let progression speed differ by disease state (for example,
#'   faster early decline). Default `0.25`.
#' @param mortality Control-arm death rates -- how quickly patients die, before
#'   any treatment effect is applied. Either a single rate used in every state,
#'   or a length-`K+1` vector giving a separate death rate for each state a
#'   patient can die from (states 0, 1, ..., K). Use a vector to let mortality
#'   rise with severity. This vector is one element longer than `progression`
#'   because there are `K` steps between the non-fatal states but `K+1` states
#'   (0 through K) from which death can occur. Default `0.07`.
#' @param dropout_rate Intensity of random (non-administrative) censoring.
#'   Default `0.05`.
#' @param seed Optional integer random seed.
#'
#' @return A data frame in long multistate format with columns `id`, `time`,
#'   `status` (`0` censored, `1..K` non-fatal states, `K+1` death) and `arm`
#'   (`0` control, `1` active treatment), suitable for [auc_fit()].
#'
#' @seealso [auc_fit()]
#'
#' @examples
#' # Balanced 1:1 trial, common hazard ratio on every transition
#' sim <- sim_cumburden(n = 200, hr = 0.8, seed = 42)
#' table(arm = tapply(sim$arm, sim$id, `[`, 1))   # exactly 100 / 100
#'
#' # 2:1 allocation; treatment slows progression more than mortality,
#' # with severity-dependent baseline hazards
#' sim2 <- sim_cumburden(
#'   n = 600, K = 4, tau = 6, ratio = 2,
#'   hr_progression = 0.75, hr_death = 0.90,
#'   progression = c(0.30, 0.25, 0.20, 0.15),
#'   mortality   = c(0.02, 0.03, 0.05, 0.08, 0.12),
#'   seed = 1
#' )
#'
#' @export
sim_cumburden <- function(n = 500, K = 4, tau = 6, ratio = 1,
                          hr = 0.8, hr_progression = hr, hr_death = hr,
                          progression = 0.25, mortality = 0.07,
                          dropout_rate = 0.05, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)
  if (K < 1) stop("`K` must be at least 1.", call. = FALSE)
  if (ratio <= 0) stop("`ratio` must be positive.", call. = FALSE)

  # Baseline (control-arm) transition intensities.
  lambda <- if (length(progression) == 1L) rep(progression, K) else progression
  if (length(lambda) != K)
    stop("`progression` must have length 1 or K.", call. = FALSE)
  mu <- if (length(mortality) == 1L) rep(mortality, K + 1L) else mortality
  if (length(mu) != K + 1L)
    stop("`mortality` must have length 1 or K+1.", call. = FALSE)

  # Fixed allocation to the requested active:control ratio, then permuted.
  n_active <- round(n * ratio / (ratio + 1))
  if (n_active < 1L || n_active > n - 1L)
    stop("`n` and `ratio` leave an arm empty; increase `n` or adjust `ratio`.",
         call. = FALSE)
  arms <- sample(c(rep(1L, n_active), rep(0L, n - n_active)))

  out <- vector("list", n)

  for (i in seq_len(n)) {
    arm <- arms[i]
    fp  <- if (arm == 1L) hr_progression else 1   # progression multiplier
    fd  <- if (arm == 1L) hr_death       else 1   # death multiplier
    cens <- min(stats::rexp(1, dropout_rate), tau)

    rows_t <- numeric(0); rows_s <- integer(0)
    clock <- 0; state <- 0L; died <- FALSE

    repeat {
      death_rate_j <- mu[state + 1L] * fd
      t_death <- stats::rexp(1, death_rate_j)
      t_prog  <- if (state < K) stats::rexp(1, lambda[state + 1L] * fp) else Inf
      step    <- min(t_death, t_prog)

      if (clock + step > cens) break                # censored before next event
      clock <- clock + step

      if (t_death <= t_prog) {                      # death
        rows_t <- c(rows_t, clock); rows_s <- c(rows_s, K + 1L)
        died <- TRUE
        break
      } else {                                      # progress one state
        state  <- state + 1L
        rows_t <- c(rows_t, clock); rows_s <- c(rows_s, state)
      }
    }

    if (!died) {                                    # add censoring row
      rows_t <- c(rows_t, cens); rows_s <- c(rows_s, 0L)
    }

    out[[i]] <- data.frame(id = i, time = rows_t, status = rows_s, arm = arm)
  }

  do.call(rbind, out)
}
