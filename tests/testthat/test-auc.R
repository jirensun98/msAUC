test_that(".km_rmst reproduces survival::survmean", {
  set.seed(1)
  for (n in c(40, 120, 300)) {
    time <- rexp(n, 0.3); ind <- rbinom(n, 1, 0.7)
    km <- survival::survfit(survival::Surv(time, ind) ~ 1)
    taus <- c(0.5, 1, 2, 4, max(time) + 1)
    ref  <- vapply(taus, function(t)
      as.numeric(survival:::survmean(km, rmean = t)[[1]]["rmean"]), numeric(1))
    mine <- msAUC:::.km_rmst(time, ind, taus)
    expect_equal(mine, ref, tolerance = 1e-12)
  }
})

test_that(".auc_influence gives the Greenwood-type RMST variance (early censoring, ties)", {
  set.seed(7)
  for (n in c(60, 150, 400)) {
    t_ev <- ceiling(rexp(n, 0.25) * 2) / 2      # half-unit grid: many tied events
    cens <- runif(n, 0.5, 6)                    # censoring well before tau
    time <- pmin(t_ev, cens); ind <- as.integer(t_ev <= cens)
    tau  <- 5
    b  <- msAUC:::.auc_influence(time, ind, tau)
    km <- survival::survfit(survival::Surv(time, ind) ~ 1)
    keep <- km$n.event > 0 & km$time <= tau
    u <- km$time[keep]; r <- km$n.risk[keep]; e <- km$n.event[keep]
    w  <- b$U - msAUC:::.km_rmst(time, ind, u)
    gw <- sum(ifelse(r > e, w^2 * e / (r * (r - e)), 0))
    expect_equal(sum(b$influence^2) / n^2, gw, tolerance = 1e-10)
    expect_equal(sum(b$influence), 0, tolerance = 1e-8)
    expect_equal(b$U, as.numeric(survival:::survmean(km, rmean = tau)[[1]]["rmean"]),
                 tolerance = 1e-12)
  }
})

test_that("auc_fit returns a well-formed object", {
  sim <- sim_cumburden(n = 300, tau = 6, seed = 3)
  fit <- auc_fit(sim, tau = 6)
  expect_s3_class(fit, "aucfit")
  expect_named(fit$ratio, c("estimate", "se", "conf.low", "conf.high", "p.value"))
  expect_named(fit$difference, c("estimate", "se", "conf.low", "conf.high", "p.value"))
  expect_equal(fit$n_states, 5L)
  expect_true(all(fit$auc > 0))
})

test_that("ratio and difference are internally consistent", {
  sim <- sim_cumburden(n = 300, tau = 6, seed = 4)
  fit <- auc_fit(sim, tau = 6)
  expect_equal(unname(fit$ratio["estimate"]),
               unname(fit$auc["treatment"] / fit$auc["control"]),
               tolerance = 1e-10)
  expect_equal(unname(fit$difference["estimate"]),
               unname(fit$auc["treatment"] - fit$auc["control"]),
               tolerance = 1e-10)
})

test_that("decomposition sums to the negative AUC difference", {
  sim <- sim_cumburden(n = 300, tau = 6, seed = 5)
  fit <- auc_fit(sim, tau = 6)
  expect_equal(sum(fit$decomposition$estimate),
               unname(-fit$difference["estimate"]),
               tolerance = 1e-10)
})

test_that("endpoint selection changes the number of states and scoring", {
  sim <- sim_cumburden(n = 300, tau = 6, seed = 6)
  fit_all <- auc_fit(sim, tau = 6)
  fit_sub <- auc_fit(sim, tau = 6, states = c(4, 5))
  expect_equal(fit_all$n_states, 5L)
  expect_equal(fit_sub$n_states, 2L)
  expect_equal(fit_sub$decomposition$score, c(1L, 2L))
  expect_equal(fit_sub$states, c(4, 5))
})

test_that("including all states reproduces the manual composite endpoints", {
  # The state-5 (death) endpoint built by auc_fit must match a direct
  # construction: time = subject follow-up max, event = died.
  sim <- sim_cumburden(n = 200, tau = 6, seed = 8)
  fit <- auc_fit(sim, tau = 6)
  death_ep <- fit$endpoints[[5]]            # highest state = death
  manual <- do.call(rbind, by(sim, sim$id, function(d) {
    data.frame(time = max(d$time), ind = as.integer(any(d$status == 5)),
               arm = d$arm[1])
  }))
  m_trt  <- manual[manual$arm == 1, ]
  got_t  <- death_ep$trt[order(death_ep$trt$time), ]
  exp_t  <- m_trt[order(m_trt$time), c("time", "ind")]
  expect_equal(sort(got_t$time), sort(exp_t$time), tolerance = 1e-10)
  expect_equal(sum(got_t$ind), sum(exp_t$ind))
})

test_that("input validation triggers informative errors", {
  sim <- sim_cumburden(n = 100, seed = 9)
  expect_error(auc_fit(sim), "tau")
  expect_error(auc_fit(sim, tau = 6, states = c(1, 99)), "subset")
  bad <- sim; bad$arm <- 1
  expect_error(auc_fit(bad, tau = 6), "two distinct")
})

test_that("sim_cumburden yields a valid multistate data frame", {
  sim <- sim_cumburden(n = 150, K = 4, seed = 10)
  expect_true(all(c("id", "time", "status", "arm") %in% names(sim)))
  expect_true(all(sim$status %in% 0:5))
  expect_true(all(sim$arm %in% c(0, 1)))
  expect_true(all(sim$time > 0))
})

test_that("auc_curve_data returns usable curves and ratio", {
  sim <- sim_cumburden(n = 200, tau = 6, seed = 11)
  fit <- auc_fit(sim, tau = 6)
  cd <- auc_curve_data(fit)
  expect_true(all(c("arm", "time", "score") %in% names(cd$curves)))
  expect_true(all(c("time", "ratio") %in% names(cd$ratio)))
  expect_true(all(cd$ratio$time <= 6))
})
