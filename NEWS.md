# msAUC 0.1.1

* Corrected the influence-function weight used for standard errors. The
  denominator is now the at-risk proportion times `1 - dLambda(u)`, i.e.
  `(Y(u) - d(u)) / n`, instead of the estimated survival function. Version
  0.1.0 understated standard errors when subjects were censored before `tau`.
  For a single state the variance now equals the Greenwood-type variance of the
  Kaplan-Meier restricted mean, including with tied event times.
* `auc_fit()` now computes variances as `sum(IF^2) / n^2`, matching the
  published formula (influence values sum to zero within arm).
* Point estimates are unchanged. Removed the unused internal
  `.cox_baseline_surv()`.

# msAUC 0.1.0

* Initial release.
