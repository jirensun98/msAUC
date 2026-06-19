# msAUC

<!-- badges: start -->
<!-- badges: end -->

**msAUC** quantifies *cumulative disease burden* in clinical trials whose
primary outcome is a progressive, ordinal multistate process with a terminal
event (death). Rather than reducing each patient to their first event, it
integrates the **mean cumulative severity score** over a restriction window
`[0, tau]` and summarises the treatment effect as an **AUC ratio** and an **AUC
difference**, with a per-state decomposition.

The method was developed to capture the totality of chronic kidney disease
progression (e.g. successive eGFR-decline thresholds, end-stage renal disease,
and death) in outcome trials, where conventional time-to-first
analyses discard information about how far and how fast patients progress.

## Installation

You can install the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("jirensun98/msAUC", build_vignettes = TRUE)
```

`build_vignettes = TRUE` makes the tutorial available inside R (see **Tutorial**
below) and needs `knitr`, `rmarkdown`, and pandoc; drop it for a faster install
if you don't need the tutorial.

## Tutorial

A step-by-step worked example is provided as a package vignette. After
installing with `build_vignettes = TRUE`, open it in R with

``` r
vignette("msAUC-tutorial")   # the rendered tutorial
browseVignettes("msAUC")      # or list all vignettes
```

You can also read the source on GitHub at
[`vignettes/msAUC-tutorial.Rmd`](https://github.com/jirensun98/msAUC/blob/main/vignettes/msAUC-tutorial.Rmd)
— GitHub shows the source text, while the command above gives the rendered
version with output and figures.

## Data format

`msAUC` uses the same long multistate layout as the **rmt** package: one or
more rows per subject, with a `status` code that is

- `0` for censoring,
- `1, ..., K` for the ordinal non-fatal states in **increasing** severity,
- `K + 1` for death (the most severe state).

Each subject's records end with either a death row or a censoring row.

| id | time | status | arm |
|---:|-----:|-------:|----:|
| 1  | 1.8  | 1      | 0   |
| 1  | 3.2  | 2      | 0   |
| 1  | 4.0  | 0      | 0   |
| 2  | 2.5  | 5      | 1   |

## Quick start

``` r
library(msAUC)

# Simulate example data: 4 ordinal states + death, treatment slows progression
sim <- sim_cumburden(n = 600, tau = 6, hr = 0.8, seed = 1)

# Fit the AUC cumulative-burden model
fit <- auc_fit(sim, tau = 6)
fit
summary(fit)      # adds the per-state decomposition

# Figure: mean cumulative score by arm + AUC ratio over time
auc_plot(fit)
```

`auc_fit()` reports the AUC ratio (treatment / control) and the AUC difference
(treatment − control). A **ratio below 1** or a **difference below 0** indicates
*lower* cumulative burden under treatment.

## Choosing which endpoints to include

The argument `states` controls which ordinal states enter the composite. The
retained states are re-scored `1, ..., m` in increasing severity (death always
carries the top score), reproducing the endpoint-definition sensitivity
analysis from the source manuscript:

``` r
auc_fit(sim, tau = 6)                       # all states
auc_fit(sim, tau = 6, states = c(2, 3, 4, 5))  # drop the mildest state
auc_fit(sim, tau = 6, states = c(4, 5))        # ESRD + death only
```

## Method summary

At any time `t`, a subject's *cumulative severity score* `N(t)` is the number of
severity thresholds they have crossed by `t`: it is 0 while they are event-free,
steps up by 1 each time they reach a more severe state, and tops out at `m` (the
number of retained states) at death. Averaging `N(t)` over the subjects in an
arm gives that arm's mean cumulative-score curve, which climbs as the cohort
progresses. The estimand is the **area under that mean curve** over `[0, tau]` —
a single number capturing both how severe events are and how long patients spend
in each state.

That area is obtained from survival curves rather than by integrating the score
directly. For each retained state `k`, `msAUC` forms the endpoint "time to
reach state `k` or worse" and estimates its restricted mean survival time
`RMST_k` — the average event-free time before that threshold is crossed, equal
to the area under the matching Kaplan–Meier curve. Summing over the `m`
thresholds,

```
AUC_a = m * tau - sum_k RMST_k(a)
```

Why the identity holds: `m * tau` is the most burden the window could possibly
hold, and each `RMST_k` is time the score stayed *below* threshold `k`;
subtracting that event-free time at every threshold leaves exactly the burden
accumulated. So a **lower** AUC means **less** burden.

The two arms are compared by the ratio `AUC_1 / AUC_0` and the difference
`AUC_1 - AUC_0` (treatment vs. control): a ratio below 1, or a difference below
0, indicates lower burden under treatment. The difference splits additively
across the thresholds, so each severity level contributes a known share of the
overall effect — the event-free time gained under treatment at that level.
Standard errors come from the martingale influence-function representation of
each RMST, with the two arms treated as independent.

## References

- Claggett B, Tian L, Fu H, et al. Quantifying the totality of treatment effect
  with multiple event-time observations in the presence of a terminal event.
  *Statistics in Medicine* 2018; 37(25): 3589–3598.
  <https://doi.org/10.1002/sim.7907>
- Sun J, et al. Improve the precision of area under the curve
  estimation for recurrent events through covariate adjustment. *Statistics in
  Medicine* 2025; 44(15–17): e70187. <https://doi.org/10.1002/sim.70187>

## License

MIT © Jiren Sun
