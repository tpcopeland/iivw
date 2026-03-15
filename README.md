# iivw

![R 4.0+](https://img.shields.io/badge/R-%E2%89%A5%204.0-brightgreen) ![MIT License](https://img.shields.io/badge/License-MIT-blue)

Inverse intensity of visit weighting for longitudinal data with informative visit processes.

## The Problem

In observational longitudinal studies, patients are not seen on a fixed schedule. Sicker patients tend to visit more often, meaning the observed data over-represents high-severity timepoints. Naive analyses of these data produce biased treatment effect estimates because the outcome influences when it is measured.

`iivw` corrects this by computing observation-level weights that account for the informative visit process, restoring the analysis to what would have been observed under a non-informative (e.g., fixed-interval) design.

## Weighting Methods

| Method | Corrects for | Model |
|--------|-------------|-------|
| **IIW** | Informative visit timing | Andersen-Gill recurrent-event Cox |
| **IPTW** | Confounding by indication | Cross-sectional logistic regression |
| **FIPTIW** | Both simultaneously | IIW x IPTW (product weight) |

**IIW** (inverse intensity weighting; Buzkova & Lumley 2007) models the visit process as a recurrent event using an Andersen-Gill counting-process Cox model. Covariates that drive visit frequency (e.g., disease severity, recent clinical events) enter the intensity model. Each observation receives weight exp(-xb), down-weighting visits that occurred because the patient was sicker and up-weighting visits from patients who visited less frequently than expected.

**IPTW** (inverse probability of treatment weighting) fits a logistic propensity score model on cross-sectional data (one row per subject) to estimate P(treatment | confounders). Weights are always stabilized using the marginal treatment prevalence in the numerator, preventing extreme weights from subjects with near-certain treatment assignment.

**FIPTIW** (Tompkins et al. 2025) combines both corrections as a simple product. This is the recommended approach when treatment assignment is non-random and visit timing is outcome-dependent, which is the typical situation in observational registry data.

Stabilized IIW weights are also supported via `stabcov`, fitting a marginal intensity model in the numerator to further reduce weight variability.

## Installation

```r
remotes::install_github("tpcopeland/iivw")
```

## Functions

| Function | Purpose |
|----------|---------|
| `iivw_weight()` | Compute IIW, IPTW, or FIPTIW weights |
| `iivw_fit()` | Fit weighted outcome model (GEE or mixed) |

`iivw_weight()` handles the full weight pipeline: counting-process construction, Cox model fitting, propensity score estimation, weight computation, optional truncation, and diagnostic reporting. It returns an S3 object with metadata so that `iivw_fit()` can automatically retrieve the weight variable, panel structure, and weight type.

`iivw_fit()` fits the weighted outcome model using either GEE-equivalent estimation (GLM with clustered sandwich SEs via `sandwich::vcovCL`, as required by IIW theory) or mixed-effects models (via `lme4`). Time can enter the model as linear, polynomial, or natural cubic spline terms. Bootstrap SEs are available for inference that accounts for the weight estimation step.

## Quick Start

### IIW only (correct for visit timing)

```r
library(iivw)

d <- iivw_weight(data, id = "id", time = "days",
                 visit_cov = c("edss", "relapse"))

fit <- iivw_fit(edss ~ relapse, data = d, timespec = "linear")
```

### FIPTIW (correct for visit timing + treatment confounding)

```r
d <- iivw_weight(data, id = "id", time = "days",
                 visit_cov = c("edss", "relapse"),
                 treat = "treated",
                 treat_cov = c("age", "sex", "edss_bl"),
                 truncate = c(1, 99))

fit <- iivw_fit(edss ~ treated + age + sex + edss_bl,
                data = d, timespec = "quadratic")
```

### FIPTIW with treatment-time interaction

```r
fit <- iivw_fit(edss ~ treated + age + sex + edss_bl,
                data = d, timespec = "ns(3)",
                interaction = "treated")
```

### Categorical treatment (3+ levels)

```r
fit <- iivw_fit(edss ~ treatment + age + sex + edss_bl,
                data = d, timespec = "ns(3)",
                categorical = "treatment",
                interaction = "treatment")
```

### Lagged covariates

```r
d <- iivw_weight(data, id = "id", time = "days",
                 visit_cov = c("edss", "relapse"),
                 lagvars = c("edss", "relapse"))
```

## Key Options

### `iivw_weight()`

| Argument | Description |
|----------|-------------|
| `data` | Data frame in long panel format (one row per visit) |
| `id` | Column name for subject identifier (required) |
| `time` | Column name for continuous visit time (required) |
| `visit_cov` | Column names for visit intensity Cox model (required) |
| `treat` | Column name for binary (0/1) treatment (triggers FIPTIW) |
| `treat_cov` | Column names for propensity score model (defaults to `visit_cov`) |
| `stabcov` | Column names for stabilized IIW numerator model |
| `lagvars` | Time-varying covariates to lag by one visit |
| `truncate` | Percentile trimming bounds, e.g. `c(1, 99)` |
| `wtype` | Override auto-detection: `"iivw"`, `"iptw"`, or `"fiptiw"` |
| `entry` | Column name for subject-specific study entry time (default: 0) |
| `prefix` | Prefix for weight column names (default: `"iivw_"`) |

### `iivw_fit()`

| Argument | Description |
|----------|-------------|
| `formula` | Outcome model formula (e.g., `y ~ treat + age`) |
| `data` | Data frame from `iivw_weight()` |
| `model` | `"gee"` (default) or `"mixed"` |
| `timespec` | `"linear"`, `"quadratic"`, `"cubic"`, `"ns(k)"`, or `"none"` |
| `family` | GLM family (default: `gaussian()`) |
| `link` | GLM link function override (default: canonical) |
| `interaction` | Covariates to interact with time terms |
| `categorical` | Variables to expand into dummy indicators |
| `basecat` | Named list of reference categories, e.g. `list(arm = 0)` |
| `bootstrap` | Number of bootstrap replicates (0 = sandwich SE) |
| `cluster` | Override clustering variable (default: id from metadata) |
| `level` | Confidence level (default: 95) |

## Return Values

### `iivw_weight()` returns:

A data frame of class `"iivw_weights"` with added weight columns (`iivw_iw`, `iivw_tw`, `iivw_weight`) and attributes containing diagnostics: N, n_ids, mean/SD/min/max/P1/P99 of weights, ESS, and n_truncated.

### `iivw_fit()` returns:

An S3 object of class `"iivw_fit"` with `print()`, `summary()`, `coef()`, and `vcov()` methods. Contains coefficients, standard errors, z-statistics, p-values, confidence intervals, the underlying model object, and metadata about the time specification, interactions, and categorical expansions used.

## Cross-Validation

The package QA suite includes cross-validation against two independent R implementations.

### Part A: IrregLong (Phenobarb dataset)

Cross-validated against the [IrregLong](https://CRAN.R-project.org/package=IrregLong) package (Pullenayegum 2020) using the Phenobarb pharmacokinetic dataset (59 neonates, 155 concentration measurements at irregular times).

| Test | Result |
|------|--------|
| Cox coefficients vs `coxph` reference | Match within 0.001 |
| IIW weight correlation with IrregLong `iiw.weights()` | r > 0.9 |

### Part B: FIPTIW Simulation (Tompkins et al. 2025)

Cross-validated against the Tompkins et al. R implementation using their simulation DGP (200 subjects, thinned Poisson visit process, known treatment effect beta = 0.5).

| Test | Result |
|------|--------|
| Cox coefficients (D, Wt, Z) vs reference | Match within 0.01 |
| IPTW weights vs reference | Max difference < 0.001 |
| FIPTIW weight correlation with reference | r > 0.75 |
| Treatment effect recovery (beta = 0.5) | Within 1.0 of truth |

### Test Summary

| Suite | Tests | Status |
|-------|-------|--------|
| Functional | 42 | All pass |
| Validation | 15 | All pass |
| Cross-validation | 9 | All pass |

## Dependencies

**Imports:** survival, sandwich, splines, stats

**Suggests:** lme4 (for mixed models), testthat, IrregLong, nlme

## Also Available

A Stata implementation is available at [tpcopeland/Stata-Tools](https://github.com/tpcopeland/Stata-Tools/tree/main/iivw).

## Author

Timothy P Copeland
Department of Clinical Neuroscience, Karolinska Institutet

## License

MIT

## References

- Buzkova P, Lumley T (2007). Longitudinal data analysis for generalized linear models with follow-up dependent on outcome-related variables. *Canadian Journal of Statistics* 35(4):485-500.
- Tompkins G, Dubin JA, Wallace M (2025). On flexible inverse probability of treatment and intensity weighting. *Statistical Methods in Medical Research*.
- Pullenayegum EM (2016). Multiple outputation for the analysis of longitudinal data subject to irregular observation. *Statistics in Medicine* 35(11):1800-1818.
- Lin H, Scharfstein DO, Rosenheck RA (2004). Analysis of longitudinal data with irregular, outcome-dependent follow-up. *Journal of the Royal Statistical Society B* 66(3):791-813.
