# iivw 1.0.0

* Initial release.
* `iivw_weight()`: Compute IIW, IPTW, and FIPTIW weights for longitudinal
  data with irregular visit patterns.
* `iivw_fit()`: Fit weighted outcome models using GEE or mixed model
  estimation with flexible time specifications.
* Supports stabilized IIW weights via `stabcov`.
* Supports lagged covariates via `lagvars`.
* Supports weight truncation at specified percentiles.
* Supports categorical predictor expansion and time-covariate interactions.
* Supports bootstrap standard errors.
