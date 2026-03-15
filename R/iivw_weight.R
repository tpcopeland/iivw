#' Compute Inverse Intensity and Visit Weights
#'
#' Computes inverse intensity of visit weights (IIW) to correct for informative
#' visit processes in longitudinal clinic-based data. Optionally computes IPTW
#' for confounding by indication and their product (FIPTIW).
#'
#' Visit intensity is modeled via an Andersen-Gill recurrent-event Cox model on
#' the counting process of visits. Weights are the inverse of the estimated
#' conditional intensity ratio (Buzkova & Lumley 2007).
#'
#' @param data A data.frame in long panel format (one row per subject-visit).
#' @param id Character string; column name for subject identifier.
#' @param time Character string; column name for continuous visit time.
#' @param visit_cov Character vector; covariate column names for the visit
#'   intensity Cox model.
#' @param treat Character string or NULL; column name for binary (0/1)
#'   time-invariant treatment indicator. If specified, IPTW weights are
#'   computed and the default weight type becomes FIPTIW.
#' @param treat_cov Character vector or NULL; covariates for the propensity
#'   score model. Defaults to `visit_cov` if NULL and `treat` is specified.
#' @param wtype Character or NULL; weight type override. One of `"iivw"`,
#'   `"iptw"`, or `"fiptiw"`. If NULL (default), auto-detected: FIPTIW when
#'   `treat` is specified, IIW otherwise.
#' @param stabcov Character vector or NULL; covariates for the stabilized IIW
#'   numerator model. When specified, stabilized weights are computed as
#'   `exp(xb_stab - xb_full)` instead of `exp(-xb_full)`.
#' @param lagvars Character vector or NULL; time-varying covariates to lag by
#'   one visit within each subject. Creates columns named `{var}_lag1`.
#' @param entry Character string or NULL; column name for subject-specific
#'   study entry time. Default start time is 0 for all subjects.
#' @param truncate Numeric vector of length 2 or NULL; percentile bounds
#'   (0-100) for weight truncation (winsorization). E.g., `c(1, 99)`.
#' @param prefix Character string; prefix for generated weight column names.
#'   Default `"iivw_"`.
#'
#' @return A data.frame of class `"iivw_weights"` with added weight columns:
#'   * `{prefix}iw`: IIW component (if weight type includes IIW)
#'   * `{prefix}tw`: IPTW component (if weight type includes IPTW)
#'   * `{prefix}weight`: Final combined weight
#'
#'   The returned object has attributes storing metadata: `iivw_id`, `iivw_time`,
#'   `iivw_weight_type`, `iivw_weight_var`, `iivw_prefix`, `iivw_treat`,
#'   `iivw_diagnostics`.
#'
#' @references
#' Buzkova P, Lumley T (2007). "Longitudinal data analysis for generalized
#' linear models with follow-up dependent on outcome-related variables."
#' *Canadian Journal of Statistics*, 35(4), 485-500.
#'
#' Tompkins JD, Dubin JA, Wallace MP (2025). "Fully inverse probability of
#' treatment and intensity weighted estimators for non-randomized longitudinal
#' studies subject to irregular visit processes."
#'
#' @export
#' @examples
#' # Simulate panel data
#' set.seed(42)
#' n_subj <- 50
#' n_vis <- 4
#' d <- data.frame(
#'   id = rep(1:n_subj, each = n_vis),
#'   time = rep((0:(n_vis - 1)) * 6, n_subj),
#'   severity = rnorm(n_subj * n_vis, 3, 1)
#' )
#'
#' # IIW weights
#' result <- iivw_weight(d, id = "id", time = "time", visit_cov = "severity")
iivw_weight <- function(data, id, time, visit_cov,
                        treat = NULL, treat_cov = NULL,
                        wtype = NULL, stabcov = NULL,
                        lagvars = NULL, entry = NULL,
                        truncate = NULL, prefix = "iivw_") {

  # =========================================================================
  # Validate inputs
  # =========================================================================

  if (!is.character(id) || length(id) != 1L) {
    stop("`id` must be a single column name string", call. = FALSE)
  }
  if (!is.character(time) || length(time) != 1L) {
    stop("`time` must be a single column name string", call. = FALSE)
  }
  if (!is.character(visit_cov) || length(visit_cov) < 1L) {
    stop("`visit_cov` must be a character vector of column names", call. = FALSE)
  }

  .validate_iivw_data(data, id, time, visit_cov,
                       treat = treat, treat_cov = treat_cov,
                       stabcov = stabcov, lagvars = lagvars,
                       entry = entry)

  # =========================================================================
  # Determine weight type
  # =========================================================================

  if (is.null(wtype)) {
    wtype <- if (!is.null(treat)) "fiptiw" else "iivw"
  }
  wtype <- match.arg(wtype, c("iivw", "iptw", "fiptiw"))

  if (wtype %in% c("iptw", "fiptiw") && is.null(treat)) {
    stop("'", wtype, "' requires `treat` argument", call. = FALSE)
  }

  # Default treat_cov to visit_cov if not specified
  if (!is.null(treat) && is.null(treat_cov)) {
    treat_cov <- visit_cov
  }

  # Validate truncation
  if (!is.null(truncate)) {
    if (length(truncate) != 2L || !is.numeric(truncate)) {
      stop("`truncate` must be a numeric vector of length 2", call. = FALSE)
    }
    if (truncate[1] >= truncate[2]) {
      stop("`truncate` lower bound must be less than upper bound", call. = FALSE)
    }
    if (truncate[1] < 0 || truncate[2] > 100) {
      stop("`truncate` values must be between 0 and 100", call. = FALSE)
    }
  }

  # =========================================================================
  # Sort data
  # =========================================================================

  data <- data[order(data[[id]], data[[time]]), ]

  # =========================================================================
  # Lag variables (if requested)
  # =========================================================================

  if (!is.null(lagvars)) {
    for (v in lagvars) {
      lagname <- paste0(v, "_lag1")
      data[[lagname]] <- .lag_by_group(data[[v]], data[[id]])
    }
  }

  # Build full covariate list (original + lagged)
  visit_covars <- visit_cov
  if (!is.null(lagvars)) {
    visit_covars <- c(visit_covars, paste0(lagvars, "_lag1"))
  }

  # =========================================================================
  # Display header
  # =========================================================================

  wtype_display <- toupper(wtype)
  ids <- unique(data[[id]])
  n_ids <- length(ids)
  N <- nrow(data)

  cat("\n")
  cat(strrep("-", 70), "\n")
  cat("iivw_weight -", wtype_display, "Weight Computation\n")
  cat(strrep("-", 70), "\n\n")
  cat("ID variable:      ", id, "\n")
  cat("Time variable:    ", time, "\n")
  if (wtype %in% c("iivw", "fiptiw")) {
    cat("Visit covariates: ", paste(visit_covars, collapse = " "), "\n")
  }
  if (!is.null(treat)) {
    cat("Treatment:        ", treat, "\n")
  }
  if (!is.null(treat_cov)) {
    cat("Treatment covars: ", paste(treat_cov, collapse = " "), "\n")
  }
  cat("Weight type:      ", wtype_display, "\n")
  if (!is.null(truncate)) {
    cat("Truncation:       ", truncate[1], "th -", truncate[2], "th percentile\n")
  }
  cat("\n")

  # =========================================================================
  # IIW COMPONENT: Visit intensity model
  # =========================================================================

  iw_col <- paste0(prefix, "iw")

  if (wtype %in% c("iivw", "fiptiw")) {
    cat("Fitting visit intensity model (Andersen-Gill Cox)...\n")

    # Build counting process
    data <- .make_counting_process(data, id, time, entry)

    # Fit Andersen-Gill Cox model
    cox_formula <- as.formula(paste(
      "survival::Surv(.start, .stop, .event) ~",
      paste(.bt(visit_covars), collapse = " + "),
      "+ cluster(", .bt(id), ")"
    ))

    cat("  Visit model: coxph(Surv(start, stop, event) ~",
        paste(visit_covars, collapse = " + "), ")\n")

    # Suppress expected warning about zero-length intervals for first obs
    cox_fit <- suppressWarnings(survival::coxph(cox_formula, data = data))

    # Compute linear predictor for ALL rows manually (X %*% beta)
    # predict() may drop rows with invalid Surv intervals (start >= stop)
    cox_coefs <- coef(cox_fit)
    xb_full <- as.numeric(as.matrix(data[, names(cox_coefs), drop = FALSE]) %*% cox_coefs)

    # Compute IIW weights
    if (!is.null(stabcov)) {
      # Stabilized: fit numerator model with stabcov only
      stab_formula <- as.formula(paste(
        "survival::Surv(.start, .stop, .event) ~",
        paste(.bt(stabcov), collapse = " + "),
        "+ cluster(", .bt(id), ")"
      ))

      cat("  Stabilization model: coxph(... ~",
          paste(stabcov, collapse = " + "), ")\n")

      stab_fit <- suppressWarnings(survival::coxph(stab_formula, data = data))
      stab_coefs <- coef(stab_fit)
      xb_stab <- as.numeric(as.matrix(data[, names(stab_coefs), drop = FALSE]) %*% stab_coefs)
      data[[iw_col]] <- exp(xb_stab - xb_full)
    } else {
      # Unstabilized: exp(-xb)
      data[[iw_col]] <- exp(-xb_full)
    }

    # First observation per subject: weight = 1
    first_obs <- c(TRUE, data[[id]][-1] != data[[id]][-nrow(data)])
    data[[iw_col]][first_obs] <- 1

    # Clean up counting process columns
    data$.start <- NULL
    data$.stop <- NULL
    data$.event <- NULL
  }

  # =========================================================================
  # IPTW COMPONENT: Treatment model
  # =========================================================================

  tw_col <- paste0(prefix, "tw")

  if (wtype %in% c("iptw", "fiptiw")) {
    cat("Fitting treatment model (logistic)...\n")

    # Cross-sectional data: first obs per subject
    first_obs <- c(TRUE, data[[id]][-1] != data[[id]][-nrow(data)])
    subj_data <- data[first_obs, ]

    # Fit propensity score model on cross-sectional data
    ps_formula <- as.formula(paste(.bt(treat), "~", paste(.bt(treat_cov), collapse = " + ")))
    cat("  Treatment model: glm(", treat, "~",
        paste(treat_cov, collapse = " + "), ")\n")

    ps_fit <- glm(ps_formula, family = binomial(link = "logit"), data = subj_data)

    # Predict for all observations
    ps <- predict(ps_fit, newdata = data, type = "response")

    # Marginal treatment probability (from cross-sectional data)
    p_treat <- mean(subj_data[[treat]], na.rm = TRUE)

    # Stabilized IPTW
    data[[tw_col]] <- ifelse(
      data[[treat]] == 1,
      p_treat / ps,
      (1 - p_treat) / (1 - ps)
    )
  }

  # =========================================================================
  # COMBINE WEIGHTS
  # =========================================================================

  w_col <- paste0(prefix, "weight")

  if (wtype == "fiptiw") {
    data[[w_col]] <- data[[iw_col]] * data[[tw_col]]
  } else if (wtype == "iivw") {
    data[[w_col]] <- data[[iw_col]]
  } else if (wtype == "iptw") {
    data[[w_col]] <- data[[tw_col]]
  }

  # Warn if missing weights
  n_miss <- sum(is.na(data[[w_col]]))
  if (n_miss > 0L) {
    cat("Note:", n_miss, "observations have missing weights (missing covariates)\n")
  }

  # =========================================================================
  # TRUNCATION
  # =========================================================================

  n_truncated <- 0L
  if (!is.null(truncate)) {
    cat("Truncating weights at", truncate[1], "th and", truncate[2],
        "th percentiles...\n")

    trunc_result <- .truncate_weights(data[[w_col]], truncate[1], truncate[2])
    data[[w_col]] <- trunc_result$weights
    n_truncated <- trunc_result$n_truncated
    cat("  Truncated", n_truncated, "observations (",
        trunc_result$n_lo, "low,", trunc_result$n_hi, "high)\n")
  }

  # =========================================================================
  # DIAGNOSTICS
  # =========================================================================

  w <- data[[w_col]]
  w_valid <- w[!is.na(w)]
  w_mean <- mean(w_valid)
  w_sd <- sd(w_valid)
  w_min <- min(w_valid)
  w_max <- max(w_valid)
  w_q <- quantile(w_valid, probs = c(0.01, 0.50, 0.99))
  w_p1 <- w_q[1]
  w_p50 <- w_q[2]
  w_p99 <- w_q[3]
  ess <- .compute_ess(w_valid)

  cat("\nWeight distribution:\n")
  cat(sprintf("  Mean:     %9.4f\n", w_mean))
  cat(sprintf("  SD:       %9.4f\n", w_sd))
  cat(sprintf("  Min:      %9.4f\n", w_min))
  cat(sprintf("  Median:   %9.4f\n", w_p50))
  cat(sprintf("  Max:      %9.4f\n", w_max))
  cat(sprintf("  P1:       %9.4f\n", w_p1))
  cat(sprintf("  P99:      %9.4f\n", w_p99))
  cat("\n")
  cat(sprintf("Observations:          %9.0f\n", N))
  cat(sprintf("Subjects:              %9.0f\n", n_ids))
  cat(sprintf("Effective sample size: %9.1f (of %d)\n", ess, N))

  if (abs(w_mean - 1) > 0.2) {
    cat("\nNote: weight mean is", sprintf("%.3f", w_mean), "\n")
    cat("  Consider checking model specification or using truncation.\n")
  }

  # List created variables
  created <- w_col
  if (wtype %in% c("iivw", "fiptiw")) created <- c(iw_col, created)
  if (wtype %in% c("iptw", "fiptiw")) created <- c(tw_col, created)
  cat("\nVariables created:", paste(created, collapse = " "), "\n")
  cat("Next step: iivw_fit() to fit weighted outcome model\n")
  cat(strrep("-", 70), "\n")

  # =========================================================================
  # SET CLASS AND ATTRIBUTES
  # =========================================================================

  class(data) <- c("iivw_weights", setdiff(class(data), "iivw_weights"))
  attr(data, "iivw_id") <- id
  attr(data, "iivw_time") <- time
  attr(data, "iivw_weight_type") <- wtype
  attr(data, "iivw_weight_var") <- w_col
  attr(data, "iivw_prefix") <- prefix
  if (!is.null(treat)) attr(data, "iivw_treat") <- treat

  diag <- list(
    N = N,
    n_ids = n_ids,
    mean_weight = w_mean,
    sd_weight = w_sd,
    min_weight = w_min,
    max_weight = w_max,
    p1_weight = unname(w_p1),
    p99_weight = unname(w_p99),
    ess = ess,
    n_truncated = n_truncated
  )
  attr(data, "iivw_diagnostics") <- diag

  data
}
