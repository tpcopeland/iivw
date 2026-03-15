#' Fit Weighted Outcome Model
#'
#' Fits a weighted outcome model using weights from [iivw_weight()]. Supports
#' GEE-style estimation (GLM with clustered robust standard errors) or mixed
#' models with random intercepts. Provides flexible time specifications,
#' categorical predictor expansion, and time-covariate interactions.
#'
#' @param formula A formula specifying the outcome model (e.g., `y ~ treat`).
#'   Time terms are added automatically based on `timespec`.
#' @param data A data.frame with class `"iivw_weights"` (output from
#'   [iivw_weight()]), or a regular data.frame if `weight_var`, `id_var`, and
#'   `time_var` are specified.
#' @param model Character; `"gee"` (default) for GLM with clustered sandwich
#'   standard errors, or `"mixed"` for mixed model with random intercept.
#' @param family A family object or string for the GLM family. Default
#'   `gaussian()`.
#' @param link Character or NULL; link function override. Default uses the
#'   canonical link for the family.
#' @param timespec Character; time specification for the model. One of
#'   `"linear"` (default), `"quadratic"`, `"cubic"`, `"none"`, or `"ns(k)"`
#'   for natural splines with `k` degrees of freedom.
#' @param interaction Character vector or NULL; covariates to interact with
#'   time terms. Creates product terms for each covariate x each time variable.
#' @param categorical Character vector or NULL; variables in the formula to
#'   expand into dummy indicators.
#' @param basecat Named list or NULL; base categories for categorical variables.
#'   E.g., `list(arm = 0)`. Default uses lowest value.
#' @param cluster Character or NULL; column name for clustering variable.
#'   Default uses the id variable from `iivw_weight()` metadata.
#' @param bootstrap Integer; number of bootstrap replicates. 0 (default) uses
#'   sandwich standard errors only.
#' @param level Numeric; confidence level (default 95).
#' @param weight_var Character or NULL; weight column name (auto-detected from
#'   `iivw_weights` class if NULL).
#' @param id_var Character or NULL; id column name (auto-detected if NULL).
#' @param time_var Character or NULL; time column name (auto-detected if NULL).
#'
#' @return An object of class `"iivw_fit"` containing:
#' \describe{
#'   \item{coefficients}{Named vector of point estimates.}
#'   \item{se}{Named vector of standard errors.}
#'   \item{vcov}{Variance-covariance matrix.}
#'   \item{z}{Z-statistics.}
#'   \item{p}{P-values (two-sided).}
#'   \item{ci_lower}{Lower confidence bounds.}
#'   \item{ci_upper}{Upper confidence bounds.}
#'   \item{model_fit}{The underlying model object (glm or lmer).}
#'   \item{N}{Number of observations used.}
#'   \item{model}{Model type used.}
#'   \item{weight_type}{Weight type from iivw_weight.}
#'   \item{timespec}{Time specification used.}
#'   \item{weight_var}{Weight variable name.}
#'   \item{level}{Confidence level.}
#'   \item{family}{Family used.}
#'   \item{formula_used}{The actual formula fit.}
#'   \item{time_vars}{Time variable names added.}
#'   \item{ix_vars}{Interaction variable names created.}
#'   \item{cat_vars}{Categorical dummy variable names created.}
#' }
#'
#' @export
#' @examples
#' set.seed(42)
#' n_subj <- 50
#' n_vis <- 4
#' d <- data.frame(
#'   id = rep(1:n_subj, each = n_vis),
#'   time = rep((0:(n_vis - 1)) * 6, n_subj),
#'   severity = rnorm(n_subj * n_vis, 3, 1),
#'   outcome = rnorm(n_subj * n_vis, 50, 5)
#' )
#' d <- iivw_weight(d, id = "id", time = "time", visit_cov = "severity")
#' fit <- iivw_fit(outcome ~ severity, data = d)
iivw_fit <- function(formula, data,
                     model = c("gee", "mixed"),
                     family = gaussian(),
                     link = NULL,
                     timespec = c("linear", "quadratic", "cubic", "none"),
                     interaction = NULL,
                     categorical = NULL,
                     basecat = NULL,
                     cluster = NULL,
                     bootstrap = 0L,
                     level = 95,
                     weight_var = NULL,
                     id_var = NULL,
                     time_var = NULL) {

  # =========================================================================
  # Retrieve metadata
  # =========================================================================

  if (inherits(data, "iivw_weights")) {
    if (is.null(weight_var)) weight_var <- attr(data, "iivw_weight_var")
    if (is.null(id_var)) id_var <- attr(data, "iivw_id")
    if (is.null(time_var)) time_var <- attr(data, "iivw_time")
    weight_type <- attr(data, "iivw_weight_type")
    prefix <- attr(data, "iivw_prefix")
  } else {
    if (is.null(weight_var) || is.null(id_var) || is.null(time_var)) {
      stop("data is not from iivw_weight(); specify weight_var, id_var, time_var",
           call. = FALSE)
    }
    weight_type <- "custom"
    prefix <- "iivw_"
  }

  if (is.null(cluster)) cluster <- id_var

  # =========================================================================
  # Validate arguments
  # =========================================================================

  model <- match.arg(model)

  # Handle timespec: could be "ns(3)" etc.
  ns_df <- NULL
  if (is.character(timespec) && length(timespec) == 1L && grepl("^ns\\(\\d+\\)$", timespec)) {
    ns_df <- as.integer(sub("^ns\\((\\d+)\\)$", "\\1", timespec))
  } else if (length(timespec) > 1L || !(timespec[1] %in% c("linear", "quadratic", "cubic", "none"))) {
    # match.arg for the standard options
    timespec <- match.arg(timespec)
  }
  ts <- if (!is.null(ns_df)) paste0("ns(", ns_df, ")") else timespec[1]

  if (!is.null(interaction) && ts == "none") {
    stop("interaction requires time variables; not compatible with timespec='none'",
         call. = FALSE)
  }

  if (!is.null(basecat) && is.null(categorical)) {
    stop("basecat requires categorical", call. = FALSE)
  }

  if (model == "mixed" && !requireNamespace("lme4", quietly = TRUE)) {
    stop("mixed model requires the 'lme4' package", call. = FALSE)
  }

  bootstrap <- as.integer(bootstrap)

  # Handle family
  if (is.character(family)) {
    family <- get(family, mode = "function")()
  }
  if (!is.null(link)) {
    family <- do.call(family$family, list(link = link))
  }

  # =========================================================================
  # Parse formula
  # =========================================================================

  formula_terms <- terms(formula)
  depvar <- as.character(formula[[2]])
  indepvars <- attr(formula_terms, "term.labels")

  # =========================================================================
  # Mark sample
  # =========================================================================

  keep <- complete.cases(data[, c(depvar, indepvars, weight_var, time_var, cluster), drop = FALSE])
  N <- sum(keep)
  if (N == 0L) {
    stop("no complete observations", call. = FALSE)
  }

  # =========================================================================
  # Display header
  # =========================================================================

  wtype_display <- toupper(weight_type)
  cat("\n")
  cat(strrep("-", 70), "\n")
  cat("iivw_fit -", wtype_display, "Weighted Outcome Model\n")
  cat(strrep("-", 70), "\n\n")
  cat("Model type:       ", model, "\n")
  cat("Outcome:          ", depvar, "\n")
  cat("Predictors:       ", paste(indepvars, collapse = " "), "\n")
  cat("Time spec:        ", ts, "\n")
  if (!is.null(interaction)) {
    cat("Interactions:     ", paste(interaction, collapse = " "), "\n")
  }
  if (!is.null(categorical)) {
    cat("Categorical:      ", paste(categorical, collapse = " "), "\n")
  }
  if (model == "gee") {
    cat("Family:           ", family$family, "\n")
    cat("Link:             ", family$link, "\n")
    cat("Estimation:        GLM with clustered robust SEs\n")
  }
  cat("Weight var:       ", weight_var, "\n")
  cat("Cluster var:      ", cluster, "\n")
  if (bootstrap > 0L) {
    cat("Bootstrap reps:   ", bootstrap, "\n")
  }
  cat("\n")

  # =========================================================================
  # Build time specification variables
  # =========================================================================

  time_vars <- character(0)
  time_vars_created <- character(0)

  if (ts != "none") {
    time_vars <- time_var

    if (ts %in% c("quadratic", "cubic")) {
      tsq_name <- paste0(prefix, "time_sq")
      data[[tsq_name]] <- data[[time_var]]^2
      time_vars <- c(time_vars, tsq_name)
      time_vars_created <- c(time_vars_created, tsq_name)
    }
    if (ts == "cubic") {
      tcu_name <- paste0(prefix, "time_cu")
      data[[tcu_name]] <- data[[time_var]]^3
      time_vars <- c(time_vars, tcu_name)
      time_vars_created <- c(time_vars_created, tcu_name)
    }
    if (!is.null(ns_df)) {
      # Natural cubic splines using splines::ns()
      t_vals <- data[[time_var]][keep]
      ns_basis <- splines::ns(t_vals, df = ns_df)
      ns_knots <- attr(ns_basis, "knots")
      ns_bknots <- attr(ns_basis, "Boundary.knots")

      # Create basis for ALL rows using the same knots
      ns_all <- splines::ns(data[[time_var]], df = ns_df,
                             knots = ns_knots,
                             Boundary.knots = ns_bknots)

      time_vars <- character(0)
      for (j in seq_len(ns_df)) {
        tns_name <- paste0(prefix, "tns", j)
        data[[tns_name]] <- ns_all[, j]
        time_vars <- c(time_vars, tns_name)
        time_vars_created <- c(time_vars_created, tns_name)
      }
    }
  }

  # =========================================================================
  # Expand categorical variables
  # =========================================================================

  expanded_indepvars <- indepvars
  cat_vars_created <- character(0)
  expanded_interaction <- interaction

  if (!is.null(categorical)) {
    for (cvar in categorical) {
      if (!(cvar %in% indepvars)) {
        stop("'", cvar, "' in categorical not found in predictor variables",
             call. = FALSE)
      }

      vals <- sort(unique(data[[cvar]][keep & !is.na(data[[cvar]])]))
      if (length(vals) < 2L) {
        stop("'", cvar, "' in categorical has fewer than 2 unique values",
             call. = FALSE)
      }

      # Determine base category
      base_val <- vals[1]
      if (!is.null(basecat) && cvar %in% names(basecat)) {
        if (basecat[[cvar]] %in% vals) {
          base_val <- basecat[[cvar]]
        } else {
          cat("note: basecat for", cvar, "not found; using lowest value\n")
        }
      }

      non_base <- vals[vals != base_val]
      dummy_list <- character(0)

      for (lev in non_base) {
        vname <- paste0(prefix, "cat_", cvar, "_", lev)
        if (nchar(vname) > 32L) vname <- substr(vname, 1, 32)

        data[[vname]] <- as.integer(data[[cvar]] == lev)
        dummy_list <- c(dummy_list, vname)
        cat_vars_created <- c(cat_vars_created, vname)
      }

      # Replace original variable in predictor list (use first match only)
      idx <- which(expanded_indepvars == cvar)[1]
      if (!is.na(idx)) {
        expanded_indepvars <- c(
          expanded_indepvars[seq_len(idx - 1)],
          dummy_list,
          if (idx < length(expanded_indepvars)) expanded_indepvars[(idx + 1):length(expanded_indepvars)]
        )
      }

      # Replace in interaction if present
      if (!is.null(expanded_interaction) && cvar %in% expanded_interaction) {
        ix_idx <- which(expanded_interaction == cvar)[1]
        expanded_interaction <- c(
          expanded_interaction[seq_len(ix_idx - 1)],
          dummy_list,
          if (ix_idx < length(expanded_interaction)) expanded_interaction[(ix_idx + 1):length(expanded_interaction)]
        )
      }
    }
  }

  # =========================================================================
  # Build interaction variables
  # =========================================================================

  ix_vars <- character(0)
  ix_vars_created <- character(0)

  if (!is.null(expanded_interaction) && length(time_vars) > 0L) {
    for (ivar in expanded_interaction) {
      # Warn if not in predictors
      if (!(ivar %in% expanded_indepvars)) {
        cat("note:", ivar, "specified in interaction but not in predictors\n")
      }

      for (tvar in time_vars) {
        # Map time variable to suffix
        if (tvar == time_var) {
          suffix <- "time"
        } else if (tvar == paste0(prefix, "time_sq")) {
          suffix <- "tsq"
        } else if (tvar == paste0(prefix, "time_cu")) {
          suffix <- "tcu"
        } else {
          # Spline basis: strip prefix
          suffix <- sub(paste0("^", gsub("([.|()\\^{}+*?\\[\\]$])", "\\\\\\1", prefix)),
                        "", tvar)
        }

        # Covariate portion: strip cat prefix for cleaner names
        cat_prefix_str <- paste0(prefix, "cat_")
        if (startsWith(ivar, cat_prefix_str)) {
          ivar_portion <- sub(paste0("^", gsub("([.|()\\^{}+*?\\[\\]$])", "\\\\\\1", cat_prefix_str)),
                              "", ivar)
        } else {
          ivar_portion <- ivar
        }

        ix_name <- paste0(prefix, "ix_", ivar_portion, "_", suffix)
        if (nchar(ix_name) > 32L) ix_name <- substr(ix_name, 1, 32)

        data[[ix_name]] <- data[[ivar]] * data[[tvar]]
        ix_vars <- c(ix_vars, ix_name)
        ix_vars_created <- c(ix_vars_created, ix_name)
      }
    }
  }

  # =========================================================================
  # Build full covariate list and formula
  # =========================================================================

  all_covars <- c(expanded_indepvars, time_vars, ix_vars)
  fit_formula <- as.formula(paste(
    .bt(depvar), "~",
    paste(.bt(all_covars), collapse = " + ")
  ))

  # =========================================================================
  # Fit model
  # =========================================================================

  if (model == "gee") {
    cat("Fitting weighted GEE model...\n\n")

    if (bootstrap > 0L) {
      # Point estimate from full data (fit first for coefficient count)
      model_fit <- glm(fit_formula, data = data, family = family,
                        weights = data[[weight_var]], subset = keep)
      b <- coef(model_fit)
      n_coef <- length(b)

      # Subject-level bootstrap
      unique_ids <- unique(data[[id_var]][keep])
      # Pre-compute row indices per subject for efficiency
      row_idx_by_id <- split(seq_len(nrow(data)), data[[id_var]])

      bs_coefs <- matrix(NA_real_, nrow = bootstrap, ncol = n_coef)
      for (r in seq_len(bootstrap)) {
        sampled_ids <- sample(unique_ids, replace = TRUE)
        bs_rows <- unlist(row_idx_by_id[as.character(sampled_ids)],
                          use.names = FALSE)
        bs_data <- data[bs_rows, ]
        bs_fit <- tryCatch(
          glm(fit_formula, data = bs_data, family = family,
              weights = bs_data[[weight_var]]),
          error = function(e) NULL
        )
        if (!is.null(bs_fit)) {
          bs_coefs[r, ] <- coef(bs_fit)
        }
      }

      se <- apply(bs_coefs, 2, sd, na.rm = TRUE)
      V <- stats::cov(bs_coefs, use = "complete.obs")
    } else {
      model_fit <- glm(fit_formula, data = data, family = family,
                        weights = data[[weight_var]], subset = keep)

      V <- sandwich::vcovCL(model_fit, cluster = data[[cluster]][keep])
      b <- coef(model_fit)
      se <- sqrt(diag(V))
    }

  } else if (model == "mixed") {
    cat("Fitting weighted mixed model...\n\n")

    mixed_formula <- as.formula(paste(
      .bt(depvar), "~", paste(.bt(all_covars), collapse = " + "),
      "+ (1 |", .bt(id_var), ")"
    ))

    if (family$family == "gaussian" && family$link == "identity") {
      model_fit <- lme4::lmer(mixed_formula, data = data,
                               weights = data[[weight_var]],
                               subset = keep)
    } else {
      model_fit <- lme4::glmer(mixed_formula, data = data,
                                family = family,
                                weights = data[[weight_var]],
                                subset = keep)
    }

    b <- lme4::fixef(model_fit)
    V <- as.matrix(stats::vcov(model_fit))
    se <- sqrt(diag(V))
  }

  # =========================================================================
  # Compute test statistics
  # =========================================================================

  z_crit <- stats::qnorm((100 + level) / 200)
  z <- b / se
  p <- 2 * stats::pnorm(-abs(z))
  ci_lo <- b - z_crit * se
  ci_hi <- b + z_crit * se

  # =========================================================================
  # Display summary
  # =========================================================================

  cat(strrep("-", 70), "\n")

  # Show first predictor's effect
  first_pred <- expanded_indepvars[1]
  if (first_pred %in% names(b) && !is.na(se[first_pred]) && se[first_pred] > 0) {
    cat(wtype_display, "-weighted effect of ", first_pred, ":\n", sep = "")
    cat(sprintf("  Coefficient: %9.4f (SE: %7.4f)\n", b[first_pred], se[first_pred]))
    cat(sprintf("  %d%% CI: %9.4f - %9.4f\n", level, ci_lo[first_pred], ci_hi[first_pred]))
    cat(sprintf("  p-value:     %9.4f\n", p[first_pred]))
  }

  cat("\n")
  cat(strrep("-", 70), "\n")

  # =========================================================================
  # Build result
  # =========================================================================

  result <- list(
    coefficients = b,
    se = se,
    vcov = V,
    z = z,
    p = p,
    ci_lower = ci_lo,
    ci_upper = ci_hi,
    model_fit = model_fit,
    N = N,
    model = model,
    weight_type = weight_type,
    timespec = ts,
    weight_var = weight_var,
    level = level,
    family = family,
    formula_used = fit_formula,
    time_vars = time_vars,
    ix_vars = ix_vars,
    cat_vars = cat_vars_created
  )
  class(result) <- "iivw_fit"
  result
}
