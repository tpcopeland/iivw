# Internal utility functions for iivw
# Not exported; prefixed with . per package convention

#' Backtick-quote a variable name if it is non-syntactic
#'
#' @param x Character vector of variable names.
#' @return Character vector with non-syntactic names backtick-quoted.
#' @noRd
.bt <- function(x) {
  ifelse(make.names(x) == x, x, paste0("`", x, "`"))
}

#' Validate data for iivw_weight
#'
#' Checks: data.frame, required columns, >=2 visits per subject, unique id-time,
#' treatment is binary and time-invariant (if specified).
#'
#' @param data A data.frame.
#' @param id Character; column name for subject identifier.
#' @param time Character; column name for visit time.
#' @param visit_cov Character vector; covariate column names.
#' @param treat Character or NULL; treatment column name.
#' @param treat_cov Character vector or NULL; treatment model covariates.
#' @param stabcov Character vector or NULL; stabilization covariates.
#' @param lagvars Character vector or NULL; variables to lag.
#' @param entry Character or NULL; entry time column.
#'
#' @return Invisible NULL on success; stops with informative error otherwise.
#' @noRd
.validate_iivw_data <- function(data, id, time, visit_cov,
                                treat = NULL, treat_cov = NULL,
                                stabcov = NULL, lagvars = NULL,
                                entry = NULL) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame, not ", class(data)[1], call. = FALSE)
  }
  if (nrow(data) == 0L) {
    stop("no observations in `data`", call. = FALSE)
  }

  # Check required columns exist
  required <- c(id, time, visit_cov)
  if (!is.null(treat)) required <- c(required, treat)
  if (!is.null(treat_cov)) required <- c(required, treat_cov)
  if (!is.null(stabcov)) required <- c(required, stabcov)
  if (!is.null(lagvars)) required <- c(required, lagvars)
  if (!is.null(entry)) required <- c(required, entry)
  required <- unique(required)

  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0L) {
    stop("columns not found in data: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  # Check time is numeric
  if (!is.numeric(data[[time]])) {
    stop("`time` column '", time, "' must be numeric", call. = FALSE)
  }

  # Check covariates are numeric
  non_numeric <- visit_cov[!vapply(data[visit_cov], is.numeric, logical(1))]
  if (length(non_numeric) > 0L) {
    stop("visit_cov columns must be numeric: ",
         paste(non_numeric, collapse = ", "), call. = FALSE)
  }

  # Check >= 2 visits per subject
  visit_counts <- tapply(data[[id]], data[[id]], length)
  if (any(visit_counts < 2L)) {
    n_single <- sum(visit_counts < 2L)
    stop(n_single, " subject(s) have fewer than 2 visits; ",
         "IIW requires at least 2 visits per subject", call. = FALSE)
  }

  # Check unique id-time combinations
  id_time <- paste(data[[id]], data[[time]], sep = "|||")
  if (anyDuplicated(id_time)) {
    stop("duplicate id-time combinations found; ",
         "each subject-visit must be uniquely identified by id and time",
         call. = FALSE)
  }

  # Validate treatment variable
  if (!is.null(treat)) {
    treat_vals <- data[[treat]]
    if (!all(treat_vals[!is.na(treat_vals)] %in% c(0, 1))) {
      stop("`treat` column '", treat, "' must be binary (0/1)", call. = FALSE)
    }

    # Check time-invariant
    treat_by_id <- tapply(treat_vals, data[[id]], function(x) length(unique(x)))
    if (any(treat_by_id > 1L)) {
      stop("`treat` must be time-invariant within subjects", call. = FALSE)
    }

    # Check both groups present
    if (all(treat_vals == 0, na.rm = TRUE) || all(treat_vals == 1, na.rm = TRUE)) {
      stop("`treat` must have observations in both groups", call. = FALSE)
    }
  }

  invisible(NULL)
}


#' Compute effective sample size
#'
#' @param w Numeric vector of weights.
#' @return Scalar ESS = (sum(w))^2 / sum(w^2).
#' @noRd
.compute_ess <- function(w) {
  w <- w[!is.na(w)]
  (sum(w))^2 / sum(w^2)
}


#' Truncate weights at specified percentiles
#'
#' @param w Numeric vector of weights.
#' @param lower Lower percentile (0-100).
#' @param upper Upper percentile (0-100).
#' @return List with truncated weights and count of truncated observations.
#' @noRd
.truncate_weights <- function(w, lower, upper) {
  valid <- !is.na(w)
  bounds <- quantile(w[valid], probs = c(lower / 100, upper / 100))
  lo_val <- bounds[1]
  hi_val <- bounds[2]

  n_lo <- sum(w[valid] < lo_val)
  n_hi <- sum(w[valid] > hi_val)

  w[valid & w < lo_val] <- lo_val
  w[valid & w > hi_val] <- hi_val

  list(weights = w, n_truncated = n_lo + n_hi, n_lo = n_lo, n_hi = n_hi)
}


#' Lag a variable within groups
#'
#' @param x Numeric vector to lag.
#' @param id Vector of group identifiers (must be sorted by id then time).
#' @return Numeric vector with within-group lag-1; NA for first obs per group.
#' @noRd
.lag_by_group <- function(x, id) {
  n <- length(x)
  result <- rep(NA_real_, n)
  if (n < 2L) return(result)

  # Assumes data is already sorted by id, time
  same_id <- id[-1] == id[-n]
  result[c(FALSE, same_id)] <- x[c(same_id, FALSE)]
  result
}


#' Build counting process structure
#'
#' Creates start/stop/event columns for Andersen-Gill recurrent event model.
#'
#' @param data A data.frame sorted by id, time.
#' @param id Character; id column name.
#' @param time Character; time column name.
#' @param entry Character or NULL; entry time column name.
#' @return data.frame with added .start, .stop, .event columns.
#' @noRd
.make_counting_process <- function(data, id, time, entry = NULL) {
  ids <- data[[id]]
  times <- data[[time]]
  n <- nrow(data)

  # Start time: previous visit time (or entry/0 for first)
  start <- rep(0, n)
  if (n > 1L) {
    same_id <- ids[-1] == ids[-n]
    start[c(FALSE, same_id)] <- times[c(same_id, FALSE)]
  }

  # Apply entry times for first obs per subject
  if (!is.null(entry)) {
    first_obs <- c(TRUE, ids[-1] != ids[-n])
    start[first_obs] <- data[[entry]][first_obs]
  }

  data$.start <- start
  data$.stop <- times
  data$.event <- 1L
  data
}
