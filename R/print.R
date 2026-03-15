#' Print iivw_weights Object
#'
#' Prints a summary of the weight diagnostics from [iivw_weight()].
#'
#' @param x An `iivw_weights` object.
#' @param ... Additional arguments (ignored).
#'
#' @return `x` invisibly.
#' @export
print.iivw_weights <- function(x, ...) {
  diag <- attr(x, "iivw_diagnostics")
  wtype <- toupper(attr(x, "iivw_weight_type"))

  cat("iivw_weights:", wtype, "weights\n")
  cat("  Observations:", diag$N, "\n")
  cat("  Subjects:    ", diag$n_ids, "\n")
  cat("  Weight mean: ", sprintf("%.4f", diag$mean_weight), "\n")
  cat("  Weight SD:   ", sprintf("%.4f", diag$sd_weight), "\n")
  cat("  ESS:         ", sprintf("%.1f", diag$ess), "\n")
  cat("  Weight var:  ", attr(x, "iivw_weight_var"), "\n")
  cat(sprintf("  [%d x %d data.frame]\n", nrow(x), ncol(x)))
  invisible(x)
}


#' Print iivw_fit Object
#'
#' Prints a coefficient table for a fitted IIVW model.
#'
#' @param x An `iivw_fit` object.
#' @param ... Additional arguments (ignored).
#'
#' @return `x` invisibly.
#' @export
print.iivw_fit <- function(x, ...) {
  wtype <- toupper(x$weight_type)

  cat("\n", wtype, "-weighted outcome model (", x$model, ")\n", sep = "")
  cat("Family:", x$family$family, ", Link:", x$family$link, "\n")
  cat("Time spec:", x$timespec, "\n")
  cat("N =", x$N, "\n\n")

  # Coefficient table
  tab <- data.frame(
    Estimate = x$coefficients,
    SE = x$se,
    z = x$z,
    `p-value` = x$p,
    check.names = FALSE
  )
  tab[[paste0(x$level, "% CI lower")]] <- x$ci_lower
  tab[[paste0(x$level, "% CI upper")]] <- x$ci_upper

  print(round(tab, 4))
  cat("\n")
  invisible(x)
}


#' Summarize iivw_fit Object
#'
#' @param object An `iivw_fit` object.
#' @param ... Additional arguments (ignored).
#'
#' @return The `iivw_fit` object invisibly.
#' @export
summary.iivw_fit <- function(object, ...) {
  print.iivw_fit(object, ...)
  if (length(object$time_vars) > 0L) {
    cat("Time variables: ", paste(object$time_vars, collapse = ", "), "\n")
  }
  if (length(object$ix_vars) > 0L) {
    cat("Interaction variables: ", paste(object$ix_vars, collapse = ", "), "\n")
  }
  if (length(object$cat_vars) > 0L) {
    cat("Categorical dummies: ", paste(object$cat_vars, collapse = ", "), "\n")
  }
  cat("Formula:", deparse(object$formula_used), "\n")
  invisible(object)
}


#' Extract Coefficients from iivw_fit
#'
#' @param object An `iivw_fit` object.
#' @param ... Additional arguments (ignored).
#'
#' @return Named numeric vector of coefficients.
#' @export
coef.iivw_fit <- function(object, ...) {
  object$coefficients
}


#' Extract Variance-Covariance Matrix from iivw_fit
#'
#' @param object An `iivw_fit` object.
#' @param ... Additional arguments (ignored).
#'
#' @return Variance-covariance matrix.
#' @export
vcov.iivw_fit <- function(object, ...) {
  object$vcov
}
