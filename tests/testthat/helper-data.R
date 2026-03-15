#' Create a basic panel dataset for testing
#'
#' @param n_subj Number of subjects.
#' @param n_vis Number of visits per subject.
#' @param seed Random seed.
#' @return A data.frame in long panel format.
make_panel_data <- function(n_subj = 20L, n_vis = 5L, seed = 42L) {
  set.seed(seed)
  N <- n_subj * n_vis
  data.frame(
    id = rep(seq_len(n_subj), each = n_vis),
    time = rep((seq_len(n_vis) - 1L) * 6, n_subj),
    severity = rnorm(N, 3, 1)
  )
}


#' Create a panel dataset with treatment for FIPTIW testing
#'
#' @param n_subj Number of subjects.
#' @param n_vis Number of visits per subject.
#' @param seed Random seed.
#' @return A data.frame with treatment and outcome.
make_treated_panel <- function(n_subj = 20L, n_vis = 5L, seed = 42L) {
  set.seed(seed)
  N <- n_subj * n_vis
  baseline_sev <- rep(rnorm(n_subj, 3, 1), each = n_vis)
  treated <- rep(rbinom(n_subj, 1, 0.5), each = n_vis)

  data.frame(
    id = rep(seq_len(n_subj), each = n_vis),
    time = rep((seq_len(n_vis) - 1L) * 6, n_subj),
    severity = baseline_sev + rnorm(N, 0, 0.3),
    baseline_sev = baseline_sev,
    treated = treated,
    outcome = 50 - 0.1 * rep((seq_len(n_vis) - 1L) * 6, n_subj) -
              baseline_sev + 0.5 * treated + rnorm(N, 0, 2)
  )
}
