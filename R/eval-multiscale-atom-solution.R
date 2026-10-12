#' Evaluate multiscale feature representation
#'
#' Calculates feature representation and conservation cost for a selected set of
#' planning units in a multiscale conservation problem.
#'
#' @param pu Planning units containing the solution to evaluate, in the same row
#'   order used to build `data`.
#' @param solution_col Name of the column indicating selected planning units.
#' @param data A `multiscale_objective_data` object returned by
#'   [build_multiscale_objective_data()].
#' @param targets Optional relative representation targets for the features.
#'
#' @return A list containing:
#' \itemize{
#'   \item `feature_held`: amount of each feature represented by the solution.
#'   \item `feature_total`: total amount of each feature in the multiscale data.
#'   \item `relative_held`: proportion of each feature represented by the solution.
#'   \item `deduplicated_cost`: total multiscale cost of the solution.
#'   \item `n_selected_pu`: number of selected planning units.
#'   \item `n_selected_atoms`: number of spatial accounting units represented.
#'   \item `selected_atoms`: indexes of those spatial accounting units.
#' }
#'
#' When `targets` are supplied, the list also contains `target_relative`,
#' `target_absolute`, and `target_met`.
#'
#' @details
#' Feature representation is calculated across the planning-unit hierarchy so
#' that protection at different spatial resolutions contributes to the same
#' representation estimate. Spatial contributions shared by overlapping
#' planning units are counted once.
#'
#' For each feature, the function reports the total available amount, the amount
#' represented by the selected planning units, and the proportion represented.
#' When `targets` are supplied, it also reports the corresponding target amounts
#' and whether each target is met.
#' The function also reports the cost of the selected area using the same
#' multiscale accounting.
#'
#' @examples
#' pu <- data.frame(
#'   h3_address = c("872a1072bffffff","882a1072b1fffff","882a1072b3fffff"),
#'   res = c(7L, 8L, 8L),
#'   area_km2 = c(14, 2, 2),
#'   msopt_A = c(3, 1, 2),
#'   msopt_B = c(6, 2, 1),
#'   selected = c(1L, 1L, 0L))
#'
#' maps <- build_h3_maps(pu, res_levels = c(7L, 8L))
#'
#' ms_data <- build_multiscale_objective_data(pu, maps,
#' feature_cols = c(A = "msopt_A", B = "msopt_B"))
#'
#' result <- eval_multiscale_atom_solution(pu,
#'   solution_col = "selected", data = ms_data, targets = c(A = 2/3, B = 1/2))
#'
#' result$feature_held
#' result$relative_held
#' result$target_met
#' result$deduplicated_cost
#'
#' @family multiscale minimum-set
#' @export
eval_multiscale_atom_solution <- function(pu, solution_col, data, targets = NULL) {
  if (!inherits(data, "multiscale_objective_data")) {
    stop("`data` must come from build_multiscale_objective_data().")
  }
  df <- .ms_drop_geometry(pu)
  if (nrow(df) != length(data$atoms$atoms_by_pu)) {
    stop(
      "`pu` must contain the same number of planning units used to build `data`.")
  }
  if (!solution_col %in% names(df)) {
    stop("`solution_col` not found in `pu`.")
  }
  sel <- df[[solution_col]]
  if (!(is.numeric(sel) || is.logical(sel))) {
    stop("`solution_col` must contain numeric or logical 0/1 values.")
  }
  if (any(!is.na(sel) & !(sel %in% c(0, 1)))) {
    stop("`solution_col` must contain only 0/1 values or `NA`.")
  }
  sel <- as.integer(sel)
  sel[is.na(sel)] <- 0L
  idx <- which(sel == 1L)

  atom_on <- rep(FALSE, length(data$atoms$atom_cost))
  for (i in idx) {
    aa <- data$atoms$atoms_by_pu[[i]]
    if (length(aa)) atom_on[aa] <- TRUE
  }

  cost <- sum(data$atoms$atom_cost[atom_on], na.rm = TRUE)
  held <- if (any(atom_on)) Matrix::colSums(data$feature_atoms$feature_atom_matrix[atom_on, , drop = FALSE]) else rep(0, length(data$feature_atoms$feature_names))
  held <- as.numeric(held)
  names(held) <- data$feature_atoms$feature_names
  total <- data$feature_atoms$feature_totals
  names(total) <- data$feature_atoms$feature_names
  rel <- held / total
  rel[!is.finite(rel)] <- NA_real_

  out <- list(
    deduplicated_cost = cost,
    n_selected_pu = length(idx),
    n_selected_atoms = sum(atom_on),
    selected_atoms = which(atom_on),
    feature_held = held,
    feature_total = total,
    relative_held = rel
  )
  if (!is.null(targets)) {
    if (!is.numeric(targets)) {
      stop("`targets` must be numeric.")
    }
    
    n_features <- length(held)
    
    if (!(length(targets) %in% c(1L, n_features))) {
      stop("`targets` must contain either one value or one value per feature.")
    }
    
    if (length(targets) == 1L) {
      targets <- rep(targets, n_features)
      names(targets) <- names(held)
      
    } else if (is.null(names(targets))) {
      names(targets) <- names(held)
      
    } else {
      if (anyDuplicated(names(targets))) {
        stop("Named `targets` must have unique feature names.")
      }
      targets <- targets[names(held)]
    }
    
    if (
      anyNA(targets) ||
      any(!is.finite(targets)) ||
      any(targets < 0) ||
      any(targets > 1)
    ) {
      stop(
        "`targets` must contain finite proportions between zero and one, ",
        "with one target available for every feature."
      )
    }
    
    out$target_relative <- targets
    out$target_absolute <- targets * total
    out$target_met <- held >= out$target_absolute - 1e-8
  }
  out
}
