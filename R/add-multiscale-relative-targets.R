#' Add multiscale relative targets
#'
#' Adds feature representation targets expressed as proportions between zero and
#' one of the total amount of each feature in the multiscale planning data.
#'
#' @param x A single-zone `prioritizr` `ConservationProblem` object.
#' @param targets A numeric value or vector of proportions between zero and one.
#'   Supply either one value for all features or one value per problem feature.
#'   A named vector is matched to the problem feature names.
#' @param data A `multiscale_objective_data` object returned by
#'   [build_multiscale_objective_data()].
#'
#' @return An updated `ConservationProblem` object with the targets added.
#'
#' @details
#' Relative targets specify the minimum proportion of each feature that should be
#' represented by the solution. Feature totals are calculated across the nested
#' planning-unit hierarchy, allowing representation secured at different
#' resolutions to contribute toward the same target.
#'
#' A single target value is applied to every feature. An unnamed vector follows
#' the feature order in the problem, and a named vector is matched by feature
#' name and must provide a value for every feature.
#'
#' @examples
#' h3_vec <- c("872a1072bffffff","882a1072b1fffff","882a1072b3fffff")
#'
#' pu <- h3jsr::cell_to_polygon(h3_vec, simple = FALSE)
#' pu <- pu[match(h3_vec, pu$h3_address), ]
#'
#' pu$res <- c(7L, 8L, 8L)
#' pu$area_km2 <- c(14, 2, 2)
#' pu$A <- c(3, 1, 2)
#' pu$B <- c(6, 2, 1)
#'
#' maps <- build_h3_maps(pu)
#'
#' ms_data <- build_multiscale_objective_data(pu, maps, feature_cols = c(A = "A", B = "B"))
#'
#' p <- prioritizr::problem(pu,features = c("A", "B"),cost_column = "area_km2")
#' p <- add_multiscale_min_set_objective(p,ms_data)
#' p <- add_multiscale_relative_targets(p,c(B = 1/2, A = 2/3),ms_data)
#'
#' @family targets
#' @family multiscale minimum-set
#' @export
add_multiscale_relative_targets <- function(x, targets, data) {
  
  if (missing(x)) {
    stop("`x` is required.")
  }
  
  if (missing(targets)) {
    stop("`targets` is required.")
  }
  
  if (missing(data)) {
    stop("`data` is required.")
  }
  
  if (!inherits(x, "ConservationProblem")) {
    stop("`x` must be a prioritizr ConservationProblem object.")
  }
  
  if (!inherits(data, "multiscale_objective_data")) {
    stop("`data` must come from build_multiscale_objective_data().")
  }
  
  if (prioritizr::number_of_zones(x) != 1L) {
    stop(
      "`add_multiscale_relative_targets()` currently supports ",
      "single-zone problems only."
    )
  }
  
  if (!is.numeric(targets)) {
    stop("`targets` must be numeric.")
  }
  
  problem_features <- x$feature_names()
  n_features <- length(problem_features)
  
  if (!(length(targets) %in% c(1L, n_features))) {
    stop(
      "`targets` must contain either one value or one value per feature."
    )
  }
  
  if (length(targets) == 1L) {
    targets <- rep(targets, n_features)
    names(targets) <- problem_features
    
  } else if (is.null(names(targets))) {
    names(targets) <- problem_features
    
  } else {
    if (anyDuplicated(names(targets))) {
      stop("Named `targets` must have unique feature names.")
    }
    targets <- targets[problem_features]
  }
  
  if (
    anyNA(targets) ||
    any(!is.finite(targets)) ||
    any(targets < 0) ||
    any(targets > 1)
  ) {
    stop(
      "`targets` must contain finite proportions between zero and one, ",
      "with one target available for every problem feature."
    )
  }
  
  atom_totals <- data$feature_atoms$feature_totals[
    problem_features
  ]
  
  if (anyNA(atom_totals)) {
    stop(
      "The features in `x` do not match the features in `data`."
    )
  }
  
  absolute_targets <- as.numeric(targets * atom_totals)
  
  # Suppress prioritizr's native-column infeasibility warning because target
  # feasibility is evaluated against atom-based feature totals.
  withCallingHandlers(
    prioritizr::add_absolute_targets(
      x = x,
      targets = absolute_targets
    ),
    warning = function(w) {
      message <- conditionMessage(w)
      
      if (
        grepl("targets", message, fixed = TRUE) &&
        grepl("infeasible values", message, fixed = TRUE)
      ) {
        invokeRestart("muffleWarning")
      }
    }
  )
}
