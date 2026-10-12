#' Add a multiscale minimum-set objective
#'
#' Sets the objective of a conservation planning problem to minimize cost while
#' meeting feature representation targets across nested H3 resolutions.
#'
#' @param x A single-zone `prioritizr` `ConservationProblem` object.
#' @param data A `multiscale_objective_data` object returned by
#'   [build_multiscale_objective_data()].
#'
#' @details
#' Feature representation is calculated across the planning-unit hierarchy, so 
#' selecting a planning unit at one resolution can contribute toward targets for
#' features represented at another resolution. Spatial contributions shared by
#' overlapping parent and descendant planning units are counted once.
#'
#' This accounting is implemented using non-overlapping spatial units, referred
#' to as atoms, for both cost and feature representation.
#'
#' The `data` object and conservation problem must use the same planning-unit row
#' order and matching feature names. The formulation supports single-zone
#' problems with binary planning-unit decisions and the compressed formulation.
#'
#' @section Mathematical formulation:
#' Let `x_i` indicate whether planning unit `i` is selected and let `z_a`
#' indicate whether atom `a` is activated. If `C_a` is the set of planning units
#' that cover atom `a`, atom activation is linked to planning-unit selection by
#' `z_a >= x_i` for each `i` in `C_a` and `z_a <= sum(x_i, i in C_a)`.
#' 
#' Planning-unit decisions are binary, with `x_i` in `{0, 1}`, while atom
#' activation variables satisfy `0 <= z_a <= 1`.
#'
#' The objective minimizes `sum(c_a * z_a)`, where `c_a` is the cost of atom
#' `a`. Feature representation targets require `sum(q_aj * z_a) >= T_j`,
#' where `q_aj` is the amount of feature `j` assigned to atom `a` and `T_j`
#' is its representation target.
#'
#' @return An updated `ConservationProblem` object.
#'
#' @examples
#' h3_vec <- c("872a1072bffffff","882a1072b1fffff","882a1072b3fffff")
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
#' ms_data <- build_multiscale_objective_data(
#'   pu, maps, feature_cols = c(A = "A", B = "B"))
#'
#' p <- prioritizr::problem(
#'   pu,features = c("A", "B"),cost_column = "area_km2")
#'
#' p <- add_multiscale_min_set_objective(p, ms_data)
#' p <- add_multiscale_relative_targets(p, c(A = 2/3, B = 1/2),ms_data)
#' p <- prioritizr::add_binary_decisions(p)
#'
#' @family objectives
#' @family multiscale minimum-set
#' @export
add_multiscale_min_set_objective <- function(x, data) {
  if (missing(x)) {
    stop("`x` is required.")
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
    stop("`add_multiscale_min_set_objective()` currently supports single-zone problems only.")
  }
  
  x$add_objective(
    R6::R6Class(
      "MultiscaleAtomMinimumSetObjective",
      inherit = prioritizr::Objective,
      public = list(
        name = "multiscale atom minimum set objective",
        has_weights = FALSE,
        has_targets = TRUE,
        data = list(ms_data = data),
        
        apply = function(x, y, weights) {
          if (!inherits(x, "OptimizationProblem")) {
            stop("Internal error: `x` must be an OptimizationProblem.")
          }
          
          if (!inherits(y, "ConservationProblem")) {
            stop("Internal error: `y` must be a ConservationProblem.")
          }
          
          if (!isTRUE(x$compressed_formulation())) {
            stop("Multiscale atom min-set currently requires the compressed formulation.")
          }
          
          ms <- self$get_data("ms_data")
          fam <- ms$feature_atoms$feature_atom_matrix
          
          # Reorder atom-feature matrix to match prioritizr feature order.
          y_features <- tryCatch(
            y$feature_names(),
            error = function(e) NULL
          )
          
          if (!is.null(y_features)) {
            pos <- match(y_features, ms$feature_atoms$feature_names)
            
            if (any(is.na(pos))) {
              stop(
                "Problem feature names do not match multiscale feature atoms. Missing: ",
                paste(y_features[is.na(pos)], collapse = ", ")
              )
            }
            
            fam <- fam[, pos, drop = FALSE]
          }
          
          invisible(
            .ms_apply_atom_backend(
              x,
              y$feature_targets(),
              ms$atoms$atom_cost,
              lapply(ms$atoms$cover_by_atom, as.integer),
              fam
            )
          )
        }
      )
    )$new()
  )
}
