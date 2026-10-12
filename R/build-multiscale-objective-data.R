#' Prepare data for multiscale minimum-set optimization
#'
#' Combines the planning-unit hierarchy, costs, and feature amounts into the
#' data used for multiscale minimum-set optimization and solution evaluation.
#'
#' @param pu A planning-unit object. Its row order must match `maps`.
#' @param maps A hierarchy object returned by [build_h3_maps()].
#' @param cost_col Character name of the spatially additive planning-unit cost
#'   column.
#' @param feature_cols A character vector or named character vector identifying
#'   the multiscale feature columns, typically the output of 
#'   [multiscale_feature_cols()].
#' @param native_cols Optional character vector identifying native-resolution
#'   feature columns. When supplied, native and atom-based totals are compared
#'   as a diagnostic.
#' @param tol Numeric relative tolerance used for hierarchy reconciliation and
#'   diagnostic checks.
#'
#' @return A `multiscale_objective_data` object for use with
#'   [add_multiscale_min_set()] and [eval_multiscale_atom_solution()].
#'
#' @details
#' The prepared data combine cost and feature representation across nested
#' planning-unit resolutions. This allows planning units selected at different
#' resolutions to contribute toward the same feature representation while
#' shared spatial contributions are accounted for once.
#'
#' The planning-unit row order used to build this object must also be used for
#' the corresponding conservation problem and solution evaluation.
#'
#' @examples
#' # Small hierarchy: a parent and two of its seven children.
#' # Feature amounts below are already prepared and additive.
#' pu <- data.frame(
#'   h3_address = c("872a1072bffffff", "882a1072b1fffff", "882a1072b3fffff"),
#'   res = c(7L, 8L, 8L), area_km2 = c(14, 2, 2),
#'   msopt_A = c(3, 1, 2), msopt_B = c(6, 2, 1),
#'   selected = c(1L, 1L, 0L))
#' maps <- build_h3_maps(pu, res_levels = c(7L, 8L))
#' ms_data <- build_multiscale_objective_data(
#'   pu = pu, maps = maps,
#'   feature_cols = c(A = "msopt_A", B = "msopt_B"))
#'
#' result <- eval_multiscale_atom_solution(pu,
#'   solution_col = "selected", data = ms_data,
#'   targets = c(A = 2/3, B = 1/2))
#'
#' result$relative_held
#' result$target_met
#'
#' @family multiscale minimum-set
#' @export
build_multiscale_objective_data <- function(
    pu,
    maps,
    cost_col = "area_km2",
    feature_cols,
    native_cols = NULL,
    tol = 1e-6
) {
  atoms <- build_multiscale_atoms(
    pu = pu,
    maps = maps,
    cost_col = cost_col,
    tol = tol
  )
  
  feature_atoms <- build_multiscale_feature_atoms(
    pu = pu,
    atoms = atoms,
    feature_cols = feature_cols,
    maps = maps,
    tol = tol
  )
  
  df <- .ms_drop_geometry(pu)
  
  checks <- list(
    n_pu = nrow(df),
    n_atoms = length(atoms$atom_cost),
    cost_reconstruction_error = atoms$reconstruction_error,
    original_cost_reconstruction_error =
      atoms$original_reconstruction_error,
    n_reconciled_cost_rows = atoms$n_reconciled_rows,
    max_abs_cost_reconciliation =
      atoms$max_abs_reconciliation,
    max_relative_cost_reconciliation =
      atoms$max_relative_reconciliation,
    n_features = length(feature_atoms$feature_names),
    feature_totals = feature_atoms$feature_totals,
    feature_reconciliation =
      feature_atoms$reconciliation_summary
  )
  
  if (!is.null(native_cols)) {
    if (
      length(native_cols) !=
      length(feature_atoms$feature_names)
    ) {
      stop(
        "`native_cols` must have the same length as `feature_cols`."
      )
    }
    
    missing_native <- setdiff(native_cols, names(df))
    
    if (length(missing_native)) {
      stop(
        "Native columns missing from `pu`: ",
        paste(missing_native, collapse = ", ")
      )
    }
    
    native_totals <- vapply(
      native_cols,
      function(column) {
        sum(
          as.numeric(df[[column]]),
          na.rm = TRUE
        )
      },
      numeric(1)
    )
    
    names(native_totals) <-
      feature_atoms$feature_names
    
    atom_totals <- feature_atoms$feature_totals
    
    abs_error <- abs(native_totals - atom_totals)
    
    error_scale <- pmax(
      abs(native_totals),
      abs(atom_totals),
      .Machine$double.eps
    )
    
    rel_error <- abs_error / error_scale
    
    checks$native_totals <- native_totals
    checks$native_vs_atom_abs_error <- abs_error
    checks$native_vs_atom_rel_error <- rel_error
    
    if (any(rel_error > tol, na.rm = TRUE)) {
      warning(
        "Native feature totals and reconciled atom totals differ. ",
        "Inspect the returned object's `checks$native_vs_atom_rel_error` values."
      )
    }
  }
  
  structure(
    list(
      atoms = atoms,
      feature_atoms = feature_atoms,
      feature_cols = feature_atoms$feature_cols,
      native_cols = native_cols,
      cost_col = cost_col,
      tolerance = tol,
      checks = checks
    ),
    class = "multiscale_objective_data"
  )
}
