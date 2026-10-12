#' Build non-overlapping cost atoms for nested planning units
#'
#' Decomposes additive planning-unit costs into non-overlapping spatial
#' contributions ("atoms") so costs can be evaluated consistently across nested H3
#' resolutions.
#'
#' @param pu A planning-unit object. Its row order must match `maps`.
#' @param maps A hierarchy object returned by [build_h3_maps()].
#' @param cost_col Character name of a finite, non-negative, spatially additive
#'   planning-unit cost column.
#' @param tol Numeric relative tolerance used when reconciling parent and child
#'   values.
#'
#' @return An object of class `multiscale_atoms` containing atom costs,
#'   planning-unit--atom incidence structures, reconciled planning-unit costs,
#'   and reconciliation diagnostics.
#'
#' @details
#' Atoms are mutually exclusive accounting units defined from the planning-unit
#' hierarchy. Terminal planning units define atoms, and parents with incomplete
#' child partitions can contribute a residual atom. Complete H3 child partitions
#' are represented by their child atoms.
#'
#' `cost_col` must represent a finite, non-negative, spatially additive quantity,
#' such as area. A complete child partition is identified when directly mapped
#' children occur at one common H3 resolution and match the expected H3 child set.
#'
#' @examples
#' pu <- data.frame(
#'   h3_address = c("872a1072bffffff", "882a1072b1fffff", "882a1072b3fffff"),
#'   res = c(7L, 8L, 8L),
#'   area_km2 = c(14, 2, 2),
#'   A = c(3, 1, 2),
#'   B = c(6, 2, 1)
#' )
#' maps <- build_h3_maps(pu)
#' atoms <- build_multiscale_atoms(pu, maps)
#' atoms$atom_cost # Parent residual 10; two children 2 each
#' atoms$atoms_by_pu # Parent covers all three atoms
#'
#' @family multiscale minimum-set
#' @export
build_multiscale_atoms <- function(
    pu,
    maps,
    cost_col = "area_km2",
    tol = 1e-6
) {
  df <- .ms_drop_geometry(pu)
  n <- nrow(df)
  
  .ms_assert_maps(maps, n)
  
  if (!cost_col %in% names(df)) {
    stop("`cost_col` was not found in `pu`.")
  }
  
  cost <- suppressWarnings(as.numeric(df[[cost_col]]))
  
  if (
    length(cost) != n ||
    any(!is.finite(cost)) ||
    any(cost < 0)
  ) {
    stop("`cost_col` must be finite and non-negative for all planning units.")
  }
  
  parent <- maps$nearest_parent_row_of
  
  partition_info <- .ms_complete_child_partitions(maps)
  
  reconciliation <- .ms_reconcile_hierarchy_values(
    values = cost,
    maps = maps,
    complete_partition = partition_info$complete,
    tol = tol
  )
  
  reconciled_cost <- reconciliation$reconciled
  exclusive_cost <- reconciliation$exclusive
  
  # Define atom rows from the hierarchy: terminal rows and parents without a
  # complete child partition contribute atoms.
  atom_rows <- which(!partition_info$complete)
  atom_cost <- exclusive_cost[atom_rows]
  
  n_atoms <- length(atom_rows)
  
  cover_by_atom <- vector("list", n_atoms)
  atoms_by_pu <- vector("list", n)
  
  ii <- integer(0)
  jj <- integer(0)
  
  for (a in seq_along(atom_rows)) {
    atom_row <- atom_rows[a]
    current_row <- atom_row
    cover <- integer(0)
    
    while (!is.na(current_row)) {
      cover <- c(cover, current_row)
      
      atoms_by_pu[[current_row]] <- c(
        atoms_by_pu[[current_row]],
        a
      )
      
      ii <- c(ii, a)
      jj <- c(jj, current_row)
      
      current_row <- parent[current_row]
    }
    
    cover_by_atom[[a]] <- as.integer(cover)
  }
  
  atoms_by_pu <- lapply(
    atoms_by_pu,
    function(x) sort(unique(as.integer(x)))
  )
  
  atom_by_pu <- Matrix::sparseMatrix(
    i = ii,
    j = jj,
    x = 1,
    dims = c(n_atoms, n)
  )
  
  reconstructed <- if (n_atoms > 0L) {
    as.numeric(Matrix::crossprod(atom_by_pu, atom_cost))
  } else {
    rep(0, n)
  }
  
  reconstruction_error <- max(
    abs(reconstructed - reconciled_cost),
    na.rm = TRUE
  )
  
  original_reconstruction_error <- max(
    abs(reconstructed - cost),
    na.rm = TRUE
  )
  
  adjusted <- reconciliation$relative_adjustment > tol
  
  structure(
    list(
      cost_col = cost_col,
      
      row_cost = cost,
      exclusive_cost = exclusive_cost,
      atom_rows = atom_rows,
      atom_cost = atom_cost,
      atom_by_pu = atom_by_pu,
      cover_by_atom = cover_by_atom,
      atoms_by_pu = atoms_by_pu,
      reconstruction_error = reconstruction_error,
      
      reconciled_row_cost = reconciled_cost,
      cost_adjustment = reconciliation$adjustment,
      cost_relative_adjustment = reconciliation$relative_adjustment,
      reconciliation_action = reconciliation$action,
      complete_partition = partition_info$complete,
      observed_child_count = partition_info$observed_n,
      expected_child_count = partition_info$expected_n,
      original_reconstruction_error = original_reconstruction_error,
      n_reconciled_rows = sum(adjusted, na.rm = TRUE),
      max_abs_reconciliation = max(
        abs(reconciliation$adjustment),
        na.rm = TRUE
      ),
      max_relative_reconciliation = max(
        reconciliation$relative_adjustment,
        na.rm = TRUE
      )
    ),
    class = "multiscale_atoms"
  )
}

#' Build multiscale feature accounting for nested planning units
#'
#' Distributes hierarchy-consistent feature amounts among non-overlapping
#' spatial units ("atoms") so feature representation can be evaluated consistently across
#' nested H3 resolutions.
#'
#' @param pu A planning-unit object containing the multiscale feature columns.
#' @param atoms A `multiscale_atoms` object returned by
#'   [build_multiscale_atoms()].
#' @param feature_cols A character vector or named character vector identifying
#'   multiscale feature columns. When named, names define feature
#'   names used by the optimization problem.
#' @param maps A hierarchy object returned by [build_h3_maps()].
#' @param tol Numeric relative tolerance used when reconciling parent and child
#'   values.
#'
#' @return An object of class `multiscale_feature_atoms` containing the sparse
#'   atom-by-feature matrix, atom-based feature totals, and reconciliation
#'   diagnostics.
#'
#' @details
#' Rows of `feature_atom_matrix` are atoms and columns are features. Each value
#' represents the feature amount assigned to that atom, and `feature_totals`
#' contains the corresponding column totals. Feature names come from
#' `names(feature_cols)`, or from the column names when `feature_cols` is unnamed.
#'
#' Feature amounts should be defined so that the amount in a parent planning unit
#' can be decomposed among its non-overlapping child and remainder areas.
#' Non-finite feature values are treated as zero and negative values cause an
#' error. Reconciliation diagnostics are returned with the result.
#'
#' @examples
#' pu <- data.frame(
#'   h3_address = c("872a1072bffffff", "882a1072b1fffff", "882a1072b3fffff"),
#'   res = c(7L, 8L, 8L),
#'   area_km2 = c(14, 2, 2),
#'   A = c(3, 1, 2),
#'   B = c(6, 2, 1)
#' )
#' maps <- build_h3_maps(pu)
#' atoms <- build_multiscale_atoms(pu, maps)
#' features <- build_multiscale_feature_atoms(
#'   pu, atoms, feature_cols = c(A = "A", B = "B"), maps = maps
#' )
#' as.matrix(features$feature_atom_matrix)
#' features$feature_totals # A = 3, B = 6
#'
#' @family multiscale minimum-set
#' @export
build_multiscale_feature_atoms <- function(
    pu,
    atoms,
    feature_cols,
    maps,
    tol = 1e-6
) {
  if (!inherits(atoms, "multiscale_atoms")) {
    stop("`atoms` must come from build_multiscale_atoms().")
  }
  
  df <- .ms_drop_geometry(pu)
  n <- nrow(df)
  
  .ms_assert_maps(maps, n)
  
  if (is.list(feature_cols) && !is.character(feature_cols)) {
    feature_cols <- unlist(feature_cols, use.names = TRUE)
  }
  
  feature_names <- names(feature_cols)
  feature_cols <- as.character(unname(feature_cols))
  
  if (is.null(feature_names)) {
    feature_names <- feature_cols
  } else {
    empty_names <- is.na(feature_names) | !nzchar(feature_names)
    feature_names[empty_names] <- feature_cols[empty_names]
  }
  if (!length(feature_cols)) {
    stop("`feature_cols` must contain at least one feature column.")
  }
  
  if (anyDuplicated(feature_names)) {
    stop("Feature names in `feature_cols` must be unique.")
  }
  missing_cols <- setdiff(feature_cols, names(df))
  
  if (length(missing_cols)) {
    stop(
      "Feature columns missing from `pu`: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  children <- maps$children_by_row
  parent_rows <- which(lengths(children) > 0L)
  
  atom_rows <- atoms$atom_rows
  
  if (!length(atom_rows)) {
    stop("No footprint atoms were created.")
  }
  
  feature_matrix <- matrix(
    0,
    nrow = length(atom_rows),
    ncol = length(feature_cols),
    dimnames = list(NULL, feature_names)
  )
  
  dropped_positive <- stats::setNames(
    numeric(length(feature_cols)),
    feature_names
  )
  
  reconciliation_summary <- vector(
    "list",
    length(feature_cols)
  )
  
  # Force zero feature residuals for parents without a corresponding atom
  force_zero_residual <- rep(FALSE, n)
  force_zero_residual[parent_rows] <-
    !parent_rows %in% atoms$atom_rows
  
  for (j in seq_along(feature_cols)) {
    opt_col <- feature_cols[j]
    feature_name <- feature_names[j]
    
    values <- suppressWarnings(as.numeric(df[[opt_col]]))
    
    values[is.na(values) | !is.finite(values)] <- 0
    
    if (any(values < 0)) {
      stop("Feature column has negative values: ", opt_col)
    }
    
    reconciliation <- .ms_reconcile_hierarchy_values(
      values = values,
      maps = maps,
      complete_partition = atoms$complete_partition,
      tol = tol,
      force_zero_residual = force_zero_residual
    )
    
    exclusive_feature <- reconciliation$exclusive
    
    feature_matrix[, j] <- exclusive_feature[atom_rows]
    
    non_atom_rows <- setdiff(
      which(exclusive_feature > 0),
      atom_rows
    )
    
    if (length(non_atom_rows)) {
      dropped_positive[feature_name] <- sum(
        exclusive_feature[non_atom_rows],
        na.rm = TRUE
      )
    }
    
    adjusted <-
      reconciliation$relative_adjustment > tol
    
    reconciliation_summary[[j]] <- data.frame(
      feature = feature_name,
      n_reconciled_rows = sum(adjusted, na.rm = TRUE),
      max_abs_reconciliation = max(
        abs(reconciliation$adjustment),
        na.rm = TRUE
      ),
      max_relative_reconciliation = max(
        reconciliation$relative_adjustment,
        na.rm = TRUE
      ),
      sum_abs_reconciliation = sum(
        abs(reconciliation$adjustment),
        na.rm = TRUE
      ),
      dropped_positive = dropped_positive[feature_name],
      stringsAsFactors = FALSE
    )
  }
  
  feature_atom_matrix <- Matrix::Matrix(
    feature_matrix,
    sparse = TRUE
  )
  
  feature_totals <- Matrix::colSums(feature_atom_matrix)
  
  feature_totals <- stats::setNames(
    as.numeric(feature_totals),
    feature_names
  )
  
  reconciliation_summary <- do.call(
    rbind,
    reconciliation_summary
  )
  
  if (any(dropped_positive > 0)) {
    warning(
      "Some positive feature amounts occurred outside the footprint atom ",
      "structure. Inspect the returned object's `dropped_positive` values."
    )
  }
  
  structure(
    list(
      feature_names = feature_names,
      feature_cols = stats::setNames(
        feature_cols,
        feature_names
      ),
      feature_atom_matrix = feature_atom_matrix,
      feature_totals = feature_totals,
      dropped_positive = dropped_positive,
      reconciliation_summary = reconciliation_summary
    ),
    class = "multiscale_feature_atoms"
  )
}
