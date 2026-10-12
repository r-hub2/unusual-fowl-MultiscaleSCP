# Native-code registration for the multiscale atom backend.
#' @useDynLib MultiscaleSCP, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL

.ms_apply_atom_backend <- function(x, targets, atom_costs, cover_by_atom, fam) {
  model <- as.list(x)
  n <- model$number_of_planning_units
  if (!identical(model$number_of_zones, 1L) && model$number_of_zones != 1)
    stop("Only single-zone problems are supported.")
  if (!isTRUE(model$compressed_formulation))
    stop("Only compressed problems are supported.")
  if (length(model$vtype) != n || any(model$vtype != "B"))
    stop("The atom objective requires binary planning-unit decisions.")
  if (length(model$rhs) || length(model$sense) || length(model$row_ids))
    stop("Unexpected existing constraints before atom objective application.")
  if (any(!is.finite(model$A_x)) || any(model$A_x < 0))
    stop("Planning-unit feature amounts must be finite and non-negative.")
  if (!is.numeric(atom_costs) || !length(atom_costs) ||
      any(!is.finite(atom_costs)) || any(atom_costs < 0))
    stop("Atom costs must be finite and non-negative, with at least one atom.")
  if (length(cover_by_atom) != length(atom_costs) ||
      any(!vapply(cover_by_atom, function(i) {
        is.numeric(i) && length(i) > 0L && !anyNA(i) &&
          all(is.finite(i)) && all(i == as.integer(i)) &&
          all(i >= 1 & i <= n) && !anyDuplicated(i)
      }, logical(1)))) stop("Invalid atom-to-planning-unit mapping.")
  fam <- methods::as(fam, "dgCMatrix")
  methods::validObject(fam)
  if (nrow(fam) != length(atom_costs) ||
      ncol(fam) != model$number_of_features ||
      any(!is.finite(fam@x)) || any(fam@x < 0))
    stop("Invalid atom feature matrix.")
  if (length(targets$value) != ncol(fam) ||
      any(!is.finite(targets$value)) ||
      length(targets$sense) != ncol(fam) ||
      anyNA(targets$sense) ||
      any(!targets$sense %in% c(">=", "<=", "=")))
    stop("Supply one finite target per feature.")
  if (!is.null(targets$feature) &&
      !identical(as.integer(targets$feature), seq_len(ncol(fam))))
    stop("Targets must contain each feature once, in problem feature order.")
  updated <- rcpp_build_multiscale_atom_min_set_model(
    model, targets, atom_costs, lapply(cover_by_atom, as.integer), fam
  )
  # Preserve minimum-set identity for prioritizr presolve checks.
  updated$obj_id <- "min_set"
  
  # Rebuild the prioritizr optimization problem and update its stored pointer.
  rebuilt <- prioritizr::optimization_problem(updated)
  x$ptr <- rebuilt$ptr
  invisible(TRUE)
}
