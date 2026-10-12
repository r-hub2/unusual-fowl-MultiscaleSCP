#' Add multiscale connectivity penalties to a prioritizr problem
#'
#' @description
#' Adds penalties to a conservation planning problem to account for
#' cross-resolution connectivity between H3 planning units. Connections with
#' larger values are given greater importance when identifying solutions with
#' coordinated selection across resolutions.
#' Connectivity data can be supplied as a matrix, `Matrix::dgCMatrix`, data
#' frame in Marxan format, or four-dimensional array, following the formats
#' supported by [prioritizr::add_connectivity_penalties()].
#'
#' @inheritParams prioritizr::add_connectivity_penalties
#' @param data A symmetric connectivity object (matrix, `Matrix::dgCMatrix`,
#'   data frame in Marxan format, or 4D array) describing cross-scale links.
#' @param normalize Either `"none"` (use data as provided) or `"sym"` to apply
#'   symmetric degree normalization before penalization.
#'
#' @importFrom methods setGeneric setMethod
#' @importFrom Matrix rowSums Diagonal isSymmetric
#' @importFrom R6 R6Class
#' @importFrom assertthat is.number
#' @importFrom utils getFromNamespace
#'
#' @return The modified `ConservationProblem` object.
#' 
#' @export
methods::setGeneric(
  "add_multiscale_connectivity_penalties",
  signature = methods::signature("x", "penalty", "zones", "data", "normalize"),
  function(x, penalty, zones = diag(prioritizr::number_of_zones(x)), data,
           normalize = c("none","sym")) {
    normalize <- match.arg(normalize)
    .pz_assert_required(x)
    .pz_assert_required(penalty)
    .pz_assert_required(zones)
    .pz_assert_required(data)
    .pz_assert(
      .pz_is_conservation_problem(x),
      .pz_is_inherits(
        data, c("dgCMatrix", "data.frame", "matrix", "Matrix", "array")
      )
    )
    standardGeneric("add_multiscale_connectivity_penalties")
  }
)

#' @describeIn add_multiscale_connectivity_penalties matrix method
methods::setMethod(
  "add_multiscale_connectivity_penalties",
  methods::signature("ANY", "ANY", "ANY", "matrix", "character"),
  function(x, penalty, zones, data, normalize = c("none","sym")) {
    normalize <- match.arg(normalize)
    add_multiscale_connectivity_penalties(
      x, penalty, zones, .pz_as_Matrix(data, "dgCMatrix"), normalize
    )
  }
)

#' @describeIn add_multiscale_connectivity_penalties Matrix method
methods::setMethod(
  "add_multiscale_connectivity_penalties",
  methods::signature("ANY", "ANY", "ANY", "Matrix", "character"),
  function(x, penalty, zones, data, normalize = c("none","sym")) {
    normalize <- match.arg(normalize)
    add_multiscale_connectivity_penalties(
      x, penalty, zones, .pz_as_Matrix(data, "dgCMatrix"), normalize
    )
  }
)

#' @describeIn add_multiscale_connectivity_penalties data.frame (Marxan) method
methods::setMethod(
  "add_multiscale_connectivity_penalties",
  methods::signature("ANY", "ANY", "ANY", "data.frame", "character"),
  function(x, penalty, zones, data, normalize = c("none","sym")) {
    normalize <- match.arg(normalize)
    .pz_assert(
      .pz_is_conservation_problem(x),
      assertthat::is.number(penalty),
      .pz_all_finite(penalty),
      is.data.frame(data)
    )
    add_multiscale_connectivity_penalties(
      x, penalty, zones,
      prioritizr::marxan_connectivity_data_to_matrix(x, data, symmetric = TRUE),
      normalize
    )
  }
)

#' @describeIn add_multiscale_connectivity_penalties dgCMatrix method
methods::setMethod(
  "add_multiscale_connectivity_penalties",
  methods::signature("ANY", "ANY", "ANY", "dgCMatrix", "character"),
  function(x, penalty, zones, data, normalize = c("none","sym")) {
    normalize <- match.arg(normalize)
    .pz_assert(
      .pz_is_conservation_problem(x),
      assertthat::is.number(penalty),
      .pz_all_finite(penalty),
      .pz_is_matrix_ish(zones),
      nrow(zones) == ncol(zones),
      .pz_is_numeric_values(zones),
      .pz_all_finite(zones),
      .pz_is_numeric_values(data),
      .pz_all_finite(data),
      ncol(data) == nrow(data),
      max(zones) <= 1, min(zones) >= -1,
      .pz_number_of_total_units(x) == ncol(data),
      prioritizr::number_of_zones(x) == ncol(zones)
    )
    .pz_assert(
      Matrix::isSymmetric(data),
      msg = paste0(
        "{.arg data} does not contain symmetric connectivity values, ",
        "use {.fn add_asym_connectivity_penalties} instead."
      )
    )
    internal_add_multiscale_connectivity_penalties(
      x, penalty, as.matrix(zones), data, normalize
    )
  }
)

#' @describeIn add_multiscale_connectivity_penalties array (4D) method
methods::setMethod(
  "add_multiscale_connectivity_penalties",
  methods::signature("ANY", "ANY", "ANY", "array", "character"),
  function(x, penalty, zones, data, normalize = c("none","sym")) {
    normalize <- match.arg(normalize)
    .pz_assert(
      .pz_is_conservation_problem(x),
      assertthat::is.number(penalty),
      .pz_all_finite(penalty),
      is.null(zones),
      is.array(data),
      length(dim(data)) == 4,
      dim(data)[1] == .pz_number_of_total_units(x),
      dim(data)[2] == .pz_number_of_total_units(x),
      dim(data)[3] == prioritizr::number_of_zones(x),
      dim(data)[4] == prioritizr::number_of_zones(x),
      .pz_all_finite(data)
    )
    internal_add_multiscale_connectivity_penalties(x, penalty, zones, data, normalize)
  }
)

internal_add_multiscale_connectivity_penalties <- function(
    x, penalty, zones, data, normalize = c("none","sym")) {
  normalize <- match.arg(normalize)
  .pz_assert(
    .pz_is_conservation_problem(x),
    assertthat::is.number(penalty),
    .pz_all_finite(penalty),
    .internal = TRUE
  )
  x$add_penalty(
    R6::R6Class(
      "MultiscaleConnectivityPenalty",
      inherit = .pz_get("Penalty"),
      public = list(
        name = "multiscale connectivity penalties",
        data = list(
          penalty   = penalty,
          zones     = zones,
          data      = data,
          normalize = normalize
        ),
        apply = function(x, y) {
          .pz_assert(
            inherits(x, "OptimizationProblem"),
            inherits(y, "ConservationProblem"),
            .internal = TRUE
          )
          d         <- self$get_data("data")
          z         <- self$get_data("zones")
          norm_mode <- self$get_data("normalize")
          indices   <- y$planning_unit_indices()
          
          norm_sym <- function(M) {
            rs <- Matrix::rowSums(M)
            sc <- 1 / sqrt(pmax(rs, 1e-12))
            Dh <- Matrix::Diagonal(x = sc)
            Dh %*% M %*% Dh
          }
          maybe_norm <- function(M) if (norm_mode == "sym") norm_sym(M) else M
          
          m <- list()
          if (inherits(d, "dgCMatrix")) {
            d <- d[indices, indices, drop = FALSE]
            d <- maybe_norm(d)
            for (z1 in seq_len(ncol(z))) {
              m[[z1]] <- list()
              for (z2 in seq_len(nrow(z))) {
                m[[z1]][[z2]] <- d * z[z1, z2]
              }
            }
          } else if (inherits(d, "array")) {
            for (z1 in seq_len(dim(d)[3])) {
              m[[z1]] <- list()
              for (z2 in seq_len(dim(d)[4])) {
                sl <- .pz_as_Matrix(
                  d[indices, indices, z1, z2],
                  "dgCMatrix"
                )
                m[[z1]][[z2]] <- maybe_norm(sl)
              }
            }
          } else {
            cli::cli_abort(
              "Failed calculations for add_multiscale_connectivity_penalties.",
              .internal = TRUE
            )
          }
          
          m <- lapply(
            m,
            function(x) lapply(
              x,
              function(y) .pz_as_Matrix(Matrix::tril(y), "dgCMatrix")
            )
          )
          .pz_rcpp_apply_conn_penalty(
            x$ptr,
            self$get_data("penalty"),
            m
          )
          invisible(TRUE)
        }
      )
    )$new()
  )
}

# ---- internal hooks into prioritizr
.pz_get <- function(name) getFromNamespace(name, ns = "prioritizr")

.pz_assert_required          <- .pz_get("assert_required")
.pz_assert                   <- .pz_get("assert")
.pz_is_conservation_problem  <- .pz_get("is_conservation_problem")
.pz_is_inherits              <- .pz_get("is_inherits")
.pz_as_Matrix                <- .pz_get("as_Matrix")
.pz_all_finite               <- .pz_get("all_finite")
.pz_is_matrix_ish            <- .pz_get("is_matrix_ish")
.pz_is_numeric_values        <- .pz_get("is_numeric_values")
.pz_number_of_total_units    <- .pz_get("number_of_total_units")
.pz_rcpp_apply_conn_penalty  <- .pz_get("rcpp_apply_connectivity_penalties")