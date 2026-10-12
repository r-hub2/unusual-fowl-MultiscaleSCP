#' Prepare a conservation feature for multiscale planning
#'
#' Extracts feature amounts from raster or vector data and prepares them for use
#' across a nested H3 planning-unit hierarchy. The returned planning-unit object
#' contains a native-resolution feature column for reporting and a multiscale
#' feature column used to represent the feature across resolutions.
#'
#' @param pu An `sf` planning-unit object. Its row order must match `maps`.
#' @param source A single-layer `SpatRaster`, `RasterLayer`, `RasterStack`,
#'   `RasterBrick`, `sf`, or `sfc` object containing the feature data.
#' @param native_res Integer H3 resolution at which the feature is represented
#'   for reporting.
#' @param maps A hierarchy object returned by [build_h3_maps()].
#' @param feature_name Optional feature name. When `NULL`, a name is inferred
#'   from `source`.
#' @param res_col Character name of the planning-unit resolution column.
#' @param opt_prefix Character prefix used for the multiscale feature column.
#' @param overwrite Logical; should existing output columns be replaced?
#' @param scale_raster_by_area Logical; should raster means be multiplied by
#'   planning-unit area?
#' @param value_col Optional character name of a value column in a vector
#'   source. When `NULL`, vector amounts are based on geometry measure alone.
#'
#' @details
#' Feature amounts are extracted from planning units at `native_res` or finer.
#' Coarser planning units receive amounts aggregated through the H3 hierarchy.
#' The native-resolution column contains the feature amount at `native_res`,
#' while the multiscale column records the amounts used to represent the feature
#' across resolutions. By default, these columns are named
#' `r<resolution>_<feature>` and `msopt_<feature>`.
#' 
#' Where finer planning units completely partition a parent, the parent amount
#' is obtained by summing child amounts. Where finer planning units cover only
#' part of a parent, the feature amount in the uncovered parent area is combined
#' with the child amounts.
#'
#' For raster data, values are extracted with
#' [exactextractr::exact_extract()] using a mean weighted by raster-cell
#' coverage fraction. By default, this mean is multiplied by the extracted
#' geometry area in square kilometres to produce an area-scaled feature amount.
#' With `scale_raster_by_area = FALSE`, the coverage-weighted mean is used
#' directly.
#'
#' Polygon features contribute area in square kilometres, line features
#' contribute length in kilometres, and point features contribute counts.
#' Single and multipart geometries can be combined within one geometry
#' family; each point in a MULTIPOINT is counted separately. Sources mixing
#' polygons, lines, or points, and GEOMETRYCOLLECTION features, must be
#' separated before loading.
#' Supplying `value_col` weights these amounts by the corresponding source
#' values.
#' 
#' @return The input `pu` object with a native-resolution feature column, a
#'   multiscale feature column, and metadata describing the feature.
#'
#' @examples
#' # Create one H3 parent and one of its children.
#' parent <- "872a1072bffffff" 
#' child <- h3jsr::get_children(parent, 8L, simple = TRUE)
#' child <- as.character(unlist(child, use.names = FALSE))[1]
#'
#' h3_vec <- c(parent, child)
#'
#' pu <- h3jsr::cell_to_polygon(h3_vec, simple = FALSE)
#' pu <- pu[match(h3_vec, pu$h3_address), ]
#' pu$res <- c(7L, 8L)
#'
#' maps <- build_h3_maps(pu,res_levels = c(7L, 8L))
#'
#' # Feature occupies the child planning unit.
#' feature <- sf::st_sf(geometry = sf::st_geometry(pu[2, ]))
#'
#' pu <- load_multiscale_feature(pu = pu,source = feature,
#'   native_res = 8L,maps = maps,feature_name = "habitat")
#'
#' sf::st_drop_geometry(pu)[,c("h3_address", "res", "r8_habitat", "msopt_habitat")]
#' multiscale_feature_cols(pu)
#'
#' @family multiscale minimum-set
#' @export
load_multiscale_feature <- function(
    pu,
    source,
    native_res,
    maps,
    feature_name = NULL,
    res_col = "res",
    opt_prefix = "msopt_",
    overwrite = FALSE,
    scale_raster_by_area = TRUE,
    value_col = NULL
) {
  df <- .ms_drop_geometry(pu)
  n <- nrow(df)
  .ms_assert_maps(maps, n)
  
  if (!inherits(pu, "sf")) {
    stop("`pu` must be an sf object.")
  }
  
  if (!res_col %in% names(df)) {
    stop("`res_col` not found in `pu`.")
  }
  
  native_res <- as.integer(native_res)[1]
  if (!native_res %in% as.integer(maps$res_levels)) {
    stop("`native_res` must be one of maps$res_levels.")
  }
  
  nm <- .ms_source_name(source, feature_name)
  native_col <- paste0("r", native_res, "_", nm)
  opt_col <- paste0(opt_prefix, nm)
  
  if (!overwrite && any(c(native_col, opt_col) %in% names(pu))) {
    stop("Output column already exists. Use `overwrite = TRUE` to replace it.")
  }
  
  source_type <- .ms_source_type(source)

  partition_info <- .ms_complete_child_partitions(maps)
  
  direct_pu <- .ms_direct_feature_geometries(
    pu = pu,
    maps = maps,
    native_res = native_res,
    complete_partition = partition_info$complete
  )
  
  direct_amount <- numeric(n)
  
  if (nrow(direct_pu) > 0L) {
    extracted <- .ms_extract_source_amount(
      pu = direct_pu,
      source = source,
      scale_raster_by_area = scale_raster_by_area,
      value_col = value_col
    )
    
    direct_amount[direct_pu$.ms_source_row] <- extracted
  }
  
  direct_amount[!is.finite(direct_amount)] <- 0
  direct_amount[direct_amount < 0] <- 0
  
  opt_amount <- direct_amount
  res_vec <- as.integer(maps$res_vec)
  res_levels <- sort(unique(as.integer(maps$res_levels)))
  children <- maps$children_by_row
  
  for (r in rev(res_levels)) {
    idx <- which(res_vec == r)
    
    for (i in idx) {
      ch <- children[[i]]
      
      if (length(ch)) {
        opt_amount[i] <- direct_amount[i] + sum(opt_amount[ch], na.rm = TRUE)
      }
    }
  }
  
  opt_amount[!is.finite(opt_amount)] <- 0
  opt_amount[opt_amount < 0] <- 0
  
  native <- rep(0, n)
  native_rows <- which(res_vec == native_res)
  native[native_rows] <- opt_amount[native_rows]
  
  pu[[native_col]] <- native
  pu[[opt_col]] <- opt_amount
  
  meta <- attr(pu, "multiscale_features")
  if (is.null(meta)) meta <- list()
  meta[[nm]] <- list(
    feature = nm,
    native_res = native_res,
    native_col = native_col,
    opt_col = opt_col,
    source_type = source_type,
    loading_rule = "finest_available_then_aggregate_upward_h3_complete_partitions",
    units = if (inherits(source, c("SpatRaster", "RasterLayer", "RasterStack", "RasterBrick")) && scale_raster_by_area) {
      "raster_mean_times_geometry_area_km2"
    } else if (source_type == "raster") {
      "raster_mean"
    } else {
      "geometry_amount"
    }
  )
  attr(pu, "multiscale_features") <- meta
  
  pu
}

#' Retrieve feature columns for multiscale optimization
#'
#' Returns the multiscale feature columns created by
#' [load_multiscale_feature()] in a named vector that can be passed directly to
#' the `feature_cols` argument of [build_multiscale_objective_data()].
#'
#' @param pu A planning-unit object containing a `multiscale_features`
#'   attribute.
#' @param names_as Character; use `"native"` to name the returned vector with
#'   native feature columns, or `"feature"` to use feature names.
#'
#' @return A named character vector mapping native feature columns or feature
#'   names to their corresponding multiscale feature columns.
#'
#' @details
#' With `names_as = "native"`, the vector is named using the native feature
#' columns, such as `r8_habitat`. With `names_as = "feature"`, it is named using
#' the feature names supplied to [load_multiscale_feature()].
#'
#' @examples
#' parent <- "872a1072bffffff" # Create one H3 parent and one of its children.
#' child <- h3jsr::get_children(parent, 8L, simple = TRUE)
#' child <- as.character(unlist(child, use.names = FALSE))[1]
#'
#' h3_vec <- c(parent, child)
#'
#' pu <- h3jsr::cell_to_polygon(h3_vec, simple = FALSE)
#' pu <- pu[match(h3_vec, pu$h3_address), ]
#' pu$res <- c(7L, 8L)
#'
#' maps <- build_h3_maps(pu,res_levels = c(7L, 8L))
#'
#' feature <- sf::st_sf(geometry = sf::st_geometry(pu[2, ]))
#'
#' pu <- load_multiscale_feature(pu = pu,source = feature,
#'   native_res = 8L,maps = maps,feature_name = "habitat")
#'
#' # Map native feature columns to optimization columns.
#' multiscale_feature_cols(pu)
#'
#' @family multiscale minimum-set
#' @export
multiscale_feature_cols <- function(pu, names_as = c("native", "feature")) {
  names_as <- match.arg(names_as)
  meta <- attr(pu, "multiscale_features")
  if (is.null(meta) || !length(meta)) stop("No `multiscale_features` metadata found.")
  vals <- vapply(meta, function(x) x$opt_col, character(1))
  if (names_as == "native") names(vals) <- vapply(meta, function(x) x$native_col, character(1))
  vals
}
