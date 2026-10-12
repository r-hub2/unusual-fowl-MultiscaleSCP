.ms_drop_geometry <- function(x) {
  tryCatch(sf::st_drop_geometry(x), error = function(e) x)
}

.ms_safe_name <- function(x) {
  x <- as.character(x)[1]
  x <- gsub("\\.[A-Za-z0-9]+$", "", basename(x))
  x <- make.names(x)
  gsub("\\.", "_", x)
}

.ms_source_name <- function(source, feature_name = NULL) {
  if (!is.null(feature_name)) return(.ms_safe_name(feature_name))
  nms <- tryCatch(names(source), error = function(e) NULL)
  if (!is.null(nms) && length(nms) && !is.na(nms[1]) && nzchar(nms[1])) {
    return(.ms_safe_name(nms[1]))
  }
  stop("Could not infer feature name from `source`; please supply `feature_name`.")
}

.ms_assert_maps <- function(maps, n = NULL) {
  req <- c("h3_vec", "res_vec", "res_levels", "nearest_parent_row_of", "children_by_row")
  miss <- setdiff(req, names(maps))
  if (length(miss)) stop("`maps` is missing: ", paste(miss, collapse = ", "))
  if (!is.null(n) && length(maps$res_vec) != n) {
    stop("`maps` length does not match nrow(pu).")
  }
  invisible(TRUE)
}

.ms_descendants <- function(i, children_by_row) {
  frontier <- children_by_row[[i]]
  if (!length(frontier)) return(integer(0))
  out <- integer(0)
  while (length(frontier)) {
    out <- c(out, frontier)
    nxt <- integer(0)
    for (j in frontier) {
      ch <- children_by_row[[j]]
      if (length(ch)) nxt <- c(nxt, ch)
    }
    frontier <- setdiff(unique(nxt), out)
  }
  unique(out)
}

.ms_terminal_rows_up_to_res <- function(maps, native_res) {
  res_vec <- as.integer(maps$res_vec)
  children <- maps$children_by_row
  n <- length(res_vec)
  out <- rep(FALSE, n)
  for (i in seq_len(n)) {
    if (res_vec[i] > native_res) next
    desc <- .ms_descendants(i, children)
    has_desc_up_to_native <- length(desc) && any(res_vec[desc] <= native_res & res_vec[desc] > res_vec[i])
    out[i] <- !has_desc_up_to_native
  }
  which(out)
}

.ms_area_km2 <- function(pu) {
  if (!inherits(pu, "sf")) {
    stop("Planning-unit area requires an `sf` object.")
  }
  
  as.numeric(
    units::set_units(
      sf::st_area(pu),
      "km^2",
      mode = "standard"
    )
  )
}

.ms_extract_raster_amount <- function(pu, source, scale_by_area = TRUE) {
  if (!inherits(pu, "sf")) stop("Raster extraction requires `pu` to be an sf object.")
  if (!requireNamespace("exactextractr", quietly = TRUE)) {
    stop("Package `exactextractr` is required for raster feature loading.")
  }
  if (inherits(source, "SpatRaster")) {
    if (terra::nlyr(source) != 1L) stop("`load_multiscale_feature()` expects a single raster layer.")
    r <- source
  } else if (inherits(source, c("RasterLayer", "RasterStack", "RasterBrick"))) {
    if (raster::nlayers(source) != 1L) stop("`load_multiscale_feature()` expects a single raster layer.")
    r <- source
  } else {
    stop("Unsupported raster class.")
  }
  v <- exactextractr::exact_extract(r, pu, fun = "mean")
  if (is.data.frame(v)) v <- v[[1]]
  v <- as.numeric(v)
  v[!is.finite(v)] <- 0
  if (isTRUE(scale_by_area)) v <- v * .ms_area_km2(pu)
  v
}

.ms_extract_sf_amount <- function(pu, source, value_col = NULL) {
  if (!inherits(pu, "sf")) stop("Vector feature loading requires `pu` to be an sf object.")
  if (!inherits(source, "sf")) {
    if (inherits(source, "sfc")) source <- sf::st_sf(value = 1, geometry = source)
    else stop("Unsupported vector source. Use an sf/sfc object.")
  }
  source <- sf::st_make_valid(sf::st_transform(source, sf::st_crs(pu)))
  pu2 <- sf::st_make_valid(pu)
  pu2$.ms_row_id <- seq_len(nrow(pu2))
  if (!is.null(value_col) && !value_col %in% names(source)) stop("`value_col` not found in source.")
  if (is.null(value_col)) source$.ms_value <- 1 else source$.ms_value <- as.numeric(source[[value_col]])
  source$.ms_value[!is.finite(source$.ms_value)] <- 0

  inter <- suppressWarnings(sf::st_intersection(pu2[, ".ms_row_id"], source[, ".ms_value"]))
  out <- numeric(nrow(pu2))
  if (!nrow(inter)) return(out)

  gtypes <- unique(as.character(sf::st_geometry_type(source)))
  if (all(gtypes %in% c("POLYGON", "MULTIPOLYGON"))) {
    measure <- as.numeric(units::set_units(
      sf::st_area(inter), "km^2", mode = "standard"
    ))
  } else if (all(gtypes %in% c("LINESTRING", "MULTILINESTRING"))) {
    measure <- as.numeric(units::set_units(
      sf::st_length(inter), "km", mode = "standard"
    ))
  } else if (all(gtypes %in% c("POINT", "MULTIPOINT"))) {
    geometries <- sf::st_geometry(inter)
    empty <- sf::st_is_empty(geometries)
    measure <- vapply(seq_along(geometries), function(i) {
      if (empty[i]) return(0)
      as.numeric(nrow(sf::st_coordinates(geometries[[i]])))
    }, numeric(1))
  } else {
    stop(
      "`source` must contain only polygons, only lines, or only points ",
      "(single or multipart). Mixed dimensions and GEOMETRYCOLLECTION ",
      "features are not supported; separate them before loading."
    )
  }
  amount <- measure * inter$.ms_value
  s <- tapply(amount, inter$.ms_row_id, sum, na.rm = TRUE)
  out[as.integer(names(s))] <- as.numeric(s)
  out[!is.finite(out)] <- 0
  out
}

# Internal hierarchy reconciliation helpers
.ms_complete_child_partitions <- function(maps) {
  .ms_assert_maps(maps)
  
  h3_vec <- as.character(maps$h3_vec)
  res_vec <- as.integer(maps$res_vec)
  children <- maps$children_by_row
  
  n <- length(h3_vec)
  
  complete <- rep(FALSE, n)
  observed_n <- lengths(children)
  expected_n <- rep(NA_integer_, n)
  
  parent_rows <- which(observed_n > 0L)
  
  for (i in parent_rows) {
    ch <- children[[i]]
    
    # A complete partition must use one common child resolution.
    child_res <- unique(res_vec[ch])
    
    if (
      length(child_res) != 1L ||
      !is.finite(child_res) ||
      child_res <= res_vec[i]
    ) {
      next
    }
    
    expected <- tryCatch(
      h3jsr::get_children(
        h3_address = h3_vec[i],
        res = child_res,
        simple = TRUE
      ),
      error = function(e) character(0)
    )
    
    # get_children() normally returns a list for one parent.
    if (is.list(expected)) {
      expected <- expected[[1]]
    }
    
    expected <- as.character(expected)
    expected_n[i] <- length(expected)
    
    complete[i] <-
      length(expected) > 0L &&
      setequal(h3_vec[ch], expected)
  }
  
  list(
    complete = complete,
    observed_n = observed_n,
    expected_n = expected_n
  )
}


.ms_reconcile_hierarchy_values <- function(
    values,
    maps,
    complete_partition,
    tol = 1e-6,
    force_zero_residual = NULL
) {
  values <- suppressWarnings(as.numeric(values))
  
  if (
    length(values) != length(maps$res_vec) ||
    any(!is.finite(values)) ||
    any(values < 0)
  ) {
    stop("Hierarchy values must be finite and non-negative.")
  }
  
  if (
    !is.numeric(tol) ||
    length(tol) != 1L ||
    !is.finite(tol) ||
    tol < 0
  ) {
    stop("`tol` must be a single finite non-negative value.")
  }
  
  n <- length(values)
  
  if (is.null(force_zero_residual)) {
    force_zero_residual <- rep(FALSE, n)
  }
  
  if (length(force_zero_residual) != n) {
    stop("`force_zero_residual` must have one value per planning unit.")
  }
  
  force_zero_residual <- as.logical(force_zero_residual)
  force_zero_residual[is.na(force_zero_residual)] <- FALSE
  
  children <- maps$children_by_row
  res_vec <- as.integer(maps$res_vec)
  res_levels <- sort(unique(res_vec))
  
  original <- values
  reconciled <- values
  exclusive <- values
  
  action <- rep("terminal", n)
  child_total <- rep(0, n)
  
  # Finest to coarsest so that child values are already reconciled.
  for (r in rev(res_levels)) {
    idx <- which(res_vec == r)
    
    for (i in idx) {
      ch <- children[[i]]
      
      if (!length(ch)) {
        next
      }
      
      ch_sum <- sum(reconciled[ch], na.rm = TRUE)
      child_total[i] <- ch_sum
      
      difference <- original[i] - ch_sum
      
      scale <- max(
        abs(original[i]),
        abs(ch_sum),
        .Machine$double.eps
      )
      
      approximately_equal <-
        abs(difference) <= tol * scale
      
      if (isTRUE(complete_partition[i])) {
        # All expected H3 children are present.
        reconciled[i] <- ch_sum
        exclusive[i] <- 0
        action[i] <- "complete_child_partition"
        
      } else if (isTRUE(force_zero_residual[i])) {
        # Set the parent residual to zero when no corresponding parent atom exists.
        reconciled[i] <- ch_sum
        exclusive[i] <- 0
        action[i] <- "matched_to_atom_structure"
        
      } else if (approximately_equal) {
        # Absorb small scale-dependent numerical discrepancies.
        reconciled[i] <- ch_sum
        exclusive[i] <- 0
        action[i] <- "within_tolerance"
        
      } else if (difference > 0) {
        # Children cover only part of the parent.
        reconciled[i] <- original[i]
        exclusive[i] <- difference
        action[i] <- "retained_parent_residual"
        
      } else {
        # Children exceed the parent value.
        # Preserve all child amounts and raise the parent total.
        reconciled[i] <- ch_sum
        exclusive[i] <- 0
        action[i] <- "raised_to_child_total"
      }
    }
  }
  
  # Clamp residual negative values introduced by numerical precision.
  exclusive[exclusive < 0] <- 0
  
  adjustment <- reconciled - original
  
  adjustment_scale <- pmax(
    abs(original),
    abs(reconciled),
    .Machine$double.eps
  )
  
  relative_adjustment <- abs(adjustment) / adjustment_scale
  
  list(
    original = original,
    reconciled = reconciled,
    exclusive = exclusive,
    child_total = child_total,
    adjustment = adjustment,
    relative_adjustment = relative_adjustment,
    action = action
  )
}


# ---- 1. feature loading ------------------------------------------------------

.ms_direct_feature_geometries <- function(
    pu,
    maps,
    native_res,
    complete_partition = NULL,
    tol_area = 0
) {
  if (!inherits(pu, "sf")) {
    stop("Feature loading with hierarchy residuals requires `pu` to be an sf object.")
  }
  
  n <- nrow(pu)
  .ms_assert_maps(maps, n)
  
  res_vec <- as.integer(maps$res_vec)
  children <- maps$children_by_row
  geom <- sf::st_geometry(pu)
  crs <- sf::st_crs(pu)

  if (is.null(complete_partition)) {
    complete_partition <- rep(FALSE, n)
  }
  if (length(complete_partition) != n) {
    stop("`complete_partition` must have one value per planning unit.")
  }
  complete_partition <- as.logical(complete_partition)
  complete_partition[is.na(complete_partition)] <- FALSE
  
  # Direct extraction is done only where the feature is defined at this
  # resolution or a finer one. Coarser parents are filled by aggregation.
  candidate_rows <- which(res_vec >= native_res)
  
  out_rows <- integer(0)
  out_geoms <- list()
  
  for (i in candidate_rows) {
    g <- geom[i]
    ch <- children[[i]]
    
    # Complete H3 child partitions are identified from child indexes. No parent
    # geometry residual is created for these partitions.
    if (length(ch) && isTRUE(complete_partition[i])) {
      next
    }

    # For partial child partitions, extract feature amounts from the uncovered
    # parent remainder and aggregate child amounts upward.
    if (length(ch)) {
      child_union <- suppressWarnings(sf::st_union(geom[ch]))
      g <- suppressWarnings(sf::st_difference(g, child_union))
    }
    
    if (!length(g) || all(sf::st_is_empty(g))) {
      next
    }
    
    g <- suppressWarnings(sf::st_make_valid(g))
    g <- suppressWarnings(sf::st_collection_extract(g, "POLYGON"))
    
    if (!length(g) || all(sf::st_is_empty(g))) {
      next
    }
    
    a <- suppressWarnings(as.numeric(sum(sf::st_area(g), na.rm = TRUE)))
    if (!is.finite(a) || a <= tol_area) {
      next
    }
    
    g <- suppressWarnings(sf::st_union(g))
    
    out_rows <- c(out_rows, i)
    out_geoms[[length(out_geoms) + 1L]] <- sf::st_geometry(g)[[1]]
  }
  
  if (!length(out_rows)) {
    out <- pu[0, , drop = FALSE]
    out$.ms_source_row <- integer(0)
    return(out[, ".ms_source_row", drop = FALSE])
  }
  
  sf::st_sf(
    .ms_source_row = out_rows,
    geometry = sf::st_sfc(out_geoms, crs = crs)
  )
}

.ms_extract_source_amount <- function(
    pu,
    source,
    scale_raster_by_area = TRUE,
    value_col = NULL
) {
  if (inherits(source, c("SpatRaster", "RasterLayer", "RasterStack", "RasterBrick"))) {
    return(
      .ms_extract_raster_amount(
        pu,
        source,
        scale_by_area = scale_raster_by_area
      )
    )
  }
  
  if (inherits(source, c("sf", "sfc"))) {
    return(
      .ms_extract_sf_amount(
        pu,
        source,
        value_col = value_col
      )
    )
  }
  
  stop("`source` must be a single-layer raster or an sf/sfc vector object.")
}

.ms_source_type <- function(source) {
  if (inherits(source, c("SpatRaster", "RasterLayer", "RasterStack", "RasterBrick"))) {
    return("raster")
  }
  if (inherits(source, c("sf", "sfc"))) {
    return("vector")
  }
  stop("`source` must be a single-layer raster or an sf/sfc vector object.")
}

# Internal hooks into prioritizr
.pz_get <- function(name) getFromNamespace(name, ns = "prioritizr")
.pz_assert_required <- function(...) .pz_get("assert_required")(...)
.pz_assert <- function(...) .pz_get("assert")(...)
.pz_is_conservation_problem <- function(x) .pz_get("is_conservation_problem")(x)
