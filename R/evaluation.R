#' Count cross-resolution overlaps in a selection
#'
#' @description
#' For a given 0/1 selection column, counts how many planning units are
#' simultaneously selected at a finer resolution and at their ancestor cell
#' at a coarser resolution (using the hierarchy in `maps`).
#'
#' For each pair of resolutions `r_low < r_high`, it reports the number of
#' co-selected ancestor–descendant pairs, plus a total across all pairs.
#'
#' @param s An `sf` or data frame of planning units. Must contain a resolution
#'   column `res` and a 0/1 selection column referenced by `flag_col`.
#' @param flag_col Name of the 0/1 selection column to analyse
#'   (e.g. `"selected_by_resolution"`).
#' @param maps A hierarchy list as returned by [build_h3_maps()], containing
#'   at least `res_levels` and `nearest_parent_row_of`.
#'
#' @return A named integer vector with one entry per resolution pair
#'   (e.g. `"5-6"`, `"6-7"`) and a `"total"` element.
#'   
#' @examples
#' parent <- "872a1072bffffff"
#' kids   <- c("882a1072b1fffff", "882a1072b3fffff")
#'
#' s <- data.frame(
#'   h3_address = c(parent, kids),
#'   res = c(7L, 8L, 8L),
#'   selected_by_resolution_final = c(1L, 1L, 0L)
#' )
#' maps <- build_h3_maps(s, res_levels = c(7L, 8L))
#'
#' compute_overlaps_by_resolution(s, "selected_by_resolution_final", maps)
#'
#' @export
compute_overlaps_by_resolution <- function(
    s,
    flag_col = "selected_by_resolution_final",  # or "selected_by_resolution_final"
    maps
){
  df   <- sf::st_drop_geometry(s)
  sel  <- dplyr::coalesce(suppressWarnings(as.integer(df[[flag_col]])), 0L)
  resv <- as.integer(df$res)
  
  res_levels <- sort(unique(as.integer(maps$res_levels)))
  parent_row <- maps$nearest_parent_row_of
  L <- length(res_levels)
  if (L < 2L) return(c(total = 0L))
  
  # Vectorized ancestor lookup for many indices at once:
  # climb via nearest-parent links until parent is at/below r_low.
  ancestor_at_vec <- function(ii, r_low) {
    a <- ii
    # move off self to parent if currently finer than r_low
    move <- resv[a] > r_low
    if (any(move)) a[move] <- parent_row[a[move]]
    
    # climb while still above r_low
    while (any(!is.na(a) & resv[a] > r_low)) {
      w <- !is.na(a) & resv[a] > r_low
      a[w] <- parent_row[a[w]]
    }
    
    # keep only those exactly at r_low
    a[is.na(a) | resv[a] != r_low] <- NA_integer_
    a
  }
  
  # Precompute row indices per resolution once
  idx_by_res <- split(seq_len(length(resv)), resv)
  
  out <- list(); tot <- 0L
  for (i_low in seq_len(L - 1L)) {
    r_low <- res_levels[i_low]
    for (i_high in seq.int(i_low + 1L, L)) {
      r_high <- res_levels[i_high]
      idx_high <- idx_by_res[[as.character(r_high)]]
      if (is.null(idx_high) || !length(idx_high)) next
      
      anc <- ancestor_at_vec(idx_high, r_low)
      keep <- which(!is.na(anc))
      co <- sum(sel[idx_high[keep]] == 1L & sel[anc[keep]] == 1L, na.rm = TRUE)
      
      out[[paste0(r_low, "-", r_high)]] <- co
      tot <- tot + co
    }
  }
  c(unlist(out), total = tot)
}

#' Cross-scale coverage summary by resolution
#'
#' @description
#' Summarize how much area is selected at each resolution in a multiscale H3 system.
#' This function adjusts for overlap between coarse and fine planning units so that
#' selected area is not double-counted across resolutions.
#'
#' The output reports, for each resolution, the total available area, the effective
#' selected area (after cross-scale adjustment), and the proportion selected.
#'
#' @param s An `sf` or data frame of planning units. Must contain columns
#'   `res` and `area_km2`, plus the selection column referenced by `flag_col`.
#' @param flag_col Name of the 0/1 selection column to summarise.
#' @param cs_idx A cross-scale index list from [build_crossscale_index()],
#'   containing at least `res_levels`, `rows_by_res`, `pos_in_res`,
#'   `anc_at_res`, and `desc_at_res`.
#'
#' @return A tibble with columns:
#'   * `res` – resolution,
#'   * `area_total_km2` – total area at that resolution,
#'   * `area_selected_km2` – cross-scale selected area, and
#'   * `coverage_prop` – proportion selected.
#'
#' @examples
#' parent <- "872a1072bffffff"
#' kids   <- c("882a1072b1fffff", "882a1072b3fffff")
#'
#' s <- data.frame(
#'   h3_address = c(parent, kids),
#'   res = c(7L, 8L, 8L),
#'   area_km2 = c(14, 2, 2),
#'   selected_by_resolution_final = c(0L, 1L, 0L)
#' )
#' maps   <- build_h3_maps(s, res_levels = c(7L, 8L))
#' cs_idx <- build_crossscale_index(maps)
#'
#' summarize_coverage_by_resolution(s, "selected_by_resolution_final", cs_idx)
#'
#' @export
summarize_coverage_by_resolution <- function(
    s,
    flag_col = "selected_by_resolution_final",
    cs_idx
){
  df <- s |>
    sf::st_drop_geometry() |>
    dplyr::mutate(.flag = dplyr::coalesce(suppressWarnings(as.integer(.data[[flag_col]])), 0L))
  
  res_vec      <- as.integer(df$res)
  res_levels   <- cs_idx$res_levels
  rows_by_res  <- cs_idx$rows_by_res
  pos_in_res   <- cs_idx$pos_in_res
  anc_at_res   <- cs_idx$anc_at_res
  desc_at_res  <- cs_idx$desc_at_res
  
  out <- lapply(res_levels, function(r0){
    idx_r0 <- rows_by_res[[as.character(r0)]]
    if (!length(idx_r0)) {
      return(tibble::tibble(
        res = r0,
        area_total_km2    = 0,
        area_selected_km2 = 0,
        coverage_prop     = NA_real_
      ))
    }
    
    cov <- numeric(length(idx_r0))
    pos_r0 <- pos_in_res[[as.character(r0)]]
    
    # (1) coarser/equal: include EQUAL properly
    coarse_eq <- which(df$.flag == 1L & res_vec <= r0)
    if (length(coarse_eq)) {
      for (i in coarse_eq) {
        if (res_vec[i] == r0) {
          p <- pos_r0[i]
          if (!is.na(p)) cov[p] <- df$area_km2[i]
        } else {
          cells <- desc_at_res[[as.character(r0)]][[i]]
          if (length(cells)) {
            p <- pos_r0[cells]
            cov[p] <- df$area_km2[cells]
          }
        }
      }
    }
    # (2) finer: top-up ancestor, capped by ancestor area
    fine <- .ms_coverage_rows(
      which(df$.flag == 1L & res_vec > r0), res_vec, cs_idx)
    if (length(fine)) {
      for (i in fine) {
        a <- anc_at_res[[as.character(r0)]][i]
        if (!is.na(a)) {
          p <- pos_r0[a]
          if (!is.na(p)) {
            rem <- df$area_km2[a] - cov[p]
            if (rem > 0) cov[p] <- cov[p] + min(df$area_km2[i], rem)
          }
        }
      }
    }
    
    denom   <- sum(df$area_km2[idx_r0], na.rm = TRUE)
    covered <- sum(cov, na.rm = TRUE)
    
    tibble::tibble(
      res               = r0,
      area_total_km2    = denom,
      area_selected_km2 = covered,
      coverage_prop     = ifelse(denom > 0, covered / denom, NA_real_)
    )
  })
  
  dplyr::bind_rows(out) |> dplyr::arrange(res)
}

#' Cross-scale coverage summary by strata
#'
#' @description
#' Summarize how much area is selected within each stratum (e.g. country or region)
#' in a multiscale H3 system. This function reports a cross-scale coverage metric
#' that avoids double-counting when selected planning units overlap across
#' resolutions.
#'
#' For each stratum, coverage is reported at a single "native" resolution chosen
#' from the planning units labeled with that stratum.
#'
#' @param s An `sf` object or data frame of planning units. Must contain columns
#'   `res`, `area_km2`, and the strata column referenced by `strata_col`.
#' @param flag_col Name of the 0/1 selection column to summarize.
#' @param cs_idx A cross-scale index list returned by
#'   [build_crossscale_index()], containing at least `anc_at_res` and `desc_at_res`.
#' @param strata_col Name of the strata column in `s` (default `"strata"`).
#'
#' @return A tibble with one row per stratum and the following columns:
#'   * `strata`
#'   * `native_res`
#'   * `area_total_native_km2`
#'   * `area_selected_km2`
#'   * `coverage_prop`
#'
#' @importFrom stats setNames
#' 
#' @examples
#' parent <- "872a1072bffffff"
#' kids   <- c("882a1072b1fffff", "882a1072b3fffff")
#'
#' s <- data.frame(
#'   h3_address = c(parent, kids),
#'   res = c(7L, 8L, 8L),
#'   strata = c("A", "B", "B"),
#'   area_km2 = c(14, 2, 2),
#'   selected_by_resolution_final = c(0L, 1L, 0L)
#' )
#' maps   <- build_h3_maps(s, res_levels = c(7L, 8L))
#' cs_idx <- build_crossscale_index(maps)
#'
#' summarize_coverage_by_strata(
#'   s,
#'   "selected_by_resolution_final",
#'   cs_idx,
#'   strata_col = "strata"
#' )
#' 
#' @export
summarize_coverage_by_strata <- function(
    s,
    flag_col = "selected_by_resolution_final",
    cs_idx,
    strata_col = "strata"
){
  # Prepare inputs
  df <- s |>
    sf::st_drop_geometry() |>
    dplyr::mutate(.flag = dplyr::coalesce(suppressWarnings(as.integer(.data[[flag_col]])), 0L))
  
  res_vec  <- as.integer(df$res)
  area_km2 <- df$area_km2
  
  # cross-scale index bits
  anc_at_res  <- cs_idx$anc_at_res
  desc_at_res <- cs_idx$desc_at_res
  
  # native resolution per strata = mode of self-labeled rows (ties -> coarser)
  tab <- df |>
    dplyr::filter(!is.na(.data[[strata_col]])) |>
    dplyr::count(!!rlang::sym(strata_col), res, name = "n")
  
  native_res_of <- {
    xs <- split(tab, tab[[strata_col]])
    out <- lapply(xs, function(d){ d <- d[order(-d$n, d$res), ]; d$res[1] })
    stats::setNames(unlist(out), names(xs))
  }
  
  strata_vals <- sort(unique(df[[strata_col]]))
  strata_vals <- strata_vals[!is.na(strata_vals)] # strata to report
  
  # precompute selected indices once, and their coarse/fine splits per r0
  sel_idx_all <- .ms_coverage_rows(which(df$.flag == 1L), res_vec, cs_idx)
  sel_res_all <- res_vec[sel_idx_all]
  
  r0_vals <- unique(as.integer(native_res_of))
  r0_vals <- r0_vals[is.finite(r0_vals)]
  
  sel_coarse_eq_by_r0 <- list()
  sel_fine_by_r0      <- list()
  for (r0 in r0_vals) {
    key <- as.character(r0)
    sel_coarse_eq_by_r0[[key]] <- sel_idx_all[sel_res_all <= r0]
    sel_fine_by_r0[[key]]      <- sel_idx_all[sel_res_all >  r0]
  }
  
  # per-stratum cross-scale coverage (lineage-based)
  out <- lapply(strata_vals, function(A){
    r0 <- native_res_of[[A]]
    if (is.null(r0) || is.na(r0)) {
      return(tibble::tibble(
        strata = A, native_res = NA_integer_,
        area_total_native_km2 = 0,
        area_selected_km2     = 0,
        coverage_prop         = NA_real_
      ))
    }
    
    # denominator: only self-labeled r0 cells for strata A
    idx_native <- which(df[[strata_col]] == A & res_vec == r0)
    if (!length(idx_native)) {
      return(tibble::tibble(
        strata = A, native_res = r0,
        area_total_native_km2 = 0,
        area_selected_km2     = 0,
        coverage_prop         = NA_real_
      ))
    }
    
    denom <- sum(area_km2[idx_native], na.rm = TRUE)
    
    # Numerator: credit any selected PU whose ancestor/descendant at r0 lies in idx_native
    pos_native <- { v <- rep(NA_integer_, nrow(df)); v[idx_native] <- seq_along(idx_native); v }
    cov <- numeric(length(idx_native))
    
    is_native <- rep(FALSE, nrow(df))
    is_native[idx_native] <- TRUE
    
    key <- as.character(r0)
    sel_coarse_eq <- sel_coarse_eq_by_r0[[key]]
    sel_fine      <- sel_fine_by_r0[[key]]
    
    desc_r0 <- desc_at_res[[key]]
    anc_r0  <- anc_at_res[[key]]
    
    # (1) coarser/equal than r0: propagate selection down to native descendants
    if (length(sel_coarse_eq)) {
      for (i in sel_coarse_eq) {
        cells <- if (res_vec[i] == r0) i else (desc_r0[[i]] %||% integer(0))
        if (length(cells)) {
          cells <- cells[is_native[cells]]
          if (length(cells)) cov[pos_native[cells]] <- area_km2[cells]
        }
      }
    }
    
    # (2) finer than r0: top-up the ancestor r0 cell if that ancestor is native to A
    if (length(sel_fine)) {
      for (i in sel_fine) {
        a <- anc_r0[i]
        if (!is.na(a) && is_native[a]) {
          p <- pos_native[a]
          rem <- area_km2[a] - cov[p]
          if (rem > 0) cov[p] <- cov[p] + min(area_km2[i], rem)
        }
      }
    }
    
    covered <- sum(cov, na.rm = TRUE)
    
    tibble::tibble(
      strata                = A,
      native_res            = r0,
      area_total_native_km2 = denom,
      area_selected_km2     = covered,
      coverage_prop         = ifelse(denom > 0, covered / denom, NA_real_)
    )
  })
  
  dplyr::bind_rows(out) |>
    dplyr::arrange(dplyr::desc(coverage_prop))
}

#' Geometry-based feature coverage evaluation from sf column features and selected PUs
#'
#' @description
#' Calculate how well features are represented by a selected set of planning units.
#'
#' This function evaluates feature coverage using feature values stored directly
#' in planning unit attributes. For each feature, it calculates the total amount
#' available, the amount held by the selected planning units, and the proportion
#' of the feature represented by the solution.
#'
#' If targets are provided, the function also reports absolute and relative
#' shortfalls from those targets.
#'
#' Feature coverage is calculated using fractional overlap between planning
#' units and the union of selected planning units.
#'
#' @param s An \code{sf} object of planning units (all resolutions).
#' @param flag_col Name of the 0/1 selection column in \code{s}.
#' @param feature_cols Character vector of feature column names in \code{s}.
#'   Feature names must start with \code{r<res>_} so the native resolution can be inferred.
#' @param res_col Name of the resolution column in \code{s} (default \code{"res"}).
#' @param targets Optional named numeric vector of RELATIVE targets (proportions in \eqn{[0,1]}).
#'   For example, \code{targets = c(r7_featA = 0.30)} means a target of 30\% of the total
#'   available amount of \code{r7_featA} within its native-resolution footprint.
#'   If unnamed, the first value is applied to all features.
#'
#' @return A tibble with one row per feature and columns:
#'   \itemize{
#'     \item \code{feature}, \code{native_res}
#'     \item \code{total} Total feature amount across all PUs at the native resolution.
#'     \item \code{held} Feature amount held by the selected area (via fractional overlap), restricted to the native-resolution footprint.
#'     \item \code{rel_held} Relative held amount (\code{held / total}).
#'     \item \code{target} Relative target proportion (if \code{targets} supplied).
#'     \item \code{abs_shortfall} Absolute shortfall in feature units (\eqn{\max(0, target * total - held)}).
#'     \item \code{rel_shortfall} Relative shortfall in proportions (\eqn{\max(0, target - rel\_held)}).
#'   }
#'   
#' @examples
#' # Minimal polygons
#' p1 <- sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0))))
#' p2 <- sf::st_polygon(list(rbind(c(1, 0), c(2, 0), c(2, 1), c(1, 1), c(1, 0))))
#'
#' s <- sf::st_sf(
#'   res = c(7L, 7L),
#'   area_km2 = c(1, 1),
#'   selected_by_resolution = c(1L, 0L),
#'   r7_featA = c(10, 5),
#'   geometry = sf::st_sfc(p1, p2, crs = 4326)
#' )
#'
#' eval_geom_feature_coverage(
#'   s,
#'   flag_col = "selected_by_resolution",
#'   feature_cols = "r7_featA",
#'   res_col = "res"
#' )
#' 
#' @section Deprecated:
#' This function is deprecated and retained for compatibility. For multiscale 
#' workflows, use [build_multiscale_objective_data()] together with
#' [eval_multiscale_atom_solution()]. The evaluator uses precomputed
#' multiscale data and matching planning-unit row order, and has different
#' arguments and accounting semantics from this geometry-based function.
#'
#' @export
eval_geom_feature_coverage <- function(
    s,
    flag_col = "selected_by_resolution",
    feature_cols,
    res_col = "res",
    targets = NULL
) {
  .Deprecated("eval_multiscale_atom_solution", package = "MultiscaleSCP")

  
  if (!inherits(s, "sf")) stop("`s` must be an sf object.")
  if (!flag_col %in% names(s)) stop("Missing flag_col: ", flag_col)
  if (!res_col %in% names(s)) stop("Missing res_col: ", res_col)
  if (missing(feature_cols) || length(feature_cols) == 0) stop("Provide feature_cols (character vector).")
  
  # geodesic areas in lon/lat
  if ("sf_use_s2" %in% getNamespaceExports("sf")) sf::sf_use_s2(TRUE)
  
  # coerce selection/res
  s[[flag_col]] <- as.integer(s[[flag_col]])
  s[[flag_col]][is.na(s[[flag_col]])] <- 0L
  s[[res_col]] <- as.integer(s[[res_col]])
  
  # Union of selected PUs (any resolution).
  sel_idx <- which(s[[flag_col]] == 1L)
  if (length(sel_idx) == 0) stop("No selected PUs found in ", flag_col, ".")
  
  geom_sel <- sf::st_make_valid(sf::st_union(s[sel_idx, ]))
  
  # helper: target
  get_target <- function(f) {
    if (is.null(targets)) return(NA_real_)
    if (!is.null(names(targets)) && nzchar(names(targets)[1])) return(as.numeric(targets[[f]]))
    as.numeric(targets[1])
  }
  
  # infer native resolution from feature name like "r7_Whatever"
  get_native_res <- function(f) {
    m <- regexpr("^r([0-9]+)_", f, perl = TRUE)
    if (m[1] == -1) NA_integer_
    else as.integer(sub("^r([0-9]+)_.*$", "\\1", f))
  }
  
  #precompute per-resolution fractional coverages
  p_by_res <- list()
  
  for (r in sort(unique(s[[res_col]]))) {
    
    idx <- which(s[[res_col]] == r)
    if (length(idx) == 0) next
    
    sf_r <- sf::st_make_valid(s[idx, ])
    sf_r$.row_id <- idx
    
    area_r <- as.numeric(sf::st_area(sf_r))
    area_r[!is.finite(area_r) | area_r <= 0] <- NA_real_
    names(area_r) <- as.character(idx)
    
    # one intersection for all polygons at this resolution
    inter <- suppressWarnings(sf::st_intersection(sf_r, geom_sel))
    
    if (nrow(inter) == 0) {
      p_by_res[[as.character(r)]] <- stats::setNames(rep(0, length(idx)), idx)
      next
    }
    
    inter_area <- as.numeric(sf::st_area(inter))
    inter_area[!is.finite(inter_area) | inter_area < 0] <- 0
    
    # sum intersection areas by original polygon id
    a_int <- tapply(inter_area, inter$.row_id, sum)
    a_int_full <- rep(0, length(idx))
    names(a_int_full) <- as.character(idx)
    a_int_full[names(a_int)] <- as.numeric(a_int)
    
    p <- a_int_full / area_r
    p[!is.finite(p) | p < 0] <- 0
    p[p > 1] <- 1
    
    p_by_res[[as.character(r)]] <- p
  }
  
  # Evaluate features using precomputed p_j
  out <- lapply(feature_cols, function(f) {
    
    if (!f %in% names(s)) stop("Feature not found in s: ", f)
    
    v <- as.numeric(s[[f]])
    v[!is.finite(v) | is.na(v)] <- 0
    
    r0 <- get_native_res(f)
    if (is.na(r0)) stop("Cannot infer native resolution from feature name: ", f,
                        "\nExpected names like 'r5_*', 'r6_*', 'r7_*'.")
    
    idx_r0 <- which(s[[res_col]] == r0)
    if (length(idx_r0) == 0) {
      return(tibble::tibble(
        feature = f, native_res = r0,
        total = 0, held = 0, rel_held = NA_real_,
        target = get_target(f),
        abs_shortfall = NA_real_, rel_shortfall = NA_real_
      ))
    }
    
    p <- p_by_res[[as.character(r0)]][as.character(idx_r0)]
    p[is.na(p)] <- 0
    
    total <- sum(v[idx_r0])
    held  <- sum(v[idx_r0] * p)
    rel_held <- if (total > 0) held / total else NA_real_
    
    tgt_rel <- get_target(f)       
    tgt_abs <- if (!is.na(tgt_rel)) tgt_rel * total else NA_real_
    
    abs_short <- if (!is.na(tgt_abs)) max(0, tgt_abs - held) else NA_real_
    rel_short <- if (!is.na(tgt_rel) && total > 0) max(0, tgt_rel - rel_held) else NA_real_
    
    tibble::tibble(
      feature = f,
      native_res = r0,
      total = total,
      held = held,
      rel_held = rel_held,
      target = tgt_rel,
      abs_shortfall = abs_short,
      rel_shortfall = rel_short
    )
  })
  
  dplyr::bind_rows(out)
}

#' Exact feature coverage from rasters and selected PUs
#'
#' @description
#' Calculate how well raster-based features are represented by a selected
#' set of planning units.
#'
#' This function evaluates feature coverage from raster layers using
#' exact area-weighted extraction. For each feature, it calculates the
#' total amount available within its native-resolution planning unit
#' footprint, the amount held by the selected planning units, and the
#' proportion of the feature represented by the solution.
#'
#' If targets are provided, the function also reports absolute and
#' relative shortfalls from those targets.
#'
#' @param spp_all A \code{SpatRaster} or \code{RasterStack}/\code{RasterBrick} of feature layers.
#'   Layer names must start with \code{r<res>_} (e.g. \code{r5_}, \code{r6_}, \code{r7_}).
#' @param pu_sf Planning units as an \code{sf} polygon object containing all resolutions.
#'   Must include the resolution column specified by \code{res_col}.
#' @param solution_col Name of the 0/1 (or logical) selection column in \code{pu_sf}
#'   indicating selected planning units.
#' @param targets Optional named numeric vector of RELATIVE targets (proportions in \eqn{[0,1]}).
#'   For example, \code{targets = c(r7_featA = 0.30)} means a target of holding 30\% of the
#'   total available amount of \code{r7_featA} within its native-resolution domain.
#'   If unnamed, the first value is applied to all features.
#' @param res_col Name of the resolution column in \code{pu_sf} (default \code{"res"}).
#'
#' @return A data frame with columns:
#'   \itemize{
#'     \item \code{feature}
#'     \item \code{native_res}
#'     \item \code{total_area} Total feature amount within the native-resolution domain (area-integral units).
#'     \item \code{held_area} Feature amount held by the selected area, restricted to that domain (area-integral units).
#'     \item \code{relative_held} Relative held amount (\code{held_area / total_area}).
#'     \item \code{target} Relative target proportion (if \code{targets} supplied).
#'     \item \code{absolute_shortfall} Absolute shortfall in feature units (\eqn{\max(0, target * total\_area - held\_area)}).
#'     \item \code{relative_shortfall} Relative shortfall in proportions (\eqn{\max(0, target - relative\_held)}).
#'   }
#'
#' @importFrom sf st_union st_intersection st_make_valid
#' @importFrom dplyr bind_rows
#' 
#' @examples
#' \donttest{
#' # Tiny raster with a single layer named with r7_ prefix
#' r <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2)
#' terra::values(r) <- 1
#' names(r) <- "r7_featA"
#'
#' # One square PU covering the raster
#' pu <- sf::st_sf(
#'   res = 7L,
#'   selected = 1L,
#'   geometry = sf::st_sfc(
#'     sf::st_polygon(list(rbind(c(0, 0), c(2, 0), c(2, 2), c(0, 2), c(0, 0)))),
#'     crs = 4326
#'   )
#' )
#'
#' eval_exact_raster_coverage(
#'   spp_all = r,
#'   pu_sf = pu,
#'   solution_col = "selected",
#'   targets = c(r7_featA = 0.3),
#'   res_col = "res"
#' )
#' }
#' 
#' @section Deprecated:
#' This function is deprecated and retained for compatibility. For multiscale
#' workflows, use [build_multiscale_objective_data()] together with
#' [eval_multiscale_atom_solution()]. The atom-based evaluator uses precomputed
#' multiscale data and matching planning-unit row order, and has different
#' arguments and accounting semantics from this raster-based function.
#'
#' @export
eval_exact_raster_coverage <- function(
    spp_all,          # SpatRaster or RasterStack/Brick
    pu_sf,            # Planning units as a single sf
    solution_col,     # name of the column with 0/1 selection
    targets = NULL,   # Optional named target vector (e.g., c("r7_feat" = 0.3))
    res_col = "res"   # resolution column in pu_sf (default "res")
) {
  .Deprecated("eval_multiscale_atom_solution", package = "MultiscaleSCP")

  # convert rasters to a RasterStack for exactextractr
  if (inherits(spp_all, "SpatRaster")) {
    rstack <- raster::stack(spp_all)
  } else if (inherits(spp_all, c("RasterStack","RasterBrick"))) {
    rstack <- spp_all
  } else {
    stop("spp_all must be a SpatRaster or RasterStack/Brick.")
  }
  feat_names <- names(rstack)
  
  # Build selection vector
  if (!solution_col %in% names(pu_sf)) stop("solution_col not found in pu_sf: ", solution_col)
  if (!res_col %in% names(pu_sf)) stop("res_col not found in pu_sf: ", res_col)
  
  sel <- pu_sf[[solution_col]]
  sel[is.na(sel)] <- 0L
  sel <- as.integer(sel)
  if (!any(sel == 1L)) stop("No selected planning units found in ", solution_col, ".")
  
  # union of selected PUs (any resolution)
  geom_sel_all <- sf::st_make_valid(sf::st_union(pu_sf[sel == 1L, ]))
  
  # helper: infer native resolution from feature name like "r7_Arbutus_unedo"
  get_native_res_from_name <- function(nm) {
    m <- regexpr("^r([0-9]+)_", nm, perl = TRUE)
    if (m[1] == -1) return(NA_integer_)
    as.integer(sub("^r([0-9]+)_.*$", "\\1", nm, perl = TRUE))
  }
  
  # helper: get target for a feature
  get_target <- function(f) {
    if (is.null(targets)) return(NA_real_)
    if (!is.null(names(targets)) && nzchar(names(targets)[1])) {
      return(as.numeric(targets[[f]]))
    }
    as.numeric(targets[1])
  }
  
  # Loop over features (domain-aware)
  res_list <- lapply(feat_names, function(f) {
    
    r <- rstack[[f]]
    native_res <- get_native_res_from_name(f)
    if (is.na(native_res)) {
      stop("Could not infer native resolution from layer name: ", f,
           "\nExpected names like 'r5_*', 'r6_*', 'r7_*'.")
    }
    
    pu_native <- pu_sf[as.integer(pu_sf[[res_col]]) == native_res, ]
    if (nrow(pu_native) == 0) {
      return(data.frame(
        feature = f,
        native_res = native_res,
        total_area = 0,
        held_area = 0,
        relative_held = NA_real_,
        target_area = get_target(f),
        absolute_shortfall = NA_real_,
        relative_shortfall = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    
    # Domain where this feature is defined (consistent with the PU feature assignment)
    geom_all_f <- sf::st_make_valid(sf::st_union(pu_native))
    # selected geometry restricted to that domain
    geom_sel_f <- sf::st_make_valid(sf::st_intersection(geom_sel_all, geom_all_f))
    
    # total feature mass in its intended domain
    total <- exactextractr::exact_extract(
      r,
      geom_all_f,
      fun = "sum",
      coverage_area = TRUE,
      progress = FALSE
    )
    
    # held feature mass = selected (any res) but only inside intended domain
    held <- exactextractr::exact_extract(
      r,
      geom_sel_f,
      fun = "sum",
      coverage_area = TRUE,
      progress = FALSE
    )
    
    total <- as.numeric(total)
    held  <- as.numeric(held)
    
    rel_held <- if (is.finite(total) && total > 0) held / total else NA_real_
    
    tgt_rel <- get_target(f)
    tgt_abs <- if (!is.na(tgt_rel)) tgt_rel * total else NA_real_
    
    abs_short <- if (!is.na(tgt_abs)) pmax(0, tgt_abs - held) else NA_real_
    rel_short <- if (!is.na(tgt_rel) && is.finite(total) && total > 0) pmax(0, tgt_rel - rel_held) else NA_real_
    
    data.frame(
      feature            = f,
      native_res         = native_res,
      total_area         = total,
      held_area          = held,
      relative_held      = rel_held,
      target             = tgt_rel,
      absolute_shortfall = abs_short,
      relative_shortfall = rel_short,
      stringsAsFactors   = FALSE
    )
  })
  
  dplyr::bind_rows(res_list)
}
