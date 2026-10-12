#' Stratified multiscale selection by resolution
#'
#' @description
#' Select planning units to meet area-based targets at each H3 resolution.
#' Selection is based on a single score column (e.g., \code{\link{compute_pu_scores}})
#' and a cross-scale index produced by \code{\link{build_crossscale_index}}.
#'
#' Resolutions are processed from finest to coarsest. At each resolution, the
#' function accounts for area already covered by selected finer units to avoid
#' double-counting. To keep selections non-overlapping, it prefers units that do
#' not have finer descendants in the input dataset; coarser units with descendants
#' are only used if needed to reach the target.
#'
#' The result is a resolution-stratified selection that meets per-resolution
#' area targets while minimizing overlap between coarse and fine planning units.
#'
#' @param s An \code{sf} object or data frame of planning units. Must contain at
#'   least the columns referenced by \code{score_col}, \code{area_col}, and
#'   \code{res_col}.
#' @param cs_idx A cross-scale index list as returned by
#'   \code{\link{build_crossscale_index}}.
#' @param target Default area-based target proportion (e.g. \code{0.3}) applied
#'   to all resolutions if \code{target_by_res} is \code{NULL}.
#' @param score_col Name of the column in \code{s} containing the planning-unit
#'   score used to rank candidates (e.g. \code{"ensemble_score"}).
#' @param area_col Name of the column in \code{s} containing planning-unit areas
#'   (e.g. \code{"area_km2"}).
#' @param res_col Name of the column in \code{s} containing integer resolution
#'   codes (e.g. H3 resolution).
#' @param out_col Name of the output column created in \code{s} containing the
#'   final 0/1 selection flag (default \code{"selected_by_resolution"}).
#' @param blocked_col Optional name of a column to store a global "has finer
#'   descendants" flag as 0/1 (1 = has finer descendants). If \code{NULL}
#'   (default), this column is not written.
#' @param target_by_res Optional named numeric vector of per-resolution targets,
#'   e.g. \code{c("5" = 0.30, "6" = 0.20, "7" = 0.15)}. Names should match the
#'   values in \code{res_col}. When supplied, these override \code{target} for
#'   the corresponding resolutions.
#'
#' @return The input \code{s} with an additional column named by \code{out_col}
#'   (0/1), and optionally a column named by \code{blocked_col} if requested.
#'
#' @examples
#' # Build a tiny multiscale set: 2 parents (r7), 2 children (r8) under parent1
#' parent1 <- "872a1072bffffff"
#' parent2 <- "872a10729ffffff"
#' kids1   <- c("882a1072b1fffff", "882a1072b3fffff")
#'
#' h3_vec  <- c(parent1, parent2, kids1)
#' res_vec <- c(7L, 7L, 8L, 8L)
#'
#' s <- data.frame(
#'   h3_address     = h3_vec,
#'   res            = res_vec,
#'   area_km2       = c(14, 14, 2, 2),
#'   ensemble_score = c(0.2, 0.9, 0.8, 0.1)
#' )
#'
#' maps   <- build_h3_maps(s, res_levels = c(7L, 8L))
#' cs_idx <- build_crossscale_index(maps)
#'
#' out <- compute_selection_by_resolution(
#'   s, cs_idx,
#'   target    = 0.5,
#'   score_col = "ensemble_score",
#'   area_col  = "area_km2",
#'   res_col   = "res",
#'   out_col   = "selected_by_resolution"
#' )
#'
#' out[, c("res", "area_km2", "ensemble_score", "selected_by_resolution")]
#'
#' @export
compute_selection_by_resolution <- function(
    s,
    cs_idx,
    target        = 0.30,
    score_col     = "ensemble_score",
    area_col      = "area_km2",
    res_col       = "res",
    out_col       = "selected_by_resolution",
    blocked_col   = NULL,  # If NULL, don't write the blocked column.
    target_by_res = NULL   # Optional per-resolution targets e.g. c("5"=0.3,"6"=0.2,"7"=0.1)
){
  # drop geometry if present
  df <- tryCatch(sf::st_drop_geometry(s), error = function(...) s)
  N  <- nrow(df)
  
  # pull vectors
  res_vec <- as.integer(df[[res_col]])
  area    <- suppressWarnings(as.numeric(df[[area_col]])); area[!is.finite(area)] <- 0
  score   <- suppressWarnings(as.numeric(df[[score_col]])); score[!is.finite(score)] <- -Inf
  
  # index bits from cs_idx
  res_levels    <- cs_idx$res_levels
  rows_by_res   <- cs_idx$rows_by_res
  pos_in_res    <- cs_idx$pos_in_res
  anc_at_res    <- cs_idx$anc_at_res
  finer_buckets <- cs_idx$finer_rows_by_r0cell
  
  # helper to get target for a given resolution r0
  get_target_for_res <- function(r0) {
    if (is.null(target_by_res)) return(target)
    
    # prefer named lookup if available
    if (!is.null(names(target_by_res))) {
      val <- target_by_res[[as.character(r0)]]
    } else {
      # fallback: assume same order as res_levels
      idx <- match(r0, res_levels)
      val <- target_by_res[idx]
    }
    if (is.null(val) || is.na(val)) target else as.numeric(val)
  }
  # precompute "blocked_for_coarser": TRUE if a row has any finer descendants in the dataset
  blocked_global <- rep(FALSE, N)
  for (r0 in res_levels) {
    idx_r0 <- rows_by_res[[as.character(r0)]]
    if (!length(idx_r0)) next
    buckets <- finer_buckets[[as.character(r0)]]
    if (length(buckets)) {
      has_desc <- vapply(buckets, function(x) length(x) > 0L, logical(1))
      if (any(has_desc)) {
        jj <- as.integer(names(buckets)[has_desc])
        blocked_global[jj] <- TRUE
      }
    }
  }
  
  sel <- integer(N)  # output 0/1 selection vector
  
  # Keep a compact list of selected indices to avoid repeated full scans of sel
  sel_idx <- integer(0)
  res_sel <- integer(0)  # res_vec[sel_idx]
  # process from finest -> coarser
  for (r0 in rev(res_levels)) {
    idx_r0 <- rows_by_res[[as.character(r0)]]
    if (!length(idx_r0)) next
    denom <- sum(area[idx_r0], na.rm = TRUE)
    if (!is.finite(denom) || denom <= 0) next
    # resolution-specific target
    t_r0        <- get_target_for_res(r0)
    target_area <- t_r0 * denom
    # coverage on native r0 from already selected finer rows
    cov_native <- numeric(length(idx_r0))
    pos_r0     <- pos_in_res[[as.character(r0)]]
    # only look at selected indices
    fine_idx <- .ms_coverage_rows(sel_idx[res_sel > r0], res_vec, cs_idx)
    if (length(fine_idx)) {
      a <- anc_at_res[[as.character(r0)]][fine_idx]
      keep <- which(!is.na(a))
      if (length(keep)) {
        a <- a[keep]
        gain_by_a <- tapply(area[fine_idx[keep]], a, sum, na.rm = TRUE)
        a_idx <- as.integer(names(gain_by_a))
        # cap by ancestor area
        capped <- pmin(as.numeric(gain_by_a), area[a_idx])
        cov_native[pos_r0[a_idx]] <- capped
      }
    }
    
    covered_now <- sum(cov_native, na.rm = TRUE)
    need <- max(0, target_area - covered_now)
    
    if (need > 0) {
      # Tier 1: native cells without finer descendants
      cand1 <- idx_r0[(sel[idx_r0] == 0L) & (!blocked_global[idx_r0])]
      if (length(cand1)) {
        g <- pmax(0, area[cand1] - cov_native[pos_r0[cand1]])
        keep <- which(g > 0)
        if (length(keep)) {
          cand1 <- cand1[keep]
          ord <- order(score[cand1], decreasing = TRUE)
          for (j in cand1[ord]) {
            if (need <= 0) break
            mg <- max(0, area[j] - cov_native[pos_r0[j]])
            if (mg <= 0) next
            sel[j] <- 1L
            sel_idx <- c(sel_idx, j)
            res_sel <- c(res_sel, res_vec[j])
            
            cov_native[pos_r0[j]] <- area[j]
            need <- need - mg
          }
        }
      }
    }
    if (need > 0) {
      # 2) fallback: allow blocked parents (those with finer descendants)
      cand2 <- idx_r0[(sel[idx_r0] == 0L) & (blocked_global[idx_r0])]
      if (length(cand2)) {
        g <- pmax(0, area[cand2] - cov_native[pos_r0[cand2]])
        keep <- which(g > 0)
        if (length(keep)) {
          cand2 <- cand2[keep]
          ord <- order(score[cand2], decreasing = TRUE)
          for (j in cand2[ord]) {
            if (need <= 0) break
            mg <- max(0, area[j] - cov_native[pos_r0[j]])
            if (mg <= 0) next
            sel[j] <- 1L
            sel_idx <- c(sel_idx, j)
            res_sel <- c(res_sel, res_vec[j])
            
            cov_native[pos_r0[j]] <- area[j]
            need <- need - mg
          }
        }
      }
    }
    # If still need > 0: infeasible at this res; move on.
  }
  # write outputs
  s[[out_col]] <- as.integer(sel)
  # Write blocked_col only when requested.
  if (!is.null(blocked_col)) {
    s[[blocked_col]] <- as.integer(blocked_global) #1=has finer descendants,0=no descendants
  }
  
  s
}

#' Stratified multiscale selection by strata and resolution
#'
#' @description
#' Select planning units to meet area-based targets for each strata (e.g. country)
#' at each H3 resolution. Selection is based on a single score column (typically
#' \code{ensemble_score}) and a cross-scale index produced by
#' \code{\link{build_crossscale_index}}.
#'
#' Resolutions are processed from finest to coarsest. At each strata and resolution,
#' the function accounts for area already covered by selected planning units at
#' other resolutions to avoid double-counting. Within each strata, it prefers
#' native-resolution units that do not have finer descendants in the input dataset,
#' using coarser units with descendants only when needed to reach the target.
#'
#' The result is a single 0/1 selection that meets per-strata area targets while
#' minimizing cross-scale overlap.
#'
#' @param s An `sf` object or data frame of planning units. Must contain at
#'   least the columns referenced by `admin_strata`, `score_col`, `area_col`,
#'   and `res_col`.
#' @param cs_idx A cross-scale index list as returned by
#'   [build_crossscale_index()], containing at least `res_levels`,
#'   `rows_by_res`, `anc_at_res`, `desc_at_res`, and `finer_rows_by_r0cell`.
#' @param strata_masks Optional named list of logical vectors (length
#'   `nrow(s)`) defining the "in-scope" units for each strata. Names should
#'   match the values in `admin_strata`. If a given strata is missing from
#'   `strata_masks`, a fallback mask `admin_strata == A` is used.
#' @param target Default area-based target proportion (e.g. `0.3`) applied to
#'   all strata if `target_by_strata` is `NULL` or does not contain an entry
#'   for that strata.
#' @param admin_strata Name of the column in `s` containing strata labels
#'   (e.g. `"admin"`).
#' @param score_col Name of the column in `s` containing the PU score to rank
#'   candidates.
#' @param area_col Name of the column in `s` with PU areas.
#' @param res_col Name of the column in `s` with integer resolution codes.
#' @param out_col Name of the output column that will be created in `s`
#'   containing the final 0/1 selection flag (default `"selected_by_admin"`).
#' @param blocked_col Optional name of a column to store the global
#'   "has finer descendants" flag as 0/1. If `NULL` (default), this column
#'   is not written.
#' @param target_by_strata Optional named numeric vector of per-strata targets,
#'   e.g. `c("Ireland" = 0.30, "Norway" = 0.25)`. Names should match the
#'   values in `admin_strata`. When present, these override `target` for the
#'   corresponding strata.
#'
#' @return The input `s` with an additional column `out_col` (0/1), and
#'   optionally `blocked_col` if requested.
#'
#' @examples
#' # Tiny multiscale example with two strata (A and B)
#' parent1 <- "872a1072bffffff"
#' kids1   <- c("882a1072b1fffff", "882a1072b3fffff")
#' parent2 <- "872a10729ffffff"
#' parent3 <- "872a10774ffffff"
#'
#' h3_vec  <- c(parent1, parent2, parent3, kids1)
#' res_vec <- c(7L, 7L, 7L, 8L, 8L)
#'
#' s <- data.frame(
#'   h3_address     = h3_vec,
#'   res            = res_vec,
#'   admin          = c("A", "B", "B", "A", "A"),
#'   area_km2       = c(14, 14, 14, 2, 2),
#'   ensemble_score = c(0.2, 0.9, 0.2, 0.8, 0.1)
#' )
#'
#' maps   <- build_h3_maps(s, res_levels = c(7L, 8L))
#' cs_idx <- build_crossscale_index(maps)
#'
#' strata_masks <- list(
#'   A = (s$admin == "A"),
#'   B = (s$admin == "B")
#' )
#'
#' out <- compute_selection_by_strata(
#'   s, cs_idx,
#'   strata_masks     = strata_masks,
#'   target           = 0.3,
#'   admin_strata     = "admin",
#'   score_col        = "ensemble_score",
#'   area_col         = "area_km2",
#'   res_col          = "res",
#'   out_col          = "selected_by_admin",
#'   target_by_strata = c(A = 0.12, B = 0.5)
#' )
#'
#' out[, c("admin", "res", "area_km2", "ensemble_score", "selected_by_admin")]
#'
#' @export
compute_selection_by_strata <- function(
    s,
    cs_idx,
    strata_masks,
    target           = 0.30,
    admin_strata     = "admin",
    score_col        = "ensemble_score",
    area_col         = "area_km2",
    res_col          = "res",
    out_col          = "selected_by_admin",
    blocked_col      = NULL,    #optional output column; if NULL, nothing is written
    target_by_strata = NULL
) {
  df <- tryCatch(sf::st_drop_geometry(s), error = function(...) s)
  N  <- nrow(df)
  
  res_vec <- as.integer(df[[res_col]])
  area    <- suppressWarnings(as.numeric(df[[area_col]]));  area[!is.finite(area)] <- 0
  score   <- suppressWarnings(as.numeric(df[[score_col]])); score[!is.finite(score)] <- -Inf
  admin_v <- df[[admin_strata]]
  
  res_levels    <- cs_idx$res_levels
  rows_by_res   <- cs_idx$rows_by_res
  anc_at_res    <- cs_idx$anc_at_res
  desc_at_res   <- cs_idx$desc_at_res
  finer_buckets <- cs_idx$finer_rows_by_r0cell
  
  # global "has finer descendants" flag (used internally only)
  blocked_global <- rep(FALSE, N)
  for (r0 in res_levels) {
    idx_r0 <- rows_by_res[[as.character(r0)]]
    if (!length(idx_r0)) next
    buckets <- finer_buckets[[as.character(r0)]]
    if (!length(buckets)) next
    has_desc <- vapply(buckets, function(x) length(x) > 0L, logical(1))
    if (any(has_desc)) {
      jj <- as.integer(names(buckets)[has_desc])
      blocked_global[jj] <- TRUE
    }
  }
  
  # target resolver for strata A
  get_target_for <- function(A) {
    if (!is.null(target_by_strata) && length(target_by_strata) > 0L) {
      nms <- names(target_by_strata)
      if (!is.null(nms) && length(nms) > 0L) {
        keyA <- as.character(A)
        if (!is.na(keyA) && keyA %in% nms) {
          val <- target_by_strata[[keyA]]
          if (!is.null(val) && is.finite(val)) return(as.numeric(val))
        }
      }
    }
    target
  }
  
  sel <- integer(N)
  
  # coverage on native r0 cells
  cov_for_admin_res <- function(r0, idx_native, pos_native) {
    cov <- numeric(length(idx_native))
    
    sel1 <- .ms_coverage_rows(which(sel == 1L), res_vec, cs_idx)
    if (!length(sel1)) return(cov)
    
    # coarser/equal selections (any strata)
    coarse_eq <- sel1[res_vec[sel1] <= r0]
    if (length(coarse_eq)) {
      desc_r0 <- desc_at_res[[as.character(r0)]]
      idx_native_area <- area[idx_native]
      for (i in coarse_eq) {
        cells <- desc_r0[[i]]
        if (!length(cells)) next
        pp <- pos_native[cells]
        pp <- pp[!is.na(pp)]
        if (!length(pp)) next
        cov[pp] <- idx_native_area[pp]
      }
    }
    
    # finer selections (any strata)
    fine <- sel1[res_vec[sel1] > r0]
    if (length(fine)) {
      anc_r0 <- anc_at_res[[as.character(r0)]]
      for (i in fine) {
        a <- anc_r0[i]
        if (is.na(a)) next
        p <- pos_native[a]
        if (is.na(p)) next
        rem <- area[a] - cov[p]
        if (rem > 0) cov[p] <- cov[p] + min(area[i], rem)
      }
    }
    cov
  }
  
  # marginal gain for candidate i (uses in_scope), but avoids intersect/%in%
  mgain_admin_res <- function(i, r0, covA, pos_native, idx_native_area, in_scope) {
    if (!in_scope[i]) return(0)
    
    if (res_vec[i] <= r0) {
      cells <- desc_at_res[[as.character(r0)]][[i]]
      if (!length(cells)) return(0)
      pp <- pos_native[cells]
      pp <- pp[!is.na(pp)]
      if (!length(pp)) return(0)
      sum(pmax(0, idx_native_area[pp] - covA[pp]), na.rm = TRUE)
    } else {
      a <- anc_at_res[[as.character(r0)]][i]
      if (is.na(a)) return(0)
      p <- pos_native[a]
      if (is.na(p)) return(0)
      rem <- idx_native_area[p] - covA[p]
      if (rem <= 0) return(0)
      selected <- .ms_coverage_rows(which(sel == 1L), res_vec, cs_idx)
      already_covered <- any(vapply(anc_at_res, function(a) {
        a[i] %in% selected
      }, logical(1)))
      if (already_covered) return(0)
      ancestors <- anc_at_res[[as.character(res_vec[i])]][selected]
      descendants <- selected[which(
        res_vec[selected] > res_vec[i] & ancestors == i
      )]
      new_area <- max(0, area[i] - sum(area[descendants], na.rm = TRUE))
      min(new_area, rem)
    }
  }
  
  # resolutions: finest -> coarsest
  for (r0 in rev(res_levels)) {
    idx_r0 <- rows_by_res[[as.character(r0)]]
    if (!length(idx_r0)) next
    
    admins_r0 <- unique(admin_v[idx_r0])
    admins_r0 <- admins_r0[!is.na(admins_r0)]
    if (!length(admins_r0)) next
    
    for (A in admins_r0) {
      in_scope <- strata_masks[[A]]
      if (is.null(in_scope)) {
        in_scope <- (admin_v == A)
      }
      if (!is.logical(in_scope) || length(in_scope) != N) next
      
      idx_native <- which((admin_v == A) & (res_vec == r0) & in_scope)
      if (!length(idx_native)) next
      
      pos_native <- rep(NA_integer_, N)
      pos_native[idx_native] <- seq_along(idx_native)
      
      idx_native_area <- area[idx_native]
      denom <- sum(idx_native_area, na.rm = TRUE)
      if (!is.finite(denom) || denom <= 0) next
      
      t_A         <- get_target_for(A)
      target_area <- t_A * denom
      
      covA        <- cov_for_admin_res(r0, idx_native, pos_native)
      covered_now <- sum(covA, na.rm = TRUE)
      need        <- max(0, target_area - covered_now)
      if (need <= 0) next
      
      ## Tier 1: native r0, not blocked
      cand1 <- idx_native[(sel[idx_native] == 0L) & (!blocked_global[idx_native])]
      if (length(cand1)) {
        g <- pmax(0, area[cand1] - covA[pos_native[cand1]])
        keep <- which(g > 0)
        if (length(keep)) {
          cand1 <- cand1[keep]
          ord   <- order(score[cand1], decreasing = TRUE)
          for (j in cand1[ord]) {
            if (need <= 0) break
            mg <- max(0, area[j] - covA[pos_native[j]])
            if (mg <= 0) next
            sel[j] <- 1L
            covA[pos_native[j]] <- area[j]
            need <- need - mg
          }
        }
      }
      if (need <= 0) next
      
      ## Tier 2: native r0, blocked
      cand2 <- idx_native[(sel[idx_native] == 0L) & (blocked_global[idx_native])]
      if (length(cand2)) {
        g <- pmax(0, area[cand2] - covA[pos_native[cand2]])
        keep <- which(g > 0)
        if (length(keep)) {
          cand2 <- cand2[keep]
          ord   <- order(score[cand2], decreasing = TRUE)
          for (j in cand2[ord]) {
            if (need <= 0) break
            mg <- max(0, area[j] - covA[pos_native[j]])
            if (mg <= 0) next
            sel[j] <- 1L
            covA[pos_native[j]] <- area[j]
            need <- need - mg
          }
        }
      }
      if (need <= 0) next
      
      ## Tier 3: any in-scope row by marginal gain
      pool <- which((sel == 0L) & in_scope)
      if (length(pool)) {
        mg <- vapply(
          pool,
          mgain_admin_res,
          numeric(1),
          r0              = r0,
          covA            = covA,
          pos_native      = pos_native,
          idx_native_area = idx_native_area,
          in_scope        = in_scope
        )
        keep <- which(mg > 0)
        if (length(keep)) {
          pool <- pool[keep]; mg <- mg[keep]
          ord <- order(score[pool], mg, decreasing = TRUE)
          for (i in pool[ord]) {
            if (need <= 0) break
            g <- mgain_admin_res(i, r0, covA, pos_native, idx_native_area, in_scope)
            if (g <= 0) next
            sel[i] <- 1L
            
            if (res_vec[i] <= r0) {
              cells <- desc_at_res[[as.character(r0)]][[i]]
              if (length(cells)) {
                pp <- pos_native[cells]
                pp <- pp[!is.na(pp)]
                if (length(pp)) covA[pp] <- idx_native_area[pp]
              }
            } else {
              a <- anc_at_res[[as.character(r0)]][i]
              if (!is.na(a)) {
                p <- pos_native[a]
                if (!is.na(p)) {
                  covA[p] <- covA[p] + g
                }
              }
            }
            need <- need - g
          }
        }
      }
    }
  }
  
  s[[out_col]] <- as.integer(sel)
  
  # only write blocked_col if asked to
  if (!is.null(blocked_col)) {
    s[[blocked_col]] <- as.integer(blocked_global)
  }
  
  s
}
