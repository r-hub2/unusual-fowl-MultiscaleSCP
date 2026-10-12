#' Build H3 hierarchy maps for multiscale planning units
#'
#' @description
#' Creates parent--child relationships among H3 planning units represented at
#' multiple resolutions.
#'
#' @param s_or_h3 An `sf` object or data frame containing `h3_address` and
#'   `res` columns, or a character vector of H3 indexes.
#' @param res_vec Optional integer vector of H3 resolutions. Required when
#'   `s_or_h3` is supplied as a character vector.
#' @param res_levels Optional integer vector of resolutions represented in the
#'   hierarchy. Defaults to the unique resolutions in the planning units.
#'
#' @details
#' The supplied `res` or `res_vec` values are used directly and should match
#' the resolutions encoded by the corresponding H3 indexes.
#'
#' @return A list with elements `h3_vec`, `res_vec`, `res_levels`,
#'   `row_idx_by_h3`, `nearest_parent_row_of`, and `children_by_row`.
#'
#' @examples
#' # Minimal example: two resolutions, parent-child relationship
#' h3_child  <- "8a2a1072b59ffff" # example H3 index
#' h3_parent <- "892a1072b5bffff"
#'
#' h3_vec <- c(h3_parent, h3_child)
#' h3_res <- h3jsr::get_res(h3_vec, simple = TRUE)
#'
#' maps <- build_h3_maps(s_or_h3 = h3_vec, res_vec = h3_res)
#'
#' str(maps, max.level = 1)
#' maps$nearest_parent_row_of
#'
#' @importFrom stats setNames
#' @importFrom h3jsr get_parent
#' @export
build_h3_maps <- function(s_or_h3, res_vec = NULL, res_levels = NULL) {
  
  if (is.null(res_vec)) {
    if (
      !is.data.frame(s_or_h3) ||
      !all(c("h3_address", "res") %in% names(s_or_h3))
    ) {
      stop(
        "`s_or_h3` must contain `h3_address` and `res` columns ",
        "when `res_vec` is not supplied."
      )
    }
    
    h3_vec <- as.character(s_or_h3$h3_address)
    res_vec <- suppressWarnings(as.integer(s_or_h3$res))
    
  } else {
    h3_vec <- as.character(s_or_h3)
    res_vec <- suppressWarnings(as.integer(res_vec))
  }
  
  N <- length(h3_vec)
  
  if (N == 0L || length(res_vec) != N) {
    stop("`res_vec` must contain one resolution for every H3 index.")
  }
  
  if (anyNA(h3_vec) || any(!nzchar(h3_vec))) {
    stop("H3 indexes must be non-missing and non-empty.")
  }
  
  if (anyDuplicated(h3_vec)) {
    stop("H3 indexes must be unique.")
  }
  
  if (
    anyNA(res_vec) ||
    any(res_vec < 0L) ||
    any(res_vec > 15L)
  ) {
    stop("H3 resolutions must be integers between 0 and 15.")
  }
  
  res_levels <- if (is.null(res_levels)) {
    sort(unique(res_vec))
  } else {
    sort(unique(as.integer(res_levels)))
  }
  
  if (
    !length(res_levels) ||
    anyNA(res_levels) ||
    any(res_levels < 0L) ||
    any(res_levels > 15L) ||
    !all(res_vec %in% res_levels)
  ) {
    stop(
      "`res_levels` must contain valid H3 resolutions and include every ",
      "resolution present in the planning units."
    )
  }
  
  row_idx_by_h3 <- setNames(seq_len(N), h3_vec)
  
  nearest_parent_row_of <- rep(NA_integer_, N)
  for (i in seq_len(N)) {
    r <- res_vec[i]
    if (r == res_levels[1]) next
    cur_h3 <- h3_vec[i]
    pos <- match(r, res_levels)
    parent_found <- NA_integer_
    for (k in seq.int(pos - 1L, 1L)) {
      r_up <- res_levels[k]
      cur_h3 <- h3jsr::get_parent(cur_h3, r_up, simple = TRUE)
      parent_found <- as.integer(row_idx_by_h3[cur_h3])
      if (!is.na(parent_found)) break
    }
    nearest_parent_row_of[i] <- parent_found
  }
  
  children_by_row <- vector("list", N)
  for (j in seq_len(N)) {
    p <- nearest_parent_row_of[j]
    if (!is.na(p)) children_by_row[[p]] <- c(children_by_row[[p]], j)
  }
  
  list(
    h3_vec = h3_vec,
    res_vec = res_vec,
    res_levels = res_levels,
    row_idx_by_h3 = row_idx_by_h3,
    nearest_parent_row_of = nearest_parent_row_of,
    children_by_row = children_by_row
  )
}

#' Build cross-scale index structures for H3-based SCP workflows
#'
#' @description
#' Extends the basic hierarchy from [build_h3_maps()] into full ancestor,
#' descendant and resolution-index mappings used by all multiscale selection
#' and evaluation functions.
#'
#' @param maps The list returned by [build_h3_maps()].
#'
#' @return A list with elements:
#' * `res_levels`
#' * `rows_by_res`
#' * `pos_in_res`
#' * `anc_at_res`
#' * `desc_at_res`
#' * `finer_rows_by_r0cell`
#'
#' @examples
#' h3_child  <- "8a2a1072b59ffff"
#' h3_parent <- "892a1072b5bffff"
#'
#' h3_vec <- c(h3_parent, h3_child)
#' h3_res <- h3jsr::get_res(h3_vec, simple = TRUE)
#'
#' maps <- build_h3_maps(s_or_h3 = h3_vec, res_vec = h3_res)
#' 
#' cs_idx <- build_crossscale_index(maps)
#'
#' names(cs_idx)
#' cs_idx$res_levels
#' 
#' @export
build_crossscale_index <- function(maps) {
  
  h3_vec  <- maps$h3_vec
  res_vec <- maps$res_vec
  res_levels <- sort(unique(as.integer(maps$res_levels)))
  nearest_parent_row_of <- maps$nearest_parent_row_of
  
  N <- length(h3_vec)
  
  rows_by_res <- setNames(lapply(res_levels, function(r) which(res_vec == r)),
                          as.character(res_levels))
  
  pos_in_res  <- setNames(vector("list", length(res_levels)), as.character(res_levels))
  for (r0 in res_levels) {
    idx <- rows_by_res[[as.character(r0)]]
    pos <- integer(N); pos[] <- NA_integer_
    if (length(idx)) pos[idx] <- seq_along(idx)
    pos_in_res[[as.character(r0)]] <- pos
  }
  
  anc_at_res <- setNames(vector("list", length(res_levels)), as.character(res_levels))
  for (r0 in res_levels) {
    anc <- rep(NA_integer_, N)
    for (i in seq_len(N)) {
      if (res_vec[i] == r0) {
        anc[i] <- i
      } else if (res_vec[i] > r0) {
        p <- nearest_parent_row_of[i]
        while (!is.na(p) && res_vec[p] > r0) p <- nearest_parent_row_of[p]
        if (!is.na(p) && res_vec[p] == r0) anc[i] <- p
      }
    }
    anc_at_res[[as.character(r0)]] <- anc
  }
  
  desc_at_res <- setNames(vector("list", length(res_levels)), as.character(res_levels))
  for (r0 in res_levels) {
    desc_list <- vector("list", N)
    idx_r0 <- rows_by_res[[as.character(r0)]]
    for (j in idx_r0) {
      p <- j
      while (!is.na(p)) {
        desc_list[[p]] <- c(desc_list[[p]], j)
        p <- nearest_parent_row_of[p]
      }
    }
    for (i in seq_len(N)) if (length(desc_list[[i]])) desc_list[[i]] <- sort(unique(desc_list[[i]]))
    desc_at_res[[as.character(r0)]] <- desc_list
  }
  
  finer_rows_by_r0cell <- setNames(vector("list", length(res_levels)), as.character(res_levels))
  for (r0 in res_levels) {
    idx_r0 <- rows_by_res[[as.character(r0)]]
    bucket <- setNames(vector("list", length(idx_r0)), as.character(idx_r0))
    anc <- anc_at_res[[as.character(r0)]]
    finer_idx <- which(res_vec > r0 & !is.na(anc))
    for (i in finer_idx) {
      a <- anc[i]
      bucket[[as.character(a)]] <- c(bucket[[as.character(a)]], i)
    }
    names(bucket) <- as.character(idx_r0)
    finer_rows_by_r0cell[[as.character(r0)]] <- bucket
  }
  
  list(
    res_levels           = res_levels,
    rows_by_res          = rows_by_res,
    pos_in_res           = pos_in_res,
    anc_at_res           = anc_at_res,
    desc_at_res          = desc_at_res,
    finer_rows_by_r0cell = finer_rows_by_r0cell
  )
}
