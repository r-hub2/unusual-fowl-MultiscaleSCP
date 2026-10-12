#' Deduplicate multiscale selections across H3 resolutions
#'
#' @description
#' Removes redundant planning-unit selections across nested H3 resolutions.
#' Because finer and coarser H3 cells overlap perfectly (parent–child
#' hierarchy), a selection at one resolution makes selections at other
#' resolutions redundant.
#'
#' This function enforces a consistent hierarchy using one of two strategies:
#'
#' * **`"coarser_first"`** – keep coarser selected cells and drop all selected
#'   descendants (finer cells).  
#'   Useful if coarse-scale representation should dominate.
#'
#' * **`"finer_first"`** – keep finer selected cells and drop selected ancestors
#'   (coarser cells).  
#'   Useful when fine-scale detail should dominate and coarse cells should not
#'   “override” them.
#'
#' @details
#' When \code{mode = "finer_first"}, removing selected coarser cells can reduce
#' the total selected area if only a subset of their descendant cells were
#' selected (i.e., partial coverage of the parent footprint).
#'
#' The function uses the precomputed parent/child
#' lists from [`build_h3_maps()`], and produces a new column
#' `paste0(sel_col, "_deduplicated")`.
#'
#' @param s An `sf` or data frame containing at least the selection column
#'   (`sel_col`).
#' @param sel_col Name of the 0/1 column to deduplicate (e.g. `"solution_1"`).
#' @param h3_vec Character vector of H3 addresses (one per PU; same order as `s`).
#' @param res_vec Integer vector of H3 resolutions (same length/order as `h3_vec`).
#' @param res_levels Vector of resolutions in hierarchical order
#'   (e.g. `c(5, 6, 7)`).
#' @param nearest_parent_row_of Integer vector where each position gives the row
#'   index of the nearest parent cell in `s` (or `NA` if none), as returned by
#'   [`build_h3_maps()`].
#' @param children_by_row List of integer vectors giving, for each PU row,
#'   the row indices of all direct children (finer-resolution descendants),
#'   as returned by [`build_h3_maps()`].
#' @param mode Either `"coarser_first"` or `"finer_first"` (default). Controls
#'   whether coarse or fine PUs are retained when duplicates occur.
#'
#' @return The input `s` with an additional column
#'   `paste0(sel_col, "_deduplicated")` containing a clean 0/1 selection.
#'
#' @examples
#' # One parent (res7) with two children (res8)
#' parent <- "872a1072bffffff"
#' kids   <- c("882a1072b1fffff", "882a1072b3fffff")
#'
#' h3_vec  <- c(parent, kids)
#' res_vec <- c(7L, 8L, 8L)
#'
#' maps <- build_h3_maps(h3_vec, res_vec = res_vec)
#'
#' s <- data.frame(solution_1 = c(1L, 1L, 1L))
#'
#' # Keep the coarser cell (drops children)
#' out_coarse <- deduplicate_h3_selections(
#'   s, sel_col = "solution_1",
#'   h3_vec = maps$h3_vec, res_vec = maps$res_vec, res_levels = maps$res_levels,
#'   nearest_parent_row_of = maps$nearest_parent_row_of,
#'   children_by_row = maps$children_by_row,
#'   mode = "coarser_first"
#' )
#'
#' # Keep the finer cells (drops parent)
#' out_fine <- deduplicate_h3_selections(
#'   s, sel_col = "solution_1",
#'   h3_vec = maps$h3_vec, res_vec = maps$res_vec, res_levels = maps$res_levels,
#'   nearest_parent_row_of = maps$nearest_parent_row_of,
#'   children_by_row = maps$children_by_row,
#'   mode = "finer_first"
#' )
#'
#' out_coarse
#' out_fine
#'
#' @export
deduplicate_h3_selections <- function(
    s,
    sel_col,
    h3_vec, res_vec, res_levels,
    nearest_parent_row_of, children_by_row,
    mode = c("coarser_first","finer_first")
) {
  mode <- match.arg(mode)
  sel <- suppressWarnings(as.integer(s[[sel_col]])); sel[is.na(sel)] <- 0L
  N <- length(sel)
  
  # Collect descendants of row i (BFS over children_by_row)
  drop_all_descendants <- function(i) {
    frontier <- children_by_row[[i]]
    if (!length(frontier)) return(integer(0))
    out <- integer(0)
    seen <- rep(FALSE, N)
    while (length(frontier)) {
      out <- c(out, frontier)
      seen[frontier] <- TRUE
      # Expand frontier to next generation of children
      #if (!length(nxt)) break
      nxt <- integer(0L)
      for (k in seq_along(frontier)) {
        ch <- children_by_row[[frontier[k]]]
        if (length(ch)) nxt <- c(nxt, ch)  # append only
      }
      frontier <- nxt[!seen[nxt]]
    }
    #unique(out)
    out
  }
  
  # helper: walk parents upward to collect ancestors of row i.
  drop_all_ancestors <- function(i) {
    out <- integer(0)
    p <- nearest_parent_row_of[i]
    while (!is.na(p)) { out <- c(out, p); p <- nearest_parent_row_of[p] }
    out
  }
  
  if (mode == "coarser_first") {
    # keep coarser; drop any selected descendants under a kept cell
    for (r in res_levels) {
      idx_r <- which(res_vec == r & sel == 1L)
      if (!length(idx_r)) next
      for (i in idx_r) {
        desc <- drop_all_descendants(i)
        if (length(desc)) sel[desc] <- 0L
      }
    }
  } else { # finer_first
    # keep finer; drop any selected ancestors of kept fine cells
    for (r in rev(res_levels)) {
      idx_r <- which(res_vec == r & sel == 1L)
      if (!length(idx_r)) next
      for (i in idx_r) {
        anc <- drop_all_ancestors(i)
        if (length(anc)) sel[anc] <- 0L
      }
    }
  }
  
  s[[paste0(sel_col, "_deduplicated")]] <- as.integer(sel)
  s
}