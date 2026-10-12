#' Compute importance scores for planning units
#'
#' @description
#' This function calculates a continuous importance score
#' (\code{ensemble_score}) for each planning unit by combining two sources
#' of information: (1) a rarity-weighted aggregation of feature values, and
#' (2) how consistently each planning unit is selected across one or more
#' solutions.
#'
#' The relative influence of these components is controlled by
#' \code{alpha_freq}. When multiple solutions are provided, selection
#' consistency reflects how often each planning unit is selected across the
#' portfolio. When a single solution is provided, it reflects whether the
#' planning unit is selected or not. If no solution columns are present, the
#' score is based entirely on the rarity-weighted feature component.
#'
#' Optionally, \code{locked_in} planning units can be forced to rank above
#' all others.
#'
#' @param s An \code{sf} object with one row per planning unit. Must contain
#'   any \code{solution_*} columns used to compute selection consistency, and
#'   optionally a logical \code{locked_in} column.
#' @param feature_mat A numeric matrix or data frame with one row per
#'   planning unit and one column per feature (e.g. proportion of each
#'   feature in each PU).
#' @param feature_weights Optional named numeric vector of user weights
#'   for the features. Names must match \code{colnames(feature_mat)}. If
#'   \code{NULL}, all features receive equal weight.
#' @param alpha_freq Numeric in \eqn{[0,1]} controlling the trade-off between
#'   feature score and selection consistency. Values near 0 emphasize features,
#'   values near 1 emphasize consistency across solutions.
#' @param freq_gamma Numeric exponent applied to selection consistency
#'   (after percentile ranking) to emphasize planning units that are selected
#'   in many solutions (\code{freq_gamma > 1} increases contrast).
#' @param lock_mode Character, either \code{"none"} (locked-in PUs are treated
#'   like any others, aside from their higher selection consistency) or
#'   \code{"top"} (locked-in PUs are always ranked above non-locked ones).
#' @param winsor_p Optional numeric in \eqn{[0,0.5)}. If provided, feature
#'   values are winsorized at the \code{winsor_p} and \code{1 - winsor_p}
#'   quantiles prior to rank/ECDF scaling to \eqn{[0,1]}. If \code{NULL} (default), no
#'   winsorization is performed.
#' @param cost_col Optional name of a numeric column in \code{s} containing
#'   planning-unit costs. If provided, the feature-based score is divided by
#'   \code{cost^cost_beta} prior to final rescaling.
#' @param cost_beta Numeric exponent applied to costs when \code{cost_col} is
#'   provided. Defaults to 1 (benefit per cost).
#'
#' @return The same \code{sf} object \code{s} with two new numeric columns:
#'   \code{selection_consistency} and \code{ensemble_score}.
#'   
#' @examples
#' # Minimal sf with 3 planning units + two solution columns
#' s <- sf::st_as_sf(
#'   data.frame(
#'     x = c(0, 1, 2),
#'     y = c(0, 0, 0),
#'     solution_1 = c(1, 0, 1),
#'     solution_2 = c(1, 1, 0),
#'     locked_in  = c(FALSE, TRUE, FALSE)
#'   ),
#'   coords = c("x", "y"), crs = 4326
#' )
#'
#' # Feature matrix must match nrow(s)
#' feature_mat <- data.frame(
#'   f1 = c(0.2, 0.0, 0.8),
#'   f2 = c(0.0, 0.5, 0.1)
#' )
#'
#' s2 <- compute_pu_scores(
#'   s, feature_mat,
#'   alpha_freq = 0.25,
#'   lock_mode  = "top"
#' )
#'
#' s2[, c("selection_consistency", "ensemble_score", "locked_in")]
#' 
#' @export
compute_pu_scores <- function(
    s, feature_mat, feature_weights = NULL,
    alpha_freq = 0.25,
    freq_gamma = 1.5,
    lock_mode  = c("top","none"), 
    winsor_p   = NULL,
    cost_col   = NULL,
    cost_beta  = 1
){
  lock_mode <- match.arg(lock_mode)
  
  # helpers:
  # ECDF/rank scaling to [0,1]
  ecdf01 <- function(x) {
    x <- as.numeric(x)
    ok <- is.finite(x)
    out <- rep(0, length(x))
    if (!any(ok)) return(out)
    xr <- x[ok]
    # percent-rank with average ranks for ties
    r <- rank(xr, ties.method = "average")
    n <- length(r)
    if (n == 1L) {
      out[ok] <- 1
    } else {
      out[ok] <- (r - 1) / (n - 1)   # in [0,1]
    }
    out
  }
  # winsorize extremes (optional), then ECDF scale to [0,1]
  winsor01 <- function(x, p = NULL) {
    x <- as.numeric(x)
    ok <- is.finite(x)
    if (!any(ok)) return(rep(0, length(x)))
    if (is.null(p)) return(ecdf01(x))
    
    if (!is.numeric(p) || length(p) != 1L || is.na(p) || p < 0 || p >= 0.5) {
      stop("winsor_p must be NULL or a single number in [0, 0.5).")
    }
    
    qs <- stats::quantile(x[ok], probs = c(p, 1 - p), na.rm = TRUE, names = FALSE, type = 7)
    x[x < qs[1]] <- qs[1]
    x[x > qs[2]] <- qs[2]
    ecdf01(x)
  }
  
  df  <- sf::st_drop_geometry(s)
  n   <- nrow(df)
  if (!n) {
    s$selection_consistency <- numeric(0)
    s$ensemble_score <- numeric(0)
    return(s)
  }
  
  # 1) selection consistency (percentile rank of selection counts)
  sol_cols <- grep("^solution_", names(df), value = TRUE)
  k <- if (length(sol_cols)) {
    rowSums(df[, sol_cols, drop = FALSE], na.rm = TRUE)
  } else {
    rep(0, n)
  }
  
  freq <- ecdf01(k)            # percentile rank in [0,1]
  freq <- freq^freq_gamma      # emphasize high-consistency PUs (gamma > 1)
  
  # 2) feature score with rarity-weighted linear combination (WLC)
  if (is.null(dim(feature_mat)) || nrow(feature_mat) != n) {
    stop("feature_mat must be a matrix/data.frame with nrow == nrow(s).")
  }
  feature_mat <- as.matrix(feature_mat)
  
  # user weights (default 1), rarity multiplier R = 1 / total availability
  if (is.null(feature_weights)) {
    w_user <- rep(1, ncol(feature_mat))
    names(w_user) <- colnames(feature_mat)
  } else {
    if (is.null(names(feature_weights))) names(feature_weights) <- colnames(feature_mat)
    w_user <- feature_weights[colnames(feature_mat)]
    w_user[is.na(w_user)] <- 1
  }
  
  total_avail <- colSums(feature_mat, na.rm = TRUE)
  R <- 1 / (total_avail + 1e-9)
  w <- w_user * R
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) {
    w[] <- 1 / length(w)
  } else {
    w <- w / sw
  }
  
  # optional winsorization per feature, then ECDF scale to [0,1]
  Fnorm <- matrix(0, nrow = n, ncol = ncol(feature_mat))
  for (j in seq_len(ncol(feature_mat))) {
    Fnorm[, j] <- winsor01(feature_mat[, j], p = winsor_p)
  }
  
  lin <- as.numeric(Fnorm %*% w)  # WLC score
  
  # optional: cost-adjust the feature score (benefit per cost)
  if (!is.null(cost_col)) {
    if (!cost_col %in% names(df)) stop("cost_col not found in s.")
    cost <- as.numeric(df[[cost_col]])
    cost[!is.finite(cost) | cost <= 0] <- NA_real_
    lin <- lin / (cost^cost_beta)
    lin[!is.finite(lin)] <- 0
  }
  
  lin <- ecdf01(lin)  # normalize to [0,1]
  
  # 3) combine
  score <- (1 - alpha_freq) * lin + alpha_freq * freq
  
  # 4) optional: prioritize locked_in
  locked <- as.logical(df$locked_in)
  locked[is.na(locked)] <- FALSE
  
  if (any(locked) && lock_mode == "top") {
    score[locked] <- 1.01
  }
  
  s$selection_consistency <- freq
  s$ensemble_score <- score
  s
}

#' Compute cross-scale planning-unit scores by resolution
#'
#' @description
#' This function computes a feature-based importance score at each resolution
#' present in a multiscale planning-unit set, and then propagates those
#' resolution-specific scores through the H3 hierarchy so every planning unit
#' receives a score "as seen from" each resolution.
#'
#' Propagation rules:
#' * **Downward (parent \eqn{\rightarrow} children):** undamped broadcast.
#'   All descendants of a given parent at resolution \code{r} inherit the same
#'   \code{score_from_r} value.
#' * **Upward (children \eqn{\rightarrow} parent):** aggregation of descendant
#'   scores at resolution \code{r} using \code{mean} (default) or \code{max}.
#'
#' @param s An \code{sf} object with one row per planning unit. Must contain a
#'   resolution column (default \code{"res"}).
#' @param feature_mat A numeric matrix/data.frame with one row per planning unit
#'   and one column per feature. Feature columns are expected to be prefixed by
#'   resolution (e.g., \code{r5_}, \code{r6_}, ...), unless \code{feature_cols_by_res}
#'   is supplied.
#' @param maps Output of [build_h3_maps()].
#' @param cs_idx Output of [build_crossscale_index()].
#' @param feature_weights Optional named numeric vector of user weights for
#'   features. Names must match \code{colnames(feature_mat)}.
#' @param winsor_p Optional winsorization level passed to the feature-only score
#'   computation. See [compute_pu_scores()].
#' @param cost_col Optional name of a numeric cost column in \code{s}. If
#'   provided, feature-based scores are divided by \code{cost^cost_beta} within
#'   each resolution before propagation.
#' @param cost_beta Numeric exponent applied to costs when \code{cost_col} is
#'   provided.
#' @param agg_mode Character, aggregation used for upward propagation
#'   (children \eqn{\rightarrow} parent). One of \code{"mean"} or \code{"max"}.
#' @param res_col Name of the resolution column in \code{s}. Defaults to \code{"res"}.
#' @param res_levels Optional integer vector of resolutions to compute. Defaults
#'   to \code{maps$res_levels}.
#' @param feature_cols_by_res Optional named list mapping each resolution (as a
#'   character) to a character vector of feature column names. If \code{NULL},
#'   columns are inferred using the prefix pattern \code{^r{res}_}.
#' @param prefix Column name prefix for outputs. Defaults to \code{"score_from_r"}.
#'
#' @return The same \code{sf} object \code{s} with one new numeric column per
#'   resolution, named \code{paste0(prefix, res)} (e.g., \code{score_from_r5}).
#'   
#' @examples
#' # Tiny 2-level setup: 2 parents (r7), and 2 children (r8) under the first parent
#' parent1 <- "872a1072bffffff"
#' parent2 <- "872a10729ffffff"
#' kids    <- c("882a1072b1fffff", "882a1072b3fffff")
#'
#' h3_vec  <- c(parent1, parent2, kids)
#' res_vec <- c(7L, 7L, 8L, 8L)
#'
#' s <- sf::st_as_sf(
#'   data.frame(
#'     h3_address = h3_vec,
#'     res        = res_vec,
#'     cost       = c(1, 1, 1, 1),
#'     x          = c(0, 1, 2, 3),
#'     y          = c(0, 0, 0, 0)
#'   ),
#'   coords = c("x", "y"), crs = 4326
#' )
#'
#' # Feature columns prefixed by resolution (r7_, r8_)
#' feature_mat <- data.frame(
#'   r7_f1 = c(1.0, 0.2, 0.0, 0.0),
#'   r8_f1 = c(0.0, 0.0, 0.6, 0.2)
#' )
#'
#' maps   <- build_h3_maps(s)
#' cs_idx <- build_crossscale_index(maps)
#'
#' s2 <- compute_pu_scores_crossscale(
#'   s, feature_mat,
#'   maps = maps, cs_idx = cs_idx,
#'   agg_mode = "mean",
#'   res_col = "res"
#' )
#'
#' s2[, c("res", "score_from_r7", "score_from_r8")]
#' 
#' @export
compute_pu_scores_crossscale <- function(
    s,
    feature_mat,
    maps,
    cs_idx,
    feature_weights = NULL,
    winsor_p = NULL,
    cost_col = NULL,
    cost_beta = 1,
    agg_mode = c("mean", "max"),
    res_col = "res",
    res_levels = NULL,
    feature_cols_by_res = NULL,
    prefix = "score_from_r"
) {
  agg_mode <- match.arg(agg_mode)
  
  df <- sf::st_drop_geometry(s)
  n <- nrow(df)
  if (!n) return(s)
  
  if (is.null(dim(feature_mat)) || nrow(feature_mat) != n) {
    stop("feature_mat must be a matrix/data.frame with nrow == nrow(s).")
  }
  feature_mat <- as.matrix(feature_mat)
  
  if (!res_col %in% names(df)) stop("res_col not found in s.")
  res_vec <- as.integer(df[[res_col]])
  
  if (is.null(res_levels)) {
    res_levels <- sort(unique(as.integer(maps$res_levels)))
  } else {
    res_levels <- sort(unique(as.integer(res_levels)))
  }
  
  # ECDF/rank scaling to [0,1]:
  ecdf01 <- function(x) {
    x <- as.numeric(x)
    ok <- is.finite(x)
    out <- rep(0, length(x))
    if (!any(ok)) return(out)
    v <- x[ok]
    if (length(unique(v)) <= 1L) {  #flat input: all 0s (stable behaviour)
      out[ok] <- 0
      return(out)
    }
    r <- rank(v, ties.method = "average")
    out[ok] <- (r - 1) / (length(v) - 1)
    out
  }
  
  winsor01 <- function(x, p = NULL) {
    x <- as.numeric(x)
    if (!any(is.finite(x))) return(rep(0, length(x)))
    if (is.null(p)) return(ecdf01(x))
    if (!is.numeric(p) || length(p) != 1L || is.na(p) || p < 0 || p >= 0.5) {
      stop("winsor_p must be NULL or a single number in [0, 0.5).")
    }
    qs <- stats::quantile(x, probs = c(p, 1 - p), na.rm = TRUE, names = FALSE, type = 7)
    x[x < qs[1]] <- qs[1]
    x[x > qs[2]] <- qs[2]
    ecdf01(x)
  }
  
  # costs (optional)
  cost_vec <- NULL
  if (!is.null(cost_col)) {
    if (!cost_col %in% names(df)) stop("cost_col not found in s.")
    cost_vec <- as.numeric(df[[cost_col]])
    cost_vec[!is.finite(cost_vec) | cost_vec <= 0] <- NA_real_
  }
  
  # feature columns per resolution
  if (is.null(feature_cols_by_res)) {
    feature_cols_by_res <- setNames(vector("list", length(res_levels)), as.character(res_levels))
    cn <- colnames(feature_mat)
    for (r in res_levels) {
      r_chr <- as.character(r)
      feature_cols_by_res[[r_chr]] <- grep(paste0("^r", r_chr, "_"), cn, value = TRUE)
    }
  }
  
  # compute base scores on native r-cells (feature-only) for each resolution
  base_by_res <- setNames(vector("list", length(res_levels)), as.character(res_levels))
  
  for (r in res_levels) {
    r_chr <- as.character(r)
    idx_r <- which(res_vec == r)
    cols  <- feature_cols_by_res[[r_chr]]
    
    base <- rep(NA_real_, n)
    
    if (length(idx_r) && length(cols)) {
      X <- feature_mat[idx_r, cols, drop = FALSE]
      
      if (is.null(feature_weights)) {
        w_user <- rep(1, ncol(X))
        names(w_user) <- colnames(X)
      } else {
        if (is.null(names(feature_weights))) names(feature_weights) <- colnames(feature_mat)
        w_user <- feature_weights[colnames(X)]
        w_user[is.na(w_user)] <- 1
      }
      
      total_avail <- colSums(X, na.rm = TRUE)
      R <- 1 / (total_avail + 1e-9)
      w <- w_user * R
      sw <- sum(w)
      if (!is.finite(sw) || sw <= 0) {
        w[] <- 1 / length(w)
      } else {
        w <- w / sw
      }
      
      Fnorm <- matrix(0, nrow = nrow(X), ncol = ncol(X))
      for (j in seq_len(ncol(X))) {
        Fnorm[, j] <- winsor01(X[, j], p = winsor_p)
      }
      
      lin <- as.numeric(Fnorm %*% w)
      
      if (!is.null(cost_vec)) {
        cst <- cost_vec[idx_r]
        lin <- lin / (cst^cost_beta)
        lin[!is.finite(lin)] <- 0
      }
      
      # Apply ECDF scaling
      lin <- ecdf01(lin)
      base[idx_r] <- lin
      
    } else if (length(idx_r)) {
      base[idx_r] <- 0
    }
    
    base_by_res[[r_chr]] <- base
  }
  
  # propagate per-resolution base scores through hierarchy
  for (r in res_levels) {
    r_chr <- as.character(r)
    base <- base_by_res[[r_chr]]
    out  <- rep(NA_real_, n)
    
    anc       <- cs_idx$anc_at_res[[r_chr]]
    desc_list <- cs_idx$desc_at_res[[r_chr]]
    
    for (i in seq_len(n)) {
      if (res_vec[i] == r) {
        out[i] <- base[i]
      } else if (res_vec[i] > r) {
        a <- anc[i]
        if (!is.na(a)) out[i] <- base[a]
      } else {
        d <- desc_list[[i]]
        if (length(d)) {
          v <- base[d]
          v <- v[is.finite(v)]
          if (length(v)) {
            out[i] <- if (agg_mode == "max") max(v) else mean(v)
          }
        }
      }
    }
    
    # rescale propagated scores via ECDF
    out <- ecdf01(out)
    s[[paste0(prefix, r)]] <- out
  }
  
  s
}
