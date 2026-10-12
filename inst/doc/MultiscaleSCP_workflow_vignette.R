## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE,
  fig.width = 7,
  fig.height = 4
)

library(MultiscaleSCP)
library(sf)
library(prioritizr)

have_highs <- requireNamespace("highs", quietly = TRUE)

## ----planning-units-----------------------------------------------------------
seed <- "87182e586ffffff"

grandparent <- h3jsr::get_parent(seed, 5L, simple = TRUE)

parents <- h3jsr::get_children(grandparent, 6L, simple = TRUE)
parents <- sort(as.character(unlist(parents, use.names = FALSE)))[1:2]

parent_full <- parents[1]
parent_partial <- parents[2]

# get the seven children for the first parent.
children_full <- h3jsr::get_children(parent_full, 7L, simple = TRUE)
children_full <- sort(as.character(unlist(children_full, use.names = FALSE)))

# Only two available children for the second parent.
children_partial_all <- h3jsr::get_children(parent_partial, 7L, simple = TRUE)
children_partial_all <- sort(as.character(unlist(children_partial_all, use.names = FALSE)))
children_partial <- children_partial_all[1:2]

h3_vec <- c(parent_full,children_full,parent_partial,children_partial)

pu <- h3jsr::cell_to_polygon(h3_vec, simple = FALSE)
pu <- pu[match(h3_vec, pu$h3_address), ]
pu <- sf::st_transform(pu, 3857)

pu$res <- c(6L,rep(7L, length(children_full)),6L,rep(7L, length(children_partial)))

# Synthetic additive costs.
pu$cost <- c(14,rep(2, length(children_full)),14,rep(2, length(children_partial)))

maps <- build_h3_maps(pu,res_levels = c(6L, 7L))

knitr::kable(sf::st_drop_geometry(pu)[, c("h3_address", "res", "cost")])

## ----features-----------------------------------------------------------------
full_parent_row <- match(parent_full, pu$h3_address)
partial_parent_row <- match(parent_partial, pu$h3_address)

full_child_rows <- match(children_full, pu$h3_address)
partial_child_rows <- match(children_partial, pu$h3_address)

# The common resolution-6 feature occurs once below each parent.
common_rows <- c(
  full_child_rows[1],
  partial_child_rows[1]
)

feature_common <- sf::st_sf(
  geometry = sf::st_point_on_surface(
    sf::st_geometry(pu)[common_rows]
  )
)

# The unique resolution-6 feature occurs only below the partial parent.
# Put it inside its second available child.
feature_unique <- sf::st_sf(
  geometry = sf::st_point_on_surface(
    sf::st_geometry(pu)[partial_child_rows[2]]
  )
)

# Fine feature: two occurrences below each parent.
fine_rows <- c(full_child_rows[1:2],partial_child_rows)

feature_fine <- sf::st_sf(geometry = sf::st_point_on_surface(
  sf::st_geometry(pu)[fine_rows]))

pu <- load_multiscale_feature(pu = pu,source = feature_common,
  native_res = 6L,maps = maps,feature_name = "FeatureA")

pu <- load_multiscale_feature(pu = pu,source = feature_unique,
  native_res = 6L,maps = maps,feature_name = "FeatureB")

pu <- load_multiscale_feature(pu = pu,source = feature_fine,
  native_res = 7L,maps = maps,feature_name = "FeatureC")

feature_cols <- multiscale_feature_cols(pu)

knitr::kable(sf::st_drop_geometry(pu)[,c("h3_address","res","cost","r6_FeatureA","msopt_FeatureA")])

## ----feature-map, echo=FALSE, fig.width=7, fig.height=4.5---------------------
plot(sf::st_geometry(pu[pu$res == 6L, ]),col = "grey95",border = "grey35",
     lwd = 2,main = "Synthetic multiscale conservation features")

plot(sf::st_geometry(pu[pu$res == 7L, ]),add = TRUE,
  col = "aliceblue",border = "grey60")

# Redraw parent boundaries above the children.
plot(sf::st_geometry(pu[pu$res == 6L, ]),add = TRUE,
  col = NA,border = "grey25",lwd = 2)

plot(sf::st_geometry(feature_common),add = TRUE,pch = 1,
  col = "darkorange",cex = 1.8,lwd = 2)

plot(sf::st_geometry(feature_unique),add = TRUE,
  pch = 4,col = "firebrick",cex = 1.8,lwd = 2)

plot(sf::st_geometry(feature_fine),add = TRUE,
  pch = 16,col = "navy",cex = 1.2)

legend("topleft",
  legend = c("FeatureA (res. 6)","FeatureB (res. 6)","FeatureC (res. 7)"),
  pch = c(1, 4, 16),col = c("darkorange", "firebrick", "navy"),bty = "n")

## ----atoms--------------------------------------------------------------------
ms_data <- build_multiscale_objective_data(pu = pu,maps = maps,feature_cols = feature_cols,
                                           native_cols = names(feature_cols),cost_col = "cost")

## ----problems-----------------------------------------------------------------
targets <- c(r6_FeatureA = 1/2, r6_FeatureB = 1, r7_FeatureC = 1/2)

# Keep the same order used by the problem and atom data.
targets <- targets[names(feature_cols)]

# Conventional minimum-set formulation.
p_classic <- prioritizr::problem(pu,features = names(feature_cols),cost_column = "cost") |>
  prioritizr::add_min_set_objective() |>
  prioritizr::add_relative_targets(targets) |>
  prioritizr::add_binary_decisions()

# Multiscale minimum-set formulation.
p_atom <- prioritizr::problem(pu,features = names(feature_cols),cost_column = "cost") |>
  add_multiscale_min_set(data = ms_data,targets = targets) |>
  prioritizr::add_binary_decisions()

## ----solve, eval=have_highs---------------------------------------------------
s_classic <- solve(prioritizr::add_highs_solver(p_classic,gap = 0,verbose = FALSE))
s_atom <- solve(prioritizr::add_highs_solver(p_atom,gap = 0,verbose = FALSE))

classic_eval <- eval_multiscale_atom_solution(s_classic,"solution_1",ms_data,targets)
atom_eval <- eval_multiscale_atom_solution(s_atom,"solution_1",ms_data,targets)

comparison <- data.frame(
  Approach = c("Conventional minimum-set","Multiscale atom minimum-set"),
  Selected_PUs = c(sum(s_classic$solution_1),sum(s_atom$solution_1)),
  Cost = c(classic_eval$deduplicated_cost,atom_eval$deduplicated_cost),
  FeatureA_held = c(classic_eval$relative_held["r6_FeatureA"],
    atom_eval$relative_held["r6_FeatureA"]),
  FeatureB_held = c(classic_eval$relative_held["r6_FeatureB"],
    atom_eval$relative_held["r6_FeatureB"]),
  FeatureC_held = c(classic_eval$relative_held["r7_FeatureC"],
    atom_eval$relative_held["r7_FeatureC"]),
  Targets_met = c(all(classic_eval$target_met),all(atom_eval$target_met)))

knitr::kable(comparison,digits = 2,
  col.names = c("Approach","Selected PUs","Cost", "Feature A held","Feature B held","Feature C held","All targets met"))

## ----solution-maps, eval=have_highs, fig.width=8, fig.height=4----------------
plot_solution <- function(x, title) {
  parent_rows <- x$res == 6L
  child_rows <- x$res == 7L
  parent_selected <- x$solution_1[parent_rows] == 1
  child_selected <- x$solution_1[child_rows] == 1
  plot(sf::st_geometry(x[parent_rows, ]),
    col = ifelse(parent_selected, "lightblue", "grey95"),border = "grey30",
    lwd = 2,main = title)
  plot(sf::st_geometry(x[child_rows, ]),
    add = TRUE,col = ifelse(child_selected, "steelblue", NA),
    border = "grey60",lwd = 0.8)
  plot(sf::st_geometry(x[parent_rows, ]),add = TRUE,col = NA,border = "grey20",lwd = 2)
}

op <- par(mfrow = c(1, 2),mar = c(0.5, 0.5, 3, 0.5))
plot_solution(s_classic,"Conventional minimum-set")
plot_solution(s_atom,"Multiscale minimum-set")

## ----connectivity-matrix------------------------------------------------------
vertical <- build_multiscale_connectivity_matrix(
  maps,
  symmetric = TRUE
)

vertical

## ----connectivity-example, eval=FALSE-----------------------------------------
# p_connected <- prioritizr::problem(pu,features = names(feature_cols),cost_column = "cost") |>
#   add_multiscale_min_set(data = ms_data,targets = targets) |>
#   add_multiscale_connectivity_penalties(penalty = 0.05, data = vertical) |>
#   prioritizr::add_binary_decisions() |>
#   prioritizr::add_highs_solver(gap = 0,verbose = FALSE)
# 
# # s_connected <- solve(p_connected)

