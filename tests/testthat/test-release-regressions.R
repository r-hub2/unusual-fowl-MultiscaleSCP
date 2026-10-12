# Regression tests for vector feature measurements and nested area accounting.

regression_square <- function(x, y, side = 1000) {
  sf::st_polygon(list(rbind(
    c(x, y), c(x + side, y), c(x + side, y + side),
    c(x, y + side), c(x, y)
  )))
}

regression_vector_fixture <- function() {
  pu <- sf::st_sf(
    h3_address = "872a1072bffffff", res = 7L,
    geometry = sf::st_sfc(regression_square(0, 0, 10000), crs = 3857)
  )
  list(pu = pu, maps = build_h3_maps(pu))
}

regression_chain_fixture <- function(sibling = FALSE) {
  parent <- "872a1072bffffff"
  children <- as.character(unlist(
    h3jsr::get_children(parent, 8L, simple = TRUE), use.names = FALSE
  ))
  grandchildren <- as.character(unlist(
    h3jsr::get_children(children[1], 9L, simple = TRUE), use.names = FALSE
  ))
  pu <- data.frame(
    h3_address = c(parent, children[1], grandchildren[1]),
    res = c(7L, 8L, 9L), area_km2 = c(100, 20, 2),
    feature = c(100, 20, 2), selected = c(0L, 1L, 1L),
    ensemble_score = c(1, 2, 3), strata = c("A", "B", "C")
  )
  if (sibling) {
    pu <- rbind(pu, data.frame(
      h3_address = children[2], res = 8L, area_km2 = 10,
      feature = 10, selected = 1L, ensemble_score = 4, strata = "B"
    ))
  }
  maps <- build_h3_maps(pu)
  list(pu = pu, maps = maps, cs_idx = build_crossscale_index(maps))
}

test_that("single and multipart polygons are measured by area", {
  f <- regression_vector_fixture()
  source <- sf::st_sf(
    value = c(2, 3),
    geometry = sf::st_sfc(
      regression_square(1000, 1000),
      sf::st_multipolygon(list(
        unclass(regression_square(3000, 1000)),
        unclass(regression_square(5000, 1000))
      )), crs = 3857
    )
  )
  x <- load_multiscale_feature(f$pu, source, 7L, f$maps, "mixed")
  expect_equal(x$msopt_mixed, 3, tolerance = 1e-8)
  x <- load_multiscale_feature(
    f$pu, source, 7L, f$maps, "weighted", value_col = "value"
  )
  expect_equal(x$msopt_weighted, 8, tolerance = 1e-8)
})

test_that("single and multipart lines are measured by length", {
  f <- regression_vector_fixture()
  segment <- function(y) rbind(c(1000, y), c(2000, y))
  source <- sf::st_sf(geometry = sf::st_sfc(
    sf::st_linestring(segment(4000)),
    sf::st_multilinestring(list(segment(5000), segment(6000))), crs = 3857
  ))
  x <- load_multiscale_feature(f$pu, source, 7L, f$maps, "lines")
  expect_equal(x$msopt_lines, 3, tolerance = 1e-8)
})

test_that("multipart points are counted individually after clipping", {
  f <- regression_vector_fixture()
  source <- sf::st_sf(value = c(1, 2), geometry = sf::st_sfc(
    sf::st_point(c(1000, 1000)),
    sf::st_multipoint(rbind(c(2000, 2000), c(3000, 3000), c(12000, 12000))),
    crs = 3857
  ))
  x <- load_multiscale_feature(f$pu, source, 7L, f$maps, "points")
  expect_equal(x$msopt_points, 3)
  x <- load_multiscale_feature(
    f$pu, source, 7L, f$maps, "weighted", value_col = "value"
  )
  expect_equal(x$msopt_weighted, 5)
})

test_that("sources with incompatible geometry dimensions fail clearly", {
  f <- regression_vector_fixture()
  source <- sf::st_sf(geometry = sf::st_sfc(
    regression_square(1000, 1000), sf::st_point(c(5000, 5000)), crs = 3857
  ))
  expect_error(
    load_multiscale_feature(f$pu, source, 7L, f$maps, "mixed"),
    "must contain only polygons"
  )
})

test_that("three-level area coverage agrees with the atom union", {
  f <- regression_chain_fixture()
  data <- build_multiscale_objective_data(f$pu, f$maps, feature_cols = "feature")
  selections <- expand.grid(parent = 0:1, child = 0:1, grandchild = 0:1)
  for (i in seq_len(nrow(selections))) {
    f$pu$selected <- as.integer(unlist(selections[i, ], use.names = FALSE))
    summary <- summarize_coverage_by_resolution(f$pu, "selected", f$cs_idx)
    strata <- summarize_coverage_by_strata(f$pu, "selected", f$cs_idx)
    atom <- eval_multiscale_atom_solution(f$pu, "selected", data)
    expected <- if (f$pu$selected[1] == 1L) 100 else
      if (f$pu$selected[2] == 1L) 20 else
        if (f$pu$selected[3] == 1L) 2 else 0
    expect_equal(summary$area_selected_km2[summary$res == 7L], expected)
    expect_equal(strata$area_selected_km2[strata$strata == "A"], expected)
    expect_equal(atom$deduplicated_cost, expected)
  }
})

test_that("distinct branches add area and reordered rows retain coverage", {
  f <- regression_chain_fixture(sibling = TRUE)
  for (order in list(seq_len(nrow(f$pu)), c(4L, 3L, 1L, 2L))) {
    pu <- f$pu[order, ]
    cs_idx <- build_crossscale_index(build_h3_maps(pu))
    summary <- summarize_coverage_by_resolution(pu, "selected", cs_idx)
    strata <- summarize_coverage_by_strata(pu, "selected", cs_idx)
    expect_equal(summary$area_selected_km2[summary$res == 7L], 30)
    expect_equal(strata$area_selected_km2[strata$strata == "A"], 30)
  }
})

test_that("area selectors do not mistake 20 distinct units for 22", {
  f <- regression_chain_fixture()
  by_resolution <- compute_selection_by_resolution(
    f$pu, f$cs_idx,
    target_by_res = c("7" = 0.21, "8" = 0.5, "9" = 1)
  )
  expect_equal(by_resolution$selected_by_resolution, rep(1L, 3))
  masks <- list(A = rep(TRUE, 3), B = rep(TRUE, 3), C = rep(TRUE, 3))
  by_strata <- compute_selection_by_strata(
    f$pu, f$cs_idx, strata_masks = masks, admin_strata = "strata",
    target_by_strata = c(A = 0.21, B = 0.5, C = 1)
  )
  expect_equal(by_strata$selected_by_admin, rep(1L, 3))
})
