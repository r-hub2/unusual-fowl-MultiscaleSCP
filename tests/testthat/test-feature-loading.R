test_that("polygon loading retains partial residuals and native columns", {
  f <- atom_spatial_fixture()
  pu <- load_multiscale_feature(f$pu, f$source, 7L, f$maps, "habitat")
  expect_equal(pu$msopt_habitat, c(7, 1, 1), tolerance = 1e-8)
  expect_equal(pu$r7_habitat, c(7, 0, 0), tolerance = 1e-8)
  expect_equal(multiscale_feature_cols(pu), c(r7_habitat = "msopt_habitat"))
  expect_equal(multiscale_feature_cols(pu, "feature"), c(habitat = "msopt_habitat"))
  expect_error(load_multiscale_feature(pu, f$source, 7L, f$maps, "habitat"),
               "already exists")
  pu2 <- load_multiscale_feature(pu, f$source, 7L, f$maps, "habitat",
                                 overwrite = TRUE)
  expect_equal(pu2$msopt_habitat, pu$msopt_habitat)
})

test_that("fine-native features use the fine domain and aggregate upward", {
  f <- atom_spatial_fixture()
  pu <- load_multiscale_feature(f$pu, f$source, 8L, f$maps, "habitat")
  expect_equal(pu$msopt_habitat, c(2, 1, 1), tolerance = 1e-8)
  expect_equal(pu$r8_habitat, c(0, 1, 1), tolerance = 1e-8)
  d <- build_multiscale_objective_data(pu, f$maps,
      feature_cols = multiscale_feature_cols(pu), native_cols = "r8_habitat")
  expect_equal(unname(d$feature_atoms$feature_totals), 2, tolerance = 1e-8)
  expect_equal(unname(d$checks$native_vs_atom_abs_error), 0, tolerance = 1e-8)
})

test_that("complete child sets avoid a separate parent feature residual", {
  f <- atom_spatial_fixture(TRUE)
  pu <- load_multiscale_feature(f$pu, f$source, 7L, f$maps, "habitat")
  expect_equal(pu$msopt_habitat, c(7, rep(1, 7)), tolerance = 1e-8)
  d <- build_multiscale_objective_data(pu, f$maps,
                                      feature_cols = multiscale_feature_cols(pu))
  expect_equal(length(d$atoms$atom_cost), 7L)
  expect_equal(unname(d$feature_atoms$feature_totals), 7, tolerance = 1e-8)
})

test_that("vector weights, lines, and simple points have known amounts", {
  f <- atom_spatial_fixture()
  f$source$value <- 3
  pu <- load_multiscale_feature(f$pu, f$source, 7L, f$maps, "weighted",
                                value_col = "value")
  expect_equal(pu$msopt_weighted, c(21, 3, 3), tolerance = 1e-8)
  line <- sf::st_sf(geometry = sf::st_sfc(
    sf::st_linestring(rbind(c(0, 500), c(7000, 500))), crs = 3857))
  pu <- load_multiscale_feature(f$pu, line, 7L, f$maps, "line")
  expect_equal(pu$msopt_line, c(7, 1, 1), tolerance = 1e-8)
  points <- sf::st_as_sf(data.frame(x = c(500, 2500), y = c(500, 500)),
                          coords = c("x", "y"), crs = 3857)
  pu <- load_multiscale_feature(f$pu, points, 7L, f$maps, "points")
  expect_equal(pu$msopt_points, c(2, 1, 0))
})

test_that("uniform raster densities integrate over partial footprints", {
  f <- atom_spatial_fixture()
  r <- terra::rast(nrows = 1, ncols = 7, xmin = 0, xmax = 7000,
                    ymin = 0, ymax = 1000, crs = "EPSG:3857")
  terra::values(r) <- 2
  pu <- load_multiscale_feature(f$pu, r, 7L, f$maps, "density")
  expect_equal(pu$msopt_density, c(14, 2, 2), tolerance = 1e-8)
  expect_equal(pu$r7_density, c(14, 0, 0), tolerance = 1e-8)
  expect_error(load_multiscale_feature(f$pu, c(r, r), 7L, f$maps, "bad"),
               "single raster layer")
})
