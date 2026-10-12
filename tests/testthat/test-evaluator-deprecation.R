test_that("deprecated geometric evaluator retains its result", {
  f <- atom_spatial_fixture()
  f$pu$r7_habitat <- c(7, 0, 0)
  f$pu$selected <- c(1L, 0L, 0L)
  old_s2 <- sf::sf_use_s2()
  on.exit(suppressMessages(sf::sf_use_s2(old_s2)), add = TRUE)
  expect_warning(
    result <- eval_geom_feature_coverage(f$pu, "selected", "r7_habitat"),
    "deprecated"
  )
  expect_equal(result$total, 7)
  expect_equal(result$held, 7)
})

test_that("deprecated raster evaluator retains full-coverage representation", {
  f <- atom_spatial_fixture()
  f$pu$selected <- c(1L, 0L, 0L)
  r <- terra::rast(nrows = 1, ncols = 7, xmin = 0, xmax = 7000,
                    ymin = 0, ymax = 1000, crs = "EPSG:3857")
  terra::values(r) <- 2
  names(r) <- "r7_density"
  expect_warning(
    result <- eval_exact_raster_coverage(r, f$pu, "selected"),
    "deprecated"
  )
  expect_equal(result$relative_held, 1)
  expect_equal(result$held_area, result$total_area)
})
