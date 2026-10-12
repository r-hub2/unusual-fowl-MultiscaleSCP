test_that("named targets follow feature order and atom totals", {
  f <- atom_target_fixture(c("B", "A"))
  p <- add_multiscale_min_set(f$problem, f$data, c(A = 2/3, B = 1/2))
  p <- prioritizr::add_binary_decisions(p)
  op <- prioritizr::compile(p)
  expect_equal(op$rhs()[op$row_ids() == "ms_feature_target"], c(3, 2))
})

test_that("scalar and positional targets are converted correctly", {
  f <- atom_target_fixture()
  p <- add_multiscale_relative_targets(f$problem, 0.5, f$data)
  expect_equal(as.numeric(p$feature_targets()$value), c(1.5, 3))
  f <- atom_target_fixture()
  p <- add_multiscale_relative_targets(f$problem, c(1, 0.5), f$data)
  expect_equal(as.numeric(p$feature_targets()$value), c(3, 3))
})

test_that("invalid proportions and missing target names are rejected", {
  f <- atom_target_fixture()
  for (bad in list(-0.1, 1.1, NA_real_, Inf, c(A = 0.5, other = 0.5))) {
    expect_error(add_multiscale_relative_targets(f$problem, bad, f$data),
                 "finite proportions")
  }
  expect_error(add_multiscale_relative_targets(f$problem, c(0.1, 0.2, 0.3), f$data),
               "one value")
  expect_error(add_multiscale_relative_targets(f$problem,c(A = 0.5, A = 0.5),f$data),
    "unique feature names")
})
