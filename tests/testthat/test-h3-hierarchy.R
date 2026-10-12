test_that("build_h3_maps validates hierarchy inputs", {
  child <- "8a2a1072b59ffff"
  parent <- h3jsr::get_parent(child, 7L, simple = TRUE)
  h3 <- c(parent, child)
  res <- h3jsr::get_res(h3, simple = TRUE)
  
  expect_silent(build_h3_maps(h3, res))
  
  expect_error(
    build_h3_maps(h3, res[1]),
    "one resolution for every H3 index"
  )
  
  expect_error(
    build_h3_maps(c(parent, parent), c(7L, 7L)),
    "must be unique"
  )
  
  expect_error(
    build_h3_maps(h3, res, res_levels = res[1]),
    "include every"
  )
})