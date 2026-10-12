test_that("partial hierarchies retain residuals and complete sets do not", {
  f <- atom_hierarchy_fixture()
  a <- build_multiscale_atoms(f$pu, f$maps)
  expect_equal(a$atom_cost, c(5, 1, 1))
  expect_equal(a$atoms_by_pu[[1]], 1:3)
  expect_equal(a$reconstruction_error, 0)
  expect_false(a$complete_partition[1])

  f <- atom_hierarchy_fixture(TRUE)
  # Complete H3 partition: child totals determine the reconciled parent.
  f$pu$area_km2[1] <- 7.5
  a <- build_multiscale_atoms(f$pu, f$maps)
  expect_true(a$complete_partition[1])
  expect_equal(a$atom_cost, rep(1, 7))
  expect_equal(a$reconciled_row_cost[1], 7)
  expect_equal(a$cost_adjustment[1], -0.5)
  expect_false(1L %in% a$atom_rows)
})

test_that("atom evaluation reuses data without double counting", {
  f <- atom_target_fixture()
  f$pu$selected <- c(1L, 1L, 0L)
  a <- eval_multiscale_atom_solution(f$pu, "selected", f$data)
  expect_equal(a$deduplicated_cost, 7)
  expect_equal(a$feature_held, c(A = 3, B = 6))
  expect_equal(a$n_selected_atoms, 3L)
  f$pu$selected <- c(0L, 1L, 1L)
  a <- eval_multiscale_atom_solution(f$pu, "selected", f$data)
  expect_equal(a$deduplicated_cost, 2)
  expect_equal(a$feature_held, c(A = 3, B = 3))
  f$pu$selected <- 0L
  a <- eval_multiscale_atom_solution(f$pu, "selected", f$data)
  expect_equal(a$deduplicated_cost, 0)
  expect_equal(a$feature_held, c(A = 0, B = 0))
})

test_that("atom evaluation validates selections and targets", {
  f <- atom_target_fixture()
  
  f$pu$selected <- c(1L, 2L, 0L)
  expect_error(
    eval_multiscale_atom_solution(f$pu, "selected", f$data),
    "only 0/1"
  )
  
  f$pu$selected <- c(1L, 0L, 0L)
  
  expect_error(
    eval_multiscale_atom_solution(
      f$pu, "selected", f$data,
      targets = c(A = 0.5, C = 0.5)
    ),
    "one target available for every feature"
  )
  
  expect_error(
    eval_multiscale_atom_solution(
      f$pu, "selected", f$data,
      targets = c(A = 0.5, B = 1.2)
    ),
    "between zero and one"
  )
  
  expect_error(
    eval_multiscale_atom_solution(
      f$pu[-1, ], "selected", f$data
    ),
    "same number of planning units"
  )
})

test_that("inconsistent cost totals are reconciled and bad costs rejected", {
  f <- atom_hierarchy_fixture()
  f$pu$area_km2[1] <- 1
  a <- build_multiscale_atoms(f$pu, f$maps)
  expect_equal(a$reconciled_row_cost[1], 2)
  expect_equal(a$reconciliation_action[1], "raised_to_child_total")
  for (bad in c(-1, NA_real_, Inf)) {
    f$pu$area_km2[1] <- bad
    expect_error(build_multiscale_atoms(f$pu, f$maps), "finite and non-negative")
  }
})

test_that("multiscale feature names must be unique", {
  f <- atom_target_fixture()
  expect_error(
    build_multiscale_objective_data(f$pu,f$maps,feature_cols = c(A = "A", A = "B")),
    "must be unique"
  )
})

test_that("feature residuals follow complete H3 partitions", {
  f <- atom_hierarchy_fixture(TRUE)
  f$pu$A <- c(10, rep(1, 7))
  d <- build_multiscale_objective_data(f$pu, f$maps, feature_cols = "A")
  expect_equal(d$feature_atoms$feature_totals, c(A = 7))
  expect_equal(d$feature_atoms$reconciliation_summary$n_reconciled_rows, 1)
})

test_that("free terminal cells retain features and absolute targets", {
  f <- atom_hierarchy_fixture()
  f$pu$area_km2 <- c(1, 0, 1)
  f$pu$A <- c(3, 2, 1)
  d <- build_multiscale_objective_data(f$pu, f$maps, feature_cols = "A")
  expect_true(2L %in% d$atoms$atom_rows)
  expect_equal(d$feature_atoms$feature_totals, c(A = 3))
  expect_equal(d$feature_atoms$dropped_positive, c(A = 0))
  f$pu$selected <- c(0L, 1L, 0L)
  result <- eval_multiscale_atom_solution(f$pu, "selected", d, targets = 2/3)
  expect_equal(result$deduplicated_cost, 0)
  expect_equal(result$feature_held, c(A = 2))
  expect_true(all(result$target_met))
  p <- prioritizr::problem(f$pu, "A", cost_column = "area_km2")
  p <- add_multiscale_min_set(p, d, targets = 2/3)
  p <- prioritizr::add_binary_decisions(p)
  op <- prioritizr::compile(p)
  expect_equal(op$rhs()[op$row_ids() == "ms_feature_target"], 2)
  z <- as.numeric(seq_along(d$atoms$atom_cost) %in% d$atoms$atoms_by_pu[[2]])
  decision <- c(0, 1, 0, z)
  lhs <- drop(as.matrix(op$A()) %*% decision)
  expect_true(all(ifelse(op$sense() == ">=", lhs >= op$rhs() - 1e-9,
                 ifelse(op$sense() == "<=", lhs <= op$rhs() + 1e-9,
                        abs(lhs - op$rhs()) < 1e-9))))
  expect_equal(sum(op$obj() * decision), 0)
})

test_that("a free partial parent remainder retains its feature amount", {
  f <- atom_hierarchy_fixture()
  f$pu$area_km2 <- c(2, 1, 1)
  f$pu$A <- c(5, 1, 1)
  d <- build_multiscale_objective_data(f$pu, f$maps, feature_cols = "A")
  expect_equal(d$atoms$atom_cost, c(0, 1, 1))
  expect_equal(as.numeric(d$feature_atoms$feature_atom_matrix[, 1]), c(3, 1, 1))
  f$pu$selected <- c(1L, 1L, 1L)
  result <- eval_multiscale_atom_solution(f$pu, "selected", d)
  expect_equal(result$deduplicated_cost, 2)
  expect_equal(result$feature_held, c(A = 5))
})

test_that("complete parents have no extra atom even with free children", {
  f <- atom_hierarchy_fixture(TRUE)
  f$pu$area_km2 <- c(6, 0, rep(1, 6))
  f$pu$A <- c(7, rep(1, 7))
  d <- build_multiscale_objective_data(f$pu, f$maps, feature_cols = "A")
  expect_false(1L %in% d$atoms$atom_rows)
  expect_equal(length(d$atoms$atom_cost), 7L)
  expect_equal(d$atoms$atom_cost, c(0, rep(1, 6)))
  f$pu$selected <- c(1L, rep(0L, 7))
  result <- eval_multiscale_atom_solution(f$pu, "selected", d)
  expect_equal(result$deduplicated_cost, 6)
  expect_equal(result$feature_held, c(A = 7))
})

test_that("an entirely free hierarchy still retains all features", {
  f <- atom_hierarchy_fixture(TRUE)
  f$pu$area_km2 <- 0
  f$pu$A <- c(7, rep(1, 7))
  d <- build_multiscale_objective_data(f$pu, f$maps, feature_cols = "A")
  expect_equal(length(d$atoms$atom_cost), 7L)
  expect_equal(d$feature_atoms$feature_totals, c(A = 7))
  f$pu$selected <- c(1L, rep(0L, 7))
  result <- eval_multiscale_atom_solution(f$pu, "selected", d, targets = 1)
  expect_equal(result$deduplicated_cost, 0)
  expect_equal(result$feature_held, c(A = 7))
  expect_true(all(result$target_met))
})
