# No licensed solver or external spatial data is needed: enumerate this small
# model exactly and verify compiled constraints against the defining equations.
fixture <- function(order = c("A", "B")) {
  q <- cbind(A = c(1, 1, 1, 0, 0, 0, 0), B = c(1, 1, 1, 1, 1, 1, 0))
  pu <- data.frame(id = seq_len(8), cost = c(7, rep(1, 7)),
                   A = c(3, q[, "A"]), B = c(6, q[, "B"]))
  d <- structure(list(
    atoms = list(atom_cost = rep(1, 7),
                 cover_by_atom = lapply(2:8, function(i) c(1L, i))),
    feature_atoms = list(feature_atom_matrix = methods::as(q, "dgCMatrix"),
                         feature_names = c("A", "B"),
                         feature_totals = c(A = 3, B = 6))
  ), class = "multiscale_objective_data")
  problem <- prioritizr::problem(pu, order, cost_column = "cost")
  problem <- MultiscaleSCP::add_multiscale_min_set(problem, d, c(A = 2/3, B = 1/2))
  problem <- prioritizr::add_binary_decisions(problem)
  list(problem = problem, data = d, q = q)
}

binary_assignments <- function(n) {
  as.matrix(expand.grid(rep(list(c(0, 1)), n)))
}

feasible <- function(op, v) {
  lhs <- drop(as.matrix(op$A()) %*% v)
  ok <- switch_sense(lhs, op$rhs(), op$sense())
  all(ok) && all(v >= op$lb() - 1e-9 & v <= op$ub() + 1e-9)
}
switch_sense <- function(lhs, rhs, sense) {
  ifelse(sense == ">=", lhs >= rhs - 1e-9,
         ifelse(sense == "<=", lhs <= rhs + 1e-9, abs(lhs - rhs) < 1e-9))
}

test_that("all 256 PU selections agree with atom accounting", {
  f <- fixture(); op <- prioritizr::compile(f$problem)
  
  expect_identical(op$obj_id(), "min_set")
  expect_true(prioritizr::presolve_check(op, warn = FALSE))
  
  expect_equal(op$ncol(), 15)
  expect_equal(op$obj(), c(rep(0, 8), rep(1, 7)))
  expect_equal(op$vtype(), c(rep("B", 8), rep("C", 7)))
  X <- binary_assignments(8)
  costs <- rep(Inf, nrow(X))
  for (k in seq_len(nrow(X))) {
    x <- X[k, ]; z <- as.numeric(x[1] == 1 | x[2:8] == 1)
    held <- drop(crossprod(f$q, z))
    expected <- all(held >= c(2, 3))
    expect_identical(feasible(op, c(x, z)), expected)
    if (expected) costs[k] <- sum(z)
  }
  expect_equal(min(costs), 3)
  expect_true(feasible(op, c(1, rep(0, 7), rep(1, 7))))
  expect_true(feasible(op, c(rep(1, 8), rep(1, 7))))
  expect_equal(sum(op$obj() * c(rep(1, 8), rep(1, 7))), 7)
})

test_that("links uniquely force the correct atoms, including zero-cost atoms", {
  f <- fixture(); f$data$atoms$atom_cost[7] <- 0
  expect_warning(
    p <- MultiscaleSCP::add_multiscale_min_set_objective(f$problem, f$data),
    "Overwriting previously defined objective"
  )
  op <- prioritizr::compile(p)
  link <- op$row_ids() %in% c("ms_atom_upper", "ms_atom_lower")
  A <- as.matrix(op$A())[link, , drop = FALSE]
  X <- binary_assignments(8); Z <- binary_assignments(7)
  for (k in seq_len(nrow(X))) {
    x <- X[k, ]; expected <- as.numeric(x[1] == 1 | x[2:8] == 1)
    values <- sweep(A[, 9:15, drop = FALSE] %*% t(Z),
                    1, drop(A[, 1:8, drop = FALSE] %*% x), "+")
    valid <- colSums(values > 1e-9) == 0
    expect_equal(sum(valid), 1)
    expect_equal(as.numeric(Z[valid, ]), expected)
  }
})

test_that("feature order and locked decisions survive reconstruction", {
  op <- prioritizr::compile(fixture(c("B", "A"))$problem)
  target <- op$row_ids() == "ms_feature_target"
  expect_equal(op$rhs()[target], c(3, 2))
  expect_equal(unname(as.matrix(op$A())[target, 9:15]),
               unname(t(fixture()$q[, c("B", "A")])))
  p <- prioritizr::add_locked_out_constraints(fixture()$problem, 1L)
  op <- prioritizr::compile(p)
  expect_equal(op$ub()[1], 0)
  expect_false(feasible(op, c(1, rep(0, 7), rep(1, 7))))
})

test_that("unsupported decisions fail explicitly", {
  expect_warning(
    p <- prioritizr::add_proportion_decisions(fixture()$problem),
    "Overwriting previously defined decision"
  )
  expect_error(prioritizr::compile(p), "binary")
  expect_error(prioritizr::compile(fixture()$problem, compressed_formulation = FALSE),
               "compressed")
})

test_that("invalid atom mapping fails before mutation", {
  f <- fixture(); f$data$atoms$cover_by_atom[[1]] <- c(1L, 99L)
  expect_warning(
    p <- MultiscaleSCP::add_multiscale_min_set_objective(f$problem, f$data),
    "Overwriting previously defined objective"
  )
  expect_error(prioritizr::compile(p), "mapping")
})

test_that("multiscale minimum-set solves with HiGHS", {
  testthat::skip_if_not_installed("highs")
  
  f <- fixture()
  
  p <- prioritizr::add_highs_solver(
    f$problem,
    gap = 0,
    verbose = FALSE
  )
  
  s <- solve(p)
  
  expect_identical(
    unname(attr(s, "status")),
    "Optimal"
  )
  
  expect_equal(
    unname(attr(s, "objective")),
    3
  )
  
  expect_equal(
    unname(attr(s, "gap")),
    0
  )
})