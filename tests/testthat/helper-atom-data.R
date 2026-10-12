# Genuine H3 relationships; idealized rectangular geometries isolate accounting
# from geodesic effects. Real H3 geometries are exercised in the help examples.
atom_hierarchy_fixture <- function(complete = FALSE) {
  parent <- "872a1072bffffff"
  children <- h3jsr::get_children(parent, res = 8L, simple = TRUE)
  if (is.list(children)) children <- children[[1]]
  children <- as.character(children)
  stopifnot(length(children) == 7L)
  if (!complete) children <- children[1:2]
  n <- length(children)
  pu <- data.frame(id = seq_len(n + 1L),
                   h3_address = c(parent, children),
                   res = c(7L, rep(8L, n)),
                   area_km2 = c(7, rep(1, n)))
  maps <- build_h3_maps(pu)
  list(pu = pu, maps = maps)
}

atom_spatial_fixture <- function(complete = FALSE) {
  f <- atom_hierarchy_fixture(complete)
  rectangle <- function(xmin, xmax) {
    sf::st_polygon(list(rbind(c(xmin, 0), c(xmax, 0),
                             c(xmax, 1000), c(xmin, 1000), c(xmin, 0))))
  }
  n <- nrow(f$pu) - 1L
  polygons <- c(list(rectangle(0, 7000)),
                lapply(seq_len(n), function(i) rectangle((i - 1) * 1000, i * 1000)))
  f$pu <- sf::st_sf(f$pu, geometry = sf::st_sfc(polygons, crs = 3857))
  f$source <- sf::st_sf(value = 1, geometry = sf::st_geometry(f$pu)[1])
  f
}

atom_target_fixture <- function(order = c("A", "B")) {
  f <- atom_hierarchy_fixture()
  f$pu$A <- c(3, 1, 2)
  f$pu$B <- c(6, 2, 1)
  f$data <- build_multiscale_objective_data(
    f$pu, f$maps, feature_cols = c(A = "A", B = "B")
  )
  f$problem <- prioritizr::problem(f$pu, order, cost_column = "area_km2")
  f
}
