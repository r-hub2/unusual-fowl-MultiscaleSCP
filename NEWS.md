# MultiscaleSCP 0.2.0

* Add a multiscale minimum-set objective for nested H3 planning units, allowing
  feature representation to be shared across resolutions while accounting for
  spatial overlap.
* Add functions for loading conservation features across H3 resolutions and
  preparing multiscale optimization data.
* Add multiscale relative representation targets and solution evaluation for
  feature representation, cost, and target achievement.
* Add and register the Rcpp backend used by the multiscale minimum-set
  formulation.
* Deprecate `eval_geom_feature_coverage()` and
  `eval_exact_raster_coverage()` in favour of
  `eval_multiscale_atom_solution()`.
* Add runnable examples, automated tests, and a reproducible synthetic workflow
  vignette.
* Prevent double counting of nested selections in area coverage summaries
  and area-based selection helpers.
* Require prioritizr (>= 9.0.1).

# MultiscaleSCP 0.1.0

* Initial CRAN submission.
