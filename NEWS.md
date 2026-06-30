# fastLISA 1.0.0

* First release.
* Local Indicators of Spatial Association computed by a plain-C backend with
  optional OpenMP multi-threading and a `xoshiro256++` permutation RNG:
  * `local_moran()` -- univariate local Moran's I
  * `local_moran_bv()` -- bivariate local Moran's I
  * `local_moran_eb()` -- Empirical-Bayes-rate local Moran's I
  * `local_geary()` -- univariate local Geary's C
  * `local_multigeary()` -- multivariate local Geary's C
  * `local_g()` / `local_gstar()` -- Getis-Ord G and G*
* Accepts any `spdep` `listw` weights object, including custom and
  non-contiguity (distance-decay) weights.
* Returns compact statistic-specific matrices containing the observed
  statistic, permutation Z-score, and pseudo-p-value, with optional
  permutation-moment columns.
* Folded pseudo-p-values are returned for Moran, G, and G* statistics;
  Geary statistics return tail-adaptive pseudo-p-values.
* `p.method` selects the permutation pseudo-p-value method: `"count"` (default)
  counts permutations at least as extreme as the observed value, while `"rank"`
  uses `spdep`'s ties-averaged rank. Both return a folded (smaller-tail) value.
* Cluster codes follow `rgeoda` conventions, including an `Isolated` category
  for observations with no neighbours.
* `n.cores` defaults to `1L`; raise it to use multiple OpenMP threads.
* The permutation RNG is re-seeded per observation, so for a fixed `iseed` the
  pseudo-p-values are identical for any `n.cores` (verified in
  `tests/multicore-reproducible.R`).
