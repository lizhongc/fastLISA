# fastLISA

Fast, reproducible **Local Indicators of Spatial Association (LISA)** with
arbitrary spatial weights.

`fastLISA` computes seven families of LISA statistics with a plain-C backend,
optional OpenMP multi-threading, and a modern `xoshiro256++` random number
generator for permutation-based inference. It accepts **any** `spdep` `listw`
spatial weights object — including custom and non-contiguity (e.g.
distance-decay) weights — and returns compact, inspectable, `spdep`-compatible
matrices.

## Why fastLISA

The two established R tools force a trade-off:

* **spdep** accepts any `listw` and integrates with the R spatial ecosystem,
  but its conditional-permutation inference runs largely in R and is slow on
  large maps.
* **rgeoda** is fast, but builds its own weights internally — so custom weight
  *values* are ignored (it falls back to binary contiguity) — and returns
  objects that need package-specific accessors.

`fastLISA` closes the gap:

* **Honours custom weight values.** It uses the actual `listw` weights, so
  distance-decay and other non-binary schemes are respected (unlike rgeoda).
* **Fast.** OpenMP-parallelised C kernels are one to two orders of magnitude
  faster than spdep and competitive with rgeoda at equal thread counts.
* **Reproducible at any thread count.** The permutation RNG is re-seeded per
  observation, so for a fixed `iseed` the pseudo-p-values are identical
  regardless of `n.cores` — a guarantee neither spdep nor rgeoda offers.
* **Inspectable output.** Results are plain matrices with `spdep`-style classes
  and cluster/quadrant attributes, not opaque pointers.

## Statistics

| Function            | Statistic                                  |
|---------------------|--------------------------------------------|
| `local_moran()`     | Univariate local Moran's *I*               |
| `local_moran_bv()`  | Bivariate local Moran's *I*                |
| `local_moran_eb()`  | Empirical-Bayes-rate local Moran's *I*     |
| `local_geary()`     | Univariate local Geary's *C*               |
| `local_multigeary()`| Multivariate local Geary's *C*             |
| `local_g()`         | Getis-Ord local *G*                        |
| `local_gstar()`     | Getis-Ord local *G\**                      |

Each returns the observed statistic, a permutation *z*-score, and a pseudo
*p*-value (folded for Moran/G/G\*, tail-adaptive for Geary), with optional
permutation-moment columns. Cluster codes follow `rgeoda` conventions,
including an `Isolated` category for observations with no neighbours.

## Installation

Install the released version from CRAN:

```r
install.packages("fastLISA")
```

Or the development version from source (requires a C99 compiler; OpenMP is used
when available):

```r
# install.packages("remotes")
remotes::install_github("lizhongc/fastLISA")
```

Or from a local clone:

```sh
R CMD INSTALL fastLISA
```

`spdep` is suggested for constructing `listw` weights and for the examples.

## Quick start

```r
library(spdep)
library(fastLISA)

nb <- cell2nb(7, 7)             # 49 cells on a 7 x 7 grid
lw <- nb2listw(nb, style = "W") # row-standardised weights
x  <- as.numeric(seq_len(49))   # a simple gradient

res <- local_moran(x, lw, nsim = 999L, iseed = 1L, n.cores = 1L)
head(res)

attr(res, "cluster")            # High-High / Low-Low / outliers / Isolated ...
```

Custom (e.g. distance-decay) weights are passed through unchanged:

```r
coords  <- as.matrix(expand.grid(x = 1:7, y = 1:7))   # 49 grid points
dnb     <- dnearneigh(coords, 0, 2)                    # neighbours within distance 2
glist   <- lapply(nbdists(dnb, coords), function(d) exp(-d))  # distance decay
lw_exp  <- nb2listw(dnb, glist = glist, style = "B")
res_exp <- local_g(x, lw_exp, nsim = 999L, iseed = 1L)
```

All functions share the same interface: `nsim` permutations, an optional integer
`iseed` for reproducibility, a significance cutoff `p.value`, `n.cores`
(default `1L`; raise it to use multiple OpenMP threads), and `p.method` to choose
the pseudo-p-value method — `"count"` (default) or `spdep`'s ties-averaged
`"rank"`.

## Reproducibility

Because the RNG is re-seeded per observation rather than per thread, the same
`iseed` yields bit-identical pseudo-p-values whether you run on 1 core or many:

```r
a <- local_moran(x, lw, nsim = 999L, iseed = 42L, n.cores = 1L)
b <- local_moran(x, lw, nsim = 999L, iseed = 42L, n.cores = 8L)
identical(c(a), c(b))   # TRUE -- same statistics and pseudo-p-values
```

(`c()` strips attributes; the full objects differ only in the recorded `call`,
which stores the `n.cores` value you passed.)

## Documentation

See the package help (`?local_moran`, `?local_g`, `?local_geary`, ...) and the package vignette:

```r
vignette("fastLISA")
```

## References

* Anselin, L. (1995). Local Indicators of Spatial Association — LISA.
  *Geographical Analysis* 27(2), 93–115.
  <https://doi.org/10.1111/j.1538-4632.1995.tb00338.x>
* Getis, A. & Ord, J. K. (1992). The Analysis of Spatial Association by Use of
  Distance Statistics. *Geographical Analysis* 24(3), 189–206.
  <https://doi.org/10.1111/j.1538-4632.1992.tb00261.x>
* Anselin, L. (2019). A Local Indicator of Multivariate Spatial Association:
  Extending Geary's *c*. *Geographical Analysis* 51(2), 133–150.
  <https://doi.org/10.1111/gean.12164>

## License

GPL-3. 
