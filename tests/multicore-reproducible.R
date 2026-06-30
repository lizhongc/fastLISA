## Per-observation RNG seeding must make results thread-independent: for a fixed
## iseed the output is bit-identical for any n.cores.  Two cores keeps the test
## within CRAN's example/test thread policy.  Self-contained (no spdep): a small
## row-standardised grid listw is built by hand.

library(fastLISA)

grid_nb <- function(nr, nc) {
  idx <- function(r, c) (r - 1L) * nc + c
  nb  <- vector("list", nr * nc)
  for (r in seq_len(nr)) for (cc in seq_len(nc)) {
    nbrs <- integer(0)
    if (r > 1L)  nbrs <- c(nbrs, idx(r - 1L, cc))
    if (r < nr)  nbrs <- c(nbrs, idx(r + 1L, cc))
    if (cc > 1L) nbrs <- c(nbrs, idx(r, cc - 1L))
    if (cc < nc) nbrs <- c(nbrs, idx(r, cc + 1L))
    nb[[idx(r, cc)]] <- sort(as.integer(nbrs))
  }
  nb
}

make_listw <- function(nb) {
  w <- lapply(nb, function(z)
    if (length(z) == 1L && z[1] == 0L) 0 else rep(1 / length(z), length(z)))
  structure(list(style = "W", neighbours = nb, weights = w), class = "listw")
}

make_grid_listw <- function(nr, nc) make_listw(grid_nb(nr, nc))

lw <- make_grid_listw(4, 4)
x  <- seq_len(16)
y  <- rev(x)

functions <- list(
  function(nc) local_moran(x, lw, nsim = 99, iseed = 7, n.cores = nc),
  function(nc) local_moran_bv(x, y, lw, nsim = 99, iseed = 7, n.cores = nc),
  function(nc) local_g(x, lw, nsim = 99, iseed = 7, n.cores = nc),
  function(nc) local_gstar(x, lw, nsim = 99, iseed = 7, n.cores = nc),
  function(nc) local_geary(x, lw, nsim = 99, iseed = 7, n.cores = nc),
  function(nc) local_multigeary(cbind(x, y), lw, nsim = 99, iseed = 7, n.cores = nc)
)

for (fun in functions) {
  one <- fun(1L)
  two <- fun(2L)
  stopifnot(identical(one, two))
  stopifnot(identical(two, fun(2L)))
}

cat("multicore-reproducible: OK\n")
