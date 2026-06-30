## Tests for the p.method = c("count", "rank") option.
##
## "count" (default) must be bit-identical to the historical behaviour; "rank"
## must return the spdep ties-averaged rank, folded to the smaller tail.
## Self-contained (no spdep): a small row-standardised grid listw is built here.

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

make_grid_listw <- function(nr, nc) {
  nb <- grid_nb(nr, nc)
  w  <- lapply(nb, function(z) rep(1 / length(z), length(z)))
  structure(list(style = "W", neighbours = nb, weights = w), class = "listw")
}

strip_call <- function(z) { attr(z, "call") <- NULL; z }

## --- p.method default is "count" and bit-identical to explicit count ---------
lw <- make_grid_listw(5, 5)
x  <- as.numeric(seq_len(25))
stopifnot(identical(
  strip_call(local_moran(x, lw, nsim = 199, iseed = 7)),
  strip_call(local_moran(x, lw, nsim = 199, iseed = 7, p.method = "count"))
))
stopifnot(identical(
  strip_call(local_g(x, lw, nsim = 199, iseed = 7)),
  strip_call(local_g(x, lw, nsim = 199, iseed = 7, p.method = "count"))
))

## --- rank p-values are folded, share the NA pattern, and differ from count ---
lw   <- make_grid_listw(6, 6)
x    <- as.numeric(seq_len(36))
x[c(4, 19)] <- NA
nsim <- 499L
funs <- list(
  moran      = function(m) local_moran(x, lw, nsim = nsim, iseed = 11, p.method = m),
  geary      = function(m) local_geary(x, lw, nsim = nsim, iseed = 11, p.method = m),
  g          = function(m) local_g(x, lw, nsim = nsim, iseed = 11, p.method = m),
  gstar      = function(m) local_gstar(x, lw, nsim = nsim, iseed = 11, p.method = m),
  multigeary = function(m) local_multigeary(cbind(x, rev(x)), lw, nsim = nsim,
                                            iseed = 11, p.method = m)
)
for (f in funs) {
  pc <- f("count")[, 3]
  pr <- f("rank")[, 3]
  stopifnot(identical(is.na(pc), is.na(pr)))
  stopifnot(max(pr, na.rm = TRUE) <= 0.5 + 1e-12)
  stopifnot(min(pr, na.rm = TRUE) >= 1 / (nsim + 1) - 1e-12)
  stopifnot(any(abs(pc - pr) > 1e-9, na.rm = TRUE))
}

## --- rank p-values are reproducible across thread counts ---------------------
lw <- make_grid_listw(4, 4)
x  <- seq_len(16)
y  <- rev(x)
funs <- list(
  function(nc) local_moran(x, lw, nsim = 99, iseed = 3, n.cores = nc, p.method = "rank"),
  function(nc) local_g(x, lw, nsim = 99, iseed = 3, n.cores = nc, p.method = "rank"),
  function(nc) local_geary(x, lw, nsim = 99, iseed = 3, n.cores = nc, p.method = "rank"),
  function(nc) local_multigeary(cbind(x, y), lw, nsim = 99, iseed = 3,
                                n.cores = nc, p.method = "rank")
)
for (fun in funs) stopifnot(identical(fun(1L), fun(2L)))

## --- C rank index reproduces R's ties-averaged rank including ties -----------
## The C kernel computes ri = as.integer(n_less + (n_eq + 2)/2); this must equal
## as.integer(rank(c(perm, obs), ties.method = "average")[nsim + 1]) for any input.
## This is a pure-R identity check (no spdep, no package call).
c_index <- function(perm, obs) {
  n_less <- sum(perm < obs)
  n_eq   <- sum(perm == obs)
  as.integer(n_less + (n_eq + 2) / 2)
}
r_index <- function(perm, obs) as.integer(rank(c(perm, obs))[length(perm) + 1L])
set.seed(99)
for (trial in seq_len(200)) {
  nsim_t <- sample(5:50, 1)
  perm   <- sample(c(rnorm(nsim_t), rep(0, nsim_t)), nsim_t, replace = TRUE)
  obs    <- sample(c(perm, rnorm(1), 0), 1)
  stopifnot(identical(c_index(perm, obs), r_index(perm, obs)))
}

## --- invalid p.method is rejected via match.arg ------------------------------
lw <- make_grid_listw(4, 4)
x  <- as.numeric(seq_len(16))
stopifnot(inherits(try(local_moran(x, lw, p.method = "folded"), silent = TRUE), "try-error"))
stopifnot(inherits(try(local_g(x, lw, p.method = "zzz"),        silent = TRUE), "try-error"))

cat("p-method: OK\n")
