test_that("seeded results are identical across thread counts", {
  skip_if_not_installed("spdep")

  lw <- spdep::nb2listw(spdep::cell2nb(4, 4))
  x <- seq_len(16)
  y <- rev(x)

  # Each statistic, run single-threaded vs two-threaded with the same seed.
  # Per-observation RNG seeding must make the result thread-independent, so
  # the two runs (and a repeat) must be bit-identical.  Two cores keeps the
  # test within CRAN's example/test thread policy.
  functions <- list(
    function(nc) local_moran(x, lw, nsim = 99, iseed = 7, n.cores = nc),
    function(nc) local_moran_bv(x, y, lw, nsim = 99, iseed = 7, n.cores = nc),
    function(nc) local_g(x, lw, nsim = 99, iseed = 7, n.cores = nc),
    function(nc) local_gstar(x, lw, nsim = 99, iseed = 7, n.cores = nc),
    function(nc) local_geary(x, lw, nsim = 99, iseed = 7, n.cores = nc),
    function(nc) local_multigeary(cbind(x, y), lw, nsim = 99, iseed = 7,
                                  n.cores = nc)
  )

  for (fun in functions) {
    one <- fun(1L)
    two <- fun(2L)
    expect_identical(one, two)
    expect_identical(two, fun(2L))
  }
})
