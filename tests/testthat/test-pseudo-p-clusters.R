test_that("public LISA APIs expose renamed arguments and pseudo-p values", {
  skip_if_not_installed("spdep")

  lw <- spdep::nb2listw(spdep::cell2nb(3, 3))
  x <- seq_len(9)
  y <- rev(x)

  functions <- list(local_moran_bv, local_moran, local_moran_eb, local_g,
                    local_gstar, local_geary, local_multigeary)
  for (fun in functions) {
    args <- names(formals(fun))
    expect_true(all(c("p.value", "n.cores") %in% args))
    expect_false(any(c("alternative", "significance_cutoff", "cpu_threads") %in% args))
  }

  bv <- local_moran_bv(x, y, lw, nsim = 19, n.cores = 1)
  expect_equal(colnames(bv),
               c("Ibvi", "Z.Ibvi", "Pr(folded) Sim"))
  expect_true(all(bv[, "Pr(folded) Sim"] <= 0.5, na.rm = TRUE))
  expect_equal(colnames(local_moran(x, lw, nsim = 19, n.cores = 1)),
               c("Ii", "Z.Ii", "Pr(folded) Sim"))
  expect_equal(colnames(local_moran_eb(x, x + 20, lw, nsim = 19, n.cores = 1)),
               c("Ii", "Z.Ii", "Pr(folded) Sim"))
  expect_equal(colnames(local_g(x, lw, nsim = 19, n.cores = 1)),
               c("Gi", "Z.Gi", "Pr(folded) Sim"))
  expect_equal(colnames(local_gstar(x, lw, nsim = 19, n.cores = 1)),
               c("G*i", "Z.G*i", "Pr(folded) Sim"))
  expect_equal(colnames(local_geary(x, lw, nsim = 19, n.cores = 1)),
               c("Ci", "Z.Ci", "Pr Sim"))
  expect_equal(colnames(local_multigeary(cbind(x, y), lw, nsim = 19, n.cores = 1)),
               c("Ci", "Z.Ci", "Pr Sim"))

  expect_error(local_moran(x, lw, significance_cutoff = 0.05), "unused argument")
  expect_error(local_moran(x, lw, cpu_threads = 1), "unused argument")
})

test_that("CSR conversion helper is internal", {
  skip_if_not_installed("spdep")

  lw <- spdep::nb2listw(spdep::cell2nb(3, 3))
  csr <- fastLISA:::.listw2csr(lw)

  expect_named(csr, c("row_ptr", "col_idx", "weights"))
  expect_false(exists("listw2csr", envir = asNamespace("fastLISA"),
                      inherits = FALSE))
  expect_false("listw2csr" %in% getNamespaceExports("fastLISA"))
})

test_that("remaining LISA outputs use optional matrix moment columns", {
  skip_if_not_installed("spdep")

  for (fun in list(local_g, local_gstar, local_geary, local_multigeary)) {
    expect_identical(formals(fun)$moments, FALSE)
  }

  lw <- spdep::nb2listw(spdep::cell2nb(3, 3))
  x <- seq_len(9)
  y <- rev(x)
  functions <- list(
    G = function(moments) local_g(x, lw, nsim = 19, iseed = 7, n.cores = 1,
                                  moments = moments),
    Gstar = function(moments) local_gstar(x, lw, nsim = 19, iseed = 7,
                                          n.cores = 1, moments = moments),
    Geary = function(moments) local_geary(x, lw, nsim = 19, iseed = 7,
                                          n.cores = 1, moments = moments),
    Multigeary = function(moments) local_multigeary(cbind(x, y), lw, nsim = 19,
                                                    iseed = 7, n.cores = 1,
                                                    moments = moments)
  )

  expected <- list(
    G = c("Gi", "Z.Gi", "Pr(folded) Sim", "E.Gi", "Var.Gi", "Skew.Gi", "Kurt.Gi"),
    Gstar = c("G*i", "Z.G*i", "Pr(folded) Sim",
              "E.G*i", "Var.G*i", "Skew.G*i", "Kurt.G*i"),
    Geary = c("Ci", "Z.Ci", "Pr Sim", "E.Ci", "Var.Ci", "Skew.Ci", "Kurt.Ci"),
    Multigeary = c("Ci", "Z.Ci", "Pr Sim", "E.Ci", "Var.Ci", "Skew.Ci", "Kurt.Ci")
  )

  for (name in names(functions)) {
    compact <- functions[[name]](FALSE)
    expanded <- functions[[name]](TRUE)
    expect_equal(colnames(expanded), expected[[name]])
    expect_equal(c(expanded[, 1:3]), c(compact), tolerance = 1e-12)
    expect_equal(attr(expanded, "cluster"), attr(compact, "cluster"))
    expect_identical(class(compact),
                     c(if (name %in% c("G", "Gstar")) "localG" else "localC",
                       "matrix", "array"))
    expect_equal(expanded[, 2],
                 (expanded[, 1] - expanded[, 4]) / sqrt(expanded[, 5]),
                 tolerance = 1e-12)
    expect_false(any(c("internals", "pseudo-p", "mean", "var", "skew", "kurt") %in%
                     names(attributes(compact))))
  }

  expect_identical(attr(functions$G(FALSE), "gstari"), FALSE)
  expect_identical(attr(functions$Gstar(FALSE), "gstari"), TRUE)
})

test_that("Moran moment columns are optional and statistic-specific", {
  skip_if_not_installed("spdep")

  for (fun in list(local_moran_bv, local_moran, local_moran_eb)) {
    expect_identical(formals(fun)$moments, FALSE)
    expect_false("empirical.dist" %in% names(formals(fun)))
  }

  lw <- spdep::nb2listw(spdep::cell2nb(3, 3))
  x <- seq_len(9)
  y <- rev(x)

  functions <- list(
    Bivariate = function(moments) {
      local_moran_bv(x, y, lw, nsim = 19, iseed = 7, n.cores = 1,
                     moments = moments)
    },
    Moran = function(moments) {
      local_moran(x, lw, nsim = 19, iseed = 7, n.cores = 1,
                  moments = moments)
    },
    EB = function(moments) {
      local_moran_eb(x, x + 20, lw, nsim = 19, iseed = 7, n.cores = 1,
                     moments = moments)
    }
  )

  for (fun in functions) {
    compact <- fun(FALSE)
    expanded <- fun(TRUE)

    expect_equal(ncol(compact), 3L)
    expect_equal(dim(expanded[, 1:3]), dim(compact))
    expect_equal(c(expanded[, 1:3]), c(compact), tolerance = 1e-12)
    expect_equal(attr(expanded, "cluster"), attr(compact, "cluster"))
    expect_equal(attr(expanded, "quadr"), attr(compact, "quadr"))
  }

  expect_equal(colnames(functions$Bivariate(TRUE)),
               c("Ibvi", "Z.Ibvi", "Pr(folded) Sim",
                 "E.Ibvi", "Var.Ibvi", "Skew.Ibvi", "Kurt.Ibvi"))
  expect_equal(colnames(functions$Moran(TRUE)),
               c("Ii", "Z.Ii", "Pr(folded) Sim",
                 "E.Ii", "Var.Ii", "Skew.Ii", "Kurt.Ii"))
  expect_equal(colnames(functions$EB(TRUE)),
               c("Ii", "Z.Ii", "Pr(folded) Sim",
                 "E.Ii", "Var.Ii", "Skew.Ii", "Kurt.Ii"))
  expect_false("cluster" %in% colnames(functions$EB(FALSE)))

  moran <- functions$Moran(TRUE)
  moran_bv <- local_moran_bv(x, x, lw, nsim = 19, iseed = 7, n.cores = 1,
                             moments = TRUE)
  expect_equal(c(moran[, 4:7]), c(moran_bv[, 4:7]), tolerance = 1e-12)

  eb <- .rate_standardize_eb(x, x + 20)
  eb_bv <- local_moran_bv(eb$results, eb$results, lw, nsim = 19, iseed = 7,
                          n.cores = 1, moments = TRUE)
  expect_equal(c(functions$EB(TRUE)[, 4:7]), c(eb_bv[, 4:7]), tolerance = 1e-12)

  expect_error(local_moran_bv(x, y, lw, empirical.dist = TRUE), "unused argument")
  expect_error(local_moran(x, lw, empirical.dist = TRUE), "unused argument")
  expect_error(local_moran_eb(x, x + 20, lw, empirical.dist = TRUE),
               "unused argument")
})

test_that("Moran moments preserve undefined and isolated rows", {
  skip_if_not_installed("spdep")

  nb <- structure(list(c(2L, 4L), 1L, 0L, 1L), class = "nb",
                  region.id = as.character(1:4), call = quote(user_nb),
                  type = "user", sym = FALSE)
  lw <- spdep::nb2listw(nb, zero.policy = TRUE)
  x <- c(1, NA, 3, 4)

  res <- local_moran(x, lw, nsim = 9, n.cores = 1, moments = TRUE)
  expect_true(all(is.na(res[2:3, c("Z.Ii", "Pr(folded) Sim",
                                   "E.Ii", "Var.Ii", "Skew.Ii", "Kurt.Ii")])))
  expect_equal(as.character(attr(res, "cluster"))[2:3],
               c("Undefined", "Isolated"))
})

test_that("significance cutoff filters C-computed clusters", {
  skip_if_not_installed("spdep")

  lw <- spdep::nb2listw(spdep::cell2nb(3, 3))
  x <- seq_len(9)
  y <- rev(x)

  functions <- list(
    function(cutoff) local_moran_bv(x, y, lw, nsim = 19, p.value = cutoff, n.cores = 1),
    function(cutoff) local_g(x, lw, nsim = 19, p.value = cutoff, n.cores = 1),
    function(cutoff) local_gstar(x, lw, nsim = 19, p.value = cutoff, n.cores = 1),
    function(cutoff) local_geary(x, lw, nsim = 19, p.value = cutoff, n.cores = 1),
    function(cutoff) local_multigeary(cbind(x, y), lw, nsim = 19, p.value = cutoff, n.cores = 1)
  )

  for (fun in functions) {
    expect_true(all(as.character(attr(fun(0), "cluster")) == "Not significant"))
    expect_true(any(as.character(attr(fun(1), "cluster")) != "Not significant"))
  }
})

test_that("G and Gstar never manufacture isolated clusters", {
  skip_if_not_installed("spdep")

  lw <- spdep::nb2listw(spdep::cell2nb(3, 3))
  x <- seq_len(9)

  for (i in seq_len(200)) {
    g <- local_g(x, lw, nsim = 9, p.value = 0, n.cores = 1)
    gs <- local_gstar(x, lw, nsim = 9, p.value = 0, n.cores = 1)

    expect_true(all(as.character(attr(g, "cluster")) == "Not significant"))
    expect_true(all(as.character(attr(gs, "cluster")) == "Not significant"))
  }
})

test_that("single-core seeded results are reproducible", {
  skip_if_not_installed("spdep")

  lw <- spdep::nb2listw(spdep::cell2nb(4, 4))
  x <- seq_len(16)
  y <- rev(x)

  functions <- list(
    function() local_moran(x, lw, nsim = 19, iseed = 7, n.cores = 1),
    function() local_moran_bv(x, y, lw, nsim = 19, iseed = 7, n.cores = 1),
    function() local_g(x, lw, nsim = 19, iseed = 7, n.cores = 1),
    function() local_gstar(x, lw, nsim = 19, iseed = 7, n.cores = 1),
    function() local_geary(x, lw, nsim = 19, iseed = 7, n.cores = 1),
    function() local_multigeary(cbind(x, y), lw, nsim = 19, iseed = 7,
                                n.cores = 1)
  )

  for (fun in functions) {
    expect_identical(fun(), fun())
  }
})

test_that("R-level scaling matches manually prepared data", {
  skip_if_not_installed("spdep")

  lw <- spdep::nb2listw(spdep::cell2nb(3, 3))
  x <- c(1:4, NA, 6:9)
  y <- c(9:3, NA, 1)
  valid <- !is.na(x) & !is.na(y)
  xs <- x
  ys <- y
  xs[valid] <- as.numeric(scale(x[valid]))
  ys[valid] <- as.numeric(scale(y[valid]))

  bv_scaled <- local_moran_bv(x, y, lw, nsim = 19, iseed = 7, n.cores = 1)
  bv_manual <- local_moran_bv(xs, ys, lw, nsim = 19, scale = FALSE,
                             iseed = 7, n.cores = 1)
  expect_equal(c(bv_scaled), c(bv_manual), tolerance = 1e-12)
  expect_equal(attr(bv_scaled, "cluster"), attr(bv_manual, "cluster"))

  gx <- c(1:4, NA, 6:9)
  gvalid <- !is.na(gx)
  gxs <- gx
  gxs[gvalid] <- as.numeric(scale(gx[gvalid]))
  geary_scaled <- local_geary(gx, lw, nsim = 19, iseed = 7, n.cores = 1)
  geary_manual <- local_geary(gxs, lw, nsim = 19, scale = FALSE,
                              iseed = 7, n.cores = 1)
  expect_equal(c(geary_scaled), c(geary_manual), tolerance = 1e-12)
  expect_equal(attr(geary_scaled, "cluster"), attr(geary_manual, "cluster"))

  df <- cbind(x, y)
  dfs <- cbind(xs, ys)
  mg_scaled <- local_multigeary(df, lw, nsim = 19, iseed = 7, n.cores = 1)
  mg_manual <- local_multigeary(dfs, lw, nsim = 19, scale = FALSE,
                                iseed = 7, n.cores = 1)
  expect_equal(c(mg_scaled), c(mg_manual), tolerance = 1e-12)
  expect_equal(attr(mg_scaled, "cluster"), attr(mg_manual, "cluster"))
})

test_that("scale FALSE passes unstandardised values to C", {
  skip_if_not_installed("spdep")

  lw <- spdep::nb2listw(spdep::cell2nb(3, 3))
  x <- seq(10, 90, by = 10)
  y <- seq(2, 18, by = 2)

  bv <- local_moran_bv(x, y, lw, nsim = 9, scale = FALSE, n.cores = 1)
  expected_bv <- x * spdep::lag.listw(lw, y)
  expect_equal(unname(bv[, "Ibvi"]), unname(expected_bv), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(
    bv[, "Ibvi"],
    local_moran_bv(x, y, lw, nsim = 9, scale = TRUE, n.cores = 1)[, "Ibvi"]
  )))

  geary <- local_geary(x, lw, nsim = 9, scale = FALSE, n.cores = 1)
  expected_geary <- vapply(seq_along(x), function(i) {
    sum(lw$weights[[i]] * (x[i] - x[lw$neighbours[[i]]])^2)
  }, numeric(1))
  expect_equal(unname(geary[, "Ci"]), expected_geary, tolerance = 1e-12)
})

test_that("Undefined and Isolated cluster labels are preserved", {
  skip_if_not_installed("spdep")

  nb <- structure(list(c(2L, 4L), 1L, 0L, 1L), class = "nb",
                  region.id = as.character(1:4), call = quote(user_nb),
                  type = "user", sym = FALSE)
  lw <- spdep::nb2listw(nb, zero.policy = TRUE)
  x <- c(1, NA, 3, 4)

  res <- local_geary(x, lw, nsim = 9, p.value = 0, n.cores = 1)
  expect_equal(as.character(attr(res, "cluster"))[2:3], c("Undefined", "Isolated"))

  expanded <- local_geary(x, lw, nsim = 9, p.value = 0, n.cores = 1,
                          moments = TRUE)
  expect_true(all(is.na(expanded[2:3, c("Z.Ci", "Pr Sim", "E.Ci", "Var.Ci",
                                       "Skew.Ci", "Kurt.Ci")])))

  remaining <- list(
    local_g(x, lw, nsim = 9, n.cores = 1, moments = TRUE),
    local_gstar(x, lw, nsim = 9, n.cores = 1, moments = TRUE),
    local_multigeary(cbind(x, x * 2), lw, nsim = 9, n.cores = 1,
                     moments = TRUE)
  )
  for (result in remaining) {
    expect_true(all(is.na(result[2:3, 2:7])))
    expect_equal(as.character(attr(result, "cluster"))[2:3],
                 c("Undefined", "Isolated"))
  }
})
