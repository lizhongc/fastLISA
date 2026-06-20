# local_moran_bv(): Bivariate Local Moran's I_bv,i statistic (Anselin 1995).
#
# Correlates a variable x at i with the spatial lag of a second variable y over
# i's neighbours. On sample (n-1) standardised z_x, z_y (when scale = TRUE), with
# row-standardised weights w*_ij:
#
#   I_bv,i = z_x,i * sum_j w*_ij z_y,j = z_x,i * lag(z_y)_i
#
# Positive I_bv,i => i's x value coincides with high lagged y nearby; negative =>
# spatial mismatch. Univariate local Moran is the x = y special case (see
# local_moran). The listw is converted to CSR internally and passed to the C
# backend; the result is a matrix of class c("localmoran", "matrix", "array").
#
# Inference: conditional permutation (nsim reps) of the neighbour y-values, focal
# fixed. Folded two-tailed pseudo-p
#   Pr = (min(g, nsim - g) + 1) / (nsim + 1),  g = #{I_perm >= I_obs}.
# Z.Ibvi = (I_bv,i - E_perm) / sqrt(Var_perm); Skew/Kurt (moments = TRUE) use the
# e1071 type-3 convention. Backend re-seeds per observation -> identical for any
# n.cores.
#
# Inputs:
#   x, y     numeric vectors of length n (y is the lagged variable).
#   listw    spdep 'listw' spatial weights (any style).
#   nsim     integer permutations (default 999L; >= 1).
#   scale    logical; if TRUE (default) sample-standardise x and y in R.
#   iseed    integer RNG seed, or NULL.
#   p.value  cluster significance cutoff (default 0.05).
#   n.cores  number of OpenMP threads (default 1L).
#   moments  logical; if TRUE append E.Ibvi, Var.Ibvi, Skew.Ibvi, Kurt.Ibvi.
#
# Output: numeric matrix of class c("localmoran", "matrix", "array"); columns
#   Ibvi, Z.Ibvi, "Pr(folded) Sim" (+ moment columns when moments = TRUE).
#   Attributes: quadr (data.frame of mean/median/pysal Moran-scatterplot
#   quadrants), cluster (factor: Not significant / High-High / Low-Low /
#   Low-High / High-Low / Undefined / Isolated), call.
#
# Internal helpers defined below: .listw_lag (spatial lag), .scale_sample (n-1
# z-score), .quadrant_factor (Moran-scatterplot quadrant factor).

local_moran_bv <- function(x, y, listw,
                           nsim    = 999L,
                           scale   = TRUE,
                           iseed   = NULL,
                           p.value = 0.05,
                           n.cores = 1L,
                           moments = FALSE)
{
  ## ------------------------------------------------------------------ ##
  ## 1.  Input validation                                                 ##
  ## ------------------------------------------------------------------ ##
  if (!is.numeric(x))
    stop("'x' must be a numeric vector.")
  if (!is.numeric(y))
    stop("'y' must be a numeric vector.")
  if (!inherits(listw, "listw"))
    stop("'listw' must be a spdep listw object.")

  n <- length(listw$neighbours)
  if (length(x) != n)
    stop("Length of 'x' does not match number of observations in 'listw'.")
  if (length(y) != n)
    stop("Length of 'y' does not match number of observations in 'listw'.")

  n.cores <- max(1L, as.integer(n.cores))
  p.value <- as.double(p.value)
  seed    <- if (is.null(iseed)) 123456789.0 else as.double(iseed)
  nsim    <- as.integer(nsim)
  if (nsim < 1L) {
    stop("nsim must be at least 1.")
  }

  ## ------------------------------------------------------------------ ##
  ## 2.  Build quadrant attribute (on original data scale, before scaling) ##
  ## ------------------------------------------------------------------ ##
  lbs <- c("Low", "High")
  ly  <- .listw_lag(y, listw)

  xx    <- mean(x,  na.rm = TRUE)
  lyy   <- mean(ly, na.rm = TRUE)
  xmed  <- median(x,  na.rm = TRUE)
  lymed <- median(ly, na.rm = TRUE)

  quadr_mean   <- .quadrant_factor(x, ly, xx,   lyy,   lbs)
  quadr_median <- .quadrant_factor(x, ly, xmed, lymed, lbs)

  ## ------------------------------------------------------------------ ##
  ## 3.  Identify Undefined and Isolates                                  ##
  ## ------------------------------------------------------------------ ##
  isolate_mask <- vapply(listw$neighbours, function(nbrs) {
    length(nbrs) == 0L || (length(nbrs) == 1L && nbrs[1] == 0L)
  }, logical(1L))
  na_mask      <- is.na(x) | is.na(y)
  undef_mask   <- na_mask | isolate_mask
  undef_flag   <- as.integer(na_mask)

  n_valid <- sum(!na_mask)
  if (n_valid < 3L) {
    stop("Too few valid observations (need at least 3).")
  }

  x_in <- .prepare_lisa_data(x, na_mask, scale)
  y_in <- .prepare_lisa_data(y, na_mask, scale)

  x_sc <- x_in
  y_sc <- y_in 
  x_sc[na_mask] <- NA_real_
  y_sc[na_mask] <- NA_real_
  ly_sc <- .listw_lag(y_sc, listw)
  quadr_pysal   <- .quadrant_factor(x_sc, ly_sc, 0, 0, lbs)

  quadr <- data.frame(
    mean   = quadr_mean,
    median = quadr_median,
    pysal  = quadr_pysal,
    stringsAsFactors = TRUE
  )

  ## ------------------------------------------------------------------ ##
  ## 4.  Call C backend (nsim >= 1 required)                            ##
  ## ------------------------------------------------------------------ ##

  csr <- .listw2csr(listw)
  raw <- .Call(
    "r_bi_localmoran",
    csr$row_ptr,
    csr$col_idx,
    csr$weights,
    x_in,
    y_in,
    undef_flag,
    nsim,
    seed,
    n.cores,
    p.value
  )
  obs_Ibvi     <- raw$bimoran
  p_folded_sim <- raw$p_val
  p_folded_sim[undef_mask] <- NA_real_
  cluster_int  <- raw$cluster

  # Empirical moments from permutation test
  E.Ibvi    <- raw$mean
  E.Ibvi[undef_mask]    <- NA_real_
  Var.Ibvi  <- raw$var
  Var.Ibvi[undef_mask]  <- NA_real_
  Skew.Ibvi <- raw$skew
  Skew.Ibvi[undef_mask] <- NA_real_
  Kurt.Ibvi <- raw$kurt
  Kurt.Ibvi[undef_mask] <- NA_real_

  Z.Ibvi <- (obs_Ibvi - E.Ibvi) / sqrt(Var.Ibvi)
  Z.Ibvi[undef_mask] <- NA_real_

  ## ------------------------------------------------------------------ ##
  ## 5.  Assemble Output Matrix                                           ##
  ## ------------------------------------------------------------------ ##
  if (moments) {
    res <- matrix(
      c(obs_Ibvi, Z.Ibvi, p_folded_sim, E.Ibvi, Var.Ibvi, Skew.Ibvi, Kurt.Ibvi),
      nrow  = n,
      ncol  = 7L,
      dimnames = list(
        attr(listw, "region.id"),
        c("Ibvi", "Z.Ibvi", "Pr(folded) Sim", "E.Ibvi", "Var.Ibvi", "Skew.Ibvi", "Kurt.Ibvi")
      )
    )
  } else {
    res <- matrix(
      c(obs_Ibvi, Z.Ibvi, p_folded_sim),
      nrow  = n,
      ncol  = 3L,
      dimnames = list(
        attr(listw, "region.id"),
        c("Ibvi", "Z.Ibvi", "Pr(folded) Sim")
      )
    )
  }

  ## ------------------------------------------------------------------ ##
  ## 6.  Attach attributes                                               ##
  ## ------------------------------------------------------------------ ##
  cluster <- factor(
    cluster_int,
    levels = 0:6,
    labels = c(
      "Not significant",
      "High-High",
      "Low-Low",
      "Low-High",
      "High-Low",
      "Undefined",
      "Isolated"
    )
  )

  attr(res, "quadr")   <- quadr
  attr(res, "cluster") <- cluster
  attr(res, "call")    <- match.call()

  class(res) <- c("localmoran", "matrix", "array")
  res
}


## -----------------------------------------------------------------------
## Internal helpers (not exported)
## -----------------------------------------------------------------------

# .listw_lag(x, listw): spatial lag sum_j w_ij x_j over i's neighbours using the
# raw (non-standardised) listw weights; returns NA for observations with no
# neighbours. Used to build the Moran scatter-plot quadrants.
.listw_lag <- function(x, listw) {
  vapply(seq_along(listw$neighbours), function(i) {
    nbrs <- listw$neighbours[[i]]
    wts  <- listw$weights[[i]]
    if (length(nbrs) == 0L || (length(nbrs) == 1L && nbrs[1] == 0L)) return(NA_real_)
    sum(wts * x[nbrs])
  }, numeric(1L))
}

# .scale_sample(x): sample (n-1) z-score standardisation (x - mean)/sd, matching
# the C backend; a zero or NA sd is treated as 1 to avoid division by zero.
.scale_sample <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x,   na.rm = TRUE)   # sd() uses n-1 by default in R
  if (is.na(s) || s == 0) s <- 1.0
  (x - m) / s
}

# .quadrant_factor(x, ly, cx, cly, lbs): Moran scatter-plot quadrant as an
# interaction factor of (x vs centre cx) by (lag ly vs centre cly), with the two
# level labels in lbs (e.g. Low/High) -> "Low-Low", "High-Low", etc.
.quadrant_factor <- function(x, ly, cx, cly, lbs) {
  interaction(
    cut(x,  c(-Inf, cx,  Inf), labels = lbs),
    cut(ly, c(-Inf, cly, Inf), labels = lbs),
    sep = "-"
  )
}
