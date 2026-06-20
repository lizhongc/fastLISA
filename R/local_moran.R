# local_moran(): Univariate Local Moran's I_i statistic (Anselin 1995).
#
# The classic LISA correlating a value with its neighbours' average. On the
# sample (n-1) standardised variable z, with row-standardised weights w*_ij:
#
#   I_i = z_i * sum_j w*_ij z_j = z_i * lag(z)_i
#
# Positive I_i => i is similar to its neighbours (High-High or Low-Low cluster);
# negative => spatial outlier (High-Low or Low-High). This is the x = y case of
# bivariate Moran, so it delegates to local_moran_bv (scale = TRUE) and relabels
# columns to Ii / Z.Ii (and E.Ii/Var.Ii/Skew.Ii/Kurt.Ii when moments = TRUE).
#
# Inference, folded two-tailed pseudo-p, Z-score, e1071 type-3 moments,
# NA->Undefined / neighbourless->Isolated handling, and per-observation RNG
# seeding (identical for any n.cores) are all as in local_moran_bv.
#
# Inputs: x (numeric n), listw, nsim, iseed, p.value, n.cores, moments.
# Output: numeric matrix of class c("localmoran", "matrix", "array"); columns Ii,
#   Z.Ii, "Pr(folded) Sim" (+ moment columns when moments = TRUE). Attributes:
#   quadr (Moran-scatterplot quadrants), cluster (Not significant / High-High /
#   Low-Low / Low-High / High-Low / Undefined / Isolated), call.
local_moran <- function(x, listw,
                        nsim    = 999L,
                        iseed   = NULL,
                        p.value = 0.05,
                        n.cores = 1L,
                        moments = FALSE) 
{
  # Univariate Local Moran's I on standardised variable x is mathematically equivalent
  # to Bivariate Local Moran's I with x = y and scale = TRUE.
  # We reuse local_moran_bv to avoid duplicate standardisation and permutation code.
  res_bv <- local_moran_bv(x = x, y = x, listw = listw, nsim = nsim,
                          scale = TRUE,
                          iseed = iseed,
                          p.value = p.value,
                          n.cores = n.cores,
                          moments = moments)

  res_mat <- unclass(res_bv)
  colnames(res_mat)[1:2] <- c("Ii", "Z.Ii")
  if (moments) {
    colnames(res_mat)[4:7] <- c("E.Ii", "Var.Ii", "Skew.Ii", "Kurt.Ii")
  }

  # Copy attributes from local_moran_bv result
  attr(res_mat, "quadr")   <- attr(res_bv, "quadr")
  attr(res_mat, "cluster") <- attr(res_bv, "cluster")
  attr(res_mat, "call")    <- match.call()

  class(res_mat) <- c("localmoran", "matrix", "array")
  res_mat
}
