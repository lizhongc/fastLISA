# local_geary(): Univariate Local Geary's C_i statistic (Anselin 1995).
#
# A squared-difference LISA measuring local dissimilarity between a unit and its
# neighbours. On the sample (n-1) standardised variable z (when scale = TRUE),
# with row-standardised weights w*_ij:
#
#   C_i = sum_j w*_ij (z_i - z_j)^2 = z_i^2 - 2 z_i * lag(z)_i + lag(z^2)_i
#
# Small C_i  => i resembles its neighbours (positive spatial association);
# large C_i  => i differs from them (negative association / spatial outlier).
#
# Inference: one-tailed conditional permutation (nsim reps). The observed C_i is
# compared with the permutation mean to choose the tail, then
#   Pr = (t + 1) / (nsim + 1),  t = #{permuted C in the chosen tail}.
# Z.Ci = (C_i - E_perm) / sqrt(Var_perm); Skew.Ci/Kurt.Ci (moments = TRUE) use the
# e1071 type-3 convention.
#
# Reproducibility/threads: C backend re-seeds per observation -> identical for any
# n.cores; ignored without OpenMP.
#
# Inputs:
#   x        numeric vector, length n.
#   listw    spdep 'listw' spatial weights.
#   nsim     integer permutations (default 999L; >= 1).
#   scale    logical; if TRUE (default) z-score-standardise x in R before the test.
#   iseed    integer RNG seed, or NULL.
#   p.value  cluster significance cutoff (default 0.05).
#   n.cores  number of OpenMP threads (default 1L).
#   moments  logical; if TRUE append E.Ci, Var.Ci, Skew.Ci, Kurt.Ci.
#
# Output: numeric matrix of class c("localC", "matrix", "array"); columns Ci,
#   Z.Ci, "Pr Sim" (+ moment columns when moments = TRUE). Attributes: cluster
#   (factor, levels Not significant / High-High / Low-Low / Other Positive /
#   Negative / Undefined / Isolated) and call.
local_geary <- function(x, listw,
                        nsim    = 999L,
                        scale   = TRUE,
                        iseed   = NULL,
                        p.value = 0.05,
                        n.cores = 1L,
                        moments = FALSE,
                        p.method = c("count", "rank")) {

  ## ------------------------------------------------------------------ ##
  ## 1. Input validation                                                ##
  ## ------------------------------------------------------------------ ##
  if (!is.numeric(x))
    stop("'x' must be a numeric vector.")
  if (!inherits(listw, "listw"))
    stop("'listw' must be a spdep listw object.")

  n <- length(listw$neighbours)
  if (length(x) != n)
    stop("Length of 'x' does not match number of observations in 'listw'.")

  n.cores  <- max(1L,as.integer(n.cores))
  p.value  <- as.double(p.value)
  p.method <- match.arg(p.method)
  seed     <- if (is.null(iseed)) 123456789.0 else as.double(iseed)
  nsim     <- as.integer(nsim)
  if (nsim < 1L) {
    stop("nsim must be at least 1.")
  }

  ## ------------------------------------------------------------------ ##
  ## 2. Identify NAs and Isolates                                       ##
  ## ------------------------------------------------------------------ ##
  isolate_mask <- vapply(listw$neighbours, function(nbrs) {
    length(nbrs) == 0L || (length(nbrs) == 1L && nbrs[1] == 0L)
  }, logical(1L))
  na_mask      <- is.na(x)
  undef_mask   <- na_mask | isolate_mask
  undef_flag   <- as.integer(na_mask)

  n_valid <- sum(!na_mask)
  if (n_valid < 3L) {
    stop("Too few valid observations (need at least 3).")
  }

  ## ------------------------------------------------------------------ ##
  ## 3. Call C backend (nsim >= 1 required)                            ##
  ## ------------------------------------------------------------------ ##
  x_in <- .prepare_lisa_data(x, na_mask, scale)

  csr <- .listw2csr(listw)
  raw <- .Call(
    "r_localgeary",
    csr$row_ptr,
    csr$col_idx,
    csr$weights,
    x_in,
    undef_flag,
    nsim,
    seed,
    n.cores,
    p.value,
    p.method == "rank"
  )
  obs_geary <- raw$geary
  p_sim     <- raw$p_val
  p_sim[undef_mask] <- NA_real_

  E.Ci   <- raw$mean
  E.Ci[undef_mask]   <- NA_real_
  Var.Ci <- raw$var
  Var.Ci[undef_mask]   <- NA_real_
  Skew.Ci <- raw$skew
  Skew.Ci[undef_mask] <- NA_real_
  Kurt.Ci <- raw$kurt
  Kurt.Ci[undef_mask] <- NA_real_

  # Z-score based on permutation moments
  Z.Ci <- (obs_geary - E.Ci) / sqrt(Var.Ci)
  Z.Ci[undef_mask] <- NA_real_

  ## ------------------------------------------------------------------ ##
  ## 4. Build Result Matrix and Cluster Factor                           ##
  ## ------------------------------------------------------------------ ##
  if (moments) {
    res <- matrix(
      c(obs_geary, Z.Ci, p_sim, E.Ci, Var.Ci, Skew.Ci, Kurt.Ci),
      nrow = n,
      dimnames = list(attr(listw, "region.id"),
                      c("Ci", "Z.Ci", "Pr Sim",
                        "E.Ci", "Var.Ci", "Skew.Ci", "Kurt.Ci"))
    )
  } else {
    res <- matrix(
      c(obs_geary, Z.Ci, p_sim),
      nrow = n,
      dimnames = list(attr(listw, "region.id"),
                      c("Ci", "Z.Ci", "Pr Sim"))
    )
  }

  cluster <- factor(raw$cluster, levels = 0:6,
                    labels = c("Not significant", "High-High", "Low-Low",
                               "Other Positive", "Negative", "Undefined",
                               "Isolated"))

  ## ------------------------------------------------------------------ ##
  ## 5. Return Value Assembly                                           ##
  ## ------------------------------------------------------------------ ##
  attr(res, "cluster")    <- cluster
  attr(res, "call")       <- match.call()

  class(res) <- c("localC", "matrix", "array")
  res
}
