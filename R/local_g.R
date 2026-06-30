# local_g(): Local Getis-Ord G_i statistic (Getis & Ord 1992; Ord & Getis 1995).
#
# Detects local clustering of high values ("hot spots") and low values ("cold
# spots"). For observation i, with row-standardised spatial weights w*_ij and the
# focal value EXCLUDED from both the lag and the denominator:
#
#   G_i = ( sum_{j != i} w*_ij x_j ) / ( sum_k x_k - x_i )
#
# A large G_i means i is surrounded by high values; a small G_i a low-value
# neighbourhood. G_i has no self term (see local_gstar for the self-inclusive G*).
#
# Inference: conditional permutation (nsim reps) with the focal x_i held fixed
# while neighbour values are permuted. Folded two-tailed pseudo p-value
#   Pr = (min(g, nsim - g) + 1) / (nsim + 1),  g = #{G_perm >= G_obs}.
# Z.Gi = (G_i - E_perm) / sqrt(Var_perm) from the permutation moments;
# Skew.Gi/Kurt.Gi (moments = TRUE) use the e1071 type-3 convention.
#
# Reproducibility/threads: the C backend re-seeds its RNG per observation, so
# results are identical for any n.cores; n.cores is ignored without OpenMP.
#
# Inputs:
#   x        numeric vector, length n (one value per spatial unit).
#   listw    spdep 'listw' spatial weights (any style).
#   nsim     integer number of permutations (default 999L; must be >= 1).
#   iseed    integer RNG seed, or NULL for the package default.
#   p.value  significance cutoff used to filter the cluster factor (default 0.05).
#   n.cores  number of OpenMP threads (default 1L).
#   moments  logical; if TRUE append E.Gi, Var.Gi, Skew.Gi, Kurt.Gi.
#
# Output: numeric matrix of class c("localG", "matrix", "array"), n rows, columns
#   Gi, Z.Gi, "Pr(folded) Sim" (plus the four moment columns when moments = TRUE).
#   Attributes:
#     cluster  factor, levels Not significant / High-High / Low-Low / Undefined /
#              Isolated (NA obs -> Undefined; neighbourless obs -> Isolated).
#     gstari   logical FALSE (this is G, not G*).
#     call     the matched call.
local_g <- function(x, listw,
                    nsim    = 999L,
                    iseed   = NULL,
                    p.value = 0.05,
                    n.cores = 1L,
                    moments = FALSE,
                    p.method = c("count", "rank"))
{
  ## ------------------------------------------------------------------ ##
  ## 1. Input validation                                                ##
  ## ------------------------------------------------------------------ ##
  if (!is.numeric(x)) {
    stop("'x' must be a numeric vector.")
  }
  if (!inherits(listw, "listw")) {
    stop("'listw' must be a spdep listw object.")
  }

  n <- length(listw$neighbours)
  if (length(x) != n) {
    stop("Length of 'x' does not match number of observations in 'listw'.")
  }

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
  na_mask    <- is.na(x)
  undef_mask <- na_mask | isolate_mask
  undef_flag <- as.integer(na_mask)

  n_valid <- sum(!na_mask)
  if (n_valid < 3L) {
    stop("Too few valid observations (need at least 3).")
  }

  ## ------------------------------------------------------------------ ##
  ## 3. Call C backend (nsim >= 1 required)                            ##
  ## ------------------------------------------------------------------ ##
  x_in <- as.double(x)
  x_in[na_mask] <- 0.0

  csr <- .listw2csr(listw)
  raw <- .Call(
    "r_localg",
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
  obs_g    <- raw$g
  p_folded <- raw$p_val
  p_folded[undef_mask] <- NA_real_

  EG_sim <- raw$mean
  EG_sim[undef_mask] <- NA_real_
  VG_sim <- raw$var
  VG_sim[undef_mask] <- NA_real_
  Skew.Gi <- raw$skew
  Skew.Gi[undef_mask] <- NA_real_
  Kurt.Gi <- raw$kurt
  Kurt.Gi[undef_mask] <- NA_real_

  # Standardised Z-score based on permutation moments
  Z.Gi <- (obs_g - EG_sim) / sqrt(VG_sim)
  Z.Gi[undef_mask] <- NA_real_

  ## ------------------------------------------------------------------ ##
  ## 4. Build Result Matrix and Cluster Factor                           ##
  ## ------------------------------------------------------------------ ##
  if (moments) {
    res <- matrix(
      c(obs_g, Z.Gi, p_folded, EG_sim, VG_sim, Skew.Gi, Kurt.Gi),
      nrow = n,
      dimnames = list(attr(listw, "region.id"),
                      c("Gi", "Z.Gi", "Pr(folded) Sim",
                        "E.Gi", "Var.Gi", "Skew.Gi", "Kurt.Gi"))
    )
  } else {
    res <- matrix(
      c(obs_g, Z.Gi, p_folded),
      nrow = n,
      dimnames = list(attr(listw, "region.id"),
                      c("Gi", "Z.Gi", "Pr(folded) Sim"))
    )
  }

  cluster <- factor(raw$cluster, levels = 0:4,
                    labels = c("Not significant", "High-High", "Low-Low",
                               "Undefined", "Isolated"))

  ## ------------------------------------------------------------------ ##
  ## 5. Return Value Assembly                                           ##
  ## ------------------------------------------------------------------ ##
  attr(res, "cluster")   <- cluster
  attr(res, "gstari")    <- FALSE
  attr(res, "call")      <- match.call()

  class(res) <- c("localG", "matrix", "array")
  res
}
