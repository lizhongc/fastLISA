# local_multigeary(): Multivariate Local Geary's C_i statistic (Anselin 2019).
#
# Extends local Geary's C to K variables: the average over variables of the
# univariate squared-difference statistic. On sample (n-1) standardised variables
# z^1..z^K (when scale = TRUE), with row-standardised weights w*_ij:
#
#   C_i = (1/K) * sum_{v=1}^{K} sum_j w*_ij (z^v_i - z^v_j)^2
#
# Small C_i  => i is multivariate-similar to its neighbours (positive
# association); large C_i  => multivariate dissimilarity / spatial outlier.
#
# Inference: one-tailed conditional permutation (nsim reps); each replicate
# applies the SAME permuted neighbour configuration across all K variables.
#   Pr = (t + 1) / (nsim + 1),  t = #{permuted C in the chosen tail}.
# Z.Ci = (C_i - E_perm) / sqrt(Var_perm); Skew.Ci/Kurt.Ci (moments = TRUE) use the
# e1071 type-3 convention.
#
# Reproducibility/threads: C backend re-seeds per observation -> identical for any
# n.cores; ignored without OpenMP.
#
# Inputs:
#   df       data.frame or matrix, n rows x K columns (one column per variable; an
#            'sf' geometry column is dropped). A row with any NA is Undefined.
#   listw    spdep 'listw' spatial weights.
#   nsim     integer permutations (default 999L; >= 1).
#   scale    logical; if TRUE (default) z-score-standardise each column in R.
#   iseed    integer RNG seed, or NULL.
#   p.value  cluster significance cutoff (default 0.05).
#   n.cores  number of OpenMP threads (default 1L).
#   moments  logical; if TRUE append E.Ci, Var.Ci, Skew.Ci, Kurt.Ci.
#   p.method pseudo p-value rule: "count" (default, standard) or "rank"
#            (spdep ties-averaged). The two differ only under exact ties.
#
# Output: numeric matrix of class c("localC", "matrix", "array"); columns Ci,
#   Z.Ci, "Pr Sim" (+ moment columns when moments = TRUE). Attributes: cluster
#   (factor, levels Not significant / Positive / Negative / Undefined / Isolated)
#   and call.
local_multigeary <- function(df, listw,
                             nsim    = 999L,
                             scale   = TRUE,
                             iseed   = NULL,
                             p.value = 0.05,
                             n.cores = 1L,
                             moments = FALSE,
                             p.method = c("count", "rank"))
{
  ## ------------------------------------------------------------------ ##
  ## 1. Input validation                                                ##
  ## ------------------------------------------------------------------ ##
  if (!is.data.frame(df) && !is.matrix(df))
    stop("'df' must be a data.frame or matrix.")
  if (inherits(df, "sf")) {
    df[[attr(df, "sf_column")]] <- NULL
    df <- as.data.frame(df)
  }
  if (!inherits(listw, "listw"))
    stop("'listw' must be a spdep listw object.")

  n <- length(listw$neighbours)
  if (nrow(df) != n)
    stop("Number of rows in 'df' does not match number of observations in 'listw'.")

  n.cores  <- max(1L,as.integer(n.cores))
  p.value  <- as.double(p.value)
  p.method <- match.arg(p.method)
  seed     <- if (is.null(iseed)) 123456789.0 else as.double(iseed)
  nsim     <- as.integer(nsim)
  if (nsim < 1L) {
    stop("nsim must be at least 1.")
  }

  # Convert to matrix of doubles
  xorig <- as.matrix(df)
  nc    <- ncol(xorig)
  if (nc < 1L)
    stop("'df' must have at least 1 column.")

  ## ------------------------------------------------------------------ ##
  ## 2. Identify NAs and Isolates                                       ##
  ## ------------------------------------------------------------------ ##
  isolate_mask <- vapply(listw$neighbours, function(nbrs) {
    length(nbrs) == 0L || (length(nbrs) == 1L && nbrs[1] == 0L)
  }, logical(1L))
  na_mask      <- rowSums(is.na(xorig)) > 0L
  undef_mask   <- na_mask | isolate_mask
  undef_flag   <- as.integer(na_mask)

  n_valid <- sum(!na_mask)
  if (n_valid < 3L) {
    stop("Too few valid observations (need at least 3).")
  }

  ## ------------------------------------------------------------------ ##
  ## 3. Call C backend (nsim >= 1 required)                            ##
  ## ------------------------------------------------------------------ ##
  # Prepare data_list: list of numeric vectors
  data_list <- lapply(seq_len(nc), function(j) {
    .prepare_lisa_data(xorig[, j], na_mask, scale)
  })

  csr <- .listw2csr(listw)
  raw <- .Call(
    "r_localmultigeary",
    csr$row_ptr,
    csr$col_idx,
    csr$weights,
    data_list,
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
  E.Ci[undef_mask]    <- NA_real_
  Var.Ci <- raw$var
  Var.Ci[undef_mask]  <- NA_real_
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

  cluster <- factor(raw$cluster, levels = 0:4,
                    labels = c("Not significant", "Positive", "Negative",
                               "Undefined", "Isolated"))

  ## ------------------------------------------------------------------ ##
  ## 5. Return Value Assembly                                           ##
  ## ------------------------------------------------------------------ ##
  attr(res, "cluster")     <- cluster
  attr(res, "call")        <- match.call()

  class(res) <- c("localC", "matrix", "array")
  res
}
