# local_gstar(): Local Getis-Ord G*_i statistic (Getis & Ord 1992; Ord & Getis 1995).
#
# Like local_g but INCLUDES the focal unit i itself (a self-neighbour with weight
# 1). With m_i valid neighbours, row-standardised neighbour weights w_ij, and the
# global total S = sum_k x_k:
#
#   G*_i = [ ((sum_{j in N_i} w_ij x_j) / (sum_j w_ij)) * m_i + x_i ] / (m_i + 1) / S
#
# i.e. the average value over the focal unit AND its neighbours, divided by the
# global sum. Large/small G*_i flag hot/cold spots (focal unit included).
#
# Inference, folded two-tailed pseudo-p, Z-score, e1071 type-3 moments,
# NA->Undefined / neighbourless->Isolated handling, and per-observation RNG
# seeding (results identical for any n.cores) are all as in local_g.
#
# Inputs: as local_g (x, listw, nsim, iseed, p.value, n.cores, moments).
# Output: numeric matrix of class c("localG", "matrix", "array"); columns G*i,
#   Z.G*i, "Pr(folded) Sim" (+ E.G*i, Var.G*i, Skew.G*i, Kurt.G*i when
#   moments = TRUE). Attributes: cluster (Not significant / High-High / Low-Low /
#   Undefined / Isolated), gstari = TRUE, call.
local_gstar <- function(x, listw,
                        nsim    = 999L,
                        iseed   = NULL,
                        p.value = 0.05,
                        n.cores = 1L,
                        moments = FALSE) 
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

  n.cores <- max(1L,as.integer(n.cores))
  p.value <- as.double(p.value)
  seed    <- if (is.null(iseed)) 123456789.0 else as.double(iseed)
  nsim    <- as.integer(nsim)
  if (nsim < 1L) {
    stop("nsim must be at least 1.")
  }

  ## ------------------------------------------------------------------ ##
  ## 2. Identify NAs and Isolates                                       ##
  ## ------------------------------------------------------------------ ##
  isolate_mask <- vapply(listw$neighbours, function(nbrs) {
    length(nbrs) == 0L || (length(nbrs) == 1L && nbrs[1] == 0L)
  }, logical(1L))
  na_mask <- is.na(x)
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
    "r_localgstar",
    csr$row_ptr,
    csr$col_idx,
    csr$weights,
    x_in,
    undef_flag,
    nsim,
    seed,
    n.cores,
    p.value
  )
  obs_gstar <- raw$gstar_val
  p_folded  <- raw$p_val
  p_folded[undef_mask] <- NA_real_

  EG_sim <- raw$mean
  EG_sim[undef_mask] <- NA_real_
  VG_sim <- raw$var
  VG_sim[undef_mask] <- NA_real_
  Skew.Gstari <- raw$skew
  Skew.Gstari[undef_mask] <- NA_real_
  Kurt.Gstari <- raw$kurt
  Kurt.Gstari[undef_mask] <- NA_real_

  # Standardised Z-score based on permutation moments
  Z.Gstari <- (obs_gstar - EG_sim) / sqrt(VG_sim)
  Z.Gstari[undef_mask] <- NA_real_

  ## ------------------------------------------------------------------ ##
  ## 4. Build Result Matrix and Cluster Factor                           ##
  ## ------------------------------------------------------------------ ##
  if (moments) {
    res <- matrix(
      c(obs_gstar, Z.Gstari, p_folded, EG_sim, VG_sim, Skew.Gstari, Kurt.Gstari),
      nrow = n,
      dimnames = list(attr(listw, "region.id"),
                      c("G*i", "Z.G*i", "Pr(folded) Sim",
                        "E.G*i", "Var.G*i", "Skew.G*i", "Kurt.G*i"))
    )
  } else {
    res <- matrix(
      c(obs_gstar, Z.Gstari, p_folded),
      nrow = n,
      dimnames = list(attr(listw, "region.id"),
                      c("G*i", "Z.G*i", "Pr(folded) Sim"))
    )
  }

  cluster <- factor(raw$cluster, levels = 0:4,
                    labels = c("Not significant", "High-High", "Low-Low",
                               "Undefined", "Isolated"))

  ## ------------------------------------------------------------------ ##
  ## 5. Return Value Assembly                                           ##
  ## ------------------------------------------------------------------ ##
  attr(res, "cluster")   <- cluster
  attr(res, "gstari")    <- TRUE
  attr(res, "call")      <- match.call()

  class(res) <- c("localG", "matrix", "array")
  res
}
