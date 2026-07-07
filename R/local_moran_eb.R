# local_moran_eb(): Local Moran's I on Empirical Bayes (EB) standardised rates
# (Assuncao & Reis 1999; Anselin 1995).
#
# For event/population rate data, raw rates from small populations are noisy. EB
# rate standardisation variance-stabilises them before computing local Moran's I.
# With p_i = event_i / base_i, global rate b = sum(event)/sum(base), and an EB
# variance component a_hat:
#
#   z_i = (p_i - b) / sqrt(a_hat + b / base_i)
#
# Local Moran's I is then computed on z. A SECOND (sample n-1 z-score)
# standardisation is applied internally because Moran's I requires it; the two
# standardisations are NOT redundant. This matches libgeoda/rgeoda EBLocalMoran.
# See .rate_standardize_eb below for the variance-component math.
#
# Inference/output mirror local_moran (delegates to local_moran_bv with x = y);
# clusters are read off the EB-rate Moran scatterplot. Undefined = NA event/base
# or base <= 0; Isolated = no neighbours.
#
# Inputs:
#   event    numeric vector of event counts, length n.
#   base     numeric vector of populations at risk, length n.
#   listw    spdep 'listw' spatial weights.
#   nsim, iseed, p.value, n.cores, moments  as in local_moran.
#   p.method pseudo p-value rule: "count" (default, standard) or "rank" ties-averaged.
#
# Output: numeric matrix of class c("local_moran_eb", "matrix"); columns Ii,
#   Z.Ii, "Pr(folded) Sim" (+ E.Ii/Var.Ii/Skew.Ii/Kurt.Ii when moments = TRUE).
#   Attributes: quadr, cluster (Not significant / High-High / Low-Low / Low-High /
#   High-Low / Undefined / Isolated), call, nsim.
local_moran_eb <- function(event, base, listw,
                           nsim    = 999L,
                           iseed   = NULL,
                           p.value = 0.05,
                           n.cores = 1L,
                           moments = FALSE,
                           p.method = c("count", "rank"))
{
  # Input sanity checks
  if (!is.numeric(event)) stop("'event' must be a numeric vector.")
  if (!is.numeric(base))  stop("'base' must be a numeric vector.")

  # Step 1: Perform Empirical Bayes (EB) rate standardisation.
  # This is a *variance-stabilising* transform z_i = (p_i - b) / se_i (see
  # .rate_standardize_eb below); it shrinks the noisy rates of small-population
  # areas toward the global rate but does NOT produce a mean-0/sd-1 variable.
  eb <- .rate_standardize_eb(event, base)

  # Step 2: Compute Local Moran's I on the EB rates.
  # The two standardisations here serve distinct purposes and are NOT redundant:
  #   - Step 1 (EB)            = variance stabilisation of the rates.
  #   - scale = TRUE (R side)  = the z-score standardisation that local Moran's I
  #                              itself requires, applied to the EB rates.
  # This two-stage sequence matches libgeoda/rgeoda's EBLocalMoran.
  # Note: Univariate Local Moran's I on a standardised variable x is mathematically
  # equivalent to Bivariate Local Moran's I with x = y and scale = TRUE, so we reuse
  # local_moran_bv to avoid duplicating the permutation engine.
  eb_input <- eb$results
  eb_input[eb$undef] <- NA_real_
  res <- local_moran_bv(x = eb_input, y = eb_input, listw = listw, nsim = nsim,
                       scale = TRUE, iseed = iseed,
                       p.value = p.value,
                       n.cores = n.cores,
                       moments = moments,
                       p.method = match.arg(p.method))

  # Step 3: Identify isolates and undefined units.
  # Isolate: node with 0 neighbors.
  # Undefined: NA event/base, or base <= 0.
  isolate_mask <- vapply(listw$neighbours, function(nbrs) {
    length(nbrs) == 0L || (length(nbrs) == 1L && nbrs[1] == 0L)
  }, logical(1L))
  undef_mask <- eb$undef | isolate_mask

  # Step 4: Scale results and lag for cluster classification.
  x_sc  <- .scale_sample(eb$results)
  ly_sc <- .listw_lag(x_sc, listw)

  # Step 5: Assign cluster codes matching the 0-6 convention.
  # HH = 1, LL = 2, LH = 3, HL = 4, Undefined = 5, Isolated = 6.
  cluster_int <- rep(0L, length(event))
  cluster_int[x_sc > 0 & ly_sc > 0] <- 1L # High-High
  cluster_int[x_sc < 0 & ly_sc < 0] <- 2L # Low-Low
  cluster_int[x_sc < 0 & ly_sc > 0] <- 3L # Low-High
  cluster_int[x_sc > 0 & ly_sc < 0] <- 4L # High-Low
  cluster_int[undef_mask]           <- 5L # Undefined
  cluster_int[isolate_mask]         <- 6L # Isolated

  # Step 6: Apply significance filtering.
  # Units with p-values above p.value are recoded to 0 (Not significant).
  p_col    <- res[, "Pr(folded) Sim"]
  sig_mask <- !is.na(p_col) & p_col <= p.value
  cluster_int[!sig_mask & !undef_mask] <- 0L

  # Step 7: Build output matrix and attach attributes.
  out <- unclass(res)
  colnames(out)[1:2]   <- c("Ii", "Z.Ii")
  if (moments) {
    colnames(out)[4:7] <- c("E.Ii", "Var.Ii", "Skew.Ii", "Kurt.Ii")
  }
  out[eb$undef, "Ii"]  <- 0.0
  out[eb$undef, "Pr(folded) Sim"] <- NA_real_

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

  attr(out, "quadr")   <- attr(res, "quadr")
  attr(out, "cluster") <- cluster
  attr(out, "call")    <- match.call()
  attr(out, "nsim")    <- nsim

  class(out) <- c("local_moran_eb", "matrix")
  out
}

# Internal helper for EB Rate standardisation, adapted from the GPL-licensed
# GeoDa/libgeoda EBLocalMoran formulation.
#
# Math behind EB Rate:
# Let event_i be case count, base_i be population. Raw rate p_i = event_i / base_i.
# Global rate b = sum(event) / sum(base).
# We estimate variance components:
#   gamma = sum_i base_i * (p_i - b)^2
#   a = (gamma / sum(base)) - (b / (sum(base) / N_valid))
#   a_hat = max(a, 0)
# Then standard error for area i:
#   se_i = sqrt(a_hat + b / base_i)
# Standardised rate is then:
#   z_i = (p_i - b) / se_i
.rate_standardize_eb <- function(event, base) {
  n <- length(event)
  undef <- is.na(event) | is.na(base) | (base <= 0)

  p <- numeric(n)
  p[!undef] <- event[!undef] / base[!undef]

  sP <- sum(base[!undef])
  sE <- sum(event[!undef])
  if (sP == 0) {
    return(list(results = numeric(n), undef = rep(TRUE, n)))
  }

  b_hat <- sE / sP

  obs_valid <- sum(!undef)
  gamma <- sum(base[!undef] * (p[!undef] - b_hat)^2)

  a <- (gamma / sP) - (b_hat / (sP / obs_valid))
  a_hat <- max(a, 0.0)

  results <- numeric(n)
  se <- sqrt(a_hat + b_hat / base[!undef])
  results[!undef] <- ifelse(se > 0, (p[!undef] - b_hat) / se, 0.0)

  list(results = results, undef = undef)
}
