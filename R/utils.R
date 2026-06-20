# utils.R – Helper utilities for fastLISA
#
# .listw2csr()  Convert a spdep listw object to 0-based CSR format
#               expected by the C backend.

.listw2csr <- function(listw) {
  if (!inherits(listw, "listw"))
    stop("'listw' must be a spdep listw object.")

  # Step 1: Compute cumulative counts of non-isolate neighbors.
  # spdep represents isolated units as a single 0L neighbor, which we treat as 0 active neighbors.
  counts <- vapply(listw$neighbours, function(nbrs) {
    if (length(nbrs) == 0L || (length(nbrs) == 1L && nbrs[1] == 0L)) 0L else length(nbrs)
  }, integer(1L))
  # row_ptr has length n + 1. The i-th row's neighbors start at col_idx[row_ptr[i] + 1]
  # and end at col_idx[row_ptr[i + 1]] (in R's 1-based indexing).
  row_ptr <- as.integer(c(0L, cumsum(counts)))

  if (sum(counts) == 0L) {
    # All observations are isolates (no spatial connections). Return empty index and weights vectors.
    return(list(
      row_ptr = row_ptr,
      col_idx = integer(0L),
      weights = double(0L)
    ))
  }

  # Step 2: Extract and clean neighbor indices.
  # We discard the 0L isolates and map R's 1-based indices to C's 0-based indices.
  clean_neighbours <- lapply(listw$neighbours, function(nbrs) {
    if (length(nbrs) == 1L && nbrs[1] == 0L) integer(0L) else nbrs
  })
  col_idx <- as.integer(unlist(clean_neighbours)) - 1L  # 1-based to 0-based

  # Step 3: Extract and row-standardise weights.
  # For each unit i, if the sum of weights is greater than 0, we divide each weight by the sum
  # so that the sum of weights for each row is exactly 1.0 (row-standardisation).
  clean_weights <- lapply(seq_along(listw$weights), function(i) {
    w <- listw$weights[[i]]
    if (counts[i] == 0L || length(w) == 0L) return(numeric(0L))
    w
  })

  norm_weights <- lapply(clean_weights, function(w) {
    if (length(w) == 0L) return(numeric(0L))
    s <- sum(w)
    if (s == 0) return(w)
    w / s
  })
  weights <- as.double(unlist(norm_weights))

  list(row_ptr = row_ptr, col_idx = col_idx, weights = weights)
}

# .prepare_lisa_data(x, undef_mask, scale): coerce x to double; if scale = TRUE,
# apply the sample (n-1) z-score over the valid (non-undefined) entries; then set
# undefined entries to 0 (the C kernels skip them via the undef flag). Returns the
# prepared numeric vector passed to the C backend.
.prepare_lisa_data <- function(x, undef_mask, scale = TRUE) {
  out <- as.double(x)
  valid <- !undef_mask

  if (scale) {
    m <- mean(out[valid])
    s <- sd(out[valid])
    if (is.na(s) || s == 0) s <- 1.0
    out[valid] <- (out[valid] - m) / s
  }

  out[undef_mask] <- 0.0
  out
}
