/* localbimoran.c
 *
 * Pure-C computation kernels for bivariate Local Moran's I.
 * No SEXP types are used here; all data is passed as plain C arrays.
 * The SEXP wrapper lives in R_export.c.
 *
 * Statistic, permutation-tail, and cluster conventions are adapted from
 * GeoDa/libgeoda's BiLocalMoran implementation (GPL-3-or-later).
 */
#include <stdlib.h>
#include <math.h>
#include "fastLISA.h"

/* ------------------------------------------------------------------
 * compute_bimoran
 *   Compute observed bivariate Moran's I for every observation.
 *
 * Inputs:
 *   n          - number of observations
 *   row_ptr    - CSR row pointers (length n+1, 0-based)
 *   col_idx    - CSR column indices (0-based)
 *   weights    - CSR weight values (row-standardised)
 *   z1, z2     - standardised variable vectors (length n)
 *   undef      - 1 if observation is NA/undefined, else 0 (length n)
 *   n_threads  - number of OpenMP threads
 *
 * Outputs (caller-allocated, length n each):
 *   bimoran_out - I_{bv,i} = z1_i * lag(z2)_i
 *   splag_out   - spatial lag of z2
 *   cluster_out - initial cluster code (before significance filtering)
 * ------------------------------------------------------------------ */
void compute_bimoran(int n, const int *row_ptr, const int *col_idx, const double *weights, const double *z1, const double *z2, const int *undef, double *bimoran_out, double *splag_out, int *cluster_out, int n_threads)
{
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(n_threads)
    #endif
    for (int i = 0; i < n; i++)
    {
        if (undef[i])
        {
            bimoran_out[i] = 0.0;
            splag_out[i]   = 0.0;
            cluster_out[i] = CLUSTER_UNDEFINED;
            continue;
        }

        int start = row_ptr[i];
        int end   = row_ptr[i + 1];

        if (start == end)
        {
            bimoran_out[i] = 0.0;
            splag_out[i]   = 0.0;
            cluster_out[i] = CLUSTER_ISOLATED;
            continue;
        }

        double lag   = 0.0;
        double w_sum = 0.0;
        for (int k = start; k < end; k++)
        {
            int j = col_idx[k];
            if (j != i && !undef[j])
            {
                lag   += weights[k] * z2[j];
                w_sum += weights[k];
            }
        }
        if (w_sum > 0.0)
            lag /= w_sum;

        splag_out[i]   = lag;
        bimoran_out[i] = z1[i] * lag;

        if      (z1[i] > 0.0 && lag > 0.0)
            cluster_out[i] = CLUSTER_HH;
        else if (z1[i] < 0.0 && lag < 0.0)
            cluster_out[i] = CLUSTER_LL;
        else if (z1[i] < 0.0 && lag > 0.0)
            cluster_out[i] = CLUSTER_LH;
        else if (z1[i] > 0.0 && lag < 0.0)
            cluster_out[i] = CLUSTER_HL;
        else
            cluster_out[i] = CLUSTER_NOT_SIG;
    }
}

/* ------------------------------------------------------------------
 * compute_bimoran_pvalues
 *   Conditional permutation test for bivariate Moran's I.
 *   For each valid observation i, the lag variable (z2) is shuffled
 *   among valid observations excluding i itself, and the fraction of
 *   permuted statistics at least as extreme as the observed value
 *   (folded, two-tailed) becomes the pseudo p-value.
 *
 * Parallelism & reproducibility (this design is shared by every *_pvalues
 * kernel in the package):
 *   - Observations are independent, so the i-loop runs in an OpenMP parallel
 *     region with schedule(dynamic, 8): per-observation cost varies with the
 *     neighbour count, and dynamic chunks keep the threads load-balanced.
 *   - Each thread owns a private perm_ws scratch workspace from plain malloc
 *     (perm_ws_alloc), NOT R_Calloc, because R's allocator is not thread-safe
 *     inside a parallel region. A barrier follows allocation so no thread starts
 *     the loop while another is still recording an allocation failure.
 *   - Errors are never raised inside the region: on allocation failure a thread
 *     sets alloc_failed and skips the loop body; error() is called once, after
 *     the region, with every per-thread buffer already freed.
 *   - Reproducibility: the RNG is re-seeded from the observation index i
 *     (base_seed + i*const) at the top of each iteration, so observation i draws
 *     the same permutation stream regardless of which thread runs it or in what
 *     order. Results are therefore bit-identical for any n_threads and any
 *     schedule (checked by tests/testthat/test-multicore-reproducible.R).
 *
 * Inputs:
 *   bimoran_obs - observed I_{bv,i} values (length n)
 *   permutations, base_seed, n_threads — as usual
 *
 * Output (caller-allocated, length n):
 *   pval_out    - pseudo p-values (NaN for undefined/isolated)
 * ------------------------------------------------------------------ */
void compute_bimoran_pvalues(int n, const int *row_ptr, const int *col_idx, const double *weights, const double *z1, const double *z2, const int *undef, const double *bimoran_obs, int permutations, uint64_t base_seed, int n_threads, int rank_pval, double *pval_out, double *mean_out, double *var_out, double *skew_out, double *kurt_out)
{
    int *pool      = R_Calloc(n, int);
    int  pool_size = 0;
    for (int i = 0; i < n; i++)
    {
        if (!undef[i])
            pool[pool_size++] = i;
    }

    /* e1071 type-3 moment corrections, constant across i: compute once.
       r^1.5 and r^2 via plain arithmetic (no pow); r = (m-1)/m. */
    double r         = (double)(permutations - 1) / permutations;
    double skew_corr = r * sqrt(r);
    double kurt_corr = r * r;

    int alloc_failed = 0;

    #ifdef _OPENMP
    #pragma omp parallel num_threads(n_threads)
    #endif
    {
        perm_ws ws;
        if (perm_ws_alloc(&ws, n, permutations) != 0)
        {
            #ifdef _OPENMP
            #pragma omp atomic write
            #endif
            alloc_failed = 1;
        }

        #ifdef _OPENMP
        #pragma omp barrier
        #pragma omp for schedule(dynamic, 8)
        #endif
        for (int i = 0; i < n; i++)
        {
            if (alloc_failed)
                continue;
            if (undef[i])
            {
                pval_out[i] = R_NaN;
                mean_out[i] = R_NaN;
                var_out[i]  = R_NaN;
                skew_out[i] = R_NaN;
                kurt_out[i] = R_NaN;
                continue;
            }

            /* Seed the RNG from the observation index for thread-independent,
               schedule-independent reproducibility. */
            uint64_t sm_state = base_seed + (uint64_t)i * 0x9e3779b97f4a7c15ULL;
            uint64_t rng_state[4];
            rng_state[0] = splitmix64(&sm_state);
            rng_state[1] = splitmix64(&sm_state);
            rng_state[2] = splitmix64(&sm_state);
            rng_state[3] = splitmix64(&sm_state);
            if (rng_state[0] == 0 && rng_state[1] == 0 && rng_state[2] == 0 && rng_state[3] == 0)
            {
                rng_state[0] = 1ULL;
            }

            /* Collect valid neighbours: count and weights in a single pass. */
            double w_sum = 0.0;
            int    nn    = 0;
            for (int k = row_ptr[i]; k < row_ptr[i + 1]; k++)
            {
                int j = col_idx[k];
                if (j != i && !undef[j])
                {
                    ws.w_valid[nn] = weights[k];
                    w_sum         += weights[k];
                    nn++;
                }
            }
            if (nn == 0)
            {
                pval_out[i] = R_NaN;
                mean_out[i] = R_NaN;
                var_out[i]  = R_NaN;
                skew_out[i] = R_NaN;
                kurt_out[i] = R_NaN;
                continue;
            }

            int local_size = 0;
            for (int p = 0; p < pool_size; p++)
            {
                if (pool[p] != i)
                    ws.draw[local_size++] = pool[p];
            }

            double   sum      = 0.0;
            uint64_t count_ge = 0;
            uint64_t count_eq = 0;

            // Step 4a: Run the permutation trials
            for (int perm = 0; perm < permutations; perm++)
            {
                // Fisher-Yates shuffle: select a random subset of neighbors from the pool
                // of valid non-focal observations.
                for (int k = 0; k < nn; k++)
                {
                    int idx = k + rng_int(rng_state, local_size - k);
                    int tmp = ws.draw[k];
                    ws.draw[k]   = ws.draw[idx];
                    ws.draw[idx] = tmp;
                }

                // Compute the spatial lag of variable z2 over the permuted neighbors
                double perm_lag = 0.0;
                for (int k = 0; k < nn; k++)
                {
                    perm_lag += ws.w_valid[k] * z2[ws.draw[k]];
                }
                if (w_sum > 0.0)
                    perm_lag /= w_sum;

                // Calculate the permuted Bivariate Moran's I statistic: I_i = z1_i * lag(z2)
                double val = z1[i] * perm_lag;
                ws.perm_vals[perm] = val;
                sum += val;
                if (val >= bimoran_obs[i])
                    count_ge++;
                if (val == bimoran_obs[i])
                    count_eq++;
            }

            // rgeoda folded pseudo-p: smaller tail, without doubling.
            uint64_t count_folded = count_ge;
            if ((uint64_t)permutations - count_ge < count_folded)
            {
                count_folded = (uint64_t)permutations - count_ge;
            }
            double p_count = ((double)count_folded + 1.0) / ((double)permutations + 1.0);

            // spdep rank-based folded pseudo-p: averaged rank of the observed value
            // among the nsim+1 values, mapped to the smaller tail (no doubling).
            uint64_t n_less = (uint64_t)permutations - count_ge;
            double xrank = (double)n_less + ((double)count_eq + 2.0) / 2.0;
            int    ri    = (int)xrank;
            double gr    = (double)ri / ((double)permutations + 1.0);
            double ls    = ((double)permutations + 2.0 - (double)ri) / ((double)permutations + 1.0);
            double p_rank = fmin(gr, ls);

            pval_out[i] = rank_pval ? p_rank : p_count;

            // Step 4c: Compute empirical moments of the permutation distribution
            // (Used by R wrapper to build expectation, variance, Z-score, and analytical p-values).
            double mean = sum / permutations;

            double sum2 = 0.0;
            double sum3 = 0.0;
            double sum4 = 0.0;
            for (int k = 0; k < permutations; k++)
            {
                double diff = ws.perm_vals[k] - mean;
                sum2 += diff * diff;
                sum3 += diff * diff * diff;
                sum4 += diff * diff * diff * diff;
            }

            double var = sum2 / (permutations - 1); // Sample variance
            double m2 = sum2 / permutations;        // 2nd central moment
            double m3 = sum3 / permutations;        // 3rd central moment
            double m4 = sum4 / permutations;        // 4th central moment

            double skew = 0.0;
            double kurt = 0.0;
            if (m2 > 0.0)
            {
                // Adjust skewness and kurtosis to match R's e1071 (type 3) standard:
                skew = (m3 / (m2 * sqrt(m2))) * skew_corr;
                kurt = (m4 / (m2 * m2)) * kurt_corr - 3.0;
            }

            mean_out[i] = mean;
            var_out[i]  = var;
            skew_out[i] = skew;
            kurt_out[i] = kurt;
        }

        perm_ws_free(&ws);
    }

    R_Free(pool);

    if (alloc_failed)
    {
        error("Memory allocation failed in OpenMP parallel region.");
    }
}
