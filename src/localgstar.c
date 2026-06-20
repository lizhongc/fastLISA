/* localgstar.c
 *
 * Pure-C computation kernels for Local Getis-Ord G*.
 * No SEXP types are used here; all data is passed as plain C arrays.
 * The SEXP wrapper lives in R_export.c.
 *
 * Statistic, permutation-tail, and cluster conventions are adapted from
 * GeoDa/libgeoda's UniGstar implementation (GPL-3-or-later).
 */
#include <stdlib.h>
#include <math.h>
#include "fastLISA.h"

/* ------------------------------------------------------------------
 * compute_localgstar
 *   Compute observed Local G* statistic for every observation.
 *
 *   G*_i = (lag_i + x_i) / (w_sum + 1) / sum_x
 *   where lag_i = sum_{j != i} W*_ij x_j  (excluding self)
 *   and the self-contribution x_i is added with weight 1.
 *
 * Inputs:
 *   z        - raw data (NA→0) (length n)
 *   sum_x    - total sum of valid x values
 *   n_threads
 *
 * Outputs (caller-allocated, length n):
 *   gstar_out   - observed G*_i values
 * ------------------------------------------------------------------ */
void compute_localgstar(int n, const int *row_ptr, const int *col_idx, const double *weights, const double *z, const int *undef, double sum_x, double *gstar_out, int *cluster_out, int n_threads)
{
    if (sum_x == 0.0) {
        for (int i = 0; i < n; i++) {
            gstar_out[i] = 0.0;
            cluster_out[i] = row_ptr[i] == row_ptr[i + 1] ? G_CLUSTER_ISOLATED : G_CLUSTER_UNDEFINED;
        }
        return;
    }

#ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
    for (int i = 0; i < n; i++) {
        cluster_out[i] = G_CLUSTER_NOT_SIG;

        if (undef[i]) {
            gstar_out[i] = 0.0;
            cluster_out[i] = G_CLUSTER_UNDEFINED;
            continue;
        }

        int start = row_ptr[i];
        int end   = row_ptr[i + 1];
        if (start == end) {
            gstar_out[i] = 0.0;
            cluster_out[i] = G_CLUSTER_ISOLATED;
            continue;
        }

        double lag   = 0.0;
        double w_sum = 0.0;
        int    nn    = 0;
        for (int k = start; k < end; k++) {
            int j = col_idx[k];
            if (j != i && !undef[j]) {
                lag   += weights[k] * z[j];
                w_sum += weights[k];
                nn++;
            }
        }

        if (nn == 0 || w_sum == 0.0) {
            gstar_out[i] = 0.0;
            cluster_out[i] = G_CLUSTER_ISOLATED;
            continue;
        }

        /* Row-standardise self-loop and neighbors together to match rgeoda/spdep */
        double lag_star = ((lag / w_sum) * nn + z[i]) / (nn + 1.0);
        gstar_out[i] = lag_star / sum_x;
    }

    double sum_g = 0.0;
    int n_g = 0;
    for (int i = 0; i < n; i++) {
        if (!undef[i] && cluster_out[i] != G_CLUSTER_ISOLATED) {
            sum_g += gstar_out[i];
            n_g++;
        }
    }
    double mean_g = n_g > 0 ? sum_g / n_g : 0.0;
    for (int i = 0; i < n; i++) {
        if (!undef[i] && cluster_out[i] != G_CLUSTER_ISOLATED) {
            cluster_out[i] = gstar_out[i] >= mean_g ? G_CLUSTER_HH : G_CLUSTER_LL;
        }
    }
}

/* ------------------------------------------------------------------
 * compute_localgstar_pvalues
 *   Conditional permutation test for Local G*.
 *   Folded two-tailed: min(count_ge, permutations-count_ge).
 *   Self-contribution x_i is fixed (not permuted).
 *
 *   Parallelism & reproducibility: same design as compute_bimoran_pvalues in
 *   bimoran.c — OpenMP dynamic schedule, per-thread perm_ws via malloc, and a
 *   per-observation RNG seed that makes results identical for any n_threads.
 *
 * Output (caller-allocated, length n):
 *   pval_out - pseudo p-values (NaN for undefined/isolated)
 * ------------------------------------------------------------------ */
void compute_localgstar_pvalues(int n, const int *row_ptr, const int *col_idx, const double *weights, const double *z, const int *undef, const double *gstar_obs, double sum_x, int permutations, uint64_t base_seed, int n_threads, double *pval_out, double *mean_out, double *var_out, double *skew_out, double *kurt_out)
{
    int *pool      = R_Calloc(n, int);
    int  pool_size = 0;
    for (int i = 0; i < n; i++) {
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
        if (perm_ws_alloc(&ws, n, permutations) != 0) {
#ifdef _OPENMP
            #pragma omp atomic write
#endif
            alloc_failed = 1;
        }

#ifdef _OPENMP
        #pragma omp barrier
        #pragma omp for schedule(dynamic, 8)
#endif
        for (int i = 0; i < n; i++) {
            if (alloc_failed)
                continue;
            if (undef[i]) {
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
            if (rng_state[0] == 0 && rng_state[1] == 0 && rng_state[2] == 0 && rng_state[3] == 0) {
                rng_state[0] = 1ULL;
            }

            /* Collect valid neighbours: count and weights in a single pass. */
            double w_sum = 0.0;
            int    nn    = 0;
            for (int k = row_ptr[i]; k < row_ptr[i + 1]; k++) {
                int j = col_idx[k];
                if (j != i && !undef[j]) {
                    ws.w_valid[nn] = weights[k];
                    w_sum         += weights[k];
                    nn++;
                }
            }
            if (nn == 0) {
                pval_out[i] = R_NaN;
                mean_out[i] = R_NaN;
                var_out[i]  = R_NaN;
                skew_out[i] = R_NaN;
                kurt_out[i] = R_NaN;
                continue;
            }

            int local_size = 0;
            for (int p = 0; p < pool_size; p++) {
                if (pool[p] != i)
                    ws.draw[local_size++] = pool[p];
            }

            double   sum      = 0.0;
            uint64_t count_ge = 0;

            // Step 4a: Run the permutation trials
            for (int perm = 0; perm < permutations; perm++) {
                // Fisher-Yates shuffle to draw a random subset of neighbors from the pool
                for (int k = 0; k < nn; k++) {
                    int idx = k + rng_int(rng_state, local_size - k);
                    int tmp = ws.draw[k];
                    ws.draw[k]   = ws.draw[idx];
                    ws.draw[idx] = tmp;
                }

                // Compute the spatial lag of the variable over the permuted neighbor configuration
                double perm_lag = 0.0;
                for (int k = 0; k < nn; k++) {
                    perm_lag += ws.w_valid[k] * z[ws.draw[k]];
                }
                // Row-standardise self-loop and neighbors together to match rgeoda/spdep
                double perm_lag_star = ((perm_lag / w_sum) * nn + z[i]) / (nn + 1.0);
                double val = perm_lag_star / sum_x;
                ws.perm_vals[perm] = val;
                sum += val;
                if (val >= gstar_obs[i])
                    count_ge++;
            }

            uint64_t count_folded = count_ge;
            if ((uint64_t)permutations - count_ge < count_folded) {
                count_folded = (uint64_t)permutations - count_ge;
            }
            pval_out[i] = ((double)count_folded + 1.0) / ((double)permutations + 1.0);

            /* Compute moments */
            double mean = sum / permutations;

            double sum2 = 0.0;
            double sum3 = 0.0;
            double sum4 = 0.0;
            for (int k = 0; k < permutations; k++) {
                double diff = ws.perm_vals[k] - mean;
                sum2 += diff * diff;
                sum3 += diff * diff * diff;
                sum4 += diff * diff * diff * diff;
            }

            double var = sum2 / (permutations - 1);
            double m2 = sum2 / permutations;
            double m3 = sum3 / permutations;
            double m4 = sum4 / permutations;

            double skew = 0.0;
            double kurt = 0.0;
            if (m2 > 0.0) {
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

    if (alloc_failed) {
        error("Memory allocation failed in OpenMP parallel region.");
    }
}
