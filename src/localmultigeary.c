/* localmultigeary.c
 *
 * Pure-C computation kernels for multivariate Local Geary's C.
 * No SEXP types are used here; all data is passed as plain C arrays.
 * The SEXP wrapper lives in R_export.c.
 *
 * Statistic, permutation-tail, and cluster conventions are adapted from
 * GeoDa/libgeoda's MultiGeary implementation (GPL-3-or-later).
 */
#include <stdlib.h>
#include <math.h>
#include "fastLISA.h"

/* ------------------------------------------------------------------
 * compute_localmultigeary
 *   Compute observed multivariate local Geary's C for every observation.
 *
 *   C_i = (1/K) * sum_{v=1}^{K} sum_j W*_ij (z^v_i - z^v_j)^2
 *
 * Inputs:
 *   num_vars - number of variables K
 *   z        - array of K pointers, each to a length-n standardised vector
 *   z_sq     - array of K pointers, each to element-wise squares (z^v)^2
 *   n_threads
 *
 * Outputs (caller-allocated, length n):
 *   geary_out   - observed C_i values
 * ------------------------------------------------------------------ */
void compute_localmultigeary(int n, int num_vars, const int *row_ptr, const int *col_idx, const double *weights, double **z, double **z_sq, const int *undef, double *geary_out, int *cluster_out, int n_threads)
{
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
    for (int i = 0; i < n; i++) {
        cluster_out[i] = MG_CLUSTER_NOT_SIG;

        if (undef[i]) {
            geary_out[i] = 0.0;
            cluster_out[i] = MG_CLUSTER_UNDEFINED;
            continue;
        }

        int start = row_ptr[i];
        int end   = row_ptr[i + 1];
        if (start == end) {
            geary_out[i] = 0.0;
            cluster_out[i] = MG_CLUSTER_ISOLATED;
            continue;
        }

        double w_sum = 0.0;
        for (int k = start; k < end; k++) {
            int j = col_idx[k];
            if (j != i && !undef[j])
                w_sum += weights[k];
        }

        if (w_sum == 0.0) {
            geary_out[i] = 0.0;
            cluster_out[i] = MG_CLUSTER_ISOLATED;
            continue;
        }

        double sum_c = 0.0;
        for (int v = 0; v < num_vars; v++) {
            double sp_lag    = 0.0;
            double sp_lag_sq = 0.0;
            for (int k = start; k < end; k++) {
                int j = col_idx[k];
                if (j != i && !undef[j]) {
                    sp_lag    += weights[k] * z[v][j];
                    sp_lag_sq += weights[k] * z_sq[v][j];
                }
            }
            sp_lag    /= w_sum;
            sp_lag_sq /= w_sum;
            sum_c += z_sq[v][i] - 2.0 * z[v][i] * sp_lag + sp_lag_sq;
        }
        geary_out[i] = sum_c / num_vars;
    }
}

/* ------------------------------------------------------------------
 * compute_localmultigeary_pvalues
 *   Conditional permutation test for multivariate Local Geary's C.
 *   One-tailed: compare observed value against permutation mean.
 *
 *   Parallelism & reproducibility: same design as compute_bimoran_pvalues in
 *   bimoran.c — OpenMP dynamic schedule, per-thread perm_ws via malloc, and a
 *   per-observation RNG seed that makes results identical for any n_threads.
 *
 * Outputs:
 *   pval_out      - pseudo p-values (NaN for undefined/isolated)
 * ------------------------------------------------------------------ */
void compute_localmultigeary_pvalues(int n, int num_vars, const int *row_ptr, const int *col_idx, const double *weights, double **z, double **z_sq, const int *undef, const double *geary_obs, int permutations, uint64_t base_seed, int n_threads, double *pval_out, double *mean_out, double *var_out, int *cluster_out, double *skew_out, double *kurt_out)
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

            double sum_perm = 0.0;

            // Step 4a: Run the permutation trials
            for (int perm = 0; perm < permutations; perm++) {
                // Fisher-Yates shuffle to draw a random configuration of neighbors
                for (int k = 0; k < nn; k++) {
                    int idx = k + rng_int(rng_state, local_size - k);
                    int tmp = ws.draw[k];
                    ws.draw[k]   = ws.draw[idx];
                    ws.draw[idx] = tmp;
                }

                // Compute the multivariate Geary statistic: average of univariate Geary statistics
                // over all num_vars variables.
                double sum_c = 0.0;
                for (int v = 0; v < num_vars; v++) {
                    double perm_lag    = 0.0;
                    double perm_lag_sq = 0.0;
                    for (int k = 0; k < nn; k++) {
                        int nb = ws.draw[k];
                        perm_lag    += ws.w_valid[k] * z[v][nb];
                        perm_lag_sq += ws.w_valid[k] * z_sq[v][nb];
                    }
                    if (w_sum > 0.0) {
                        perm_lag    /= w_sum;
                        perm_lag_sq /= w_sum;
                    }
                    sum_c += z_sq[v][i] - 2.0 * z[v][i] * perm_lag + perm_lag_sq;
                }
                double val = sum_c / num_vars;
                ws.perm_vals[perm] = val;
                sum_perm += val;
            }

            // Step 4b: empirical mean of the simulated stats (accumulated above)
            double mean_perm = sum_perm / permutations;

            uint64_t count_tail = 0;
            if (geary_obs[i] <= mean_perm) {
                for (int p = 0; p < permutations; p++) {
                    if (ws.perm_vals[p] <= geary_obs[i])
                        count_tail++;
                }
                cluster_out[i] = MG_CLUSTER_POSITIVE;
            } else {
                for (int p = 0; p < permutations; p++) {
                    if (ws.perm_vals[p] > geary_obs[i])
                        count_tail++;
                }
                cluster_out[i] = MG_CLUSTER_NEGATIVE;
            }
            pval_out[i] = ((double)count_tail + 1.0) / ((double)permutations + 1.0);

            /* Compute moments */
            double sum2 = 0.0;
            double sum3 = 0.0;
            double sum4 = 0.0;
            for (int k = 0; k < permutations; k++) {
                double diff = ws.perm_vals[k] - mean_perm;
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

            mean_out[i] = mean_perm;
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
