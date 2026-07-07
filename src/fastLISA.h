#ifndef FASTLISA_H
#define FASTLISA_H

#include <stdint.h>

/* omp.h MUST be included before R's headers: Rinternals.h defines the remap
 * macro `match` (-> Rf_match), which otherwise corrupts omp.h's
 * `#pragma omp begin declare variant match(...)` clause under clang. */
#ifdef _OPENMP
  #include <omp.h>
#endif

#include <R.h>
#include <Rinternals.h>

/* Cluster codes */
#define CLUSTER_NOT_SIG      0
#define CLUSTER_HH           1
#define CLUSTER_LL           2
#define CLUSTER_LH           3
#define CLUSTER_HL           4
#define CLUSTER_UNDEFINED    5
#define CLUSTER_ISOLATED     6

/* G/G* cluster codes */
#define G_CLUSTER_NOT_SIG   0
#define G_CLUSTER_HH        1
#define G_CLUSTER_LL        2
#define G_CLUSTER_UNDEFINED 3
#define G_CLUSTER_ISOLATED  4

/* Multivariate Geary cluster codes */
#define MG_CLUSTER_NOT_SIG   0
#define MG_CLUSTER_POSITIVE  1
#define MG_CLUSTER_NEGATIVE  2
#define MG_CLUSTER_UNDEFINED 3
#define MG_CLUSTER_ISOLATED  4

/* -----------------------------------------------------------------------
 * RNG Functions (rng.c)
 * --------------------------------------------------------------------- */
uint64_t splitmix64(uint64_t *);
uint64_t xoshiro256plusplus(uint64_t [4]);
int      rng_int(uint64_t [4], int);

/* -----------------------------------------------------------------------
 * Per-thread permutation workspace (workspace.c)
 *   Bundles the three scratch buffers every permutation kernel needs so
 *   they are obtained in a single allocation per thread (reused across all
 *   observations) instead of buffer-by-buffer.  Plain malloc/free is used,
 *   not R_Calloc, because the workspace is allocated inside OpenMP regions
 *   where R_Calloc is not thread-safe.
 * --------------------------------------------------------------------- */
typedef struct
{
    int    *draw;       /* candidate neighbour indices (len n)            */
    double *w_valid;    /* neighbour weights for current obs (len n)      */
    double *perm_vals;  /* permuted statistics (len permutations)         */
    void   *block;      /* single backing allocation for the three above  */
} perm_ws;

/* Returns 0 on success, non-zero if allocation failed. */
int  perm_ws_alloc(perm_ws *, int, int);
void perm_ws_free(perm_ws *);

/* -----------------------------------------------------------------------
 * Pure computation functions — no SEXP types anywhere.
 * Each function takes plain C arrays and writes results into caller-
 * allocated output arrays.  The SEXP wrappers (R_export.c) unpack R objects,
 * call these, then pack the results back into R lists.
 * --------------------------------------------------------------------- */

/* localbimoran.c */
void compute_bimoran(int, const int *, const int *, const double *, const double *, const double *, const int *, double *, double *, int *, int);

void compute_bimoran_pvalues(int, const int *, const int *, const double *, const double *, const double *, const int *, const double *, int, uint64_t, int, int, double *, double *, double *, double *, double *);

/* localgeary.c */
void compute_localgeary(int, const int *, const int *, const double *, const double *, const int *, double *, int *, int);

void compute_localgeary_pvalues(int, const int *, const int *, const double *, const double *, const int *, const double *, int, uint64_t, int, int, double *, double *, double *, int *, double *, double *);

/* localmultigeary.c */
void compute_localmultigeary(int, int, const int *, const int *, const double *, double **, const int *, double *, int *, int);

void compute_localmultigeary_pvalues(int, int, const int *, const int *, const double *, double **, const int *, const double *, int, uint64_t, int, int, double *, double *, double *, int *, double *, double *);

/* localg.c */
void compute_localg(int, const int *, const int *, const double *, const double *, const int *, double *, int *, int);

void compute_localg_pvalues(int, const int *, const int *, const double *, const double *, const int *, const double *, int, uint64_t, int, int, double *, double *, double *, double *, double *);

/* localgstar.c */
void compute_localgstar(int, const int *, const int *, const double *, const double *, const int *, double *, int *, int);

void compute_localgstar_pvalues(int, const int *, const int *, const double *, const double *, const int *, const double *, int, uint64_t, int, int, double *, double *, double *, double *, double *);

#endif /* FASTLISA_H */
