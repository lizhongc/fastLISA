/* init.c
 *
 * R routine registration for fastLISA.
 *
 * Each r_* SEXP wrapper is declared here so the linker resolves the symbol,
 * then all wrappers are registered in CallMethods[].  The actual computation
 * lives in the per-statistic kernels (localbimoran.c, localgeary.c, ...); the
 * SEXP shims are all defined together in R_export.c.
 */
#include "fastLISA.h"        /* first: pulls omp.h before R's Rinternals.h remap macros */
#include <R_ext/Rdynload.h>

/* Forward-declare every SEXP entry point (defined in R_export.c) */
SEXP r_bi_localmoran(SEXP r_row_ptr, SEXP r_col_idx, SEXP r_weights, SEXP r_data1, SEXP r_data2, SEXP r_undef, SEXP r_permutations, SEXP r_seed, SEXP r_n_threads, SEXP r_sig_cutoff, SEXP r_rank_pval);

SEXP r_localgeary(SEXP r_row_ptr, SEXP r_col_idx, SEXP r_weights, SEXP r_data, SEXP r_undef, SEXP r_permutations, SEXP r_seed, SEXP r_n_threads, SEXP r_sig_cutoff, SEXP r_rank_pval);

SEXP r_localmultigeary(SEXP r_row_ptr, SEXP r_col_idx, SEXP r_weights, SEXP r_data_list, SEXP r_undef, SEXP r_permutations, SEXP r_seed, SEXP r_n_threads, SEXP r_sig_cutoff, SEXP r_rank_pval);

SEXP r_localg(SEXP r_row_ptr, SEXP r_col_idx, SEXP r_weights, SEXP r_data, SEXP r_undef, SEXP r_permutations, SEXP r_seed, SEXP r_n_threads, SEXP r_sig_cutoff, SEXP r_rank_pval);

SEXP r_localgstar(SEXP r_row_ptr, SEXP r_col_idx, SEXP r_weights, SEXP r_data, SEXP r_undef, SEXP r_permutations, SEXP r_seed, SEXP r_n_threads, SEXP r_sig_cutoff, SEXP r_rank_pval);

static const R_CallMethodDef CallMethods[] = {
    {"r_bi_localmoran",   (DL_FUNC)&r_bi_localmoran,   11},
    {"r_localgeary",      (DL_FUNC)&r_localgeary,      10},
    {"r_localmultigeary", (DL_FUNC)&r_localmultigeary, 10},
    {"r_localg",          (DL_FUNC)&r_localg,          10},
    {"r_localgstar",      (DL_FUNC)&r_localgstar,      10},
    {NULL, NULL, 0}
};

void R_init_fastLISA(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
