/* R_export.c
 *
 * SEXP wrappers for every fastLISA statistic, consolidated into one file.
 * Each function unpacks R objects -> calls the matching pure-C kernel
 * (compute_*.c) -> packs the results into an R named list.  The shims are
 * deliberately thin: no numerics live here.
 *
 *   r_bi_localmoran    -> localbimoran.c
 *   r_localgeary       -> localgeary.c
 *   r_localmultigeary  -> localmultigeary.c
 *   r_localg           -> localg.c
 *   r_localgstar       -> localgstar.c
 *
 * Registration lives in init.c.
 */
#include <stdlib.h>
#include "fastLISA.h"

/* Map any NaN emitted by the kernels (R_NaN) to R's NA so it prints/tests as
 * NA rather than NaN.  Necessary: the kernels emit R_NaN for undefined or
 * neighbourless observations, and some of those (e.g. nodes whose neighbours
 * are all NA) are not caught by the R-side undef mask. */
static void nan_to_na(double *x, int n)
{
    for (int i = 0; i < n; i++)
        if (ISNAN(x[i]))
            x[i] = NA_REAL;
}

/* ==================================================================
 * Bivariate Local Moran's I
 * ================================================================== */
SEXP r_bi_localmoran(SEXP r_row_ptr, SEXP r_col_idx, SEXP r_weights, SEXP r_data1, SEXP r_data2, SEXP r_undef, SEXP r_permutations, SEXP r_seed, SEXP r_n_threads, SEXP r_sig_cutoff, SEXP r_rank_pval)
{
    /* --- Unpack scalars ------------------------------------------------- */
    int      n            = length(r_data1);
    int      permutations = INTEGER(r_permutations)[0];
    uint64_t seed         = (uint64_t)REAL(r_seed)[0];
    int      n_threads    = INTEGER(r_n_threads)[0];
    double   cutoff       = REAL(r_sig_cutoff)[0];
    int      rank_pval    = (asLogical(r_rank_pval) == TRUE);
    if (seed == 0)
        seed = 123456789ULL;
    if (n_threads < 1)
        n_threads = 1;

    /* --- Unpack array pointers ------------------------------------------ */
    int    *row_ptr = INTEGER(r_row_ptr);
    int    *col_idx = INTEGER(r_col_idx);
    double *weights = REAL(r_weights);
    int    *undef   = INTEGER(r_undef);

    /* --- Prepared input arrays (scaled and NA-filled by R) --------------- */
    double *z1 = REAL(r_data1);
    double *z2 = REAL(r_data2);

    /* --- Allocate R result vectors (kernels write straight into them) --- */
    SEXP out       = PROTECT(allocVector(VECSXP,  8));
    SEXP names     = PROTECT(allocVector(STRSXP,  8));
    SEXP r_bimoran = PROTECT(allocVector(REALSXP, n));
    SEXP r_splag   = PROTECT(allocVector(REALSXP, n));
    SEXP r_pval    = PROTECT(allocVector(REALSXP, n));
    SEXP r_cluster = PROTECT(allocVector(INTSXP,  n));
    SEXP r_mean    = PROTECT(allocVector(REALSXP, n));
    SEXP r_var     = PROTECT(allocVector(REALSXP, n));
    SEXP r_skew    = PROTECT(allocVector(REALSXP, n));
    SEXP r_kurt    = PROTECT(allocVector(REALSXP, n));

    /* --- Dispatch to pure-C kernels -------------------------------------- */
    compute_bimoran(n, row_ptr, col_idx, weights, z1, z2, undef, REAL(r_bimoran), REAL(r_splag), INTEGER(r_cluster), n_threads);

    compute_bimoran_pvalues(n, row_ptr, col_idx, weights, z1, z2, undef, REAL(r_bimoran), permutations, seed, n_threads, rank_pval, REAL(r_pval), REAL(r_mean), REAL(r_var), REAL(r_skew), REAL(r_kurt));

    /* --- Apply significance cutoff -------------------------------------- */
    for (int i = 0; i < n; i++)
    {
        if (INTEGER(r_cluster)[i] >= CLUSTER_HH && INTEGER(r_cluster)[i] <= CLUSTER_HL)
        {
            if (ISNAN(REAL(r_pval)[i]) || REAL(r_pval)[i] > cutoff)
            {
                INTEGER(r_cluster)[i] = CLUSTER_NOT_SIG;
            }
        }
    }

    /* --- NaN → NA ------------------------------------------------------- */
    nan_to_na(REAL(r_pval), n);
    nan_to_na(REAL(r_mean), n);
    nan_to_na(REAL(r_var), n);
    nan_to_na(REAL(r_skew), n);
    nan_to_na(REAL(r_kurt), n);

    /* --- Assemble named list -------------------------------------------- */
    SET_VECTOR_ELT(out, 0, r_bimoran);
    SET_VECTOR_ELT(out, 1, r_splag);
    SET_VECTOR_ELT(out, 2, r_pval);
    SET_VECTOR_ELT(out, 3, r_cluster);
    SET_VECTOR_ELT(out, 4, r_mean);
    SET_VECTOR_ELT(out, 5, r_var);
    SET_VECTOR_ELT(out, 6, r_skew);
    SET_VECTOR_ELT(out, 7, r_kurt);

    SET_STRING_ELT(names, 0, mkChar("bimoran"));
    SET_STRING_ELT(names, 1, mkChar("sp_lag"));
    SET_STRING_ELT(names, 2, mkChar("p_val"));
    SET_STRING_ELT(names, 3, mkChar("cluster"));
    SET_STRING_ELT(names, 4, mkChar("mean"));
    SET_STRING_ELT(names, 5, mkChar("var"));
    SET_STRING_ELT(names, 6, mkChar("skew"));
    SET_STRING_ELT(names, 7, mkChar("kurt"));
    setAttrib(out, R_NamesSymbol, names);

    UNPROTECT(10);
    return out;
}

/* ==================================================================
 * Univariate Local Geary's C
 * ================================================================== */
SEXP r_localgeary(SEXP r_row_ptr, SEXP r_col_idx, SEXP r_weights, SEXP r_data, SEXP r_undef, SEXP r_permutations, SEXP r_seed, SEXP r_n_threads, SEXP r_sig_cutoff, SEXP r_rank_pval)
{
    /* --- Unpack scalars ------------------------------------------------- */
    int      n            = length(r_data);
    int      permutations = INTEGER(r_permutations)[0];
    uint64_t seed         = (uint64_t)REAL(r_seed)[0];
    int      n_threads    = INTEGER(r_n_threads)[0];
    double   cutoff       = REAL(r_sig_cutoff)[0];
    int      rank_pval    = (asLogical(r_rank_pval) == TRUE);
    if (seed == 0)
        seed = 123456789ULL;
    if (n_threads < 1)
        n_threads = 1;

    /* --- Unpack array pointers ------------------------------------------ */
    int    *row_ptr = INTEGER(r_row_ptr);
    int    *col_idx = INTEGER(r_col_idx);
    double *weights = REAL(r_weights);
    int    *undef   = INTEGER(r_undef);

    /* --- Data prepared by R (alias) ------------------------------------- */
    const double *z = REAL(r_data);

    /* --- Allocate R result vectors (kernels write straight into them) --- */
    SEXP out       = PROTECT(allocVector(VECSXP,  7));
    SEXP names     = PROTECT(allocVector(STRSXP,  7));
    SEXP r_geary   = PROTECT(allocVector(REALSXP, n));
    SEXP r_pval    = PROTECT(allocVector(REALSXP, n));
    SEXP r_mean    = PROTECT(allocVector(REALSXP, n));
    SEXP r_var     = PROTECT(allocVector(REALSXP, n));
    SEXP r_skew    = PROTECT(allocVector(REALSXP, n));
    SEXP r_kurt    = PROTECT(allocVector(REALSXP, n));
    SEXP r_cluster = PROTECT(allocVector(INTSXP,  n));

    /* --- Dispatch to pure-C kernels -------------------------------------- */
    compute_localgeary(n, row_ptr, col_idx, weights, z, undef, REAL(r_geary), INTEGER(r_cluster), n_threads);

    compute_localgeary_pvalues(n, row_ptr, col_idx, weights, z, undef, REAL(r_geary), permutations, seed, n_threads, rank_pval, REAL(r_pval), REAL(r_mean), REAL(r_var), INTEGER(r_cluster), REAL(r_skew), REAL(r_kurt));

    /* --- Apply significance cutoff -------------------------------------- */
    for (int i = 0; i < n; i++)
    {
        if (INTEGER(r_cluster)[i] >= CLUSTER_HH && INTEGER(r_cluster)[i] <= CLUSTER_HL)
        {
            if (ISNAN(REAL(r_pval)[i]) || REAL(r_pval)[i] > cutoff)
                INTEGER(r_cluster)[i] = CLUSTER_NOT_SIG;
        }
    }

    /* --- NaN → NA ------------------------------------------------------- */
    nan_to_na(REAL(r_pval), n);
    nan_to_na(REAL(r_mean), n);
    nan_to_na(REAL(r_var), n);
    nan_to_na(REAL(r_skew), n);
    nan_to_na(REAL(r_kurt), n);

    /* --- Assemble named list -------------------------------------------- */
    SET_VECTOR_ELT(out, 0, r_geary);
    SET_VECTOR_ELT(out, 1, r_pval);
    SET_VECTOR_ELT(out, 2, r_mean);
    SET_VECTOR_ELT(out, 3, r_var);
    SET_VECTOR_ELT(out, 4, r_skew);
    SET_VECTOR_ELT(out, 5, r_kurt);
    SET_VECTOR_ELT(out, 6, r_cluster);

    SET_STRING_ELT(names, 0, mkChar("geary"));
    SET_STRING_ELT(names, 1, mkChar("p_val"));
    SET_STRING_ELT(names, 2, mkChar("mean"));
    SET_STRING_ELT(names, 3, mkChar("var"));
    SET_STRING_ELT(names, 4, mkChar("skew"));
    SET_STRING_ELT(names, 5, mkChar("kurt"));
    SET_STRING_ELT(names, 6, mkChar("cluster"));
    setAttrib(out, R_NamesSymbol, names);

    UNPROTECT(9);
    return out;
}

/* ==================================================================
 * Multivariate Local Geary's C
 * ================================================================== */
SEXP r_localmultigeary(SEXP r_row_ptr, SEXP r_col_idx, SEXP r_weights, SEXP r_data_list, SEXP r_undef, SEXP r_permutations, SEXP r_seed, SEXP r_n_threads, SEXP r_sig_cutoff, SEXP r_rank_pval)
{
    /* --- Unpack scalars ------------------------------------------------- */
    int      n            = length(r_undef);
    int      num_vars     = length(r_data_list);
    int      permutations = INTEGER(r_permutations)[0];
    uint64_t seed         = (uint64_t)REAL(r_seed)[0];
    int      n_threads    = INTEGER(r_n_threads)[0];
    double   cutoff       = REAL(r_sig_cutoff)[0];
    int      rank_pval    = (asLogical(r_rank_pval) == TRUE);
    if (seed == 0)
        seed = 123456789ULL;
    if (n_threads < 1)
        n_threads = 1;

    /* --- Unpack array pointers ------------------------------------------ */
    int    *row_ptr = INTEGER(r_row_ptr);
    int    *col_idx = INTEGER(r_col_idx);
    double *weights = REAL(r_weights);
    int    *undef   = INTEGER(r_undef);

    /* --- Variables prepared by R: gather the per-variable pointers ------ */
    double **z = R_Calloc(num_vars, double*);
    for (int v = 0; v < num_vars; v++)
    {
        z[v] = REAL(VECTOR_ELT(r_data_list, v));
    }

    /* --- Allocate R result vectors (kernels write straight into them) --- */
    SEXP out       = PROTECT(allocVector(VECSXP,  7));
    SEXP names     = PROTECT(allocVector(STRSXP,  7));
    SEXP r_geary   = PROTECT(allocVector(REALSXP, n));
    SEXP r_pval    = PROTECT(allocVector(REALSXP, n));
    SEXP r_mean    = PROTECT(allocVector(REALSXP, n));
    SEXP r_var     = PROTECT(allocVector(REALSXP, n));
    SEXP r_skew    = PROTECT(allocVector(REALSXP, n));
    SEXP r_kurt    = PROTECT(allocVector(REALSXP, n));
    SEXP r_cluster = PROTECT(allocVector(INTSXP,  n));

    /* --- Dispatch to pure-C kernels -------------------------------------- */
    compute_localmultigeary(n, num_vars, row_ptr, col_idx, weights, z, undef, REAL(r_geary), INTEGER(r_cluster), n_threads);

    compute_localmultigeary_pvalues(n, num_vars, row_ptr, col_idx, weights, z, undef, REAL(r_geary), permutations, seed, n_threads, rank_pval, REAL(r_pval), REAL(r_mean), REAL(r_var), INTEGER(r_cluster), REAL(r_skew), REAL(r_kurt));

    R_Free(z);

    /* --- Apply significance cutoff -------------------------------------- */
    for (int i = 0; i < n; i++)
    {
        if (INTEGER(r_cluster)[i] == MG_CLUSTER_POSITIVE || INTEGER(r_cluster)[i] == MG_CLUSTER_NEGATIVE)
        {
            if (ISNAN(REAL(r_pval)[i]) || REAL(r_pval)[i] > cutoff)
                INTEGER(r_cluster)[i] = MG_CLUSTER_NOT_SIG;
        }
    }

    /* --- NaN → NA ------------------------------------------------------- */
    nan_to_na(REAL(r_pval), n);
    nan_to_na(REAL(r_mean), n);
    nan_to_na(REAL(r_var), n);
    nan_to_na(REAL(r_skew), n);
    nan_to_na(REAL(r_kurt), n);

    /* --- Assemble named list -------------------------------------------- */
    SET_VECTOR_ELT(out, 0, r_geary);
    SET_VECTOR_ELT(out, 1, r_pval);
    SET_VECTOR_ELT(out, 2, r_mean);
    SET_VECTOR_ELT(out, 3, r_var);
    SET_VECTOR_ELT(out, 4, r_skew);
    SET_VECTOR_ELT(out, 5, r_kurt);
    SET_VECTOR_ELT(out, 6, r_cluster);

    SET_STRING_ELT(names, 0, mkChar("geary"));
    SET_STRING_ELT(names, 1, mkChar("p_val"));
    SET_STRING_ELT(names, 2, mkChar("mean"));
    SET_STRING_ELT(names, 3, mkChar("var"));
    SET_STRING_ELT(names, 4, mkChar("skew"));
    SET_STRING_ELT(names, 5, mkChar("kurt"));
    SET_STRING_ELT(names, 6, mkChar("cluster"));
    setAttrib(out, R_NamesSymbol, names);

    UNPROTECT(9);
    return out;
}

/* ==================================================================
 * Local Getis-Ord G
 * ================================================================== */
SEXP r_localg(SEXP r_row_ptr, SEXP r_col_idx, SEXP r_weights, SEXP r_data, SEXP r_undef, SEXP r_permutations, SEXP r_seed, SEXP r_n_threads, SEXP r_sig_cutoff, SEXP r_rank_pval)
{
    /* --- Unpack scalars ------------------------------------------------- */
    int      n            = length(r_data);
    int      permutations = INTEGER(r_permutations)[0];
    uint64_t seed         = (uint64_t)REAL(r_seed)[0];
    int      n_threads    = INTEGER(r_n_threads)[0];
    double   cutoff       = REAL(r_sig_cutoff)[0];
    int      rank_pval    = (asLogical(r_rank_pval) == TRUE);
    if (seed == 0)
        seed = 123456789ULL;
    if (n_threads < 1)
        n_threads = 1;

    /* --- Unpack array pointers ------------------------------------------ */
    int    *row_ptr = INTEGER(r_row_ptr);
    int    *col_idx = INTEGER(r_col_idx);
    double *weights = REAL(r_weights);
    int    *undef   = INTEGER(r_undef);

    /* --- Data prepared by R (NA already set to 0); alias --------------- */
    const double *z = REAL(r_data);

    /* --- Allocate R result vectors (kernels write straight into them) --- */
    SEXP out       = PROTECT(allocVector(VECSXP,  7));
    SEXP names     = PROTECT(allocVector(STRSXP,  7));
    SEXP r_gval    = PROTECT(allocVector(REALSXP, n));
    SEXP r_pval    = PROTECT(allocVector(REALSXP, n));
    SEXP r_mean    = PROTECT(allocVector(REALSXP, n));
    SEXP r_var     = PROTECT(allocVector(REALSXP, n));
    SEXP r_skew    = PROTECT(allocVector(REALSXP, n));
    SEXP r_kurt    = PROTECT(allocVector(REALSXP, n));
    SEXP r_cluster = PROTECT(allocVector(INTSXP,  n));

    /* --- Dispatch to pure-C kernels -------------------------------------- */
    compute_localg(n, row_ptr, col_idx, weights, z, undef, REAL(r_gval), INTEGER(r_cluster), n_threads);

    compute_localg_pvalues(n, row_ptr, col_idx, weights, z, undef, REAL(r_gval), permutations, seed, n_threads, rank_pval, REAL(r_pval), REAL(r_mean), REAL(r_var), REAL(r_skew), REAL(r_kurt));

    /* --- Apply significance cutoff -------------------------------------- */
    for (int i = 0; i < n; i++)
    {
        if (INTEGER(r_cluster)[i] == G_CLUSTER_HH || INTEGER(r_cluster)[i] == G_CLUSTER_LL)
        {
            if (ISNAN(REAL(r_pval)[i]) || REAL(r_pval)[i] > cutoff)
                INTEGER(r_cluster)[i] = G_CLUSTER_NOT_SIG;
        }
    }

    /* --- NaN → NA ------------------------------------------------------- */
    nan_to_na(REAL(r_pval), n);
    nan_to_na(REAL(r_mean), n);
    nan_to_na(REAL(r_var), n);
    nan_to_na(REAL(r_skew), n);
    nan_to_na(REAL(r_kurt), n);

    /* --- Assemble named list -------------------------------------------- */
    SET_VECTOR_ELT(out, 0, r_gval);
    SET_VECTOR_ELT(out, 1, r_pval);
    SET_VECTOR_ELT(out, 2, r_mean);
    SET_VECTOR_ELT(out, 3, r_var);
    SET_VECTOR_ELT(out, 4, r_skew);
    SET_VECTOR_ELT(out, 5, r_kurt);
    SET_VECTOR_ELT(out, 6, r_cluster);

    SET_STRING_ELT(names, 0, mkChar("g_val"));
    SET_STRING_ELT(names, 1, mkChar("p_val"));
    SET_STRING_ELT(names, 2, mkChar("mean"));
    SET_STRING_ELT(names, 3, mkChar("var"));
    SET_STRING_ELT(names, 4, mkChar("skew"));
    SET_STRING_ELT(names, 5, mkChar("kurt"));
    SET_STRING_ELT(names, 6, mkChar("cluster"));
    setAttrib(out, R_NamesSymbol, names);

    UNPROTECT(9);
    return out;
}

/* ==================================================================
 * Local Getis-Ord G*
 * ================================================================== */
SEXP r_localgstar(SEXP r_row_ptr, SEXP r_col_idx, SEXP r_weights, SEXP r_data, SEXP r_undef, SEXP r_permutations, SEXP r_seed, SEXP r_n_threads, SEXP r_sig_cutoff, SEXP r_rank_pval)
{
    /* --- Unpack scalars ------------------------------------------------- */
    int      n            = length(r_data);
    int      permutations = INTEGER(r_permutations)[0];
    uint64_t seed         = (uint64_t)REAL(r_seed)[0];
    int      n_threads    = INTEGER(r_n_threads)[0];
    double   cutoff       = REAL(r_sig_cutoff)[0];
    int      rank_pval    = (asLogical(r_rank_pval) == TRUE);
    if (seed == 0)
        seed = 123456789ULL;
    if (n_threads < 1)
        n_threads = 1;

    /* --- Unpack array pointers ------------------------------------------ */
    int    *row_ptr = INTEGER(r_row_ptr);
    int    *col_idx = INTEGER(r_col_idx);
    double *weights = REAL(r_weights);
    int    *undef   = INTEGER(r_undef);

    /* --- Data prepared by R (NA already set to 0); alias --------------- */
    const double *z = REAL(r_data);

    /* --- Allocate R result vectors (kernels write straight into them) --- */
    SEXP out       = PROTECT(allocVector(VECSXP,  7));
    SEXP names     = PROTECT(allocVector(STRSXP,  7));
    SEXP r_gstar   = PROTECT(allocVector(REALSXP, n));
    SEXP r_pval    = PROTECT(allocVector(REALSXP, n));
    SEXP r_mean    = PROTECT(allocVector(REALSXP, n));
    SEXP r_var     = PROTECT(allocVector(REALSXP, n));
    SEXP r_skew    = PROTECT(allocVector(REALSXP, n));
    SEXP r_kurt    = PROTECT(allocVector(REALSXP, n));
    SEXP r_cluster = PROTECT(allocVector(INTSXP,  n));

    /* --- Dispatch to pure-C kernels (kernels handle sum_x == 0) --------- */
    compute_localgstar(n, row_ptr, col_idx, weights, z, undef, REAL(r_gstar), INTEGER(r_cluster), n_threads);

    compute_localgstar_pvalues(n, row_ptr, col_idx, weights, z, undef, REAL(r_gstar), permutations, seed, n_threads, rank_pval, REAL(r_pval), REAL(r_mean), REAL(r_var), REAL(r_skew), REAL(r_kurt));

    /* --- Apply significance cutoff -------------------------------------- */
    for (int i = 0; i < n; i++)
    {
        if (INTEGER(r_cluster)[i] == G_CLUSTER_HH || INTEGER(r_cluster)[i] == G_CLUSTER_LL)
        {
            if (ISNAN(REAL(r_pval)[i]) || REAL(r_pval)[i] > cutoff)
                INTEGER(r_cluster)[i] = G_CLUSTER_NOT_SIG;
        }
    }

    /* --- NaN → NA ------------------------------------------------------- */
    nan_to_na(REAL(r_pval), n);
    nan_to_na(REAL(r_mean), n);
    nan_to_na(REAL(r_var), n);
    nan_to_na(REAL(r_skew), n);
    nan_to_na(REAL(r_kurt), n);

    /* --- Assemble named list -------------------------------------------- */
    SET_VECTOR_ELT(out, 0, r_gstar);
    SET_VECTOR_ELT(out, 1, r_pval);
    SET_VECTOR_ELT(out, 2, r_mean);
    SET_VECTOR_ELT(out, 3, r_var);
    SET_VECTOR_ELT(out, 4, r_skew);
    SET_VECTOR_ELT(out, 5, r_kurt);
    SET_VECTOR_ELT(out, 6, r_cluster);

    SET_STRING_ELT(names, 0, mkChar("gstar_val"));
    SET_STRING_ELT(names, 1, mkChar("p_val"));
    SET_STRING_ELT(names, 2, mkChar("mean"));
    SET_STRING_ELT(names, 3, mkChar("var"));
    SET_STRING_ELT(names, 4, mkChar("skew"));
    SET_STRING_ELT(names, 5, mkChar("kurt"));
    SET_STRING_ELT(names, 6, mkChar("cluster"));
    setAttrib(out, R_NamesSymbol, names);

    UNPROTECT(9);
    return out;
}
