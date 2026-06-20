/* workspace.c
 *
 * Per-thread permutation workspace shared by all permutation kernels.
 * The three scratch buffers (draw, w_valid, perm_vals) are carved out of a
 * single backing allocation so each thread allocates once and reuses the
 * workspace across every observation.
 *
 * Plain malloc/free is used (not R_Calloc): the workspace is allocated inside
 * OpenMP parallel regions, where R's allocator is not thread-safe.  No zeroing
 * is needed because every region is fully written before it is read on each
 * observation.
 */
#include <stdlib.h>
#include "fastLISA.h"

/* Allocate draw[n], w_valid[n], perm_vals[permutations] from one block.
 * Doubles are placed first so the int region stays naturally aligned.
 * Returns 0 on success, 1 on allocation failure. */
int perm_ws_alloc(perm_ws *ws, int n, int permutations)
{
    size_t bytes = (size_t)permutations * sizeof(double) + (size_t)n * sizeof(double) + (size_t)n * sizeof(int);

    char *block = (char *)malloc(bytes);
    if (!block) {
        ws->block     = NULL;
        ws->perm_vals = NULL;
        ws->w_valid   = NULL;
        ws->draw      = NULL;
        return 1;
    }

    ws->block     = block;
    ws->perm_vals = (double *)block;
    ws->w_valid   = (double *)(block + (size_t)permutations * sizeof(double));
    ws->draw      = (int    *)(block + (size_t)permutations * sizeof(double) + (size_t)n * sizeof(double));
    return 0;
}

void perm_ws_free(perm_ws *ws)
{
    free(ws->block);
    ws->block     = NULL;
    ws->perm_vals = NULL;
    ws->w_valid   = NULL;
    ws->draw      = NULL;
}
