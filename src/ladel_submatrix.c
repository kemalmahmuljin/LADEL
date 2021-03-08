#include "ladel_types.h"
#include "ladel_global.h"

ladel_sparse_matrix *ladel_column_submatrix(ladel_sparse_matrix *M, ladel_int* cols, ladel_int nb_cols)
{
    if (!M) return NULL;
    if (!cols) return ladel_sparse_alloc_empty(M->nrow, M->ncol, M->symmetry, M->values, FALSE);

    ladel_int index, col, nnz = 0, index_M;
    for (index = 0; index < nb_cols; index++)
    {
        col = cols[index];
        nnz += (M->nz) ? M->nz[col] : M->p[col+1] - M->p[col]; 
    }

    if (nnz == 0) return ladel_sparse_alloc_empty(M->nrow, M->ncol, M->symmetry, M->values, FALSE);

    ladel_sparse_matrix *M_sub = ladel_sparse_alloc(M->nrow, nb_cols, nnz, M->symmetry, M->values, FALSE);
    nnz = 0;
    M_sub->p[0] = 0;
    for (index = 0; index < nb_cols; index++)
    {
        col = cols[index];
        LADEL_FOR(index_M, M, col)
        {
            M_sub->i[nnz] = M->i[index_M];
            M_sub->x[nnz] = M->x[index_M];
            nnz++; 
        }
        M_sub->p[index+1] = nnz;
    }
    
    return M_sub;
}
ladel_sparse_matrix *ladel_leading_principal_submatrix(ladel_sparse_matrix *M, ladel_int p, ladel_int pattern)
{
    if (!M) return NULL;
    ladel_int index;
    ladel_int ndim = M->ncol - p, nnzmax = M->p[ndim], nnz = nnzmax ;
    for (index = 0; index < nnzmax; index++)
    {
        if (M->i[index] > ndim - 1)
        {
            nnz--;            
        }
    }
    ladel_sparse_matrix *M_sub = ladel_sparse_alloc(ndim, ndim, nnz, M->symmetry, M->values, FALSE);

    M_sub->p[0] = 0;

    int it  = 0, diff = 0;
    for (index = 1; index < ndim+1; index++)
    {
        for (int j=M->p[index-1]; j<M->p[index]; j++){
            if (M->i[j] < ndim)
            {
                M_sub->i[it] = M->i[j];
                if (!pattern)
                {
                    M_sub->x[it] = M->x[j];
                }
                it++;
            }
            else{
            diff = M->p[index] - it; 
            break;
            }
        }
        M_sub->p[index] = M->p[index]-diff;
    }
    return M_sub;
}