/**
 * @file ladel_matmat.h
 * @author Ben Hermans
 * @brief Routines to compute matrix matrix products. For now only @f$AA^T@f$ and @f$A\Sigma A^T@f$, 
 * with @f$\Sigma@f$ a diagonal matrix, are supported.
 */

#ifndef LADEL_MATMAT_H
#define LADEL_MATMAT_H

#include "ladel_types.h"

ladel_sparse_matrix *ladel_mat_mat_transpose(ladel_sparse_matrix *M, ladel_sparse_matrix *M_transpose, ladel_work *work);

ladel_sparse_matrix *ladel_mat_mat_transpose_pattern(ladel_sparse_matrix *M, ladel_sparse_matrix *M_transpose, ladel_work *work);

ladel_sparse_matrix *ladel_mat_diag_mat_transpose(ladel_sparse_matrix *M, ladel_sparse_matrix *M_transpose, ladel_double *diag, ladel_work *work);

ladel_sparse_matrix *ladel_mat_mat_transpose_advanced(ladel_sparse_matrix *M, ladel_sparse_matrix *M_transpose, ladel_double *diag, int values, ladel_work *work);

#endif /* LADEL_MATMAT_H */