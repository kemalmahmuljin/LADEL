/**
 * @file ladel_ldl_numeric.h
 * @author Ben Hermans
 * @brief The numerical part of the factorization.
 */

#ifndef LDL_NUMERIC_H
#define LDL_NUMERIC_H

#include "ladel_types.h"

/**
 * Numerical part of the factorization of @f$M + \alpha \begin{bmatrix}I_{n} & \\ & 0\end{bmatrix}@f$.
 * 
 * @param Mpp   The matrix M or its symmetric permutation
 * @param d     Diagonal parameters @f$\alpha@f$ and @f$n@f$
 * @param sym   Symbolics of the factorization
 * @param LD    Output factor struct
 * @param work  LADEL workspace
 * @return      Status
 */
ladel_int ladel_ldl_numeric_with_diag(ladel_sparse_matrix *Mpp, 
                                      ladel_diag          d, 
                                      ladel_symbolics     *sym, 
                                      ladel_factor        *LD, 
                                      ladel_work          *work);


/**
 * Numerical part of the factorization of @f$M@f$. Computing M = LDLT + E 
 * 
 * @param Mpp   The matrix M or its symmetric permutation
 * @param sym   Symbolics of the factorization
 * @param LD    Output factor struct
 * @param work  LADEL workspace
 * @param error_array Diagonal matrix D containing the errors in the factorization
 * @param beta  minimally allowed matrix dependent diagonal value
 * @return      Status
 * 
 */
ladel_int ladel_ldl_numeric_with_modification(ladel_sparse_matrix *Mpp, 
                                              ladel_diag          d, 
                                              ladel_symbolics     *sym, 
                                              ladel_factor        *LD, 
                                              ladel_work          *work,
                                              ladel_double        *error_array,
                                              ladel_double        beta);

/**
 * Numerical part of the factorization of @f$M@f$.
 * 
 * @param Mpp   The matrix M or its symmetric permutation
 * @param sym   Symbolics of the factorization
 * @param LD    Output factor struct
 * @param work  LADEL workspace
 * @param error_array Diagonal matrix D containing the errors in the factorization
 * @param beta  minimally allowed matrix dependent diagonal value
 * @return      Status
 * 
 */
ladel_int ladel_ldl_numeric(ladel_sparse_matrix *Mpp, 
                            ladel_symbolics     *sym, 
                            ladel_factor        *LD, 
                            ladel_work          *work,
                            ladel_double        *error_array,
                            ladel_double        beta);

#endif /*LDL_NUMERIC_H*/