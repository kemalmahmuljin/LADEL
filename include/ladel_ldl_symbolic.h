/**
 * @file ladel_ldl_symbolic.h
 * @author Ben Hermans
 * @brief The symbolic part of the factorization.
 */

#ifndef LADEL_LDL_SYMBOLIC_H
#define LADEL_LDL_SYMBOLIC_H

#include "ladel_types.h"


void quicksort(ladel_int *element_list, 
               ladel_int *ind_list,
               ladel_int low,
               ladel_int high);
               
/**
 * Symbolic part of the factorization.
 * 
 * @param M                 The matrix
 * @param sym               Symbolics of the factorization
 * @param ordering_method   Indicator to choose the ordering method: @a AMD, @a NO_ORDERING or @a GIVEN_ORDERING
 * @param Mpp               Output symmetric permutation of M if ordering is requested
 * @param work              LADEL workspace
 * @return                  Status
 */
ladel_int ladel_ldl_symbolic(   ladel_sparse_matrix *M, 
                                ladel_symbolics     *sym, 
                                ladel_int           ordering_method, 
                                ladel_int           num_fixed,
                                ladel_sparse_matrix *Mpp, 
                                ladel_work          *work);



#endif /*LADEL_LDL_SYMBOLIC_H*/
