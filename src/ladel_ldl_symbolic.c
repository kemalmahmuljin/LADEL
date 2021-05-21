#include "ladel_types.h"
#include "ladel_constants.h"
#include "ladel_global.h"
#include "ladel_permutation.h"
#include "ladel_etree.h"
#include "ladel_postorder.h"
#include "ladel_col_counts.h"
#include "ladel_debug_print.h"
#include "ladel_submatrix.h"

#ifdef DAMD
#include "amd.h"
#endif /*DAMD*/
void quick_sort(ladel_int* element_list, ladel_int* ind_list, ladel_int low, ladel_int high)
{
	ladel_int pivot, value1, value2, temp, temp_ind;
	if (low < high)
	{
		pivot = low;
		value1 = low;
		value2 = high;
		while (value1 < value2)
		{
			while (element_list[value1] <= element_list[pivot] && value1 <= high)
				value1++;
			while (element_list[value2] > element_list[pivot] && value2 >= low)
				value2--;
			if (value1 < value2)
			{
				temp = element_list[value1];
                temp_ind = ind_list[value1];
				
                element_list[value1] = element_list[value2];
				element_list[value2] = temp;

                ind_list[value1] = ind_list[value2];
                ind_list[value2] = temp_ind;
            }
		}
		temp = element_list[value2];
		temp_ind = ind_list[value2];

        element_list[value2] = element_list[pivot];
		element_list[pivot] = temp;
        ind_list[value2] = ind_list[pivot];
        ind_list[pivot] = temp_ind;

		quick_sort(element_list, ind_list, low, value2 - 1);
		quick_sort(element_list, ind_list, value2 + 1, high);
	}
}
ladel_int ladel_ldl_symbolic(ladel_sparse_matrix *M, ladel_symbolics *sym, ladel_int ordering_method, ladel_int num_fixed, ladel_sparse_matrix *Mpp, ladel_work* work)
{
    if (!M || !sym || !Mpp || !work) return FAIL;

    ladel_sparse_matrix *Mwork = M;
    if (ordering_method == AMD)
    {
        #ifdef DAMD
        ladel_int status;
        double Info [AMD_INFO];

        #ifdef DLONG
        status = amd_l_order(M->ncol, M->p, M->i, sym->p, NULL, Info);
        #else /*DLONG*/
        status = amd_order(M->ncol, M->p, M->i, sym->p, NULL, Info);
        #endif
        if (status != AMD_OK) return FAIL;
        
        #else /*DAMD*/
        sym->p = ladel_free(sym->p);
        #endif
    } else if (ordering_method > GIVEN_ORDERING)
    {   
        // This means AMD but, num_fixed last cols and rows unchanged
        ladel_int status, pattern = 1;
        double Info [AMD_INFO];

        ladel_sparse_matrix *M_sub = ladel_leading_principal_submatrix(M, num_fixed, pattern); 
        status = amd_l_order(M_sub->ncol, M_sub->p, M_sub->i, sym->p, NULL, Info);
        ladel_int nz[num_fixed], perm[num_fixed];
        
        for (int i=M->ncol-num_fixed; i<M->ncol; i++)
            perm[i-M->ncol+num_fixed] = i-M->ncol+num_fixed; 
        
        if (ordering_method == 3)
        {
            for (int i=M->ncol-num_fixed; i<M->ncol; i++)
            {
                if (M->nz) nz[i-M->ncol+num_fixed] = M->nz[i];
                else nz[i-M->ncol+num_fixed] = M->p[i+1] - M->p[i];
            }
            quick_sort(nz, perm, 0, num_fixed-1);
        }
        for (int i=0; i< num_fixed; i++)
        {
            sym->p[M_sub->ncol + i] = M_sub->ncol + perm[i];
        }
    } else if (ordering_method == NO_ORDERING)
    {
        sym->p = ladel_free(sym->p);
    }
    
    if (sym->p)
    {
        ladel_permute_symmetric_matrix(M, sym->p, Mpp, work);
        Mwork = Mpp;
        ladel_invert_permutation_vector(sym->p, sym->pinv, M->ncol);
    }
    #ifdef SIMPLE_COL_COUNTS
    ladel_etree_and_col_counts(Mwork, sym, work);
    #else
    ladel_etree(Mwork, sym, work);
    ladel_postorder(Mwork, sym, work);
    ladel_col_counts(Mwork, sym, work);
    //ladel_print_dense_int_vector_matlab(sym->col_counts, sym->ncol);
    #endif /* SIMPLE_COL_COUNTS */

    return SUCCESS;
}

