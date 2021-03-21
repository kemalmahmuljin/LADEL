#include "ladel_types.h"
#include "ladel_global.h"
#include "ladel_constants.h"
#include "ladel_pattern.h"
#include "ladel_transpose.h"
#include "ladel_debug_print.h"

ladel_int ladel_ldl_numeric_with_diag(ladel_sparse_matrix *Mpp, ladel_diag d, ladel_symbolics *sym, ladel_factor *LD, ladel_work* work)
{
    if (!Mpp || !sym || !LD || !work) return FAIL;
    
    ladel_int row, col, index, ncol = Mpp->ncol, start, index_in_pattern;
    ladel_int *pattern = sym->pattern;
    ladel_double diag_elem, temp, L_elem;
    
    ladel_sparse_matrix *L = LD->L;
    ladel_double *D = LD->D, *Dinv = LD->Dinv;
    ladel_double *rhs = work->array_double_all_zeros_ncol1;
    ladel_int *col_pointers = work->array_int_ncol1;

    L->p[0] = col_pointers[0] = 0;
    for (index = 1; index < ncol; index++) 
        L->p[index] = col_pointers[index] = sym->col_counts[index-1];

    L->p[ncol] = sym->col_counts[ncol-1];

    for (col = 0; col < ncol; col++)
    {
        LADEL_FOR(index, Mpp, col)
            rhs[Mpp->i[index]] = Mpp->x[index];
        diag_elem = rhs[col];
        if ((LD->p && LD->p[col] < d.diag_size) || (!LD->p && col < d.diag_size)) diag_elem += d.diag_elem;
        rhs[col] = 0;

        start = ladel_nonzero_pattern_of_row_in_L(Mpp, sym, col);
        for (index_in_pattern = start; index_in_pattern < ncol; index_in_pattern++)
        {
            row = pattern[index_in_pattern];
            temp = rhs[row];
            L_elem = Dinv[row]*temp; /*L(col,row) = rhs(row) / L(row,row) / D(row,row)*/
            diag_elem -= L_elem*temp; /*D(col, col) -= L(col,row)*D(row,row)*L(col,row)*/
            rhs[row] = 0;
            
            /* Gaussian elimination */  
            for (index = L->p[row]; index < col_pointers[row]; index++)
                rhs[L->i[index]] -= L->x[index]*temp;
            
            index = col_pointers[row];
            col_pointers[row]++;
            L->i[index] = col;
            L->x[index] = L_elem;
        }
        /*Return FAIL if eigenvalue (close to) zero*/
        if (LADEL_ABS(diag_elem) < 1e-15) 
        {
            ladel_print("LADEL ERROR: MATRIX (POSSIBLY) NOT FULL RANK (diagonal element of %le)\n", diag_elem);
            return FAIL; 
        }

        D[col] = diag_elem;
        Dinv[col] = 1/diag_elem;
    }

    for (index = 0; index < ncol; index++) L->nz[index] = col_pointers[index] - L->p[index];
    return SUCCESS;
}



ladel_int ladel_ldl_numeric_with_modification(ladel_sparse_matrix *Mpp, ladel_diag d, ladel_symbolics *sym, ladel_factor *LD, ladel_work* work, ladel_double* error_array, ladel_double beta)
{
    if (!Mpp || !sym || !LD || !work) return FAIL;
    
    ladel_int row, col, index, ncol = Mpp->ncol, start, index_in_pattern;
    ladel_int *pattern = sym->pattern;
    ladel_double diag_elem, temp, L_elem;
    ladel_int *added = ladel_malloc(ncol, sizeof(ladel_int));;
    ladel_double theta;

    ladel_sparse_matrix *L = LD->L;
    ladel_double *D = LD->D, *Dinv = LD->Dinv;
    ladel_double *rhs = work->array_double_all_zeros_ncol1;
    ladel_int *col_pointers = work->array_int_ncol1;
    
    ladel_sparse_matrix *C = ladel_sparse_allocate_and_copy(L);

    ladel_sparse_matrix *Mppt = ladel_transpose(Mpp, TRUE, work);
    L->p[0] = C->p[0] = col_pointers[0] = 0;
    for (index = 1; index < ncol; index++) 
    {
        L->p[index] = col_pointers[index] = sym->col_counts[index-1];
    }
    L->p[ncol] = sym->col_counts[ncol-1];
     
    for (col = 0; col < ncol; col++)
    {
        LADEL_FOR(index, Mpp, col)
            rhs[Mpp->i[index]] = Mpp->x[index];
        diag_elem = rhs[col];
        if ((LD->p && LD->p[col] < d.diag_size) || (!LD->p && col < d.diag_size)) diag_elem += d.diag_elem;
        rhs[col] = 0;

        start = ladel_nonzero_pattern_of_row_in_L(Mpp, sym, col);
        for (index_in_pattern = start; index_in_pattern < ncol; index_in_pattern++)
        {
            row = pattern[index_in_pattern];
            temp = rhs[row];
            L_elem = Dinv[row]*temp; /*L(col,row) = rhs(row) / L(row,row) / D(row,row)*/
            diag_elem -= L_elem*temp; /*D(col, col) -= L(col,row)*D(row,row)*L(col,row)*/
            rhs[row] = 0;
            
            /* Gaussian elimination */  
            for (index = L->p[row]; index < col_pointers[row]; index++) // only satisfied when fill-in is introduced
            {
                rhs[L->i[index]] -= L->x[index]*temp;
            }
            index = col_pointers[row];
            col_pointers[row]++;
            L->i[index] = col;
            L->x[index] = L_elem;
        }

        /* Lets Try something */
                
        /*
        1. -> loop over every possible addition and mark non-zero pattern
        */

        for(index = Mppt->p[col]+1; index < Mppt->p[col+1]; index++)
            added[Mppt->i[index]] = 1;
        for (index_in_pattern = start; index_in_pattern < ncol; index_in_pattern++)
        {
            row =  col_pointers[pattern[index_in_pattern]];
            for(index = C->p[pattern[index_in_pattern]]; index < C->p[pattern[index_in_pattern]+1]; index++)
            {
                if (C->i[index] > col+1)
                    added[C->i[index]] = 1;
            }
        }

        C->p[col+1] = C->p[col];
        row = 0;
        for(index = 0; index < ncol; index++){
            if (added[index])
            {
                C->i[C->p[col] + row] = index;
                C->x[C->p[col] + row] = 0.0;
                C->p[col+1] +=  1;
                added[index] = 0;
                row ++;
            }
        }
        /*
        2 . -> Assigning values to these indices
        */ 
        row = 0;
        for(index = Mppt->p[col]+1; index < Mppt->p[col+1]; index++)
        {
            while ( C->i[C->p[col] + row] < Mppt->i[index] )
                row++;
            C->x[C->p[col] + row] += Mppt->x[index];            
        }
        for (index_in_pattern = start; index_in_pattern < ncol; index_in_pattern++)
        {
            row = 0;
            for(index = C->p[pattern[index_in_pattern]]; index < C->p[pattern[index_in_pattern]+1]; index++)
            {
                if (C->i[index] > col+1){
                    while (C->i[C->p[col] + row] < C->i[index])
                        row++;
                    C->x[C->p[col] + row] -= L->x[col_pointers[pattern[index_in_pattern]]-1]*C->x[index];
                }
            }
        }
        theta = 0.0;
        for (index = C->p[col]; index < C->p[col+1]; index++)
            // if (LADEL_SIGN(theta) == 1)
                theta = LADEL_MAX(theta, LADEL_ABS(C->x[index]));
        
        /*Return FAIL if eigenvalue (close to) zero*/
        error_array[col] = Mppt->x[Mppt->p[col]] /*diag_elem*/;
        diag_elem = LADEL_SIGN(Mppt->x[Mppt->p[col]])*LADEL_MAX(LADEL_ABS(Mppt->x[Mppt->p[col]]), LADEL_MAX(theta*theta/beta, 1e-15));
        // ladel_print("\nCorrect diagonal value %f, Concurring value %f.\n", Mppt->x[Mppt->p[col]], theta*theta/beta);
        // diag_elem = LADEL_MAX(LADEL_ABS(Mppt->x[Mppt->p[col]]), LADEL_MAX(theta*theta/beta, 1e-15));
        error_array[col] -= diag_elem;
        if (LADEL_ABS(diag_elem) < 1e-15) 
        {
            ladel_print("LADEL ERROR: MATRIX (POSSIBLY) NOT FULL RANK (diagonal element of %le)\n", diag_elem);
            return FAIL; 
        }

        row = col+1;
        for(index = C->p[col]; index < C->p[col+1]; index++)
        {
            if (C->i[index] == row)
            {
                Mppt->x[Mppt->p[row]] -= C->x[index]*C->x[index]/diag_elem; 
                row ++;
            }
        }

        D[col] = diag_elem;
        Dinv[col] = 1/diag_elem;
    }
    for (index = 0; index < ncol; index++) L->nz[index] = col_pointers[index] - L->p[index];
    return SUCCESS;
}


ladel_int ladel_ldl_numeric(ladel_sparse_matrix *Mpp, ladel_symbolics *sym, ladel_factor *LD, ladel_work* work, ladel_double* error_array, ladel_double beta)
{
    ladel_diag d;
    d.diag_elem = 0;
    d.diag_size = 0;
    if (error_array) return ladel_ldl_numeric_with_modification(Mpp, d, sym, LD, work, error_array, beta);
    else return ladel_ldl_numeric_with_diag(Mpp, d, sym, LD, work);
}