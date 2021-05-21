#include "ladel_types.h"
#include "ladel_global.h"
#include "ladel_constants.h"
#include "ladel_pattern.h"
#include "ladel_transpose.h"
#include "ladel_permutation.h"
#include "ladel_debug_print.h"
#include <stdio.h>

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

ladel_int ladel_ldl_numeric_with_modification(ladel_sparse_matrix *Mpp, ladel_diag d, ladel_symbolics *sym, ladel_factor *LD, ladel_work* work, ladel_double beta, ladel_int n)
{
     if (!Mpp || !sym || !LD || !work) return FAIL;
    
    ladel_int row, col, index, ncol = Mpp->ncol, start, index_in_pattern, active;
    ladel_int *pattern = sym->pattern;
    ladel_double diag_elem, temp, L_elem;
    
    ladel_sparse_matrix *L = LD->L;
    ladel_double *D = LD->D, *Dinv = LD->Dinv, *error_array = LD->E;
    ladel_double *rhs = work->array_double_all_zeros_ncol1;
    ladel_int *col_pointers = work->array_int_ncol1;

    L->p[0] = col_pointers[0] = 0;
    for (index = 1; index < ncol; index++) 
        L->p[index] = col_pointers[index] = sym->col_counts[index-1];

    L->p[ncol] = sym->col_counts[ncol-1];
    for (col = 0; col < ncol; col++)
    {
        active = 1;
        LADEL_FOR(index, Mpp, col){
            rhs[Mpp->i[index]] = Mpp->x[index];
        }
        diag_elem = rhs[col];
        if (LD->p && LD->p[col] >= n && LADEL_ABS(diag_elem - 1) < LADEL_E_MACH) active = 0;
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

        if ((LD->p) && (2*(LD->p[col] < n) - 1) * diag_elem < LADEL_E_MACH*LADEL_MAX(beta, 1) && active /*LADEL_E_MACH *LADEL_MAX(beta, 1)*/)
        {
        // if ( (LADEL_ABS(diag_elem) < LADEL_E_MACH*LADEL_MAX(beta, 1)) \
        //       || (diag_elem < 0 && LD->p && LD->p[col] < n)   \
        //       || (diag_elem > 0 && LD->p && LD->p[col] >= n) ){
    
            // ladel_print("\nReg %le, orig %le, pos %ld and col %ld\n",LD->DEF_REG, diag_elem, sym->pinv[col], col);
            // fprintf(fp,"At col %ld (original is %ld) (pinv is %ld), subtract from row %ld (original is %ld) (pinv is %ld) , current value %.16le, next value %.16le ;\n", col, sym->p[col], sym->pinv[col], row, sym->p[row], sym->pinv[row], diag_elem, diag_elem - L_elem*temp);
            // ladel_print("\nWe have a problem at col %d with true ind %d old value is %le\n ", col, LD->p[col]+1, diag_elem);


            // if (diag_elem < 0 && LD->p[col] < n)
            //     ladel_print("\nShould be positive at col %d with true ind %d old value is %le\n ", col, LD->p[col]+1, diag_elem);
            // if (diag_elem > 0 && LD->p[col] >= n)
            //     ladel_print("\nShould be negative at col %d with true ind %d old value is %le\n ", col, LD->p[col]+1, diag_elem);
            
            // if (col == 453){
                // ladel_print("\n\n\n\nCol %d with true ind %d old value is %le\n ", col, LD->p[col]+1, diag_elem);
            // }
            
            // diag_elem = LADEL_SIGN(Mppt->x[Mppt->p[col]])*LADEL_MAX(LADEL_ABS(Mppt->x[Mppt->p[col]]), LADEL_MAX(theta*theta/beta, 1e-9));
            if(LD->p[col] >= n){

                //error_array[col] = -diag_elem;
                // if (diag_elem > 0){
                //     ladel_print("\nWe have a problem at col %d %d\n ", col, sym->pinv[col]);   // shouldnt be an issue
                // }
                // ladel_print("\nshould belong to m/p region col %ld and original %ld with diag val being %le", col, LD->p[col], deletet);
                // ladel_print("\n start %d ncol %d\n", start, ncol);
                // if (/* && */diag_elem > 0)
                // {
                    // ladel_print_dense_int_vector_matlab(sym->p, ncol);
                    error_array[LD->p[col]] = -(LD->REG + diag_elem);
                    diag_elem = -LD->REG/*LADEL_SQRT_E_MACH*/;
                    // ladel_print("\nWrong sign col %d converted %d and difference is %d\n ", col, LD->p[col],  ncol-start);   // shouldnt be an issue
                // }
                // else{
                //     error_array[LD->p[col]] = -(LD->REG + diag_elem);
                //     diag_elem = -LD->REG/*LADEL_SQRT_E_MACH*/;
                // }
                // error_array[sym->p[col]] = -LD->DEF_REG/*LADEL_SQRT_E_MACH*/; // to get positive diagonal E
            }
            else{
                // 
                // ladel_print("\nCol\t%ld\t perm \t%ld\tand diag \t%le\t", col, sym->p[col], diag_elem);

                error_array[LD->p[col]] = LD->REG - diag_elem;//LD->DEF_REG;
                diag_elem = LD->REG;//LD->DEF_REG;

                
            }
        
            // ladel_print("\nTrue value of alternative is: %le, value of beta %le and value of theta %le \n", theta*theta/beta, beta, theta);
            // if (LADEL_ABS(Mppt->x[Mppt->p[col]]) < theta*theta/beta){
            //     ladel_print("\nAlternative was finally chosen !!!!!\n");
            // }
            // ladel_print("\nChosen diagelem %le, and spec lim %le prev %le", diag_elem, theta*theta/beta, Mppt->x[Mppt->p[col]]);
        }
        // ladel_print("\n\n\n\nCol %d with true ind %d old value is %le\n ", col, sym->p[col]+1, diag_elem);
        /*Return FAIL if eigenvalue (close to) zero*/
        if (LADEL_ABS(diag_elem) < 1e-15) 
        {
            // ladel_print("\nAt column %d\n", sym->p[col]);
            ladel_print("\nWe have a problem at col %d %d\n ", col, LD->p[col]);
            ladel_print("LADEL ERROR: MATRIX (POSSIBLY) NOT FULL RANK (diagonal element of %le)\n", diag_elem);
            return FAIL; 
        }


// ladel_print("\n\n\n\nCol %d with true ind %d old value is %le\n ", col, LD->p[col]+1, diag_elem);
        D[col] = diag_elem;
        Dinv[col] = 1/diag_elem;
    }
    // ladel_print("\nFinal value %le" , D[955]);
    // fclose(fp);

    for (index = 0; index < ncol; index++) L->nz[index] = col_pointers[index] - L->p[index];
    

    //   FILE *fp; 
    //     fp = fopen("../symp.m","w");
        
    //     ladel_print_dense_int_vector_matlab_to_file(fp,sym->p, sym->ncol);
    //     ladel_print("yu|||||||||||||");
    //     fclose(fp);
    //     fp = fopen("../sympinv.m","w");
    //     ladel_print_dense_int_vector_matlab_to_file(fp,sym->pinv, sym->ncol);
    //     fclose(fp);


    
    // ladel_double *local_err = (ladel_double *) ladel_malloc(ncol, sizeof(ladel_double)); 

    // ladel_permute_vector(error_array, sym->p, sym->ncol, local_err);
    // error_array = local_err;
    // ladel_print_dense_vector_matlab(error_array, ncol);
    // ladel_free(local_err);
    // ladel_print("\nCount %ld and dCount %ld\n", count, dcount);
    return SUCCESS;
}

// ladel_int ladel_ldl_numeric_with_modification(ladel_sparse_matrix *Mpp, ladel_diag d, ladel_symbolics *sym, ladel_factor *LD, ladel_work* work, ladel_double* error_array, ladel_double beta, ladel_int n)
// {
//     if (!Mpp || !sym || !LD || !work) return FAIL;
    
//     ladel_int row, col, index, ncol = Mpp->ncol, start, index_in_pattern;
//     ladel_int *pattern = sym->pattern;
//     ladel_double diag_elem, temp, L_elem;
//     ladel_int *added = (ladel_int *) ladel_calloc(ncol, sizeof(ladel_int)); // should be initialised with 0's
//     ladel_int *lastrow = (ladel_int *) ladel_calloc(ncol, sizeof(ladel_int)); 
//     ladel_int *lastrowindex = (ladel_int *) ladel_calloc(ncol, sizeof(ladel_int));
//     ladel_double theta;
//     ladel_sparse_matrix *L = LD->L;
//     ladel_double *D = LD->D, *Dinv = LD->Dinv;
//     ladel_double *rhs = work->array_double_all_zeros_ncol1;
//     ladel_int *col_pointers = work->array_int_ncol1;
//     ladel_sparse_matrix *Mppt = ladel_transpose(Mpp, TRUE, work);

// // ladel_print_sparse_matrix_matlab(Mppt);
// // ladel_print_dense_int_vector_matlab(sym->p, ncol);
// // ladel_print_dense_int_vector_matlab(sym->pinv, ncol);
// /*
//  print matrices to file
// */

//    FILE *fptrM;

// //    fptrM = fopen("../Matr.m","w");
// //    ladel_print_sparse_matrix_matlab_to_file(fptrM, Mpp);
// //    fclose(fptrM);
// //    fptrM = fopen("../MatrT.m","w");
// //    ladel_print_sparse_matrix_matlab_to_file(fptrM, Mppt);
// //    fclose(fptrM);


//     // copy column information  -- made using Mbasis --
//     L->p[0] = col_pointers[0] = 0;
//     for (index = 1; index < ncol; index++) 
//     {
//         L->p[index] = col_pointers[index] = sym->col_counts[index-1];
//     }
//     L->p[ncol] = sym->col_counts[ncol-1];

// //    fptrM = fopen("../Lp.m","w");
// //    ladel_print_dense_int_vector_matlab_to_file(fptrM, L->p, ncol);
// //    fclose(fptrM);

//     // ladel_print("\ninitial value\n");
//     // ladel_print_dense_int_vector_matlab(col_pointers, ncol);
//     for (col = 0; col < ncol; col++)
//     {
//      //***/ladel_print("\n#######################  %d  ##########################\n", col);
//         start = ladel_nonzero_pattern_of_row_in_L(Mpp, sym, col);
//         /* Lets Try something */
                
//         /*
//         1. -> loop over every possible addition and mark non-zero pattern
//         */
        
//         for(index = Mppt->p[col]+1; index < Mppt->p[col+1]; index++)
//         {
//             added[Mppt->i[index]] = 1;
//           //***/ladel_print("\nRow %d and col %d, total col elems %d, value is %f\n", Mppt->i[index], col, Mppt->p[col+1] -(Mppt->p[col]+1), Mppt->x[index]);
//         }
//         for (index_in_pattern = start; index_in_pattern < ncol; index_in_pattern++) // corresponds to number of nonzero L's in row
//         {
//             row =  pattern[index_in_pattern];
//           //***/ladel_print("\nEntered index %d, row %d associated L value %f and start %d\n", col_pointers[row], row, L->x[col_pointers[row]], start);
//             // for(index = col_pointers[row]+1; index <  L->p[row+1]; index++) //old loop
//             for(index = col_pointers[row]+1; index <  lastrowindex[row]+1; index++)
//             {   
//                 //if (L->i[index] > L->i[index-1] && L->i[index] <= lastrow[row]/*+1*/) /* only here for the colcounts extra zeros < dont know why they're there, something about symbolic factorization > */
//                 //{
//                   //**/ladel_print("\nRow %d and col %d and value %f , index %d??, %d\n", L->i[index], col, L->x[index], index, L->i[index] > col);
//                     added[L->i[index]] = 2;
//                 //}
//             }
//         }
//         row = 0;
//         for(index = col+1; index < ncol; index++){
//             if (added[index])
//             {
//                 L->i[L->p[col] + row] = index;
//                 // ladel_print("\nIndex: %d \n", index);
//                 L->x[L->p[col] + row] = 0.0;
//                 if (added[index] == 1)
//                 {
//                     added[index] = 0;       // -----> try
//                 }
//                 lastrow[col] = index;
//                 lastrowindex[col] = L->p[col] + row;
                
//                 //***/ladel_print("\nColumn %d and row %d of L need to have nonzero\n ", col, index);
//                 row ++;
//             }
//         }
//         //2 . -> Assigning values to these indices
//         for (index_in_pattern = start; index_in_pattern < ncol; index_in_pattern++) // number of nonzero L's = ncol - start
//         {
//             //***/ladel_print("\nSOME VALUES: L index: %d and its value %f, Diag index %d and its value %f\n", L->i[col_pointers[pattern[index_in_pattern]]],L->x[col_pointers[pattern[index_in_pattern]]], Mppt->p[pattern[index_in_pattern]], Mppt->x[Mppt->p[pattern[index_in_pattern]]]);
//             L->x[col_pointers[pattern[index_in_pattern]]] /= Mppt->x[Mppt->p[pattern[index_in_pattern]]];
//             row = 0;
//             // for(index = col_pointers[pattern[index_in_pattern]]+1; index <  L->p[pattern[index_in_pattern]+1]; index++) //old loop
//             for(index = col_pointers[pattern[index_in_pattern]]+1; index <  lastrowindex[pattern[index_in_pattern]]+1; index++)
//             {
//                 // /***/ladel_print("\nL row %d with col %d\n",L->i[index], col);
//                 // if (L->i[index] > L->i[index-1] && L->i[index] <= lastrow[pattern[index_in_pattern]]/*+1*/)//if (L->i[index] > col/*+1*/)
//                 // {
//                     // if (added[L->i[index]]==2) // not neccesary
//                     // {
//                         while (L->i[L->p[col] + row] < L->i[index])
//                             row++;
//                         added[L->i[index]]=0; //------>>> try
//                         // ladel_print("\nAffected row is %d \n", L->i[index]);
                   
//                     //ladel_print("\nL value to multiply with %f, colptr %d\n", L->x[col_pointers[pattern[index_in_pattern]]],col_pointers[pattern[index_in_pattern]]);
//                     L->x[L->p[col] + row] -= L->x[ col_pointers[pattern[index_in_pattern]] ] * L->x[index];
//                     row++;
//                     // }
//                     // else
//                     // {
//                     //     ladel_print("\nHow even...\n");
//                     // }
//                 // }
//             //***/ladel_print("\nColpointer associated with column %d is %d\nOld row %d, new row %d\nIndex previous %d, next %d\n\n", pattern[index_in_pattern], index, L->i[index], L->i[index+1], L->p[col], L->p[col+1]);
//             }
//             col_pointers[pattern[index_in_pattern]]++;
//         }

//         row = 0;
//         for(index = Mppt->p[col]+1; index < Mppt->p[col+1]; index++)
//         {
//             while ( L->i[L->p[col] + row] < Mppt->i[index] )
//                 row++;
//             L->x[L->p[col] + row] += Mppt->x[index];
//             row++;
//         }
// //    fptrM = fopen("../col_pointers.m","a");
// //    ladel_print_dense_int_vector_matlab_to_file(fptrM, L->i, L->p[col]);
// //    fclose(fptrM);

// // //    -- old diagonal adjustment

//         // theta = LADEL_ABS(L->x[L->p[col+1] -1]); // last entry
//         // //***/ladel_print("\nDumbass row %d, value %f\n", L->i[L->p[col+1] -1], L->x[L->p[col+1] - 1]);
//         // for (index = L->p[col+1]-2; index > L->p[col]-1; index--)
//         // {
//         //     // if (LADEL_SIGN(theta) == 1)
//         //     // ladel_print("\nDumbass row %d, value %f\n", L->i[index], L->x[index]);
//         //         theta = LADEL_MAX(theta, LADEL_ABS(L->x[index]));
//         // }

// theta = 0;
// index = L->p[col]+1;
// // ladel_print("\nRow values %d vs %d \n", L->i[index], lastrow[col]);
// // while (L->i[index] <= lastrow[col])
//         //ladel_print("\nCOL %d, row %d, index %d value %le , lastrow index %d and row %d \n", col, L->i[index], index, L->x[index], lastrowindex[col], L->i[lastrowindex[col]]);

// for (index = L->p[col]+1; index < lastrowindex[col]+1; index++)
// {
//     if((L->i[index] > ncol) ||(L->i[index] < 0))
//     {
//         ladel_print("\nCOL %d, row %d, value %le , lastrow index %d and row %d first  \n", col, L->i[index], L->x[index], lastrowindex[col], L->i[lastrowindex[col]]);
//         return FAIL;
//     }
//     // ladel_print("\nCOL %d, row %d, value %le \n", col, L->i[index], L->x[index]);
//     theta = LADEL_MAX(theta, LADEL_ABS(L->x[index]));
// }
//     // ladel_print("\nFirst index %d, last index %d with rows %d and %d\n", L->p[col]+1, lastrowindex[col], L->i[L->p[col]+1], L->i[lastrowindex[col]]);
// // while (index < lastrowindex[col]+1)
// // {
// //     theta = LADEL_MAX(theta, LADEL_ABS(L->x[index]));
// //     index++;
// // }
//         if (col == ncol-1)
//         { 
//             theta = 0;
//         }
//         // ladel_print("\nValue of theta %f\n", theta);
//         /*Return FAIL if eigenvalue (close to) zero*/
 
//         error_array[col] = -Mppt->x[Mppt->p[col]];
//         /***/ladel_print("\nCurrent diag elem %f \n", Mppt->x[Mppt->p[col]]);
 
//    fptrM = fopen("../errors.txt","a");
//    fprintf(fptrM, "BEFORE(%ld)= %le;",col,  Mppt->x[Mppt->p[col]]);
//    fclose(fptrM);
//         if ( (LADEL_ABS(Mppt->x[Mppt->p[col]]) < 1e-15) /*|| (SUCCESS)*/ ){
//             ladel_print("\ndd %f \n",Mppt->x[Mppt->p[col]]);
//             diag_elem = LADEL_SIGN(Mppt->x[Mppt->p[col]])*LADEL_MAX(LADEL_ABS(Mppt->x[Mppt->p[col]]), LADEL_MAX(theta*theta/beta, 1e-9));
//             // if(sym->p[col] >= n){
//             //     Mppt->x[Mppt->p[col]] -= 1e-7;
//             //     // ladel_print("\nyo\n");
//             //     error_array[col] = -1e-7;
//             // }
//             // else{
//             //     Mppt->x[Mppt->p[col]] += 5.6e-12;
//             //     // ladel_print("\nno\n");
//             //     error_array[col] = -5e-12;
//             // }
//             // ladel_print("\nTrue value of alternative is: %le, value of beta %le and value of theta %le \n", theta*theta/beta, beta, theta);
//             // if (LADEL_ABS(Mppt->x[Mppt->p[col]]) < theta*theta/beta){
//             //     ladel_print("\nAlternative was finally chosen !!!!!\n");
//             // }
//             ladel_print("\nChosen diagelem %le, and spec lim %le prev %le", diag_elem, theta*theta/beta, Mppt->x[Mppt->p[col]]);
//         }
//    fptrM = fopen("../errors.txt","a");
//    fprintf(fptrM, "\tBEFORE(%ld)= %le;\n",col,  Mppt->x[Mppt->p[col]]);
// //    ladel_print_dense_int_vector_matlab_to_file(fptrM, Mppt->x[Mppt->p[col]], 1);
//    fclose(fptrM);
//         // else
//         // {
//         //     diag_elem = Mppt->x[Mppt->p[col]];
//         //     ladel_print("\nDEFAULT Chosen diagelem %f", diag_elem);

//         // }
//         // 
//         error_array[col] += Mppt->x[Mppt->p[col]];
//         // if(error_array[col] > 1e-10)
//         //     ladel_print("\ne-value %f", error_array[col]);
//         // ladel_print("%d %f\n", col, Mppt->x[Mppt->p[col]]);

//         diag_elem = Mppt->x[Mppt->p[col]] ;
 
// // //  end old diagonal adjustement */


 
//         if (LADEL_ABS(diag_elem) < 1e-15) 
//         {
//             ladel_print("LADEL ERROR: MATRIX (POSSIBLY) NOT FULL RANK (diagonal element of %le)\n", diag_elem);
//             return FAIL; 
//         }
//         //***/ladel_print("\nSHow me the col starts %d associated row %d, %d\n", L->p[col], L->i[L->p[col]], L->i[0]);
//         /////////////////// Temporarely commented out //////////
//         // Updating diagonal entries // 
//         //ladel_print_dense_int_vector_matlab(L->p, ncol);
//         // for(index = L->p[col]; index < L->p[col+1]; index++) //old loop
//         for(index = L->p[col]; index < lastrowindex[col]+1; index++)
//         {
//             // ladel_print("\n%d, %f, %d", index, L->x[index], L->i[index]);
//             // if (L->i[index] < 0 || L->i[index] >ncol)
//             // {
//             //     // for (int i=0; i<ncol;i++){
//             //     //     if (i%5 == 0)
//             //     //     {
//             //     //         ladel_print("\n");
//             //     //     }
//             //     //     ladel_print("%d -", L->p[index]);
//             //     // }
//             //     ladel_print("\nLadel index %d \n", index);
//             //     ladel_print("\nPossibly error found %d, start is %d while end is %d and col %d of total col %d\n", L->i[index],L->p[col],L->p[col+1], col+1, ncol);
//             //     return FAIL; 
//             // }
//             // if (L->i[index] > col && L->i[index] < lastrow[col])
//             // {
//                 //***/ladel_print("\nRow of L %d, and corresponding value %f", L->i[index],L->x[index]);
//                 Mppt->x[Mppt->p[L->i[index]]] -= L->x[index]*L->x[index]/diag_elem; 
//             // }
//         }
//         /////////////////// Temporarely commented out //////////

//         // row = 0;
//         // for(index = col+1; index < ncol; index++)
//         // {

//         //     if (added[index])
//         //     {
//         //         while ( L->i[L->p[col] + row] < index )
//         //             row++;

//         //         //if (L->i[L->p[col] + row] < 0 || L->i[L->p[col] + row] >ncol)
//         //         ladel_print("\nUFF row used: %d added row: %d", L->i[L->p[col] + row], index);
//         //         Mppt->x[Mppt->p[index]] -= L->x[L->p[col] + row]*L->x[L->p[col] + row]/diag_elem;
//         //         row++;
//         //         added[index] = 0;
//         //     }
//         // }
//         D[col] = diag_elem;
//         Dinv[col] = 1/diag_elem;
//     }
//     // ladel_print_dense_int_vector_matlab(L->i, L->p[ncol]);
//     for (index = 0; index < ncol; index++) L->nz[index] = col_pointers[index] - L->p[index];
//     // ladel_print_dense_int_vector_matlab(sym->pinv, ncol);
//     // ladel_print_dense_int_vector_matlab(lastrow, ncol);
// ladel_free(added); 
// ladel_free(lastrow);
// ladel_free(lastrowindex);
// ladel_sparse_free(Mppt);
    
//     //ladel_print("\n\n\n%d, %d ,%d, %d, %f\n\n\n\n", L->ncol, L->nrow, L->p[ncol], L->i[4], L->x[4]);
//     // ladel_double *hhh = L->x;
//     // ladel_int *iii = L->i;
//     // ladel_print_dense_vector_matlab(hhh, 6);
//     // ladel_print_dense_int_vector_matlab(iii, 6);
//     // ladel_print_dense_int_vector_matlab(L->i, L->p[ncol]);
    
//     return SUCCESS;
// }

ladel_int ladel_ldl_numeric(ladel_sparse_matrix *Mpp, ladel_symbolics *sym, ladel_factor *LD, ladel_work* work , ladel_double beta, ladel_int n)
{
    ladel_diag d;
    d.diag_elem = 0;
    d.diag_size = 0;
    if (LD->E) return ladel_ldl_numeric_with_modification(Mpp, d, sym, LD, work, beta, n);
    else return ladel_ldl_numeric_with_diag(Mpp, d, sym, LD, work);
}




// ladel_int ladel_ldl_numeric_with_modification(ladel_sparse_matrix *Mpp, ladel_diag d, ladel_symbolics *sym, ladel_factor *LD, ladel_work* work, ladel_double* error_array, ladel_double beta, ladel_int n)
// {
//     if (!Mpp || !sym || !LD || !work) return FAIL;
    
//     ladel_int row, col, index, ncol = Mpp->ncol, start, index_in_pattern;
//     ladel_int *pattern = sym->pattern;
//     ladel_double diag_elem, temp, L_elem;
//     ladel_int *added = (ladel_int *) ladel_calloc(ncol, sizeof(ladel_int)); // should be initialised with 0's
//     ladel_int *lastrowindex = (ladel_int *) ladel_calloc(ncol, sizeof(ladel_int));
//     ladel_double theta;
//     ladel_sparse_matrix *L = LD->L;
//     ladel_double *D = LD->D, *Dinv = LD->Dinv;
//     ladel_int *col_pointers = work->array_int_ncol1;
//     ladel_sparse_matrix *Mppt = ladel_transpose(Mpp, TRUE, work);

// /*
//  print matrices to file
// */
//         ladel_print("\n col: %d\n", col);

// //    FILE *fptrM;

// //    fptrM = fopen("../Matr.m","w");
// //    ladel_print_sparse_matrix_matlab_to_file(fptrM, Mpp);
// //    fclose(fptrM);
// //    fptrM = fopen("../MatrT.m","w");
// //    ladel_print_sparse_matrix_matlab_to_file(fptrM, Mppt);
// //    fclose(fptrM);


//     // copy column information  -- made using Mbasis --
//     L->p[0] = col_pointers[0] = 0;
//     for (index = 1; index < ncol; index++) 
//     {
//         L->p[index] = col_pointers[index] = sym->col_counts[index-1];
//     }
//     L->p[ncol] = sym->col_counts[ncol-1];


//     for (col = 0; col < ncol; col++)
//     {
//      //***/ladel_print("\n#######################  %d  ##########################\n", col);
//         start = ladel_nonzero_pattern_of_row_in_L(Mpp, sym, col);
//         /* Lets Try something */
                
//         /*
//         1. -> loop over every possible addition and mark non-zero pattern
//         */
        
//         for(index = Mppt->p[col]+1; index < Mppt->p[col+1]; index++)
//         {
//             added[Mppt->i[index]] = 1;
//           //***/ladel_print("\nRow %d and col %d, total col elems %d, value is %f\n", Mppt->i[index], col, Mppt->p[col+1] -(Mppt->p[col]+1), Mppt->x[index]);
//         }
//         for (index_in_pattern = start; index_in_pattern < ncol; index_in_pattern++) // corresponds to number of nonzero L's in row
//         {
//             row =  pattern[index_in_pattern];
//           //***/ladel_print("\nEntered index %d, row %d associated L value %f and start %d\n", col_pointers[row], row, L->x[col_pointers[row]], start);
//             // for(index = col_pointers[row]+1; index <  L->p[row+1]; index++) //old loop
//             for(index = col_pointers[row]+1; index <  lastrowindex[row]+1; index++)
//             {   
//                 added[L->i[index]] = 2;
//             }
//         }
//         row = 0;
//         for(index = col+1; index < ncol; index++){
//             if (added[index])
//             {
//                 L->i[L->p[col] + row] = index;
//                 L->x[L->p[col] + row] = 0.0;
//                 added[index] = 0;       // -----> try
//                 //***/ladel_print("\nColumn %d and row %d of L need to have nonzero\n ", col, index);
//                 lastrowindex[col] = L->p[col] + row;                

//                 row ++;
//             }
//         }
//         // lastrowindex[col] = L->p[col] + row-1;                

//         //2 . -> Assigning values to these indices
//         for (index_in_pattern = start; index_in_pattern < ncol; index_in_pattern++) // number of nonzero L's = ncol - start
//         {
//             //***/ladel_print("\nSOME VALUES: L index: %d and its value %f, Diag index %d and its value %f\n", L->i[col_pointers[pattern[index_in_pattern]]],L->x[col_pointers[pattern[index_in_pattern]]], Mppt->p[pattern[index_in_pattern]], Mppt->x[Mppt->p[pattern[index_in_pattern]]]);
//             L->x[col_pointers[pattern[index_in_pattern]]] /= Mppt->x[Mppt->p[pattern[index_in_pattern]]];
//             row = 0;
//             // for(index = col_pointers[pattern[index_in_pattern]]+1; index <  L->p[pattern[index_in_pattern]+1]; index++) //old loop
//             for(index = col_pointers[pattern[index_in_pattern]]+1; index <  lastrowindex[pattern[index_in_pattern]]+1; index++)
//             {
//                 while (L->i[L->p[col] + row] < L->i[index])
//                     row++;
//                 L->x[L->p[col] + row] -= L->x[ col_pointers[pattern[index_in_pattern]] ] * L->x[index];
//                 row++;
//             }
//             col_pointers[pattern[index_in_pattern]]++;
//         }

//         row = 0;
//         for(index = Mppt->p[col]+1; index < Mppt->p[col+1]; index++)
//         {
//             while ( L->i[L->p[col] + row] < Mppt->i[index] )
//                 row++;
//             L->x[L->p[col] + row] += Mppt->x[index];
//             row++;
//         }

//         theta = 0;
//         // index = L->p[col]+1;

//         for (index = L->p[col]+1; index < lastrowindex[col]+1; index++)
//         {
//             if((L->i[index] > ncol) ||(L->i[index] < 0))
//             {
//                 ladel_print("\nCOL %d, row %d, value %le , lastrow index %d and row %d first  \n", col, L->i[index], L->x[index], lastrowindex[col], L->i[lastrowindex[col]]);
//                 return FAIL;
//             }
//             // ladel_print("\nCOL %d, row %d, value %le \n", col, L->i[index], L->x[index]);
//             theta = LADEL_MAX(theta, LADEL_ABS(L->x[index]));
//         }
//         if (col == ncol-1)
//         { 
//             theta = 0;
//         }

//         if ( (LADEL_ABS(Mppt->x[Mppt->p[col]]) < 1e-15) /*|| (SUCCESS)*/ )
//         {
//             error_array[col] = -Mppt->x[Mppt->p[col]];
//             diag_elem = LADEL_SIGN(Mppt->x[Mppt->p[col]])*LADEL_MAX(LADEL_ABS(Mppt->x[Mppt->p[col]]), LADEL_MAX(theta*theta/beta, 1e-14));
//             ladel_print("\nOriginal %le while other %le\n",  Mppt->x[Mppt->p[col]], theta*theta/beta);
//             Mppt->x[Mppt->p[col]] = diag_elem;

//             // if(sym->p[col] >= n){
//             //     Mppt->x[Mppt->p[col]] -= 1e-7;
//             //     // ladel_print("\nyo\n");
//             //     error_array[col] = -1e-7;
//             // }
//             // else{
//             //     Mppt->x[Mppt->p[col]] += 5.6e-12;
//             //     // ladel_print("\nno\n");
//             //     error_array[col] = 5.6e-12;
//             // }
//         }
//         // diag_elem = Mppt->x[Mppt->p[col]];
//         // Mppt->x[Mppt->p[col]] = diag_elem;
//         error_array[col] += Mppt->x[Mppt->p[col]];

 
//         if (LADEL_ABS(diag_elem) < 1e-15) 
//         {
//             ladel_print("LADEL ERROR: MATRIX (POSSIBLY) NOT FULL RANK (diagonal element of %le)\n", diag_elem);
//             return FAIL; 
//         }
//         for(index = L->p[col]; index < lastrowindex[col]+1; index++)
//         {
//             Mppt->x[Mppt->p[L->i[index]]] -= L->x[index]*L->x[index]/diag_elem; 
//         }
//         D[col] = diag_elem;
//         Dinv[col] = 1/diag_elem;
//     }
//     for (index = 0; index < ncol; index++) L->nz[index] = col_pointers[index] - L->p[index];

//     // Free allocated variables
//     ladel_free(added); 
//     ladel_free(lastrowindex);
//     ladel_sparse_free(Mppt);

//     return SUCCESS;
// }


