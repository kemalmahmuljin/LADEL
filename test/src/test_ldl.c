#include "minunit.h"
#include "ladel_types.h"
#include "ladel_global.h"
#include "ladel_constants.h"
#include "ladel.h"
#include "ladel_debug_print.h"
#include "ladel_ldl_symbolic.h"


#define NCOL 5
#define NROW 5
#define NZMAX 10
#define TOL 1e-10

#define NCOL2 4
#define NROW2 4
#define NZMAX2 5

#define NCOL3 8
#define NROW3 8
#define NZMAX3 19
#define REG1 1E-16

static ladel_work *work;
static ladel_work *work2;
static ladel_work *work3;
static ladel_sparse_matrix *M, *Mbasis;
static ladel_sparse_matrix *M2, *Mnz;
static ladel_sparse_matrix *M3, *Mbasis3;
static ladel_symbolics *sym;
static ladel_factor *LD;
static ladel_symbolics *sym2;
static ladel_factor *LD2;
static ladel_symbolics *sym3;
static ladel_factor *LD3;
void ldl_suite_setup(void)
{
    work = ladel_workspace_allocate(NCOL);
    work2 = ladel_workspace_allocate(NCOL2);
    work3 = ladel_workspace_allocate(NCOL3);
    M = ladel_sparse_alloc(NROW, NCOL, NZMAX, UPPER, TRUE, FALSE);
    M->p[0] = 0; M->p[1] = 1; M->p[2] = 3; M->p[3] = 4; M->p[4] = 7; M->p[5] = 10;
    M->i[0] = 0; M->i[1] = 0; M->i[2] = 1; M->i[3] = 2; M->i[4] = 1; M->i[5] = 2; M->i[6] = 3; M->i[7] = 0; M->i[8] = 3; M->i[9] = 4;  
    M->x[0] = 1; M->x[1] = 10; M->x[2] = 2; M->x[3] = -3; M->x[4] = 11; M->x[5] = 12; M->x[6] = 4; M->x[7] = -3; M->x[8] = 2; M->x[9] = -5;
    Mbasis = ladel_sparse_alloc(NROW, NCOL, 15, UPPER, FALSE, FALSE); /*full matrix*/
    Mbasis->p[0] = 0; Mbasis->p[1] = 1; Mbasis->p[2] = 3; Mbasis->p[3] = 6; Mbasis->p[4] = 10; Mbasis->p[5] = 15;
    Mbasis->i[0] = 0; Mbasis->i[1] = 0; Mbasis->i[2] = 1; Mbasis->i[3] = 0; Mbasis->i[4] = 1; Mbasis->i[5] = 2; Mbasis->i[6] = 0; Mbasis->i[7] = 1; 
    Mbasis->i[8] = 2; Mbasis->i[9] = 3; Mbasis->i[10] = 0; Mbasis->i[11] = 1; Mbasis->i[12] = 2; Mbasis->i[13] = 3; Mbasis->i[14] = 4;

    M2 = ladel_sparse_alloc(NCOL2, NCOL2, NZMAX2, UPPER, TRUE, FALSE);
    M2->p[0] = 0; M2->p[1] = 1; M2->p[2] = 2; M2->p[3] = 3; M2->p[4] = 5; 
    M2->i[0] = 0; M2->i[1] = 1; M2->i[2] = 2; M2->i[3] = 0; M2->i[4] = 3; 
    M2->x[0] = 2; M2->x[1] = 1; M2->x[2] = 3; M2->x[3] = 1; M2->x[4] = 2; 

    /* Adding one value in the col_pointers (from col 3 onwards) */
    Mnz = ladel_sparse_alloc(NROW, NCOL, NZMAX+1, UPPER, TRUE, TRUE);
    Mnz->p[0] = 0; Mnz->p[1] = 1; Mnz->p[2] = 3; Mnz->p[3] = 5; Mnz->p[4] = 8; Mnz->p[5] = 11;
    Mnz->nz[0] = 1; Mnz->nz[1] = 2; Mnz->nz[2] = 1; Mnz->nz[3] = 3; Mnz->nz[4] = 3; 

    Mnz->i[0] = 0; Mnz->i[1] = 0; Mnz->i[2] = 1; Mnz->i[3] = 2; Mnz->i[4] = 52342304; Mnz->i[5] = 1; Mnz->i[6] = 2; Mnz->i[7] = 3; Mnz->i[8] = 0; Mnz->i[9] = 3; Mnz->i[10] = 4;  
    Mnz->x[0] = 1; Mnz->x[1] = 10; Mnz->x[2] = 2; Mnz->x[3] = -3; Mnz->x[4] = 454359; Mnz->x[5] = 11; Mnz->x[6] = 12; Mnz->x[7] = 4; Mnz->x[8] = -3; Mnz->x[9] = 2; Mnz->x[10] = -5;


    /* KKT matrix with 5 primal variables, 2 equality constraints and 1 inequality constraint */
    M3 = ladel_sparse_alloc(NROW3, NCOL3, NZMAX3, UPPER, TRUE, FALSE);
    M3->p[0] = 0; M3->p[1] = 1; M3->p[2] = 3; M3->p[3] = 4; M3->p[4] = 7; M3->p[5] = 10;
    M3->p[6] = 13; M3->p[7] = 16; M3->p[8] = 19;

    M3->i[0] = 0; M3->i[1] = 0; M3->i[2] = 1; M3->i[3] = 2; M3->i[4] = 1; M3->i[5] = 2; M3->i[6] = 3; M3->i[7] = 0; M3->i[8] = 3; M3->i[9] = 4;  
    M3->i[10] = 0; M3->i[11] = 4; M3->i[12] = 5; M3->i[13] = 1; M3->i[14] = 2; M3->i[15] = 6; M3->i[16] = 2; M3->i[17] = 4; M3->i[18] = 7; 
    
    M3->x[0] = 1; M3->x[1] = 10; M3->x[2] = 2; M3->x[3] = 3; M3->x[4] = 11; M3->x[5] = 12; M3->x[6] = 4; M3->x[7] = -3; M3->x[8] = 2; M3->x[9] = -5;
    M3->x[10] = 1; M3->x[11] = 3; M3->x[12] = -REG1; M3->x[13] = 2; M3->x[14] = 4; M3->x[15] = -REG1; M3->x[16] = 1; M3->x[17] = 1; M3->x[18] = -1; 


    Mbasis3 = M3;
    // Mbasis3 = ladel_sparse_alloc(NROW, NCOL, 15, UPPER, FALSE, FALSE); /*full matrix*/
    // Mbasis3->p[0] = 0; Mbasis3->p[1] = 1; Mbasis3->p[2] = 3; Mbasis3->p[3] = 6; Mbasis3->p[4] = 10; Mbasis3->p[5] = 15;
    // Mbasis3->i[0] = 0; Mbasis3->i[1] = 0; Mbasis3->i[2] = 1; Mbasis3->i[3] = 0; Mbasis3->i[4] = 1; Mbasis3->i[5] = 2; Mbasis3->i[6] = 0; Mbasis3->i[7] = 1; 
    // Mbasis3->i[8] = 2; Mbasis3->i[9] = 3; Mbasis3->i[10] = 0; Mbasis3->i[11] = 1; Mbasis3->i[12] = 2; Mbasis3->i[13] = 3; Mbasis3->i[14] = 4;


}

void ldl_suite_teardown(void)
{
    ladel_workspace_free(work);
    ladel_workspace_free(work2);
    ladel_workspace_free(work3);
    ladel_sparse_free(M);
    ladel_sparse_free(Mbasis);
    ladel_sparse_free(M2);
    ladel_sparse_free(Mnz);
    ladel_sparse_free(M3);
    ladel_sparse_free(Mbasis3);
    
}

void ldl_test_setup(void)
{
    sym = ladel_symbolics_alloc(NCOL);  
    sym2 = ladel_symbolics_alloc(NCOL2);
    sym3 = ladel_symbolics_alloc(NCOL3);  
}

void ldl_test_teardown(void)
{
    sym = ladel_symbolics_free(sym);
    sym2 = ladel_symbolics_free(sym2);
    sym3 = ladel_symbolics_free(sym3);
    LD = ladel_factor_free(LD);
    LD2 = ladel_factor_free(LD2);
    LD3 = ladel_factor_free(LD3);
}

MU_TEST(test_simple_ldl)
{
    ladel_double x[NCOL] = {1, 2, 3, 4, 5};
    ladel_double y[NCOL], y_ref[NCOL] = {-1.738103756708408e-01, -1.081932021466905e-01, 4.379964221824687e-01, 3.594991055456172e-01, -7.519141323792485e-01};

    ladel_int status = ladel_factorize(M, sym, NO_ORDERING, 0, &LD, work);
    
    mu_assert_long_eq(status, SUCCESS);
    ladel_dense_solve(LD, x, y, work);
    ladel_int index;

    for (index = 0; index < NCOL; index++)
        mu_assert_double_eq(y[index], y_ref[index], TOL);
}

MU_TEST(test_simple_ldl2)
{
    ladel_double x[NCOL2] = {1, 2, 3, 4};
    ladel_double y[NCOL2], y_ref[NCOL2] = {-2.0/3.0, 2, 1, 7.0/3};

    ladel_int status = ladel_factorize(M2, sym2, NO_ORDERING, 0, &LD2, work2);
    
    mu_assert_long_eq(status, SUCCESS);
    ladel_dense_solve(LD2, x, y, work2);
    ladel_int index;

    for (index = 0; index < NCOL2; index++)
        mu_assert_double_eq(y[index], y_ref[index], TOL);
}

MU_TEST(test_simple_ldl_nz)
{
    ladel_double x[NCOL] = {1, 2, 3, 4, 5};
    ladel_double y[NCOL], y_ref[NCOL] = {-1.738103756708408e-01, -1.081932021466905e-01, 4.379964221824687e-01, 3.594991055456172e-01, -7.519141323792485e-01};

    ladel_int status = ladel_factorize(Mnz, sym, NO_ORDERING, 0, &LD, work);

    mu_assert_long_eq(status, SUCCESS);
    ladel_dense_solve(LD, x, y, work);
    ladel_int index;

    for (index = 0; index < NCOL; index++)
        mu_assert_double_eq(y[index], y_ref[index], TOL);
}

MU_TEST(test_simple_ldl_with_diag)
{
    ladel_sparse_matrix *Q = ladel_sparse_alloc(5, 5, 0, UPPER, TRUE, FALSE);
    Q->p[0] = 0; Q->p[1] = 0; Q->p[2] = 0; Q->p[3] = 0; Q->p[4] = 0; Q->p[5] = 0;
    ladel_diag d;
    d.diag_elem = 2;
    d.diag_size = 5;

    ladel_double x[NCOL] = {1, 2, 3, 4, 5};
    ladel_double y[NCOL];

    ladel_int status = ladel_factorize_with_diag(Q, d, sym, NO_ORDERING, 0, &LD, work);
    mu_assert_long_eq(status, SUCCESS);
    ladel_dense_solve(LD, x, y, work);
    ladel_int index;

    for (index = 0; index < NCOL; index++)
        mu_assert_double_eq(y[index], x[index]/d.diag_elem, TOL);

    Q = ladel_sparse_free(Q);
}

MU_TEST(test_simple_ldl_with_partial_diag)
{
    ladel_diag d;
    d.diag_elem = 2;
    d.diag_size = 2;

    ladel_double x[NCOL] = {1, 2, 3, 4, 5};
    ladel_double y[NCOL], y_ref[NCOL] = {-1.570512820512821e-01,-8.208255159474673e-02,4.177611006879300e-01,3.544402751719825e-01,-7.639931207004377e-01};

    ladel_int status = ladel_factorize_with_diag(M, d, sym, NO_ORDERING, 0, &LD, work);
    mu_assert_long_eq(status, SUCCESS);
    ladel_dense_solve(LD, x, y, work);
    ladel_int index;

    for (index = 0; index < NCOL; index++)
        mu_assert_double_eq(y[index], y_ref[index], TOL);

}

#ifdef DAMD
MU_TEST(test_simple_ldl_with_amd)
{
    ladel_double x[NCOL] = {1, 2, 3, 4, 5};
    ladel_double y[NCOL], y_ref[NCOL] = {-1.738103756708408e-01, -1.081932021466905e-01, 4.379964221824687e-01, 3.594991055456172e-01, -7.519141323792485e-01};

    ladel_int status = ladel_factorize(M, sym, AMD, 0, &LD, work);
    //ladel_print_dense_int_vector_matlab(sym->p, sym->ncol);
    mu_assert_long_eq(status, SUCCESS);
    ladel_dense_solve(LD, x, y, work);
    ladel_int index;

    for (index = 0; index < NCOL; index++)
        mu_assert_double_eq(y[index], y_ref[index], TOL);
    
}
#endif

MU_TEST(test_advanced_ldl)
{
    ladel_double x[NCOL] = {1, 2, 3, 4, 5};
    ladel_double y[NCOL], y_ref[NCOL] = {-1.738103756708408e-01, -1.081932021466905e-01, 4.379964221824687e-01, 3.594991055456172e-01, -7.519141323792485e-01};

    ladel_int status = ladel_factorize_advanced(M, sym, NO_ORDERING, 0, &LD, Mbasis, work, 0, 0);
    
    mu_assert_long_eq(status, SUCCESS);
    ladel_dense_solve(LD, x, y, work);
    ladel_int index;

    for (index = 0; index < NCOL; index++)
        mu_assert_double_eq(y[index], y_ref[index], TOL);
}

MU_TEST(test_advanced_ldl_with_modification)
{
    ladel_double beta = 1, reg = 1e-7;
    
    ladel_int status = ladel_factorize_advanced_with_reg(M3, sym3, GIVEN_ORDERING, 0, &LD3, Mbasis3, work3, beta, 4, reg);

    mu_assert_long_eq(status, SUCCESS);
    // ladel_print_dense_int_vector_matlab(sym3->pinv, sym3->ncol);
    // ladel_print_dense_int_vector_matlab(sym3->p, sym3->ncol);
    // ladel_print("First print changed KKT matrix \n");
    // ladel_print_sparse_matrix_matlab(M3);


    ladel_print("\n\nResulting prints: \n\n");
    ladel_print("Original matrix \n");
    ladel_print_sparse_matrix_matlab(M3);
    ladel_print("Diagonal\n");
    ladel_print_dense_vector_matlab(LD3->D, sym3->ncol);
    ladel_print("Lower triangular\n");
    // ladel_print_sparse_matrix_matlab(LD3->L);
    // ladel_print("Diagonal error\n");
    // ladel_print_dense_vector_matlab(LD3->E, sym3->ncol);
    // ladel_print("\nEnd of prints\n");

}

MU_TEST(test_fixed_part)
{    

    ladel_int num_fixed = 0;
    ladel_print_sparse_matrix_matlab(M3);

    ladel_int status = ladel_factorize_advanced(M3, sym3, GIVEN_ORDERING , num_fixed, &LD3, Mbasis3, work3, 0, 0);
    // mu_assert_long_eq(status, SUCCESS);
    // for (int index = 0; index < NCOL; index++)
    //     mu_assert_int_eq(sym3->p[index], index);
}

MU_TEST_SUITE(suite_ldl)
{
    MU_SUITE_CONFIGURE(ldl_suite_setup, ldl_suite_teardown, ldl_test_setup, ldl_test_teardown);
    // MU_RUN_TEST(test_simple_ldl);
    // MU_RUN_TEST(test_simple_ldl2);
    // MU_RUN_TEST(test_simple_ldl_nz);
    // MU_RUN_TEST(test_simple_ldl_with_diag);
    // MU_RUN_TEST(test_simple_ldl_with_partial_diag);
    // #ifdef DAMD
    // MU_RUN_TEST(test_simple_ldl_with_amd);
    // #endif
    // MU_RUN_TEST(test_advanced_ldl);
    // MU_RUN_TEST(test_fixed_part);
    MU_RUN_TEST(test_advanced_ldl_with_modification);
}