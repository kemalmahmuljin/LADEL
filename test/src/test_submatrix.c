#include "minunit.h"
#include "ladel.h"

#define NROW 4
#define NCOL 5
#define NZMAX 9
#define TOL 1e-8

static ladel_sparse_matrix *M;
static ladel_sparse_matrix *M2;

/*
M = [1.2  0  0    0    0.5;
     0   -2  1.1  0    0;
     3.6  0  1.5  0  0;
     0   -3  0    1.7 -0.5;]
*/

/*
M2 = [          M           ;
      0,5  0    0    -0,7 0 ]

*/

void submatrix_suite_setup(void)
{
    M = ladel_sparse_alloc(NROW, NCOL, NZMAX, UNSYMMETRIC, TRUE, FALSE);
    M->p[0] = 0; M->p[1] = 2; M->p[2] = 4; M->p[3] = 6; M->p[4] = 7; M->p[5] = 9;
    M->i[0] = 0; M->i[1] = 2; M->i[2] = 1; M->i[3] = 3; M->i[4] = 1; M->i[5] = 2; M->i[6] = 3; M->i[7] = 0; M->i[8] = 3;
    M->x[0] = 1.2; M->x[1] = 3.6; M->x[2] = -2; M->x[3] = -3; M->x[4] = 1.1; M->x[5] = 1.5; M->x[6] = 1.7; M->x[7] = 0.5; M->x[8] = -0.5;


    M2 = ladel_sparse_alloc(NCOL, NCOL, NZMAX+2, UNSYMMETRIC, TRUE, FALSE);
    M2->p[0] = 0; M2->p[1] = 3; M2->p[2] = 5; M2->p[3] = 7; M2->p[4] = 9; M2->p[5] = 11;
    M2->i[0] = 0; M2->i[1] = 2; M2->i[2] = 4; 
    M2->i[3] = 1; M2->i[4] = 3; 
    M2->i[5] = 1; M2->i[6] = 2; 
    M2->i[7] = 3; M2->i[8] = 4; M2->i[9] = 0; M2->i[10] = 3;

    M2->x[0] = 1.2; M2->x[1] = 3.6; M2->x[2] = 0.5; 
    M2->x[3] = -2; M2->x[4] = -3; 
    M2->x[5] = 1.1; M2->x[6] = 1.5; 
    M2->x[7] = 1.7; M2->x[8] = -0.7; M2->x[9] = 0.5; M2->x[10] = -0.5;
}

void submatrix_suite_teardown(void)
{
    ladel_sparse_free(M);
    ladel_sparse_free(M2);
    
}

MU_TEST(test_column_submatrix)
{
    ladel_int cols[3] = {0, 2, 4};
    ladel_int nb_cols = 3, pattern = 0;
    ladel_sparse_matrix *M_sub = ladel_column_submatrix(M, cols, nb_cols);

    ladel_int p = 2;
    ladel_sparse_matrix *M_sub2 = ladel_leading_principal_submatrix(M2, p, pattern);





    mu_assert_long_eq(M_sub2->ncol, 3);
    mu_assert_long_eq(M_sub2->nrow, 3);
    
    
    mu_assert_long_eq(M_sub->ncol, 3);
    mu_assert_long_eq(M_sub->nrow, NROW);
    mu_assert_long_eq(M_sub->nzmax, 6);

    mu_assert_long_eq(M_sub->p[0], 0);
    mu_assert_long_eq(M_sub->p[1], 2);
    mu_assert_long_eq(M_sub->p[2], 4);
    mu_assert_long_eq(M_sub->p[3], 6);

    mu_assert_long_eq(M_sub->i[0], 0);
    mu_assert_long_eq(M_sub->i[1], 2);
    mu_assert_long_eq(M_sub->i[2], 1);
    mu_assert_long_eq(M_sub->i[3], 2);
    mu_assert_long_eq(M_sub->i[4], 0);
    mu_assert_long_eq(M_sub->i[5], 3);

    mu_assert_double_eq(M_sub->x[0], 1.2, TOL);
    mu_assert_double_eq(M_sub->x[1], 3.6, TOL);
    mu_assert_double_eq(M_sub->x[2], 1.1, TOL);
    mu_assert_double_eq(M_sub->x[3], 1.5, TOL);
    mu_assert_double_eq(M_sub->x[4], 0.5, TOL);
    mu_assert_double_eq(M_sub->x[5], -0.5, TOL);

    ladel_sparse_free(M_sub);
}

MU_TEST_SUITE(suite_submatrix)
{
    MU_SUITE_CONFIGURE(submatrix_suite_setup, submatrix_suite_teardown, NULL, NULL);
    MU_RUN_TEST(test_column_submatrix);
}