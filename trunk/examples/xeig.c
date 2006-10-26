/*
  Несимметричная проблема собственных чисел
*/

#include "nl.h"

int main()
{
    size_t n = 5;
    double **A, **Q, *scal, *wr, *wi;
    size_t *perm, *iter;
    size_t rc;
    size_t low, high;

    A = nl_dmatrix_create(n, n);
    Q = nl_dmatrix_create(n, n);
    scal = nl_dvector_create(n);
    wr = nl_dvector_create(n);
    wi = nl_dvector_create(n);
    perm = nl_xvector_create(n);
    iter = nl_xvector_create(n);

    A[0][0] = .11; A[0][1] = 1.0; A[0][2] = 18; A[0][3] = 11; A[0][4] = 21; 
    A[1][0] = .10; A[1][1] = 1.4; A[1][2] = 41; A[1][3] = 12; A[1][4] = 13; 
    A[2][0] = .18; A[2][1] = 4.1; A[2][2] = 17; A[2][3] = 23; A[2][4] = 23; 
    A[3][0] = .11; A[3][1] = 1.2; A[3][2] = 23; A[3][3] = 11; A[3][4] = 11; 
    A[4][0] = .21; A[4][1] = 1.3; A[4][2] = 23; A[4][3] = 17; A[4][4] = 17; 

    printf("\nМатрица A\n");
    nl_dmatrix_print(A, n, n, NULL);

    eig_balance(A, n, &low, &high, scal);

    printf("\nСбалансированная матрица:\n");
    nl_dmatrix_print(A, n, n, NULL);

    printf("\nКоэффиценты:\n");
    nl_dvector_print(scal, n, NULL);

    printf("\nlow = %d, \t high = %d\n", low, high);

    eig_hess_reduction(A, n, low, high, perm);

    printf("\nФорма Хессенберга\n");
    nl_dmatrix_print(A, n, n, NULL);

    printf("\nperm:\n");
    nl_xvector_print(perm, n, NULL);

    eig_hess_transform_matrix (A, n, low, high, perm, Q);  // Если нужны собственные векторы

    printf("\nМатрица перехода к форме Хессенберга\n");
    nl_dmatrix_print(Q, n, n, NULL);

    eig_hess(A, n, low, high, wr, wi, 1, Q, iter, &rc);

    printf("\nВещественные части собственных чисел:\n");
    nl_dvector_print(wr, n, NULL);
    printf("\nМнимые части собственных чисел:\n");
    nl_dvector_print(wi, n, NULL);
    printf("\nСобственные векторы для сбалансированной матрицы\n");
    nl_dmatrix_print(Q, n, n, NULL);

    eig_balance_inverse(Q, n, low, high, scal);  // Если нужны собственные векторы
    printf("\nСобственные векторы для исходной матрицы\n");
    nl_dmatrix_print(Q, n, n, NULL);

    eig_norm_Inf(Q, n, wi); // Если нужно нормировать собственные векторы
    printf("\nНормированные (чебышева норма) собственные векторы \n");
    nl_dmatrix_print(Q, n, n, NULL);


    nl_dmatrix_free(A, n);
    nl_dmatrix_free(Q, n);
    nl_dvector_free(scal);
    nl_dvector_free(wr);
    nl_dvector_free(wi);
    nl_xvector_free(perm);
    
    return 0;
}
