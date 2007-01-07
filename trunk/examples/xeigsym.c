/*
  Симметричная проблема собственных значений.
  $QR$-алгоритм.
*/

#include "nl.h"

int main()
{
    size_t n = 5;
    size_t rc;
    double *A, *d, *e;

    A = nl_dmatrix_create(n, n);
    d = nl_dvector_create(n);
    e = nl_dvector_create(n);
    
    A[0]  = 11;
    A[5]  = 10; A[6]  = 14; 
    A[10] = 18; A[11] = 41; A[12] = 17; 
    A[15] = 11; A[16] = 12; A[17] = 23; A[18] = 11; 
    A[20] = 21; A[21] = 13; A[22] = 23; A[23] = 17; A[24] = 17; 
    
    printf("\nМатрица A (нижняя треугольная часть):\n");
    nl_dmatrix_print(A, n, n, NULL);

    eig_tridiag_reduction(A, n, 1, d, e);

    printf("\nТрехдиагональный вид\n");
    printf("Диагональ:\n");
    nl_dvector_print(d, n, NULL);

    printf("\nПоддиагональ = наддиагональ:\n");
    nl_dvector_print(e, n, NULL);

    printf("\nМатрица перехода Q:\n");
    nl_dmatrix_print(A, n, n, NULL);

    eig_tridiag(d, e, n, 1, A, &rc);

    printf("\nСобственные числа:\n");
    nl_dvector_print(d, n, NULL);

    printf("\nСобственные векторы:\n");
    nl_dmatrix_print(A, n, n, NULL);

    nl_dmatrix_free(A);
    nl_dvector_free(d);
    nl_dvector_free(e);
    
    return 0;
}
