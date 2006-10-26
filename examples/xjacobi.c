/*
  Симметричная проблема собственных значений.
  Метод Якоби
*/

#include "nl.h"

int main()
{
    size_t n = 5;
    double **A, **V, *d;
    int nrot, rc;

    A = nl_dmatrix_create(n, n);
    V = nl_dmatrix_create(n, n);
    d = nl_dvector_create(n);
  
    A[0][0] = 11; A[0][1] = 10; A[0][2] = 18; A[0][3] = 11; A[0][4] = 21; 
    A[1][0] = 10; A[1][1] = 14; A[1][2] = 41; A[1][3] = 12; A[1][4] = 13; 
    A[2][0] = 18; A[2][1] = 41; A[2][2] = 17; A[2][3] = 23; A[2][4] = 23; 
    A[3][0] = 11; A[3][1] = 12; A[3][2] = 23; A[3][3] = 11; A[3][4] = 17; 
    A[4][0] = 21; A[4][1] = 13; A[4][2] = 23; A[4][3] = 17; A[4][4] = 17; 
  
    printf("\nМатрица A (верхняя треугольная часть):\n");
    nl_dmatrix_print(A, n, n, NULL);

    eig_jacobi(A, n, d, 1, V, &nrot, &rc);

    if (rc)
    {
      printf("\nЧисло итераций превысило 50\n");
    }
    else
    {
      printf("\nСобственные числа:\n");
      nl_dvector_print(d, n, NULL);

      printf("\nСобственные векторы:\n");
      nl_dmatrix_print(V, n, n, NULL);

      printf("\nКоличество выполненных вращений = %d\n", nrot);
    }

    nl_dmatrix_free(A, n);
    nl_dmatrix_free(V, n);
    nl_dvector_free(d);
    
    return 0;
}
