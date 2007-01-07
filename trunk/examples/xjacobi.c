/*
  Симметричная проблема собственных значений.
  Метод Якоби
*/

#include "nl.h"

int main()
{
    size_t n = 5;
    double *A, *V, *d, *work;
    int nrot, rc;

    A = nl_dmatrix_create(n, n);
    V = nl_dmatrix_create(n, n);
    d = nl_dvector_create(n);
    work = nl_dvector_create(2*n);
  
    A[0]  = 11; A[1]  = 10; A[2]  = 18; A[3]  = 11; A[4]  = 21; 
    A[5]  = 10; A[6]  = 14; A[7]  = 41; A[8]  = 12; A[9]  = 13; 
    A[10] = 18; A[11] = 41; A[12] = 17; A[13] = 23; A[14] = 23; 
    A[15] = 11; A[16] = 12; A[17] = 23; A[18] = 11; A[19] = 17; 
    A[20] = 21; A[21] = 13; A[22] = 23; A[23] = 17; A[24] = 17; 
  
    printf("\nМатрица A (верхняя треугольная часть):\n");
    nl_dmatrix_print(A, n, n, NULL);

    eig_jacobi(A, n, d, 1, V, &nrot, &rc, work);

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

    nl_dmatrix_free(A);
    nl_dmatrix_free(V);
    nl_dvector_free(d);
    nl_dvector_free(work);
    
    return 0;
}
