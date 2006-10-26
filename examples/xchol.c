/*    
  Пример использования функций из модуля @chol.h@
  Решение симметричной положительно определенной 
  системы $Ax=b$ с помощью разложения Холецкого,
  $A =   \left(  \begin{array}{rrr}     0.1 & 0.2 & 0.3 \\       0.2 & 0.7 & 0.9 \\      0.3 & 0.9 & 1.3 \\   \end{array}  \right)  $, $\quad b =    \left(   \begin{array}{r}      0.6 \\        1.8 \\       2.5 \\    \end{array}   \right)  $.
*/

#include "nl.h"

int main()
{
  double **A, *b;
  size_t n = 3;

  A = nl_dmatrix_create(n, n);
  b = nl_dvector_create(n);

  A[0][0] = .1; A[0][1] = .2; A[0][2] =  .3;   b[0] =  .6;
  A[1][0] = .2; A[1][1] = .7; A[1][2] =  .9;   b[1] = 1.8;
  A[2][0] = .3; A[2][1] = .9; A[2][2] = 1.3;   b[2] = 2.5;

  printf("Матрица A:\n");
  nl_dmatrix_print(A, n, n, NULL);

  chol_decomp(A, n); 

  printf("\nРазложение Холецкого:\n");
  printf("\n(в нижней треугольной части):\n");
  nl_dmatrix_print(A, n, n, NULL);

  printf("\nВектор b:\n");
  nl_dvector_print(b, n, NULL);

  chol_solve(A, n, b);

  printf("\nРешение системы Ax = b:\n");
  nl_dvector_print(b, n, NULL);

  nl_dmatrix_free(A, n);
  nl_dvector_free(b);

  return 0;
}
