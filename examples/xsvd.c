/*    
  Пример использования функций из модуля @svd.h@
  Решение задачи наименьших квадратов $Ax=b$
  с помощью сингулярного разложения
  $A =   \left(  \begin{array}{rrr}     1 & 5 & 9 \\       2 & 6 & 10 \\      3 & 7 & 11 \\      4 & 8 & 12 \\   \end{array}  \right)  $, $\quad b =    \left(   \begin{array}{r}      1 \\        1 \\       1 \\       2 \\    \end{array}   \right)  $
*/

#include "nl.h"

int main()
{
  double **A, **U, **V, *w, *b, *x;
  size_t ierr;
  size_t m = 4;
  size_t n = 3;

  A = nl_dmatrix_create(m, n);
  U = nl_dmatrix_create(m, n);
  V = nl_dmatrix_create(n, n);
  w = nl_dvector_create(n);
  b = nl_dvector_create(m);
  x = nl_dvector_create(n);

  A[0][0] = 1; A[0][1] = 5; A[0][2] =  9;   b[0] = 1;
  A[1][0] = 2; A[1][1] = 6; A[1][2] = 10;   b[1] = 1;
  A[2][0] = 3; A[2][1] = 7; A[2][2] = 11;   b[2] = 1;
  A[3][0] = 4; A[3][1] = 8; A[3][2] = 12;   b[3] = 2;

  svd_decomp(A, m, n, w, 1, U, 1, V, &ierr);
  if (ierr)
  {
    printf("Сингулярные числа не были найдены ");
    printf("за 30 итераций\n");
    return 1;
  }

  printf("Матрица A:\n");
  nl_dmatrix_print(A, m, n, NULL);

  printf("\nМатрица U:\n");
  nl_dmatrix_print(U, m, n, NULL);

  printf("\nМатрица V:\n");
  nl_dmatrix_print(V, n, n, NULL);

  printf("\nСингулярные числа w:\n");
  nl_dvector_print(w, n, NULL);

  svd_correct(w, n, 1e-16);
  svd_least_squares(U, w, V, m, n, b, x);

  printf("\nПравая часть системы:\n");
  nl_dvector_print(b, m, NULL);

  printf("\nНормальное псевдорешение:\n");
  nl_dvector_print(x, n, NULL);

  nl_dmatrix_free(A, m);
  nl_dmatrix_free(U, m);
  nl_dmatrix_free(V, n);
  nl_dvector_free(w);
  nl_dvector_free(b);
  nl_dvector_free(x);
  
  return 0;
}
