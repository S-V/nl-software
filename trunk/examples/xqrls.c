/*    
  Пример использования функций из модуля @qr.h@
  Решение задачи наименьших квадратов с помощью 
  QR-разложения
  $A =   \left(  \begin{array}{rrr}     1 & 5 & 1 \\       2 & 6 & 10 \\      3 & 7 & 11 \\      4 & 8 & 12 \\   \end{array}  \right)  $, $\quad b =    \left(   \begin{array}{r}      1 \\        1 \\       1 \\       2 \\    \end{array}   \right)  $
*/

#include "nl.h"

int main()
{
  double *A, *b, *t, *r, *work;
  size_t m = 4;
  size_t n = 3;

  A = nl_dmatrix_create(m, n);
  b = nl_dvector_create(m);
  t = nl_dvector_create(n);
  r = nl_dvector_create(m);
  work = nl_dvector_create(m);

  A[0] = 1; A[1]  = 5; A[2]  =  1;   b[0] = 1;
  A[3] = 2; A[4]  = 6; A[5]  = 10;   b[1] = 1;
  A[6] = 3; A[7]  = 7; A[8]  = 11;   b[2] = 1;
  A[9] = 4; A[10] = 8; A[11] = 12;   b[3] = 2;

  printf("Матрица A:\n");
  nl_dmatrix_print(A, m, n, NULL);

  qr_decomp(A, m, n, t);

  printf("\nQR-разложение:\n");
  nl_dmatrix_print(A, m, n, NULL);

  printf("\nВектор t:\n");
  nl_dvector_print(t, n, NULL);

  printf("\nВектор b:\n");
  nl_dvector_print(b, m, NULL);

  qr_least_squares(A, m, n, t, b, r, work);

  printf("\nПсевдорешение системы Ax = b:\n");
  nl_dvector_print(b, n, NULL);

  printf("\nНевязки:\n");
  nl_dvector_print(r, m, NULL);

  nl_dmatrix_free(A);
  nl_dvector_free(b);
  nl_dvector_free(t);
  nl_dvector_free(r);
  nl_dvector_free(work);

  return 0;
}
