/*
  Пример использования функций из модуля @sparse.h@
  Метод бисопряженных градиентов
  $  A=  \left(  \begin{array}{rrrrr}    31 & 1 &   &   &   \\     1 & 1 & 1 &   &   \\       & 1 & 1 & 1 &   \\       &   & 1 & 1 & 1 \\       &   &   & 1 & 1 \\  \end{array}  \right)  ,\quad  b=  \left(  \begin{array}{r}     32 \\     3  \\     3  \\     3  \\     2  \\  \end{array}  \right)  $
*/
#include <stdlib.h>
#include <string.h> 

#include "nl.h"

int main()
{
  size_t n = 5;
  size_t nz = 13;
  double A[] = {31, 24, 31, 41, 31, 4, 4, 31, 1, 21, 1, 1, 1};
  size_t I[] = {0, 1, 2, 3, 4, 0, 1, 2, 3, 1, 2, 3, 4};
  size_t J[] = {0, 1, 2, 3, 4, 1, 2, 3, 4, 0, 1, 2, 3};
  size_t *IA, *JA;
  double *AN;
  double *b;
  double x[] = {0, 0, 0, 0, 0};
  int it;

  double xx[] = {1, 1, 1, 1, 1};

  b = nl_dvector_create(n);

  sp_create(n, nz, &IA, &JA, &AN);
  sp_convert(nz, A, I, J, n, IA, JA, AN);

  sp_mult_col(IA, JA, AN, xx, n, b);

  printf("Метод бисопряженных градиентов\n");
  printf("Матрица A:\n");
  sp_print_list(IA, JA, AN, n, n, 0, 0);

  printf("\nВектор b:\n");
  nl_dvector_print(b, n, 0);

  it = sp_biconj(IA, JA, AN, b, n, 1e-3, 20, x);
  
  printf("\nРешение системы Ax = b:\n");
  nl_dvector_print(x, n, 0);

  printf("\nЧисло итераций = %i \n", it);

  sp_free(IA, JA, AN);
  nl_dvector_free(b);

  return 0;
}
