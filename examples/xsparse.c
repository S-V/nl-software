/*
  Пример использования функций из модуля @sparse.h@
  Конвертация разреженной матрицы в представление 
  в виде списка ненулевых элементов и обратно
  $A=  \left(  \begin{array}{rrrrr}    1 & -1 &    &    &    \\   -5 &  2 & -2 &    &    \\      & -6 &  3 & -3 &    \\      &    & -7 &  4 & -4 \\      &    &    & -8 &  5 \\  \end{array}  \right).  $
*/
#include <string.h> 

#include "nl.h"

int main()
{
  size_t n = 5;
  size_t nz = 13;
  double A[] = {1, 2, 3, 4, 5,-1,-2,-3,-4,-5,-6,-7,-8};
  size_t I[] = {0, 1, 2, 3, 4, 0, 1, 2, 3, 1, 2, 3, 4};
  size_t J[] = {0, 1, 2, 3, 4, 1, 2, 3, 4, 0, 1, 2, 3};
  size_t *IA, *JA;
  double *AN, *AF;

  sp_create(n, nz, &IA, &JA, &AN);
  sp_convert(nz, A, I, J, n, IA, JA, AN);

  printf("\nПредставление в виде ");
  printf("списка ненулевых элементов:\n");
  sp_print_list(IA, JA, AN, n, n, NULL, NULL);

  printf("Разреженное представление:\n");
  sp_print(IA, JA, AN, n, n, NULL);

  printf("\nПлотное представление:\n");
  AF = nl_dmatrix_create(n, n);
  sp_full(IA, JA, AN, n, n, AF);
  nl_dmatrix_print(AF, n, n, NULL);

  sp_free(IA, JA, AN);
  nl_dmatrix_free(AF);

  return 0;
}
