/*
  Пример использования функций из модуля @sparse.h@
  Решение разреженной положительно определенной системы 
  $Ax=b$ с помощью $LDL\transpose$-разложения (Холецкого).
  $  A=  \left(  \begin{array}{rrrrr}     3 & 1 &   &   &   \\     1 & 3 & 1 &   &   \\       & 1 & 3 & 1 &   \\       &   & 1 & 3 & 1 \\       &   &   & 1 & 3 \\  \end{array}  \right)  ,\quad  b=  \left(  \begin{array}{r}     4  \\     5  \\     5  \\     5  \\     4  \\  \end{array}  \right)  $.
*/
#include <string.h> 

#include "nl.h"

int main()
{
  double A[] = {1, 1, 1, 1};
  size_t I[] = {0, 1, 2, 3};
  size_t J[] = {4, 4, 4, 4};
  double D[] = {1, 1, 1, 1, 5};
  double DINV[5];
  size_t *IA, *JA, *IU, *JU; 
  double *AN, *UN;
  double b[] = {2, 2, 2, 2, 9};
  double x[5];
  size_t *xwork;

  sp_create(5, 4, &IA, &JA, &AN);
  sp_convert(4, A, I, J, 5, IA, JA, AN);

  xwork = nl_xvector_create(10);

  printf("Разложение A = U'*D*U\n");
  printf("\nВерхнетреугольная часть матрицы A:\n");
  sp_print_list_sym(IA, JA, AN, D, 5, 0, 0);

  sp_create(5, 4, &IU, &JU, &UN);
  sp_chol_symb(IA, JA, 5, IU, JU, 5, xwork);
  sp_chol_num(IA, JA, AN, D, IU, JU, 5, UN, DINV, xwork);

  printf("\nВерхнетреугольная часть матрицы U:\n");
  sp_print_list(IU, JU, UN, 5, 5, 0, 0);

  printf("\nЭлементы, обратные диагональным элементам ");
  printf("матрицы D:\n");
  nl_dvector_print(DINV, 5, 0);

  sp_chol_solve(IU, JU, UN, DINV, b, 5, x);

  printf("\nВектор b:\n");
  nl_dvector_print(b, 5, 0);

  printf("\nРешение системы Ax = b:\n");
  nl_dvector_print(x, 5, 0);

  sp_free(IA, JA, AN);
  sp_free(IU, JU, UN);

  nl_xvector_free(xwork);

  return 0;
}
