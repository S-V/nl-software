#include "nl.h"

/*
  Пример использования функций из модуля @mda.h@
  Метод минимальной степени
*/

int main()
{

  double *A, *SN, *SD, *UN, *UD, *x, *b;  
  size_t *I, *J, *IA, *JA, *IS, *JS, *D, *P, *IP, *M, *L, *IU, *JU, *xwork;
  size_t n, nz, k, unz;  

  n = 5;
  nz = n - 1;

  A = nl_dvector_create(nz);
  I = nl_xvector_create(nz);
  J = nl_xvector_create(nz);
  xwork = nl_xvector_create(2*n);

  /* Генерируем верхнетреугольную часть матрицы A
     - стреловидный портрет */

  for (k = 0; k < nz; k++)
    A[k] = 1.;
  for (k = 0; k < nz; k++)
  {
    I[k] = 0;
    J[k] = k + 1;
  }

  sp_create_sym(n, nz, &IS, &JS, &SN, &SD);
  sp_convert(nz, A, I, J, n, IS, JS, SN);
  sp_order(IS, JS, SN, n);

  SD[0] = n;
  for(k = 1; k < n; k++) 
    SD[k] = 1.;

  printf("\nНаддиагональные элементы до разложения Холецкого:\n");
  sp_print_list_sym(IS, JS, SN, SD, n, NULL, NULL);

  /* Конвертируем представление в полное и вызываем алгоритм mda */

  mda_create(n, nz, &IA, &JA, &D, &P, &IP, &M, &L);
  mda_convert(n, IS, JS, IA, JA, D, P, IP);
  unz = mda_order(n, IA, JA, M, L, D, P, IP);

  printf("\nПерестановка, найденная методом mda:\n");
  nl_xvector_print(P, n, NULL);
  printf("\nОбратная перестановка:\n");
  nl_xvector_print(IP, n, NULL);

  /* Разложение Холецкого */

  sp_create_sym(n, unz, &IU, &JU, &UN, &UD);
  mda_chol_symb(n, unz, IS, JS, IU, JU, M, L, D, IP);

  // sp_order_profile(IU, JU, n); //?
  mda_chol_num(n, IS, JS, IU, JU, IP, SN, SD, UN, UD, xwork);

  printf("\nНаддиагональные элементы после числ. разложения Холецкого:\n");
  sp_print_list_sym(IU, JU, UN, UD, n, NULL, NULL);

  /* Генерируем решение и правую часть системы */

  x = nl_dvector_create(n);
  for(k = 0; k < n; k++) 
    x[k] = (k + 1)*11;

  b = nl_dvector_create(n);
  sp_mult_col_sym(IS, JS, SN, SD, x, n, b);
  printf("\nПравая часть системы:\n");
  nl_dvector_print(b, n, NULL);

  /* Решаем систему */

  mda_chol_solve(n, IU, JU, P, UN, UD, b, x);

  printf("\nРешение:\n");
  nl_dvector_print(x, n, NULL);


  /* Освобождаем память */

  nl_dvector_free(A);
  nl_xvector_free(I);
  nl_xvector_free(J);
  nl_xvector_free(xwork);
  sp_free_sym(IS, JS, SN, SD);
  mda_free(IA, JA, D, P, IP, M, L);
  sp_free_sym(IU, JU, UN, UD);
  nl_dvector_free(x);
  nl_dvector_free(b);

  return 0;

}
