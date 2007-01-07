/*    
  Пример использования функции @gall_laplace@
*/

#include "nl.h"

int main()
{
  size_t *IA, *JA, m, n, N, nz;
  double *AN, *AD, *A;

  m = n = 3;
  N = m*n;

  printf("5-точечная аппроксимация уравнения Пуассона\n");
  printf("Решетка 3x3\n");
  printf("\nРазреженное полное представление:\n");

  gall_laplace_size(m, n, &N, &nz);
  sp_create(N, nz, &IA, &JA, &AN);
  gall_laplace(m, n, IA, JA, AN);
  sp_print(IA, JA, AN, N, N, NULL);

  A = nl_dmatrix_create(N, N);
  sp_full(IA, JA, AN, N, N, A);

  printf("\nПлотное представление:\n");
  nl_dmatrix_print(A, N, N, NULL);

  sp_free(IA, JA, AN);

  printf("\nРазреженное неполное представление\n");
  printf("(верхнетреугольная часть):\n");

  gall_laplace_sym_size(m, n, &N, &nz);
  sp_create_sym(N, nz, &IA, &JA, &AN, &AD);
  gall_laplace_sym(m, n, IA, JA, AN, AD);
  sp_print(IA, JA, AN, N, N, NULL);

  A = nl_dmatrix_create(N, N);
  sp_full(IA, JA, AN, N, N, A);

  printf("\nПлотное представление:\n");
  nl_dmatrix_print(A, N, N, NULL);

  printf("\nДиагональные элементы:\n");
  nl_dvector_print(AD, N, NULL);

  sp_free_sym(IA, JA, AN, AD);
  nl_dmatrix_free(A);

  return 0;
}
