/*
  Пример использования функций из модуля @band.h@
  Решение ленточной системы
  $    \left\{    \begin{narrowarray}{rcrcrcrcrcrcrccr}       4 x\sb 1 & + &  x\sb 2 &   &         &   &         &   &         &   &         &   &         & ~=~ & 5, \\       2 x\sb 1 & + & 4x\sb 2 & + &  x\sb 3 &   &         &   &         &   &         &   &         & ~=~ & 7, \\         x\sb 1 & + & 2x\sb 2 & + & 4x\sb 3 & + &  x\sb 4 &   &         &   &         &   &         & ~=~ & 8, \\                &   &  x\sb 2 & + & 2x\sb 3 & + & 4x\sb 4 & + &  x\sb 5 &   &         &   &         & ~=~ & 8, \\                &   &         &   &  x\sb 3 & + & 2x\sb 4 & + & 4x\sb 5 & + &  x\sb 6 &   &         & ~=~ & 8, \\                &   &         &   &         &   &  x\sb 4 & + & 2x\sb 5 & + & 4x\sb 6 & + &  x\sb 7 & ~=~ & 8, \\                &   &         &   &         &   &         &   &  x\sb 5 & + & 2x\sb 6 & + & 4x\sb 7 & ~=~ & 7. \\    \end{narrowarray}      \right.   $.
*/

#include "nl.h"

int main()
{
  size_t n = 7;
  size_t m1 = 2;
  size_t m2 = 1;
  size_t m = m1 + m2 + 1;

  size_t *p;
  double **C, **L, *b, *x;
  int sgn;

  C = nl_dmatrix_create(n, m);
  L = nl_dmatrix_create(n, m1);
  b = nl_dvector_create(n);
  x = nl_dvector_create(n);
  p = nl_xvector_create(n);

  C[0][0] = 0; C[0][1] = 0; C[0][2] = 4; C[0][3] = 1;
  C[1][0] = 0; C[1][1] = 2; C[1][2] = 4; C[1][3] = 1;
  C[2][0] = 1; C[2][1] = 2; C[2][2] = 4; C[2][3] = 1;
  C[3][0] = 1; C[3][1] = 2; C[3][2] = 4; C[3][3] = 1;
  C[4][0] = 1; C[4][1] = 2; C[4][2] = 4; C[4][3] = 1;
  C[5][0] = 1; C[5][1] = 2; C[5][2] = 4; C[5][3] = 1;
  C[6][0] = 1; C[6][1] = 2; C[6][2] = 4; C[6][3] = 0;

  x[0] = 1; x[1] = 1; x[2] = 1; x[3] = 1;
  x[4] = 1; x[5] = 1; x[6] = 1; 
 
  printf("\nМатрица A:\n");
  nl_dmatrix_print(C, n, m, NULL);

  band_mult_col(C, n, m1, m2, x, b);
   
  printf("\nВектор b:\n");
  nl_dvector_print(b, n, NULL);

  band_decomp(C, n, m1, m2, L, p, &sgn);
  band_solve(C, n, m1, m2, L, p, b);

  printf("\nРешение системы Ax = b:\n");
  nl_dvector_print(b, n, NULL);

  nl_dmatrix_free(C, n);
  nl_dmatrix_free(L, n);
  nl_dvector_free(b);
  nl_xvector_free(p);
  
  return 0;
}
