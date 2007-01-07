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
  double *C, *L, *b, *x;
  int sgn;

  C = nl_dmatrix_create(n, m);
  L = nl_dmatrix_create(n, m1);
  b = nl_dvector_create(n);
  x = nl_dvector_create(n);
  p = nl_xvector_create(n);

  C[0]  = 0;  C[1]  = 0;  C[2]  = 4;  C[3]  = 1;
  C[4]  = 0;  C[5]  = 2;  C[6]  = 4;  C[7]  = 1;
  C[8]  = 1;  C[9]  = 2;  C[10] = 4;  C[11] = 1;
  C[12] = 1;  C[13] = 2;  C[14] = 4;  C[15] = 1;
  C[16] = 1;  C[17] = 2;  C[18] = 4;  C[19] = 1;
  C[20] = 1;  C[21] = 2;  C[22] = 4;  C[23] = 1;
  C[24] = 1;  C[25] = 2;  C[26] = 4;  C[27] = 0;

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

  nl_dmatrix_free(C);
  nl_dmatrix_free(L);
  nl_dvector_free(b);
  nl_xvector_free(p);
  
  return 0;
}
