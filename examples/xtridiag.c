/*
  Пример использования функций из модуля @band.h@
  Решение трехдиагональной системы
  %%$  \left\{    \begin{narrowarray}{rcrcrcr}       3x\sb 1 & + & 6x\sb 2 &   &      & = &  9,\\        x\sb 1 & + & 4x\sb 2 & + & 7x\sb 3 & = & 12,\\            &   & 2x\sb 2 & + & 5x\sb 3 & = &  7.\\    \end{narrowarray}  \right.  $%%
*/

#include "nl.h"

int main()
{
  size_t n = 3;
  double a[] = {1, 2};
  double d[] = {3, 4, 5};
  double c[] = {6, 7};
  double b[] = {9, 12, 7};
  double *x;

  x = nl_dvector_create(n);

  band_tridiag(a, d, c, b, x, n);
  nl_dvector_print(x, n, NULL);

  nl_dvector_free(x);

  return 0;
}
