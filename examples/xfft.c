/*
  Пример использования функций из модуля @fft.h@
  Быстрое преобразование Фурье
  $x=(0, 1+i, 2+2i, 3+3i, 4+4i, 5+5i, 6+6i, 7+7i)$
*/

#include "nl.h"

int main()
{
  size_t n;

  double *re, *im;

  n = 8;

  re = nl_dvector_create(n);
  im = nl_dvector_create(n);

  re[0] = 0; im[0] = 0;
  re[1] = 1; im[1] = 1;
  re[2] = 2; im[2] = 2;
  re[3] = 3; im[3] = 3;
  re[4] = 4; im[4] = 4;
  re[5] = 5; im[5] = 5;
  re[6] = 6; im[6] = 6;
  re[7] = 7; im[7] = 7;

  printf("Исходный вектор:\n");
  nl_dvector_print(re, n, NULL);
  nl_dvector_print(im, n, NULL);

  fft_transform(re, im, n);

  printf("\nfft:\n");
  nl_dvector_print(re, n, NULL);
  nl_dvector_print(im, n, NULL);

  fft_inverse(re, im, n);

  printf("\nifft:\n");
  nl_dvector_print(re, n, NULL);
  nl_dvector_print(im, n, NULL);

  nl_dvector_free(re);
  nl_dvector_free(im);

  return 0;
}
