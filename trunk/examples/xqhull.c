#include <time.h>
#include "nl.h"

/*
  Алгоритм Quick hull построения выпуклой оболочки точек на плоскости
*/

int main()
{

  double *x, *y;
  size_t n, *k, *stack, nk, j;

  n = 10;

  x = nl_dvector_create(n);
  y = nl_dvector_create(n);
  k = nl_xvector_create(n);
  stack = nl_xvector_create(4*ceil(log(n)/log(2)) + 1);

  x[0] = 0.9501; y[0] = 0.6154;
  x[1] = 0.2311; y[1] = 0.7919;
  x[2] = 0.6068; y[2] = 0.9218;
  x[3] = 0.4860; y[3] = 0.7382;
  x[4] = 0.8913; y[4] = 0.1763;
  x[5] = 0.7621; y[5] = 0.4057;
  x[6] = 0.4565; y[6] = 0.9355;
  x[7] = 0.0185; y[7] = 0.9169;
  x[8] = 0.8214; y[8] = 0.4103;
  x[9] = 0.4447; y[9] = 0.8936;

  conv_qhull(x, y, n, k, &nk, stack);

  printf("Заданные точки:\n");
  for (j = 0; j < n; j++) 
    printf("  %f  %f \n", x[j], y[j]);

  printf("\nВершины выпуклой оболочки:\n");
  for (j = 0; j < nk; j++) 
    printf("  %f  %f \n", x[k[j]], y[k[j]]);

  nl_dvector_free(x);
  nl_dvector_free(y);
  nl_xvector_free(k);
  nl_xvector_free(stack);

  return 0;
}
