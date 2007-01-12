#include <stdio.h>
#include <math.h>
#include "nl.h"

double fun(double *x)
{

  double tmp1, tmp2;

  tmp1 = x[1] - x[0]*x[0];
  tmp2 = 1 - x[0];

  return tmp1*tmp1 + tmp2*tmp2;
}

int main(void)
{
  size_t n;
  double *x0, f, *work;
  double f0;
  int rc, nfunevals, niter;

  n = 2;

  x0 = nl_dvector_create(n);
  work = nl_dvector_create(2*n);

  x0[0] = -6;
  x0[1] = 10;

  rc = optim_hooke_jeeves(fun, n, x0, 0.5, 1e-6, 1e-6, 500, 200, &f, &nfunevals, &niter, work);

  if (rc < 0)
  {
    printf("\nЧисло итераций или количество вычисленных значений функции\n");
    printf("превысило максимально допустимое!\n" );
  }
  else
    printf("\nУспешное завершение\n");

  printf("\nКоличество вычисленных значений функции: %d\n", nfunevals);
  printf("Число итераций: %d\n",niter);

  printf("\nВычисленная точка минимума:\n");
  nl_dvector_print(x0, n, "  %12.6e");

  printf("\nЗначения функции: %e\n", f);
  
  nl_dvector_free(x0);
  nl_dvector_free(work);

  return 0;
}


