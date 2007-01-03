#include <stdio.h>
#include <math.h>
#include "nl.h"

/*
  Пример использования функций из модуля @interp.h@.

  Интегрирование таблично заданной функции с помощью
  эрмитова кубического интерполянта и кубического сплайна
*/
int main()
{
   size_t j, n;
   double *x, *y, *b, *c, *d, xx, *work;

   n = 5;

   x = nl_dvector_create(n);
   y = nl_dvector_create(n);
   b = nl_dvector_create(n - 1);
   c = nl_dvector_create(n - 1);
   d = nl_dvector_create(n);
   work = nl_dvector_create(2*n);

   x[0] = 0.00;
   x[1] = 0.25;
   x[2] = 0.50;
   x[3] = 0.75;
   x[4] = 1.00;

   for (j = 0; j < n; j++)
      y[j] = sqrt(x[j]);

   interp_pchip(x, y, n, b, c, d);
   printf("\nИнтеграл на основе interp_pchip  = %f\n", 
     interp_quad(x, y, b, c, d, n, 0, 1));

   interp_spline(x, y, n, b, c, d, work);
   printf("\nИнтеграл на основе interp_spline = %f\n", 
     interp_quad(x, y, b, c, d, n, 0, 1));

   nl_dvector_free(x);
   nl_dvector_free(y);
   nl_dvector_free(b);
   nl_dvector_free(c);
   nl_dvector_free(d);
   nl_dvector_free(work);

   return 0;
}