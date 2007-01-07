#include <stdio.h>
#include <math.h>
#include "nl.h"

/* Использование функций из модуля @interp.h@
   
   Построение и вычисление значений кубического эрмитова интерполянта 
   и кубического сплайна
*/

int main()
{
   size_t n;
   double *x, *y, *b, *c, *d, xx, *work;

   n = 5;

   x = nl_dvector_create(n);
   y = nl_dvector_create(n);
   b = nl_dvector_create(n - 1);
   c = nl_dvector_create(n - 1);
   d = nl_dvector_create(n);
   work = nl_dvector_create(3*n);

   x[0] = -2; y[0] = 0;
   x[1] = -1; y[1] = 0;
   x[2] =  0; y[2] = 1;
   x[3] =  1; y[3] = 0;
   x[4] =  2; y[4] = 0;

   printf("\nЭрмитов интерполянт\n\n");

   interp_pchip(x, y, n, b, c, d);

   printf("     x       y    \n");
   printf(" ---------------- \n");

   for (xx = -2; xx <= 2; xx += .25)
   {
      printf("  %5.2f  %7.4f\n", xx, interp_eval(x, y, b, c, d, n, xx));
   }

   printf("\nКубический сплайн\n\n");

   interp_spline(x, y, n, b, c, d, work);

   printf("     x       y    \n");
   printf(" ---------------- \n");

   for (xx = -2; xx <= 2; xx += .25)
   {
      printf("  %5.2f  %7.4f\n", xx, interp_eval(x, y, b, c, d, n, xx));
   }

   nl_dvector_free(x);
   nl_dvector_free(y);
   nl_dvector_free(b);
   nl_dvector_free(c);
   nl_dvector_free(d);
   nl_dvector_free(work);

   return 0;
}

