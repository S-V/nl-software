#include <stdio.h>
#include <math.h>
#include "nl.h"

/*
   Пример использования функций из модуля @ode.h@ 
   Решение задачи Коши
   методом Рунге-Кутты-Фельберга.
   Решение периодично. Период равен $8$.   
*/

void fun(double t, double* y, double* ydot)
{
   double alpha, r;

   alpha = 3.14159/4;
   r = y[0]*y[0] + y[2]*y[2];
   r = r*sqrt(r)/alpha/alpha;

   ydot[0] = y[1];
   ydot[1] = -y[0]/r;
   ydot[2] = y[3];
   ydot[3] = -y[2]/r;
}

int main()
{
   size_t n;
   double t, t0, tout, h;
   double *y0, *y, *ydot, *work;
   int rc, nfunevals, totalnfunevals;

   n = 4;

   t0 = 0;
   tout = 12;
   h = 12;

   y0 = nl_dvector_create(n);
   y = nl_dvector_create(n);
   ydot = nl_dvector_create(n);
   work = nl_dvector_create(5*n);

   y0[0] = 1 - 0.25;
   y0[1] = 0;
   y0[2] = 0;
   y0[3] = 3.14159/4 * sqrt((1 + 0.25)/(1 - 0.25));

   printf("   t         y[0]         y[2]   \n");
   printf(" ------------------------------- \n");
   printf(" %5.2f    %9.6f    %9.6f \n", t0, y0[0], y0[2]);

   totalnfunevals = 0;

   for (t = 0; t < tout; t += h)
   {
      ode_rkf45(&fun, n, t, t + h, y0, 1e-3, 1e-6, 600, &rc, y, ydot, &nfunevals, work);

      totalnfunevals += nfunevals;
      
      switch (rc)
      {
          case 0:
              printf("Число обращений к функции, вычисляющей правую часть,\n");
              printf("превысило допустимое максимальное значение\n");
              return 0;
          case -1:
              printf("Размер шага меньше допустимого минимального значения\n");
              return -1;
      }
      printf(" %5.2f    %9.6f    %9.6f \n", t + h, y[0], y[2]);
   }

   printf("Общее количество обращений к функции,\n");
   printf("вычисляющей правую часть, = %d", totalnfunevals);

   nl_dvector_free(y0);
   nl_dvector_free(y);
   nl_dvector_free(ydot);
   nl_dvector_free(work);

   return 1;

}