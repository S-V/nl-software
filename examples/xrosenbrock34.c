#include <stdio.h>
#include <math.h>
#include "nl.h"

void jac(double t, double* y, double* fprimet, double** fprimey)
{
   fprimet[0] = 0;
   fprimet[1] = 0;
   fprimet[2] = 0;

   fprimey[0][0] = -0.013 - 1000*y[2];
   fprimey[0][1] = 0;
   fprimey[0][2] = -1000*y[0];
   fprimey[1][0] = 0;
   fprimey[1][1] = -2500*y[2];
   fprimey[1][2] = -2500*y[1];
   fprimey[2][0] = -0.013 - 1000*y[2];
   fprimey[2][1] = -2500*y[2];
   fprimey[2][2] = -1000*y[0] - 2500*y[1];
}

void fun(double t, double* y, double* ydot)
{
   ydot[0] = -0.013*y[0] - 1000*y[0]*y[2];
   ydot[1] = -2500*y[1]*y[2];
   ydot[2] = -0.013*y[0] - 1000*y[0]*y[2] - 2500*y[1]*y[2];
}


int main()
{
   size_t n, *p;
   double t, t0, tout, h;
   double *y0, *ydot, *y, *fprimet, *work, **a, **fprimey;
   int rc, nfunevals, totalnfunevals;

   n = 3;

   y0 = nl_dvector_create(n);
   ydot = nl_dvector_create(n);
   y = nl_dvector_create(n);
   fprimet = nl_dvector_create(n);
   work = nl_dvector_create(5*n);
   p = nl_xvector_create(n);

   a = nl_dmatrix_create(n, n);
   fprimey = nl_dmatrix_create(n, n);

   t0 = 0;
   tout = 50;
   h = 1;

   y0[0] = 1;
   y0[1] = 1;
   y0[2] = 0;

   printf("   t         y[0]         y[1]         y[2]   \n");
   printf(" -------------------------------------------- \n");
   printf(" %5.2f    %9.6f    %9.6f    %9.6f \n", t0, y0[0], y0[1], y0[2]);

   totalnfunevals = 0;

   for (t = 0; t < tout; t += h)
   {
      ode_rosenbrock34(&fun, &jac, n, t, t + h, y0, 1e-3, 1e-6, 600, 
         &rc, y, ydot, fprimet, fprimey, &nfunevals, a, work, p);

      totalnfunevals += nfunevals;
      
      switch (rc)
      {
          case 0:
              printf("„исло обращений к функции, вычисл€ющей правую часть,\n");
              printf("превысило допустимое максимальное значение\n");
              return 0;
          case -1:
              printf("–азмер шага меньше допустимого минимального значени€\n");
              return -1;
      }
      printf(" %5.2f    %9.6f    %9.6f    %9.6f \n", t + h, y[0], y[1], y[2]);
   }

   printf("ќбщее количество обращений к функции,\n");
   printf("вычисл€ющей правую часть, = %d", totalnfunevals);

   nl_dvector_free(y0);
   nl_dvector_free(ydot);
   nl_dvector_free(y);
   nl_dvector_free(fprimet);
   nl_dvector_free(work);

   nl_dmatrix_free(a, n);
   nl_dmatrix_free(fprimey, n);

   nl_xvector_free(p);

   return 1;

}