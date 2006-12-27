#include <math.h>
#include "nl.h"

/* 
  Один шаг метода Рунге-Кутты-Фельберга.
*/

void ode_rkf_step (
  void (*f)(double, double*, double*), 
  size_t n, 
  double t0, 
  double tout, 
  double* y0, 
  double* ydot, 
  double* y, 
  double* err, 
  double* work)
{
   double h, err_j;
   double *k1, *k2, *k3, *k4, *k5, *k6;
   size_t j;

   k1 = ydot;
   k2 = work;
   k3 = work + n;
   k4 = work + 1*n;
   k5 = work + 2*n;
   k6 = work + 3*n;
                    
   h = tout - t0;
       
   k1 = ydot;
   for (j = 0; j < n; j++)
      y[j] = y0[j] + h*k1[j]/4;

   (*f)(t0 + h/4, y, k2);
   for (j = 0; j < n; j++)
      y[j] = y0[j] + h*(3*k1[j] + 9*k2[j])/32;

   (*f)(t0 + 3*h/8, y, k3);     
   for (j = 0; j < n; j++)
      y[j] = y0[j] + h*(1932*k1[j] - 7200*k2[j] + 7296*k3[j])/2197;

   (*f)(t0 + 12*h/13, y, k4);
   for (j = 0; j < n; j++)
      y[j] = y0[j] + h*(439*k1[j]/216 - 8*k2[j] + 3680*k3[j]/513 - 845*k4[j]/4104);

   (*f)(tout, y, k5);
   for (j = 0; j < n; j++)
      y[j] = y0[j] + h*(-8*k1[j]/27 + 2*k2[j] - 3544*k3[j]/2565 + 1859*k4[j]/4104 - 11*k5[j]/40);

   (*f)(t0 + h/2, y, k6);
   for (j = 0; j < n; j++)
      y[j] = y0[j] + h*(16*k1[j]/135 + 6656*k3[j]/12825 + 28561*k4[j]/56430 - 9*k5[j]/50 + 2*k6[j]/55);


   *err = 0;

   for (j = 0; j < n; j++)
   {
      err_j = fabs(h*(k1[j]/360 - 128*k3[j]/4275 - 2197*k4[j]/75240 + k5[j]/50 + 2*k6[j]/55));
      if (err_j > *err)
        *err = err_j;
   }

}



void ode_rkf(
  void (*f)(double, double*, double*), 
  size_t n, 
  double t0, 
  double tout, 
  double *y0, 
  double reltol, 
  double abstol, 
  int maxfun, 
  int *rc, 
  double *y, 
  double *ydot,
  double *work)
{

   double eps, epsp1, sqrt_eps, err, t, h, hmin, delta;
   int fun_evals;

   eps = 1.0;

   do
   {
     eps = eps/2;
     epsp1 = 1.0 + eps;
   }
   while (epsp1 > 1.0);

   sqrt_eps = sqrt(eps);

   if (fabs(t0) > fabs(tout))
      hmin = sqrt_eps*fabs(t0);
   else 
      hmin = sqrt_eps*fabs(tout); 

   t = t0;
   h = tout - t0;

   (*f)(t, y0, ydot);
   fun_evals = 1;

   /* Основной цикл */

   while (tout > t0 && t < tout || tout < t0 && t > tout) 
   {
       if (tout > t0 && t + h > tout || tout < 0 && t + h < tout)
          h = tout - t;

       ode_rkf_step(f, n, t, t + h, y0, ydot, y, &err, work);
       fun_evals += 6;

       if (fun_evals > maxfun)
       {
          *rc = 0;
          return;
       }

       if(err <= abstol || err <= reltol*nl_dvector_norm_inf(y, n)) 
       {
          t += h;
          nl_dvector_copy(y0, y, n);

          (*f)(t, y0, ydot);
       }
       else
       {
          /* Выбор длины нового шага */

	  if (err == 0)
             delta = 4;
          else 
          {
             if (reltol < eps)
                delta = 0.8*pow(eps/err, 0.2);
             else
                delta = 0.8*pow(reltol/err, 0.2);
             if (delta < 0.1)
               delta = 0.1;
             else if (delta > 4)
               delta = 4;
          }
          h *= delta;

          if (h < hmin)
          {
             *rc = -1;
             return;
          }
       }
   }

   *rc = 1;
}


void ode_rosenbrock_step(
  void (*f)(double, double*, double*), 
  size_t n, 
  double t0, 
  double tout, 
  double *y0, 
  double *ydot0,
  double *fprimet0, 
  double **fprimey0, 
  double *y, 
  double *err,
  double **a, 
  double *work,
  size_t *p)
{
   double h, err_j, *g1, *g2, *g3, *g4, *ydot;
   size_t i, j;
   int sgn;

   g1 = work;
   g2 = work + n;
   g3 = work + 2*n;
   g4 = work + 3*n;
   ydot = work + 4*n;

   h = tout - t0;

   for (i = 0; i < n; i++) 
   {
      for (j = 0; j < n; j++) 
            a[i][j] = -fprimey0[i][j];
      a[i][i] += 2/h;
   }

   lu_decomp(a, n, p, &sgn);

   for (j = 0; j < n; j++)
      g1[j] = ydot0[j] + h*fprimet0[j]/2;

   lu_solve(a, n, p, g1);

   for (j = 0; j < n; j++)
      y[j] = y0[j] + 2*g1[j];

   (*f)(t0 + h, y, ydot);

   for (j = 0; j < n; j++)
      g2[j] = ydot[j] - h*3*fprimet0[j]/2 - 8*g1[j]/h;

   lu_solve(a, n, p, g2);

   for (j = 0; j < n; j++)
      y[j] = y0[j] + 48*g1[j]/25 + 6*g2[j]/25;

   (*f)(t0 + 3*h/5, y, ydot);

   for (j = 0; j < n; j++)
      g3[j] = ydot[j] + h*121*fprimet0[j]/50 + (372*g1[j]/25 + 12*g2[j]/5)/h;

   lu_solve(a, n, p, g3);

   for (j = 0; j < n; j++)
      g4[j] = ydot[j] + h*29*fprimet0[j]/250 + (-112*g1[j]/125 - 54*g2[j]/125 - 2*g3[j]/5)/h;

   lu_solve(a, n, p, g4);

   for (j = 0; j < n; j++) 
      y[j] = y0[j] + 19*g1[j]/9 + g2[j]/2 + 25*g3[j]/108 + 125*g4[j]/108;

   *err = 0;

   for (j = 0; j < n; j++) 
   {
      err_j = fabs(17*g1[j]/54 + 7*g2[j]/36 + 125*g4[j]/108);

      if (err_j > *err)
         *err = err_j;
   }
}

/**
    Метод Розенброка 3(4) порядка с константами Шампэня 
    для решения жестких задач

    \todo
      - предоставить возможность вычислять матрицу Якоби с помощью конечных разностей
      - задавать минимальный, максимальный и начальный шаг
*/

void ode_rosenbrock(
  void (*f)(double, double*, double*), 
  void (*jacobian)(double, double*, double*, double**), 
  size_t n, 
  double t0, 
  double tout, 
  double *y0, 
  double reltol, 
  double abstol, 
  int maxfun, 
  int *rc, 
  double *y, 
  double *ydot,
  double *fprimet, 
  double **fprimey, 
  double **a, 
  double *work,
  size_t *p)
{
   double eps, epsp1, sqrt_eps, err, t, h, hmin, delta;
   int fun_evals;

   eps = 1.0;

   do
   {
     eps = eps/2;
     epsp1 = 1.0 + eps;
   }
   while (epsp1 > 1.0);

   sqrt_eps = sqrt(eps);

   if (fabs(t0) > fabs(tout))
      hmin = sqrt_eps*fabs(t0);
   else 
      hmin = sqrt_eps*fabs(tout); 

   t = t0;
   h = tout - t0;

   /* Основной цикл */

   (*f)(t, y0, ydot);
   (*jacobian)(t, y0, fprimet, fprimey);

   fun_evals = 1;

   while (tout > t0 && t < tout || tout < t0 && t > tout) 
   {
       if (tout > t0 && t + h > tout || tout < 0 && t + h < tout)
          h = tout - t;

       ode_rosenbrock_step(f, n, t, t + h, y0, ydot, fprimet, fprimey, y, &err, a, work, p);

       fun_evals += 2;
       
       if (fun_evals > maxfun)
       {
          *rc = 0;
          return;
       }
       
       if (err <= abstol || err <= reltol*nl_dvector_norm_inf(y, n)) 
       {
          t += h;
          nl_dvector_copy(y0, y, n);

          (*f)(t, y0, ydot);
          (*jacobian)(t, y0, fprimet, fprimey);
         
          fun_evals += 1;
       }
       else
       {
          /* Выбор длины нового шага */

	  if (err == 0)
             delta = 4;
          else 
          {
             if (reltol < eps)
                delta = 0.8*pow(eps/err, 0.25);
             else
                delta = 0.8*pow(reltol/err, 0.25);
             if (delta < 0.1)
               delta = 0.1;
             else if (delta > 4)
               delta = 4;
          }
          h *= delta;

          if (h < hmin)
          {
             *rc = -1;
             return;
          }
       }
   }

   *rc = 1;
}

   
