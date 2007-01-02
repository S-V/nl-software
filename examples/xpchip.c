#include <stdio.h>
#include <math.h>
#include "nl.h"

/**
  \f$n \ge 3\f$
  x[j] < x[j + 1]

Трудоемкость \f$25n + O(1)$\f
*/

void pchip(double *x, double *y, size_t n, double *b, double *c, double *d)
{
   double w1, w2, h, delta, Delta, Delta1;
   size_t j;

   /* Вычисляем d[j] (производные) */

   Delta = y[1] - y[0];

   for (j = 1; j < n - 1; j++)
   {
      Delta1 = y[j + 1] - y[j];
      if (NL_SIGN(Delta1) == NL_SIGN(Delta))
      {
         w1 = 2*x[j + 1] - x[j] - x[j - 1];
         w2 = x[j + 1] + x[j] - 2*x[j - 1];

         d[j] = (w1 + w2)/(w1/Delta + w2/Delta1);
      }
      else
      {
         d[j] = 0;
      }
      Delta = Delta1;
   }

   /* Вычисляем d[0] и d[n - 1] */

   d[0] = ((x[1] - 2*x[0] + x[2])*(y[1] - y[0]) - (x[1] - x[0])*(y[2] - y[1]))/(x[2]-x[0]);

   if (NL_SIGN(d[0]) != NL_SIGN(y[1] - y[0]))
      d[0] = 0;
   else if (NL_SIGN(y[1] - y[0]) != NL_SIGN(y[2] - y[1]) && fabs(d[0]) > fabs(3*(y[1] - y[0])))
      d[0] = 3*(y[1] - y[0]);

   d[n - 1] = ((2*x[n - 1] - x[n - 2] - x[n - 3])*(y[n - 1] - y[n - 2]) - (x[n - 1] - x[n - 2])*(y[n - 2] - y[n - 3]))/(x[n - 1]-x[n - 3]);

   if (NL_SIGN(d[n - 1]) != NL_SIGN(y[n - 1] - y[n - 2]))
      d[n - 1] = 0;
   else if (NL_SIGN(y[n - 1] - y[n - 2]) != NL_SIGN(y[n - 2] - y[n - 3]) && fabs(d[n - 1]) > fabs(3*(y[n - 1] - y[n - 2])))
      d[n - 1] = 3*(y[n - 1] - y[n - 2]);

   /* Вычисляем b[j] и c[j] */

   for (j = 0; j < n - 1; j++)
   {
      h = x[j + 1] - x[j];
      delta = (y[j + 1] - y[j])/h;

      b[j] = (d[j] - 2*delta + d[j + 1])/h/h;
      c[j] = (3*delta - 2*d[j] - d[j + 1])/h;
   }
   
}

/**
  \f$n \ge 4\f$
  x[j] < x[j + 1]

  Трудоемкость: \f$32n + O(1)$\f
*/

void spline(double *x, double *y, size_t n, double *b, double *c, double *d, double *work)
{
   double *a, *r, delta, delta1, h, h0, h1, hlst, h2ndlst;
   size_t j;
  
   a = work;
   r = work + n;

   /* Вычисляем d[j] (производные) */

   h = x[1] - x[0];
   delta = (y[1] - y[0])/h;

   for (j = 1; j < n - 1; j++)
   {
      c[j] = h;				/* над-диагональные элементы */
      h = b[j - 1] = x[j + 1] - x[j];	/* под-диагональные элементы */
      a[j] = 2*(b[j - 1] + c[j]);	/* диагональ */

      delta1 = (y[j + 1] - y[j])/b[j - 1];
      r[j] = 3*(b[j - 1]*delta + c[j]*delta1);	/* правая часть */
      delta = delta1;
   }

   a[0] = x[2] - x[1];
   c[0] = x[2] - x[0];
   b[n - 2] = x[n - 1] - x[n - 3];
   a[n - 1] = x[n - 2] - x[n - 3];

   /* Учитываем концевые условия 
      (непрерывность второй производной во второй и предпоследней точках) */

   r[0] = 0;
   r[n - 1] = 0;

   h0 = x[1] - x[0];
   h1 = x[2] - x[1];

   r[0] = ((h0 + 2*c[0])*h1*(y[1] - y[0])/h0 + h0*h0*(y[2] - y[1])/h1)/c[0];

   hlst = x[n - 1] - x[n - 2];
   h2ndlst = x[n - 2] - x[n - 3];

   r[n - 1] = (hlst*hlst*(y[n - 2] - y[n - 3])/h2ndlst
            + (2*b[n - 2] + hlst)*h2ndlst*(y[n - 1] - y[n - 2])/hlst)/b[n - 2];

   band_tridiag(b, a, c, r, d, n);

   /* Вычисляем b[j] и c[j] */

   for (j = 0; j < n - 1; j++)
   {
      h = x[j + 1] - x[j];

      b[j] = (d[j] - 2*(y[j + 1] - y[j])/h + d[j + 1])/h/h;
      c[j] = (3*(y[j + 1] - y[j])/h - 2*d[j] - d[j + 1])/h;
   }
   
}

double eval(double *x, double *y, double *b, double *c, double *d, double xx, size_t n)
{
   size_t low, up, j;
   double s;

   low = 0;
   up = n;

   do
   {
      j = (low + up)/2;
      if (x[j] <= xx)
        low = j;
      else
        up = j;
   }
   while (up > low + 1);

   j = low;
   s = xx - x[j];

   return y[j] + s*(d[j] + s*(c[j] + s*b[j]));
}

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
   work = nl_dvector_create(2*n);

   x[0] = -2; y[0] = 0;
   x[1] = -1; y[1] = 0;
   x[2] =  0; y[2] = 1;
   x[3] =  1; y[3] = 0;
   x[4] =  2; y[4] = 0;

   //pchip(x, y, n, b, c, d);
   spline(x, y, n, b, c, d, work);

   printf("  x[j] = ");
   nl_dvector_print(x, n, NULL);
   printf("  y[j] = ");
   nl_dvector_print(y, n, NULL);
   printf("  b[j] = ");
   nl_dvector_print(b, n - 1, NULL);
   printf("  c[j] = ");
   nl_dvector_print(c, n - 1, NULL);
   printf("  d[j] = ");
   nl_dvector_print(d, n, NULL);

   printf("     x       y    \n");
   printf(" ---------------- \n");

   for (xx = -2; xx <= 2; xx += .25)
   {
      printf("  %5.2f  %7.4f\n", xx, eval(x, y, b, c, d, xx, n));
   }

   nl_dvector_free(x);
   nl_dvector_free(y);
   nl_dvector_free(b);
   nl_dvector_free(c);
   nl_dvector_free(d);
   nl_dvector_free(work);

   return 0;
}