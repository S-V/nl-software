#include <math.h>
#include "nl.h"

void interp_pchip(double *x, double *y, size_t n, double *b, double *c, double *d)
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

void interp_spline(double *x, double *y, size_t n, double *b, double *c, double *d, double *work)
{
   double *a, *r, delta, delta1, h, h0, h1, hlst, h2ndlst, *triwork;
   size_t j;
  
   a = work;
   r = work + n;
   triwork = work + 2*n;

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

   band_tridiag(b, a, c, r, d, n, triwork);

   /* Вычисляем b[j] и c[j] */

   for (j = 0; j < n - 1; j++)
   {
      h = x[j + 1] - x[j];

      b[j] = (d[j] - 2*(y[j + 1] - y[j])/h + d[j + 1])/h/h;
      c[j] = (3*(y[j + 1] - y[j])/h - 2*d[j] - d[j + 1])/h;
   }
   
}

size_t find_interval(double *x, size_t n, double xx)
{
   size_t low, up, j;

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

   return low;
}

double interp_eval(double *x, double *y, double *b, double *c, double *d, size_t n, double xx)
{
   size_t j;
   double s;

   j = find_interval(x, n, xx);
   s = xx - x[j];

   return y[j] + s*(d[j] + s*(c[j] + s*b[j]));
}

double interp_quad(double *x, double *y, double *b, double *c, double *d, size_t n, double x1, double x2)
{
   double q, s;
   size_t j, j1, j2;

   j1 = find_interval(x, n, x1);
   j2 = find_interval(x, n, x2);

   if (j1 == j2)
   {
      j = j1;
      s = x2 - x[j];
      q = s*(y[j] + s*(d[j]/2 + s*(c[j]/3 + s*b[j]/4)));
      s = x1 - x[j];
      q -= s*(y[j] + s*(d[j]/2 + s*(c[j]/3 + s*b[j]/4)));
   }
   else
   {
      j = j1;
      s = x[j + 1] - x[j];
      q = s*(y[j] + s*(d[j]/2 + s*(c[j]/3 + s*b[j]/4)));
      s = x1 - x[j];
      q -= s*(y[j] + s*(d[j]/2 + s*(c[j]/3 + s*b[j]/4)));

      for (j = j1 + 1; j < j2; j++)
      {
         s = x[j + 1] - x[j];
         q += s*(y[j] + s*(d[j]/2 + s*(c[j]/3 + s*b[j]/4)));
      }

      j = j2;
      s = x2 - x[j];
      q += s*(y[j] + s*(d[j]/2 + s*(c[j]/3 + s*b[j]/4)));
   }

   return q;
}

