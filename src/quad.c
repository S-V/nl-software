#include <math.h>
#include <float.h>
#include "nl.h"

double quad_simpson_step(double (*f)(double), double a, double b, 
   double fa, double fm, double fb, double is, int maxfunevals, int *rc, int *nfunevals)
{
   double m, h, flm, frm, q1, q2;

   m = (a + b)/2; 
   h = (b - a)/4;

   flm = (*f)(a + h);
   frm = (*f)(b - h);

   *nfunevals += 2;

   q1 = h/1.5 * (fa + 4*fm + fb);
   q2 = h/3 * (fa + 4*(flm + frm) + 2*fm + fb);
   q1 = (16*q2 - q1)/15;

   if (m <= a || b <= m) 
   {
      /* “екущий отрезок не содержит чисел, представимых в машине. 
         “ребуема€ точность не может быть достигнута */
 
      *rc = -1;
   }
   else if (*nfunevals > maxfunevals) 
   {
      /*  оличество вычислений значений подынтегральной функции
         превысило максимально допустимое */
 
      *rc = -2;
   }
   else if (is + (q1 - q2) == is)
   {
      *rc = 0;
   } 
   else
   {
      q2 = quad_simpson_step(f, a, m, fa, flm, fm, is, maxfunevals, rc, nfunevals);
      if (*rc == 0)
      {
         q2 += quad_simpson_step(f, m, b, fm, frm, fb, is, maxfunevals, rc, nfunevals);
         if (*rc == 0)
            q1 = q2;
      }
    
   }

   return q1;
}


double quad_simpson(double (*f)(double), double a, double b, double reltol, int maxfunevals, 
   int *rc, int *nfunevals)
{
   double fa, fm, fb, f1, f2, f3, f4, f5, h, is; 

   reltol = NL_MAX(reltol, DBL_EPSILON);

   fa = (*f)(a);
   fm = (*f)((a + b)/2);
   fb = (*f)(b);

   h = b - a;

   f1 = (*f)(a + 0.9501*h);
   f2 = (*f)(a + 0.2311*h);
   f3 = (*f)(a + 0.6068*h);
   f4 = (*f)(a + 0.4860*h);
   f5 = (*f)(a + 0.8913*h);

   *nfunevals = 8;

   is = h/8*(fa + fm + fb + f1 + f2 + f3 + f4 + f5);

   if (is == 0)
     is = h;

   is = is*reltol/DBL_EPSILON;

   return quad_simpson_step(f, a, b, fa, fm, fb, is, maxfunevals, rc, nfunevals);
}


double quad_lk47_step(double (*f)(double), double a, double b, 
   double fa, double fb, double is, int maxfunevals, int *rc, int *nfunevals)
{
   double h, ll, l, m, r, rr, fll, fl, fm, fr, frr, q1, q2;
   
   const double alpha = sqrt(2./3.); 
   const double beta = 1/sqrt(5.);

   h = (b - a)/2; 
   m = (a + b)/2;

   ll = m - alpha*h;
   l = m - beta*h;
   r = m + beta*h;
   rr = m + alpha*h;

   fll = (*f)(ll);
   fl = (*f)(l);
   fm = (*f)(m);
   fr = (*f)(r);
   frr = (*f)(rr);

   *nfunevals += 5;

   q1 = (h/1470)*(77*(fa + fb) + 432*(fll + frr) + 625*(fl + fr) + 672*fm);
   q2 = (h/6)*(fa + fb + 5*(fl + fr));

   if (m <= a || b <= m)
   {
      /* “екущий отрезок не содержит чисел, представимых в машине. 
         “ребуема€ точность не может быть достигнута */
 
      *rc = -1;
   }
   else if (*nfunevals > maxfunevals) 
   {
      /*  оличество вычислений значений подынтегральной функции
         превысило максимально допустимое */
 
      *rc = -2;
   }
   else if (is + (q1 - q2) == is)
   {
      *rc = 0;
   }
   else
   {
      q2 = quad_lk47_step(f, a, ll, fa, fll, is, maxfunevals, rc, nfunevals);
      if (*rc == 0)
      {
         q2 += quad_lk47_step(f, ll, l, fll, fl, is, maxfunevals, rc, nfunevals);
         if (*rc == 0)
         {
            q2 += quad_lk47_step(f, l, m, fl, fm, is, maxfunevals, rc, nfunevals);
            if (*rc == 0)
            {
               q2 += quad_lk47_step(f, m, r, fm, fr, is, maxfunevals, rc, nfunevals);
               if (*rc == 0)
               {
                  q2 += quad_lk47_step(f, r, rr, fr, frr, is, maxfunevals, rc, nfunevals);
                  if (*rc == 0)
                  {
                     q2 += quad_lk47_step(f, rr, b, frr, fb, is, maxfunevals, rc, nfunevals);
                     if (*rc == 0)
                        q1 = q2;
                  }
               } 
            } 
         } 
      } 
   }

   return q1;
}


double quad_lk47(double (*f)(double), double a, double b, double reltol, int maxfunevals,
   int *rc, int *nfunevals)
{
   double m, h, x1, x2, x3, fa, fb, q1, q2, is, err_q1, err_q2, s, r;
   double f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13;

   const double alpha = sqrt(2./3.); 
   const double beta = 1/sqrt(5.);

   reltol = NL_MAX(reltol, DBL_EPSILON);

   m = (a + b)/2; 
   h = (b - a)/2;

   x1 = .942882415695480; 
   x2 = .641853342345781;
   x3 = .236383199662150;

   f1 = (*f)(a);
   f2 = (*f)(m - x1*h);
   f3 = (*f)(m - alpha*h);
   f4 = (*f)(m - x2*h);
   f5 = (*f)(m - beta*h);
   f6 = (*f)(m - x3*h);
   f7 = (*f)(m);
   f8 = (*f)(m + x3*h);
   f9 = (*f)(m + beta*h);
   f10 = (*f)(m + x2*h);
   f11 = (*f)(m + alpha*h);
   f12 = (*f)(m + x1*h);
   f13 = (*f)(b);

   *nfunevals = 13;

   fa = f1; 
   fb = f13;

   q1 = (h/1470)*(77*(f1 + f13) + 432*(f3 + f11) + 625*(f5 + f9) + 672*f7);
   q2 = (h/6)*(f1 + f13 + 5*(f5 + f9));
   is = h*(0.0158271919734802*(f1 + f13) 
         + 0.0942738402188500*(f2 + f12) 
         + 0.155071987336585*(f3 + f11)
         + 0.188821573960182*(f4 + f10) 
         + 0.199773405226859*(f5 + f9) 
         + 0.224926465333340*(f6 + f8)
         + 0.242611071901408*f7);

   s = NL_SIGN(is); 
   if(s == 0) 
      s = 1; 

   err_q1 = fabs(q1 - is);
   err_q2 = fabs(q2 - is);

   if (err_q2 != 0) 
      r = err_q1/err_q2;
   else
      r = 1;

   if (r > 0 && r < 1)
      reltol = reltol/r;

   is = s*fabs(is)*reltol/DBL_EPSILON;

   if (is == 0) 
      is = b - a;

   return quad_lk47_step(f, a, b, fa, fb, is, maxfunevals, rc, nfunevals);
}

