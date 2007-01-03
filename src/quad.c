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


double quad_gk715_step(double (*f)(double), double a, double b, 
   double is, int maxfunevals, int *rc, int *nfunevals)
{
   double h, m, q1, q2;
   double x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15;
   double f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15;
   
   h = (b - a)/2; 
   m = (a + b)/2;

   x1  = m - 0.9914553711208126*h;
   x2  = m - 0.9491079123427585*h;
   x3  = m - 0.8648644233597691*h;
   x4  = m - 0.7415311855993944*h;
   x5  = m - 0.5860872354676911*h;
   x6  = m - 0.4058451513773972*h;
   x7  = m - 0.2077849550078985*h;
   x8  = m;
   x9  = m + 0.2077849550078985*h;
   x10 = m + 0.4058451513773972*h;
   x11 = m + 0.5860872354676911*h;
   x12 = m + 0.7415311855993944*h;
   x13 = m + 0.8648644233597691*h;
   x14 = m + 0.9491079123427585*h;
   x15 = m + 0.9914553711208126*h;

   f1  = (*f)(x1);
   f2  = (*f)(x2);
   f3  = (*f)(x3);
   f4  = (*f)(x4);
   f5  = (*f)(x5);
   f6  = (*f)(x6);
   f7  = (*f)(x7);
   f8  = (*f)(x8);
   f9  = (*f)(x9);
   f10 = (*f)(x10);
   f11 = (*f)(x11);
   f12 = (*f)(x12);
   f13 = (*f)(x13);
   f14 = (*f)(x14);
   f15 = (*f)(x15);

   *nfunevals += 15;

   q1 = h*(0.02293532201052922*(f1 + f15) + 0.06309209262997855*(f2 + f14)
         + 0.1047900103222502*(f3 + f13) + 0.1406532597155259*(f4 + f12)
         + 0.1690047266392679*(f5 + f11) + 0.1903505780647854*(f6 + f10)
         + 0.2044329400752989*(f7 + f9)  + 0.2094821410847278*f8);

   q2 = h*(0.1294849661688697*(f2 + f14) + 0.2797053914892767*(f4 + f12)
         + 0.3818300505051189*(f6 + f10) + 0.4179591836734694*f8);

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
      q2 = quad_gk715_step(f, a, m, is, maxfunevals, rc, nfunevals);
      if (*rc == 0)
      {
          q2 += quad_gk715_step(f, m, b, is, maxfunevals, rc, nfunevals);
          if (*rc == 0)
             q1 = q2;
      }
   }

   return q1;
}


double quad_gk715(double (*f)(double), double a, double b, double reltol, int maxfunevals,
   int *rc, int *nfunevals)
{
   double m, h, is, q1, q2, err_q1, err_q2, r, s;
   double f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14;
   double f15, f16, f17, f18, f19, f20, f21, f22, f23, f24, f25;
   double x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14;
   double x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25;

   reltol = NL_MAX(reltol, DBL_EPSILON);

   m = (a + b)/2; 
   h = (b - a)/2;

   x1  = m - 0.9914553711208126*h;
   x2  = m - 0.9501292851471813*h;
   x3  = m - 0.9491079123427585*h;
   x4  = m - 0.8912989661489030*h;
   x5  = m - 0.8648644233597691*h;
   x6  = m - 0.7415311855993944*h;
   x7  = m - 0.6068425835417951*h;
   x8  = m - 0.5860872354676911*h;
   x9  = m - 0.4859824687093055*h;
   x10 = m - 0.4058451513773972*h;
   x11 = m - 0.2311385135742917*h;
   x12 = m - 0.2077849550078985*h;
   x13 = m;                       
   x14 = m + 0.2077849550078985*h;
   x15 = m + 0.2311385135742917*h;
   x16 = m + 0.4058451513773972*h;
   x17 = m + 0.4859824687093055*h;
   x18 = m + 0.5860872354676911*h;
   x19 = m + 0.6068425835417951*h;
   x20 = m + 0.7415311855993944*h;
   x21 = m + 0.8648644233597691*h;
   x22 = m + 0.8912989661489030*h;
   x23 = m + 0.9491079123427585*h;
   x24 = m + 0.9501292851471813*h;
   x25 = m + 0.9914553711208126*h;

   f1  = (*f)(x1);
   f2  = (*f)(x2);
   f3  = (*f)(x3);
   f4  = (*f)(x4);
   f5  = (*f)(x5);
   f6  = (*f)(x6);
   f7  = (*f)(x7);
   f8  = (*f)(x8);
   f9  = (*f)(x9);
   f10 = (*f)(x10);
   f11 = (*f)(x11);
   f12 = (*f)(x12);
   f13 = (*f)(x13);
   f14 = (*f)(x14);
   f15 = (*f)(x15);
   f16 = (*f)(x16);
   f17 = (*f)(x17);
   f18 = (*f)(x18);
   f19 = (*f)(x19);
   f20 = (*f)(x20);
   f21 = (*f)(x21);
   f22 = (*f)(x22);
   f23 = (*f)(x23);
   f24 = (*f)(x24);
   f25 = (*f)(x25);

   *nfunevals = 25;

   q1 = h*(0.02293532201052922*(f1 + f25) + 0.06309209262997855*(f3 + f23)
         + 0.1047900103222502*(f5  + f21) + 0.1406532597155259*(f6 + f20)
         + 0.1690047266392679*(f5  + f11) + 0.1903505780647854*(f8 + f18)
         + 0.2044329400752989*(f12 + f14) + 0.2094821410847278*f13);

   q2 = h*(0.1294849661688697*(f3 + f23) + 0.2797053914892767*(f6 + f20)
         + 0.3818300505051189*(f8 + f18) + 0.4179591836734694*f13);

   is = (b-a)/25*(f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8 + f9 + f10 
                + f11 + f12 + f13 + f14 + f15 + f16 + f17 + f18 
                + f19 + f20 + f21 + f22 + f23 + f24 + f25);

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

   return quad_gk715_step(f, a, b, is, maxfunevals, rc, nfunevals);
}


