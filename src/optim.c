#include <math.h>
#include <float.h>
#include "util.h"
#include "optim.h"

double optim_min(double (*fun)(double), double a, double b, double abstol) 
{

  double eps, tol1, tol2, c, d, e, x, u, v, w, fx, fu, fv, fw, xm, p, q, r;


  c = 0.5*(3 - sqrt(5.0));

  eps = sqrt(DBL_EPSILON);

  v = a + c*(b - a);
  w = v;
  x = v;
  e = 0;
  fx = (*fun)(x);
  fv = fx;
  fw = fx;

  /* Основной цикл */

  while (1)
  {
      xm = 0.5*(a + b);
      tol1 = eps*fabs(x) + abstol/3;
      tol2 = 2*tol1;

      if (fabs(x - xm) <= tol2 - 0.5*(b - a)) 
        break;

      if (fabs(e) <= tol1)
      {
         /* Шаг метода золотого сечения */

         if (x >= xm) 
           e = a - x;
         else 
           e = b - x;

         d = c*e;
      }
      else
      {
        /* Построение параболической интерполянты */

        r = (x - w)*(fx - fv);
        q = (x - v)*(fx - fw);
        p = (x - v)*q - (x - w)*r;
        q = 2*(q - r);

        if (q > 0) 
          p = -p;

        q = fabs(q);
        r = e;
        e = d;

        /* Приемлема ли парабола? */

        if (fabs(p) >= fabs(0.5*q*r) || p <= q*(a - x) || p >= q*(b - x)) 
        {
           /* Все-таки метод золотого сечения */

           if (x >= xm) 
             e = a - x;
           else 
             e = b - x;

           d = c*e;
        }
        else
        {
          /* Шаг параболической интерполяции */

          d = p/q;
          u = x + d;

          /* f не должна вычисляться слишком близко к a или b */

          if ((u - a) < tol2)
          {
            if (xm >= x)
              d = tol1;
            else
              d = -tol1;
          }
          else if ((b - u) < tol2)
          { 
            if (xm >= x)
              d = tol1;
            else
              d = -tol1;
          }
        }
      }

      /* f не должна вычисляться слишком близко к x */

      if (fabs(d) >= tol1) 
         u = x + d;
      else
      {
         if (d >= 0) 
           u = x + tol1;
         else
           u = x - tol1;
      }

      fu = (*fun)(u);

      /* Пересчет a, b, v, w, x */

      if (fu <= fx) 
      {
        if (u >= x) 
          a = x;
        else 
          b = x;

        v = w;
        fv = fw;
        w = x;
        fw = fx;
        x = u;
        fx = fu;
      }
      else
      {
        if (u < x) 
          a = u;
        else
          b = u;

        if (fu <= fw || w == x)
        {
          v = w;  
          fv = fw;
          w = u;  
          fw = fu;
        }
        else if (fu <= fv || v == x || v == w) 
        {
          v = u;  
          fv = fu;
        }
      }
   }

   return x;

}

void optim_nelder_mead(
  double (*fun)(double*),
  size_t n,  
  double *x0,
  double *f0,
  int initsimplex,
  double *x, /* (n + 1)*n matrix: vertices of the current simplex */
  double *f, /* n + 1 values of the function in the vertices */
  double tolf, double tolx,
  int maxfunevals, int maxiter,
  int *rc,
  int *nfunevals, int *niter,
  double *work)
{
  double fr, fe, fc;
  double d, diam;
  double *xbar, *xr, *xe, *xc;
  size_t imin, imax, imax2; /* the indices of the highest (worst), 
    next-highest, and lowest (best) vertices of the simplex */
  size_t i, j;

  const double alpha = 1.0;
  const double beta  = 2.0;
  const double gamma = 0.5;
  const double delta = 0.5;
  const double rel_init_simplex_quota = 0.05; /* Params for the initial simplex */
  const double abs_init_simplex_quota = 0.00025;

  xbar = work;
  xr = work + n;
  xe = work + 2*n;
  xc = work + 3*n;

  /* Compute an initial simplex */

  if (initsimplex)
  {
    for (j = 0; j < n; j++)
      x0[j] = x[j]; /* x(0, j) */
  }
  else
  {
    for (j = 0; j < n; j++)
      x[j] = x0[j];

    for (i = 1; i <= n; i++)
      x[i*n + i - 1] = x0[i - 1] != 0.0 ? 
        (1 + rel_init_simplex_quota) * x0[i - 1] : abs_init_simplex_quota;
  }

  for (i = 0; i <= n; i++)
    f[i] = (*fun)(x + i*n);
  
  *rc = 1;
  *nfunevals = n + 1;
  *niter = 0;

  while (*nfunevals < maxfunevals && *niter < maxiter)
  {
    /* Determine which vertex of the current simplex 
       is the highest (worst), next-highest, and lowest (best) */
            // n > 2 !
    if (f[0] < f[1])
    {
      imin = 0;
      imax = 1;
      imax2 = 0;
    }
    else
    {
      imin = 1;
      imax = 0;
      imax2 = 1;
    }
    for (i = 2; i <= n; i++) 
    {
      if (f[i] < f[imin]) 
        imin = i;
      else if (f[i] > f[imax]) 
      {
        imax2 = imax;
        imax = i;
      } 
      else if (f[i] >= f[imax2]) 
        imax2 = i;
    }

    /* Find diameter of the simplex (in infinity-norm) */

    diam = 0.0; 
    for (i = 0; i <= n; i++)
      if (i != imin)
        for (j = 0; j < n; j++)
        {
          d = fabs(x[i*n + j] - x[imin*n + j]);
          if (d > diam)
            diam = d;
        }

    if (f[imax] - f[imin] <= tolf && diam <= tolx)
    {
      *rc = 0;
      break;
    }
    
    /* Compute xbar = average of the n best points */

    for (j = 0; j < n; j++)
      xbar[j] = 0.0;

    for (i = 0; i <= n; i++)
      if (i != imax)
        for (j = 0; j < n; j++)
          xbar[j] += x[i*n + j];

    for (j = 0; j < n; j++)
      xbar[j] /= (double)n;
    
    /* Compute the reflection point xr */
    
    for (j = 0; j < n; j++)
      xr[j] = (1 + alpha) * xbar[j] - alpha * x[imax*n + j];

    fr = (*fun)(xr);
    (*nfunevals)++;
    
    if (fr < f[imin])
    {
      /* Compute the expansion point xe */

      for (j = 0; j < n; j++)
        xe[j] = beta * xr[j] + (1 - beta) * xbar[j];

      fe = (*fun)(xe);
      (*nfunevals)++;

      if (fe < f[imin]) /* fe < f[imin] Ё«Ё fr ??? */
      {
        /* Expand the simplex */

        for (j = 0; j < n; j++)
          x[imax*n + j] = xe[j];
        f[imax] = fe;
      }
      else
      {
        /* Reflect the simplex */

        for (j = 0; j < n; j++)
          x[imax*n + j] = xr[j];
        f[imax] = fr;
      }
    }
    else /* fr >= f[imin] */
    {
      if (fr < f[imax2])
      {
        /* Reflect the simplex */

        for (j = 0; j < n; j++)
          x[imax*n + j] = xr[j];
        f[imax] = fr;
      }
      else /* fr >= f[imax2] */
      {
        /* Perform contraction */

        if (fr < f[imax])
        {
          /* Perform an outside contraction */
          for (j = 0; j < n; j++)
            xc[j] = gamma * xr[j] + (1 - gamma) * xbar[j];
        }
        else
        {  
          /* Perform an inside contraction */
          for (j = 0; j < n; j++)
            xc[j] = gamma * x[imax*n + j] + (1 - gamma) * xbar[j];
        }

        fc = (*fun)(xc);
        (*nfunevals)++;
  
        if (fc < fr && fc < f[imax])  /* ... corresponingly */
        {
          /* Contract outside */
          for (j = 0; j < n; j++)
            x[imax*n + j] = xc[j];
          f[imax] = fc;
        }
        else
        {
          /* Perform a shrink */
          for (i = 0; i <= n; i++)
            if (i != imin)
            {
              for (j = 0; j < n; j++)
                x[i*n + j] = x[imin*n + j] + delta * (x[i*n + j] - x[imin*n + j]);
              f[i] = (*fun)(x + i*n);
            }
          (*nfunevals) += n;
        }
      }
    }

    (*niter)++;
  }

  for (j = 0; j < n; j++)
    x0[j] = x[imin*n + j];

  *f0 = f[imin];

}

