#include <math.h>
#include "util.h"
#include "optim.h"

void opt_nelder_mead(
  double (*fun)(double*),
  size_t n,  
  double *x0,
  double *f0,
  int initsimplex,
  double **x, /* (n + 1)*n matrix: vertices of the current simplex */
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
      x0[j] = x[0][j];
  }
  else
  {
    for (j = 0; j < n; j++)
      x[0][j] = x0[j];

    for (i = 1; i <= n; i++)
      x[i][i - 1] = x0[i - 1] != 0.0 ? 
        (1 + rel_init_simplex_quota) * x0[i - 1] : abs_init_simplex_quota;
  }

  for (i = 0; i <= n; i++)
    f[i] = (*fun)(x[i]);
  
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
          d = fabs(x[i][j] - x[imin][j]);
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
          xbar[j] += x[i][j];

    for (j = 0; j < n; j++)
      xbar[j] /= (double)n;
    
    /* Compute the reflection point xr */
    
    for (j = 0; j < n; j++)
      xr[j] = (1 + alpha) * xbar[j] - alpha * x[imax][j];

    fr = (*fun)(xr);
    (*nfunevals)++;
    
    if (fr < f[imin])
    {
      /* Compute the expansion point xe */

      for (j = 0; j < n; j++)
        xe[j] = beta * xr[j] + (1 - beta) * xbar[j];

      fe = (*fun)(xe);
      (*nfunevals)++;

      if (fe < f[imin]) /* fe < f[imin] или fr ??? */
      {
        /* Expand the simplex */

        for (j = 0; j < n; j++)
          x[imax][j] = xe[j];
        f[imax] = fe;
      }
      else
      {
        /* Reflect the simplex */

        for (j = 0; j < n; j++)
          x[imax][j] = xr[j];
        f[imax] = fr;
      }
    }
    else /* fr >= f[imin] */
    {
      if (fr < f[imax2])
      {
        /* Reflect the simplex */

        for (j = 0; j < n; j++)
          x[imax][j] = xr[j];
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
            xc[j] = gamma * x[imax][j] + (1 - gamma) * xbar[j];
        }

        fc = (*fun)(xc);
        (*nfunevals)++;
  
        if (fc < fr && fc < f[imax])  /* ... corresponingly */
        {
          /* Contract outside */
          for (j = 0; j < n; j++)
            x[imax][j] = xc[j];
          f[imax] = fc;
        }
        else
        {
          /* Perform a shrink */
          for (i = 0; i <= n; i++)
            if (i != imin)
            {
              for (j = 0; j < n; j++)
                x[i][j] = x[imin][j] + delta * (x[i][j] - x[imin][j]);
              f[i] = (*fun)(x[i]);
            }
          (*nfunevals) += n;
        }
      }
    }

    (*niter)++;
  }

  for (j = 0; j < n; j++)
    x0[j] = x[imin][j];

  *f0 = f[imin];

}
