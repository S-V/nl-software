#include "util.h"
#include "eig.h"

static void RSWAP(double *a, double *b)
{
  double t;
  t = *a; *a = *b; *b = t;
}

static double MACH_EPS(void)
{
  double mach_eps = 1;
  do
    mach_eps /= 2;
  while((mach_eps + 1) > 1);
  return mach_eps * 2;
}

static void COMDIV(double ar, double ai, double br, double bi, double *cr, double *ci)  
{
  double tmp;

  if (fabs(br) > fabs(bi))
  {
    tmp = bi / br;
    br  = tmp * bi + br;
    *cr  = (ar + tmp * ai) / br;
    *ci  = (ai - tmp * ar) / br;
  }
  else
  {
    tmp = br / bi;
    bi  = tmp * br + bi;
    *cr  = (tmp * ar + ai) / bi;
    *ci  = (tmp * ai - ar) / bi;
  }
}

static double COMABS(double ar, double ai)
{
  if (ar == 0.0 && ai == 0.0) 
    return 0.0;

  ar = fabs(ar);
  ai = fabs(ai);

  if (ai > ar)
    RSWAP(&ai, &ar);

  if (ai == 0.0)
    return ar;
  else
    return ar * sqrt(1.0 + ai / ar * ai / ar);
}

static double PYTHAG(double a, double b)
{
    a = fabs(a);
    b = fabs(b);
    if (a > b) return a * sqrt(1.0 + (b / a) * (b / a));
    else return (b == 0.0 ? 0.0 : b * sqrt(1.0 + (a / b) * (a / b)));
}


void eig_balance(double **A, size_t n, size_t *low, size_t *high, double *scal)
{
  double b2, r, c, f, g, s, BASIS;
  size_t m, k, j, i;
  int iter;

  BASIS = 2;
  b2 = BASIS * BASIS;
  m = 0;
  k = n - 1;
  iter = 1;
  while (iter == 1)
  {
    iter = 0;
    for (j = k + 1; (j--) > 0;)
    {
      r = 0;
      for (i = 0; i <= k; i++)
        if (i != j)  r += fabs(A[j][i]);
      if (r == 0) 
      {
        iter = 1;
        scal[k] = j;
        if (j != k)
        {
          for (i = 0; i <= k; i++) 
	     RSWAP(&(A[i][j]), &(A[i][k]));
          for (i = m; i < n; i++) 
             RSWAP(&(A[j][i]), &(A[k][i]));
        }
        k--;
      }
    }
  }

  iter = 1;
  while (iter == 1)
  {
    iter = 0;
    for (j = m; j <= k; j++)
    {
      c = 0.0;
      for (i = m; i <= k; i++)
        if (i != j)  c += fabs(A[i][j]);
      if (c == 0.0) 
      {
        iter = 1;
        scal[m] = j;
        if (j != m) 
        {
          for (i = 0; i <= k; i++) 
             RSWAP(&(A[i][j]), &(A[i][m]));
          for (i = m; i < n; i++) 
             RSWAP(&(A[j][i]), &(A[m][i]));
        }
        m++;
      }
    }
  }

  (*low) = m;
  (*high) = k;
  for (i = m; i <= k; i++) 
    scal[i] = 1.0;

  iter = 1;
  while (iter == 1)
  {
    iter = 0;
    for (i = m; i <= k; i++)
    {
      c = r = 0.0;
      for (j = m; j <= k; j++)
      {
        if (j != i) 
        {
          c += fabs(A[j][i]);
          r += fabs(A[i][j]);
        }
      }
      g = r / BASIS;
      f = 1.0;
      s = c + r;

      while (c < g)
      {
        f *= BASIS; /* f = r * BASIS; */
        c *= b2;
      }

      g = r * BASIS;
      while (c >= g)
      {
        f /= BASIS;
        c /= b2;
      }

      if ((c + r) / f < 0.95 * s) 
      {
        iter = 1;
        g = 1.0 / f;
        scal[i] *= f;
        for (j = m; j < n; j++) 
          A[i][j] *= g;
        for (j = 0; j <= k; j++)  
          A[j][i] *= f;
      }
    }
  }
}

void eig_hess_reduction(double **A, size_t n, size_t low, size_t high, size_t *perm)
{
  double  x, y;
  size_t m, i, j;

  for(m = low + 1; m < high; m++)
  {
    i = m;
    x = 0.0;
    for(j = m; j <= high; j++)
    {
      if (fabs(A[j][m-1]) > fabs(x)) 
      {
        x = A[j][m-1];
        i = j;
      }
    }

    perm[m] = i;
    if (i != m) 
    {
      for(j = m - 1; j < n; j++)
          RSWAP(&(A[i][j]), &(A[m][j]));
      for(j = 0; j <= high; j++)
          RSWAP(&(A[j][i]), &(A[j][m]));
    }

    if (x != 0.0) 
    {
      for(i = m + 1; i <= high; i++)
      {
        y = A[i][m-1];
        if (y != 0.0) 
        {
          y /= x;
          A[i][m-1] = y;
          for(j = m; j < n; j++)
              A[i][j] = A[i][j] - y * A[m][j];
          for(j = 0; j <= high; j++)
              A[j][m] = A[j][m] + y * A[j][i];
        }
      }
    }
  }
}

void eig_hess_transform_matrix(double **A, size_t n, size_t low, size_t high, 
                               size_t *perm, double **Q)
{
  size_t i, k, j;

  for (i = 0; i < n; i++)
  {
    for (k = 0; k < n; k++)
       Q[i][k] = 0.0;
    Q[i][i] = 1.0;
  }

  for (i = high - 1; i > low; i--)
  {
    j = perm[i];
    for (k = i + 1; k <= high; k++) 
      Q[k][i] = A[k][i - 1];
    if (i != j)
    {
      for (k = i; k <= high; k++)
      {
        Q[i][k] = Q[j][k];
        Q[j][k] = 0.0;
      }
      Q[j][i] = 1.0;
    }
  }
}


// complex division
void eig_vectors(double **A, size_t n, size_t low, size_t high, 
             double *wr, double *wi, double **Q)
{
  size_t i, j, k, l, m, en, na;
  double p, q, r, s, t, w, x, y, z, ra, sa, vr, vi, norm, temp;
  double mach_eps;

  mach_eps = MACH_EPS();

  m = r = s = z = norm = 0.0;

  for (i = 0; i < n; i++)  /* find norm of A */
    for (j = i; j < n; j++)
      norm += fabs(A[i][j]);

  /* if (norm == 0.0) 0.0 matrix */

  for (en = n; (en--) > 0; ) /* transform back */
  {
    p = wr[en];
    q = wi[en];
    na = en - 1;
    if (q == 0.0)
    {
      m = en;
      A[en][en] = 1.0;
      for (i = en; (i--) > 0; )
      {
        w = A[i][i] - p;
        r = A[i][en];
        for (j = m; j <= na; j++)
          r += A[i][j] * A[j][en];
        if (wi[i] < 0.0) 
        {
          z = w;
          s = r;
        }
        else
        {
          m = i;
          if (wi[i] == 0.0)
          {
            if (w != 0.0) 
              temp = w; 
            else 
              temp = mach_eps * norm;
            A[i][en] = -r/temp;   
          }         
          else
          {
            /* Solve the linear system:
              | w   x |  | A[i][en]   |   | -r |
              |       |  |            | = |    |
              | y   z |  | A[i+1][en] |   | -s | */
            x = A[i][i + 1];
            y = A[i + 1][i];
            q = (wr[i] - p)*(wr[i] - p) + wi[i]*wi[i];
            A[i][en] = (x * s - z * r) / q;
            t = A[i][en];
            if (fabs(x) > fabs(z)) 
               temp = (-r -w * t) / x; 
            else 
               temp = (-s -y * t) / z;
            A[i+1][en] = temp;
          }
        } /* wi[i] < 0 */
      } /* i loop */
    }
    else if (q < 0.0) 
    {
      m = na;
      if (fabs(A[en][na]) > fabs(A[na][en]))
      {
        A[na][na] = - (A[en][en] - p) / A[en][na];
        A[na][en] = - q / A[en][na];
      }
      else
      {
        COMDIV(-A[na][en], 0.0, A[na][na]-p, q, &A[na][na], &A[na][en]);
      }
      A[en][na] = 1.0;
      A[en][en] = 0.0;
      for (i = na; (i--) > 0; )
      {
        w = A[i][i] - p;
        ra = A[i][en];
        sa = 0.0;
        for (j = m; j <= na; j++)
        {
          ra = ra + A[i][j] * A[j][na];
          sa = sa + A[i][j] * A[j][en];
        }

        if (wi[i] < 0.0) 
        {
          z = w;
          r = ra;
          s = sa;
        }
        else
        {
          m = i;
          if (wi[i] == 0.0)
            COMDIV(-ra, -sa, w, q, &A[i][na], &A[i][en]);
          else
          {
            /* Solve complex linear system:
               | w+i*q     x | | A[i][na] + i*A[i][en]  |   | -ra+i*sa |
               |             | |                        | = |          |
               |   y    z+i*q| | A[i+1][na]+i*A[i+1][en]|   | -r+i*s   | */
            x = A[i][i + 1];
            y = A[i + 1][i];
            vr = (wr[i] - p) * (wr[i] - p) + wi[i] * wi[i] - q * q;
            vi = 2.0 * q * (wr[i] - p);
            if (vr == 0.0 && vi == 0.0)
              vr = mach_eps * norm * (fabs(w) + fabs(q) + fabs(x) + fabs(y) + fabs(z));

            COMDIV (x * r - z * ra + q * sa, x * s - z * sa -q * ra, vr, vi, &A[i][na], &A[i][en]);
            if (fabs(x) > fabs(z) + fabs(q))
            {
              A[i+1][na] = (-ra - w * A[i][na] + q * A[i][en]) / x;
              A[i+1][en] = (-sa - w * A[i][en] - q * A[i][na]) / x;
            }
            else
            {
              COMDIV (-r - y * A[i][na], -s - y * A[i][en], z, q, &A[i + 1][na], &A[i+1][en]);
            }
          } /* wi[i] = 0 */
        } /* wi[i] < 0 */
      } /* i loop */
    } /* else if q < 0 */
  } /* en loop */

  for (i = 0; i < n; i++) /* Eigenvectors for the evalues for) */
  {
    if (i < low || i > high) /* rows < low and rows > high */
    {
      for (k = i + 1; k < n; k++);
        Q[i][k] = A[i][k];
    }
  }

  j = n;
  while ((j--) > low)
  {
    if (j <= high)
      m = j; 
    else 
      j = high;

    if (wi[j] < 0.0)
    {
      l = j - 1;
      for (i = low; i <= high; i++)
      {
        y = z = 0.0;
        for (k = low; k <= m; k++)
        {
          y += Q[i][k] * A[k][l];
          z += Q[i][k] * A[k][j];
        }
        Q[i][l] = y;
        Q[i][j] = z;
      }
    }
    else
    {
      if (wi[j] == 0.0)
      {
        for (i = low; i <= high; i++)
        {
          z = 0.0;
          for (k = low; k <= m; k++)
            z += Q[i][k] * A[k][j];
          Q[i][j] = z;
        }
      }
    }
  } /* j loop */
}



void eig_hess(double **A, size_t n, size_t low, size_t high,
         double *wr, double *wi, int matq, double **Q, size_t *iter, size_t *rc)
{
  size_t i, j, m, k, l, na, ll, en;
  int cur_iter, flag, MAXIT;
  double p, q, r, s, t, w, x, y, z;
  double mach_eps;

  mach_eps = MACH_EPS();

  MAXIT = 30;
  
  p = q = r = 0.0; 

  for(i = 0; i < n; i++)
  {
    if (i < low || i > high)
    {
      wr[i] = A[i][i];
      wi[i] = 0.0;
      iter[i] = 0;
    }
  }

  en = high;
  t = 0.0;

  flag = 1;
  while (en >= low && flag)
  {
    cur_iter = 0;
    na = en - 1;

    while(1)
    {
      ll = 0;                          
      for(l = en; l > low; l--) /* search for small subdiagonal element */
      {
        if (fabs(A[l][l-1]) <= mach_eps * (fabs(A[l-1][l-1]) + fabs(A[l][l])))
        {
          ll = l; /* save current index */
          break;
        }
      }
      l = ll; /* restore l */

      x = A[en][en];
      if (l == en) /* found one evalue */
      {
        wr[en] = x + t;
        A[en][en] = x + t;
        wi[en] = 0.0;
        iter[en] = cur_iter;
        if (en > 0) 
          en--;
        else
          flag = 0;
        break; /* exit from loop while(1) */
      }

      y = A[na][na];
      w = A[en][na] * A[na][en];

      if (l == na) /* found two evalues */
      {
        p = (y - x) * 0.5;
        q = p * p + w;
        z = sqrt(fabs(q));
        x = x + t;
        A[en][en] = x + t;
        A[na][na] = y + t;
        iter[en] = -cur_iter;
        iter[na] = cur_iter;
        if (q >= 0.0) /* real eigenvalues */
        {
          if (p < 0.0) 
            z = p - z; 
	  else 
            z = p + z;
          wr[na] = x + z;
          wr[en] = x - w / z;
          wi[na] = 0.0;
          wi[en] = 0.0;
          x = A[en][na];
          r = sqrt(x * x + z * z);

          if (matq) 
          {
            p = x / r;
            q = z / r;
            for(j = na; j < n; j++)
            {
              z = A[na][j];
              A[na][j] = q * z + p * A[en][j];
              A[en][j] = q * A[en][j] - p * z;
            }

            for(i = 0; i <= en; i++)
            {
              z = A[i][na];
              A[i][na] = q * z + p * A[i][en];
              A[i][en] = q * A[i][en] - p * z;
            }

            if (matq) for(i = low; i <= high; i++)
            {
              z = Q[i][na];
              Q[i][na] = q * z + p * Q[i][en];
              Q[i][en] = q * Q[i][en] - p * z;
            }
          } 
        }
        else
        {  /* pair of complex */
          wr[na] = x + p;
          wr[en] = x + p;
          wi[na] =   z;
          wi[en] = - z;
        } 

        if (en > 1) 
          en -= 2;
        else
          flag = 0;
        break; /* exit while(1) */
      } /* if l = na */

      if (cur_iter >= MAXIT) 
      {
        iter[en] = MAXIT + 1;
        *rc = en;
        return;
      }                                  

      if (cur_iter == 10 || cur_iter == 20) 
      {
        t += x;
        for(i = low; i <= en; i++)
           A[i][i] -= x;
        s = fabs(A[en][na]) + fabs(A[na][en-2]);
        x = 0.75 * s; 
        y = x;
        w = -0.4375 * s * s;
      }

      cur_iter++;

      for(m = en - 1; m > l; m--)
      {
        z = A[m - 1][m - 1];
        r = x - z;
        s = y - z;
        p = ( r * s - w ) / A[m][m - 1] + A[m - 1][m];
        q = A[m][m] - z - r - s;
        r = A[m + 1][m];
        s = fabs(p) + fabs(q) + fabs (r);
        p = p / s;
        q = q / s;
        r = r / s;
        if (m == l + 1) break;
        if (fabs(A[m - 1][m - 2]) * (fabs(q) + fabs(r)) <= mach_eps * fabs(p)
                 * (fabs(A[m - 2][m - 2]) + fabs(z) + fabs(A[m][m]))) 
          break;
      }
   
      for(i = m + 1; i <= en; i++)
        A[i][i - 2] = 0.0;
      for(i = m + 2; i <= en; i++)
        A[i][i - 3] = 0.0;

      for(k = m - 1; k <= na; k++)
      {
        if (k != m - 1) /* double QR step, for rows l to en */
        {               /* and columns m to en */
          p = A[k][k - 1];
          q = A[k + 1][k - 1];
          if (k != na)
            r = A[k + 2][k - 1]; 
          else 
            r = 0.0;
          x = fabs(p) + fabs(q) + fabs(r);
          if (x == 0.0) 
            continue; /* next k */
          p = p / x;
          q = q / x;
          r = r / x;
        }
        s = sqrt(p*p + q*q + r*r);
        if (p < 0.0) 
          s = -s;

        if (k != m - 1) 
          A[k][k - 1] = -s * x;
        else if (l != m - 1) 
          A[k][k - 1] = -A[k][k - 1];

        p = p + s;
        x = p / s;
        y = q / s;
        z = r / s;
        q = q / p;
        r = r / p;

        for(j = k; j < n; j++) /* modify rows */
        {
          p = A[k][j] + q * A[k + 1][j];
          if (k != na) 
          {
            p = p + r * A[k + 2][j];
            A[k + 2][j] = A[k + 2][j] - p * z;
          }
          A[k + 1][j] = A[k + 1][j] - p * y;
          A[k][j]   = A[k][j] - p * x;
        }

        if (k + 3 < en)  
           j = k + 3; 
        else 
           j = en;
        for(i = 0; i <= j; i++) /* modify columns */
        {
          p = x * A[i][k] + y * A[i][k+1];
          if (k != na) 
          {
            p = p + z * A[i][k+2];
            A[i][k+2] = A[i][k+2] - p * r;
          }
          A[i][k+1] = A[i][k+1] - p * q;
          A[i][k]   = A[i][k] - p;
        }

        if (matq) /* if eigenvectors are needed */
        {
          for(i = low; i <= high; i++)
          {
            p = x * Q[i][k] + y * Q[i][k+1];
            if (k != na)
            {
              p = p + z * Q[i][k+2];
              Q[i][k+2] = Q[i][k+2] - p * r;
            }
            Q[i][k+1] = Q[i][k+1] - p * q;
            Q[i][k]   = Q[i][k] - p;
          }
        }
      } /* k loop */

    } /* while (1<2) */

 } /* while en >= low All evalues found */


  if (matq) /* transform evectors back */
    eig_vectors(A, n, low, high, wr, wi, Q);

  *rc = 0;
}



void eig_balance_inverse(double **Q, size_t n, size_t low, size_t high, double *scal)
{
  size_t  i, j, k;
  double s;

  for (i = low; i <= high; i++)
  {
    s = scal[i];
    for (j = 0; j < n; j++)  
      Q[i][j] *= s;
  }

  for (i = low; (i--) > 0; )
  {
    k = (int)scal[i];    /* round? */
    if (k != i) 
      for (j = 0; j < n; j++)
        RSWAP(&Q[i][j], &Q[k][j]);
  }

  for (i = high + 1; i < n; i++)
  {
    k = (int)scal[i];    /* round? */
    if (k != i)
      for (j = 0; j < n; j++)
        RSWAP(&Q[i][j], &Q[k][j]);
  }
}





/* Complex absolute value */
void eig_norm_Inf(double **Q, size_t n, double *wi)
{
  double maxi, tr, ti;
  size_t i, j;

  j = 0;
  while (j < n)
  {
    if (wi[j] == 0.0)
    {
      maxi = Q[0][j];
      for (i = 1; i < n; i++)
      {
        if (fabs(Q[i][j]) > fabs(maxi))  
          maxi = Q[i][j];
      }

      if (maxi != 0.0)
      {
        maxi = 1.0 / maxi;
        for (i = 0; i < n; i++) 
           Q[i][j] *= maxi;
      }
    }
    else
    {
      tr = Q[0][j];
      ti = Q[0][j + 1];
      for (i = 1; i < n; i++)
      {
        if (COMABS(Q[i][j], Q[i][j + 1]) > COMABS(tr, ti))
        { 
          tr = Q[i][j];
          ti = Q[i][j + 1];
        }
      }

      if (tr != 0.0 || ti != 0.0)
      {
        for (i = 0; i < n; i++)
          COMDIV(Q[i][j], Q[i][j+1], tr, ti, &(Q[i][j]), &(Q[i][j + 1]));
      }
      j++; /* raise j by two */
    }
    j++;
  } 
}



void eig_tridiag_reduction(double **A, size_t n, int matq, double *d, double *a)
{
    size_t l, k, j, i;
    double scale, hh, h, g, f;

    for (i = n - 1; i >= 1; i--) {
        l = i - 1;
        h = scale = 0.0;
        if (l > 0) {
            for (k = 0; k <= l; k++)
                scale += fabs(A[i][k]);
            if (scale == 0.0)
                a[i] = A[i][l];
            else {
                for (k = 0; k <= l; k++) {
                    A[i][k] /= scale;
                    h += A[i][k] * A[i][k];
                }
                f = A[i][l];
                g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
                a[i] = scale * g;
                h -= f * g;
                A[i][l] = f - g;
                f = 0.0;
                for (j = 0; j <= l; j++) {
                    A[j][i] = A[i][j]/h;
                    g = 0.0;
                    for (k = 0; k <= j; k++)
                        g += A[j][k] * A[i][k];
                    for (k = j + 1; k <= l; k++)
                        g += A[k][j] * A[i][k];
                    a[j] = g / h;
                    f += a[j] * A[i][j];
                }
                hh = f / (h + h);
                for (j = 0; j <= l; j++) {
                    f = A[i][j];
                    a[j] = g = a[j] - hh * f;
                    for (k = 0; k <= j; k++)
                        A[j][k] -= (f * a[k] + g * A[i][k]);
                }
            }
        } else
            a[i] = A[i][l];
        d[i] = h;
    }
    d[0] = 0.0;
    a[0] = 0.0;

    if (matq)
      for (i = 0; i < n; i++) 
      {
          l = i;
          if (d[i]) {
              for (j = 0; j < l; j++) {
                  g = 0.0;
                  for (k = 0; k < l; k++)
                      g += A[i][k] * A[k][j];
                  for (k = 0; k < l; k++)
                      A[k][j] -= g * A[k][i];
              }
          }
          d[i] = A[i][i];
          A[i][i] = 1.0;
          for (j = 0; j < l; j++) 
              A[j][i] = A[i][j] = 0.0;
      }
    else
      for (i = 0; i < n; i++) 
      {
          d[i] = A[i][i];
      }
}

void eig_tridiag(double *d, double *a, size_t n, int matq, double **Q, size_t *rc)
{
    size_t m, l, i, k;
    int iter;
    double s, r, p, g, f, dd, c, b;

    for (i = 1; i < n; i++) 
        a[i - 1] = a[i];
    a[n - 1] = 0.0;
    for (l = 0; l < n; l++) {
        iter = 0;
        do {
            for (m = l; m < n - 1; m++) {
                dd = fabs(d[m]) + fabs(d[m + 1]);
                if ((double)(fabs(a[m]) + dd) == dd) break;
            }
            if (m != l) {
                if (iter++ == 30)
                {
                   *rc = l+1;
                   return; /* Too many iterations */
                }
                g = (d[l + 1] - d[l]) / (2.0 * a[l]);
                r = PYTHAG(g, 1.0);
                g = d[m] - d[l] + a[l] / (g + (g >= 0.0 ? fabs(r) : -fabs(r)));
                s = c = 1.0;
                p = 0.0;
                for (i = m; i > l; i--) 
                {
                    f = s * a[i - 1];
                    b = c * a[i - 1];
                    a[i] = (r = PYTHAG(f, g));
                    if (r == 0.0) 
                    {
                        d[i] -= p;
                        a[m] = 0.0;
                        break;
                    }
                    s = f / r;
                    c = g / r;
                    g = d[i] - p;
                    r = (d[i - 1] - g) * s + 2.0 * c * b;
                    d[i] = g + (p = s * r);
                    g = c * r - b;
                  
                    if (matq)
                      for (k = 0; k < n; k++) 
                      {
                          f = Q[k][i];
                          Q[k][i] = s * Q[k][i - 1] + c * f;
                          Q[k][i - 1] = c * Q[k][i - 1] - s * f;
                      }
                }
                if (r == 0.0 && i > l) continue;
                d[l] -= p;
                a[l] = g;
                a[m] = 0.0;
            }
        } while (m != l);
    }
    *rc = 0;
}

/*

#define RADIX 2.0

void eig_balance_old(double **a, size_t n)
{
    size_t j, i;
    double s, r, g, f, c, sqrdx;
    int last;

    sqrdx = RADIX * RADIX;
    last = 0;
    while (last == 0) 
    {
        last = 1;
        for (i = 0; i < n; i++) 
        {
            r = c = 0.0;
            for (j = 0; j < n; j++)
                if (j != i) 
                {
                    c += fabs(a[j][i]);
                    r += fabs(a[i][j]);
                }
            if (c && r) 
            {
                g = r / RADIX;
                f = 1.0;
                s = c + r;
                while (c < g) 
                {
                    f *= RADIX;
                    c *= sqrdx;
                }
                g = r * RADIX;
                while (c > g) 
                {
                    f /= RADIX;
                    c /= sqrdx;
                }
                if ((c + r)  / f < 0.95 * s) 
                {
                    last = 0;
                    g = 1.0 / f;

                    for (j = 0; j < n; j++) 
                        a[i][j] *= g;

                    for (j = 0; j < n; j++) 
                        a[j][i] *= f;
                }
            }
        }
    }
}
#undef RADIX



#define SWAP(g,h) {y=(g);(g)=(h);(h)=y;}

void eig_hess_reduction_old(double **a, size_t n)
{
    size_t m, j, i;
    double y, x;

    for (m = 1; m < n - 1; m++) 
    {
        x = 0.0;
        i = m;
        for (j = m; j < n; j++) 
        {
            if (fabs(a[j][m - 1]) > fabs(x)) 
            {
                x = a[j][m - 1];
                i = j;
            }
        }
        if (i != m) 
        {
            for (j = m - 1; j < n; j++) SWAP(a[i][j],a[m][j])
            for (j = 0; j < n; j++) SWAP(a[j][i],a[j][m])
        }
        if (x) 
        {
            for (i = m + 1; i < n; i++) 
            {
                if ((y = a[i][m - 1]) != 0.0) 
                {
                    y /= x;
                    a[i][m - 1] = y;
                    for (j = m; j < n; j++)
                        a[i][j] -= y * a[m][j];
                    for (j = 0; j < n; j++)
                        a[j][m] += y * a[j][i];
                }
            }
        }
    }
}
#undef SWAP
    

void eig_hess_old(double **a, size_t n, double *wr, double *wi)
{
    size_t nn, l, m, k, j, i, mmin;
    int its;
    double z, y, x, w, v, u, t, s, r, q, p, anorm;

    anorm = fabs(a[0][0]);
    for (i = 1; i < n; i++)
        for (j = (i - 1); j < n; j++)
            anorm += fabs(a[i][j]);
    nn = n;
    t = 0.0;
    while (nn > 0) {
        its = 0;
        do 
        {
            for (l = nn - 1; l > 0; l--) 
            {
                s = fabs(a[l - 1][l - 1]) + fabs(a[l][l]);
                if (s == 0.0) s = anorm;
                if ((double)(fabs(a[l][l - 1]) + s) == s) break;
            }
            x = a[nn - 1][nn - 1];
            if (l == nn - 1) 
            {
                nn--;
                wr[nn] = x + t;
                wi[nn] = 0.0;
            } else {
                y = a[nn - 2][nn - 2];
                w = a[nn - 1][nn - 2] * a[nn - 2][nn - 1];
                if (l == (nn - 2)) {
                    p = 0.5 * (y - x);
                    q = p * p + w;
                    z = sqrt(fabs(q));
                    x += t;
                    if (q >= 0.0) 
                    {
                        z = p + (p >= 0.0 ? z : -z);
                        wr[nn - 2] = wr[nn - 1] = x + z;
                        if (z) wr[nn - 1] = x - w / z;
                        wi[nn - 2] = wi[nn - 1] = 0.0;
                    } 
                    else 
                    {
                        wr[nn - 2] = wr[nn - 1] = x + p;
                        wi[nn - 2] = -(wi[nn - 1] = z);
                    }
                    nn -= 2;
                } 
                else 
                {
                    if (its == 30); // nrerror("Too many iterations in hqr");
                    if (its == 10 || its == 20) 
                    {
                        t += x;
                        for (i = 0; i < nn; i++) a[i][i] -= x;
                        s = fabs(a[nn - 1][nn - 2]) + fabs(a[nn - 2][nn - 3]);
                        y = x = 0.75 * s;
                        w = -0.4375 * s * s;
                    }
                    ++its;
                    for (m = (nn - 3); m >= l; m--) 
                    {
                        z = a[m][m];
                        r = x - z;
                        s = y - z;
                        p = (r * s - w) / a[m + 1][m] + a[m][m + 1];
                        q = a[m + 1][m + 1] - z - r - s;
                        r = a[m + 2][m + 1];
                        s = fabs(p) + fabs(q) + fabs(r);
                        p /= s;
                        q /= s;
                        r /= s;
                        if (m == l) break;
                        u = fabs(a[m][m - 1]) * (fabs(q) + fabs(r));
                        v = fabs(p) * (fabs(a[m - 1][m - 1]) + fabs(z) + fabs(a[m + 1][m + 1]));
                        if ((double)(u + v) == v) break;
                    }
                    for (i = m + 2; i < nn; i++) 
                    {
                        a[i][i - 2] = 0.0;
                        if (i != (m + 2)) a[i][i - 3] = 0.0;
                    }
                    for (k = m; k < nn - 1; k++) 
                    {
                        if (k != m) 
                        {
                            p = a[k][k - 1];
                            q = a[k + 1][k - 1];
                            r = 0.0;
                            if (k != (nn - 2)) r = a[k + 2][k - 1];
                            if ((x = fabs(p) + fabs(q) + fabs(r)) != 0.0) 
                            {
                                p /= x;
                                q /= x;
                                r /= x;
                            }
                        }
                        if ((s = (p >= 0.0 ? sqrt(p * p + q * q + r * r) : 
                                            -sqrt(p * p + q * q + r * r)) ) != 0.0) 
                        {
                            if (k == m) 
                            {
                                if (l != m)
                                a[k][k - 1] = -a[k][k - 1];
                            } 
                            else
                                a[k][k - 1] = -s * x;
                            p += s;
                            x = p/s;
                            y = q/s;
                            z = r/s;
                            q /= p;
                            r /= p;
                            for (j = k; j < nn; j++) 
                            {
                                p = a[k][j] + q * a[k + 1][j];
                                if (k != (nn - 2)) 
                                {
                                    p += r * a[k + 2][j];
                                    a[k + 2][j] -= p * z;
                                }
                                a[k + 1][j] -= p * y;
                                a[k][j] -= p * x;
                            }
                            mmin = nn < k + 4 ? nn : k + 4;
                            for (i = l; i < mmin; i++) 
                            {
                                p = x * a[i][k] + y * a[i][k + 1];
                                if (k != (nn - 2)) 
                                {
                                    p += z * a[i][k + 2];
                                    a[i][k + 2] -= p * r;
                                }
                                a[i][k + 1] -= p * q;
                                a[i][k] -= p;
                            }
                        }
                    }
                }
            }
        } 
        while (l + 2 < nn);
    }
}
*/

#define EIG_ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);

void eig_jacobi(double **A, size_t n, double *w, int matq, double **Q, int *nrot, int *rc)
{
// Функция написана на основе NR jacobi
  size_t i, j, ii;
  int it;
  double tresh, theta, tau, t, sm, s, h, g, c, *b, *z;

  b = nl_dvector_create(n);
  z = nl_dvector_create(n);

  if (matq) 
     /* initialize Q to identity matrix */
     for (i = 0; i < n; i++) 
     {
       for (j = 0; j < n; j++)
         if (i == j) 
           Q[i][j] = 1.0;
         else
           Q[i][j] = 0.0;
     }
  for (i = 0; i < n; i++) 
  {
    b[i] = w[i] = A[i][i];
    z[i] = 0.0;
  }
  *nrot = 0;
  for (it = 0; it < 50; it++) 
  {
    sm = 0.0;
    /* sum of off-diagonal elements */
    for (i = 0; i < n - 1; i++) 
    {
      for (j = i + 1; j < n; j++)
        sm += fabs(A[i][j]);
    }
    if (sm == 0.0) 
    {
      nl_dvector_free(z);
      nl_dvector_free(b);
      *rc = 0;
      return;
    }
    if (it < 4)
      tresh = 0.2*sm / (n*n);
    else
      tresh = 0.0;
    for (i = 0; i < n - 1; i++) 
    {
      for (j = i + 1; j < n; j++) 
      {
        g = 100.0 * fabs(A[i][j]);
        /* after 4 iterations, skip the rotation if the off-diagonal element is small */
        if (it > 4 && (double)(fabs(w[i]) + g) == (double)fabs(w[i])
          && (double)(fabs(w[j]) + g) == (double)fabs(w[j]))
          A[i][j] = 0.0;
        else if (fabs(A[i][j]) > tresh) 
        {
          h = w[j] - w[i];
          if ((double)(fabs(h) + g) == (double)fabs(h))
            t = A[i][j] / h;
          else 
          {
            theta = 0.5*h / A[i][j];
            t = 1.0 / (fabs(theta) + sqrt(1.0 + theta*theta));
            if (theta < 0.0) t = -t;
          }
          c = 1.0 / sqrt(1 + t*t);
          s = t*c;
          tau = s / (1.0 + c);
          h = t*A[i][j];
          z[i] -= h;
          z[j] += h;
          w[i] -= h;
          w[j] += h;
          A[i][j] = 0.0;
          for (ii = 0; ii < i; ii++) 
          {
            EIG_ROTATE(A, ii, i, ii, j)
          }
          for (ii = i + 1; ii < j; ii++) 
          {
            EIG_ROTATE(A, i, ii, ii, j)
          }
          for (ii = j + 1; ii < n; ii++) 
          {
            EIG_ROTATE(A, i, ii, j, ii)
          }
          if (matq)
            for (ii = 0; ii < n; ii++) 
            {
              EIG_ROTATE(Q, ii, i, ii, j)
            }
          ++(*nrot);
        }
      }
    }
    for (i = 0; i < n; i++) 
    {
      b[i] += z[i];
      w[i] = b[i];
      z[i] = 0.0;
    }
  }

  nl_dvector_free(z);
  nl_dvector_free(b);
  *rc = 1;  // Too many iterations
}
#undef EIG_ROTATE

