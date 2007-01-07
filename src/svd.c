#include <math.h>
#include <assert.h>
#include "util.h"
#include "svd.h"

#define SVD_SIGN(a, b) (((b) >= 0) ? fabs(a):-fabs(a))

void svd_decomp(double* A, size_t m, size_t n, double* w, int matu, double* U, 
  int matv, double* V, size_t *ierr, double *work)
{
  size_t i, j, k, l, ii, i1, kk, k1, ll, l1, mn, its;
  double c, f, g, h, s, x, y, z, scale, anorm, tmp;

  (*ierr) = 0;

    for (i = 0; i < m; i++)
      for (j = 0; j < n; j++)
        U[i*n + j] = A[i*n + j];

  // 1 этап
  // ѕриведение к двухдиагональному виду отражени€ми ’аусхолдера
  
  g = 0.0;
  scale = 0.0;
  anorm = 0.0;

  for (i = 0; i < n; i++)
  {
     l = i + 1;
     work[i] = scale * g;
     g = 0.0;
     s = 0.0;
     scale = 0.0;
     if (i < m) 
     {
       for (k = i; k < m; k++)
         scale += fabs(U[k*n + i]);

       if (scale != 0.0) 
       {

          for (k = i; k < m; k++)
          {
             U[k*n + i] /= scale;
             s += U[k*n + i] * U[k*n + i];
          }
       

          f = U[i*n + i];
          g = -SVD_SIGN(sqrt(s),f);
          h = f * g - s;
          U[i*n + i] = f - g;
          if (i != n - 1) 
          {
            for (j = l; j < n; j++)
            {
               s = 0.0;

               for (k = i; k < m; k++)
                  s += U[k*n + i] * U[k*n + j];

               f = s / h;

               for (k = i; k < m; k++)
                  U[k*n + j] += f * U[k*n + i];
            }
          }

          for (k = i; k < m; k++)
            U[k*n + i] = scale * U[k*n + i];

       }
     }
     
     w[i] = scale * g;
     g = 0.0;
     s = 0.0;
     scale = 0.0;

     if (i < m && i != n - 1) 
     {
       for (k = l; k < n; k++)
         scale += fabs(U[i*n + k]);

       if (scale != 0.0)
       {

         for (k = l; k < n; k++)
         {
            U[i*n + k] /= scale;
            s += U[i*n + k] * U[i*n + k];
         }
         

         f = U[i*n + l];
         g = -SVD_SIGN(sqrt(s),f);
         h = f * g - s;
         U[i*n + l] = f - g;

         for (k = l; k < n; k++)
           work[k] = U[i*n + k] / h;

         if (i != m - 1) 
           for (j = l; j < m; j++)
           {
              s = 0.0;

              for (k = l; k < n; k++)
                s += U[j*n + k] * U[i*n + k];

              for (k = l; k < n; k++)
                 U[j*n + k] += s * work[k];
           }

         for (k = l; k < n; k++)
           U[i*n + k] *= scale;
       }
    }
    if (anorm < (tmp = fabs(w[i]) + fabs(work[i])))
      anorm = tmp;
  }
  
  // Ќакопление правосторонних преобразований
  
  if (matv) 
  {

    for (ii = n; ii > 0; ii--)
    {
       i = ii - 1;

       if (i != n - 1)
       {
         if (g != 0.0)
         {
           for (j = l; j < n; j++)
             V[j*n + i] = (U[i*n + j] / U[i*n + l]) / g; 
                // двойное деление обходит возможный машинный ноль

           for (j = l; j < n; j++)
           {
              s = 0.0;

              for (k = l; k < n; k++)
                s += U[i*n + k] * V[k*n + j];

              for (k = l; k < n; k++)
                 V[k*n + j] += s * V[k*n + i];
           }
           
         }

         for (j = l; j < n; j++)
         {
            V[i*n + j] = 0.0;
            V[j*n + i] = 0.0;
         }
         
       }
       V[i*n + i] = 1.0;
       g = work[i];
       l = i;
    }

  }

  // Ќакопление левосторонних преобразований

  if (matu) 
  {

    mn = n;
    if (m < n) mn = m; // mn = min(m, n)

    for (l = mn; l > 0; l--)
    {
       i = l - 1;
       g = w[i];
       if (i != n - 1) 
         for (j = l; j < n; j++)
           U[i*n + j] = 0.0;

       if (g == 0.0) 
       {
         for (j = i; j < m; j++)
           U[j*n + i] = 0.0;
       }
       else
       {
         if (i != mn - 1)
         {
           for (j = l; j < n; j++)
           {
              s = 0.0;

              for (k = l; k < m; k++)
                s += U[k*n + i] * U[k*n + j];
              f = (s / U[i*n + i]) / g; 
                // двойное деление обходит возможный машинный ноль

              for (k = i; k < m; k++)
                 U[k*n + j] += f * U[k*n + i];
           }
         }
         for (j = i; j < m; j++)
           U[j*n + i] /= g;

       }

       U[i*n + i] += 1.0;
    }

  }
     
 // 2 этап
 // ѕриведение двухдиагональной матрицы к диагональному виду 

 for (kk = n; kk > 0; kk--)
 {
   k = kk - 1;
   k1 = k - 1;
   its = 0;

 // ѕроверка возможности расщеплени€

splitting_test:

   for (ll = k + 1; ll > 0; ll--)
   {   
      l = ll - 1;
      l1 = l - 1;
      if (fabs(work[l]) + anorm == anorm) goto convergence_test;
        // work[0] всегда равно нулю, поэтому
        // выхода через конец цикла не будет 
      if (fabs(w[l1]) + anorm == anorm) break;
   }
   // ≈сли l больше чем 0, то work[l] присваиваетс€ нулевое значение

   c = 0.0;
   s = 1.0;

   for (i = l; i <= k; i++)
   {
      f = s * work[i];
      work[i] *= c;

      if (fabs(f) + anorm == anorm) goto convergence_test;

      g = w[i];
      h = sqrt(f * f + g * g);
      assert(h > 0);
      w[i] = h;
      c = g / h;
      s = -f / h;
      if (matu) 
        for (j = 0; j < m; j++)
        {
           y = U[j*n + l1];
           z = U[j*n + i];
           U[j*n + l1] = y * c + z * s;
           U[j*n + i] = -y * s + z * c;
        }
   }


// ѕроверка сходимости

convergence_test:   

   z = w[k];

   if (l != k) 
   {
     // сдвиг выбираетс€ из нижнего углового минора пор€дка 2
     if (its == 30) 
     {
        (*ierr) = k;  
        break;
     }
     its++;
     x = w[l];
     y = w[k1];
     g = work[k1];
     h = work[k];
     f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
     g = sqrt(f * f + 1.0);
     f = ((x - z) * (x + z) + h * (y / (f + SVD_SIGN(g, f)) - h)) / x;
     
     // QR-преобразование
     c = 1.0;
     s = 1.0;

     for (i1 = l; i1 <= k1; i1++)
     {
        i = i1 + 1;
        g = work[i];
        y = w[i];
        h = s * g;
        g = c * g;
        z = sqrt(f * f + h * h);
        work[i1] = z;
        c = f / z;
        s = h / z;
        f = x * c + g * s;
        g = -x * s + g * c;
        h = y * s;
        y = y * c;
        if (matv) 
          for (j = 0; j < n; j++)
          {
             x = V[j*n + i1];
             z = V[j*n + i];
             V[j*n + i1] = x * c + z * s;
             V[j*n + i] = -x * s + z * c;
          }

        z = sqrt(f * f + h * h);
        assert(z > 0);
        w[i1] = z;

        // вращение может быть произвольным, если z равно нулю

        if (z != 0.0)
        {
          c = f / z;
          s = h / z;
        }

        f = c * g + s * y;
        x = -s * g + c * y;
        
        if (matu)
          
          for (j = 0; j < m; j++)
          {
             y = U[j*n + i1];
             z = U[j*n + i];
             U[j*n + i1] = y * c + z * s;
             U[j*n + i] = -y * s + z * c;
          }

     }

     work[l] = 0.0;
     work[k] = f;
     w[k] = x;
     goto splitting_test;
   }

   // —ходимость

   if (z < 0.0) 
   {
     // w[k] делаетс€ неотрицательным
     w[k] = -z;
     if (matv)
       for (j = 0; j < n; j++)
         V[j*n + k] = -V[j*n + k];
   }
 }

}

double svd_cond(double *w, size_t n)
{
  double sigma_max, sigma_min;
  size_t j;
  
  sigma_max = sigma_min = w[0];
  
  for (j = 1; j < n; j++)
  {
    if (w[j] > sigma_max)
      sigma_max = w[j];
    if (w[j] < sigma_min)
      sigma_min = w[j];
  }

  return (sigma_max / sigma_min);
 
}

void svd_correct(double *w, size_t n, double rel_err)
{
  double sigma_max, mach_eps, tau;
  size_t j;
  
  sigma_max = 0;
  
  for (j = 0; j < n; j++)
    if (w[j] > sigma_max)
      sigma_max = w[j];

  if (!rel_err)
  {
    mach_eps = 1;
    do
      mach_eps /= 2;
    while((mach_eps + 1) > 1);
    mach_eps *= 2;
    rel_err = mach_eps * n;
  }

  tau = rel_err * sigma_max;
  
  for (j = 0; j < n; j++)
    if (w[j] < tau)
      w[j] = 0;
}

void svd_least_squares(double *U, double *w, double *V, size_t m, size_t n, double *b, double *x,
  double *work)
{
  size_t jj, j, i;
  double s;

  for (j = 0; j < n; j++) 
  {
    s = 0.0;
    if (w[j]) 
    {
      for (i = 0; i < m; i++) 
        s += U[i*n + j] * b[i];

      s /= w[j];
    }
    work[j] = s;
  }

  for (j = 0; j < n; j++) 
  {
    s = 0.0;
    for (jj = 0; jj < n; jj++) 
      s += V[j*n + jj] * work[jj];
 
    x[j] = s;
  }

}

#undef SVD_SIGN
