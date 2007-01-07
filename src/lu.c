#include "nl.h"

#define LU_SWAP(a,b,temp) {temp = (a); (a)=(b); (b)=temp;}
#define LU_SWAPD(a,b) {double temp = (a); (a)=(b); (b)=temp;}

double lu_decomp(double* A, size_t n, size_t* p, int *sgn, double *work)
{
  double cond = 0;
  *sgn = 1;
  if(n != 1)
  {
    size_t i, j, k, ip;
    size_t y, Pivot;
    double norm_A = 0.;
    double Elem, Pr, t;
    double* vy = work;
    double ek, ynorm, znorm;

    // Вычисление нормы матрицы A
    for(j = 0; j < n; j++)
    {
      t = 0.;
      for(i = 0; i < n; i++)
        t += fabs(A[i*n + j]);
      if(norm_A < t)
        norm_A = t;
    };

    for(y = 0; y < n-1; y++)
    {
      // поиск максимального элемента в столбце y
      double Max = fabs(A[y*n + y]);
      Pivot = y;
      for(i = y+1; i < n; i++)
      {
        Elem = fabs(A[i*n + y]);
        if(Elem >= Max)
        {
          Max = Elem;
          Pivot = i;
        };
      };
      if(y != Pivot)
      {
        // перестановка строк Pivot и y  cblas_dswap(CblasRowMajor, A + pivot*n, 1, A + y*n, 1)
        for (j = 0; j < n; j++)
        {
           LU_SWAPD(A[Pivot*n + j], A[y*n + j])
        }
        *sgn = -*sgn;
      }
      p[y] = Pivot;
      if(Max > 0)
      {
        Pr = 1./A[y*n + y];
        for(j = y+1; j < n; j++)
        {
          A[j*n + y] *= Pr;
          t = A[j*n + y];
          for(i = y+1; i < n; i++)
            A[j*n + i] -= t*A[y*n + i];
        }
      };
    }

    // решение системы A^t*y = E
    for(k=0; k < n; k++)
    {
      t = 0.;
      for(i = 0; i<k; i++)
        t += A[i*n + k]*vy[i];
      ek = (t >= 0.)? 1. : -1.;
      /*
      if(A[k*n + k] == 0.)
      {
        return MAX_DOUBLE;
      }
      */
      vy[k] = -(ek+t)/A[k*n + k];
    };
    for(k = n-1; ; k--)
    {
      t = 0.;
      for(i = k+1; i<n; i++)
        t += A[i*n + k]*vy[i];
      vy[k] -= t;
      if(k == 0) break;
    };
    for(i = 0; i < n - 1; i++)
    {
      ip    = p[i];
      if(ip != i) LU_SWAPD(vy[ip], vy[i])
    }

    ynorm = 0.;
    for(i = 0; i < n; i++) ynorm += fabs(vy[i]);

    // решение системы A*z = y
    lu_solve(A, n, p, vy);

    znorm = 0.;
    for(i = 0; i < n; i++) znorm += fabs(vy[i]);

    cond = norm_A*znorm/ynorm;
    if (cond < 1.) cond = 1.;
    
  }
  p[n - 1] = 0;
  return cond;
}


void lu_solve(double* LU, size_t n, size_t* p, double* b)
{
  size_t i, ip, j;
  double t;

  if(n > 1)
  {
    for(i = 0; i < n - 1; i++)
    {
      ip = p[i];
      LU_SWAPD(b[ip], b[i])
    }

    for(i = 0; i < n - 1; i++)
    {
      t = b[i];
      for(j = i+1; j < n; j++)
        b[j] -= LU[j*n + i]*t;
    }
    for(i = n-1; i >= 1; i--)
    {
      b[i] /= LU[i*n + i];
      t = -b[i];
      for(j = 0; j < i; j++)
        b[j] += LU[j*n + i]*t;
    };
  }
  b[0] /= LU[0];
}

double lu_det(double* LU, size_t n, int sgn)
{
  size_t i;
  double det = LU[0]*sgn;
  for(i = 1; i < n; i++)
    det *= LU[i*n + i];
  return det;
}


void lu_invert(double* LU, size_t n, size_t* p, double* B, double *work)
{
  size_t i, j;
  double* col = work;
  for(j = 0; j < n; j++)
  {
    for(i = 0; i<n; i++)
      col[i] = 0.;
    col[j] = 1.;
    lu_solve(LU, n, p, col);
    for(i = 0; i<n; i++)
      B[i*n + j] = col[i];
  };
}

void lu_improve(double* LU, size_t n, size_t* p, double* b, double* x, double* r)
{
  lu_mult_col(LU, n, p, x, r);
  cblas_daxpy(n, -1, b, 1, r, 1);
  lu_solve(LU, n, p, r);
  cblas_daxpy(n, -1, r, 1, x, 1);

  lu_mult_col(LU, n, p, x, r);
  cblas_daxpy(n, -1, b, 1, r, 1);
}

void lu_mult_col(double* LU, size_t n, size_t* p, double* v, double* res)
{
  size_t i, j, ip;
  double t;

  // res = U*v
  for(j = 0; j < n; j++)
  {
    t = 0.;
    for(i = j; i < n; i++)
      t+= LU[j*n + i]*v[i];
    res[j] = t;
  };

  // res = L*res
  for(j = n; j > 0; j--)
  {
    t = res[j - 1];
    for(i = 0; i < j - 1; i++)
      t += LU[(j - 1)*n + i]*res[j - 1];
    res[j - 1] = t;
  };

  // перестановка элементов в res в соответсвии с p
  for(i = 0; i < n - 1; i++)
  {
    if(i != (ip = p[i]))
      LU_SWAPD(res[ip], res[i])
  }
}
