#include "util.h"
#include "chol.h"

int chol_decomp(double* A, size_t n)
{
  double diag, sum;
  size_t i, j, k;

  if((diag = A[0]) <= 0)
    return -1;

  A[0] = sqrt(diag);

  for(i = 1; i < n; i++)
  {
    diag = A[i*n + i];
    for(j = 0; j < i; j++)
    {
      sum = 0;
      for(k = 0; k < j; k++)
        sum += A[i*n + k]*A[j*n + k];
      A[i*n + j] = (A[i*n + j] - sum)/A[j*n + j];
    }
    for(k = 0; k < i; k++)
    {
      diag -= A[i*n + k]*A[i*n + k];
    }
    if(diag <= 0)
      return -i - 1;
    A[i*n + i] = sqrt(diag);
  }

  return 0;
}

void chol_solve(double* LL, size_t n, double* b)
{
  size_t i, j;
  double t;

  if(n > 1)
  {
    for(i = 0; i < n; i++)
    {
      t = b[i] /= LL[i*n + i];
      for(j = i + 1; j < n; j++)
        b[j] -= LL[j*n + i]*t;
    }
    for(i = n - 1; i >= 1; i--)
    {
      b[i] /= LL[i*n + i];
      t = -b[i];
      for(j = 0; j < i; j++)
        b[j] += LL[i*n + j]*t;
    };
  }
  b[0] /= LL[0];
}
