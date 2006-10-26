#include "util.h"
#include "chol.h"

void chol_decomp(double** A, size_t n)
{
  double diag, sum;
  size_t i, j, k;

  if((diag = A[0][0]) <= 0)
    nl_error(nl_err_matrix_is_not_positive_definite, 1);
  A[0][0] = sqrt(diag);

  for(i = 1; i < n; i++)
  {
    diag = A[i][i];
    for(j = 0; j < i; j++)
    {
      sum = 0;
      for(k = 0; k < j; k++)
        sum += A[i][k]*A[j][k];
      A[i][j] = (A[i][j] - sum)/A[j][j];
    };
    for(k = 0; k < i; k++)
    {
      diag -= A[i][k]*A[i][k];
    }
    if(diag <= 0)
      nl_error(nl_err_matrix_is_not_positive_definite, 1);
    A[i][i] = sqrt(diag);
  }
}

void chol_solve(double** LL, size_t n, double* b)
{
  size_t i, j;
  double t;

  if(n > 1)
  {
    for(i = 0; i < n; i++)
    {
      t = b[i] /= LL[i][i];
      for(j = i + 1; j < n; j++)
        b[j] -= LL[j][i]*t;
    }
    for(i = n - 1; i >= 1; i--)
    {
      b[i] /= LL[i][i];
      t = -b[i];
      for(j = 0; j < i; j++)
        b[j] += LL[i][j]*t;
    };
  }
  b[0] /= LL[0][0];
}
