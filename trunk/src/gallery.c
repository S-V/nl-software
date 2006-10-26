#include <stdlib.h>
#include <memory.h>
#include <assert.h>
#include "util.h"
#include "sparse.h"

#define GALLERY_RANDOM(num) (rand()%(int)((num) + 1))

size_t gall_rand_sym(size_t n, size_t nz, size_t *IS, size_t *JS, size_t *M)
{
  size_t i, j, k, t, m = 0;

  memset(IS, 0, (n + 1)*sizeof(size_t));
  memset(M, 0, n*sizeof(size_t));

  for(t = 0; t < nz; t++)
  {
    do
    {
      i = GALLERY_RANDOM(n - 1);
      j = GALLERY_RANDOM(n);
      k = 1 + (i < j ? i:j);
    }
    while(IS[k] == n - k);
    IS[k]++;
  }

  for(i = j = 0; i < n; i++)
  {
    IS[i + 1] += j;
    m++;
    while(j < IS[i + 1])
    {
      do
      {
        k = GALLERY_RANDOM(n - i - 1);
      }
      while(M[k] == m);
      M[k] = m;
      JS[j++] = k + i + 1;
    }
  }
  return m;
}


void gall_laplace(size_t m, size_t n, size_t *IA, size_t *JA, double *AN)
{
  size_t i, j, k, l, N, nz;

  N = n * m;
  nz = 5 * N - 2 * (m + n); // nz = 5*(m-2)*(n-2) + 8*(m-2) + 8*(n-2) + 12;

  // i,j - индексы точек решетки
  // k   - счетчик по уравнениям
  // l   - счетчик по ненулевым элементам
  for (i = 0, k = 0, l = 0; i < m; i++)
    for (j = 0; j < n; j++, k++)
    {
      IA[k] = l;
      if (i > 0) 
      {
        AN[l] = -1.;
        JA[l] = k - n;
        l++;
      }
      if (j > 0) 
      {
        AN[l] = -1;
        JA[l] = k-1;
        l++;
      }
      AN[l] = 4.;
      JA[l] = k;
      l++;
      if (j < n - 1) 
      {
        AN[l] = -1.;
        JA[l] = k + 1;
        l++;
      }
      if (i < m - 1) 
      {
        AN[l] = -1.;
        JA[l] = k + n;
        l++;
      }
    }
  assert(nz == l);
  IA[N] = nz;
}

void gall_laplace_size(size_t m, size_t n, size_t *N, size_t *nz)
{
  *N = m*n;
  *nz = 5*(*N) - 2*(m + n);
}

void gall_laplace_sym(size_t m, size_t n, size_t *IA, size_t *JA, double *AN, double *AD)
{
  size_t i, j, k, l, N, nz;

  N = n * m;
  nz = 2 * N - (m + n); // nz = 2*(m-2)*(n-2) + 4*(m-2) + 4*(n-2) + 12;

  // i,j - индексы точек решетки
  // k   - счетчик по уравнениям
  // l   - счетчик по ненулевым элементам
  for (i = 0, k = 0, l = 0; i < m; i++)
    for (j = 0; j < n; j++, k++)
    {
      IA[k] = l;
      AD[k] = 4.;
      if (j < n - 1) 
      {
        AN[l] = -1.;
        JA[l] = k + 1;
        l++;
      }
      if (i < m - 1) 
      {
        AN[l] = -1.;
        JA[l] = k + n;
        l++;
      }
    }
  assert(nz == l);  
  IA[N] = nz;
}

void gall_laplace_sym_size(size_t m, size_t n, size_t *N, size_t *nz)
{
  *N = m*n;
  *nz = 2*(*N) - (m + n);
}


/*
void laplace_sym(
  size_t gridm,
  size_t gridn,
  size_t *IS,
  size_t *JS,
  double *SN,
  double *SD)
{
  size_t i, j, k, l;

  // i, j - индексы точек решетки
  // k    - счетчик по уравнениям
  // l    - счетчик по ненулевым элементам
  for (i = 0, k = 0, l = 0; i < gridm; i++)
    for (j = 0; j < gridn; j++, k++)
    {
      IS[k] = l;
      SD[k] = 4.;
      if (j < gridn - 1)
      {
        SN[l] = -1.;
        JS[l] = k + 1;
        l++;
      }
      if (i < gridm - 1)
      {
        SN[l] = -1.;
        JS[l] = k + gridn;
        l++;
      }
    }
  IS[k] = l;
}
*/
