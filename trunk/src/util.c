#include <stdlib.h>
#include <memory.h>
#include <assert.h>
#include "util.h"

const char double_format[] = "%8.3f ";
const char integer_format[] = "%5d ";

double* nl_dvector_create(size_t n)
{
  double* v = (double*)malloc(sizeof(double)*n);
  return v;
}

size_t* nl_xvector_create(size_t n)
{
  size_t* v = (size_t*)malloc(sizeof(size_t)*n);
  return v;
}

void nl_dvector_free(double* v)
{
  free(v);
}

void nl_xvector_free(size_t* v)
{
  free(v);
}

double* nl_dmatrix_create(size_t m, size_t n)
{
  double* A = (double*)malloc(sizeof(double)*n*m);
  return A;
}

void nl_dmatrix_free(double* A)
{
  free(A);
}



void nl_dvector_print(double* vec, size_t n, const char* format)
{
  size_t i;
  const char* f = (format)? format : double_format;
  for(i = 0; i < n; i++)
  {
    printf(f, vec[i]);
  }
  printf("\n");
}

void nl_xvector_print(size_t* vec, size_t n, const char* format)
{
  size_t i;
  const char* f = (format)? format : integer_format;
  for(i = 0; i < n; i++)
  {
    printf(f, vec[i]);
  }
  printf("\n");
}

void nl_dvector_fprint(FILE* file, const double* v, size_t n, const char* format)
{
  size_t i;
  const char* f = (format)? format : double_format;
  for(i = 0; i < n; i++)
  {
    fprintf(file, f, v[i]);
  }
  fprintf(file, "\n");
}

void nl_xvector_fprint(FILE* file, const size_t* v, size_t n, const char* format)
{
  size_t i;
  const char* f = (format)? format : integer_format;
  for(i = 0; i < n; i++)
  {
    fprintf(file, f, v[i]);
  }
  fprintf(file, "\n");
}

int nl_dvector_fwrite(const char* filename, double* v, size_t n, const char* format)
{
  FILE *file = fopen(filename, "w");
  if (file == NULL)
    return 0;
  nl_dvector_fprint(file, v, n, format);
  fclose(file);
  return 1;
}

int nl_xvector_fwrite(const char* filename, size_t* v, size_t n, const char* format)
{
  FILE *file = fopen(filename, "w");
  if (file == NULL)
    return 0;
  nl_xvector_fprint(file, v, n, format);
  fclose(file);
  return 1;
}

int nl_dvector_scan(double* a, size_t n)
{
  size_t j;

  for(j = 0; j < n; j++)
  {
    if(scanf("%lf", &a[j]) == 0)
    {
       return 0;
    }
  }
  return 1;
}

int nl_xvector_scan(size_t* a, size_t n)
{
  size_t j;

  for(j = 0; j < n; j++)
  {
    if(scanf("%lf", &a[j]) == 0)
    {
       return 0;
    }
  }
  return 1;
}

int nl_dvector_fscan(FILE* file, double* a, size_t n)
{
  size_t j;

  for(j = 0; j < n; j++)
  {
    if(fscanf(file, "%lf", a + j) == 0)
    {
       return 0;
    }
  }
  return 1;
}

int nl_xvector_fscan(FILE* file, size_t* a, size_t n)
{
  size_t j;

  for(j = 0; j < n; j++)
  {
    if(fscanf(file, "%lf", &a[j]) == 0)
    {
       return 0;
    }
  }
  return 1;
}

int nl_dvector_fread(const char* filename, double* a, size_t n)
{
  FILE *file = fopen(filename, "r");

  if (file == NULL)
    return 0;

  return nl_dvector_fscan(file, a, n);
}

int nl_xvector_fread(const char* filename, size_t* a, size_t n)
{
  FILE *file = fopen(filename, "r");
  int rc;

  if (file == NULL)
    return 0;

  rc = nl_xvector_fscan(file, a, n); 

  fclose(file);

  return rc;
}


void nl_dmatrix_print(double *A, size_t m, size_t n, const char* format)
{
  size_t i, j, k;
  const char* f = (format)? format : double_format;

  k = 0;

  for(i = 0; i < m; i++)
  {
    for(j = 0; j < n; j++)
    {
      printf(f, A[k++]);
    }
    printf("\n");
  }
}

void nl_dmatrix_fprint(FILE* file, double *A, size_t m, size_t n, const char* format)
{
  size_t i, j, k;
  const char* f = (format)? format : double_format;

  k = 0;

  for(i = 0; i < m; i++)
  {
    for(j = 0; j < n; j++)
    {
      fprintf(file, f, A[k++]);
    }
    fprintf(file, "\n");
  }
}

int nl_dmatrix_fwrite(const char* filename, double *A, size_t m, size_t n, const char* format)
{
  FILE *file = fopen(filename, "w");
  if (file == NULL)
    return 0;
  nl_dmatrix_fprint(file, A, m, n, format);
  fclose(file);
  return 1;
}

void nl_dvector_permute(const double *a, const size_t *P, size_t n, double *b)
{
  size_t k;
  for(k = 0; k < n; k++)
    b[k] = a[P[k]];
}
