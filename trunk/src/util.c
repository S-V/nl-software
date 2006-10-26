#include <stdlib.h>
#include <memory.h>
#include <assert.h>
#include "util.h"

const char nl_err_str_OK[]                              = "OK";
const char nl_err_str_out_of_memory[]                   = "Out of memory";
const char nl_err_str_inconsistent_size[]               = "Inconsistent size";
const char nl_err_str_incorrect_index[]                 = "Incorrect size_t";
const char nl_err_str_matrix_must_be_square[]           = "Matrix must be square";
const char nl_err_str_matrix_is_singular[]              = "Matrix is singular";
const char nl_err_str_IO[]                              = "Input/output error";
const char nl_err_str_matrix_is_not_positive_definite[] = "Matrix is not positive definite";
const char nl_err_str_diag_elem_must_be_non_zero[]      = "Diagonal elements must be non zero";
const char nl_err_str_unknown_error[]                   = "Uknown error";
const char nl_err_str_runtime_error[]                   = "Runtime error: %s\n";


/*
const char nl_err_str_OK[]                              = "Ошибок нет";
const char nl_err_str_out_of_memory[]                   = "Недостаточно памяти";
const char nl_err_str_inconsistent_size[]               = "Неверный размер";
const char nl_err_str_incorrect_index[]                 = "Неверный индекс";
const char nl_err_str_matrix_must_be_square[]           = "Матрица должна быть квадратной";
const char nl_err_str_matrix_is_singular[]              = "Матрица вырождена";
const char nl_err_str_IO[]                              = "Ошибка при вводе/выводе";
const char nl_err_str_matrix_is_not_positive_definite[] = "Матрица не положительно определенная";
const char nl_err_str_diag_elem_must_be_non_zero[]      = "Диагональные элементы должны быть ненулевыми";
const char nl_err_str_unknown_error[]                   = "Неизвестая ошибка";
const char nl_err_str_runtime_error[]                   = "Ошибка во время исполнения: %s\n";
*/

const char* nl_err_str[nl_err_count] =
{
  nl_err_str_OK,        
  nl_err_str_out_of_memory,
  nl_err_str_inconsistent_size,    
  nl_err_str_incorrect_index,
  nl_err_str_matrix_must_be_square,  
  nl_err_str_matrix_is_singular,
  nl_err_str_matrix_is_not_positive_definite,
  nl_err_str_IO,        
  nl_err_str_diag_elem_must_be_non_zero
};

const char double_format[] = "%8.3f ";
const char integer_format[] = "%5d ";

int nl_err_exit_when_error = 1;

void nl_error_default(nl_err_enum error, int critical)
{
  const char *s;

  if(error < nl_err_count)
    s = nl_err_str[error];
  else
    s = nl_err_str_unknown_error;

  printf(nl_err_str_runtime_error, s);

  if(nl_err_exit_when_error || critical)
  {
    exit(1);
  }
}

void (*nl_error)(nl_err_enum error, int critical) = nl_error_default;

double* nl_dvector_create(size_t n)
{
  double* v = (double*)malloc(sizeof(double)*n);
  if(v == NULL)
    nl_error(nl_err_out_of_memory, 1);
  return v;
}

size_t* nl_xvector_create(size_t n)
{
  size_t* v = (size_t*)malloc(sizeof(size_t)*n);
  if(v == NULL)
    nl_error(nl_err_out_of_memory, 1);
  return v;
}


double* nl_dvector_create_copy(const double* v, size_t n)
{
  double* v_new = nl_dvector_create(n);
  memcpy(v_new, v, n*sizeof(double));
  return v_new;
}


size_t* nl_xvector_create_copy(size_t* v, size_t n)
{
  size_t* v_new = nl_xvector_create(n);
  memcpy(v_new, v, n*sizeof(size_t));
  return v_new;
}


void nl_dvector_free(double* v)
{
  free(v);
}


void nl_xvector_free(size_t* v)
{
  free(v);
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


void nl_dvector_fwrite(const char* filename, double* v, size_t n, const char* format)
{
  FILE *file = fopen(filename, "w");
  if(file == NULL)
    nl_error(nl_err_IO, 1);
  nl_dvector_fprint(file, v, n, format);
  fclose(file);
}


void nl_xvector_fwrite(const char* filename, size_t* v, size_t n, const char* format)
{
  FILE *file = fopen(filename, "w");
  if(file == NULL)
    nl_error(nl_err_IO, 1);
  nl_xvector_fprint(file, v, n, format);
  fclose(file);
}


double* nl_dvector_scan(size_t* n)
{
  size_t y;
  double* v;
  if(scanf("%u", n) != 1)
    nl_error(nl_err_IO, 1);
  v = nl_dvector_create(*n);
  for(y = 0; y < *n; y++)
  {
    if(scanf("%lf", &v[y]) == 0)
      nl_error(nl_err_IO, 1);
  }
  return v;
}


size_t* nl_xvector_scan(size_t* n)
{
  size_t y;
  size_t* v;
  if(scanf("%u", n) != 1)
    nl_error(nl_err_IO, 1);
  v = nl_xvector_create(*n);
  for(y = 0; y < *n; y++)
  {
    if(scanf("%u", &v[y]) == 0)
      nl_error(nl_err_IO, 1);
  }
  return v;
}


double* nl_dvector_fscan(FILE* file, size_t* n)
{
  size_t y;
  double* v;
  if(fscanf(file, "%u", n) != 1)
    nl_error(nl_err_IO, 1);
  v = nl_dvector_create(*n);
  for(y = 0; y < *n; y++)
  {
    if(fscanf(file, "%lf", &v[y]) == 0)
      nl_error(nl_err_IO, 1);
  }
  return v;
}


size_t* nl_xvector_fscan(FILE* file, size_t* n)
{
  size_t y;
  size_t* v;
  if(fscanf(file, "%u", n) != 1)
    nl_error(nl_err_IO, 1);
  v = nl_xvector_create(*n);
  for(y = 0; y < *n; y++)
  {
    if(fscanf(file, "%u", &v[y]) == 0)
      nl_error(nl_err_IO, 1);
  }
  return v;
}


double* nl_dvector_fread(const char* filename, size_t* n)
{
  double* v;
  FILE *file = fopen(filename, "r");
  if(file == NULL)
    nl_error(nl_err_IO, 1);
  v = nl_dvector_fscan(file, n);
  fclose(file);
  return v;
}


size_t* nl_xvector_fread(const char* filename, size_t* n)
{
  size_t* v;
  FILE *file = fopen(filename, "r");
  if(file == NULL)
    nl_error(nl_err_IO, 1);
  v = nl_xvector_fscan(file, n);
  fclose(file);
  return v;
}


//--------------------------------------------------------------------------------
//          Matrix operation
//--------------------------------------------------------------------------------

double** nl_dmatrix_create(size_t m, size_t n)
{
  double** Matrix = (double**)malloc(m*sizeof(double*));
  size_t i;
  if(Matrix == NULL)
    nl_error(nl_err_out_of_memory, 1);
  for (i = 0; i < m; i++)
  {
    Matrix[i] = nl_dvector_create(n);
  };
  return Matrix;
}


double** nl_dmatrix_create_copy(const double** A, size_t m, size_t n)
{
  double** Matrix = (double**)malloc(m*sizeof(double*));
  size_t i;
  for (i = 0; i < m; i++)
  {
      Matrix[i] = nl_dvector_create_copy(A[i], n);
  };
  return Matrix;
}

void nl_dmatrix_copy(double** A, const double** B, size_t m, size_t n)
{
  size_t i;
  size_t  nn = n * sizeof(double);
  for (i = 0; i < m; i++)
  {
    memcpy(A[i], B[i], nn);
  };
}



void nl_dmatrix_free(double** A, size_t m)
{
  size_t i;
  for(i = 0; i < m; i++)
  {
    nl_dvector_free(A[i]);
  };
  free(A);
}


double** nl_dmatrix_convert(double* A, size_t m, size_t n)
{
  double** Matrix = (double**)malloc(m*sizeof(double*));
  size_t i;
  if(Matrix == NULL)
    nl_error(nl_err_out_of_memory, 1);
  for(i = 0; i < m; i++)
    Matrix[i] = A+n*i;
  return Matrix;
}


void nl_dmatrix_print(double **mat, size_t m, size_t n, const char* format)
{
  size_t i,j;
  const char* f = (format)? format : double_format;
  for(j = 0; j < m; j++)
  {
    for(i = 0; i < n; i++)
    {
      printf(f, mat[j][i]);
    }
    printf("\n");
  }
}


void nl_dmatrix_fprint(FILE* file, double** A, size_t m, size_t n, const char* format)
{
  size_t i,j;
  const char* f = (format)? format : double_format;
  for(j = 0; j < m; j++)
  {
    for(i = 0; i < n; i++)
    {
      fprintf(file, f, A[j][i]);
    }
    fprintf(file, "\n");
  }
}


void nl_dmatrix_fwrite(const char* filename, double** A, size_t m, size_t n, const char* format)
{
  FILE *file = fopen(filename, "w");
  if(file == NULL)
    nl_error(nl_err_IO, 1);
  nl_dmatrix_fprint(file, A, m, n, format);
  fclose(file);
}


double** nl_dmatrix_scan(size_t* m, size_t* n)
{
  size_t x,y;
  double** A;
  if(scanf("%u %u\n", m, n) != 2)
    nl_error(nl_err_IO, 1);
  A = nl_dmatrix_create(*m, *n);
  for(y = 0; y < *m; y++)
  {
    for(x = 0; x < *n; x++)
    {
      if(scanf("%lf", &A[y][x]) == 0)
        nl_error(nl_err_IO, 1);
    };
  }
  return A;
}


double** nl_dmatrix_fscan(FILE* file, size_t* m, size_t* n)
{
  size_t x,y;
  double** A;
  if(fscanf(file, "%u %u\n", m, n) != 2)
    nl_error(nl_err_IO, 1);
  A = nl_dmatrix_create(*m, *n);
  for(y = 0; y < *m; y++)
  {
    for(x = 0; x < *n; x++)
    {
      if(fscanf(file, "%lf", &A[y][x]) == 0)
        nl_error(nl_err_IO, 1);
    };
  }
  return A;
}


double** nl_dmatrix_fread(const char* filename, size_t* m, size_t* n)
{
  double** A;
  FILE *file = fopen(filename, "r");
  if(file == NULL)
    nl_error(nl_err_IO, 1);
  A = nl_dmatrix_fscan(file, m, n);
  fclose(file);
  return A;
}


void nl_dmatrix_mult_dvector(double** A, size_t m, size_t n, const double* b, double* res)
{
  size_t x,y;
  double sum;
  for(y = 0; y < m; y++)
  {
    sum = 0;
    for(x = 0; x < n; x++)
    {
      sum += A[y][x]*b[x];
    };
    res[y] = sum;
  };
}

double* nl_dvector_mult_dmatrix(double* b, double** A, size_t m, size_t n)
{
  size_t x,y;
  double* v = nl_dvector_create(n);
  double sum;
  for(x = 0; x < n; x++)
  {
    sum = 0;
    for(y = 0; y < m; y++)
    {
      sum += A[y][x]*b[y];
    };
    v[x] = sum;
  };
  return v;
}


void nl_dvector_sub(double* a, const double* b, size_t n)
{
  double* a_end = a + n;
  while(a < a_end)
    *a++ -= *b++;
}


void nl_dvector_add(double* a, const double* b, size_t n)
{
  double* a_end = a + n;
  while(a < a_end)
    *a++ += *b++;
}


void nl_dmatrix_mult(size_t m, size_t n, size_t k, double** A, double** B, double** Mul)
{
  size_t i,j,l;
  double sum;
  for(j = 0; j < m; j++)
    for(i = 0; i < k; i++)
    {
      sum = 0.;
      for(l = 0; l < n; l++)
        sum += A[j][l]*B[l][i];
      Mul[j][i] = sum;
    };
}

void nl_dmatrix_sub(size_t m, size_t n, double** A, double** B)
{
  size_t i,j;
  for(j = 0; j < m; j++)
    for(i = 0; i < n; i++)
      A[j][i] -= B[j][i];
}


void nl_dmatrix_add(size_t m, size_t n, double** A, double** B)
{
  size_t i,j;
  for(j = 0; j < m; j++)
    for(i = 0; i < n; i++)
      A[j][i] += B[j][i];
}


double nl_dvector_diff_norm1(const double* a, const double* b, size_t n)
{
  size_t i;
  double d = 0.;
  for(i = 0; i < n; i++)
    d += fabs(*a++ - *b++);
  return d;
}


double nl_dvector_diff_norm2(const double* a, const double* b, size_t n)
{
  size_t i;
  double d = 0.,q;
  for(i = 0; i < n; i++)
  {
    q = *a++ - *b++;
    d += q*q;
  }
  return sqrt(d);
}


double nl_dmatrix_diff_norm_inf(double** A, double** B, size_t m, size_t n)
{
  size_t i,j;
  double d = 0.;
  const double *pA, *pB;
  for(j = 0; j < m; j++)
  {
    pA = *A++; // указатель на j-ю строку матрицы A
    pB = *B++; // указатель на j-ю строку матрицы B
    for(i = 0; i < n; i++)
      d += fabs(*pA++ - *pB++);
  }
  return d;
}


double nl_dvector_norm1(const double* a, size_t n)
{
  double norm  = 0.;
  const double* a_end= a + n;
  while(a < a_end)
  {
    norm += fabs(*a++);
  }
  return norm;
}


double nl_dvector_norm2(const double* a, size_t n)
{
  double norm  = 0.;
  const double* a_end = a + n;
  while(a < a_end)
  {
    norm += (*a) * (*a++);
  }
  return sqrt(norm);
}



double nl_dvector_norm_inf(const double* a, size_t n)
{
  double norm  = 0.;
  double abs;
  const double* a_end= a + n;
  while(a < a_end)
  {
    if( (abs = fabs(*a++)) > norm)
      norm = abs;
  }
  return norm;
}



double nl_dvector_dot(const double* a, const double* b, size_t n)
{
  double dot = 0.;
  const double* a_end= a + n;
  while(a < a_end)
  {
    dot += *a++ * *b++;
  }
  return dot;
}



void nl_dvector_copy(double* a, const double* b, size_t n)
{
  memcpy(a, b, n*sizeof(double));
}



void nl_xvector_copy(size_t* a, const size_t* b, size_t n)
{
  memcpy(a, b, n*sizeof(size_t));
}


void nl_dvector_axpy(double alpha, const double *x, const double *y, size_t n, double *z)
{
  const double *x_end = x + n;
  while (x < x_end)
  {
    *z++ = alpha * *x++ + *y++;
  }
}


void nl_dvector_permute(const double *a, const size_t *P, size_t n, double *b)
{
  size_t k;
  for(k = 0; k < n; k++)
    b[k] = a[P[k]];
}
