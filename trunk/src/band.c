#include <assert.h>
#include "util.h"
#include "band.h"

void band_tridiag(double *a, double *d, double *c, double *b, double *x, size_t n, double *work)
{

  size_t j;
  double bet;

  assert(n > 1);

  work[0] = c[0] / d[0];
  x[0] = b[0] / d[0];

  for(j = 1; j < n - 1; j++)
  {
    bet = d[j] - a[j - 1] * work[j - 1];
    work[j] = c[j] / bet;
    x[j] = (b[j] - a[j - 1] * x[j - 1]) / bet;
  }

  x[n-1] = (b[n - 1] - a[n - 2] * x[n - 2]) / (d[n - 1] - a[n - 2] * work[n-2]);

  for(j = n - 1; j > 0; j--)
  {
    x[j - 1] -= work[j - 1] * x[j];
  }
}

void band_mult_col(double *A, size_t n, size_t m1, size_t m2, double *b, double *c)
{
  size_t i, j, k, m, begin, end;

  m = m1 + m2 + 1;
  
  for (i = 0; i < n; i++) 
  {
    begin = i > m1 ? 0 : m1 - i;
    end = i < n - m2 ? m1 + m2 : n + m1 - i;
    k = i - m1;
    c[i] = 0.0;
    for (j = begin; j <= end; j++) 
      c[i] += A[i*m + j] * b[j + k];
  }
}

void band_decomp(double *A, size_t n, size_t m1, size_t m2, double *L, size_t *p, int *sgn)
{
  size_t i, j, k, l, m, s;
  double dum;

  m = m1 + m2 + 1;                                                 	

  // Перестраиваем матрицу                                                        
  l = m1;                                                          	
  for (i = 0; i < m1; i++)                                         	
  {   
    for (j = m1 - i; j < m; j++)                          
      A[i*m + j - l] = A[i*m + j];                                  
    l--;                                                           	
    for (j = m - l - 1; j < m; j++)                                	
      A[i*m + j] = 0.0;                                               	
  }     
  
  *sgn = 1;                                                   
  l = m1;                                                     
  
  // Для каждой строки...                                                       
  for (k = 0; k < n; k++)                                     
  {
    dum = A[k*m]; /* A(k, 0) */                                            
    s = k;                                                    
    if (l < n) l++;                                          
    
    // Ищем главный элемент
    for (j = k + 1; j < l; j++)                                   
    {                                                         
      if (fabs(A[j*m]) > fabs(dum)) /* A(j, 0) */                         
      {                                                       
        dum = A[j*m];                                        
        s = j;                                                
      }                                                       
    }                                                         
    p[k] = s;                                                 

    // Переставляем строки
    if (s != k)                                               
    {                                                         
      *sgn = -(*sgn);                                         
      for (j = 0; j < m; j++) 
      {
        dum = A[k*m + j];
        A[k*m + j] = A[s*m + j];
        A[s*m + j] = dum;
      }
    }                                                                                                                                                       
    
    // Гауссово исключение
    for (i = k + 1; i < l; i++)                                                                                                                            
    {
      dum = A[i*m] / A[k*m];
      L[k*m + i - k - 1] = dum;
      for (j = 1; j < m; j++) 
        A[i*m + j - 1] = A[i*m + j] - dum * A[k*m + j];
      A[i*m + m - 1] = 0.0;
    }
  }
}

void band_solve(double *A, size_t n, size_t m1, size_t m2, double *L, size_t *p, double *b)
{
  size_t i, k, l, m;
  double dum;

  m = m1 + m2 + 1;
  l = m1;

  // Прямая подстановка
  for (k = 0; k < n; k++) 
  {           
    i = p[k];
    if (i != k) 
    { 
      dum = b[k];
      b[k] = b[i];
      b[i] = dum; 
    }  
    if (l < n) l++;
    for (i = k + 1; i < l; i++) 
      b[i] -= L[k*m + i - k - 1] * b[k];
  }

  // Обратная подстановка
  l = 1;
  for (i = n; i > 0; i--) 
  {
    dum = b[i - 1];
    for (k = 1; k < l; k++) 
      dum -= A[(i - 1)*m + k] * b[k + i - 1];
    b[i - 1] = dum / A[(i - 1)*m];  /* A(i - 1, 0) */
    if (l < m) l++;
  }
}

