#include <assert.h>
#include "util.h"
#include "band.h"

void band_tridiag(double *a, double *d, double *c, double *b, double *x, size_t n)
{

  size_t j;
  double bet, *f;

  assert(n > 1);

  f = nl_dvector_create(n);

  f[0] = c[0] / d[0];
  x[0] = b[0] / d[0];

  for(j = 1; j < n - 1; j++)
  {
    bet = d[j] - a[j] * f[j - 1];
    f[j] = c[j] / bet;
    x[j] = (b[j] - a[j] * x[j - 1]) / bet;
  }

  x[n-1] = (b[n-1] - a[n-1] * x[n-2]) / (d[n-1] - a[n-1] * f[n-2]);

  for(j = n - 1; j > 0; j--)
  {
    x[j - 1] -= f[j - 1] * x[j];
  }

  nl_dvector_free(f);
}

void band_mult_col(double **A, size_t n, size_t m1, size_t m2, double *b, double *c)
{
  size_t i, j, k, begin, end;
  
  for (i = 0; i < n; i++) 
  {
    begin = i > m1 ? 0 : m1 - i;
    end = i < n - m2 ? m1 + m2 : n + m1 - i;
    k = i - m1;
    c[i] = 0.0;
    for (j = begin; j <= end; j++) c[i] += A[i][j] * b[j + k];
  }
}

void band_decomp(double **A, size_t n, size_t m1, size_t m2, double **L, size_t *p, int *sgn)
{
  size_t i, j, k, l, m, s;
  double dum;

  m = m1 + m2 + 1;                                                 	

  // Перестраиваем матрицу                                                        
  l = m1;                                                          	
  for (i = 0; i < m1; i++)                                         	
  {   
    for (j = m1 - i; j < m; j++)                          
      A[i][j - l] = A[i][j];                                  
    l--;                                                           	
    for (j = m - l - 1; j < m; j++)                                	
      A[i][j] = 0.0;                                               	
  }     
  
  *sgn = 1;                                                   
  l = m1;                                                     
  
  // Для каждой строки...                                                       
  for (k = 0; k < n; k++)                                     
  {
    dum = A[k][0];                                            
    s = k;                                                    
    if (l < n) l++;                                          
    
    // Ищем главный элемент
    for (j = k + 1; j < l; j++)                                   
    {                                                         
      if (fabs(A[j][0]) > fabs(dum))                          
      {                                                       
        dum = A[j][0];                                        
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
        dum = A[k][j];
        A[k][j] = A[s][j];
        A[s][j] = dum;
      }
    }                                                                                                                                                       
    
    // Гауссово исключение
    for (i = k + 1; i < l; i++)                                                                                                                            
    {
      dum = A[i][0] / A[k][0];
      L[k][i - k - 1] = dum;
      for (j = 1; j < m; j++) 
        A[i][j - 1] = A[i][j] - dum * A[k][j];
      A[i][m - 1] = 0.0;
    }
  }
}

void band_solve(double **A, size_t n, size_t m1, size_t m2, double **L, size_t *p, double *b)
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
      b[i] -= L[k][i - k - 1] * b[k];
  }

  // Обратная подстановка
  l = 1;
  for (i = n; i > 0; i--) 
  {
    dum = b[i - 1];
    for (k = 1; k < l; k++) 
      dum -= A[i - 1][k] * b[k + i - 1];
    b[i - 1] = dum / A[i - 1][0];
    if (l < m) l++;
  }
}

/*
void band_decomp_2(double **a, size_t n, size_t m1, size_t m2, double **al,size_t *indx, int *d)
{
	size_t i,j,k,l;
	size_t mm;
	double dum;

	mm=m1+m2+1;
	l=m1;
	for (i=1;i<=m1;i++) {
		for (j=m1+2-i;j<=mm;j++) a[i-1][j-l-1]=a[i-1][j-1];
		l--;
		for (j=mm-l;j<=mm;j++) a[i-1][j-1]=0.0;
	}
	*d=1;
	l=m1;
	for (k=1;k<=n;k++) {
		dum=a[k-1][0];
		i=k;
		if (l < n) l++;
		for (j=k+1;j<=l;j++) {
			if (fabs(a[j-1][1-1]) > fabs(dum)) {
				dum=a[j-1][0];
				i=j;
			}
		}
		indx[k-1]=i;
		if (i != k) {
			*d = -(*d);
			for (j=1;j<=mm;j++) 
			{
                               dum = a[k-1][j-1];
                               a[k-1][j-1] = a[i-1][j-1];
                               a[i-1][j-1] = dum;
			}
		}
		for (i=k+1;i<=l;i++) {
			dum=a[i-1][1-1]/a[k-1][0];
			al[k-1][i-k-1]=dum;
			for (j=2;j<=mm;j++) a[i-1][j-2]=a[i-1][j-1]-dum*a[k-1][j-1];
			a[i-1][mm-1]=0.0;
		}
	}
}


void band_solve_2(double **A, size_t n, size_t m1, size_t m2, double **L, size_t *p, double *b)
{
  size_t i, k, l, m;
  double dum;

  m = m1 + m2 + 1;
  l = m1;

  // Прямая подстановка
  for (k = 1; k <= n; k++) 
  {
    i=p[k - 1];
    if (i != k) 
    { 
      dum = b[k - 1];
      b[k - 1] = b[i - 1];
      b[i - 1] = dum; 
    }  
    if (l < n) l++;
    for (i = k + 1; i <= l; i++) 
      b[i - 1] -= L[k - 1][i - k - 1] * b[k - 1];
  }

  // Обратная подстановка
  l = 1;
  for (i = n; i >= 1; i--) 
  {
    dum = b[i - 1];
    for (k = 2; k <= l; k++) 
      dum -= A[i - 1][k - 1] * b[k + i - 1 - 1];
    b[i - 1] = dum / A[i - 1][1 - 1];
    if (l < m) l++;
  }
}

*/
