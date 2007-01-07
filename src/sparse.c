#include <memory.h>
#include <assert.h>
#include "nl.h"


/*
  Функция сортирует целочисленный массив по возрастанию.  

  - Вход:  
	- \f$v\f$  - сортируемый массив
	- \f$first\f$ - первый элемент, с которого надо начать сортировку
	- \f$last\f$  - последний сортируемый элемент
  - Выход: 
	- \f$v\f$ - отсортированный массив
      
*/
void nl_xvector_qsort(size_t *v, size_t first, size_t last);

/*
  Функция сортируют массив \f$v\f$, а также переставляет
  соответствующие элементы в массиве \f$x\f$. В отличие от #nl_xvector_qsort_x
  элементы массива \f$x\f$ имеют тип double.

  - Вход:  
	- \f$v\f$  - сортируемый массив
	- \f$x\f$  - массив, в котором надо переставлять элементы
	- \f$first\f$ - первый элемет, с которого надо начать сортировку
	- \f$last\f$  - последний сортируемый элемент
  - Выход: 
	- \f$v\f$ - отсортированный массив
	- \f$x\f$ - массив с переставленными элементами
*/
void nl_xvector_qsort_d(size_t *v, double *x, size_t first, size_t last);

/*
  Функция сортируют массив \f$v\f$, а также переставляет
  соответствующие элементы в массиве \f$x\f$. В отличие от #nl_xvector_qsort_d
  элементы массива \f$x\f$ имеют тип #size_t.

  - Вход:  
	- \f$v\f$  - сортируемый массив
	- \f$x\f$  - массив, в котором надо переставлять элементы
	- \f$first\f$ - первый элемет, с которого надо начать сортировку
	- \f$last\f$  - последний сортируемый элемент
  - Выход: 
	- \f$v\f$ - отсортированный массив
	- \f$x\f$ - массив с переставленными элементами
*/
void nl_xvector_qsort_x(size_t *v, size_t *x, size_t first, size_t last);


#define SP_IS_ZERO(val, EPS) ((val > EPS || val < -EPS)? 0:1)


void sp_create(size_t m, size_t nz, size_t **IA, size_t **JA, double **AN)
{
  /* Выделяем память */
  *AN = nl_dvector_create(nz);
  *IA = nl_xvector_create(m + 1);
  *JA = nl_xvector_create(nz);
}

void sp_create_sym(size_t n, size_t nz, size_t **IA, size_t **JA, double **AN, double **AD)
{
  /* выделяем память */
  sp_create(n, nz, IA, JA, AN);
  *AD = nl_dvector_create(n);
}

void sp_sparse(double *A, size_t m, size_t n, size_t** IA, size_t** JA, double** AN, double eps)
{
  double  *pAN;
  size_t  *pIA, *pJA;
  size_t  nzr;
  size_t  i, j, k;

  /* заполнение массивов */
  pIA  = *IA;
  *pIA = 0;
  pJA  = *JA;
  pAN  = *AN;
  k = 0;
  nzr = 0;
  for(i = 0; i < m; i++)
  {
    for(j = 0; j < n; j++)
    {
      if(!SP_IS_ZERO(A[k++], eps))
      {
        *pAN++ = A[k];
        *pJA++ = j;
        nzr++;
      }
    }
    *(++pIA) = nzr;
  }
}



void sp_sparse_sym(double* A, size_t n, size_t** IA, size_t** JA, double** AN, double** AD, double eps)
{
  double  *pAN, *pAD;
  size_t  *pIA, *pJA;
  size_t  nz, nzr;
  size_t  i, j, k;

  /* заполнение массивов */

  pIA  = *IA;
  *pIA = 0;
  pJA  = *JA;
  pAN  = *AN;
  pAD  = *AD;
  k = 0;
  nzr = 0;

  for(i = 0; i < n; i++)
  {
    k += i;
    *pAD = A[k++];

    for(j = i + 1; j < n; j++)
    {
      if(!SP_IS_ZERO(A[k++], eps))
      {
        *pAN++ = A[k];
        *pJA++ = j;
        nzr++;
      }
    }
    *(++pIA) = nzr;
  }
}



void sp_full(size_t* IA, size_t* JA, double* AN, size_t m, size_t n, double* A)
{
  size_t i, j, k;

  /* зануляем все элементы в A */

  for (k = 0; k < m*n; k++)
      A[k] = 0;

  /* заполняем все ненулевые элементы */

  k = 0;
  for (i = 0; i < m; i++)
  {

    for(j = IA[i]; j < IA[i + 1]; j++)
    {
      A[k + JA[j]] = AN[j];
    }

    k += n;
  }
}



void sp_full_sym(size_t* IA, size_t* JA, double* AN, double* AD, size_t n, double* A)
{
  size_t i, j, k;

  /* зануляем все элементы в A */

  for (k = 0; k < n*n; k++)
      A[k] = 0;

  /* заполняем все ненулевые элементы 
     на диагонали и над ней */

  k = 0;
  for (i = 0; i < n; i++)
  {
    A[k + i] = AD[i];

    for(j = IA[i]; j < IA[i + 1]; j++)
    {
      A[k + JA[j]] = AN[j];
    }

    k += n;
  }

  /* заполняем все ненулевые элементы 
     под диагональю */

  for (i = 0; i < n; i++)
    for (j = 0; j < i; j++)
      A[i*n + j] = A[j*n + i];
}



void sp_free(size_t* IA, size_t* JA, double* AN)
{
  nl_dvector_free(AN);
  nl_xvector_free(IA);
  nl_xvector_free(JA);
}



void sp_free_sym(size_t* IA, size_t* JA, double* AN, double* AD)
{
  sp_free(IA, JA, AN);
  nl_dvector_free(AD);
}



size_t sp_nz(size_t m, size_t* IA)
{
  return IA[m];
}

size_t sp_nz_full(double *A, size_t m, size_t n, double eps)
{
  size_t nz, k;

  nz = 0;
  for(k = 0; k < m*n; k++)
  {
     if(!SP_IS_ZERO(A[k], eps))
        nz++;
  }

  return nz;
}

size_t sp_nz_up(double *A, size_t n, double eps)
{
  size_t i, j, k, nz;
  
  nz = 0;
  k = 0;
  for (i = 0; i < n; i++)
  {
    k += i + 1;
    for (j = i + 1; j < n; j++)
    {
      if(!SP_IS_ZERO(A[k++], eps))
        nz++;
    }
  }

  return nz;
}

extern size_t sp_nz_low(double *A, size_t n, double eps)
{
  size_t i, j, k, nz;
  
  nz = 0;
  k = 0;

  for (i = 0; i < n; i++)
  {
    for (j = 0; j < i; j++)
    {
      if(!SP_IS_ZERO(A[k++], eps))
        nz++;
    }
    k += n - i;
  }

  return nz;
}


void sp_print(size_t* IA, size_t* JA, double* AN, size_t m, size_t n, const char* format)
{
  size_t nz = IA[m];
  printf("%u %u\n", m, n);
  nl_xvector_print(IA, m + 1, format);
  nl_xvector_print(JA, nz, format);
  nl_dvector_print(AN, nz, format);
}



void sp_print_list(size_t* IA, size_t* JA, double* AN, size_t m, size_t n, const char* xformat, const char* dformat)
{
  size_t i,j;
  size_t  *pi = IA, *pj = JA;
  size_t  nz = IA[m];
  const char* df  = (dformat)? dformat : double_format;
  const char* xf  = (xformat)? xformat : integer_format;
  double  *pN = AN;
  printf("%u %u %u\n", m, n, nz);
  for(j = 0; j < m; j++)
  {
    for(i = *pi++; i < *pi; i++)
    {
      printf(xf, j); 
      printf(xf, *pj++); 
      printf(df, *pN++);
      printf("\n");
    }
  }
}



void sp_fprint(FILE* file, size_t* IA, size_t* JA, double* AN, size_t m, size_t n, const char* format)
{
  size_t nz = IA[m];
  fprintf(file, "%u %u\n", m, n);
  nl_xvector_fprint(file, IA, m+1, format);
  nl_xvector_fprint(file, JA, nz, format);
  nl_dvector_fprint(file, AN, nz, format);
}



void sp_fprint_list(FILE* file, size_t* IA, size_t* JA, double* AN, size_t m, size_t n, const char* xformat, const char* dformat)
{
  size_t i,j;
  size_t  *pi = IA, *pj = JA;
  size_t  nz = IA[m];
  const char* df  = (dformat)? dformat : double_format;
  const char* xf  = (xformat)? xformat : integer_format;
  double  *pN = AN;
  fprintf(file, "%u %u %u\n", m, n, nz);
  for(j = 0; j < m; j++)
  {
    for(i = *pi++; i < *pi; i++)
    {
      fprintf(file, xf, j);
      fprintf(file, xf, *pj++);
      fprintf(file, df, *pN++);
      fprintf(file, "\n");
    }
  }
}



void sp_print_sym(size_t* IA, size_t* JA, double* AN, double *AD, size_t n, const char* format)
{
  size_t nz = IA[n];
  printf("%u\n", n);
  nl_xvector_print(IA, n+1, 0);
  nl_xvector_print(JA, nz, 0);
  nl_dvector_print(AN, nz, format);
  nl_dvector_print(AD,  n, format);
}



void sp_print_list_sym(size_t* IA, size_t* JA, double* AN, double* AD, size_t n, const char* xformat, const char* dformat)
{
  size_t i,j;
  size_t  *pi = IA, *pj = JA;
  const char* df  = (dformat)? dformat : double_format;
  const char* xf  = (xformat)? xformat : integer_format;
  double  *pN = AN;
  printf("%u %u\n", n, IA[n]);
  for(j = 0; j < n; j++)
  {
    printf(xf, j);
    printf(xf, j);
    printf(df, *AD++);
    printf("\n");
    for(i = *pi++; i < *pi; i++)
    {
      printf(xf, j);
      printf(xf, *pj++);
      printf(df, *pN++);
      printf("\n");
    }
  }
}



void sp_fprint_sym(FILE* file, size_t* IA, size_t* JA, double* AN, double* AD, size_t n, const char* format)
{
  size_t nz = IA[n];
  fprintf(file, "%u\n", n);
  nl_xvector_fprint(file, IA, n+1, 0);
  nl_xvector_fprint(file, JA, nz, 0);
  nl_dvector_fprint(file, AN, nz, format);
  nl_dvector_fprint(file, AD,  n, format);
}



void sp_fprint_list_sym(FILE* file, size_t* IA, size_t* JA, double* AN, double* AD, size_t n, const char* xformat, const char* dformat)
{
  size_t i,j;
  size_t  *pi = IA, *pj = JA;
  const char* df = (dformat)? dformat : double_format;
  const char* xf = (xformat)? xformat : integer_format;
  double  *pN = AN;
  fprintf(file, "%u %u\n", n, IA[n]);
  for(j = 0; j < n; j++)
  {
    printf(xf, j);
    printf(xf, j);
    printf(df, *AD++);
    printf("\n");
    for(i = *pi++; i < *pi; i++)
    {
      fprintf(file, xf, j);
      fprintf(file, xf, *pj++);
      fprintf(file, df, *pN++);
      fprintf(file, "\n");
    }
  }
}



size_t sp_list(size_t* IA, size_t* JA, double* AN, size_t m, double* A, size_t* I, size_t* J)
{
  size_t i,j;
  size_t  *pi = IA;
  for(j = 0; j < m; j++)
  {
    for(i = *pi++; i < *pi; i++)
    {
      *A++  = *AN++;
      *I++  = j;
      *J++  = *JA++;
    }
  }
  return *pi;
}




void sp_convert(size_t nz, double* A, size_t* I, size_t* J, size_t m, size_t* IA, size_t* JA, double* AN)
{
/* I, J - индексы строк и столбцов ненулевых элементов матрицы */
  size_t i,j,k;
  size_t  *pI = I, *pJ = J;
  double  *pA = A;
  size_t  *pi;

  pi  = IA;

  memset(IA, 0, sizeof(size_t)*(m+1));
  for(i = 0; i < nz; i++)
  {
    j = *pI++ + 2;
    if(j <= m)
      IA[j]++;
  }

  pi += 2;
  for(i = 3; i <= m; i++, pi++)
    *(pi+1) += *pi;

  pI = I;
  for(i = 0; i < nz; i++)
  {
    j  = *pI++ + 1;
    k  = IA[j]++;
    AN[k]  = *pA++;
    JA[k]  = *pJ++;
  }
}



void sp_transpose(size_t* IA, size_t* JA, double* AN, size_t m, size_t n, size_t* IAT, size_t* JAT, double* ATN)
{
  size_t i,j,jj,k;
  size_t  *pIAT;
  size_t  *pIA = IA, *pJA = JA;
  double  *pAN = AN;

  size_t  nz = IA[m];
  pIAT  = IAT;

  memset(IAT, 0, sizeof(size_t)*(n+1));
  for(i = 0; i < nz; i++)
  {
    j = *pJA++ + 2;
    if(j <= n)
      IAT[j]++;
  }

  pIAT += 2;
  for(i = 3; i <= n; i++, pIAT++)
    *(pIAT+1) += *pIAT;

  pJA = JA;
  for(j = 0; j < m; j++)
  {
    size_t  z = *pIA++;
    for(i = z; i < *pIA; i++)
    {
      jj  = *pJA++ + 1;
      k  = IAT[jj]++;
      ATN[k]  = *pAN++;
      JAT[k]  = j;
    }
  }
}



void sp_order(size_t* IA, size_t* JA, double* AN, size_t m)
{
  size_t i;

  for(i = 0; i < m; i++)
  {
    if (IA[i + 1] > IA[i])
      nl_xvector_qsort_d(JA, AN, IA[i], IA[i + 1] - 1);
  }
}


void sp_order_m(size_t* IA, size_t* JA, size_t m)
{
  size_t i,z;
  for(i = 0; i < m; i++)
  {
    z = *IA++;
    if(z < *IA)
      nl_xvector_qsort(JA, z, *IA - 1);
  }
}



int sp_add_symb(
  size_t *IA, 
  size_t *JA, 
  size_t *IB, 
  size_t *JB, 
  size_t m, 
  size_t n, 
  size_t *IC, 
  size_t *JC, 
  size_t size_C,
  size_t *xwork)
{
  size_t i, j, nzC = 0, z, q;
  size_t*  pIA = IA;
  size_t*  pJA = JA;
  size_t*  pIB = IB;
  size_t*  pJB = JB;
  size_t*  pIC = IC;
  size_t*  pJC = JC;
  memset(xwork, -1, sizeof(*xwork)*n);

  *pIC++ = 0;
  for(i = 0; i < m; i++)
  {
    z = *pIA++;
    for(j = z; j < *pIA; j++)
    {
      if(++nzC > size_C)
      {
        return 0;
      }
      q  = *pJA++;
      xwork[q]  = i;
      *pJC++  = q;
    }

    z = *pIB++;
    for(j = z; j < *pIB; j++)
    {
      q  = *pJB++;
      if(xwork[q] != i)
      {
        if(++nzC > size_C)
        {
          return 0;
        }
        *pJC++ = q;
      }
    }
    assert(nzC == (size_t)(pJC - JC));
    *pIC++ = pJC - JC;
  }
  return 1;
}



void sp_add_num(
  size_t* IA, 
  size_t* JA, 
  double* AN,
  size_t* IB, 
  size_t* JB, 
  double* BN, 
  size_t m, 
  size_t n,
  size_t* IC, 
  size_t* JC, 
  double* CN,
  double *work)
{
  size_t i, j;
  size_t*  pJC, *pJC_;
  size_t*  pIA = IA;
  size_t*  pJA = JA;
  size_t*  pIB = IB;
  size_t*  pJB = JB;
  double* pAN = AN;
  double* pBN = BN;
  double* pCN = CN;
  size_t*  pIC  = IC;
  size_t  z, zC;
  pJC_ = pJC = JC;

  for(i = 0; i < m; i++)
  {
    zC  = *pIC++;
    for(j = zC; j < *pIC; j++)
      work[*pJC++] = 0.;

    z  = *pIA++;
    for(j = z; j < *pIA; j++)
      work[*pJA++] = *pAN++;

    z  = *pIB++;
    for(j = z; j < *pIB; j++)
      work[*pJB++] += *pBN++;

    for(j = zC; j < *pIC; j++)
      *pCN++ = work[*pJC_++];
  }
}

void sp_mult_col(size_t* IA, size_t* JA, double* AN, double* b, size_t m, double* c)
{
  size_t i, j;
  double sum;
  for(i = 0; i < m; i++)
  {
    sum = 0;
    for(j = *IA++; j < *IA; j++)
    {
      sum += (*AN++) * b[*JA++];
    }
    *c++ = sum;
  }
}



void sp_mult_row(size_t* IA, size_t* JA, double* AN, double* b, size_t m, size_t n, double* c)
{
  size_t i,j, z;
  memset(c, 0, sizeof(*c)*n);
  for(i = 0; i < m; i++)
  {
    for(j = *IA++; j < *IA; j++)
    {
      z = *JA++;
      c[z] += (*AN++) * b[i];
    }
  }
}



void sp_mult_col_sym(size_t* IA, size_t* JA, double* AN, double* AD, double* b, size_t n, double* c)
{
  size_t i,j;
  double* pc = c;
  double* pb = b;

  for(i = 0; i < n; i++)
    *pc++ = *AD++ * (*pb++);

  pc = c;
  pb = b;
  for(i = 0; i < n; i++)
  {
    double u = 0;
    for(j = *IA++; j < *IA; j++)
    {
      u += (*AN) * b[*JA];
      c[*JA++] += (*AN++) * (*pb);
    }
    *pc++ += u;
    pb++;
  }
}



int sp_mult_symb(
  size_t* IA, 
  size_t* JA, 
  size_t* IB, 
  size_t* JB, 
  size_t m, 
  size_t k, 
  size_t* IC, 
  size_t* JC, 
  size_t size_C,
  size_t *xwork)
{
  size_t*  pIC = IC;
  size_t  nzC = 0;
  size_t  i,j,l,z,s;
  size_t*  pIA = IA;
  size_t*  pJA = JA;
  size_t*  pJC = JC;

  memset(xwork, -1, sizeof(size_t)*k);

  *pIC++ = 0;
  for(i = 0; i < m; i++)
  {
    z = *pIA++;
    for(j = z; j < *pIA; j++)
    {
      for(l = IB[*pJA]; l < IB[*pJA + 1]; l++)
      {
        s = JB[l];
        if(xwork[s] != i)
        {
          if(++nzC > size_C)
          {
            /*nl_error(nl_err_inconsistent_size, 0);*/
            return 0;
          }
          *pJC++ = s;
          xwork[s] = i;
        }
      }
      pJA++;
    }
    *pIC++ = nzC;
  }
  return 1;
}



void sp_mult_num(size_t* IA, size_t* JA, double* AN, size_t* IB, size_t* JB, double* BN, size_t* IC, size_t* JC, size_t m, size_t k, double* CN, double *work)
{
  size_t  i,j,l,zC;
  size_t  *pJC1;
  size_t  *pIC = IC;
  size_t  *pJC = pJC1 = JC;
  size_t  *pJA, *pJA_end;
  double  *pAN;
  size_t  *pJB, *pJB_end;
  double  *pBN;

  double  d;

  for(i = 0; i < m; i++)
  {
    zC = *pIC++;
    if(zC != *pIC)
    {
      pJA_end = JC + *pIC;
      while(pJC < pJA_end)
        work[*pJC++] = 0.0;

      j  = IA[i];
      pJA  = JA + j;
      pAN = AN + j;
      pJA_end  = JA + IA[i+1];
      while(pJA < pJA_end)
      {
        d  = *pAN++;
        l  = IB[*pJA];
        pJB  = JB + l;
        pBN  = BN + l;
        pJB_end = JB + IB[*pJA+1];
        while(pJB < pJB_end)
        {
          work[*pJB++] += d*(*pBN++);
        }
        pJA++;
      }

      for(j = zC; j < *pIC; j++)
        CN[j] = work[*pJC1++];
    }
  }
}

void sp_permute_sym(
	size_t n,
	size_t *IA,
	size_t *JA,
	size_t *IB,
	size_t *JB,
	double *AN,
	double *AD,
	double *BN,
	double *BD,
	size_t *IP)
{
	size_t i, j, k, imin, imax;

	memset(IB, 0, n*sizeof(size_t));

	for(i = j = 0; i < n; i++)
	{
		while(j < IA[i + 1])
		{
			k = JA[j];
			IB[NL_MIN(IP[i], IP[k]) + 1]++;
			j++;
		}
	}

	for(i = 1; i < n; i++) IB[i] += IB[i - 1];

	for(i = j = 0; i < n; i++)
	{
		BD[IP[i]] = AD[i];
		while(j < IA[i + 1])
		{
			k = JA[j];
			imin = NL_MIN(IP[i], IP[k]);
			imax = NL_MAX(IP[i], IP[k]);
			JB[IB[imin]] = imax;
			BN[IB[imin]] = AN[j];
			IB[imin]++;
			j++;
		}
	}

	for(i = n; i > 0; i--) IB[i] = IB[i - 1];
	IB[i] = 0;
}




int sp_gauss_seidel(size_t* IA, size_t* JA, double* AN, double* AD, double* b, size_t n, double f, double eps, int max_iter, double* x)
{
  size_t  i,j;
  int    iter;
  double  norm_b  = 0.;
  double  *px    = x;
  double  *pb    = b;
  double  *pAD  = AD;
  size_t  *pIA;
  size_t  *pJA;
  double  *pAN;
  double  e,u;

  /* начальное приближение к решению */
  for(i = 0; i < n; i++)
  {
    norm_b += *pb * *pb;
    *px ++ = *pb++ / *pAD++;
  }
  eps *= eps * norm_b;

  iter = 0;
  do
  {
    pAD = AD;
    pIA = IA;
    pb  = b;
    pJA  = JA;
    px  = x;
    pAN  = AN;
    e  = 0;
    for(i = 0; i < n; i++)
    {
      u  = *pb++;
      for(j = *pIA++; j < *pIA; j++)
      {
        u -= *pAN++ * x[*pJA++];
      }
      u = u/(*pAD++) - *px;
      e += u*u;
      *px++ += f*u;
    }
    iter++;
  } while(iter < max_iter && e >= eps);
  return iter;
}



int sp_gauss_seidel_m(size_t* IA, size_t* JA, double* AN, double* b, size_t n, double f, double eps, int max_iter, double* x,
  double *work)
{
  size_t  i,j;
  int    iter;
  double  norm_b  = 0;
  double  *px    = x;
  double  *pb    = b;
  double  *AD = work;
  double  *pAD  = AD;
  size_t  *pIA  = IA;
  size_t  *pJA;
  double  *pAN;
  double  e,u;

  /* заполним массив диагональными элементами */
  for(i = 0; i < n; i++)
  {
    pJA    = JA + *pIA;
    for(j = *pIA++; j < *pIA; j++)
    {
      if(*pJA++ == i)
      {
        *pAD++ = AN[j];
        goto label;
      }
    }
    /*nl_error(nl_err_diag_elem_must_be_non_zero, 0);*/
    return -1;
label:;
  }

  /* начальное приближение к решению */
  pAD = AD;
  for(i = 0; i < n; i++)
  {
    norm_b += *pb * *pb;
    *px ++ = *pb++ / *pAD++;
  }
  eps *= eps * norm_b;

  iter = 0;
  do
  {
    pAD = AD;
    pIA = IA;
    pb  = b;
    pJA  = JA;
    px  = x;
    pAN  = AN;
    e  = 0;
    for(i = 0; i < n; i++)
    {
      u  = *pb++;
      for(j = *pIA++; j < *pIA; j++)
      {
        u -= *pAN++ * x[*pJA++];
      }
      u = u/(*pAD++);
      e += u*u;
      *px++ += f*u;
    }
    iter++;
  } 
  while(iter < max_iter && e >= eps);

  return iter;
}





int sp_conj(size_t *IA, size_t *JA, double *AN, double *b, size_t n, double eps, int max_iter, double *x, double *work)
{
  double  *pr, *p_end, *pb, *pp, *pq, *px;
  double  rho, alpha, beta, rho2, eps_ = 0.;
  double  *r;
  double  *p;
  double  *q;
  int    itr = 0;

  r = work;
  p = work + n;
  q = work + 2*n;

  /* r = A*x */
  sp_mult_col(IA, JA, AN, x, n, r);
  /* r = b - A*x */
  p_end = r + n;
  pr  = r;
  pb  = b;
  while(pr < p_end)
  {
    eps_+= *pb * *pb;
    *pr  = *pb++ - *pr;
    pr++;
  }
  eps_ *= eps * eps;

  /* главный цикл */
  while(itr < max_iter)
  {
    /* rho = ||r||2 */
    rho = cblas_ddot(n, r, 1, r, 1);
    if(rho <= eps_)
    {
      break;
    }
    if (itr)
    {
      beta = rho/rho2;
      /* p = r + beta * p; */
      pp    = p;
      p_end  = p + n;
      pr    = r;
      while(pp < p_end)
      {
        *pp  = *pr++ + beta * (*pp);
        pp++;
      }
    } 
    else
    {
      cblas_dcopy(n, r, 1, p, 1);
    }

    /* q = A*p */
    sp_mult_col(IA, JA, AN, p, n, q);

    alpha  = rho / cblas_ddot(n, p, 1, q, 1);

    /* x += alpha * p */
    /* r -= alpha * q; */
    pq = q;
    px = x;
    pr = r;
    pp = p;
    p_end = x + n;
    while(px < p_end)
    {
      *px++ += alpha * *pp++;
      *pr++ -= alpha * *pq++;
    }
    rho2 = rho;
    itr++;
    }

  return itr;
}



int sp_conj_sym(size_t *IA, size_t *JA, double *AN, double *AD, double *b, size_t n, double eps, int max_iter, double *x, double *work)
{
  double  *pr,*p_end,*pb,*pp,*pq,*px;
  double  rho, alpha, beta, rho2, eps_ = 0.;
  double  *r;
  double  *p;
  double  *q;
  int    itr = 0;

  r = work;
  p = work + n;
  q = work + 2*n;

  /* r = A*x */
  sp_mult_col_sym(IA, JA, AN, AD, x, n, r);
  /* r = b - A*x */
  p_end = r + n;
  pr  = r;
  pb  = b;
  while(pr < p_end)
  {
    eps_ += *pb * *pb;
    *pr  = *pb++ - *pr;
    pr++;
  }
  eps_ *= eps * eps;

  /* главный цикл */
  while(itr < max_iter)
  {
    /* rho = ||r||2 */
    rho = cblas_ddot(n, r, 1, r, 1);
    if(rho <= eps_)
    {
      break;
    }
    if (itr)
    {
      beta = rho/rho2;
      /* p = r + beta * p; */
      pp    = p;
      p_end  = p + n;
      pr    = r;
      while(pp < p_end)
      {
        *pp  = *pr++ + beta * (*pp);
        pp++;
      }
    } else
    {
      cblas_dcopy(n, r, 1, p, 1);
    }

    /* q = A*p */
    sp_mult_col_sym(IA, JA, AN, AD, p, n, q);

    alpha  = rho / cblas_ddot(n, p, 1, q, 1);

    /* x += alpha * p */
    /* r -= alpha * q; */
    pq = q;
    px = x;
    pr = r;
    pp = p;
    p_end = x + n;
    while(px < p_end)
    {
      *px++ += alpha * *pp++;
      *pr++ -= alpha * *pq++;
    }
    rho2 = rho;
    itr++;
    }

  return itr;
}


int sp_chol_symb(size_t* IA, size_t* JA, size_t n, size_t* IU, size_t* JU, size_t size_U, size_t *xwork)
{
  size_t  *IP = xwork;
  size_t  i,j;
  size_t  *pIA = IA;
  size_t  *pJA = JA;
  size_t  *pJU = JU;
  size_t  min, last, nz = 0, nz_, max_in_i, l, je, t;

  /* засылка начальных условий */
  memset(IP, -1, sizeof(*IP)*n);
  memset(IU, -1, sizeof(*IU)*n);

  /* внешний цикл по строкам матрицы A */
  for(i = 0; i < n-1; i++)
  {
    nz_    = nz;
    min    = n;
    max_in_i= nz + n - i;

    /* копирование столбцовых индексов, соответствующих ненулевой строке i */
    /* матрицы A, из массива JA в массив JU */
    for(j = *pIA++; j < *pIA; j++)
    {
      if(++nz > size_U)
      {
        /*nl_error(nl_err_inconsistent_size, 0);*/
        return 0;
      }
      t = *pJU++ = *pJA++;
      if(t < min)
        min = t;
      IU[t] = i;
    }

    /* выделение номера последней строки, связанной со столбцом i */
    last = IP[i];
    if(last != -1)
    {
      l = last;
      do
      {
        l = IP[l];
        je = (l+1 == i)? nz_ : IU[l+1];
        IU[i] = i;

        /* обработка строки l матрицы U (копирование индексов, */
        /* если такого индекса еще не было) */
        for(j = IU[l]; j < je; j++)
        {
          t = JU[j];
          if(IU[t] != i)
          {
            if(++nz > size_U)
            {
              /*nl_error(nl_err_inconsistent_size, 0);*/
              return 0;
            }
            *pJU++  = t;
            if(t < min)
              min = t;
            IU[t]  = i;
          }
        }
      } while(nz < max_in_i && l != last);
    }
    /* если строка оказалась не пустой */
    if(min != n)
    {
      l = IP[min];
      if(l != -1)
      {
        IP[i]  = IP[l];
        IP[l]  = i;
      } else
      {
        IP[i]  = IP[min] = i;
      }
    }
    IU[i] = nz_;
  }
  IU[n] = IU[n - 1] = nz;

  return 1;
}


void sp_chol_num(size_t* IA, size_t* JA, double* AN, double* AD, size_t* IU, size_t* JU, size_t n, 
double* UN, double* DINV, size_t *xwork)
{
  size_t  *IP  = xwork;
  size_t  *IUP = xwork + n;
  size_t  i,j;
  size_t  *pIU = IU;
  size_t  *pJU;
  double  *pD   = DINV;
  double  *pUN;
  double  UM;
  size_t  zU, last, ln, l, IUC, IUD, jj;

  /* засылка начальных условий */
  memset(IP, -1, sizeof(*IP)*n);


  /* главный цикл обработки строк матриц A и U */
  for(i = 0; i < n; i++)
  {
    zU = *pIU++;

    if(zU < *pIU)
    {
      /* засылка нулей в компоненты массива DINV, которые */
      /* соответствуют ненулевым элементам строки i матрицы U. */
      pJU = JU + zU;
      for(j = zU; j < *pIU; j++)
        DINV[*pJU++] = 0.;

      for(j = *IA++;j < *IA; j++)
        DINV[*JA++] = *AN++;
    }
    *pD = *AD++;

    /* выявление строки матрицы A, имеющей ненулевой элемент */
    /* в столбце i */
    last  = IP[i];
    if(last != -1)
    {
      ln = IP[last];
      do
      {
        l  = ln;
        ln  = IP[l];
        IUC = IUP[l];
        IUD = IU[l+1];

        /* умножение соответствующих ненулевых элементов строки */
        /* l матрицы A на коэффициент и вычитаение */
        UM  = UN[IUC]*DINV[l];
        pJU = JU + IUC;
        pUN = UN + IUC;
        for(j = IUC; j < IUD; j++)
          DINV[*pJU++] -= *pUN++ * UM;

        UN[IUC] = UM;
        IUP[l]  = ++IUC;

        if(IUC != IUD)
        {
          j = JU[IUC];
          jj = IP[j];
          if(jj == -1)
          {
            IP[l] = IP[j] = l;
          } else
          {
            IP[l] = IP[jj];
            IP[jj]= l;
          }
        }
      } while(l != last);
    }

    /* все строки, связанные со столбцом i, уже обработаны */
    *pD = 1./ (*pD);
    pD++;

    if(zU < *pIU)
    {
      pJU = JU + zU;
      pUN = UN + zU;
      for(j = zU; j < *pIU; j++)
        *pUN++ = DINV[*pJU++];

      j = JU[zU];
      jj = IP[j];
      if(jj != -1)
      {
        IP[i] = IP[jj];
        IP[jj]= i;
      } else
      {
        IP[j] = IP[i] = i;
      }
    }
    IUP[i] = zU;

  }

}



void sp_chol_solve(size_t* IU, size_t* JU, double* UN, double* DINV, double* b, size_t n, double* x)
{
  size_t  i, j;
  double *px = x;

  double *pDINV = DINV;

  size_t *pIU = IU;
  size_t *pJU = JU;
  size_t *p_end;
  double *pUN = UN;
  double xx;

  /* засылка начальных значений в промежуточный вектор Z */
  memcpy(x, b, sizeof(*x)*n);

  /* цикл вычисления компонент промежуточных векторов Z и W (прямой ход) */
  for(i = 0; i < n-1; i++)
  {
    xx = *px;
    for(j = *pIU++; j < *pIU; j++)
      x[*pJU++] -= *pUN++ * xx;

    *px++ *=  *pDINV++;
  }
  *px *= *pDINV;

  /* цикл вычисления компонент вектора решения X (обратная подстановка) */
  pJU  = JU + IU[n-1]-1;
  pUN = UN + IU[n-1]-1;
  px  = x  + n-2;
  for(i = n-1; i != 0; i--)
  {
    xx = *px;
    p_end = JU + IU[i-1];
    while(pJU >= p_end)
    {
      xx -= *pUN-- * x[*pJU--];
    }
    *px--  = xx;
  }
}



void sp_row_nz_sym(size_t* IA, size_t* JA, size_t n, size_t* nz)
{
  size_t i,z;
  size_t *pIA = IA;

  memset(nz, 0, n*sizeof(*nz));

  for(i = 0; i < n; i++)
  {
    z  = *pIA++;
    nz[i] += *pIA - z + 1; /* число ненулевых элементов справа от диагонали */
    for(; z < *pIA; z++)
      nz[*JA++]++;
  }
}


void sp_colperm_sym(size_t* IA, size_t* JA, size_t n, size_t* p, size_t *xwork)
{
  size_t i;
  size_t* ps = xwork;
  size_t* s = ps;
  sp_row_nz_sym(IA, JA, n, p);

  for(i = 0; i < n; i++)
    s[i] = i;

  nl_xvector_qsort_x(p, s, 0, n-1);

  for(i = 0; i < n; i++)
    p[*ps++] = i;
}

void sp_sym_to_complete(size_t *IS, size_t *JS, double* SN, double* SD, size_t n, size_t *IA, size_t *JA, double *AN)
{
  size_t i,j,c,ind;
  /* size_t nz = 2*sp_nz(n,IS) + n; */
  size_t  *pIA,*pIS,*pJS;
  double  *pSN;

  memset(IA, 0, sizeof(*IA)*(n+1));

  pIA = IA+2;
  pIS = IS;
  pJS = JS;
  for(i = 0; i < n-1; i++)
  {
    *pIA++ += *(pIS+1) - *pIS + 1;
    for(j = *pIS++; j < *pIS; j++)
    {
      if(*pJS < n - 1)
        IA[*pJS + 2]++;
      pJS++;
    }
  }
  if(n >= 2)
  {
    pIA = IA+2;
    c = *pIA;
    for(i = 2; i < n; i++)
    {
      c  += *(++pIA);
      *pIA = c;
    }
  }

  pIS = IS;
  pJS = JS;
  pIA = IA + 1;
  pSN = SN;
  for(i = 0; i < n; i++)
  {
    JA[*pIA] = i;
    AN[*pIA]  = *SD++;
    (*pIA)++;
    for(j = *pIS++; j < *pIS; j++, pJS++, pSN++)
    {
      JA[*pIA] = *pJS;
      AN[*pIA]  = *pSN;
      (*pIA)++;

      ind = IA[*pJS+1];
      JA[ind] = i;
      AN[ind] = *pSN;
      IA[*pJS+1]++;
    }

    pIA++;
  }
}

// Быстрая сортировка

#define MAX_THRESH 16 /* число подбирается опытном путем. 
  может улучшить сортировку в разы. на Duron 750 было получено 1024*/

typedef struct 
{
  size_t lo;
  size_t hi;
} stack_node;

/* эмуляция стека */
#define STACK_SIZE 32
#define PUSH(LOW,HIGH) do {top->lo = LOW;top++->hi = HIGH;} while (0)
#define POP(LOW,HIGH)  do {LOW = (--top)->lo;HIGH = top->hi;} while (0)
#define STACK_NOT_EMPTY (stack < top)

#define SWAP_DOUBLE(a,b) {swap_double = x[a]; x[a] = x[b]; x[b] = swap_double; \
            swap_int = v[a]; v[a] = v[b]; v[b] = swap_int; }

#define SWAP_(a,b) {swap_int = v[a]; v[a] = v[b]; v[b] = swap_int; }

#define SWAP_INDEX(a,b) {swap_index = x[a]; x[a] = x[b]; x[b] = swap_index; \
            swap_int = v[a]; v[a] = v[b]; v[b] = swap_int; }


void nl_xvector_qsort_d(size_t *v, double *x, size_t first, size_t last)
{
  size_t pivot_buffer;
  size_t  swap_int;
  double  swap_double;
  
  if (last - first > MAX_THRESH) 
  {
    size_t lo  = first;
    size_t hi  = last;
    
    stack_node stack[STACK_SIZE];
    stack_node *top = stack + 1;
    
printf("I'm here!\n");

    while (STACK_NOT_EMPTY) 
    {
      size_t *left_ptr;
      size_t *right_ptr;

      size_t mid = (hi + lo) / 2;
printf("11\n");
      pivot_buffer = v[mid];
      
      if(pivot_buffer < v[lo])
      {
printf("12\n");
        SWAP_DOUBLE(mid, lo);
        pivot_buffer = v[lo];
      }
      if(v[hi] < pivot_buffer) 
      {
printf("13\n");
        SWAP_DOUBLE(hi, mid);
        pivot_buffer = v[hi];

        if(pivot_buffer < v[lo])
        {
printf("14\n");
          SWAP_DOUBLE(mid, lo);
        }
      } 
      
      left_ptr  = v + lo + 1;
      right_ptr = v + hi - 1;
    
printf("15\n");
      do 
      {
printf("16\n");
        while(*left_ptr < pivot_buffer)
          left_ptr++;
      
printf("16.1\n");
        while(pivot_buffer < *right_ptr)
          right_ptr--;
      
printf("16.2\n");
        if (left_ptr < right_ptr) 
        {
          size_t q = left_ptr - v;
          size_t w = right_ptr - v;
printf("17\n");
          SWAP_DOUBLE(q, w);
          left_ptr++;
          right_ptr--;
        } else
        {
          if (left_ptr == right_ptr) 
          {
            left_ptr++;
            right_ptr--;
            break;
          }
        }
      } while (left_ptr <= right_ptr);
      
      if (((right_ptr - v)- lo) <= MAX_THRESH) 
      {
printf("18\n");
        if ((hi - (left_ptr - v)) <= MAX_THRESH)
          POP (lo, hi); 
        else
          lo = left_ptr - v;
      } else 
        if ((hi - (left_ptr - v)) <= MAX_THRESH)
          hi = right_ptr - v;
        else 
          if (((right_ptr - v) - lo) > (hi - (left_ptr - v)))
          {
            PUSH (lo, right_ptr - v);
            lo = left_ptr - v;
          } else 
          {
            PUSH (left_ptr - v, hi);
            hi = right_ptr - v;
          }
printf("19\n");
    }
  }
  {
    size_t end_ptr = last;
    size_t run_ptr;
    size_t tmp_ptr = first;
    size_t thresh = NL_MIN (end_ptr, first + MAX_THRESH);
    
printf("2 I'm here!\n");

    for (run_ptr = tmp_ptr + 1; run_ptr <= thresh; run_ptr++)
      if(v[run_ptr] < v[tmp_ptr])
        tmp_ptr = run_ptr;
    
    if(tmp_ptr != first) 
    {
      SWAP_DOUBLE(tmp_ptr, first);
    }
    
    for (run_ptr = first + 1; (tmp_ptr = run_ptr += 1) <= end_ptr;) 
    {
      while(v[run_ptr] < v[--tmp_ptr]);
      
      if ((++tmp_ptr) != run_ptr) 
      {
        size_t trav;    
        for (trav = run_ptr + 1; --trav >= run_ptr;) 
        {
          size_t  c;
          double d;
          size_t hi, lo;
          c = v[trav];
          d = x[trav];
      
          hi = trav;
          for (lo = trav-1; lo >= tmp_ptr; ) 
          {
            v[hi] = v[lo];
            x[hi--] = x[lo--];
          }
          v[hi] = c;
          x[hi] = d;
        }

      }
    }
  }
}


void nl_xvector_qsort_x(size_t *v, size_t *x, size_t first, size_t last)
{
  size_t pivot_buffer;
  size_t swap_int;
  size_t swap_index;
  
  if (last - first > MAX_THRESH) 
  {
    size_t lo  = first;
    size_t hi  = last;
    
    stack_node stack[STACK_SIZE];
    stack_node *top = stack + 1;
    
    while (STACK_NOT_EMPTY) 
    {
      size_t *left_ptr;
      size_t *right_ptr;

      size_t mid = (hi + lo) / 2;
      pivot_buffer = v[mid];
      
      if(pivot_buffer < v[lo])
      {
        SWAP_INDEX(mid, lo);
        pivot_buffer = v[lo];
      }
      if(v[hi] < pivot_buffer) 
      {
        SWAP_INDEX(hi, mid);
        pivot_buffer = v[hi];

        if(pivot_buffer < v[lo])
        {
          SWAP_INDEX(mid, lo);
        }
      } 
      
      left_ptr  = v + lo + 1;
      right_ptr = v + hi - 1;
    
      do 
      {
        while(*left_ptr < pivot_buffer)
          left_ptr++;
      
        while(pivot_buffer < *right_ptr)
          right_ptr--;
      
        if (left_ptr < right_ptr) 
        {
          size_t q = left_ptr - v;
          size_t w = right_ptr - v;
          SWAP_INDEX(q, w);
          left_ptr++;
          right_ptr--;
        } else
        {
          if (left_ptr == right_ptr) 
          {
            left_ptr++;
            right_ptr--;
            break;
          }
        }
      } while (left_ptr <= right_ptr);
      
      if (((right_ptr - v)- lo) <= MAX_THRESH) 
      {
        if ((hi - (left_ptr - v)) <= MAX_THRESH)
          POP (lo, hi); 
        else
          lo = left_ptr - v;
      } else 
        if ((hi - (left_ptr - v)) <= MAX_THRESH)
          hi = right_ptr - v;
        else 
          if (((right_ptr - v) - lo) > (hi - (left_ptr - v)))
          {
            PUSH (lo, right_ptr - v);
            lo = left_ptr - v;
          } else 
          {
            PUSH (left_ptr - v, hi);
            hi = right_ptr - v;
          }
    }
  }
  {
    size_t end_ptr = last;
    size_t run_ptr;
    size_t tmp_ptr = first;
    size_t thresh = NL_MIN (end_ptr, first + MAX_THRESH);
    
    for (run_ptr = tmp_ptr + 1; run_ptr <= thresh; run_ptr++)
      if(v[run_ptr] < v[tmp_ptr])
        tmp_ptr = run_ptr;
    
    if(tmp_ptr != first) 
    {
      SWAP_INDEX(tmp_ptr, first);
    }
    
    for (run_ptr = first + 1; (tmp_ptr = run_ptr += 1) <= end_ptr;) 
    {
      while(v[run_ptr] < v[--tmp_ptr]);
      
      if ((++tmp_ptr) != run_ptr) 
      {
        size_t trav;    
        for (trav = run_ptr + 1; --trav >= run_ptr;) 
        {
          size_t  c;
          size_t d;
          size_t hi, lo;
          c = v[trav];
          d = x[trav];
      
          hi = trav;
          for (lo = trav-1; lo >= tmp_ptr; ) 
          {
            v[hi] = v[lo];
            x[hi--] = x[lo--];
          }
          v[hi] = c;
          x[hi] = d;
        }

      }
    }
  }
}


void nl_xvector_qsort(size_t *v, size_t first, size_t last)
{
  size_t pivot_buffer;
  size_t  swap_int;
  
  if (last - first > MAX_THRESH) 
  {
    size_t lo  = first;
    size_t hi  = last;
    
    stack_node stack[STACK_SIZE];
    stack_node *top = stack + 1;
    
    while (STACK_NOT_EMPTY) 
    {
      size_t *left_ptr;
      size_t *right_ptr;

      size_t mid = (hi + lo) / 2;
      pivot_buffer = v[mid];
      
      if(pivot_buffer < v[lo])
      {
        SWAP_(mid, lo);
        pivot_buffer = v[lo];
      }
      if(v[hi] < pivot_buffer) 
      {
        SWAP_(hi, mid);
        pivot_buffer = v[hi];

        if(pivot_buffer < v[lo])
        {
          SWAP_(mid, lo);
        }
      } 
      
      left_ptr  = v + lo + 1;
      right_ptr = v + hi - 1;
    
      do 
      {
        while(*left_ptr < pivot_buffer)
          left_ptr++;
      
        while(pivot_buffer < *right_ptr)
          right_ptr--;
      
        if (left_ptr < right_ptr) 
        {
          size_t q = left_ptr - v;
          size_t w = right_ptr - v;
          SWAP_(q, w);
          left_ptr++;
          right_ptr--;
        } else
        {
          if (left_ptr == right_ptr) 
          {
            left_ptr++;
            right_ptr--;
            break;
          }
        }
      } while (left_ptr <= right_ptr);
      
      if (((right_ptr - v)- lo) <= MAX_THRESH) 
      {
        if ((hi - (left_ptr - v)) <= MAX_THRESH)
          POP (lo, hi); 
        else
          lo = left_ptr - v;
      } else 
        if ((hi - (left_ptr - v)) <= MAX_THRESH)
          hi = right_ptr - v;
        else 
          if (((right_ptr - v) - lo) > (hi - (left_ptr - v)))
          {
            PUSH (lo, right_ptr - v);
            lo = left_ptr - v;
          } else 
          {
            PUSH (left_ptr - v, hi);
            hi = right_ptr - v;
          }
    }
  }
  {
    size_t end_ptr = last;
    size_t run_ptr;
    size_t tmp_ptr = first;
    size_t thresh = NL_MIN (end_ptr, first + MAX_THRESH);
    
    for (run_ptr = tmp_ptr + 1; run_ptr <= thresh; run_ptr++)
      if(v[run_ptr] < v[tmp_ptr])
        tmp_ptr = run_ptr;
    
    if(tmp_ptr != first) 
    {
      SWAP_(tmp_ptr, first);
    }
    
    for (run_ptr = first + 1; (tmp_ptr = run_ptr += 1) <= end_ptr;) 
    {
      while(v[run_ptr] < v[--tmp_ptr]);
      
      if ((++tmp_ptr) != run_ptr) 
      {
        size_t trav;    
        for (trav = run_ptr + 1; --trav >= run_ptr;) 
        {
          size_t  c;
          size_t hi, lo;
          c = v[trav];
      
          hi = trav;
          for (lo = trav-1; lo >= tmp_ptr; ) 
          {
            v[hi--] = v[lo--];
          }
          v[hi] = c;
        }

      }
    }
  }
}

#undef MAX_THRESH
#undef STACK_SIZE
#undef PUSH
#undef POP
#undef STACK_NOT_EMPTY
#undef SWAP_DOUBLE
#undef SWAP_
#undef SWAP_INDEX
