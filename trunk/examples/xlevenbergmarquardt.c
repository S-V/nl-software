#include <stdio.h>
#include <math.h>
#include "nl.h"

double max_diag(double *A, size_t n)
{
   double x;
   size_t i;

   x = A[0];

   for (i = n + 1; i < n*n; i += n + 1)
   {
      if (A[i] > x)
         x = A[i];
   }

   return x;
}

void A_plus_diag(double *A, size_t n, double alpha)
{
   size_t i;
   double x;

   x = A[0];

   for (i = n + 1; i < n*n; i += n + 1)
   {
      A[i] += alpha;
   }
}

double squares_sum(double *f, size_t m)
{
   double x;
   size_t j;

   x = 0;

   for (j = 0; j < m; j++)
     x += f[j]*f[j];

   return x;
}


void optim_levenbergmarquardt(
   void (*fun)(double*, double*), 
   void (*jac)(double*, double*), 
   size_t n, 
   size_t m, 
   double* x0, 
   double ftol, 
   double xtol, 
   int maxiter, 
   int maxfunevals, 
   double *f, 
   double *J, 
   double *A, 
   double *g, 
   int *rc, 
   int *niter, 
   int *nfunevals, 
   double *work)  
{
   double alpha, beta, *x, *xnew, *h, *fnew, rho, mu;

   h = work;
   xnew = work + n;
   fnew = work + 2*n; //fnew - m

   *niter = 0;
   *nfunevals = 0;
   beta = 2;
   x = x0;

   (*fun)(x, f);
   (*jac)(x, J);

   cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, n, m, 
      1, J, n, J, n, 0, A, n);  /* A = J'*J */
   cblas_dgemv(CblasRowMajor, CblasTrans, m, n, 1, J, n, f, 1, 0, g, 1); /* g = J'*f */
   cblas_dscal(n, -1, g, 1);    /* g = -g */

   if (fabs(g[cblas_idamax(n, g, 1)]) < xtol)
   {
      *rc = 0; 
      return;
   }

   alpha = 0.5*max_diag(A, n);


   while(1)
   {
      if (*niter > maxiter)
      {
         *rc = -1;
         break;
      }
      
      (*niter)++;
   
      A_plus_diag(A, n, alpha);	/* A = A + alpha*I */

      chol_decomp(A, n);
      cblas_dcopy(n, g, 1, h, 1);
      chol_solve(A, n, h);	/* решение записано в h */

      if (cblas_dnrm2(n, h, 1) <= xtol*(cblas_dnrm2(n, x, 1) + xtol))
      {
         *rc = 1; 
         return;
      }

      cblas_dcopy(n, x, 1, xnew, 1);
      cblas_daxpy(n, 1, h, 1, xnew, 1);		/* xnew = x + h */

      cblas_daxpy(n, alpha, h, 1, g, 1);	/* g += alpha*h */

      (*fun)(xnew, fnew);

      rho = 2*(squares_sum(f, m) - squares_sum(fnew, m))/cblas_ddot(n, h, 1, g, 1);

      if (rho > 0)
      {
         cblas_dcopy(n, xnew, 1, x, 1);

         (*fun)(x, f);
         (*jac)(x, J);
         
         cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, n, m, 
            1, J, n, J, n, 0, A, n);  /* A = J'*J */
         cblas_dgemv(CblasRowMajor, CblasTrans, m, n, 1, J, n, f, 1, 0, g, 1); /* g = J'*f */
         cblas_dscal(n, -1, g, 1);    /* g = -g */

         if (fabs(g[cblas_idamax(n, g, 1)]) < xtol)
         {
            *rc = 1; 
            return;
         }

         mu = 2*rho - 1;
         mu = 1 - mu*mu*mu;

         alpha = alpha*NL_MAX(0.3, mu);
         beta = 2;
      }
      else
      {
         alpha *= beta;
         beta *= 2;
      }
      
   }

}

void fn(double *x, double *f)
{
  f[0] = x[1] - x[0]*x[0];
  f[1] = 1 - x[0];
}

void jc(double *x, double *J)
{
  J[0] = -2*x[0]; J[1] = 1;
  J[2] = -1;      J[3] = 0;
}

int main()
{
  size_t n, m;
  double *x0, *f, *work, *J, *A, *g;
  double f0;
  double tolf, tolx;
  int maxfunevals, maxiter, rc, nfunevals, niter;

  n = 2;
  m = 2;

  x0 = nl_dvector_create(n);
  f = nl_dvector_create(m);
  J = nl_dmatrix_create(m, n);
  A = nl_dmatrix_create(n, n);
  g = nl_dvector_create(n);
  work = nl_dvector_create(2*n + m);

  x0[0] = -1.2;
  x0[1] = 1;
  tolf = 1.0e-6;
  tolx = 1.0e-6;
  maxfunevals = 200;
  maxiter = 50;

  optim_levenbergmarquardt(
   fn,
   jc,
   n,
   m,
   x0, 
   tolf, 
   tolx, 
   maxiter, 
   maxfunevals, 
   f,
   J,
   A,
   g,
   &rc,
   &niter,
   &nfunevals,
   work);

  if (rc == 1)
    printf("\nУспешное завершение\n");
  else
  {
    printf("\nЧисло итераций или количество вычисленных значений функции\n");
    printf("превысило максимально допустимое!\n" );
  }

  printf("\nКоличество вычисленных значений функции: %d\n", nfunevals);
  printf("Число итераций: %d\n",niter);

  printf("\nВычисленная точка минимума:\n");
  nl_dvector_print(x0, n, "  %12.6e");

  printf("\n?Значения функции: %e\n", f);
  
  nl_dvector_free(x0);
  nl_dvector_free(work);

  return 0;
}







