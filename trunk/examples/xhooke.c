#include <stdio.h>
#include <math.h>
#include "nl.h"

/* Nonlinear Optimization using the algorithm of Hooke and Jeeves  */
/*	12 February 1994	author: Mark G. Johnson 	   */


/* Find a point X where the nonlinear function f(X) has a local    */
/* minimum.  X is an n-vector and f(X) is a scalar.  In mathe-	   */
/* matical notation  f: R^n -> R^1.  The objective function f()    */
/* is not required to be continuous.  Nor does f() need to be	   */
/* differentiable.  The program does not use or require 	   */
/* derivatives of f().						   */

/* The software user supplies three things: a subroutine that	   */
/* computes f(X), an initial "starting guess" of the minimum point */
/* X, and values for the algorithm convergence parameters.  Then   */
/* the program searches for a local minimum, beginning from the    */
/* starting guess, using the Direct Search algorithm of Hooke and  */
/* Jeeves.							   */

/* This C program is adapted from the Algol pseudocode found in    */
/* "Algorithm 178: Direct Search" by Arthur F. Kaupe Jr., Commun-  */
/* ications of the ACM, Vol 6. p.313 (June 1963).  It includes the */
/* improvements suggested by Bell and Pike (CACM v.9, p. 684, Sept */
/* 1966) and those of Tomlin and Smith, "Remark on Algorithm 178"  */
/* (CACM v.12).  The original paper, which I don't recommend as    */
/* highly as the one by A. Kaupe, is:  R. Hooke and T. A. Jeeves,  */
/* "Direct Search Solution of Numerical and Statistical Problems", */
/* Journal of the ACM, Vol. 8, April 1961, pp. 212-229. 	   */

/* Calling sequence:						   */
/*  int hooke(n, startpt, endpt, rho, epsilon, itermax)	   */
/*								   */
/*     n	   {an integer}  This is the number of dimensions  */
/*		   in the domain of f().  It is the number of	   */
/*		   coordinates of the starting point (and the	   */
/*		   minimum point.)				   */
/*     startpt	   {an array of doubles}  This is the user-	   */
/*		   supplied guess at the minimum.		   */
/*     endpt	   {an array of doubles}  This is the location of  */
/*		   the local minimum, calculated by the program    */
/*     rho	   {a double}  This is a user-supplied convergence */
/*		   parameter (more detail below), which should be  */
/*		   set to a value between 0.0 and 1.0.	Larger	   */
/*		   values of rho give greater probability of	   */
/*		   convergence on highly nonlinear functions, at a */
/*		   cost of more function evaluations.  Smaller	   */
/*		   values of rho reduces the number of evaluations */
/*		   (and the program running time), but increases   */
/*		   the risk of nonconvergence.	See below.	   */
/*     epsilon	   {a double}  This is the criterion for halting   */
/*		   the search for a minimum.  When the algorithm   */
/*		   begins to make less and less progress on each   */
/*		   iteration, it checks the halting criterion: if  */
/*		   the stepsize is below epsilon, terminate the    */
/*		   iteration and return the current best estimate  */
/*		   of the minimum.  Larger values of epsilon (such */
/*		   as 1.0e-4) give quicker running time, but a	   */
/*		   less accurate estimate of the minimum.  Smaller */
/*		   values of epsilon (such as 1.0e-7) give longer  */
/*		   running time, but a more accurate estimate of   */
/*		   the minimum. 				   */
/*     itermax	   {an integer}  A second, rarely used, halting    */
/*		   criterion.  If the algorithm uses >= itermax    */
/*		   iterations, halt.				   */


/* The user-supplied objective function f(x,n) should return a C   */
/* "double".  Its  arguments are  x -- an array of doubles, and    */
/* n -- an integer.  x is the point at which f(x) should be	   */
/* evaluated, and n is the number of coordinates of x.	That is,   */
/* n is the number of coefficients being fitted.		   */

/* rho, the algorithm convergence control			   */
/*	The algorithm works by taking "steps" from one estimate of */
/*    a minimum, to another (hopefully better) estimate.  Taking   */
/*    big steps gets to the minimum more quickly, at the risk of   */
/*    "stepping right over" an excellent point.  The stepsize is   */
/*    controlled by a user supplied parameter called rho.  At each */
/*    iteration, the stepsize is multiplied by rho  (0 < rho < 1), */
/*    so the stepsize is successively reduced.			   */
/*	Small values of rho correspond to big stepsize changes,    */
/*    which make the algorithm run more quickly.  However, there   */
/*    is a chance (especially with highly nonlinear functions)	   */
/*    that these big changes will accidentally overlook a	   */
/*    promising search vector, leading to nonconvergence.	   */
/*	Large values of rho correspond to small stepsize changes,  */
/*    which force the algorithm to carefully examine nearby points */
/*    instead of optimistically forging ahead.	This improves the  */
/*    probability of convergence.				   */
/*	The stepsize is reduced until it is equal to (or smaller   */
/*    than) epsilon.  So the number of iterations performed by	   */
/*    Hooke-Jeeves is determined by rho and epsilon:		   */
/*	    rho**(number_of_iterations) = epsilon		   */
/*	In general it is a good idea to set rho to an aggressively */
/*    small value like 0.5 (hoping for fast convergence).  Then,   */
/*    if the user suspects that the reported minimum is incorrect  */
/*    (or perhaps not accurate enough), the program can be run	   */
/*    again with a larger value of rho such as 0.85, using the	   */
/*    result of the first minimization as the starting guess to    */
/*    begin the second minimization.				   */

/* Normal use: (1) Code your function f() in the C language	   */
/*	       (2) Install your starting guess {or read it in}	   */
/*	       (3) Run the program				   */
/*	       (4) {for the skeptical}: Use the computed minimum   */
/*		      as the starting point for another run	   */

/* Data Fitting:						   */
/*	Code your function f() to be the sum of the squares of the */
/*	errors (differences) between the computed values and the   */
/*	measured values.  Then minimize f() using Hooke-Jeeves.    */
/*	EXAMPLE: you have 20 datapoints (ti, yi) and you want to   */
/*	find A,B,C such that  (A*t*t) + (B*exp(t)) + (C*tan(t))    */
/*	fits the data as closely as possible.  Then f() is just    */
/*	f(x) = SUM (measured_y[i] - ((A*t[i]*t[i]) + (B*exp(t[i])) */
/*				  + (C*tan(t[i]))))^2		   */
/*	where x[] is a 3-vector consisting of {A, B, C}.	   */

/*								   */
/*  Функция написана на основе hooke.c by M.G. Johnson             */ 
/*  Следующий абзац взят из комментариев к hooke.c                 */
/*								   */
/*  The author of this software is M.G. Johnson.		   */
/*  Permission to use, copy, modify, and distribute this software  */
/*  for any purpose without fee is hereby granted, provided that   */
/*  this entire notice is included in all copies of any software   */
/*  which is or includes a copy or modification of this software   */
/*  and in all copies of the supporting documentation for such	   */
/*  software.  THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT    */
/*  ANY EXPRESS OR IMPLIED WARRANTY.  IN PARTICULAR, NEITHER THE   */
/*  AUTHOR NOR AT&T MAKE ANY REPRESENTATION OR WARRANTY OF ANY	   */
/*  KIND CONCERNING THE MERCHANTABILITY OF THIS SOFTWARE OR ITS    */
/*  FITNESS FOR ANY PARTICULAR PURPOSE. 			   */
/*								   */



/* given a point, look for a better one nearby, one coord at a time */

void
search_neighbourhood(
   double (*fun)(double*), 
   double *x, 
   double f, 
   double *delta, 
   size_t n, 
   int *nfunevals,
   double *fnew)
{
   double fmin;
   size_t j;

   fmin = f;

   for (j = 0; j < n; j++) 
   {
      x[j] += delta[j];

      f = (*fun)(x);
      (*nfunevals)++;

      if (f < fmin)
         fmin = f;
      else 
      {
         delta[j] = -delta[j];
         x[j] += 2*delta[j];

         f = (*fun)(x);
         (*nfunevals)++;

         if (f < fmin)
            fmin = f;
         else
            x[j] -= delta[j];
      }
   }

   *fnew = fmin;
}

/* Остерегаемся слишком маленьких шагов;
   fnew < f может явиться следствием ошибок округления */

int is_step_too_small(double* x, double* xnew, double* delta, size_t n)
{
   size_t j;

   for (j = 0; j < n; j++) 
   {
      if (fabs(xnew[j] - x[j]) > (0.5 * fabs(delta[j])))
         return 0;
   }
   return 1;
}


void optim_hooke(double (*fun)(double*), size_t n, double* x0, double rho, double ftol, double xtol, 
   int maxiter, int maxfunevals, double *f, int *rc, 
   int *niter, int *nfunevals, double *work)  
{
  double fnew, h, xx;
  double *x, *xnew, *delta;

  size_t j;

  x = x0;
  xnew = work;
  delta = work + n;

  for (j = 0; j < n; j++) 
  {
     xnew[j] = x[j];

     delta[j] = fabs(x[j] * rho);

     if (delta[j] == 0)
        delta[j] = rho;
  }

  h = rho;

  *f = fnew = (*fun)(xnew);

  *nfunevals = 1;
  *niter = 0;

  while (h > xtol && *nfunevals < maxfunevals && *niter < maxiter)
  {
     (*niter)++;

     cblas_dcopy(n, x, 1, xnew, 1);

     search_neighbourhood(fun, xnew, *f, delta, n, nfunevals, &fnew);

     /* Двигаемся в заданном направлении,
        пока значение целевой функции уменьшается */

     while (fnew < *f && !is_step_too_small(x, xnew, delta, n) 
        && (*nfunevals) <= maxfunevals) 
     {
        for (j = 0; j < n; j++) 
        {
           /* вычисляем знак delta[j] */

           if (xnew[j] <= x[j])
              delta[j] = -fabs(delta[j]);
           else
              delta[j] = fabs(delta[j]);

           /* двигаемся в заданном направлении */

           xx = x[j];
           x[j] = xnew[j];
           xnew[j] += xnew[j] - xx;
        }

        *f = fnew;
        search_neighbourhood(fun, xnew, *f, delta, n, nfunevals, &fnew);
     }

     if (h >= xtol && fnew >= *f) 
     {
        h *= rho;

        for (j = 0; j < n; j++) 
     	   delta[j] *= rho;
     }
  }

  if ((*nfunevals) > maxfunevals || (*niter) > maxiter)
     *rc = 1;
  else
     *rc = 0;
}

double fun(double *x)
{

  double tmp1, tmp2;

  tmp1 = x[1] - x[0]*x[0];
  tmp2 = 1 - x[0];

  return tmp1*tmp1 + tmp2*tmp2;
}

int main(void)
{
  size_t n;
  double *x0, f, *work;
  double f0;
  double tolf, tolx;
  int maxfun, maxiter, rc, nfun, niter;

  n = 2;

  x0 = nl_dvector_create(n);
  work = nl_dvector_create(2*n);

  x0[0] = -1.2;
  x0[1] = 1;
  tolf = 1.0e-6;
  tolx = 1.0e-6;
  maxfun = 200;
  maxiter = 50;

  optim_hooke(fun, n, x0, 0.6, tolf, tolx, maxiter, maxfun, 
     &f, &rc, &niter, &nfun, work);

  if (rc)
  {
    printf("\nЧисло итераций или количество вычисленных значений функции\n");
    printf("превысило максимально допустимое!\n" );
  }
  else
    printf("\nУспешное завершение\n");

  printf("\nКоличество вычисленных значений функции: %d\n", nfun);
  printf("Число итераций: %d\n",niter);

  printf("\nВычисленная точка минимума:\n");
  nl_dvector_print(x0, n, "  %12.6e");

  printf("\nЗначения функции: %e\n", f);
  
  nl_dvector_free(x0);
  nl_dvector_free(work);

  return 0;
}


