#include <math.h>
#include "nl.h"

#define QHULL_SWAP(i, j) tmp = k[i]; k[i] = k[j]; k[j] = tmp;
#define QHULL_STACK_INIT stack_pos = 0;
#define QHULL_STACK_IS_NOT_EMPTY stack_pos
#define QHULL_STACK_PUSH(p, np, r, l) stack[stack_pos++] = p;  \
                                      stack[stack_pos++] = np; \
                                      stack[stack_pos++] = r;  \
                                      stack[stack_pos++] = l;
#define QHULL_STACK_POP(p, np, r, l)  l  = stack[--stack_pos]; \
                                      r  = stack[--stack_pos]; \
                                      np = stack[--stack_pos]; \
                                      p  = stack[--stack_pos];

static void find_points(double* x, double* y, size_t* k, size_t n, 
  size_t r, size_t v, size_t l, 
  size_t cur, size_t ncur, size_t pv, 
  size_t* pr, size_t* npr, size_t* pl, size_t* npl)
{
/*
  find points "above" the line r-v and "above" the line v-l
  allocate points above r-v to the begin of the array p==pr 
  then allocate v
  then allocate points above v-l
*/
  size_t pp, ppe; /* iterators */
  double a, b, c;
  size_t tmp, j;

  QHULL_SWAP(pv, cur + ncur - 1) /* swap v to the end of the array */

  a = y[v] - y[r];
  b = x[r] - x[v];
  c = a*x[v] + b*y[v];

  *pr = pp = cur;
  ppe = cur + ncur - 1;
  (*npr) = 0;

  while (pp != ppe)
  {
    j = k[pp];
    if (a*x[j] + b*y[j] > c)
    {
      pp++;
      (*npr)++;
    }
    else
    {
      ppe--;
      QHULL_SWAP(pp, ppe)
    }
  }

  QHULL_SWAP(pp, cur + ncur - 1) /* swap v to the middle points above r-v and v-l */

  a = y[l] - y[v];
  b = x[v] - x[l];
  c = a*x[v] + b*y[v];

  pp++;
  *pl = pp;
  ppe = cur + ncur;
  (*npl) = 0;
      
  while (pp != ppe)
  {
    j = k[pp];
    if (a*x[j] + b*y[j] > c)
    {
      pp++;
      (*npl)++;
    }
    else
    {
      ppe--;
      k[pp] = k[ppe];
      k[ppe] = n; 
    }
  }
}


void conv_qhull(double* x, double* y, size_t n, size_t* k, size_t* nk, size_t* stack)
{
  size_t cur, ncur, r, v, l, pr, pl, npr, npl, pv, j, jj, tmp, stack_pos;
  double a, b, f, fmax;

  QHULL_STACK_INIT

  if (n == 0)
  {
    *nk = 0;
    return;
  }

  if (n == 1)
  {
    *nk = 1;
    *k = 0;
    return;
  }

  r = 0;
  for (j = 1; j < n; j++)
  {
    if (x[j] > x[r] || (x[j] == x[r] && y[j] > y[r]))
      r = j;
  }
 
  l = 0;
  for (j = 1; j < n; j++)
  {
    if (x[j] < x[l] || (x[j] == x[l] && y[j] < y[l]))
      l = j;
  }
 
  if (r == l)
  {
    *nk = 1;
    *k = r;
    return;
  }

  for (j = 0; j < n; j++)
    k[j] = j;

  k[0] = r; /* swap r to the begin of the array */
  k[r] = 0;

  if (l) 
    pv = l; /* position of l */
  else
    pv = r;

  find_points(x, y, k, n, r, l, r, 1, n - 1, pv, &pr, &npr, &pl, &npl);

  if (npr < npl)
  {
    QHULL_STACK_PUSH(pl, npl, l, r)
    cur = pr;
    ncur = npr;
  }
  else
  {
    QHULL_STACK_PUSH(pr, npr, r, l)
    cur = pl;
    ncur = npl;
    tmp = r; r = l; l = tmp;
  }

  while (1)
  {
    /* Consider all intervals in the stack */

    while (ncur)
    {
      /* Consider current interval */

      a = y[l] - y[r];
      b = x[r] - x[l];
      /*
        Among all points maximizing a*x + b*y we choose the dominated point v
        i.e. x = max, y = max
      */
      
      pv = cur;
      v = k[cur];
      fmax = a*x[v] + b*y[v];

      for (j = cur + 1; j < cur + ncur; j++)
      {
        jj = k[j];
        f = a*x[jj] + b*y[jj];
        if (f > fmax || (f == fmax && 
            (x[jj] > x[v] || (x[jj] == x[v] && y[jj] > y[v]))))
        {
          pv = j;
          v = jj;
          fmax = f;
        }
      }

      find_points(x, y, k, n, r, v, l, cur, ncur, pv, &pr, &npr, &pl, &npl);

      if (npr < npl)
      {
        QHULL_STACK_PUSH(pl, npl, v, l)
        cur = pr;
        ncur = npr;
        l = v;
      }
      else
      {
        QHULL_STACK_PUSH(pr, npr, r, v)
        cur = pl;
        ncur = npl;
        r = v;
      }

    }

    if (QHULL_STACK_IS_NOT_EMPTY)
    {
      QHULL_STACK_POP(cur, ncur, r, l)
    }
    else
      break;
  }


  *nk = 0;
  for (j = 0; j < n; j++)
    if (k[j] != n)
    {
      k[*nk] = k[j];
      (*nk)++;
    }

}


#undef QHULL_SWAP
#undef QHULL_STACK_INIT
#undef QHULL_STACK_IS_NOT_EMPTY
#undef QHULL_STACK_PUSH
#undef QHULL_STACK_POP

#define TORIGHT(a, b, c) (p[b].x - p[a].x)*(p[c].y - p[b].y) <= \
              (p[b].y - p[a].y)*(p[c].x - p[b].x)

static int compare_above(const void *tpi, const void *tpj)
{
  conv_point *pi, *pj;

  pi = (conv_point*)tpi;
  pj = (conv_point*)tpj;

  if (pi->x > pj->x) 
    return -1;
  if (pi->x < pj->x) 
    return 1;
  if (pi->y > pj->y) 
    return -1;
  if (pi->y < pj->y) 
    return 1;
  return 0;  
}

static int compare_under(const void *tpi, const void *tpj)
{
  conv_point *pi, *pj;

  pi = (conv_point*)tpi;
  pj = (conv_point*)tpj;

  if (pi->x < pj->x) 
    return -1;
  if (pi->x > pj->x) 
    return 1;
  if (pi->y < pj->y) 
    return -1;
  if (pi->y > pj->y) 
    return 1;
  return 0;  
}

void conv_andrew(const double* x, const double* y, size_t n, conv_point* p, size_t* nk)
{
  double xl, xr, yl, yr, a, b, c, cj;
  size_t j;
  size_t above, under;
 
  if (n == 0)
  {
    *nk = 0;
    return;
  }

  if (n == 1)
  {
    *nk = 1;
    p[0].x = x[0];
    p[0].y = y[0];
    p[0].k = 0;
    return;
  }

  xl = x[0];
  yl = y[0];
  xr = x[0];
  yr = y[0];

  for(j = 1; j < n; j++)
  {
    if (x[j] < xl || (x[j] == xl && y[j] < yl))
    {
      xl = x[j];
      yl = y[j];
    }
    if (x[j] > xr || (x[j] == xr && y[j] > yr))
    {
      xr = x[j];
      yr = y[j];
    }
  }

  a = yl - yr;
  b = xr - xl;
  c = a*xr + b*yr;

  above = 0;
  under = n;

  for (j = 0; j < n; j++)
  {
    cj = a*x[j] + b*y[j]; 
      /* It is important to assign this to some variable,
         otherwise one must insert in the following
         also the check  || j == r  */
    if (cj >= c)  
    {
      p[above].x = x[j];
      p[above].y = y[j];
      p[above].k = j;
      above++;
    }
    else 
    {
      under--;
      p[under].x = x[j];
      p[under].y = y[j];
      p[under].k = j;
    }
  }


  qsort(p, above, sizeof(conv_point), compare_above);

  qsort(p + above, n - above, sizeof(conv_point), compare_under);


  *nk = 2;

  for (j = 2; j < n; j++)
  {
      while (*nk >= 2 && TORIGHT(*nk - 2, *nk - 1, j))
        (*nk)--;

      p[*nk] = p[j];
      (*nk)++;
  }
  
  if (*nk == 2)
  {
    if (p[0].x == p[1].x && p[0].y == p[1].y)
      *nk = 1;
  }
  else if (*nk > 2)
  {
    while (*nk >= 2 && TORIGHT(*nk - 2, *nk - 1, 0))
      (*nk)--;
  }

}

#undef TORIGHT
