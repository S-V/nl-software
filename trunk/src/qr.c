#include <memory.h>
#include "nl.h"

void qr_decomp(double *A, size_t m, size_t n, double* t)
{
  size_t j, i, k;
  double norm, d, sum, alpha, beta, tau;
  size_t min = NL_MIN(n,m);
  for (j = 0; j < min; j++)
  {
    /* Вычисление преобразований Хаусхолдера, чтобы привести j-й столбец
       матрицы к множителю j-го единичного вектора. */

    norm = 0.;
    for(i = j + 1; i < m; i++)
    {
      if( (d = fabs(A[i*n + j])) > norm)
        norm = d;
    }
    sum = 0.;
    for(i = j+1; i < m; i++)
    {
      d = A[i*n + j]/norm;
      sum += d*d;
    };
    norm = sqrt(sum)*norm;


    alpha = A[j*n + j];
    beta = -(alpha >= 0.0 ? +1.0 : -1.0) * hypot(alpha, norm) ;
    /*if(beta == 0.)
    {
      err_matrix_is_singular
    };*/
    t[j] = tau = (beta-alpha)/beta;
    d   = (alpha-beta);
    for(i = j+1; i < m; i++)
    {
      A[i*n + j] /= d;
    };
    A[j*n + j] = beta;
    for(k = j+1; k < n; k++)
    {
      sum = A[j*n + k];
      for(i = j+1; i < m; i++)
        sum += A[i*n + k]*A[i*n + j];
      A[j*n + k] -= tau*sum;
      for(i = j+1; i < m; i++)
        A[i*n + k] -= tau*sum*A[i*n + j];
    }
  }
}

void GetQk(size_t k, double *QR, size_t n, size_t m, double* t, double *Q_k)
{
  size_t i,j;
  double tk = t[k];
  for(j = 0; j < k; j++)
  {
    for(i = 0; i < k; i++)
    {
      Q_k[i*n + j] = Q_k[j*n + i] = 0.;
    };
    Q_k[j*n + j] = 1.;
    for(i = k; i < m; i++)
    {
      Q_k[i*n + j] = Q_k[j*n + i] = 0.;
    };
  }
  Q_k[k*n + k] = 1-tk;
  for(i = k+1; i < m; i++)
  {
    Q_k[k*n + i] = Q_k[i*n + k] = -QR[i*n + k]*t[k];
  };
  for(j = k+1; j < m; j++)
  {
    for(i = k+1; i < j; i++)
      Q_k[i*n + j] = Q_k[j*n + i] = -QR[i*n + k]*QR[j*n + k]*t[k];
    Q_k[j*n + j] = 1 - QR[i*n + k]*QR[j*n + k]*t[k];
  }
}

void qr_unpack(double *QR, size_t m, size_t n, double* t, double *Q, double *R, double *work)
{
  size_t i,j,k;
  double *Q_k, *M;
  double min;

  for(j = 0; j < m; j++)
  {
    min = NL_MIN(j,n);
    for(i = 0; i < min; i++)
      R[j*n + i] = 0;
    for(i = j; i < n; i++)
      R[j*n + i] = QR[j*n + i];
  };

  Q_k = work;
  M  = work + m*m;
  GetQk(0, QR, n, m, t, Q);

  for(k = 1; k < n; k++)
  {
    GetQk(k, QR, n, m, t, Q_k);
    for (j = 0; j < m; j++)
      for (i = 0; i < m; i++)
        M[i*m + j] = 0;

    //nl_dmatrix_mult(m,m,m,Q,Q_k,M);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, m, m, 1, Q, n, Q_k, m, 1, M, m);

    for(j = 0; j < m; j++)
      for(i = 0; i < m; i++)
        Q[j*n + i] = M[j*n + i];
  }
}

void qr_solve(double *QR, size_t n, double* t, double* b)
{
  size_t i,j,k;
  double tau,d;
  // Q'b
  for(k = 0; k < n; k++)
  {

    if((tau = t[k]) != 0.)
    {
      /* d = v'b */
      d = b[k];
      for (i = k+1; i < n; i++)
        d += QR[i*n + k] * b[i];

      /* b = b - tau (v) (v'b) */
      d  *= tau;
      b[k] -= d;
      for (i = k+1; i < n; i++)
      {
        b[i] -= QR[i*n + k]*d;
      }
    }
  }
  // Solve Rx = b
  for(i = n-1; i >= 1; i--)
  {
    d = b[i] /= QR[i*n + i];
    for(j = 0; j < i; j++)
      b[j] -= QR[j*n + i]*d;
  };
  b[0] /= QR[0]; /* QR(0, 0) */
}

void qr_qt_mul_b(double *QR, size_t m, size_t n, double* t, double* b)
{
  size_t i,k;
  double tau,d;
  // b = Q'*b
  for(k = 0; k < n; k++)
  {
    if((tau = t[k]) != 0.)
    {
      /* d = v'b */
      d = b[k];
      for (i = k+1; i < m; i++)
        d += QR[i*n + k] * b[i];

      /* b = b - tau (v) (v'b) */
      b[k] -= tau*d;
      for (i = k+1; i < m; i++)
      {
        b[i] -= tau*QR[i*n + k]*d;
      }
    }
  }
}


void qr_q_mul_b(double *QR, size_t m, size_t n, double* t, double* b)
{
  size_t i,k;
  double tau,d;
  // b = Q*b
  for(k = n-1; ; k--)
  {
    if((tau = t[k]) != 0.)
    {
      /* d = v'b */
      d = b[k];
      for (i = k+1; i < m; i++)
        d += QR[i*n + k] * b[i];

      /* b = b - tau (v) (v'b) */
      b[k] -= tau*d;
      for (i = k+1; i < m; i++)
      {
        b[i] -= tau*QR[i*n + k]*d;
      }
    }
    if(k == 0) break;
  }
}

void qr_least_squares(double *QR, size_t m, size_t n, double* t, double* b, double* r, double *work)
{
  size_t i,j;
  double d;
  double* bb = work;

  cblas_dcopy(m, b, 1, bb, 1);

  // b = Q'b
  qr_qt_mul_b(QR, m,n,t,b);

  // Solve Rx = b
  for(i = n-1; i >= 1; i--)
  {
    d = b[i] /= QR[i*n + i];
    for(j = 0; j < i; j++)
      b[j] -= QR[j*n + i]*d;
  };
  b[0] /= QR[0]; /* QR(0, 0) */

  // r = R*b
  for(j = 0; j < n; j++)
  {
    d = 0;
    for(i = j; i < n; i++)
      d += QR[j*n + i]*b[i];
    r[j] = d;
  }
  for(j = n; j < m; j++)
    r[j] = 0.;
  // r = Q*R*b
  qr_q_mul_b(QR,m,n,t,r);
  //nl_dvector_sub(r, bb, m);
  cblas_daxpy(m, -1, bb, 1, r, 1);
}

