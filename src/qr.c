#include <memory.h>
#include "util.h"
#include "qr.h"

#define QR_MIN(a,b) (a < b)? a:b

void qr_decomp(double** A, size_t m, size_t n, double* t)
{
  size_t j, i, k;
  double norm, d, sum, alpha, beta, tau;
  size_t min = QR_MIN(n,m);
  for (j = 0; j < min; j++)
  {
    /* Вычисление преобразований Хаусхолдера, чтобы привести j-ый столбец
       матрицы к множителю j-ого единичного вектора. */

    norm = 0.;
    for(i = j+1; i < m; i++)
    {
      if( (d = fabs(A[i][j])) > norm)
        norm = d;
    }
    sum = 0.;
    for(i = j+1; i < m; i++)
    {
      d = A[i][j]/norm;
      sum += d*d;
    };
    norm = sqrt(sum)*norm;


    alpha = A[j][j];
    beta = -(alpha >= 0.0 ? +1.0 : -1.0) * hypot(alpha, norm) ;
    if(beta == 0.)
    {
      nl_error(nl_err_matrix_is_singular, 1);
    };
    t[j] = tau = (beta-alpha)/beta;
    d   = (alpha-beta);
    for(i = j+1; i < m; i++)
    {
      A[i][j] /= d;
    };
    A[j][j] = beta;
    for(k = j+1; k < n; k++)
    {
      sum = A[j][k];
      for(i = j+1; i < m; i++)
        sum += A[i][k]*A[i][j];
      A[j][k] -= tau*sum;
      for(i = j+1; i < m; i++)
        A[i][k] -= tau*sum*A[i][j];
    }
    }
}

void qr_decomp_t(double** At, size_t m, size_t n, double* t)
{
  size_t j,k;
  double norm,d,sum,alpha,beta,tau;
  double *pAt, *pAtj;
  double *pd,  *p_end;
  size_t min = QR_MIN(n,m);

  for (j = 0; j < min; j++)
  {
    /* Вычисление преобразований Хаусхолдера, чтобы привести j-ый столбец
       матрицы к множителю j-ого единичного вектора. */
    norm = 0.;
    pAt = pAtj = At[j] + j+1;
    p_end = At[j] + n;
    while(pAt < p_end)
    {
      if( (d = fabs(*pAt++)) > norm)
        norm = d;
    }
    sum = 0.;
    pAt = pAtj;
    while(pAt < p_end)
    {
      d = *pAt++/norm;
      sum += d*d;
    };
    norm = sqrt(sum)*norm;

    alpha = At[j][j];
    beta = -(alpha >= 0.0 ? +1.0 : -1.0) * hypot(alpha, norm) ;
    if(beta == 0.)
    {
      nl_error(nl_err_matrix_is_singular, 1);
    };
    t[j] = tau = (beta-alpha)/beta;
    d   = (alpha-beta);
    pAt = pAtj;
    while(pAt < p_end)
    {
      *pAt++ /= d;
    };
    At[j][j] = beta;
    for(k = j+1; k < m; k++)
    {
      sum = At[k][j];
      pd  = At[k] + j+1;
      pAt  = pAtj;
      while(pAt < p_end)
        sum += (*pd++)*(*pAt++);
      At[k][j] -= tau*sum;
      pd  = At[k] + j+1;
      pAt = pAtj;
      while(pAt < p_end)
        *pd++ -= tau*sum*(*pAt++);

    }
    }
}



void GetQk(size_t k, double** QR, size_t m, double* t, double** Q_k)
{
  size_t i,j;
  double tk = t[k];
  for(j = 0; j < k; j++)
  {
    for(i = 0; i < k; i++)
    {
      Q_k[i][j] = Q_k[j][i] = 0.;
    };
    Q_k[j][j] = 1.;
    for(i = k; i < m; i++)
    {
      Q_k[i][j] = Q_k[j][i] = 0.;
    };
  }
  Q_k[k][k] = 1-tk;
  for(i = k+1; i < m; i++)
  {
    Q_k[k][i] = Q_k[i][k] = -QR[i][k]*t[k];
  };
  for(j = k+1; j < m; j++)
  {
    for(i = k+1; i < j; i++)
      Q_k[i][j] = Q_k[j][i] = -QR[i][k]*QR[j][k]*t[k];
    Q_k[j][j] = 1 - QR[i][k]*QR[j][k]*t[k];
  }

//  nl_dmatrix_print(Q_k, m,m, "%5.2f\t"); printf("\n\n");
}


void GetQk_t(size_t k, double** QRt, size_t m, double* t, double** Q_k)
{
  size_t i,j;
  double tk = t[k];
  for(j = 0; j < k; j++)
  {
    for(i = 0; i < k; i++)
    {
      Q_k[j][i] = Q_k[i][j] = 0.;
    };
    Q_k[j][j] = 1.;
    for(i = k; i < m; i++)
    {
      Q_k[i][j] = Q_k[j][i] = 0.;
    };
  }
  Q_k[k][k] = 1-tk;
  for(i = k+1; i < m; i++)
  {
    Q_k[k][i] = Q_k[i][k] = -QRt[k][i]*t[k];
  };
  for(j = k+1; j < m; j++)
  {
    for(i = k+1; i < j; i++)
      Q_k[i][j] = Q_k[j][i] = -QRt[k][i]*QRt[k][j]*t[k];
    Q_k[j][j] = 1 - QRt[k][i]*QRt[k][j]*t[k];
  }
//  nl_dmatrix_print(Q_k, m,m, "%5.2f\t"); printf("\n\n");
}



void qr_least_unpack(double** QR, size_t m, size_t n, double* t, double** Q, double** R)
{
  size_t i,j,k;
  double** Q_k, **M;
  double min;

  for(j = 0; j < m; j++)
  {
    min = QR_MIN(j,n);
    for(i = 0; i < min; i++)
      R[j][i] = 0;
    for(i = j; i < n; i++)
      R[j][i] = QR[j][i];
  };

  Q_k = nl_dmatrix_create(m,m);
  M  = nl_dmatrix_create(m,m);
  GetQk(0, QR, m, t, Q);

  for(k = 1; k < n; k++)
  {
    GetQk(k, QR, m, t, Q_k);
    nl_dmatrix_mult(m,m,m,Q,Q_k,M);

    for(j = 0; j < m; j++)
      for(i = 0; i < m; i++)
        Q[j][i] = M[j][i];
  }
}




void qr_least_unpack_t(double** QRt, size_t m, size_t n, double* t, double** Q, double** R)
{
  size_t i,j,k;
  double** Q_k, **M;
  double min;

  for(j = 0; j < n; j++)
  {
    min = QR_MIN(j,m);
    for(i = 0; i < min; i++)
      R[j][i] = 0;
    for(i = j; i < m; i++)
      R[j][i] = QRt[i][j];
  };

  Q_k = nl_dmatrix_create(n,n);
  M  = nl_dmatrix_create(n,n);
  GetQk_t(0, QRt, n, t, Q);

  for(k = 1; k < m; k++)
  {
    GetQk_t(k, QRt, n, t, Q_k);
    nl_dmatrix_mult(n,n,n,Q,Q_k,M);

    for(j = 0; j < n; j++)
      for(i = 0; i < n; i++)
        Q[j][i] = M[j][i];
  }
}




void qr_solve(double** QR, size_t n, double* t, double* b)
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
        d += QR[i][k] * b[i];

      /* b = b - tau (v) (v'b) */
      d  *= tau;
      b[k] -= d;
      for (i = k+1; i < n; i++)
      {
        b[i] -= QR[i][k]*d;
      }
    }
  }
  // Solve Rx = b
  for(i = n-1; i >= 1; i--)
  {
    d = b[i] /= QR[i][i];
    for(j = 0; j < i; j++)
      b[j] -= QR[j][i]*d;
  };
  b[0] /= QR[0][0];
}



void qr_solve_t(double** QRt, size_t n, double* t, double* b)
{

  size_t k;
  double tau,d;
  double *pb, *pQR, *p_end;
  // Q'b
  for(k = 0; k < n; k++)
  {

    if((tau = t[k]) != 0.)
    {
      /* d = v'b */
      pb  = b + k;
      d  = *pb++;
      p_end = b + n;
      pQR = QRt[k] + k+1;
      while(pb < p_end)
        d += *pQR++ * *pb++;

      /* b = b - tau (v) (v'b) */
      pb  = b + k;
      d  *= tau;
      *pb++ -= d;
      pQR = QRt[k] + k+1;
      while(pb < p_end)
      {
        *pb++ -= *pQR++ * d;
      }
    }
  }
  // Solve Rx = b
  for(k = n-1; k >= 1; k--)
  {
    d  = b[k] /= QRt[k][k];
    pb  = b;
    p_end = b + k;
    pQR  = QRt[k];
    while(pb < p_end)
      *pb++ -= *pQR++ * d;
  };
  b[0] /= QRt[0][0];
}


void qr_qt_mul_b(double** QR, size_t m, size_t n, double* t, double* b)
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
        d += QR[i][k] * b[i];

      /* b = b - tau (v) (v'b) */
      b[k] -= tau*d;
      for (i = k+1; i < m; i++)
      {
        b[i] -= tau*QR[i][k]*d;
      }
    }
  }
}



void qr_qt_mul_b_t(double** QRt, size_t m, size_t n, double* t, double* b)
{
  size_t i,k;
  double tau,d;
  // b = Q'*b
  for(k = 0; k < m; k++)
  {
    if((tau = t[k]) != 0.)
    {
      /* d = v'b */
      d = b[k];
      for (i = k+1; i < n; i++)
        d += QRt[k][i] * b[i];

      /* b = b - tau (v) (v'b) */
      b[k] -= tau*d;
      for (i = k+1; i < n; i++)
      {
        b[i] -= tau*QRt[k][i]*d;
      }
    }
  }
}



void qr_q_mul_b(double** QR, size_t m, size_t n, double* t, double* b)
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
        d += QR[i][k] * b[i];

      /* b = b - tau (v) (v'b) */
      b[k] -= tau*d;
      for (i = k+1; i < m; i++)
      {
        b[i] -= tau*QR[i][k]*d;
      }
    }
    if(k == 0) break;
  }
}



void qr_q_mul_b_t(double** QRt, size_t m, size_t n, double* t, double* b)
{
  size_t i,k;
  double tau,d;
  // b = Q*b
  for(k = m-1; ; k--)
  {
    if((tau = t[k]) != 0.)
    {
      /* d = v'b */
      d = b[k];
      for (i = k+1; i < n; i++)
        d += QRt[k][i] * b[i];

      /* b = b - tau (v) (v'b) */
      b[k] -= tau*d;
      for (i = k+1; i < n; i++)
      {
        b[i] -= tau*QRt[k][i]*d;
      }
    }
    if(k == 0) break;
  }
}

void qr_least_squares(double** QR, size_t m, size_t n, double* t, double* b, double* r)
{
  size_t i,j;
  double d;
  double* bb = nl_dvector_create_copy(b,m);

  // b = Q'b
  qr_qt_mul_b(QR, m,n,t,b);

  // Solve Rx = b
  for(i = n-1; i >= 1; i--)
  {
    d = b[i] /= QR[i][i];
    for(j = 0; j < i; j++)
      b[j] -= QR[j][i]*d;
  };
  b[0] /= QR[0][0];

  // r = R*b
  for(j = 0; j < n; j++)
  {
    d = 0;
    for(i = j; i < n; i++)
      d += QR[j][i]*b[i];
    r[j] = d;
  }
  for(j = n; j < m; j++)
    r[j] = 0.;
  // r = Q*R*b
  qr_q_mul_b(QR,m,n,t,r);
  nl_dvector_sub(r, bb, m);
  nl_dvector_free(bb);
}



void qr_least_squares_t(double** QRt, size_t m, size_t n, double* t, double* b, double* r)
{
  size_t i,j;
  double d;
  double* bb = nl_dvector_create_copy(b,n);

  // b = Q'b
  qr_qt_mul_b_t(QRt, m,n,t,b);

  // Solve Rx = b
  for(i = m-1; i >= 1; i--)
  {
    d = b[i] /= QRt[i][i];
    for(j = 0; j < i; j++)
      b[j] -= QRt[i][j]*d;
  };
  b[0] /= QRt[0][0];

  // r = R*b
  for(j = 0; j < m; j++)
  {
    d = 0;
    for(i = j; i < m; i++)
      d += QRt[i][j]*b[i];
    r[j] = d;
  }
  for(j = m; j < n; j++)
    r[j] = 0.;
  // r = Q*R*b
  qr_q_mul_b_t(QRt,n,m,t,r);
  nl_dvector_sub(r, bb, n);
  nl_dvector_free(bb);
}
