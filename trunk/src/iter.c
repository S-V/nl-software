#include "nl.h"

int iter_bicg(
  size_t *IA, 
  size_t *JA, 
  double *AN, 
  double *b, 
  size_t n, 
  double tol, 
  int max_iter, 
  double *x,
  size_t *IM, 
  size_t *JM, 
  double *MN,
  size_t *IK, 
  size_t *JK, 
  double *KN,
  double *work)
{
  double *r, *r_tilde, *p, *p_tilde, *z, *z_tilde; 
  double r_tilde_z, r_tilde_z_new, norm_b, rel_tol, alpha, beta;   
  size_t j;
  int it;

  r = work;
  r_tilde = r + n;
  p = r_tilde + n;
  p_tilde = p + n;
  z = p_tilde + n;
  z_tilde = z + n;

  sp_mult_col(IA, JA, AN, x, n, r);
  cblas_daxpy(n, -1, b, 1, r, 1);
  cblas_dscal(n, -1, r, 1);

  cblas_dcopy(n, r, 1, r_tilde, 1);

  cblas_dcopy(n, r, 1, p, 1); /* нужно: p = M\r; */
  cblas_dcopy(n, r_tilde, 1, p_tilde, 1); /* ружно: p_tilde = r_tilde/M; */

  cblas_dcopy(n, p, 1, z, 1);
  cblas_dcopy(n, p_tilde, 1, z_tilde, 1);
  r_tilde_z = cblas_ddot(n, r_tilde, 1, z, 1);

  norm_b = cblas_dnrm2(n, b, 1);
  if (norm_b == 0)
  {
    for (j = 0; j < n; j++)
      x[j] = 0;
    return 0;
  }

  rel_tol = norm_b*tol;

  it = 0;

  while (cblas_dnrm2(n, r, 1) > rel_tol && ++it <= max_iter)
  {
    sp_mult_col(IA, JA, AN, p, n, z_tilde); /* here z_tilde is used as temp value A*p */

    alpha = r_tilde_z / cblas_ddot(n, p_tilde, 1, z_tilde, 1);
    cblas_daxpy(n, -alpha, z_tilde, 1, r, 1);

    sp_mult_row(IA, JA, AN, p_tilde, n, n, z_tilde);
    cblas_daxpy(n, -alpha, z_tilde, 1, r_tilde, 1);

    cblas_daxpy(n, alpha, p, 1, x, 1);

    cblas_dcopy(n, r, 1, z, 1); /* нужно: z = M\r; */
    cblas_dcopy(n, r_tilde, 1, z_tilde, 1); /* нужно: z_tilde = r_tilde/M; */
    r_tilde_z_new = cblas_ddot(n, r_tilde, 1, z, 1);

    beta = r_tilde_z_new/r_tilde_z;
    cblas_dscal(n, beta, p, 1);
    cblas_daxpy(n, 1, z, 1, p, 1); 

    cblas_dscal(n, beta, p_tilde, 1);
    cblas_daxpy(n, 1, z_tilde, 1, p_tilde, 1); 

    r_tilde_z = r_tilde_z_new;
  }

  return it;
}
