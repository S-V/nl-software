#include "util.h"
#include "sparse.h"

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
  nl_dvector_axpy(-1, r, b, n, r);

  nl_dvector_copy(r_tilde, r, n);

  nl_dvector_copy(p, r, n); /* p = M\r; */
  nl_dvector_copy(p_tilde, r_tilde, n); /* p_tilde = r_tilde/M; */

  nl_dvector_copy(z, p, n);
  nl_dvector_copy(z_tilde, p_tilde, n);
  r_tilde_z = nl_dvector_dot(r_tilde, z, n);

  norm_b = nl_dvector_norm2(b, n);
  if (norm_b == 0)
  {
    for (j = 0; j < n; j++)
      x[j] = 0;
    return 0;
  }

  rel_tol = norm_b*tol;

  it = 0;

  while (nl_dvector_norm2(r, n) > rel_tol && ++it <= max_iter)
  {
    sp_mult_col(IA, JA, AN, p, n, z_tilde); /* here z_tilde is used as temp value A*p */

    alpha = r_tilde_z / nl_dvector_dot(p_tilde, z_tilde, n);
    nl_dvector_axpy(-alpha, z_tilde, r, n, r);

    sp_mult_row(IA, JA, AN, p_tilde, n, n, z_tilde);
    nl_dvector_axpy(-alpha, z_tilde, r_tilde, n, r_tilde);

    nl_dvector_axpy(alpha, p, x, n, x);

    nl_dvector_copy(z, r, n); /* z = M\r; */
    nl_dvector_copy(z_tilde, r_tilde, n); /* z_tilde = r_tilde/M; */
    r_tilde_z_new = nl_dvector_dot(r_tilde, z, n);

    beta = r_tilde_z_new/r_tilde_z;
    nl_dvector_axpy(beta, p, z, n, p);

    nl_dvector_axpy(beta, p_tilde, z_tilde, n, p_tilde);

    r_tilde_z = r_tilde_z_new;
  }

  return it;
}
