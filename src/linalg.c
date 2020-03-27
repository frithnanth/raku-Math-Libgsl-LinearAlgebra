#include <gsl/gsl_complex.h>
#include <gsl/gsl_linalg.h>

void mgsl_linalg_complex_LU_det(gsl_matrix_complex *LU, int signum, gsl_complex *res)
{
  gsl_complex ret = gsl_linalg_complex_LU_det(LU, signum);
  *res = ret;
}

void mgsl_linalg_complex_LU_sgndet(gsl_matrix_complex *LU, int signum, gsl_complex *res)
{
  gsl_complex ret = gsl_linalg_complex_LU_sgndet(LU, signum);
  *res = ret;
}

void mgsl_linalg_complex_householder_transform(gsl_vector_complex *w, gsl_complex *res)
{
  gsl_complex ret = gsl_linalg_complex_householder_transform(w);
  *res = ret;
}
