
//#include <mkl_lapack.h>
//#include <clapack.h>
#include "Matrix.h"

extern "C"
void dsysvx_(const char*, const char*, int*, int*, const double*, int*, double*,
	     int*, int*, const double*, int*, double*, int*, double*,
	     double*, double*, double*, int*, int*, int*);

void Matrix::SetIdentity()
{
#ifdef _DEBUG_
  assert(m_nRow == m_nCol);
#endif
  for(int i = 0; i < m_nRow; ++i ) {
    for(int j = 0; j < m_nCol; ++j ) {
      (*this)[i][j] = ( i == j ) ? 1.0 : 0.0;
    }
  }
}

Matrix Matrix::Trans() const
{
  Matrix mat(m_nCol, m_nRow);
  for(int i = 0; i < m_nRow; ++i) {
    for(int j = 0; j < m_nCol; ++j) {
      mat[j][i] = (*this)[i][j];
    }
  }
  return mat;
}

int Matrix::SolveLinearEquation(const Vector& p, Vector& q) const
{
#ifdef _DEBUG_
  assert(m_nRow == m_nCol);
  assert(m_nRow == p.size());
  assert(m_nRow == q.size());
#endif
  int n = p.size();

  int info     = 0;
  int nrhs     = 1;
  Matrix inv(n, n);
  vector<int> ipiv(n);
  double rcond = 0.0;
  double ferr  = 0.0;
  double berr  = 0.0;
  int lwork    = 64 * n;
  vector<double> work(lwork);
  vector<int>    iwork(n);


  dsysvx_("N", "U", &n, &nrhs, &(m_elem[0]), &n, &(inv[0][0]), 
         &n, &(ipiv[0]), &(p[0]), &n, &(q[0]), &n, &rcond, 
         &ferr, &berr, &(work[0]), &lwork, &(iwork[0]), &info);

  return info;
}
