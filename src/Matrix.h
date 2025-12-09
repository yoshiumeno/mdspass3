#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <math.h>
#include <iostream>
#include "Vector.h"

class Matrix
{
public:
  Matrix();
  Matrix(int m, int n) : m_nRow(m), m_nCol(n)
  {
    m_elem.resize( m * n, 0.0 );
  };
  Matrix(const Matrix& mat) 
  : m_nRow(mat.m_nRow), m_nCol(mat.m_nCol), m_elem(mat.m_elem) {}
  ~Matrix() {}

  Vector operator*(const Vector& v) const
  {
#ifdef _DEBUG_
    assert(m_nCol == v.size());
#endif
    Vector vRet(m_nRow);
    for(int i = 0; i < m_nRow; ++i) {
      double d = 0.0;
      for(int j = 0; j < m_nCol; ++j) {
        d += (*this)[i][j] * v[j];
      }
      vRet[i] = d;
    }
    return vRet;
  }
  const double* operator[](int i) const
  {
#ifdef _DEBUG_
    assert(-1 < i && i < m_nRow);
#endif
    return &(m_elem[i * m_nCol]);
  }
  double* operator[](int i)
  {
#ifdef _DEBUG_
    assert(-1 < i && i < m_nRow);
#endif
    return &(m_elem[i * m_nCol]);
  }
  Matrix operator-() const
  {
    Matrix mat(*this);
    for(int i = 0; i < m_elem.size(); ++i) {
      mat.m_elem[i] = - mat.m_elem[i];
    }
    return mat;
  }
  Matrix operator*(const Matrix& mat) const
  {
#ifdef _DEBUG_
    assert(m_nCol == mat.m_nRow);
#endif
    int nRow = m_nRow;
    int nCol = mat.m_nCol;
    Matrix matRet(nRow, nCol);
    for(int i = 0; i < nRow; ++i) {
      for(int j = 0; j < nCol; ++j) {
        double d = 0.0;
        for(int k = 0; k < m_nCol; ++k) {
          d += (*this)[i][k] * mat[k][j];
        }
        matRet[i][j] = d;
      }
    }
    return matRet;
  }

  void SetIdentity();
  int SolveLinearEquation(const Vector& p, Vector& q) const;
  Matrix Trans() const;

public:
  vector<double> m_elem;
  int m_nRow;
  int m_nCol;
};

#endif // _MATRIX_H_

