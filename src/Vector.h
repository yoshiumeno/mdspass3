#ifndef _VECTOR_H_
#define _VECTOR_H_

#include <math.h>
#include <vector>
#include <assert.h>

using namespace std;

class Vector : public vector<double>
{
public:
  Vector() {}
  Vector(int n) : vector<double>(n, 0.0) {}
  Vector(int n, double d) : vector<double>(n, d) {}
  Vector(const Vector& v) : vector<double>(v) {}
  Vector(const vector<double>& vd) : vector<double>(vd) {}
  ~Vector() {}

  Vector operator+(const Vector& v) const
  {
#ifdef _DEBUG_
    assert(size() == v.size());
#endif
    Vector vRet(*this);
    for(int i = 0; i < size(); ++i) {
      vRet[i] += v[i];
    }
    return vRet;
  }
  Vector operator-(const Vector& v) const
  {
#ifdef _DEBUG_
    assert(size() == v.size());
#endif
    Vector vRet(*this);
    for(int i = 0; i < size(); ++i) {
      vRet[i] -= v[i];
    }
    return vRet;
  }
  double operator*(const Vector& v) const
  {
#ifdef _DEBUG_
    assert(size() == v.size());
#endif
    double d = 0.0;
    for(int i = 0; i < size(); ++i) {
      d += (*this)[i] * v[i];
    }
    return d;
  }
  Vector operator*(double d) const
  {
    Vector v(*this);
    for(int i = 0; i < size(); ++i) {
      v[i] *= d;
    }
    return v;
  }
  friend Vector operator*(double d, const Vector& v)
  {
    return v * d;
  }
  Vector operator/(double d) const
  {
    d = 1.0 / d;
    Vector vRet(*this);
    for(int i = 0; i < size(); ++i) {
      vRet[i] *= d;
    }
    return vRet;
  }
  Vector& operator+=(const Vector& v)
  {
#ifdef _DEBUG_
    assert(size() == v.size());
#endif
    *this = *this + v;
    return *this;
  }
  Vector& operator-=(const Vector& v)
  {
#ifdef _DEBUG_
    assert(size() == v.size());
#endif
    *this = *this - v;
    return *this;
  }
  Vector& operator*=(double d)
  {
    *this = *this * d;
    return *this;
  }
  Vector& operator/=(double d)
  {
    *this = *this / d;
    return *this;
  }
  Vector& operator-()
  {
    for(int i = 0; i < size(); ++i) {
      (*this)[i] = - (*this)[i];
    }
    return *this;
  }
  double Norm() const {
    double d = 0.0;
    for(int i = 0; i < size(); ++i) {
      d += ((*this)[i])*((*this)[i]);
    }
    return sqrt(d);
  }
  void Normalize() {
    double d = Norm();
    (*this) /= d;
  }
};

#endif // _VECTOR_H_

