#ifndef _VECTOR3_H_
#define _VECTOR3_H_

#include <math.h>
#include <iostream>

using namespace std;

class Vector3
{
public:
  Vector3() : x(0.0), y(0.0), z(0.0) {}
  Vector3(double xx, double yy, double zz)
    : x(xx), y(yy), z(zz) {}
  Vector3(const Vector3& v)
    : x(v.x), y(v.y), z(v.z) {}
  ~Vector3() {}

  Vector3& operator=(const Vector3& v)
  {
    x = v.x;
    y = v.y;
    z = v.z;
    return *this;
  }
  Vector3 operator+(const Vector3& v) const
  {
    double xx = x + v.x;
    double yy = y + v.y;
    double zz = z + v.z;
    return Vector3(xx, yy, zz);
  }
  Vector3 operator-(const Vector3& v) const
  {
    double xx = x - v.x;
    double yy = y - v.y;
    double zz = z - v.z;
    return Vector3(xx, yy, zz);
  }
  Vector3 operator-() const
  {
    return Vector3(-x, -y, -z);
  }
  double operator*(const Vector3& v) const
  {
    return x * v.x + y * v.y + z * v.z;
  }
  Vector3 operator*(double d) const
  {
    Vector3 v;
    v.x = x * d;
    v.y = y * d;
    v.z = z * d;
    return v;
  }
  Vector3 operator/(double d) const
  {
#ifdef __DEBUG__
    assert(d != 0.0);
#endif
    d = 1.0 / d;
    return (*this) * d;
  }
  friend Vector3 operator*(double d, const Vector3& v)
  {
    return v * d;
  }
  Vector3 operator%(const Vector3& v) const
  {
    double xx = y * v.z - z * v.y;
    double yy = z * v.x - x * v.z;
    double zz = x * v.y - y * v.x;
    return Vector3(xx, yy, zz);
  }
  Vector3& operator+=(const Vector3& v)
  {
    x += v.x;
    y += v.y;
    z += v.z;
    return *this;
  }
  Vector3& operator-=(const Vector3& v)
  {
    x -= v.x;
    y -= v.y;
    z -= v.z;
    return *this;
  }
  Vector3& operator*=(double d)
  {
    x *= d;
    y *= d;
    z *= d;
    return *this;
  }
  Vector3& operator/=(double d)
  {
#ifdef __DEBUG__
    assert(d != 0.0);
#endif
    (*this) = (*this) / d;
    return *this;
  }
  friend ostream& operator<<(ostream& os, const Vector3& v)
  {
    os << v.x << " " << v.y << " " << v.z;
    return os;
  }
  Vector3& Normalize()
  {
    double d = Norm();
#ifdef __DEBUG__
    assert(d != 0.0);
#endif
    (*this) /= d;
    return *this;
  }
  double Norm() const { return sqrt(Norm2()); }
  double Norm2() const { return x * x + y * y + z * z; }

public:
  double x;
  double y;
  double z;
};

namespace {
  const static Vector3 VECTOR3_ZERO = Vector3(0.0, 0.0, 0.0);
}

#endif // _VECTOR3_H_

