
#ifndef _UTILS_H_
#define _UTILS_H_

#include <stdlib.h>
#include <math.h>
#include "String.h"
#include "PFMPI.h"
#include <time.h>

namespace {
  const double PI               = 3.1415926535897932384642;
  const double DEGREE_TO_RADIAN = PI / 180.0;
  const double COULOMB_C        = 332.06371;
  const double KCALMOL_TO_EV    = 4.33641039e-02;
  const double EV_TO_KCALMOL    = 1.0 / KCALMOL_TO_EV;

  inline double SQR(double d) { return d * d; }
  inline double CUBE(double d) { return d * d * d; }
  inline double POWI(double d, int n)
  {
    double dRet = d;
    for(int i = 0; i < n - 1; ++i) {
      dRet *= d;
    }
    return dRet;
  }
  inline int MIN(int i1, int i2) { return (i1 < i2) ? i1 : i2; }
  inline int MAX(int i1, int i2) { return (i1 > i2) ? i1 : i2; }
}

class Utils
{
public:

  static void Error(const String& msg, int iLine = -1,
                    const String& line = "");

  static void SetRandomSeed(int n)
  {
    if( n < 0 ) {
      n = time(NULL);
      PFMPI::Bcast(n, 1);
    }
    srand(n);
  }

  // return uniform random number between [0, 1)
  static double UniformRandom()
  {
    return (double)rand()/((double)RAND_MAX+1.0);
  };

  static double NormalDistribution()
  {
    return sqrt( -2.0 * log( UniformRandom() ) ) 
      * sin( 2.0 * PI * UniformRandom() );
  }
  static int RandomInt(int n)
  {
    return rand() % n;
  }

};

#endif //_UTILS_H_
