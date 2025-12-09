#ifndef _POWELL_MINIMIZATION_H_
#define _POWELL_MINIMIZATION_H_

#include "ReaxPotential.h"
#include "Configuration.h"
#include "Input.h"
#include "Matrix.h"
#include "Chromosome.h"

using namespace std;

class LinePoint {
public:
  LinePoint() : x(0.0), f(0.0) {}
  LinePoint(double _x, double _f, const Vector& _v) 
    : x(_x), f(_f), v(_v) {}

public:
  double x; // point in line
  double f; // cost
  Vector v; // energy and force vector
};


class PowellMinimization
{
public:
  static bool Run(ReaxPotential& pot, const Input& input);

private:
  static bool InitGamma(
    const Chromosome& chrom, const Vector& vef, 
    const vector<double>& vMin, const vector<double>& vMax,
    Matrix& gamma, Matrix& d);
  static void SetupLinearEquationSystem(
    const Matrix& gamma, const Vector& vef, Matrix& les,
    Vector& p);
  static double LineMinimization(
    Chromosome& chrom, const Vector& delta, Vector& vef1,
    LinePoint& min1, LinePoint& min2);
  static bool UpdateGamma(
    const LinePoint& min1, const LinePoint& min2, int j,
    Matrix& gamma, Vector& delta);
  static bool Bracket(
    const Chromosome& chrom0, const Vector& delta,
    LinePoint& l, LinePoint& c, LinePoint& r);
  static double Brent(
    const LinePoint& l, const LinePoint& c, const LinePoint& r,
    const Chromosome& chrom0, const Vector& delta, 
    LinePoint& min1, LinePoint& min2);
};

#endif // _POWELL_MINIMIZATION_H_



