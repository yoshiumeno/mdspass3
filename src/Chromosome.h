
#ifndef _CHROMOSOME_H_
#define _CHROMOSOME_H_

#include "ReaxPotential.h"
#include "Utils.h"
#include "Vector.h"

using namespace std;

class Chromosome : public Vector
{
public:
  Chromosome() : m_dCost(-1.0), m_dF(0.0), m_dCR(0.0) {}
  Chromosome(int n) : Vector(n), m_dCost(-1.0), m_dF(0.0), m_dCR(0.0) {}
  ~Chromosome() {}

  Chromosome& operator=(const Vector& v) 
  {
    Vector::operator=(v);
    return *this;
  }

  void RandomInit(const ReaxPotential& pot);
  void AssignTo(ReaxPotential& pot) const;
  void ExtractFrom(const ReaxPotential& pot);
  void SetRandomF()
  {
    m_dF = F_LOWER + F_UPPER * Utils::UniformRandom();
  }

public:
  double m_dF;
  double m_dCR;
  double m_dCost;

  static const double F_LOWER; //= 0.1; //2022.01.15
  static const double F_UPPER; //= 0.9; //to fix compile error on Win
};

#endif // _CROMOSOME_H_

