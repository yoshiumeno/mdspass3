
#include <iostream>
#include <iomanip>
#include <math.h>
#include "Chromosome.h"
#include "Utils.h"

void Chromosome::RandomInit(const ReaxPotential& pot)
{
  int nV = pot.GetVariableNum();

  const vector<int>& viVar = pot.GetVariableIndex();

  resize(nV);
  for(int i = 0; i < nV; ++i) {
    const PotentialParam& param = pot[viVar[i]];
    double val = param.m_dValue;
    double min = param.m_dMin;
    double max = param.m_dMax;
    double d = Utils::NormalDistribution() / 3.0;
    d = fmax(d, -1.0);
    d = fmin(d,  1.0);
    double diff = ( 0.0 < d ) ? (max - val) : (val - min);
    val += d * diff;
    (*this)[i] = val;
  }
  
  SetRandomF();
  m_dCR = Utils::UniformRandom();
}

void Chromosome::AssignTo(ReaxPotential& pot) const
{
  const vector<int> viVars = pot.GetVariableIndex();
  int nVar = pot.GetVariableNum();
  for(int i = 0; i < nVar; ++i) {
    pot[viVars[i]].m_dValue = (*this)[i];
  }
}

void Chromosome::ExtractFrom(const ReaxPotential& pot)
{
  const vector<int> viVars = pot.GetVariableIndex();
  int nVar = pot.GetVariableNum();
  resize(nVar, 0.0);
  for(int i = 0; i < nVar; ++i) {
    (*this)[i] = pot[viVars[i]].m_dValue;
  }
}

