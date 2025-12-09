#ifndef _DIFFERENTIAL_EVOLUTION_H_
#define _DIFFERENTIAL_EVOLUTION_H_

#include "ReaxPotential.h"
#include "Configuration.h"
#include "Input.h"
#include "Population.h"

using namespace std;

class DifferentialEvolution
{
public:
  static bool Run(ReaxPotential& pot, const Input& input);

private:
  static void CostStatus(
    const Population& pop, int iIter, double& dCostMin, 
    double& dCostMax, double& dCostAve, double& dCostDif);
};

#endif // _DIFFERENTIAL_EVOLUTION_H_



