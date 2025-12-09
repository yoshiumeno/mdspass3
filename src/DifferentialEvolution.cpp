
#include <iostream>
#include <iomanip>
#include "DifferentialEvolution.h"
#include "Population.h"
#include "ReaxFF.h"
#include "Utils.h"

using namespace std;

static const double TAU_1 = 0.1;   // probability for changing F
static const double TAU_2 = 0.1;   // probability for changing CR
static const double JR    = 0.6;   // jumping rate for opposite algorithm

bool DifferentialEvolution::Run(
  ReaxPotential& pot, const Input& input)
{
  if( PFMPI::s_master ) {
    cout << endl;
    cout << "=== Start Differntial Evolution ===" << endl;
    cout << endl;
  }

  Population pop;
  pop.RandomInit(pot);

  Chromosome chromTrial;
  chromTrial.ExtractFrom(pot);

  int nP = pop.GetChromNum();
  int nV = pop.GetChromSize();

  int iMin = 0;
  double dMin = 0.0;
  for(int i = 1; i < nP; ++i) {
    if( pop[i].m_dCost < dMin ) {
      iMin = i;
      dMin = pop[i].m_dCost;
    }
  }
  Chromosome chromBest(pop[iMin]);

  if( PFMPI::s_master ) {
    cout << endl;
    cout << "Iteration   Min Cost        Max Cost        " 
         << "Ave Cost         Max-Min"
         << endl;
  }

  double dCostMin = 0.0;
  double dCostMax = 0.0;
  double dCostAve = 0.0;
  double dCostDif = 0.0;
  CostStatus(pop, 0, dCostMin, dCostMax, dCostAve, dCostDif);

  const vector<int> viVars = pot.GetVariableIndex();
  Population popNew(pop);
  int    iJumpStep = 0;
  double dJumpRate = JR;

  double dThresh = input.m_dThreshold;

  int iIter = 0;
  for(; iIter < input.m_nMaxIteration; ++iIter ) {
    if( dCostDif < dThresh || dCostMin < dThresh ) {
      break;
    }

    for(int i = 0; i < nP; ++i) {

      int i1, i2;
      do
        i1 = Utils::RandomInt(nP);
      while (i1 == i);

      do
        i2 = Utils::RandomInt(nP);
      while (i2 == i || i2 == i1);

      if( Utils::UniformRandom() < TAU_1 ) {
        chromTrial.SetRandomF();
      }
      else {
        chromTrial.m_dF = pop[i].m_dF;
      }

      if( Utils::UniformRandom() < TAU_2 ) {
        chromTrial.m_dCR = Utils::UniformRandom();
      }
      else {
        chromTrial.m_dCR = pop[i].m_dCR;
      }

      int k = Utils::RandomInt(nV);
      double dF  = chromTrial.m_dF;
      double dCR = chromTrial.m_dCR;

      for( int j = 0; j < nV; ++j )
      {
        if( Utils::UniformRandom() < dCR || j == k ) {
          double d = chromBest[j] + dF * (pop[i1][j] - pop[i2][j]);
          d = fmax(d, pot[viVars[j]].m_dMin);
          d = fmin(d, pot[viVars[j]].m_dMax);
          chromTrial[j] = d;
        }
        else {
          chromTrial[j] = pop[i][j];
        }
      }

      Vector energy_force;
      ReaxFF::CalcCost(chromTrial, energy_force);

      if( chromTrial.m_dCost < dCostMin )
      {
        chromBest = chromTrial;
        dCostMin = chromTrial.m_dCost;
      }

      popNew[i] = ( chromTrial.m_dCost < pop[i].m_dCost ) ?
        chromTrial : pop[i];
    }

    if( Utils::UniformRandom() < dJumpRate ) {
      popNew.CheckOpposite(pot, false);

      ++iJumpStep;
      if( 10 < iJumpStep ) {
        dJumpRate *= 0.9;
        iJumpStep  = 0;
      }
    }

    pop = popNew;

    CostStatus(pop, iIter + 1, dCostMin, dCostMax, dCostAve, dCostDif);
  }

  if( PFMPI::s_master ) {
    cout << endl;
  }

  chromBest.AssignTo(pot);

  if( input.m_nMaxIteration <= iIter ) {
    if( PFMPI::s_master ) {
      cout << "Maximum iteration reached." << endl;
    }
    return false;
  }
  else {
    if( PFMPI::s_master ) {
      cout << "Differential evolution has been converged "
           << "successfully"
           << endl;
    }
    return true;
  }
}

void DifferentialEvolution::CostStatus(
  const Population& pop, int iIter, double& dCostMin, 
  double& dCostMax, double& dCostAve, double& dCostDif)
{
  int nP = pop.GetChromNum();

  dCostMin = pop[0].m_dCost;
  dCostMax = pop[0].m_dCost;
  dCostAve = pop[0].m_dCost;
  for(int i = 1; i < nP; ++i) {
    dCostMin = fmin(dCostMin, pop[i].m_dCost);
    dCostMax = fmax(dCostMax, pop[i].m_dCost);
    dCostAve += pop[i].m_dCost;
  }
  dCostDif = dCostMax - dCostMin;
  dCostAve /= nP;

  if( PFMPI::s_master ) {
    cout << setw(5)  << iIter << " ";
    cout << fixed << setprecision(5);
    cout << setw(15) << right << dCostMin << " ";
    cout << setw(15) << right << dCostMax << " ";
    cout << setw(15) << right << dCostAve << " ";
    cout << setw(15) << right << dCostDif << endl;
  }
}

