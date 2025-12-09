
#include <math.h>
#include <map>
#include <iostream>
#include <algorithm>
#include "Population.h"
#include "Utils.h"
#include "ReaxFF.h"


void Population::RandomInit(const ReaxPotential& pot)
{
  if( PFMPI::s_master ) {
    cout << "Initializing population ... ";
  }

  int nV = pot.GetVariableNum();
  int nP = 15 * (nV + 2);

  m_chroms.resize(nP);
  for(int i = 0; i < nP; ++i) {
    m_chroms[i].RandomInit(pot);
  }

  Vector energy_force;
  for( int i = 0; i < nP; ++i ) {
    ReaxFF::CalcCost(m_chroms[i], energy_force);
  }

  CheckOpposite(pot, true);

  if( PFMPI::s_master ) {
    cout << "done" << endl;
  }
}

void Population::CheckOpposite(const ReaxPotential& pot, bool bInit)
{
  return;

  int nP = GetChromNum();
  int nV = GetChromSize();

  if( nP == 0 ) return;

  vector<double> vMin(nV);
  vector<double> vMax(nV);
  if( bInit ) {
    const vector<int>& viVar = pot.GetVariableIndex();
    for(int j = 0; j < nV; ++j) {
      vMin[j] = pot[viVar[j]].m_dMin;
      vMax[j] = pot[viVar[j]].m_dMax;
    }
  }
  else {
    for(int j = 0; j < nV; ++j) {
      vMin[j] = m_chroms[0][j];
      vMax[j] = m_chroms[0][j];
    }
    for(int i = 1; i < nP; ++i) {
      for(int j = 0; j < nV; ++j) {
        vMin[j] = fmin(vMin[j], m_chroms[i][j]);
        vMax[j] = fmax(vMax[j], m_chroms[i][j]);
      }
    }
  }

  Population opp(*this);
  for(int i = 0; i < nP; ++i) {
    for(int j = 0; j < nV; ++j) {
      opp[i][j] = vMin[j] + vMax[j] - m_chroms[i][j];
    }
  }

  Vector energy_force;
  for( int i = 0; i < nP; ++i ) {
    ReaxFF::CalcCost(opp[i], energy_force);
  }

  m_chroms.insert(m_chroms.end(), opp.m_chroms.begin(),
                  opp.m_chroms.end());
  sort(m_chroms.begin(), m_chroms.end());
  m_chroms.resize(nP);
}

