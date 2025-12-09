#ifndef _REAX_FF_H_
#define _REAX_FF_H_

#include "ReaxPotential.h"
#include "Configuration.h"
#include "ReaxParams.h"
#include "Input.h"
#include "Chromosome.h"

using namespace std;

class ReaxAtom : public RAtom
{
public:
  ReaxAtom() 
    : RAtom(), m_iOrg(0), m_iX(0), m_iY(0), m_iZ(0),
      m_totalBO(0.0), m_DeltaLP(0.0), m_nlp(0.0),
      m_dEdDelta(0.0), m_dDeltaLPdDelta(0.0) {}
  ReaxAtom(const RAtom& atom, int iOrg) 
    : RAtom(atom), m_iOrg(iOrg), m_iX(0), m_iY(0), m_iZ(0),
      m_totalBO(0.0), m_DeltaLP(0.0), m_nlp(0.0), 
      m_dEdDelta(0.0), m_dDeltaLPdDelta(0.0) {}

public:
  int    m_iOrg;
  int    m_iX;
  int    m_iY;
  int    m_iZ;
  double m_totalBO;
  double m_DeltaLP;
  double m_nlp;

  // for force calculation
  double  m_dEdDelta;
  double  m_dDeltaLPdDelta;
};

class ReaxBond
{
public:
  ReaxBond() 
    : m_iAtom1(-1), m_iAtom2(-1), m_iSymmetry(-1),
      m_BO(0.0), m_BO_Sg(0.0), m_BO_Pi(0.0), m_BO_PP(0.0),
      m_dEdBO(0.0), m_dEdBO_Pi(0.0), m_dEdBO_PP(0.0), 
      m_dBOPdr(0.0), m_dBOP_Sgdr(0.0), m_dBOP_Pidr(0.0),
      m_dBOP_PPdr(0.0), m_dBOdBOP(0.0), m_dBOdDeltaP1(0.0),
      m_dBO_PidBOP_Pi(0.0), m_dBO_PidBOP(0.0),
      m_dBO_PidDeltaP1(0.0), m_dBO_PPdBOP_PP(0.0), 
      m_dBO_PPdBOP(0.0), m_dBO_PPdDeltaP1(0.0)
  {}
  
  ReaxBond(int i1, int i2)
    : m_iAtom1(i1), m_iAtom2(i2), m_iSymmetry(-1),
      m_BO(0.0), m_BO_Sg(0.0), m_BO_Pi(0.0), m_BO_PP(0.0),
      m_dEdBO(0.0), m_dEdBO_Pi(0.0), m_dEdBO_PP(0.0), 
      m_dBOPdr(0.0), m_dBOP_Sgdr(0.0), m_dBOP_Pidr(0.0),
      m_dBOP_PPdr(0.0), m_dBOdBOP(0.0), m_dBOdDeltaP1(0.0),
      m_dBO_PidBOP_Pi(0.0), m_dBO_PidBOP(0.0),
      m_dBO_PidDeltaP1(0.0), m_dBO_PPdBOP_PP(0.0),
      m_dBO_PPdBOP(0.0), m_dBO_PPdDeltaP1(0.0)
  {}

public:
  int m_iAtom1;    // atoms in the original cell
  int m_iAtom2;    // atoms in neighbor cells
  int m_iSymmetry; // index of symmetric bond

  double m_BO;    // Bond Order (might be Bond Order Prime)
  double m_BO_Sg; // Bond Order Sigma
  double m_BO_Pi; // Bond Order Pi
  double m_BO_PP; // Bond Order Pi-Pi

  double m_dEdBO;    // dE/dBO
  double m_dEdBO_Pi; // dE/dBO_Pi
  double m_dEdBO_PP; // dE/dBO_PP

  double m_dBOPdr;    // dBOP/dr
  double m_dBOP_Sgdr; // dBOP_Sg/dr
  double m_dBOP_Pidr; // dBOP_Pi/dr
  double m_dBOP_PPdr; // dBOP_PP/dr

  double m_dBOdBOP;    // dBO/dBOP
  double m_dBOdDeltaP1; // dBO/dDeltaP1

  double m_dBO_PidBOP_Pi;  // dBO_Pi/dBOP_Pi
  double m_dBO_PidBOP;     // dBO_Pi/dBOP
  double m_dBO_PidDeltaP1; // dBO_Pi/dDeltaP1

  double m_dBO_PPdBOP_PP;  // dBO_PP/dBOPPP
  double m_dBO_PPdBOP;     // dBO_PP/dBOP
  double m_dBO_PPdDeltaP1; // dBO_PP/dDeltaP1
};

class ReaxFF
{
public:
  static bool Initialize(
    ReaxPotential& pot, const ConfigurationSet& confSet,
    const Input& input);
  //    const ConfigurationSet& confSetRef, const Input& input);
  static double CalcCost(Chromosome& chrome,
                         Vector& energy_force);
  static void ShowForce(const int i);
  static void GetForce(const int i, double& fx, double& fy, double& fz);
  static void GetStress(const int i, double& fx, double& fy, double& fz);
  static void GetDStress(const int i, double& fx, double& fy, double& fz);
  static void CalcAll();
  static void CalcAll(double& etot);
  static int GetNumberOfAtoms();
  static void SetRAtomPos(int i, double x, double y, double z);
  static void SetLattice(double c11, double c12, double c13,
			 double c21, double c22, double c23,
			 double c31, double c32, double c33);
  
private:
  ReaxFF();
  ReaxFF(const Configuration& conf) : m_conf(conf) {}
  
  double CalcEnergy();
  double CalcEnergyBond();
  double CalcEnergyLonePair();
  double CalcEnergyOverUndercoordination();
  double CalcEnergyValenceAngle();
  double CalcEnergyTorsion();
  double CalcEnergyHydrogenBond();
  double CalcEnergyC2Correction();
  double CalcEnergyTripleBond();
  double CalcEnergyVanDerWaalsCoulomb();
  double CalcEnergyPolarization();

  double GetAngle(
    int i1, int i2, int i3, double& cost, double& sint,
    Vector3& dtdv1, Vector3& dtdv2, Vector3& dtdv3) const;
  double GetTorsion(
    int i1, int i2, int i3, int i4, double& coso, double& sino,
    Vector3& domegadv1, Vector3& domegadv2, 
    Vector3& domegadv3, Vector3& domegadv4 ) const;
  bool IsPeriod12(int it) const;
  void MakeAtomsNeigh();
  void MakeBonds();
  void BondOrderPrime(ReaxBond& bond) const;
  void CorrectBondOrder(ReaxBond& bond) const;
  void TotalBondOrder();
  void TotalForce();
  void CalcChargeQEq();
  void SetBondSymmetryIndex();

private:
  Configuration        m_conf;

  // atoms in neighbor cells
  // first nAtom atoms are same as m_conf.m_atoms
  vector<ReaxAtom>     m_atomsNeigh;

  // indices of bonds for each atom
  vector<vector<int> > m_viBondsOnAtom;

  // bond list
  // atom 1 is in the cell
  // atom 2 is in the cell or neighbor cells
  // this list includes both (i, j) and (j, i) 
  // when i and j are in the cell
  vector<ReaxBond>     m_bonds;

  vector<Vector3>      m_force;
  vector<Vector3>      m_stress;
  vector<Vector3>      d_stress;
  
  static ReaxPotential      s_pot;
  static ReaxParams     s_params;
  static double         s_dEnergyWeight;
  static vector<ReaxFF> s_reax;
  //static vector<ReaxFF> s_reaxRef;
  static bool           s_bInitialized;
};


#endif // _REAX_FF_H_



