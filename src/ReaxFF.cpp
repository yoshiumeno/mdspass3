#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include "ReaxFF.h"
#include "Utils.h"
#include "Matrix.h"
#include "myheader.h"

const double CUT_BOND        = 5.0;
const double CUT_HBOND       = 7.5;
const double CUT_BO          = 0.001;
const double CUT_BO_ANGLE    = 0.001;
const double HBOND_THRESHOLD = 0.01;

ReaxPotential      ReaxFF::s_pot;
ReaxParams     ReaxFF::s_params;
double         ReaxFF::s_dEnergyWeight = 1.0;
vector<ReaxFF> ReaxFF::s_reax;
//vector<ReaxFF> ReaxFF::s_reaxRef;
bool           ReaxFF::s_bInitialized = false;

void ReaxFF::SetRAtomPos(int i, double x, double y, double z)
{
  int ii = i - 1;
  Vector3 f;
  f.x = x; f.y = y; f.z = z;
  s_reax[0].m_conf.m_atoms[ii].m_vPos = f;
}
void ReaxFF::SetLattice(double c11, double c12, double c13,
			double c21, double c22, double c23,
			double c31, double c32, double c33)
{
  Vector3 f1, f2, f3;
  f1.x = c11; f1.y = c12, f1.z = c13;
  f2.x = c21; f2.y = c22, f2.z = c23;
  f3.x = c31; f3.y = c32, f3.z = c33;
  s_reax[0].m_conf.m_lattice.m_vX = f1;
  s_reax[0].m_conf.m_lattice.m_vY = f2;
  s_reax[0].m_conf.m_lattice.m_vZ = f3;
}

int ReaxFF::GetNumberOfAtoms()
{
  //Vector3 f = s_reax[0].m_force[1];
  //printf("SUCCESS %e\n",f.x);
  int na = 0;
  if (&s_reax[0].m_conf != NULL) {
    na = s_reax[0].m_conf.GetRAtomNum();
  }
  return na;
}

void ReaxFF::GetStress(const int i, double& sxx, double& syy, double& szz)
{
  int ii = i - 1;
  Vector3 f = s_reax[0].m_stress[ii] * KCALMOL_TO_EV;
  sxx = f.x; syy = f.y; szz = f.z;
}
void ReaxFF::GetDStress(const int i, double& sxy, double& syz, double& szx)
{
  int ii = i - 1;
  Vector3 d = s_reax[0].d_stress[ii] * KCALMOL_TO_EV;
  sxy = d.x; syz = d.y; szx = d.z;
}

void ReaxFF::GetForce(const int i, double& fx, double& fy, double& fz)
{
  //int na = s_reax[0].m_conf.GetRAtomNum();
  //if (i <= na) {
  int ii = i - 1;
  Vector3 f = s_reax[0].m_force[ii] * KCALMOL_TO_EV;
  fx = f.x; fy = f.y; fz = f.z;
  //} else {
  //  printf("Force of %d-th atom is requested, but I have only %d atoms\n",i,na);
  //}
}

void ReaxFF::ShowForce(const int i)
{
  int na = s_reax[0].m_conf.GetRAtomNum();
  if (i == 0) {
    printf("Number of atoms: %d\n", na);
    for (int ii=0; ii<na; ii++) {
      printf("Force of %d-th atom:  ",ii+1);
      //printf("%d\n",s_reax.size());
      Vector3 f = s_reax[0].m_force[ii];
      printf("%e %e %e\n",f.x,f.y,f.z);
    }
  } else if (i <= na) {
    int ii = i - 1;
      printf("Force of %d-th atom:  ",ii+1);
      Vector3 f = s_reax[0].m_force[ii];
      printf("%e %e %e\n",f.x,f.y,f.z);
  } else {
    printf("Force of %d-th atom is requested, but I have only %d atoms\n",i,na);
  }
}

void ReaxFF::CalcAll()
{
  //cout << "Calculating Energy and Forces ... " << flush;
  s_reax[0].CalcEnergy();
  //cout << "done" << endl;
}
void ReaxFF::CalcAll(double& etot)
{
  //cout << "Calculating Energy and Forces (and getting Etot) ... " << flush;
  etot = s_reax[0].CalcEnergy() * eV;
  //cout << "done" << endl;
}

bool ReaxFF::Initialize(
  ReaxPotential& pot, const ConfigurationSet& confSet,
  const Input& input)
//  const ConfigurationSet& confSetRef, const Input& input)
{
  //if( confSet.size() != confSetRef.size() ) {
  //  Utils::Error("The numbers of configurations and references are diffrent");
  //  return false;
  //}
  if( confSet.size() == 0 ) {
    Utils::Error("No configuration found");
    return false;
  }

  s_pot           = pot;
  s_dEnergyWeight = input.m_dEnergyWeight;

  int nConf = confSet.size();
  cout << "Number of Confs = " << nConf << endl;

  printf("! ReaxFF initialization..\n");
  // For the time being, only one s_reax[] is necessary...
  if (s_reax.size() > 0) {
    printf("! s_reax was created before, so it is cleared.\n");
    s_reax.clear();
  }
  
  for(int i = 0; i < nConf; ++i) {
    s_reax   .push_back( ReaxFF(confSet   [i]) );
    //s_reaxRef.push_back( ReaxFF(confSetRef[i]) );
  }
  for (int ii=0; ii<s_reax.size(); ii++) {
    printf("! s_reax[%d] has %d atoms.\n",ii,s_reax[ii].m_conf.GetRAtomNum());
  }
  
  s_bInitialized = true;

  cout << "Setting up parameters ... " << flush;
  int nty = s_reax[0].m_conf.m_atomTypes.size();
  s_params.Setup(s_pot, nty);
  //s_reax[0].CalcEnergy();
  cout << "done" << endl;

  /*
  if( PFMPI::s_master ) {
    cout << "Checking potential variables ... " << flush;
  }

  // Check dependence of variables
  Chromosome chrom0;
  chrom0.ExtractFrom(pot);

  Vector vef;
  printf("1st calccost\n");
  CalcCost(chrom0, vef);
  */
  /*
  const vector<int>& viVars = pot.GetVariableIndex();
  vector<int> viVarsNew;
  for( int i = 0; i < viVars.size(); ++i ){
    int j = viVars[i];
    double dVal = pot[j].m_dValue;
    double dMax = pot[j].m_dMax;
    double dMin = pot[j].m_dMin;
    double d = 0.001 * (dMax - dMin);
    if( dMax < dVal + d ) {
      d = -d;
    }
    Chromosome chrom(chrom0);
    chrom[i] += d;
    printf("%d-th calccost\n",i);
    CalcCost(chrom, vef);
    if( chrom.m_dCost != chrom0.m_dCost ) {
      viVarsNew.push_back(j);
    }
  }
  pot.SetVariableIndex(viVarsNew);
  */

  s_pot = pot;

  if( PFMPI::s_master ) {
    cout << "done" << endl;
  }

  return true;
}

double ReaxFF::CalcCost(Chromosome& chrom, 
                        Vector& energy_force)
{
  if( !s_bInitialized ) {
    Utils::Error("Initialize ReaxFF before cost calculation.");
    return 0.0;
  }
  
  int nType = s_reax[0].m_conf.m_atomTypes.size();
  
  chrom.AssignTo(s_pot);
  s_params.Setup(s_pot, nType);

  s_reax[0].CalcEnergy();
  energy_force.clear();

  int nConf = s_reax.size();
  int nMyConf = nConf / PFMPI::s_size;
  int nAdd = nConf - nMyConf * PFMPI::s_size;
  int iConfS;
  if( PFMPI::s_rank < nAdd ) {
    ++nMyConf;
    iConfS = nMyConf * PFMPI::s_rank;
  }
  else {
    iConfS = nMyConf * PFMPI::s_rank + nAdd;
  }
  int iConfE = iConfS + nMyConf;
  
  double cost = 0.0;
  for( int i = iConfS; i < iConfE; ++i ){
    const Configuration& conf = s_reax[i].m_conf;
    int nAtom = conf.GetRAtomNum();

    double dEnergy    = 0.0;
    //double dEnergyRef = 0.0;
    dEnergy    = s_reax   [i].CalcEnergy();
    //dEnergyRef = s_reaxRef[i].CalcEnergy();

    //double eCost = (dEnergy - dEnergyRef) / nAtom 
    double eCost = dEnergy / nAtom 
      - conf.m_dEnergy0;

    cost += conf.m_dWeight * s_dEnergyWeight
      * eCost * eCost;
    energy_force.push_back(eCost);
      
    if( conf.m_bUseForce ){
      for(int j = 0; j < nAtom; ++j) {
        Vector3 f = s_reax[i].m_force[j] 
          - conf.m_force0[i];
	//- s_reaxRef[i].m_force[j] - conf.m_force0[i];
        cost += conf.m_dWeight * f.Norm2();
        energy_force.push_back(f.x);
        energy_force.push_back(f.y);
        energy_force.push_back(f.z);
      }
    }
  }

  PFMPI::AllreduceSum(cost);
  PFMPI::AllgathervVector(energy_force);

  chrom.m_dCost = cost;

  return cost;
}

double ReaxFF::CalcEnergy()
{
  //printf("### CalcEnergy is called! ###\n");
  MakeAtomsNeigh();
  if( !m_conf.m_bHasCharge ) {
    CalcChargeQEq();
  }
  MakeBonds();

  double energy = 0.0;
  if( m_conf.m_bUseForce ) {
    m_force.clear();
    m_force.resize(m_atomsNeigh.size(), VECTOR3_ZERO);
    m_stress.clear();
    m_stress.resize(m_atomsNeigh.size(), VECTOR3_ZERO);
    d_stress.clear();
    d_stress.resize(m_atomsNeigh.size(), VECTOR3_ZERO);
  }
  energy += CalcEnergyBond();
  energy += CalcEnergyLonePair();
  energy += CalcEnergyOverUndercoordination();
  energy += CalcEnergyValenceAngle();
  energy += CalcEnergyTorsion();
  energy += CalcEnergyHydrogenBond();
  energy += CalcEnergyC2Correction();
  energy += CalcEnergyTripleBond();
  energy += CalcEnergyVanDerWaalsCoulomb();
  energy += CalcEnergyPolarization();

  if( m_conf.m_bUseForce ) {
    TotalForce();
  }

  return energy * KCALMOL_TO_EV;
}

void ReaxFF::TotalForce()
{
  int nAtom = m_conf.GetRAtomNum();

  // calculate dDeltaP/dv
  vector<Vector3> dDeltaPdv(nAtom, VECTOR3_ZERO);
  for(int i = 0; i < m_bonds.size(); ++i) {
    const ReaxBond& bond = m_bonds[i];
    int i1 = bond.m_iAtom1;
    int i2 = bond.m_iAtom2;
    const ReaxAtom& atom1 = m_atomsNeigh[i1];
    const ReaxAtom& atom2 = m_atomsNeigh[i2];
    const Vector3& v1 = atom1.m_vPos;
    const Vector3& v2 = atom2.m_vPos;
    Vector3 v = v1 - v2;
    Vector3 drdv = v / v.Norm();
    dDeltaPdv[i1] += bond.m_dBOPdr * drdv;
  }

  // gather atom.m_dEdDelta
  for(int i = nAtom; i < m_atomsNeigh.size(); ++i) {
    const ReaxAtom& atom = m_atomsNeigh[i];
    ReaxAtom& atomOrg = m_atomsNeigh[atom.m_iOrg];
    atomOrg.m_dEdDelta += atom.m_dEdDelta;
  }
  for(int i = 0; i < m_atomsNeigh.size(); ++i) {
    ReaxAtom& atom = m_atomsNeigh[i];
    const ReaxAtom& atomOrg = m_atomsNeigh[atom.m_iOrg];
    atom.m_dEdDelta = atomOrg.m_dEdDelta;
  }

  for(int iBond = 0; iBond < m_bonds.size(); ++iBond) {
    const ReaxBond& bond = m_bonds[iBond];
    int i1 = bond.m_iAtom1;
    int i2 = bond.m_iAtom2;
    const ReaxAtom& atom1 = m_atomsNeigh[i1];
    const ReaxAtom& atom2 = m_atomsNeigh[i2];
    const Vector3& v1 = atom1.m_vPos;
    const Vector3& v2 = atom2.m_vPos;
    Vector3 v = v1 - v2;
    Vector3 drdv1 = v / v.Norm();

    const ReaxBond& bondSym = m_bonds[bond.m_iSymmetry];
    double dEdBO    = bond.m_dEdBO    + bondSym.m_dEdBO;
    double dEdBO_Pi = bond.m_dEdBO_Pi + bondSym.m_dEdBO_Pi;
    double dEdBO_PP = bond.m_dEdBO_PP + bondSym.m_dEdBO_PP;

    double dEdBOP = 
      dEdBO    * bond.m_dBOdBOP    + 
      dEdBO_Pi * bond.m_dBO_PidBOP + 
      dEdBO_PP * bond.m_dBO_PPdBOP;

    double dEdDeltaP1 = 
      dEdBO    * bond.m_dBOdDeltaP1    +
      dEdBO_Pi * bond.m_dBO_PidDeltaP1 +
      dEdBO_PP * bond.m_dBO_PPdDeltaP1;

    Vector3 dEdv1;
    // v1 -> r -> BO' -> E
    dEdv1 += dEdBOP * bond.m_dBOPdr * drdv1;

    // v1 -> Delta1' -> E
    dEdv1 += dEdDeltaP1 * dDeltaPdv[i1];

    // v1 -> r -> BO'^pi -> BO^pi -> E
    dEdv1 += dEdBO_Pi * bond.m_dBO_PidBOP_Pi * 
      bond.m_dBOP_Pidr * drdv1;

    // v1 -> r -> BO'^pipi -> BO^pipi -> E
    dEdv1 += dEdBO_PP * bond.m_dBO_PPdBOP_PP * 
      bond.m_dBOP_PPdr * drdv1;

    double dEdDelta12 = atom1.m_dEdDelta + atom2.m_dEdDelta;

    // v1 -> r -> BO' -> BO -> Delta1 & Delta2 -> E
    dEdv1 += dEdDelta12 *
      bond.m_dBOdBOP * bond.m_dBOPdr * drdv1;

    // v1 -> DeltaP1 -> BO -> delta1 & Delta2 -> E
    dEdv1 += dEdDelta12 *
      bond.m_dBOdDeltaP1 * dDeltaPdv[i1];

    m_force[i1] -= dEdv1;
    m_stress[i1].x -= -dEdv1.x * v.x / 2.0;
    m_stress[i1].y -= -dEdv1.y * v.y / 2.0;
    m_stress[i1].z -= -dEdv1.z * v.z / 2.0;
    d_stress[i1].x -= -dEdv1.x * v.y / 2.0;
    d_stress[i1].y -= -dEdv1.y * v.z / 2.0;
    d_stress[i1].z -= -dEdv1.z * v.x / 2.0;
    //m_stress[i2].x -= dEdv1.x * v.x / 2.0;
    //m_stress[i2].y -= dEdv1.y * v.y / 2.0;
    //m_stress[i2].z -= dEdv1.z * v.z / 2.0;

    for(int k = 0; k < m_viBondsOnAtom[i1].size(); ++k) {
      int iBond13 = m_viBondsOnAtom[i1][k];
      const ReaxBond& bond13 = m_bonds[iBond13];
      int i3 = bond13.m_iAtom2;
      const ReaxAtom& atom3 = m_atomsNeigh[i3];
      const Vector3& v3 = atom3.m_vPos;
      Vector3 v = v3 - v1;
      Vector3 drdv3 = v / v.Norm();

      Vector3 dEdv3;

      // v3 -> r13 -> BOP13 -> DeltaP1 -> E
      dEdv3 += dEdDeltaP1 * bond13.m_dBOPdr * drdv3;

      // v3 -> r13 -> BOP13 -> DeltaP1 -> BO12 -> 
      // Delta1 & Delta2 -> E
      dEdv3 += dEdDelta12 * bond.m_dBOdDeltaP1 *
        bond13.m_dBOPdr * drdv3;

      m_force[i3] -= dEdv3;
      m_stress[i3].x -= -dEdv3.x * v.x / 2.0;
      m_stress[i3].y -= -dEdv3.y * v.y / 2.0;
      m_stress[i3].z -= -dEdv3.z * v.z / 2.0;
      m_stress[i1].x -= -dEdv3.x * v.x / 2.0;
      m_stress[i1].y -= -dEdv3.y * v.y / 2.0;
      m_stress[i1].z -= -dEdv3.z * v.z / 2.0;
      d_stress[i3].x -= -dEdv3.x * v.y / 2.0;
      d_stress[i3].y -= -dEdv3.y * v.z / 2.0;
      d_stress[i3].z -= -dEdv3.z * v.x / 2.0;
      d_stress[i1].x -= -dEdv3.x * v.y / 2.0;
      d_stress[i1].y -= -dEdv3.y * v.z / 2.0;
      d_stress[i1].z -= -dEdv3.z * v.x / 2.0;
    }
  }

  // gather force
  for(int i = nAtom; i < m_atomsNeigh.size(); ++i) {
    int iOrg = m_atomsNeigh[i].m_iOrg;
    m_force[iOrg] += m_force[i];
    m_stress[iOrg] += m_stress[i];
    d_stress[iOrg] += d_stress[i];
  }
  m_force.resize(nAtom);
  m_stress.resize(nAtom);
  d_stress.resize(nAtom);
}

double ReaxFF::CalcEnergyBond()
{
  double eBond = 0.0;

  for(int i = 0; i < m_bonds.size(); ++i) {
    ReaxBond& bond    = m_bonds[i];
    if( bond.m_iSymmetry < i ) continue;
    int i1 = bond.m_iAtom1;
    int i2 = bond.m_iAtom2;
    const ReaxAtom& atom1 = m_atomsNeigh[i1];
    const ReaxAtom& atom2 = m_atomsNeigh[i2];
    int it1 = atom1.m_iType;
    int it2 = atom2.m_iType;
    double pbe1 = s_params.pbe1[it1][it2];
    double pbe2 = s_params.pbe2[it1][it2];
    double DeSg = s_params.DeSg[it1][it2];
    double DePi = s_params.DePi[it1][it2];
    double DePP = s_params.DePP[it1][it2];

    double powSg = pow(bond.m_BO_Sg, pbe2);
    double expSg = DeSg * exp( pbe1 * ( 1.0 - powSg ) );
      
    eBond -= 
      expSg * bond.m_BO_Sg + 
      DePi  * bond.m_BO_Pi + 
      DePP  * bond.m_BO_PP;

    if( m_conf.m_bUseForce ) {
      double dEdBO = - expSg * ( 1.0 - pbe1 * pbe2 * powSg );
      bond.m_dEdBO    += dEdBO;
      bond.m_dEdBO_Pi -= dEdBO + DePi;
      bond.m_dEdBO_PP -= dEdBO + DePP;
    }
  }
  
  return eBond;
}

// LonePair must be called before OverUndercoodination
// and ValenceAngle
double ReaxFF::CalcEnergyLonePair()
{
  double eLP = 0.0;

  double plp1 = s_params.plp1;
  for(int i = 0; i < m_conf.GetRAtomNum(); ++i) {
    if( m_viBondsOnAtom[i].empty() ) continue;
    ReaxAtom& atom = m_atomsNeigh[i];
    int it = atom.m_iType;
    double plp2 = s_params.plp2[it];
    double valE = s_params.valE[it];
    double val  = s_params.val [it];

    double DeltaE  = atom.m_totalBO - valE;
    double d       = 2.0 + DeltaE  - 2.0 * int(0.5 * DeltaE);
    double exp1    = exp(-plp1 * d * d);
    double nlp     = - int(0.5 * DeltaE) + exp1;
    double nlp_opt = 0.5 * ( valE - val );
    double DeltaLP = nlp_opt - nlp;

    double expLP   = exp(-75.0 * DeltaLP);
    double expLP1  = 1.0 + expLP;
    double expLP1i = 1.0 / expLP1;
    double e       = plp2 * DeltaLP * expLP1i;

    eLP += e;

    atom.m_DeltaLP = DeltaLP; // reuse in Overcoodination
    atom.m_nlp     = nlp;     // reuse in ValenceAngle

    if( m_conf.m_bUseForce ) {
      double dEdDeltaLP = (plp2 + e * 75.0 * expLP) * expLP1i;
      double dDeltaLPdDeltaE = plp1 * exp1 * 2.0 * d;
      atom.m_dEdDelta += dEdDeltaLP * dDeltaLPdDeltaE;
 
      // reuse in Overcoordination and ValenceAngle
      atom.m_dDeltaLPdDelta = dDeltaLPdDeltaE;
    }
  }
  
  for(int i = 0; i < m_atomsNeigh.size(); ++i) {
    ReaxAtom& atom = m_atomsNeigh[i];
    atom.m_DeltaLP = m_atomsNeigh[atom.m_iOrg].m_DeltaLP;
    atom.m_nlp     = m_atomsNeigh[atom.m_iOrg].m_nlp;
  }

  if( m_conf.m_bUseForce ) {
    for(int i = 0; i < m_atomsNeigh.size(); ++i) {
      ReaxAtom& atom    = m_atomsNeigh[i];
      ReaxAtom& atomOrg = m_atomsNeigh[atom.m_iOrg];
      atom.m_dDeltaLPdDelta = atomOrg.m_dDeltaLPdDelta;
    }
  }

  return eLP;
}

double ReaxFF::CalcEnergyOverUndercoordination()
{
  double eOver  = 0.0;
  double eUnder = 0.0;
  
  double povun3 = s_params.povun3;
  double povun4 = s_params.povun4;
  double povun6 = s_params.povun6;
  double povun7 = s_params.povun7;
  double povun8 = s_params.povun8;

  int nAtom = m_conf.GetRAtomNum();
  vector<double> sumDBO    (nAtom, 0.0);
  vector<double> sumDeltaBO(nAtom, 0.0);
  for(int i = 0; i < m_bonds.size(); ++i) {
    int i1         = m_bonds[i].m_iAtom1;
    int i2         = m_bonds[i].m_iAtom2;
    double BO      = m_bonds[i].m_BO;
    double BO_Pi   = m_bonds[i].m_BO_Pi;
    double BO_PP   = m_bonds[i].m_BO_PP;
    const ReaxAtom& atom1 = m_atomsNeigh[i1];
    const ReaxAtom& atom2 = m_atomsNeigh[i2];
    int it1 = atom1.m_iType;
    int it2 = atom2.m_iType;
    double povun1   = s_params.povun1[it1][it2];
    double DeSg     = s_params.DeSg  [it1][it2];
    double val2     = s_params.val[it2];
    double DeltaLP2 = atom2.m_DeltaLP;
    double Delta2   = atom2.m_totalBO - val2;

    // only for period 1 and 2 element
    double DeltaLP2mask = IsPeriod12(it2) ? 1.0 : 0.0;

    sumDBO    [i1] += povun1 * DeSg * BO;
    sumDeltaBO[i1] += (Delta2 - DeltaLP2mask * DeltaLP2) * 
      (BO_Pi + BO_PP);
  }

  for(int i1 = 0; i1 < nAtom; ++i1) {
    ReaxAtom& atom1 = m_atomsNeigh[i1];
    int it1         = atom1.m_iType;
    double povun2   = s_params.povun2[it1];
    double povun5   = s_params.povun5[it1];
    double val1     = s_params.val   [it1];
    double DeltaLP1 = atom1.m_DeltaLP;
    double Delta1   = atom1.m_totalBO - val1;

    // only for period 1 and 2 element
    double DeltaLP1mask = IsPeriod12(it1) ? 1.0 : 0.0;

    // Over Coordination
    double expp4 = povun3 * exp(povun4 * sumDeltaBO[i1]);
    double expp4i = 1.0 / (1.0 + expp4);
    double DeltaLPC1 = Delta1 - DeltaLP1mask * DeltaLP1 * expp4i;

    double DeltaLPC1Vali = 1.0 / (DeltaLPC1 + val1 + 1.0e-8);
    double expp2  = exp(povun2 * DeltaLPC1);
    double expp2i = 1.0 / (1.0 + expp2);
    double eov = sumDBO[i1] * DeltaLPC1Vali * DeltaLPC1 * expp2i;
    eOver += eov;

    // Under Coordnation
    bool bUnder = !m_viBondsOnAtom[i1].empty();

    double expp2un  = 1.0 / expp2;
    double expp2uni = 1.0 / (1.0 + expp2un);
    double expp6    = exp( povun6 * DeltaLPC1);
    double expp8    = exp( povun8 * sumDeltaBO[i1]);
    double expp7i   = 1.0 / ( 1.0 + povun7 * expp8 );
    double eun      = - povun5 * ( 1.0 - expp6 ) * expp2uni 
      * expp7i;
    if( bUnder ) eUnder += eun;

    if( !m_conf.m_bUseForce ) continue;

    double dDeltaLPC1dDelta1 = 1.0 - 
      DeltaLP1mask * expp4i * atom1.m_dDeltaLPdDelta;

    // dEover/dDelta1
    double dEoverdDeltaLPC1 = - eov * DeltaLPC1Vali
      + sumDBO[i1] * DeltaLPC1Vali * expp2i
      - eov * expp2i * expp2 * povun2;
    atom1.m_dEdDelta += dEoverdDeltaLPC1 * dDeltaLPC1dDelta1;

    // dEunder/dDelta1
    double dEunderdDeltaLPC1 = 0.0;
    if( bUnder ) {
      dEunderdDeltaLPC1 = povun5 * expp7i * povun6 * expp6 
        * expp2uni + eun * povun2 * expp2un * expp2uni;
      atom1.m_dEdDelta += dEunderdDeltaLPC1 
        * dDeltaLPC1dDelta1;
    }

    for(int j = 0; j < m_viBondsOnAtom[i1].size(); ++j) {
      int iBond = m_viBondsOnAtom[i1][j];
      ReaxBond& bond12    = m_bonds[iBond];
      int i2 = bond12.m_iAtom2;
      ReaxAtom& atom2 = m_atomsNeigh[i2];
      int it2 = atom2.m_iType;

      // dEover/dDelta2
      double BO_PiPP = bond12.m_BO_Pi + bond12.m_BO_PP;
      double DeltaLP2mask = IsPeriod12(it2) ? 1.0 : 0.0;
      double dDeltaLPC1dDelta2 = DeltaLP1mask * 
        DeltaLP1 * SQR(expp4i) * expp4 * povun4 * 
        (1.0 - DeltaLP2mask * atom2.m_dDeltaLPdDelta) * BO_PiPP;
      atom2.m_dEdDelta += dEoverdDeltaLPC1 
        * dDeltaLPC1dDelta2;

      // dEover/dBO
      double povun1 = s_params.povun1[it1][it2];
      double DeSg   = s_params.DeSg  [it1][it2];
      double dEoverdBO  = povun1 * DeSg * DeltaLPC1Vali 
        * DeltaLPC1 * expp2i;
      bond12.m_dEdBO += dEoverdBO;
        
      // dEover/dBO_Pi, dEover/dBO_PP
      double val2     = s_params.val[it2];
      double DeltaLP2 = atom2.m_DeltaLP;
      double Delta2   = atom2.m_totalBO - val2;
      double dDeltaLPC1dBO_Pi = DeltaLP1mask * 
        DeltaLP1 * SQR(expp4i) * expp4 * povun4 * 
        (Delta2 - DeltaLP2mask * DeltaLP2);

      double dEoverdBO_Pi = dEoverdDeltaLPC1 
        * dDeltaLPC1dBO_Pi;
      bond12.m_dEdBO_Pi += dEoverdBO_Pi;
      bond12.m_dEdBO_PP += dEoverdBO_Pi;

      // dEunder/dDelta2
      double dEunderdDelta2 = eun * (expp7i - 1.0) * povun8 * 
        (1.0 - DeltaLP2mask * atom2.m_dDeltaLPdDelta) * BO_PiPP;
      dEunderdDelta2 += dEunderdDeltaLPC1 * dDeltaLPC1dDelta2;
      atom2.m_dEdDelta += dEunderdDelta2;

      // dEunder/dBO_Pi, dEunder/dBO_PP
      double dEunderdBO_Pi = eun * (expp7i - 1.0) * povun8 * 
        (Delta2 - DeltaLP2mask * DeltaLP2);
      dEunderdBO_Pi += dEunderdDeltaLPC1 * dDeltaLPC1dBO_Pi;
      bond12.m_dEdBO_Pi += dEunderdBO_Pi;
      bond12.m_dEdBO_PP += dEunderdBO_Pi;
    }
  }

  return eOver + eUnder;
}


double ReaxFF::CalcEnergyValenceAngle()
{
  const double CUT_BO_ANGLE_SQR = CUT_BO_ANGLE * CUT_BO_ANGLE;

  double eAngle     = 0.0;
  double ePenalty   = 0.0;
  double eCoalition = 0.0;

  double pval6  = s_params.pval6;
  double pval8  = s_params.pval8;
  double pval9  = s_params.pval9;
  double pval10 = s_params.pval10;
  double ppen2  = s_params.ppen2;
  double ppen3  = s_params.ppen3;
  double ppen4  = s_params.ppen4;
  double pcoa2  = s_params.pcoa2;
  double pcoa3  = s_params.pcoa3;
  double pcoa4  = s_params.pcoa4;


  int nAtom = m_conf.GetRAtomNum();
  for(int i2 = 0; i2 < nAtom; ++i2) {
    ReaxAtom& atom2 = m_atomsNeigh[i2];
    int it2 = atom2.m_iType;

    double pval3     = s_params.pval3 [it2];
    double pval5     = s_params.pval5 [it2];
    double val2      = s_params.val   [it2];
    double valang2   = s_params.valang[it2];
    double valE2     = s_params.valE  [it2];
    double Delta2    = atom2.m_totalBO - val2;
    double Deltaang2 = atom2.m_totalBO - valang2;
    double DeltaE2   = atom2.m_totalBO - valE2;
    double nlp2      = atom2.m_nlp;
    bool   bnlp0     = 0.0 <= DeltaE2 - 2.0 * int(0.5 * DeltaE2);
    if( bnlp0 ) {
      nlp2 = 0.0;
    }
    
    double totalBOP = 0.0;
    double totalBO8 = 0.0;
    for(int i = 0; i < m_viBondsOnAtom[i2].size(); ++i) {
      int iBond = m_viBondsOnAtom[i2][i];
      const ReaxBond& bond = m_bonds[iBond];
      totalBOP += bond.m_BO_Pi + bond.m_BO_PP;
      totalBO8 += POWI(bond.m_BO, 8);
    }
    double expBO8 = exp(-totalBO8);
    double SBO = totalBOP - (1.0 - expBO8) * 
      (Deltaang2 + pval8 * nlp2);

    double SBO2 = 0.0;
    if( SBO < 0.0 ) {
      SBO2 = 0.0;
    }
    else if( SBO < 1.0 ) {
      SBO2 = pow(SBO, pval9);
    }
    else if( SBO < 2.0 ) {
      SBO2 = 2.0 - pow(2.0 - SBO, pval9);
    }
    else {
      SBO2 = 2.0;
    }

    double dSBOdDelta2  = 0.0;
    double dSBO2dSBO    = 0.0;
    if( m_conf.m_bUseForce && 0.0 < SBO2 ) {
      dSBOdDelta2 = expBO8 - 1.0;
      if( !bnlp0 ) {
        double dnlpdDelta2 = - atom2.m_dDeltaLPdDelta;
        dSBOdDelta2 *= (1.0 + pval8 * dnlpdDelta2);
      }
      if( SBO < 1.0 ) {
        dSBO2dSBO =   pval9 * SBO2 / SBO;
      }
      else {
        dSBO2dSBO = - pval9 * (2.0 - SBO2) / (2.0 - SBO);
      }
    }

    // for penalty energy
    double expp3   = exp( - ppen3 * Delta2 );
    double expp4   = exp(   ppen4 * Delta2 );
    double expp34i = 1.0 / (1.0 + expp3 + expp4);
    double f9     = (2.0 + expp3) * expp34i;

    // for coalition energy
    double valval2   = s_params.valval[it2];
    double Deltaval2 = atom2.m_totalBO - valval2;
    double expp2     = exp( pcoa2 * Deltaval2 );
    double expp2i    = 1.0 / ( 1.0 + expp2 );

    int nBonds = m_viBondsOnAtom[i2].size();
    for(int j1 = 0; j1 < nBonds; ++j1) {
      int iBond21 = m_viBondsOnAtom[i2][j1];
      ReaxBond& bond21 = m_bonds[iBond21];
      double bo21 = bond21.m_BO;
      double BOA21 = bo21 - CUT_BO_ANGLE;
      if( BOA21 <= 0.0 ) continue;
      int i1 = bond21.m_iAtom2;
      ReaxAtom& atom1 = m_atomsNeigh[i1];
      int it1 = atom1.m_iType;

      for(int j3 = j1 + 1; j3 < nBonds; ++j3) {
        int iBond23 = m_viBondsOnAtom[i2][j3];
        ReaxBond& bond23 = m_bonds[iBond23];
        double BO23 = bond23.m_BO;
        double BOA23 = BO23 - CUT_BO_ANGLE;
        if( BOA23 <= 0.0 ) continue;
        if( bo21 * BO23 <= CUT_BO_ANGLE_SQR ) continue;  
        int i3 = bond23.m_iAtom2;
        ReaxAtom& atom3 = m_atomsNeigh[i3];
        int it3 = atom3.m_iType;

        // === Valence Angle Energy ===
        double pval1   = s_params.pval1  [it1][it2][it3];
        double pval2   = s_params.pval2  [it1][it2][it3];
        double pval4   = s_params.pval4  [it1][it2][it3];
        double pval7   = s_params.pval7  [it1][it2][it3];
        double Theta00 = s_params.Theta00[it1][it2][it3];
        Theta00 *= DEGREE_TO_RADIAN;

        double pow21p4 = pow(BOA21, pval4);
        double pow23p4 = pow(BOA23, pval4);
        double expp3_21 = exp( -pval3 * pow21p4);
        double expp3_23 = exp( -pval3 * pow23p4);
        double f7_21 = 1.0 - expp3_21;
        double f7_23 = 1.0 - expp3_23;

        double expp6   = exp(  pval6 * Deltaang2);
        double expp7   = exp(- pval7 * Deltaang2);
        double expp67i = 1.0 / (1.0 + expp6 + expp7);
        double f8 = pval5 - (pval5 - 1.0) * (2.0 + expp6) 
          * expp67i;

        double expp10 = exp(-pval10 * (2.0 - SBO2) );
        double Theta0 = PI - Theta00 * (1.0 - expp10 );
        
        double cosTheta = 1.0;
        double sinTheta = 0.0;
        Vector3 dThetadv1;
        Vector3 dThetadv2;
        Vector3 dThetadv3;
        double Theta = GetAngle(i1, i2, i3, cosTheta, sinTheta,
                                dThetadv1, dThetadv2, dThetadv3);
        double dTheta2 = SQR(Theta0 - Theta);

        double p1expp2    = pval1 * exp(-pval2 * dTheta2 );
        double p1expp2_p1 = - p1expp2;
        if( 0 <= pval1 ){
          p1expp2_p1 += pval1;
        }

        double eang = f7_21 * f7_23 * f8 * p1expp2_p1;
        eAngle += eang;

        // === Penalty Energy ===
        double ppen1 = s_params.ppen1[it1][it2][it3];

        double d21 = BOA21 - 2.0;
        double d23 = BOA23 - 2.0;
        double exp21 = exp( - ppen2 * d21 * d21 );
        double exp23 = exp( - ppen2 * d23 * d23 );

        double epen = ppen1 * f9 * exp21 * exp23;

        ePenalty += epen;

        // === Coalition Energy ===
        double pcoa1 = s_params.pcoa1[it1][it2][it3];

        double d21p3 = BOA21 - atom1.m_totalBO;
        double d23p3 = BOA23 - atom3.m_totalBO;
        double exppc3_21 = exp( - pcoa3 * SQR(d21p3) );
        double exppc3_23 = exp( - pcoa3 * SQR(d23p3) );
        double expp3 = exppc3_21 * exppc3_23;

        double d21p4 = BOA21 - 1.5;
        double d23p4 = BOA23 - 1.5;
        double expp4_21 = exp( - pcoa4 * SQR(d21p4) );
        double expp4_23 = exp( - pcoa4 * SQR(d23p4) );
        double expp4 = expp4_21 * expp4_23;

        double ecoa = pcoa1 * expp2i * expp3 * expp4;

        eCoalition += ecoa;

        if( !m_conf.m_bUseForce ) continue;

        // === Valence Angle Force ===

        // dEang/dBO21, dEang/dBO23
        double d = pval3 * pval4;
        double df7dBO21 = d * pow21p4 / BOA21 * expp3_21;
        double df7dBO23 = d * pow23p4 / BOA23 * expp3_23;
        double dEangdBO21 = df7dBO21 * f7_23 * f8 * p1expp2_p1;
        double dEangdBO23 = f7_21 * df7dBO23 * f8 * p1expp2_p1;
        bond21.m_dEdBO += dEangdBO21;
        bond23.m_dEdBO += dEangdBO23;

        // dEang/dDelta2
        double dEangdDelta2 = 0.0;
          
        double df8dDelta2 = (1.0 - pval5) * expp67i 
          * (pval6 * expp6 - (2.0 + expp6) * 
             (pval6 * expp6 - pval7 * expp7) * expp67i);
        dEangdDelta2 += f7_21 * f7_23 * df8dDelta2 
          * p1expp2_p1;
          
        double dEangdTheta0 = - 2.0 * f7_21 * f7_23 * f8 
          * p1expp2 * pval2 * (Theta0 - Theta);
        double dTheta0dSBO2 = - Theta00 * expp10 * pval10;
        double dEangdSBO2 = dEangdTheta0 * dTheta0dSBO2;
        double dEangdSBO  = dEangdSBO2 * dSBO2dSBO;
        dEangdDelta2 += dEangdSBO * dSBOdDelta2;
          
        atom2.m_dEdDelta += dEangdDelta2;

        for(int j4 = 0; j4 < nBonds; ++j4) {
          int iBond24 = m_viBondsOnAtom[i2][j4];
          ReaxBond& bond24 = m_bonds[iBond24];

          // dEang/dBO24
          double dSBOdBO24 =  - 8.0 * POWI(bond24.m_BO, 7)
            * expBO8 * (Deltaang2 + pval8 * nlp2);
          double dEangdBO24 = dEangdSBO * dSBOdBO24;
          bond24.m_dEdBO += dEangdBO24;

          // dEang/dBO_Pi24, dEang/dBO_PP24
          bond24.m_dEdBO_Pi += dEangdSBO;
          bond24.m_dEdBO_PP += dEangdSBO;
        }

        // dEang/dv1, dEang/dv2, dEang/dv3 via Theta
        //double dEangdTheta = - dEangdTheta0;
        double dEangdTheta = dEangdTheta0; // is this ok?
        Vector3 dEangdv1 = dEangdTheta * dThetadv1;
        Vector3 dEangdv2 = dEangdTheta * dThetadv2;
        Vector3 dEangdv3 = dEangdTheta * dThetadv3;
        m_force[i1] -= dEangdv1;
        m_force[i2] -= dEangdv2;
	m_force[i3] -= dEangdv3;
	Vector3 v = atom1.m_vPos - atom2.m_vPos;////////
	m_stress[i1].x -= dEangdv1.x * v.x / 2.0*2;
	m_stress[i1].y -= dEangdv1.y * v.y / 2.0*2;
	m_stress[i1].z -= dEangdv1.z * v.z / 2.0*2;
	d_stress[i1].x -= dEangdv1.x * v.y / 2.0*2;
	d_stress[i1].y -= dEangdv1.y * v.z / 2.0*2;
	d_stress[i1].z -= dEangdv1.z * v.x / 2.0*2;
	v = atom2.m_vPos - atom1.m_vPos;
	v*=0;
	//m_stress[i2].x -= dEangdv2.x * v.x / 2.0;
	//m_stress[i2].y -= dEangdv2.y * v.y / 2.0;
	//m_stress[i2].z -= dEangdv2.z * v.z / 2.0;
	//d_stress[i2].x -= dEangdv2.x * v.y / 2.0;
	//d_stress[i2].y -= dEangdv2.y * v.z / 2.0;
	//d_stress[i2].z -= dEangdv2.z * v.x / 2.0;
	v = atom3.m_vPos - atom2.m_vPos;
	m_stress[i3].x -= dEangdv3.x * v.x / 2.0*2;
	m_stress[i3].y -= dEangdv3.y * v.y / 2.0*2;
	m_stress[i3].z -= dEangdv3.z * v.z / 2.0*2;
	d_stress[i3].x -= dEangdv3.x * v.y / 2.0*2;
	d_stress[i3].y -= dEangdv3.y * v.z / 2.0*2;
	d_stress[i3].z -= dEangdv3.z * v.x / 2.0*2;
	
        // === Penalty Force ===

        // dEpen/dBO21, dEpen/dBO23
        double dEpendBO21 = - 2.0 * ppen2 * d21 * epen;
        double dEpendBO23 = - 2.0 * ppen2 * d23 * epen;
        bond21.m_dEdBO += dEpendBO21;
        bond23.m_dEdBO += dEpendBO23;

        // dEpen/dDelta2
        double df9dDelta2 = - ppen3 * expp3
          + (2.0 + expp3) * (ppen3 * expp3 - ppen4 * expp4) 
          * expp34i;
        df9dDelta2 *= expp34i;
        double dEpendDelta2 = ppen1 * df9dDelta2 * exp21 * exp23;
        atom2.m_dEdDelta += dEpendDelta2;

        // === Coalition Force ===

        // dEcoa/dBO21, dEcoa/dBO23
        double ecoa2 = -2.0 * ecoa;
        double dEcoadBO21 = ecoa2 * (pcoa3 * d21p3 + pcoa4 * d21p4);
        double dEcoadBO23 = ecoa2 * (pcoa3 * d23p3 + pcoa4 * d23p4);
        bond21.m_dEdBO += dEcoadBO21;
        bond23.m_dEdBO += dEcoadBO23;

        // dEcoa/dDelta2
        double dEcoadDelta2 = - pcoa2 * expp2i * expp2 * ecoa;
        atom2.m_dEdDelta += dEcoadDelta2;

        // dEcoa/dDelta1, dEcoa/dDelta3
        double dEcoadDelta1 = - ecoa2 * pcoa3 * d21p3;
        double dEcoadDelta3 = - ecoa2 * pcoa3 * d23p3;
        atom1.m_dEdDelta += dEcoadDelta1;
        atom3.m_dEdDelta += dEcoadDelta3;
      }
    }
  }

  return eAngle + ePenalty + eCoalition;
}

double ReaxFF::CalcEnergyTorsion()
{
  const double CUT_BO_ANGLE_SQR = CUT_BO_ANGLE * CUT_BO_ANGLE;

  double eTorsion   = 0.0;
  double eConjugate = 0.0;

  double ptor2 = s_params.ptor2;
  double ptor3 = s_params.ptor3;
  double ptor4 = s_params.ptor4;
  double pcot2 = s_params.pcot2;

  int nAtom = m_conf.GetRAtomNum();

  // only atom2 is in the cell
  // atom1, atom3 and atom4 might be outside of the cell
  for(int iBond23 = 0; iBond23 < m_bonds.size(); ++iBond23) {

    ReaxBond& bond23 = m_bonds[iBond23];
    if( bond23.m_iSymmetry < iBond23 ) continue;
    int i2 = bond23.m_iAtom1;
    int i3 = bond23.m_iAtom2;
    ReaxAtom& atom2 = m_atomsNeigh[i2];
    ReaxAtom& atom3 = m_atomsNeigh[i3];
    int it2 = atom2.m_iType;
    int it3 = atom3.m_iType;
    double valang2    = s_params.valang[it2];
    double valang3    = s_params.valang[it3];
    double Deltaang2  = atom2.m_totalBO - valang2;
    double Deltaang3  = atom3.m_totalBO - valang3;
    double Deltaang23 = Deltaang2 + Deltaang3;
    double BO23       = bond23.m_BO;
    double BO_Pi23    = bond23.m_BO_Pi;
    double BOA23      = BO23 - CUT_BO_ANGLE;
    if( BOA23 <= 0.0 ) continue;

    int iOrg3 = atom3.m_iOrg;

    double expp3 = exp( - ptor3 * Deltaang23 );
    double expp4 = exp(   ptor4 * Deltaang23 );
    double expp34i = 1.0 / (1.0 + expp3 + expp4);
    double f11  = (2.0 + expp3 ) * expp34i;

    for(int j2 = 0; j2 < m_viBondsOnAtom[i2].size(); ++j2) {
      int iBond12 = m_viBondsOnAtom[i2][j2];
      ReaxBond& bond12 = m_bonds[iBond12];
      int i1 = bond12.m_iAtom2;
      if(i1 == i3) continue;
      double BO12  = bond12.m_BO;
      double BOA12 = BO12 - CUT_BO_ANGLE;
      if( BOA12 <= 0.0 ) continue;
      if( BO12 * BO23 <= CUT_BO_ANGLE_SQR ) continue;

      double cos123 = 1.0;
      double sin123 = 0.0;
      Vector3 dt123dv1;
      Vector3 dt123dv2;
      Vector3 dt123dv3;
      double Theta123 = GetAngle(i1, i2, i3, cos123, sin123,
                                 dt123dv1, dt123dv2, dt123dv3);

      if( sin123 < 1.0e-10 ) continue;

      const ReaxAtom& atom1 = m_atomsNeigh[i1];
      int iOrg1 = atom1.m_iOrg;
      int it1 = atom1.m_iType;

      for(int j3 = 0; j3 < m_viBondsOnAtom[iOrg3].size(); ++j3) {
        int iBond34 = m_viBondsOnAtom[iOrg3][j3];
        ReaxBond& bond34 = m_bonds[iBond34];
        int i4 = bond34.m_iAtom2;
        const ReaxAtom& atom4 = m_atomsNeigh[i4];
        int iOrg4 = atom4.m_iOrg;
        // check atom4 is not same of atom2
        int iX34 = atom3.m_iX + atom4.m_iX;
        int iY34 = atom3.m_iY + atom4.m_iY;
        int iZ34 = atom3.m_iZ + atom4.m_iZ;
        if( iX34 == 0 && iY34 == 0 && iZ34 == 0 && 
            iOrg4 == i2 ) continue;
        // check atom4 is not same of atom1
        if( iX34 == atom1.m_iX &&
            iY34 == atom1.m_iY && 
            iZ34 == atom1.m_iZ && 
            iOrg4 == iOrg1 ) continue;

        double BO34 = bond34.m_BO;
        double BOA34 = BO34 - CUT_BO_ANGLE;
        if( BOA34 <= 0.0 ) continue;
        if( BO23 * BO34 <= CUT_BO_ANGLE_SQR ) continue;
        if( BO12 * BO23 * BO34 <= CUT_BO_ANGLE ) continue;
        int it4 = atom4.m_iType;

        double cos234 = 1.0;
        double sin234 = 0.0;
        Vector3 dt234dv2;
        Vector3 dt234dv3;
        Vector3 dt234dv4;
        double Theta234 = GetAngle(i2, i3, i4, cos234, sin234,
                                   dt234dv2, dt234dv3, dt234dv4);
        if( sin234 < 1.0e-10 ) continue;

        double cosomega = -1.0;
        double sinomega =  0.0;
        Vector3 domegadv1;
        Vector3 domegadv2;
        Vector3 domegadv3;
        Vector3 domegadv4;
        double omega = GetTorsion(
          i1, i2, i3, i4, cosomega, sinomega,
          domegadv1, domegadv2, domegadv3, domegadv4);

        // === Torsion Energy ===
        double ptor1 = s_params.ptor1[it1][it2][it3][it4];
        double v1    = s_params.v1   [it1][it2][it3][it4];
        double v2    = s_params.v2   [it1][it2][it3][it4];
        double v3    = s_params.v3   [it1][it2][it3][it4];

        double expp2_12 = exp(- ptor2 * BOA12);
        double expp2_23 = exp(- ptor2 * BOA23);
        double expp2_34 = exp(- ptor2 * BOA34);
        double f10 = ( 1.0 - expp2_12 ) * ( 1.0 - expp2_23 )
          * ( 1.0 - expp2_34 );

        double vo1 = 0.5 * v1 * (1.0 + cos(      omega));
        double vo2 = 0.5 * v2 * (1.0 - cos(2.0 * omega));
        double vo3 = 0.5 * v3 * (1.0 + cos(3.0 * omega));
        double dPi = BO_Pi23 + f11 - 2.0;
        double expp1 = exp(ptor1 * dPi * dPi);

        double vo = vo1 + vo2 * expp1 + vo3;
        double etor = f10 * sin123 * sin234 * vo;

        eTorsion += etor;

        // === Conjugate Energy ===
        double pcot1 = s_params.pcot1[it1][it2][it3][it4];

        double BOAC12 = BOA12 - 1.5;
        double BOAC23 = BOA23 - 1.5;
        double BOAC34 = BOA34 - 1.5;
        double exppc12 = exp( - pcot2 * SQR(BOAC12) );
        double exppc23 = exp( - pcot2 * SQR(BOAC23) );
        double exppc34 = exp( - pcot2 * SQR(BOAC34) );
        double f12 = exppc12 * exppc23 * exppc34;

        double oss = 1.0 + 
          (SQR(cosomega) - 1.0) * sin123 * sin234 ;
        double econj = f12 * pcot1 * oss;

        eConjugate += econj;

        if( !m_conf.m_bUseForce ) continue;

        // === Torsion Force ===

        // dEtor/dBO_Pi23
        double dEtordBO_Pi23 = f10 * sin123 * sin234
          * vo2 * expp1 * 2.0 * ptor1 * dPi;
        bond23.m_dEdBO_Pi += dEtordBO_Pi23;

        // dEtor/dDelta2, dEtor/dDelta3
        double dEtordf11 = dEtordBO_Pi23;
        double df11dDelta2 = - ptor3 * expp3
          + (2.0 + expp3) * (ptor3 * expp3 - ptor4 * expp4) 
          * expp34i;
        df11dDelta2 *= expp34i;
        double dEtordDelta2 = dEtordf11 * df11dDelta2;
        atom2.m_dEdDelta += dEtordDelta2;
        atom3.m_dEdDelta += dEtordDelta2;

        // dEtor/dBO12, dEtor/dBO23, dEtor/dBO34
        double df10dBO12 = ptor2 * expp2_12 * (1.0 - expp2_23)
          * (1.0 - expp2_34);
        double df10dBO23 = ptor2 * expp2_23 * (1.0 - expp2_34)
          * (1.0 - expp2_12);
        double df10dBO34 = ptor2 * expp2_34 * (1.0 - expp2_12)
          * (1.0 - expp2_23);
        double dEtordf10 = sin123 * sin234 * vo;
        double dEtordBO12 = dEtordf10 * df10dBO12;
        double dEtordBO23 = dEtordf10 * df10dBO23;
        double dEtordBO34 = dEtordf10 * df10dBO34;
        bond12.m_dEdBO += dEtordBO12;
        bond23.m_dEdBO += dEtordBO23;
        bond34.m_dEdBO += dEtordBO34;

        // dEtor/dv1, dtor/dv2 and dEtor/dv3 via Theta123
        double dEtordTheta123 = f10 * cos123 * sin234 * vo;
        m_force[i1] -= dEtordTheta123 * dt123dv1;
        m_force[i2] -= dEtordTheta123 * dt123dv2;
        m_force[i3] -= dEtordTheta123 * dt123dv3;
	Vector3 v = atom1.m_vPos - atom2.m_vPos;/////////
	m_stress[i1].x -= dEtordTheta123 * dt123dv1.x * v.x / 2.0;
	m_stress[i1].y -= dEtordTheta123 * dt123dv1.y * v.y / 2.0;
	m_stress[i1].z -= dEtordTheta123 * dt123dv1.z * v.z / 2.0;
	d_stress[i1].x -= dEtordTheta123 * dt123dv1.x * v.y / 2.0;
	d_stress[i1].y -= dEtordTheta123 * dt123dv1.y * v.z / 2.0;
	d_stress[i1].z -= dEtordTheta123 * dt123dv1.z * v.x / 2.0;
	v = atom2.m_vPos - atom1.m_vPos;/////////
	v*=0;///
	m_stress[i2].x -= dEtordTheta123 * dt123dv2.x * v.x / 2.0;
	m_stress[i2].y -= dEtordTheta123 * dt123dv2.y * v.y / 2.0;
	m_stress[i2].z -= dEtordTheta123 * dt123dv2.z * v.z / 2.0;
	d_stress[i2].x -= dEtordTheta123 * dt123dv2.x * v.y / 2.0;
	d_stress[i2].y -= dEtordTheta123 * dt123dv2.y * v.z / 2.0;
	d_stress[i2].z -= dEtordTheta123 * dt123dv2.z * v.x / 2.0;
	v = atom3.m_vPos - atom2.m_vPos;/////////
	m_stress[i3].x -= dEtordTheta123 * dt123dv3.x * v.x / 2.0;
	m_stress[i3].y -= dEtordTheta123 * dt123dv3.y * v.y / 2.0;
	m_stress[i3].z -= dEtordTheta123 * dt123dv3.z * v.z / 2.0;
	d_stress[i3].x -= dEtordTheta123 * dt123dv3.x * v.y / 2.0;
	d_stress[i3].y -= dEtordTheta123 * dt123dv3.y * v.z / 2.0;
	d_stress[i3].z -= dEtordTheta123 * dt123dv3.z * v.x / 2.0;

        // dEtor/dv2, dEtor/dv3 and dEtor/dv4 via Theta234
        double dEtordTheta234 = f10 * sin123 * cos234 * vo;
        m_force[i2] -= dEtordTheta234 * dt234dv2;
        m_force[i3] -= dEtordTheta234 * dt234dv3;
        m_force[i4] -= dEtordTheta234 * dt234dv4;
	v = atom2.m_vPos - atom3.m_vPos;/////////
	//v*=0;///
	m_stress[i2].x -= dEtordTheta234 * dt234dv2.x * v.x / 2.0;
	m_stress[i2].y -= dEtordTheta234 * dt234dv2.y * v.y / 2.0;
	m_stress[i2].z -= dEtordTheta234 * dt234dv2.z * v.z / 2.0;
	d_stress[i2].x -= dEtordTheta234 * dt234dv2.x * v.y / 2.0;
	d_stress[i2].y -= dEtordTheta234 * dt234dv2.y * v.z / 2.0;
	d_stress[i2].z -= dEtordTheta234 * dt234dv2.z * v.x / 2.0;
	v = atom3.m_vPos - atom2.m_vPos;/////////
	m_stress[i3].x -= dEtordTheta234 * dt234dv3.x * v.x / 2.0;
	m_stress[i3].y -= dEtordTheta234 * dt234dv3.y * v.y / 2.0;
	m_stress[i3].z -= dEtordTheta234 * dt234dv3.z * v.z / 2.0;
	d_stress[i3].x -= dEtordTheta234 * dt234dv3.x * v.y / 2.0;
	d_stress[i3].y -= dEtordTheta234 * dt234dv3.y * v.z / 2.0;
	d_stress[i3].z -= dEtordTheta234 * dt234dv3.z * v.x / 2.0;
	v = atom4.m_vPos - atom2.m_vPos;/////////
	m_stress[i4].x -= dEtordTheta234 * dt234dv4.x * v.x / 2.0;
	m_stress[i4].y -= dEtordTheta234 * dt234dv4.y * v.y / 2.0;
	m_stress[i4].z -= dEtordTheta234 * dt234dv4.z * v.z / 2.0;
	d_stress[i4].x -= dEtordTheta234 * dt234dv4.x * v.y / 2.0;
	d_stress[i4].y -= dEtordTheta234 * dt234dv4.y * v.z / 2.0;
	d_stress[i4].z -= dEtordTheta234 * dt234dv4.z * v.x / 2.0;

        // dEtor/dv1, dEtor/dv2 dEtor/dv3 and dEtor/dv4 via omega
        double sin2omega = 2.0 * sinomega * cosomega;
        double sin3omega = sinomega * 
          (4.0 * SQR(cosomega) - 1.0);
        double dvo1domega = - 0.5 * v1 * sinomega;
        double dvo2domega =   0.5 * v2 * sin2omega * 2.0;
        double dvo3domega = - 0.5 * v3 * sin3omega * 3.0;
        double dvodomega = dvo1domega + dvo2domega * expp1
          + dvo3domega;
        double dEtordomega = f10 * sin123 * sin234 * dvodomega;
        m_force[i1] -= dEtordomega * domegadv1;
        m_force[i2] -= dEtordomega * domegadv2;
        m_force[i3] -= dEtordomega * domegadv3;
        m_force[i4] -= dEtordomega * domegadv4;
	v = atom1.m_vPos - atom2.m_vPos;/////////
	m_stress[i1].x -= dEtordomega * domegadv1.x * v.x / 2.0;
	m_stress[i1].y -= dEtordomega * domegadv1.y * v.y / 2.0;
	m_stress[i1].z -= dEtordomega * domegadv1.z * v.z / 2.0;
	d_stress[i1].x -= dEtordomega * domegadv1.x * v.y / 2.0;
	d_stress[i1].y -= dEtordomega * domegadv1.y * v.z / 2.0;
	d_stress[i1].z -= dEtordomega * domegadv1.z * v.x / 2.0;
	v = atom2.m_vPos - atom1.m_vPos;/////////
	v*=0;///
	m_stress[i2].x -= dEtordomega * domegadv2.x * v.x / 2.0;
	m_stress[i2].y -= dEtordomega * domegadv2.y * v.y / 2.0;
	m_stress[i2].z -= dEtordomega * domegadv2.z * v.z / 2.0;
	d_stress[i2].x -= dEtordomega * domegadv2.x * v.y / 2.0;
	d_stress[i2].y -= dEtordomega * domegadv2.y * v.z / 2.0;
	d_stress[i2].z -= dEtordomega * domegadv2.z * v.x / 2.0;
	v = atom3.m_vPos - atom2.m_vPos;/////////
	m_stress[i3].x -= dEtordomega * domegadv3.x * v.x / 2.0;
	m_stress[i3].y -= dEtordomega * domegadv3.y * v.y / 2.0;
	m_stress[i3].z -= dEtordomega * domegadv3.z * v.z / 2.0;
	d_stress[i3].x -= dEtordomega * domegadv3.x * v.y / 2.0;
	d_stress[i3].y -= dEtordomega * domegadv3.y * v.z / 2.0;
	d_stress[i3].z -= dEtordomega * domegadv3.z * v.x / 2.0;
	v = atom4.m_vPos - atom2.m_vPos;/////////
	m_stress[i4].x -= dEtordomega * domegadv4.x * v.x / 2.0;
	m_stress[i4].y -= dEtordomega * domegadv4.y * v.y / 2.0;
	m_stress[i4].z -= dEtordomega * domegadv4.z * v.z / 2.0;
	d_stress[i4].x -= dEtordomega * domegadv4.x * v.y / 2.0;
	d_stress[i4].y -= dEtordomega * domegadv4.y * v.z / 2.0;
	d_stress[i4].z -= dEtordomega * domegadv4.z * v.x / 2.0;

        // === Conjugation Force ===

        // dEconj/dBO12, dEconj/dBO23, dEconj/dBO34
        double df12dBO12 = - 2.0 * pcot2 * BOAC12 * f12;
        double df12dBO23 = - 2.0 * pcot2 * BOAC23 * f12;
        double df12dBO34 = - 2.0 * pcot2 * BOAC34 * f12;
        double dEconjdBO12 = df12dBO12 * pcot1 * oss;
        double dEconjdBO23 = df12dBO23 * pcot1 * oss;
        double dEconjdBO34 = df12dBO34 * pcot1 * oss;
        bond12.m_dEdBO += dEconjdBO12;
        bond23.m_dEdBO += dEconjdBO23;
        bond34.m_dEdBO += dEconjdBO34;

        // dEconj/dv1, dconj/dv2 and dEconj/dv3 via Theta123
        double dEconjdTheta123 = f12 * pcot1 
          * (SQR(cosomega) - 1.0 ) * cos123 * sin234;
        m_force[i1] -= dEconjdTheta123 * dt123dv1;
        m_force[i2] -= dEconjdTheta123 * dt123dv2;
        m_force[i3] -= dEconjdTheta123 * dt123dv3;
	v = atom1.m_vPos - atom2.m_vPos;/////////
	m_stress[i1].x -= dEconjdTheta123 * dt123dv1.x * v.x / 2.0;
	m_stress[i1].y -= dEconjdTheta123 * dt123dv1.y * v.y / 2.0;
	m_stress[i1].z -= dEconjdTheta123 * dt123dv1.z * v.z / 2.0;
	d_stress[i1].x -= dEconjdTheta123 * dt123dv1.x * v.y / 2.0;
	d_stress[i1].y -= dEconjdTheta123 * dt123dv1.y * v.z / 2.0;
	d_stress[i1].z -= dEconjdTheta123 * dt123dv1.z * v.x / 2.0;
	v = atom2.m_vPos - atom1.m_vPos;/////////
	v*=0;///
	m_stress[i2].x -= dEconjdTheta123 * dt123dv2.x * v.x / 2.0;
	m_stress[i2].y -= dEconjdTheta123 * dt123dv2.y * v.y / 2.0;
	m_stress[i2].z -= dEconjdTheta123 * dt123dv2.z * v.z / 2.0;
	d_stress[i2].x -= dEconjdTheta123 * dt123dv2.x * v.y / 2.0;
	d_stress[i2].y -= dEconjdTheta123 * dt123dv2.y * v.z / 2.0;
	d_stress[i2].z -= dEconjdTheta123 * dt123dv2.z * v.x / 2.0;
	v = atom3.m_vPos - atom2.m_vPos;/////////
	m_stress[i3].x -= dEconjdTheta123 * dt123dv3.x * v.x / 2.0;
	m_stress[i3].y -= dEconjdTheta123 * dt123dv3.y * v.y / 2.0;
	m_stress[i3].z -= dEconjdTheta123 * dt123dv3.z * v.z / 2.0;
	d_stress[i3].x -= dEconjdTheta123 * dt123dv3.x * v.y / 2.0;
	d_stress[i3].y -= dEconjdTheta123 * dt123dv3.y * v.z / 2.0;
	d_stress[i3].z -= dEconjdTheta123 * dt123dv3.z * v.x / 2.0;

        // dEconj/dv1, dconj/dv2 and dEconj/dv3 via Theta234
        double dEconjdTheta234 = f12 * pcot1 
          * (SQR(cosomega) - 1.0 ) * sin123 * cos234;
        m_force[i2] -= dEconjdTheta234 * dt234dv2;
        m_force[i3] -= dEconjdTheta234 * dt234dv3;
        m_force[i4] -= dEconjdTheta234 * dt234dv4;
	v = atom2.m_vPos - atom3.m_vPos;/////////
	v*=0;///
	m_stress[i2].x -= dEconjdTheta234 * dt234dv2.x * v.x / 2.0;
	m_stress[i2].y -= dEconjdTheta234 * dt234dv2.y * v.y / 2.0;
	m_stress[i2].z -= dEconjdTheta234 * dt234dv2.z * v.z / 2.0;
	d_stress[i2].x -= dEconjdTheta234 * dt234dv2.x * v.y / 2.0;
	d_stress[i2].y -= dEconjdTheta234 * dt234dv2.y * v.z / 2.0;
	d_stress[i2].z -= dEconjdTheta234 * dt234dv2.z * v.x / 2.0;
	v = atom3.m_vPos - atom2.m_vPos;/////////
	m_stress[i3].x -= dEconjdTheta234 * dt234dv3.x * v.x / 2.0;
	m_stress[i3].y -= dEconjdTheta234 * dt234dv3.y * v.y / 2.0;
	m_stress[i3].z -= dEconjdTheta234 * dt234dv3.z * v.z / 2.0;
	d_stress[i3].x -= dEconjdTheta234 * dt234dv3.x * v.y / 2.0;
	d_stress[i3].y -= dEconjdTheta234 * dt234dv3.y * v.z / 2.0;
	d_stress[i3].z -= dEconjdTheta234 * dt234dv3.z * v.x / 2.0;
	v = atom4.m_vPos - atom2.m_vPos;/////////
	m_stress[i4].x -= dEconjdTheta234 * dt234dv4.x * v.x / 2.0;
	m_stress[i4].y -= dEconjdTheta234 * dt234dv4.y * v.y / 2.0;
	m_stress[i4].z -= dEconjdTheta234 * dt234dv4.z * v.z / 2.0;
	d_stress[i4].x -= dEconjdTheta234 * dt234dv4.x * v.y / 2.0;
	d_stress[i4].y -= dEconjdTheta234 * dt234dv4.y * v.z / 2.0;
	d_stress[i4].z -= dEconjdTheta234 * dt234dv4.z * v.x / 2.0;

        // dEconj/dv1, dconj/dv2, dconj/dv3 and dEconj/dv4
        // via omega
        double dEconjdomega = - f12 * pcot1 
          * sin2omega * sin123 * sin234;
        m_force[i1] -= dEconjdomega * domegadv1;
        m_force[i2] -= dEconjdomega * domegadv2;
        m_force[i3] -= dEconjdomega * domegadv3;
        m_force[i4] -= dEconjdomega * domegadv4;
	v = atom1.m_vPos - atom2.m_vPos;/////////
	m_stress[i1].x -= dEconjdomega * domegadv1.x * v.x / 2.0;
	m_stress[i1].y -= dEconjdomega * domegadv1.y * v.y / 2.0;
	m_stress[i1].z -= dEconjdomega * domegadv1.z * v.z / 2.0;
	d_stress[i1].x -= dEconjdomega * domegadv1.x * v.y / 2.0;
	d_stress[i1].y -= dEconjdomega * domegadv1.y * v.z / 2.0;
	d_stress[i1].z -= dEconjdomega * domegadv1.z * v.x / 2.0;
	v = atom2.m_vPos - atom1.m_vPos;/////////
	v*=0;///
	m_stress[i2].x -= dEconjdomega * domegadv2.x * v.x / 2.0;
	m_stress[i2].y -= dEconjdomega * domegadv2.y * v.y / 2.0;
	m_stress[i2].z -= dEconjdomega * domegadv2.z * v.z / 2.0;
	d_stress[i2].x -= dEconjdomega * domegadv2.x * v.y / 2.0;
	d_stress[i2].y -= dEconjdomega * domegadv2.y * v.z / 2.0;
	d_stress[i2].z -= dEconjdomega * domegadv2.z * v.x / 2.0;
	v = atom3.m_vPos - atom2.m_vPos;/////////
	m_stress[i3].x -= dEconjdomega * domegadv3.x * v.x / 2.0;
	m_stress[i3].y -= dEconjdomega * domegadv3.y * v.y / 2.0;
	m_stress[i3].z -= dEconjdomega * domegadv3.z * v.z / 2.0;
	d_stress[i3].x -= dEconjdomega * domegadv3.x * v.y / 2.0;
	d_stress[i3].y -= dEconjdomega * domegadv3.y * v.z / 2.0;
	d_stress[i3].z -= dEconjdomega * domegadv3.z * v.x / 2.0;
	v = atom4.m_vPos - atom2.m_vPos;/////////
	m_stress[i4].x -= dEconjdomega * domegadv4.x * v.x / 2.0;
	m_stress[i4].y -= dEconjdomega * domegadv4.y * v.y / 2.0;
	m_stress[i4].z -= dEconjdomega * domegadv4.z * v.z / 2.0;
	d_stress[i4].x -= dEconjdomega * domegadv4.x * v.y / 2.0;
	d_stress[i4].y -= dEconjdomega * domegadv4.y * v.z / 2.0;
	d_stress[i4].z -= dEconjdomega * domegadv4.z * v.x / 2.0;
      }
    } 
  }

  return eTorsion + eConjugate;
}

double ReaxFF::CalcEnergyHydrogenBond()
{
  const double CUTOFF     = fmax(CUT_HBOND, CUT_BOND);
  const double CUTOFF_SQR = CUTOFF * CUTOFF;

  double eHB = 0.0;

  int nAtom = m_conf.GetRAtomNum();
  for(int i2 = 0; i2 < nAtom; ++i2) {
    const ReaxAtom& atom2 = m_atomsNeigh[i2];
    int it2 = atom2.m_iType;
    int ihb2 = (int)round(s_params.phbond[it2]);
    if( ihb2 != 1 ) continue;
    const Vector3& v2 = atom2.m_vPos;

    vector<int> viAtom3;
    for(int i3 = 0; i3 < m_atomsNeigh.size(); ++i3) {
      const ReaxAtom& atom3 = m_atomsNeigh[i3];
      int it3 = atom3.m_iType;
      int ihb3 = (int)round(s_params.phbond[it3]);
      if( ihb3 != 2 ) continue;
      const Vector3& v3 = atom3.m_vPos;
      if( CUTOFF_SQR < (v2 - v3).Norm2() ) continue;
      viAtom3.push_back(i3);
    }

    for(int j1 = 0; j1 < m_viBondsOnAtom[i2].size(); ++j1) {
      int iBond12 = m_viBondsOnAtom[i2][j1];
      ReaxBond& bond12 = m_bonds[iBond12];
      double BO12 = bond12.m_BO;
      if( BO12 < HBOND_THRESHOLD ) continue;
      int i1 = bond12.m_iAtom2;
      const ReaxAtom& atom1 = m_atomsNeigh[i1];
      int it1 = atom1.m_iType;
      int ihb1 = (int)round(s_params.phbond[it1]);
      if( ihb1 != 2 ) continue;
      const Vector3& v1 = atom1.m_vPos;
      double r12 = (v1 - v2).Norm();

      for(int j3 = 0; j3 < viAtom3.size(); ++j3) {
        int i3 = viAtom3[j3];
        if( i3 == i1 ) continue;
        const ReaxAtom& atom3 = m_atomsNeigh[i3];
        if( atom1.m_iOrg == atom3.m_iOrg ) continue;
        double cosTheta = 1.0;
        double sinTheta = 0.0;
        Vector3 dThetadv1;
        Vector3 dThetadv2;
        Vector3 dThetadv3;
        double Theta = GetAngle(i1, i2, i3, cosTheta, sinTheta,
                                dThetadv1, dThetadv2, dThetadv3);
        double sin2 = (1.0 - cosTheta) * 0.5; // sin(Theta/2)^2
        if( sin2 < 1.0e-10 ) continue;
        double sin4  = sin2 * sin2;
        int it3 = atom3.m_iType;

        double r0hb = s_params.r0hb[it1][it2][it3];
        double phb1 = s_params.phb1[it1][it2][it3];
        double phb2 = s_params.phb2[it1][it2][it3];
        double phb3 = s_params.phb3[it1][it2][it3];

        const Vector3& v3 = atom3.m_vPos;
        Vector3 v23  = v3 - v2;
        double r23   = v23.Norm();
        double r23i  = 1.0 / r23;
        double r0hbi = 1.0 / r0hb;
        double expp2 = exp(- phb2 * BO12);
        double expp3 = 
          exp( - phb3 * (r0hb * r23i + r23 * r0hbi - 2.0) );
        double e = phb1 * (1.0 - expp2) * expp3 * sin4;

        eHB += e;

        if( !m_conf.m_bUseForce ) continue;

        // dE/dBO12
        double dEdBO12 = phb1 * phb2 * expp2 * expp3 * sin4;
        bond12.m_dEdBO += dEdBO12;

        // dE/dv1, dE/dv2 and dE/dv3 via Theta
        {
          double dEdTheta = phb1 * (1.0 - expp2) * expp3 
            * sin2 * sinTheta;
          Vector3 dEdv1 = dEdTheta * dThetadv1;
          Vector3 dEdv2 = dEdTheta * dThetadv2;
          Vector3 dEdv3 = dEdTheta * dThetadv3;
          m_force[i1] -= dEdv1;
          m_force[i2] -= dEdv2;
          m_force[i3] -= dEdv3;
	  Vector3 v = atom1.m_vPos - atom2.m_vPos;/////////
	  m_stress[i1].x -= dEdv1.x * v.x / 2.0;
	  m_stress[i1].y -= dEdv1.y * v.y / 2.0;
	  m_stress[i1].z -= dEdv1.z * v.z / 2.0;
	  d_stress[i1].x -= dEdv1.x * v.y / 2.0;
	  d_stress[i1].y -= dEdv1.y * v.z / 2.0;
	  d_stress[i1].z -= dEdv1.z * v.x / 2.0;
	  v = atom2.m_vPos - atom1.m_vPos;/////////
	  v*=0;///
	  m_stress[i2].x -= dEdv2.x * v.x / 2.0;
	  m_stress[i2].y -= dEdv2.y * v.y / 2.0;
	  m_stress[i2].z -= dEdv2.z * v.z / 2.0;
	  d_stress[i2].x -= dEdv2.x * v.y / 2.0;
	  d_stress[i2].y -= dEdv2.y * v.z / 2.0;
	  d_stress[i2].z -= dEdv2.z * v.x / 2.0;
	  v = atom3.m_vPos - atom2.m_vPos;/////////
	  m_stress[i3].x -= dEdv3.x * v.x / 2.0;
	  m_stress[i3].y -= dEdv3.y * v.y / 2.0;
	  m_stress[i3].z -= dEdv3.z * v.z / 2.0;
	  d_stress[i3].x -= dEdv3.x * v.y / 2.0;
	  d_stress[i3].y -= dEdv3.y * v.z / 2.0;
	  d_stress[i3].z -= dEdv3.z * v.x / 2.0;
        }

        // dE/dv2 and dE/dv3 via r23
        {
          double dEdr23 = e * phb3 * (r0hb * SQR(r23i) - r0hbi);
          Vector3 dr23dv3 = v23 * r23i;
          Vector3 dEdv3 = dEdr23 * dr23dv3;
          m_force[i2] -= - dEdv3;
          m_force[i3] -=   dEdv3;
	  Vector3 v = atom2.m_vPos - atom3.m_vPos;/////////
	  m_stress[i2].x -= -dEdv3.x * v.x / 2.0;
	  m_stress[i2].y -= -dEdv3.y * v.y / 2.0;
	  m_stress[i2].z -= -dEdv3.z * v.z / 2.0;
	  d_stress[i2].x -= -dEdv3.x * v.y / 2.0;
	  d_stress[i2].y -= -dEdv3.y * v.z / 2.0;
	  d_stress[i2].z -= -dEdv3.z * v.x / 2.0;
	  v = atom3.m_vPos - atom2.m_vPos;/////////
	  m_stress[i3].x -=  dEdv3.x * v.x / 2.0;
	  m_stress[i3].y -=  dEdv3.y * v.y / 2.0;
	  m_stress[i3].z -=  dEdv3.z * v.z / 2.0;
	  d_stress[i3].x -=  dEdv3.x * v.y / 2.0;
	  d_stress[i3].y -=  dEdv3.y * v.z / 2.0;
	  d_stress[i3].z -=  dEdv3.z * v.x / 2.0;
        }
      }
    }
  }

  return eHB;
}

double ReaxFF::CalcEnergyC2Correction()
{
  double eC2  = 0.0;

  double kc2 = s_params.kc2;

  int nAtom = m_conf.GetRAtomNum();
  for(int i1 = 0; i1 < nAtom; ++i1) {
    ReaxAtom& atom1 = m_atomsNeigh[i1];
    int it1 = atom1.m_iType;
    if( (int)round(s_params.mass[it1]) != 12 ) continue;
    double val1   = s_params.val[it1];
    double Delta1 = atom1.m_totalBO - val1;
    double d4 = SQR(SQR(Delta1)) * 0.04 + Delta1 + 3.0;
    
    for(int j2 = 0; j2 < m_viBondsOnAtom[i1].size(); ++j2) {
      int iBond12 = m_viBondsOnAtom[i1][j2];
      ReaxBond& bond12 = m_bonds[iBond12];
      int i2 = bond12.m_iAtom1;
      if( i1 == i2 ) i2 = bond12.m_iAtom2;
      const ReaxAtom& atom2 = m_atomsNeigh[i2];
      int it2 = atom2.m_iType; 
      if( (int)round(s_params.mass[it2]) != 12 ) continue;
      double bo = bond12.m_BO;
      if( bo <= d4 ) continue;
      double bod4 = bo - d4;
      double e = kc2 * SQR(bod4);
      eC2 += e;

      if( !m_conf.m_bUseForce ) continue;
        
      // dEdBO
      double dEdBO = 2.0 * kc2 * bod4;
      bond12.m_dEdBO += dEdBO;

      // dE/dDelta1
      double dEdDelta1 = - dEdBO * ( 1.0 + 0.16 * CUBE(Delta1) );
      atom1.m_dEdDelta += dEdDelta1;
    }
  }

  return eC2;
}

double ReaxFF::CalcEnergyTripleBond()
{
  bool bAll = (int)round(s_params.vtri) == 2;

  double eTriple  = 0.0;

  double ptri1 = s_params.ptri1;
  double ptri2 = s_params.ptri2;
  double ptri3 = s_params.ptri3;
  double ptri4 = s_params.ptri4;

  for(int i = 0; i < m_bonds.size(); ++i) {
    ReaxBond& bond = m_bonds[i];
    if( bond.m_iSymmetry < i ) continue;
    double bo = bond.m_BO;
    if( bo < 1.0 ) continue;
    int i1 = bond.m_iAtom1;
    int i2 = bond.m_iAtom2;
    ReaxAtom& atom1 = m_atomsNeigh[i1];
    ReaxAtom& atom2 = m_atomsNeigh[i2];
    int it1 = atom1.m_iType;
    int it2 = atom2.m_iType;
    int iMass1 = (int)round(s_params.mass[it1]);
    int iMass2 = (int)round(s_params.mass[it2]);
    if( !bAll ) {
      if( !(iMass1 == 12 && iMass2 == 16) &&
          !(iMass1 == 16 && iMass2 == 12) ) continue;
    }
    double val1 = s_params.val[it1];
    double val2 = s_params.val[it2];
    double totbo1 = atom1.m_totalBO;
    double totbo2 = atom2.m_totalBO;
    double Delta1 = totbo1 - val1;
    double Delta2 = totbo2 - val2;

    double expp2   = exp( - ptri2 * SQR(bo - 2.5) );
    double expp3   = exp(   ptri3 * (Delta1 + Delta2) );
    double expp4_1 = exp( - ptri4 * (totbo1 - bo) );
    double expp4_2 = exp( - ptri4 * (totbo2 - bo) );

    double expp3i = 1.0 / (1.0 + 25.0 * expp3);
    double e = ptri1 * expp2 * (expp4_1 + expp4_2)
      * expp3i;

    eTriple += e;

    if( !m_conf.m_bUseForce ) continue;

    // dE/dBO
    double dEdBO = ( ptri4 - 2.0 * ptri2 * (bo - 2.5) ) * e;
    bond.m_dEdBO += dEdBO;

    // dEdDelta1, dEdDelta2
    double dEdDelta  = - 25.0 * expp3 * ptri3 * expp3i * e;
    double d         = ptri4 * ptri1 * expp2 * expp3i;
    double dEdDelta1 = dEdDelta - d * expp4_1;
    double dEdDelta2 = dEdDelta - d * expp4_2;
    atom1.m_dEdDelta += dEdDelta1;
    atom2.m_dEdDelta += dEdDelta2;
  }

  return eTriple;
}

double ReaxFF::CalcEnergyVanDerWaalsCoulomb()
{
  double eVDW     = 0.0;
  double eCoulomb = 0.0;

  double pvdw1    = s_params.pvdw1;
  double pvdw1i   = 1.0 / pvdw1;
  double rlow     = s_params.nonbond_low;
  double rcut     = s_params.nonbond_cut;
  double rcut_sqr = rcut * rcut;

  int nAtom = m_conf.GetRAtomNum();
  for(int i1 = 0; i1 < nAtom; ++i1) {
    const ReaxAtom& atom1 = m_atomsNeigh[i1];
    const Vector3&  v1    = atom1.m_vPos;
    int             it1   = atom1.m_iType;
    double          q1    = atom1.m_dCharge;
    double dvdw1   = s_params.dvdw1 [it1];
    double rvdw1   = s_params.rvdw1 [it1];
    double alpha1  = s_params.alpha1[it1];
    double gammaw1 = s_params.gammaw[it1];
    double gamma1  = s_params.gamma[it1];
    
    for(int i2 = 0; i2 < m_atomsNeigh.size(); ++i2) {
      if( i1 == i2 ) continue;
      const ReaxAtom& atom2 = m_atomsNeigh[i2];
      int iOrg2 = atom2.m_iOrg;
      if( iOrg2 <  i1 ) continue;
      if( iOrg2 == i1 ) {
        if( atom2.m_iX < 0 ) continue;
        if( atom2.m_iX == 0 ) {
          if( atom2.m_iY < 0 ) continue;
          if( atom2.m_iY == 0 ) {
            if( atom2.m_iZ < 0 ) continue;
          }
        }
      }
      const Vector3&  v2 = atom2.m_vPos;
      Vector3 v12  = v2 - v1;
      double r_sqr = v12.Norm2();
      if( rcut_sqr <= r_sqr ) continue;
      double r = sqrt(r_sqr);
      int it2 = atom2.m_iType;

      // === van der Waals Energy ===
      double dvdw2   = s_params.dvdw1 [it2];
      double rvdw2   = s_params.rvdw1 [it2];
      double alpha2  = s_params.alpha1[it2];
      double gammaw2 = s_params.gammaw[it2];

      double dvdw  = s_params.dvdw [it1][it2];
      double rvdw  = s_params.rvdw [it1][it2];
      double alpha = s_params.alpha[it1][it2];
      if( dvdw  == 0.0 ) dvdw  = sqrt(dvdw1  * dvdw2 );
      if( rvdw  == 0.0 ) rvdw  = sqrt(rvdw1  * rvdw2 );
      if( alpha == 0.0 ) alpha = sqrt(alpha1 * alpha2);
      rvdw *= 2.0;

      double gammaw = sqrt(gammaw1 * gammaw2);

      double rtap = r - rlow;
      double tap = s_params.taps[7];
      for(int i = 0; i < 7; ++i) {
        tap = tap * rtap + s_params.taps[6-i];
      }
      
      double powr = pow(r     ,  pvdw1);
      double powg = pow(gammaw, -pvdw1);
      
      double f13 = pow( powr + powg, pvdw1i );
      double rvdwi = 1.0 / rvdw;
      double d = alpha * (1.0 - f13 * rvdwi);
      double expd  = exp(d);
      double expdh = exp(0.5 * d);
      double evdw0 = dvdw * ( expd - 2.0 * expdh );
      double evdw  = tap * evdw0;

      eVDW += evdw;

      // === Coulomb Energy ===
      double q2 = atom2.m_dCharge;
      double gamma2 = s_params.gamma[it2];
      double g = 1.0 / sqrt(gamma1 * gamma2);

      double rgi = pow(CUBE(r) + CUBE(g), - 1.0 / 3.0);
      double ecoul0 = COULOMB_C * q1 * q2 * rgi;
      double ecoul  = tap * ecoul0;

      eCoulomb += ecoul;

      if( !m_conf.m_bUseForce ) continue;

      Vector3 drdv2 = v12 / r;
      double dtapdr = 7.0 * s_params.taps[7];
      for(int i = 0; i < 6; ++i) {
        dtapdr = dtapdr * rtap + (6.0 - i) * s_params.taps[6-i];
      }

      // === van der Waals force ===
      double df13dr = f13 * powr / ( (powr + powg) * r );
      double devdw0dr  = - dvdw * ( expd - expdh )
        * alpha * rvdwi * df13dr;
      double  dEvdwdr  = dtapdr * evdw0 + tap * devdw0dr;
      Vector3 dEvdwdv2 = dEvdwdr * drdv2;
      m_force[i1] -= - dEvdwdv2;
      m_force[i2] -=   dEvdwdv2;
      m_stress[i1].x -= -dEvdwdv2.x * v12.x / 2.0;
      m_stress[i1].y -= -dEvdwdv2.y * v12.y / 2.0;
      m_stress[i1].z -= -dEvdwdv2.z * v12.z / 2.0;
      m_stress[i2].x -= -dEvdwdv2.x * v12.x / 2.0;
      m_stress[i2].y -= -dEvdwdv2.y * v12.y / 2.0;
      m_stress[i2].z -= -dEvdwdv2.z * v12.z / 2.0;
      d_stress[i1].x -= -dEvdwdv2.x * v12.y / 2.0;
      d_stress[i1].y -= -dEvdwdv2.y * v12.z / 2.0;
      d_stress[i1].z -= -dEvdwdv2.z * v12.x / 2.0;
      d_stress[i2].x -= -dEvdwdv2.x * v12.y / 2.0;
      d_stress[i2].y -= -dEvdwdv2.y * v12.z / 2.0;
      d_stress[i2].z -= -dEvdwdv2.z * v12.x / 2.0;

      // === Coulomb force ===
      double  dEcouldr  = dtapdr * ecoul0 - ecoul * CUBE(rgi) * r * r;
      Vector3 dEcouldv2 = dEcouldr * drdv2;
      m_force[i1] -= - dEcouldv2;
      m_force[i2] -=   dEcouldv2;
      m_stress[i1].x -= -dEcouldv2.x * v12.x / 2.0;
      m_stress[i1].y -= -dEcouldv2.y * v12.y / 2.0;
      m_stress[i1].z -= -dEcouldv2.z * v12.z / 2.0;
      m_stress[i2].x -= -dEcouldv2.x * v12.x / 2.0;
      m_stress[i2].y -= -dEcouldv2.y * v12.y / 2.0;
      m_stress[i2].z -= -dEcouldv2.z * v12.z / 2.0;
      d_stress[i1].x -= -dEcouldv2.x * v12.y / 2.0;
      d_stress[i1].y -= -dEcouldv2.y * v12.z / 2.0;
      d_stress[i1].z -= -dEcouldv2.z * v12.x / 2.0;
      d_stress[i2].x -= -dEcouldv2.x * v12.y / 2.0;
      d_stress[i2].y -= -dEcouldv2.y * v12.z / 2.0;
      d_stress[i2].z -= -dEcouldv2.z * v12.x / 2.0;
    }
  }
  return eVDW + eCoulomb;
}


double ReaxFF::CalcEnergyPolarization()
{
  double ePolar = 0.0;

  int nAtom = m_conf.GetRAtomNum();
  for(int i = 0; i < nAtom; ++i) {
    const ReaxAtom& atom = m_atomsNeigh[i];
    double q   = atom.m_dCharge;
    int    it  = atom.m_iType;
    double chi = s_params.chi[it];
    double eta = s_params.eta[it];
    
    double e = (chi * q + eta * q * q)
      * EV_TO_KCALMOL;

    ePolar += e;
  }

  return ePolar;
}

void ReaxFF::MakeAtomsNeigh()
{
  const Vector3& vx = m_conf.m_lattice.m_vX;
  const Vector3& vy = m_conf.m_lattice.m_vY;
  const Vector3& vz = m_conf.m_lattice.m_vZ;

  double dXY = (vx % vy).Normalize() * vz;
  double dYZ = (vy % vz).Normalize() * vx;
  double dZX = (vz % vx).Normalize() * vy;
  double dMin = fmin(fmin(dXY, dYZ), dZX);

  int nNeighCell = (int)ceil(s_params.nonbond_cut / dMin);
  int nAtom      = m_conf.GetRAtomNum();

  m_atomsNeigh.clear();
  for(int i = 0; i < nAtom; ++i) {
    m_atomsNeigh.push_back(ReaxAtom(m_conf.m_atoms[i], i));
  }
  for(int ix = -nNeighCell; ix <= nNeighCell; ++ix) {
    for(int iy = -nNeighCell; iy <= nNeighCell; ++iy) {
      for(int iz = -nNeighCell; iz <= nNeighCell; ++iz) {
        if( ix == 0 && iy == 0 && iz == 0 ) continue;
        for(int i = 0; i < nAtom; ++i) {
          ReaxAtom atom = m_atomsNeigh[i];
          atom.m_vPos += ix * vx + iy * vy + iz * vz;
          atom.m_iX = ix;
          atom.m_iY = iy;
          atom.m_iZ = iz;
          m_atomsNeigh.push_back(atom);
        }
      }
    }
  }
}

void ReaxFF::MakeBonds()
{
  const double CUT_BOND_SQR = CUT_BOND * CUT_BOND;

  int nAtom = m_conf.GetRAtomNum();

  m_bonds.clear();
  m_viBondsOnAtom.clear();
  m_viBondsOnAtom.resize(nAtom);

  // List up bonds and calculate BO'
  int iBond = 0;
  for(int i1 = 0; i1 < nAtom; ++i1) {
    const Vector3& v1  = m_atomsNeigh[i1].m_vPos;
    for(int i2 = 0; i2 < m_atomsNeigh.size(); ++i2) {
      if( i1 == i2 ) continue;
      const Vector3& v2 = m_atomsNeigh[i2].m_vPos;
      if( CUT_BOND_SQR < (v1-v2).Norm2() ) continue;
      ReaxBond bond(i1, i2);
      BondOrderPrime(bond);
      if( 0.0 < bond.m_BO ) {
        m_bonds.push_back(bond);
        m_viBondsOnAtom[i1].push_back(iBond);
        ++iBond;
      }
    }
  }

  TotalBondOrder();

  for(int i = 0; i < m_bonds.size(); ++i) {
    CorrectBondOrder(m_bonds[i]);
  }

  TotalBondOrder();

  SetBondSymmetryIndex();
}

void ReaxFF::SetBondSymmetryIndex()
{
  int iXMin = 0;
  int iYMin = 0;
  int iZMin = 0;
  int iXMax = 0;
  int iYMax = 0;
  int iZMax = 0;
  for(int i = 0; i < m_atomsNeigh.size(); ++i) {
    const ReaxAtom& atom = m_atomsNeigh[i];
    iXMin = MIN(iXMin, atom.m_iX);
    iYMin = MIN(iYMin, atom.m_iY);
    iZMin = MIN(iZMin, atom.m_iZ);
    iXMax = MAX(iXMax, atom.m_iX);
    iYMax = MAX(iYMax, atom.m_iY);
    iZMax = MAX(iZMax, atom.m_iZ);
  }

  int nX = iXMax - iXMin + 1;
  int nY = iYMax - iYMin + 1;
  int nZ = iZMax - iZMin + 1;
  int nXYZ = nX * nY * nZ;

  int nAtom = m_conf.GetRAtomNum();
  vector<int> dictBond(nAtom*nAtom*nXYZ);
  
  for(int i = 0; i < m_bonds.size(); ++i) {
    int i1 = m_bonds[i].m_iAtom1;
    int i2 = m_bonds[i].m_iAtom2;
    const ReaxAtom& atom1 = m_atomsNeigh[i1];
    const ReaxAtom& atom2 = m_atomsNeigh[i2];
    int iOrg2 = atom2.m_iOrg;
    int iX    = atom2.m_iX - iXMin;
    int iY    = atom2.m_iY - iYMin;
    int iZ    = atom2.m_iZ - iZMin;
    int iXYZ  = (iX * nY + iY) * nZ + iZ;
    int idx   = (i1 * nAtom + iOrg2) * nXYZ + iXYZ;
    dictBond[idx] = i;
  }
  for(int i = 0; i < m_bonds.size(); ++i) {
    int i1 = m_bonds[i].m_iAtom1;
    int i2 = m_bonds[i].m_iAtom2;
    const ReaxAtom& atom1 = m_atomsNeigh[i1];
    const ReaxAtom& atom2 = m_atomsNeigh[i2];
    int iOrg2 = atom2.m_iOrg;
    int iX    = - atom2.m_iX - iXMin;
    int iY    = - atom2.m_iY - iYMin;
    int iZ    = - atom2.m_iZ - iZMin;
    int iXYZ  = (iX * nY + iY) * nZ + iZ;
    int idx   = (iOrg2 * nAtom + i1) * nXYZ + iXYZ;
    m_bonds[i].m_iSymmetry = dictBond[idx];
  }
}


void ReaxFF::BondOrderPrime(ReaxBond& bond) const
{
  bond.m_BO    = 0.0;
  bond.m_BO_Sg = 0.0;
  bond.m_BO_Pi = 0.0;
  bond.m_BO_PP = 0.0;
  if( m_conf.m_bUseForce ) {
    bond.m_dBOP_Sgdr = 0.0;
    bond.m_dBOP_Pidr = 0.0;
    bond.m_dBOP_PPdr = 0.0;
    bond.m_dBOPdr    = 0.0;
  }

  int i1 = bond.m_iAtom1;
  int i2 = bond.m_iAtom2;
  const ReaxAtom& atom1 = m_atomsNeigh[i1];
  const ReaxAtom& atom2 = m_atomsNeigh[i2];
  const Vector3& v1 = atom1.m_vPos;
  const Vector3& v2 = atom2.m_vPos;
  int it1 = atom1.m_iType;
  int it2 = atom2.m_iType;
  double rSg1 = s_params.rSg[it1];
  double rSg2 = s_params.rSg[it2];
  double rPi1 = s_params.rPi[it1];
  double rPi2 = s_params.rPi[it2];
  double rPP1 = s_params.rPP[it1];
  double rPP2 = s_params.rPP[it2];
  double pbo1 = s_params.pbo1[it1][it2];
  double pbo2 = s_params.pbo2[it1][it2];
  double pbo3 = s_params.pbo3[it1][it2];
  double pbo4 = s_params.pbo4[it1][it2];
  double pbo5 = s_params.pbo5[it1][it2];
  double pbo6 = s_params.pbo6[it1][it2];
  double rSg  = s_params.rSg2[it1][it2];
  double rPi  = s_params.rPi2[it1][it2];
  double rPP  = s_params.rPP2[it1][it2];
  if( rSg <= 0.0 ) rSg = (rSg1 + rSg2) * 0.5;
  if( rPi <= 0.0 ) rPi = (rPi1 + rPi2) * 0.5;
  if( rPP <= 0.0 ) rPP = (rPP1 + rPP2) * 0.5;

  double r     = (v1-v2).Norm();
  double BO_Sg = 0.0;
  double BO_Pi = 0.0;
  double BO_PP = 0.0;
  double powSg = 0.0;
  double powPi = 0.0;
  double powPP = 0.0;
  if( 0.0 < rSg1 && 0.0 < rSg2 ) {
    powSg = pow(r / rSg, pbo2);
    BO_Sg = exp( pbo1 * powSg ) * (1.0 + CUT_BO);
  }
  if( 0.0 < rPi1 && 0.0 < rPi2 ) {
    powPi = pow(r / rPi, pbo4);
    BO_Pi = exp( pbo3 * powPi );
  }
  if( 0.0 < rPP1 && 0.0 < rPP2 ) {
    powPP = pow(r / rPP, pbo6);
    BO_PP = exp( pbo5 * powPP );
  }

  double BOP = BO_Sg + BO_Pi + BO_PP;

  if( BOP < CUT_BO ) {
    return;
  }

  if( m_conf.m_bUseForce ) {
    double ri = 1.0 / r;
    if( 0.0 < rSg1 && 0.0 < rSg2 ) {
      bond.m_dBOP_Sgdr = BO_Sg * pbo1 * pbo2 * powSg * ri;
    }
    if( 0.0 < rPi1 && 0.0 < rPi2 ) {
      bond.m_dBOP_Pidr = BO_Pi * pbo3 * pbo4 * powPi * ri;
    }
    if( 0.0 < rPP1 && 0.0 < rPP2 ) {
      bond.m_dBOP_PPdr = BO_PP * pbo5 * pbo6 * powPP * ri;
    }
    bond.m_dBOPdr = 
      bond.m_dBOP_Sgdr + 
      bond.m_dBOP_Pidr + 
      bond.m_dBOP_PPdr;
  }

  BOP   -= CUT_BO; // is this OK?
  BO_Sg -= CUT_BO;

  bond.m_BO    = BOP;
  bond.m_BO_Sg = BO_Sg;
  bond.m_BO_Pi = BO_Pi;
  bond.m_BO_PP = BO_PP;

  return;
}

void ReaxFF::CorrectBondOrder(ReaxBond& bond) const
{
  int i1 = bond.m_iAtom1;
  int i2 = bond.m_iAtom2;
  const ReaxAtom& atom1 = m_atomsNeigh[i1];
  const ReaxAtom& atom2 = m_atomsNeigh[i2];
  int it1 = atom1.m_iType;
  int it2 = atom2.m_iType;

  double ovc    = s_params.ovc   [it1][it2];
  double v13cor = s_params.v13cor[it1][it2];

  if( ovc < 0.001 && v13cor < 0.001 ) {
    if( m_conf.m_bUseForce ) {
      bond.m_dBOdBOP     = 1.0;
      bond.m_dBOdDeltaP1 = 0.0;
      
      bond.m_dBO_PidBOP_Pi  = bond.m_BO_Pi;
      bond.m_dBO_PidBOP     = 0.0;
      bond.m_dBO_PidDeltaP1 = 0.0;
      
      bond.m_dBO_PPdBOP_PP  = bond.m_BO_PP;
      bond.m_dBO_PPdBOP     = 0.0;
      bond.m_dBO_PPdDeltaP1 = 0.0;
    }
    return;
  }

  double BOP    = bond.m_BO;
  double BOP_Sg = bond.m_BO_Sg;
  double BOP_Pi = bond.m_BO_Pi;
  double BOP_PP = bond.m_BO_PP;

  double f1 = 1.0;
  double f2 = 0.0;
  double f3 = 0.0;
  double val1f2f3i = 0.0;
  double val2f2f3i = 0.0;
  double pboc1 = 0.0;
  double val1 = 0.0;
  double val2 = 0.0;
  double expf2_1 = 0.0;
  double expf2_2 = 0.0;
  double expf3_1 = 0.0;
  double expf3_2 = 0.0;
  if( 0.001 <= ovc ) {
    val1    = s_params.val[it1];
    val2    = s_params.val[it2];
    double DeltaP1 = atom1.m_totalBO - val1;
    double DeltaP2 = atom2.m_totalBO - val2;
    pboc1 = s_params.pboc1;
    double pboc2 = s_params.pboc2;
    expf2_1 = exp(- pboc1 * DeltaP1);
    expf2_2 = exp(- pboc1 * DeltaP2);
    expf3_1 = exp(- pboc2 * DeltaP1);
    expf3_2 = exp(- pboc2 * DeltaP2);
    f2 = expf2_1 + expf2_2;
    f3 = expf3_1 + expf3_2;
    f3 = - log( 0.5 * f3 ) / pboc2;
    val1f2f3i = 1.0 / ( val1 + f2 + f3 );
    val2f2f3i = 1.0 / ( val2 + f2 + f3 );
    f1 = 0.5 * ( ( val1 + f2 ) * val1f2f3i +
                 ( val2 + f2 ) * val2f2f3i );
  }

  double f4 = 1.0;
  double f5 = 1.0;
  double pboc3 = 0.0;
  double pboc4 = 0.0;
  if( 0.001 <= v13cor ) {
    double pboc3_1 = s_params.pboc3[it1];
    double pboc3_2 = s_params.pboc3[it2];
    double pboc4_1 = s_params.pboc4[it1];
    double pboc4_2 = s_params.pboc4[it2];
    double pboc5_1 = s_params.pboc5[it1];
    double pboc5_2 = s_params.pboc5[it2];
    pboc3   = sqrt(pboc3_1 * pboc3_2);
    pboc4   = sqrt(pboc4_1 * pboc4_2);
    double pboc5   = sqrt(pboc5_1 * pboc5_2);
    double valval1 = s_params.valval[it1];
    double valval2 = s_params.valval[it2];
    double DeltaP_BOC1  = atom1.m_totalBO - valval1;
    double DeltaP_BOC2  = atom2.m_totalBO - valval2;
    double d = pboc4 * SQR(BOP);
    double expf4 = exp(- pboc3 * (d - DeltaP_BOC1 ) + pboc5);
    double expf5 = exp(- pboc3 * (d - DeltaP_BOC2 ) + pboc5);
    f4 = 1.0 / (1.0 + expf4);
    f5 = 1.0 / (1.0 + expf5);
  }

  double f1f4f5   = f1*f4*f5;
  double f1f1f4f5 = f1*f1f4f5;

  bond.m_BO    = BOP    * f1f4f5;
  bond.m_BO_Pi = BOP_Pi * f1f1f4f5;
  bond.m_BO_PP = BOP_PP * f1f1f4f5;
  bond.m_BO_Sg = bond.m_BO - bond.m_BO_Pi - bond.m_BO_PP;

  // ignore tiny bonds
  if( bond.m_BO    < 1.0e-10 ) bond.m_BO    = 0.0;
  if( bond.m_BO_Sg < 1.0e-10 ) bond.m_BO_Sg = 0.0;
  if( bond.m_BO_Pi < 1.0e-10 ) bond.m_BO_Pi = 0.0;
  if( bond.m_BO_PP < 1.0e-10 ) bond.m_BO_PP = 0.0;
  
  if( !m_conf.m_bUseForce ) return;

  double df1dDeltaP1 = 0.0;
  double df1dDeltaP2 = 0.0;
  if( 0.001 <= ovc ) {
    double df1df2 =   0.5 * f3 * ( SQR(val1f2f3i) + 
                                   SQR(val2f2f3i) );
    double df1df3 = - 0.5 * ( (val1 + f2) * SQR(val1f2f3i) + 
                              (val2 + f2) * SQR(val2f2f3i) );
    double df2dDeltaP1 = - pboc1 * expf2_1;
    double df2dDeltaP2 = - pboc1 * expf2_2;
    double expf3i = 1.0 / (expf3_1 + expf3_2);
    double df3dDeltaP1 = expf3_1 * expf3i;
    double df3dDeltaP2 = expf3_2 * expf3i;
    df1dDeltaP1 = df1df2 * df2dDeltaP1 + df1df3 * df3dDeltaP1;
    df1dDeltaP2 = df1df2 * df2dDeltaP2 + df1df3 * df3dDeltaP2;
  }

  double df4dBOP = 0.0;
  double df5dBOP = 0.0;
  double df4dDeltaP1 = 0.0;
  double df5dDeltaP2 = 0.0;
  if( 0.001 <= v13cor ) {
    double d = 2.0 * pboc3 * pboc4 * BOP;
    df4dBOP = (f4 - f4 * f4) * d;
    df5dBOP = (f5 - f5 * f5) * d;
    df4dDeltaP1 = (f4 * f4 - f4) * pboc3;
    df5dDeltaP2 = (f5 * f5 - f5) * pboc3;
  }

  double df4f5dBOP = df4dBOP * f5 + f4 * df5dBOP;

  double df1f4f5dDeltaP1 = f5 
    * (df1dDeltaP1 * f4 + f1 * df4dDeltaP1);
  double df1f4f5dDeltaP2 = f4
    * (df1dDeltaP2 * f5 + f1 * df5dDeltaP2);

  double df1f1f4f5dDeltaP1 = df1dDeltaP1 * f1f4f5 
    + f1 * df1f4f5dDeltaP1;
  double df1f1f4f5dDeltaP2 = df1dDeltaP2 * f1f4f5 
    + f1 * df1f4f5dDeltaP2;

  bond.m_dBOdBOP     = f1f4f5 + BOP * f1 * df4f5dBOP;
  bond.m_dBOdDeltaP1 = BOP * df1f4f5dDeltaP1;

  bond.m_dBO_PidBOP_Pi  = f1f1f4f5;
  bond.m_dBO_PidBOP     = BOP_Pi * f1 * f1 * df4f5dBOP;
  bond.m_dBO_PidDeltaP1 = BOP_Pi * df1f1f4f5dDeltaP1;

  bond.m_dBO_PPdBOP_PP  = f1f1f4f5;
  bond.m_dBO_PPdBOP     = BOP_PP * f1 * f1 * df4f5dBOP;
  bond.m_dBO_PPdDeltaP1 = BOP_PP * df1f1f4f5dDeltaP1;
}

void ReaxFF::TotalBondOrder()
{
  // Calculate total bond order
  for(int i = 0; i < m_conf.GetRAtomNum(); ++i) {
    m_atomsNeigh[i].m_totalBO  = 0.0;
  }
  for(int i = 0; i < m_bonds.size(); ++i) {
    int i1 = m_bonds[i].m_iAtom1;
    m_atomsNeigh[i1].m_totalBO += m_bonds[i].m_BO;
  }
  for(int i = 0; i < m_atomsNeigh.size(); ++i) {
    int iOrg = m_atomsNeigh[i].m_iOrg;
    m_atomsNeigh[i].m_totalBO = m_atomsNeigh[iOrg].m_totalBO;
  }
}

// i2 might be outside of the box
// in that case i3 is connected to iOrg2
double ReaxFF::GetAngle(
  int i1, int i2, int i3, double& cost, double& sint,
  Vector3& dtdv1, Vector3& dtdv2, Vector3& dtdv3) const
{
  const ReaxAtom& atom1 = m_atomsNeigh[i1];
  const ReaxAtom& atom2 = m_atomsNeigh[i2];
  const ReaxAtom& atom3 = m_atomsNeigh[i3];
  Vector3 v1 = atom1.m_vPos;
  Vector3 v2 = atom2.m_vPos;
  Vector3 v3 = atom3.m_vPos;
  int iOrg2 = atom2.m_iOrg;
  if( iOrg2 != i2 ) {
    v3 += v2 - m_atomsNeigh[iOrg2].m_vPos;
  }
  v1 -= v2;
  v3 -= v2;
  double r1i = 1.0 / v1.Norm();
  double r3i = 1.0 / v3.Norm();
  Vector3 v1n = v1 * r1i;
  Vector3 v3n = v3 * r3i;
  cost = v1n * v3n;
  sint = (v1n % v3n).Norm();
  cost = fmin(cost,  1.0);
  cost = fmax(cost, -1.0);
  sint = fmin(sint,  1.0);
  sint = fmax(sint, -1.0);
  double Theta = atan2(sint, cost);

  if( m_conf.m_bUseForce ) {
    if( sint < 1.0e-10 ) {
      dtdv1 = VECTOR3_ZERO;
      dtdv2 = VECTOR3_ZERO;
      dtdv3 = VECTOR3_ZERO;
    }
    else {
      double sinti = 1.0 / sint;
      dtdv1 = - (v3n - cost * v1n) * r1i * sinti;
      dtdv3 = - (v1n - cost * v3n) * r3i * sinti;
      dtdv2 = - dtdv1 - dtdv3;
    }
  }

  return Theta;
}

// i3 might be outside of the box
// in that case i4 is connected to iOrg3
double ReaxFF::GetTorsion(
  int i1, int i2, int i3, int i4, double& coso, double& sino,
  Vector3& domegadv1, Vector3& domegadv2, 
  Vector3& domegadv3, Vector3& domegadv4 ) const
{
  double omega = PI;
  coso = -1.0;
  sino =  0.0;
  domegadv1 = VECTOR3_ZERO;
  domegadv2 = VECTOR3_ZERO;
  domegadv3 = VECTOR3_ZERO;
  domegadv4 = VECTOR3_ZERO;

  const ReaxAtom& atom1 = m_atomsNeigh[i1];
  const ReaxAtom& atom2 = m_atomsNeigh[i2];
  const ReaxAtom& atom3 = m_atomsNeigh[i3];
  const ReaxAtom& atom4 = m_atomsNeigh[i4];
  const Vector3& v1 = atom1.m_vPos;
  const Vector3& v2 = atom2.m_vPos;
  const Vector3& v3 = atom3.m_vPos;
  Vector3        v4 = atom4.m_vPos;
  int iOrg3 = atom3.m_iOrg;
  if( iOrg3 != i3 ) {
    v4 += v3 - m_atomsNeigh[iOrg3].m_vPos;
  }
  Vector3 v12 = v2 - v1;
  Vector3 v23 = v3 - v2;
  Vector3 v34 = v4 - v3;
  Vector3 v123 = v12 % v23;
  Vector3 v234 = v23 % v34;
  double d123 = v123.Norm();
  if( d123 == 0.0 ) return omega;
  double d234 = v234.Norm();
  if( d234 == 0.0 ) return omega;

  double d123i = 1.0 / d123;
  double d234i = 1.0 / d234;
  Vector3 v123n = v123 * d123i;
  Vector3 v234n = v234 * d234i;
  coso = v123n * v234n;
  sino = (v123n % v234n).Norm();
  coso = fmin(coso,  1.0);
  coso = fmax(coso, -1.0);
  sino = fmin(sino,  1.0);
  sino = fmax(sino, -1.0);
  omega = atan2(sino, coso);

  if( m_conf.m_bUseForce ) {
    if( sino < 1.0e-10 ) return omega;
    double sinoi = 1.0 / sino;
    
    Vector3 dcosodv1 = (v234n - coso * v123n) % v23 * d123i;
    domegadv1 = - dcosodv1 * sinoi;
    
    Vector3 dcosodv4 = (v123n - coso * v234n) % v23 * d234i;
    domegadv4 = - dcosodv4 * sinoi;

    Vector3 v13 = v3 - v1;
    Vector3 dcosodv2 = - (v234n - coso * v123n) % v13 * d123i
      + (v123n - coso * v234n) % v34 * d234i;
    domegadv2 = - dcosodv2 * sinoi;

    domegadv3 = - domegadv1 - domegadv2 - domegadv4;
  }
  
  return omega;
}

bool ReaxFF::IsPeriod12(int it) const
{
  return s_params.mass[it] < 21.0;
}

void ReaxFF::CalcChargeQEq()
{
  double rlow = s_params.nonbond_low;
  double rcut = s_params.nonbond_cut;
  double rcut_sqr = rcut * rcut;

  int nAtom  = m_conf.GetRAtomNum();
  int nNeigh = m_atomsNeigh.size();

  Matrix matH(nAtom, nAtom);
  for(int i1 = 0; i1 < nAtom; ++i1) {
    const ReaxAtom& atom1 = m_atomsNeigh[i1];
    const Vector3& v1 = atom1.m_vPos;
    int it1 = atom1.m_iType;
    double gamma1 = s_params.gamma[it1];
    double eta1   = s_params.eta  [it1];

    matH[i1][i1] = 2.0 * eta1;

    for(int i2 = 0; i2 < nNeigh; ++i2) {
      if( i1 == i2 ) continue;
      const ReaxAtom& atom2 = m_atomsNeigh[i2];
      const Vector3& v2 = atom2.m_vPos;

      double r_sqr = (v1-v2).Norm2();
      if( rcut_sqr < r_sqr ) continue;

      int it2 = atom2.m_iType;
      double gamma2 = s_params.gamma[it2];
      double g = 1.0 / sqrt(gamma1 * gamma2);

      double r = sqrt(r_sqr);
      double rtap = r - rlow;
      double tap = s_params.taps[7];
      for(int i = 0; i < 7; ++i) {
        tap = tap * rtap + s_params.taps[6-i];
      }
      
      double d = tap * 14.4 *
        pow(CUBE(r) + CUBE(g), - 1.0 / 3.0);
      
      matH[i1][atom2.m_iOrg] += d;
    }
  }

  Vector vs(nAtom);
  {
    Vector b(nAtom);
    for(int i = 0; i < nAtom; ++i) {
      int it = m_atomsNeigh[i].m_iType;
      b[i] = - s_params.chi[it];
    }
    matH.SolveLinearEquation(b, vs);
  }

  Vector vt(nAtom);
  {
    Vector b(nAtom, -1.0);
    matH.SolveLinearEquation(b, vt);
  }

  double ssum = 0.0;
  double tsum = 0.0;
  for(int i = 0; i < nAtom; ++i) {
    ssum += vs[i];
    tsum += vt[i];
  }
  Vector vq = vs - ssum / tsum * vt;

  for(int i = 0; i < nNeigh; ++i) {
    ReaxAtom& atom = m_atomsNeigh[i];
    atom.m_dCharge = vq[atom.m_iOrg];
  }  
}
