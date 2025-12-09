#ifndef _REAX_PARAMS_H_
#define _REAX_PARAMS_H_

#include "ReaxPotential.h"

using namespace std;

typedef double ReaxParam0D;

class ReaxParam1D : public vector<ReaxParam0D>
{
public:
  ReaxParam1D() {}
  void Resize(int n) { resize(n, 0.0); }
};

class ReaxParam2D : public vector<ReaxParam1D>
{
public:
  ReaxParam2D() {}
  void Resize(int n) 
  {
    if( n == size() ) return;
    resize(n);
    for(int i = 0; i < n; ++i) {
      (*this)[i].Resize(n);
    }
  }
};

class ReaxParam3D : public vector<ReaxParam2D>
{
public:
  ReaxParam3D() {}
  void Resize(int n) 
  {
    if( n == size() ) return;
    resize(n);
    for(int i = 0; i < n; ++i) {
      (*this)[i].Resize(n);
    }
  }
};

class ReaxParam4D : public vector<ReaxParam3D>
{
public:
  ReaxParam4D() {}
  void Resize(int n) 
  {
    if( n == size() ) return;
    resize(n);
    for(int i = 0; i < n; ++i) {
      (*this)[i].Resize(n);
    }
  }
};

class ReaxParams
{
public:
  ReaxParams() {}
  ~ReaxParams() {}
  
  void Setup(const ReaxPotential& pot, int nType);

private:
  void SetParam(
    const String& s123, double val, ReaxParam0D& param);
  void SetParam(
    const String& s123, double val, ReaxParam1D& param);
  void SetParam(
    const String& s123, double val, ReaxParam2D& param);
  void SetParam(
    const String& s123, double val, ReaxParam3D& param, bool bHB = false);
  void SetParam(
    const String& s123, double val, ReaxParam4D& param);
  
public:
  int m_nType;
  
  // Mass
  ReaxParam1D mass;

  // Bond Order Prime
  ReaxParam1D rSg;
  ReaxParam1D rPi;
  ReaxParam1D rPP;
  ReaxParam2D pbo1;
  ReaxParam2D pbo2;
  ReaxParam2D pbo3;
  ReaxParam2D pbo4;
  ReaxParam2D pbo5;
  ReaxParam2D pbo6;
  ReaxParam2D rSg2;
  ReaxParam2D rPi2;
  ReaxParam2D rPP2;

  // Correct Bond Order
  ReaxParam0D pboc1;
  ReaxParam0D pboc2;
  ReaxParam1D pboc3;
  ReaxParam1D pboc4;
  ReaxParam1D pboc5;
  ReaxParam1D val;
  ReaxParam1D valboc;
  ReaxParam2D ovc;
  ReaxParam2D v13cor;

  // Bond
  ReaxParam2D pbe1;
  ReaxParam2D pbe2;
  ReaxParam2D DeSg;
  ReaxParam2D DePi;
  ReaxParam2D DePP;

  // Lone Pair
  ReaxParam0D plp1;
  ReaxParam1D plp2;
  ReaxParam1D valE;

  // Overcoordination
  ReaxParam2D povun1;
  ReaxParam1D povun2;
  ReaxParam0D povun3;
  ReaxParam0D povun4;

  // Undercoordination
  ReaxParam1D povun5;
  ReaxParam0D povun6;
  ReaxParam0D povun7;
  ReaxParam0D povun8;

  // Valence Angle
  ReaxParam3D pval1;
  ReaxParam3D pval2;
  ReaxParam1D pval3;
  ReaxParam3D pval4;
  ReaxParam1D pval5;
  ReaxParam0D pval6;
  ReaxParam3D pval7;
  ReaxParam0D pval8;
  ReaxParam0D pval9;
  ReaxParam0D pval10;
  ReaxParam1D valang;
  ReaxParam3D Theta00;

  // Penalty
  ReaxParam3D ppen1;
  ReaxParam0D ppen2;
  ReaxParam0D ppen3;
  ReaxParam0D ppen4;

  // Coalition
  ReaxParam3D pcoa1;
  ReaxParam0D pcoa2;
  ReaxParam0D pcoa3;
  ReaxParam0D pcoa4;
  ReaxParam1D valval;

  // Torsion
  ReaxParam4D v1;
  ReaxParam4D v2;
  ReaxParam4D v3;
  ReaxParam4D ptor1;
  ReaxParam0D ptor2;
  ReaxParam0D ptor3;
  ReaxParam0D ptor4;

  // Four Body Conjugate
  ReaxParam4D pcot1;
  ReaxParam0D pcot2;

  // Hydrogen Bond
  ReaxParam1D phbond;
  ReaxParam3D r0hb;
  ReaxParam3D phb1;
  ReaxParam3D phb2;
  ReaxParam3D phb3;
  
  // C2 Correction
  ReaxParam0D kc2;

  // Triple Bond Correction
  ReaxParam0D vtri;
  ReaxParam0D ptri1;
  ReaxParam0D ptri2;
  ReaxParam0D ptri3;
  ReaxParam0D ptri4;

  // van der Waals
  ReaxParam0D nonbond_low;
  ReaxParam0D nonbond_cut;
  ReaxParam0D pvdw1;
  ReaxParam2D dvdw;
  ReaxParam2D rvdw;
  ReaxParam2D alpha;
  ReaxParam1D dvdw1;
  ReaxParam1D rvdw1;
  ReaxParam1D gammaw1;
  ReaxParam1D alpha1;
  ReaxParam1D gammaw;
  double taps[8];

  // Coulomb
  ReaxParam1D gamma;

  // Polarization
  ReaxParam1D chi;
  ReaxParam1D eta;
};

#endif // _REAX_PARAMS_H_


