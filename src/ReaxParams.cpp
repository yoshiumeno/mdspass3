
#include "ReaxParams.h"
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#define SET_PARAM(X, Y)                         \
  if( key == (X) ) { SetParam(s123, v, (Y)); }
#define SET_PARAM_HB(X, Y)                              \
  if( key == (X) ) { SetParam(s123, v, (Y), true); }

void ReaxParams::Setup(const ReaxPotential& pot, int nType)
{
  m_nType = nType;

  for(int i = 0; i < pot.GetParamNum(); ++i) {
    const PotentialParam& prm = pot[i];

    vector<String> key_i = prm.m_strName.Split("_");
    String key  = key_i[0];
    String s123;
    if( 1 < key_i.size() ) {
      s123 = key_i[1];
    }

    double v = prm.m_dValue;

    // Mass
    SET_PARAM("amas", mass);

    // Bond Order Prime
    SET_PARAM("rat" , rSg );
    SET_PARAM("rapt", rPi );
    SET_PARAM("vnq" , rPP );
    SET_PARAM("bop1", pbo1);
    SET_PARAM("bop2", pbo2);
    SET_PARAM("pdp" , pbo3);
    SET_PARAM("ptp" , pbo4);
    SET_PARAM("pdo" , pbo5);
    SET_PARAM("popi", pbo6);
    SET_PARAM("rsig", rSg2);
    SET_PARAM("rpi" , rPi2);
    SET_PARAM("rpi2", rPP2);

    // Correct Bond Order
    SET_PARAM("ovc"   , ovc   );
    SET_PARAM("v13cor", v13cor);
    SET_PARAM("aval"  , val   );
    SET_PARAM("vpar1" , pboc1 );
    SET_PARAM("vpar2" , pboc2 );
    SET_PARAM("bo132" , pboc3 );
    SET_PARAM("bo131" , pboc4 );
    SET_PARAM("bo133" , pboc5 );
    SET_PARAM("valf"  , valboc);

    // Bond
    SET_PARAM("psi", pbe1);
    SET_PARAM("psp", pbe2);
    SET_PARAM("de1", DeSg);
    SET_PARAM("de2", DePi);
    SET_PARAM("de3", DePP);
    
    // Lone Pair
    SET_PARAM("vpar16", plp1);
    SET_PARAM("vlp1"  , plp2);
    SET_PARAM("stlp"  , valE);

    // Overcoordination
    SET_PARAM("vover" , povun1);
    SET_PARAM("vovun" , povun2);
    SET_PARAM("vpar33", povun3);
    SET_PARAM("vpar32", povun4);

    // Undercoordination
    SET_PARAM("valp1" , povun5);
    SET_PARAM("vpar7" , povun6);
    SET_PARAM("vpar9" , povun7);
    SET_PARAM("vpar10", povun8);

    // Valence Angle
    SET_PARAM("vka"   , pval1  );
    SET_PARAM("vka3"  , pval2  );
    SET_PARAM("vval1" , pval3  );
    SET_PARAM("vval2" , pval4  );
    SET_PARAM("vval4" , pval5  );
    SET_PARAM("vpar15", pval6  );
    SET_PARAM("vkac"  , pval7  );
    SET_PARAM("vpar34", pval8  );
    SET_PARAM("vpar17", pval9  );
    SET_PARAM("vpar18", pval10 );
    SET_PARAM("valf"  , valang );
    SET_PARAM("th0"   , Theta00);

    // Penalty
    SET_PARAM("vkap"  , ppen1);
    SET_PARAM("vpar20", ppen2);
    SET_PARAM("vpar21", ppen3);
    SET_PARAM("vpar22", ppen4);

    // Coalition
    SET_PARAM("vka8"  , pcoa1);
    SET_PARAM("vpar3" , pcoa2);
    SET_PARAM("vpar39", pcoa3);
    SET_PARAM("vpar31", pcoa4);
    SET_PARAM("vval3" , valval);

    // Torsion
    SET_PARAM("v1"    , v1   );
    SET_PARAM("v2"    , v2   );
    SET_PARAM("v3"    , v3   );
    SET_PARAM("v4"    , ptor1);
    SET_PARAM("vpar24", ptor2);
    SET_PARAM("vpar25", ptor3);
    SET_PARAM("vpar26", ptor4);

    // Four Body Conjugate
    SET_PARAM("vconj" , pcot1);
    SET_PARAM("vpar28", pcot2);

    // Hydrogen Bond
    SET_PARAM   ("vnphb", phbond);
    SET_PARAM_HB("rhb"  , r0hb  );
    SET_PARAM_HB("dehb" , phb1  );
    SET_PARAM_HB("vhb1" , phb2  );
    SET_PARAM_HB("vhb2" , phb3  );

    // C2 Correction
    SET_PARAM("vpar6", kc2);

    // Triple Bond Correction
    SET_PARAM("vpar38", vtri );
    SET_PARAM("vpar11", ptri1);
    SET_PARAM("vpar8" , ptri2);
    SET_PARAM("vpar5" , ptri3);
    SET_PARAM("vpar4" , ptri4);

    // van del Waals
    SET_PARAM("vpar12", nonbond_low);
    SET_PARAM("vpar13", nonbond_cut);
    SET_PARAM("vpar29", pvdw1      );
    SET_PARAM("deodmh", dvdw       );
    SET_PARAM("rodmh" , rvdw       );
    SET_PARAM("godmh" , alpha      );
    SET_PARAM("eps"   , dvdw1      );
    SET_PARAM("rvdw"  , rvdw1      );
    SET_PARAM("alf"   , alpha1     );
    SET_PARAM("vop"   , gammaw     );

    // Coulomb
    SET_PARAM("gam", gamma);

    // Polarization
    SET_PARAM("chi", chi);
    SET_PARAM("eta", eta);
  }

  // valency_val must be same of valency_boc 
  // for period 1 and 2 elements
  for(int i = 0; i < m_nType; ++i) {
    if( mass[i] < 21.0 ) {
      valval[i] = valboc[i];
    }
  }

  // constants for toaper correction
  double rcut1 = nonbond_cut - nonbond_low;
  double rcut2 = rcut1 * rcut1;
  double rcut3 = rcut1 * rcut2;
  double rcut6 = rcut3 * rcut3;
  double rcut7 = rcut1 * rcut6;
  double rcut7i = 1.0 / rcut7;
  
  taps[7] =  20.0 * rcut7i;
  taps[6] = -70.0 * rcut7i * rcut1;
  taps[5] =  84.0 * rcut7i * rcut2;
  taps[4] = -35.0 * rcut7i * rcut3;
  taps[3] =  0.0;
  taps[2] =  0.0;
  taps[1] =  0.0;
  taps[0] =  1.0;
}

void ReaxParams::SetParam(
  const String& s123, double val, ReaxParam0D& param)
{
  param = val;
}

void ReaxParams::SetParam(
  const String& s123, double val, ReaxParam1D& param)
{
#ifdef _DEBUG_
  {
    assert( s123.length() == 1 );
    int i1 = s123.substr(0, 1).ToIntDef(0) - 1;
    assert(-1 < i1 && i1 < m_nType);
  }
#endif
  int i1 = s123.substr(0, 1).ToInt() - 1;
  param.Resize(m_nType);
  param[i1] = val;
}

void ReaxParams::SetParam(
  const String& s123, double val, ReaxParam2D& param)
{
#ifdef _DEBUG_
  {
    assert( s123.length() == 2 );
    int i1 = s123.substr(0, 1).ToIntDef(0) - 1;
    int i2 = s123.substr(1, 1).ToIntDef(0) - 1;
    assert(-1 < i1 && i1 < m_nType);
    assert(-1 < i2 && i2 < m_nType);
  }
#endif
  int i1 = s123.substr(0, 1).ToInt() - 1;
  int i2 = s123.substr(1, 1).ToInt() - 1;
  param.Resize(m_nType);
  param[i1][i2] = val;
  param[i2][i1] = val;
}

void ReaxParams::SetParam(
  const String& s123, double val, ReaxParam3D& param, bool bHB)
{
#ifdef _DEBUG_
  {
    assert( s123.length() == 3 );
    int i1 = s123.substr(0, 1).ToIntDef(0) - 1;
    int i2 = s123.substr(1, 1).ToIntDef(0) - 1;
    int i3 = s123.substr(2, 1).ToIntDef(0) - 1;
    assert(-1 < i1 && i1 < m_nType);
    assert(-1 < i2 && i2 < m_nType);
    assert(-1 < i3 && i3 < m_nType);
  }
#endif
  int i1 = s123.substr(0, 1).ToInt() - 1;
  int i2 = s123.substr(1, 1).ToInt() - 1;
  int i3 = s123.substr(2, 1).ToInt() - 1;
  param.Resize(m_nType);
  param[i1][i2][i3] = val;
  if( !bHB ) {
    param[i3][i2][i1] = val;
  }
}

void ReaxParams::SetParam(
  const String& s123, double val, ReaxParam4D& param)
{
#ifdef _DEBUG_
  {
    assert( s123.length() == 4 );
    int i1 = s123.substr(0, 1).ToIntDef(0) - 1;
    int i2 = s123.substr(1, 1).ToIntDef(0) - 1;
    int i3 = s123.substr(2, 1).ToIntDef(0) - 1;
    int i4 = s123.substr(3, 1).ToIntDef(0) - 1;
    assert(-1 < i1 && i1 < m_nType);
    assert(-1 < i2 && i2 < m_nType);
    assert(-1 < i3 && i3 < m_nType);
    assert(-1 < i4 && i4 < m_nType);
  }
#endif
  int i1 = s123.substr(0, 1).ToInt() - 1;
  int i2 = s123.substr(1, 1).ToInt() - 1;
  int i3 = s123.substr(2, 1).ToInt() - 1;
  int i4 = s123.substr(3, 1).ToInt() - 1;
  param.Resize(m_nType);
  param[i1][i2][i3][i4] = val;
  param[i4][i3][i2][i1] = val;
}

