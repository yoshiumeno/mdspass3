
#ifndef _CONFIGURATION_H_
#define _CONFIGURATION_H_

#include <vector>
#include "Vector3.h"
#include "String.h"

using namespace std;

class Lattice
{
public:
  Lattice() {}
  Lattice(const Vector3 vX, const Vector3 vY, const Vector3 vZ) 
    : m_vX(vX), m_vY(vY), m_vZ(vZ) {}
  ~Lattice() {}

  Vector3 m_vX;
  Vector3 m_vY;
  Vector3 m_vZ;
};

class RAtom
{
public:
  RAtom() : m_iType(-1) {}
  RAtom(int i, const Vector3 v) 
    : m_iType(i), m_vPos(v) {}
  ~RAtom() {}

  int     m_iType;
  Vector3 m_vPos;
  double  m_dCharge;
};

class Configuration
{
public:
  Configuration() : m_bUseForce(false), m_dWeight(1.0),
                    m_dEnergy0(0.0), m_bHasCharge(false) {}
  ~Configuration() {}

  void SetRAtomNum(int n) 
  {
    m_atoms.resize(n);
    m_force0.resize(n);
  };

  int GetRAtomNum() const 
  {
    return m_atoms.size();
  };

  bool Bcast();

public:
  vector<String> m_atomTypes;
  vector<RAtom> m_atoms;
  bool m_bUseForce;
  double m_dWeight;
  double m_dEnergy0;        // reference energy
  vector<Vector3> m_force0; // reference force
  Lattice m_lattice;
  bool m_bHasCharge;
};

class ConfigurationSet : public vector<Configuration>
{
public:
  ConfigurationSet() {}
  ConfigurationSet(int n) 
    : vector<Configuration>(n) {}
  ~ConfigurationSet() {}
 
  bool ReadFile(const String& strFile);
  bool ModConf(int n, int ntyp,
	       double h11, double h12, double h13,
	       double h21, double h22, double h23,
	       double h31, double h32, double h33,
	       double *rx, double *ry, double *rz, char **asp);

private:
  bool ReadFileMaster(const String& strFile);
};

#endif // _CONFIGURATION_H_

