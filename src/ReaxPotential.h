
#ifndef _REAXPOTENTIAL_H_
#define _REAXPOTENTIAL_H_

#include <vector>
#include "String.h"

using namespace std;

class PotentialParam
{
public:
  String m_strName;
  double m_dValue;
  double m_dMin;
  double m_dMax;
};

class ReaxPotential
{
public:
  ReaxPotential() {}
  ReaxPotential(int n) : m_params(n) {}
  ~ReaxPotential() {}

  const PotentialParam& operator[](int i) const
  {
    return m_params[i];
  }
  PotentialParam& operator[](int i)
  {
    return m_params[i];
  }

  bool ReadFile(const string& strPotFile);
  bool WriteFile(const string& strPotFile) const;
  int GetVariableNum() const { return m_iVars.size(); }
  const vector<int>& GetVariableIndex() const { return m_iVars; }
  void SetVariableIndex(const vector<int>& viVars) 
  { m_iVars = viVars; }
  int GetParamNum() const { return m_params.size(); }
  vector<String> asp;

private:
  bool ReadFileMaster(const string& strPotFile);
  bool WriteFileMaster(const string& strPotFile) const;

private:
  vector<PotentialParam> m_params;
  vector<int>            m_iVars;
  vector<String>         m_header;
};

#endif // _REAXPOTENTIAL_H_

