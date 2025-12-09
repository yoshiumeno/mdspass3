
#ifndef _INPUT_H_
#define _INPUT_H_

#include "String.h"

class Input
{
public:
  Input() { Init(); }
  ~Input() {}

  void Init();
  bool ReadFile(const String& strInputFile);

private:
  bool ReadFileMaster(const String& strInputFile);
  
public:
  int    m_nType;
  int    m_nRandomSeed;
  double m_dEnergyWeight;
  double m_dThreshold;
  int    m_nMaxIteration;
  double m_dPowellPrecision;
  String m_strConfigFile;
  String m_strReferenceFile;
  String m_strInitialPotentialFile;
  String m_strOutputPotentialFile;
};

#endif // _INPUT_H_

