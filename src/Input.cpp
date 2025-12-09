
#include <fstream>
#include <iostream>
#include "Input.h"
#include "Utils.h"

void Input::Init()
{
  m_nType                   = 0;
  m_nRandomSeed             = -1;
  m_dEnergyWeight           = 1.0;
  m_dThreshold              = 0.0;
  m_nMaxIteration           = 100;
  m_dPowellPrecision        = 0.0;
  m_strConfigFile           = "";
  m_strReferenceFile        = "";
  m_strInitialPotentialFile = "";
  m_strOutputPotentialFile  = "";
};

bool Input::ReadFile(const String& strInputFile)
{
  bool b;
  if( PFMPI::s_master ) {
    b = ReadFileMaster( strInputFile );
  }
  if( !PFMPI::Bcast(b, 1) ) return false;
  if( !b ) return false;

  if( !PFMPI::Bcast(m_nType           , 1) ) return false;
  if( !PFMPI::Bcast(m_nRandomSeed     , 1) ) return false;
  if( !PFMPI::Bcast(m_dEnergyWeight   , 1) ) return false;
  if( !PFMPI::Bcast(m_dThreshold      , 1) ) return false;
  if( !PFMPI::Bcast(m_nMaxIteration   , 1) ) return false;
  if( !PFMPI::Bcast(m_dPowellPrecision, 1) ) return false;

  return true;
}

bool Input::ReadFileMaster(const String& strInputFile)
{
  ifstream ifs( strInputFile.c_str() );
  if( !ifs ) {
    Utils::Error("Faield to open " + strInputFile);
    return false;
  }

  cout << "Reading " << strInputFile << " ... " << flush;

  Init();

  int iLine = 0;
  String line;
  while( getline(ifs, line) ) {
    ++iLine;

    String str = line.Trim();
    if( str.empty()   ) continue;
    if( str[0] == '#' ) continue;

    vector<String> tokens = str.Split();

    if( tokens.size() < 3 ) {
      Utils::Error("Insufficient number of tokens", iLine, line);
      return false;
    }

    if( tokens[1] != "=" ) {
      Utils::Error("Insufficient number of tokens", iLine, line);
      return false;
    }

    String key   = tokens[0].ToUpper();
    String value = tokens[2];

    try{
      if( key == "NTYPE" ) {
        m_nType = value.ToInt();
      }
      else if( key == "RANDOMSEED" ) {
        m_nRandomSeed = value.ToInt();
      }
      else if( key == "ENERGYWEIGHT" ) {
        m_dEnergyWeight = value.ToDouble();
      }
      else if( key == "THRESHOLD" ) {
        m_dThreshold = value.ToDouble();
      }
      else if( key == "MAXITERATION" ) {
        m_nMaxIteration = value.ToInt();
      }
      else if( key == "POWELLPRECISION" ) {
        m_dPowellPrecision = value.ToDouble();
      }
      else if( key == "CONFIGURATION" ) {
        m_strConfigFile = value;
      }
      else if( key == "REFERENCE" ) {
        m_strReferenceFile = value;
      }
      else if( key == "INITIALPOTENTIAL" ) {
        m_strInitialPotentialFile = value;
      }
      else if( key == "OUTPUTPOTENTIAL" ) {
        m_strOutputPotentialFile = value;
      }
      else {
        Utils::Error("Unknown keyword found", iLine, line);
        return false;
      }
    }
    catch(...) {
      Utils::Error("Failed to convert string to number", iLine, line);
      return false;
    }
  }
  
  if( m_strConfigFile == "" ) {
    Utils::Error("Specify a configuration file");
    return false;    
  }
  if( m_strReferenceFile == "" ) {
    Utils::Error("Specify a reference file");
    return false;    
  }
  if( m_strInitialPotentialFile == "" ) {
    Utils::Error("Specify a initial potential file");
    return false;    
  }
  if( m_strOutputPotentialFile == "" ) {
    Utils::Error("Specify a output potential file");
    return false;    
  }

  cout << "done." << endl;

  return true;
}
