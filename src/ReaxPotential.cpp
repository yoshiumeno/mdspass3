
#include <fstream>
#include <iostream>
#include <iomanip>
#include "ReaxPotential.h"
#include "Utils.h"

bool ReaxPotential::ReadFile(const string& strPotFile)
{
  bool b;
  if( PFMPI::s_master ) {
    b = ReadFileMaster( strPotFile );
  }
  if( !PFMPI::Bcast(b, 1) ) return false;
  if( !b ) return false;

  if( !PFMPI::BcastVector(m_header) ) return false;
  if( !PFMPI::BcastVector(m_params) ) return false;
  
  m_iVars.clear();
  for(int i = 0; i < m_params.size(); ++i) {
    if( m_params[i].m_dMin < m_params[i].m_dMax ) {
      m_iVars.push_back(i);
    }
  }

  return true;
}

bool ReaxPotential::ReadFileMaster(const string& strPotFile)
{
  ifstream ifs( strPotFile.c_str() );
  if( !ifs ) {
    Utils::Error("Failed to open " + strPotFile);
    return false;
  }

  cout << "Reading " << strPotFile << " ... " << flush;

  m_params.clear();

  int iLine = 0;
  String line;
  String str;
  while( getline(ifs, line) ) {
    ++iLine;

    // Addition by YU to extract atomtype info
    if (line[0] == '#' && line[1] == 'C') {
      printf("\n extracting supported atom species.. \n");
      String str = line.substr(2);
      vector<String> tokens = str.Split();
      int nToken = tokens.size();
      printf(" This ReaxFF supports %d atom species: ",nToken);
      asp.clear();
      for (int i = 0; i < nToken; i++) {
	printf(" %s ",tokens[i].c_str());
	asp.push_back(tokens[i]);
      }
      printf("\n");
    }
    //

    str = line.Trim().ToUpper();
    if( 0 < str.length() && str[0] == '#' ) {
      m_header.push_back(line);
    }
    if( str == "REAXFF" ) break;
  }

  if( str != "REAXFF" ) {
    Utils::Error("Could not find reaxff line");
    return false;
  }

  while( getline(ifs, line) ) {
    ++iLine;

    String str = line.Trim();
    if( str.empty() ) continue;
    vector<String> tokens = str.Split();

    if( tokens.size() < 4 ) {
      Utils::Error("Insufficient number of tokens", iLine, line);
      return false;
    }

    PotentialParam potpar;
    
    try {
      potpar.m_strName = tokens[0];
      potpar.m_dValue  = tokens[1].ToDouble();
      potpar.m_dMin    = tokens[2].ToDouble();
      potpar.m_dMax    = tokens[3].ToDouble();
    }
    catch(...) {
      Utils::Error("Failed to convert strings to number",
                   iLine, line);
      return false;
    }

    m_params.push_back(potpar);
  }

  ifs.close();

  cout << "done" << endl;

  return true;
}

bool ReaxPotential::WriteFile(const string& strPotFile) const
{
  bool b;
  if( PFMPI::s_master ) {
    b = WriteFileMaster( strPotFile );
  }
  if( !PFMPI::Bcast(b, 1) ) return false;
  return b;
}

bool ReaxPotential::WriteFileMaster(const string& strPotFile) const
{
  ofstream ofs( strPotFile.c_str() );
  if( !ofs ) {
    Utils::Error("Failed to open " + strPotFile);
    return false;
  }

  cout << "Writing " << strPotFile << " ... " << flush;

  for(int i = 0; i < m_header.size(); ++i) {
    ofs << m_header[i] << endl;
  }
  
  ofs << endl;
  ofs << "reaxff" << endl;

  for(int i = 0; i < m_params.size(); ++i) {
    const PotentialParam& potpar = m_params[i];
    ofs << fixed << setprecision(6)
        << potpar.m_strName << "\t "
        << potpar.m_dValue  << "\t "
        << potpar.m_dMin    << "\t "
        << potpar.m_dMax    << endl;
  }
  
  ofs.close();

  cout << "done" << endl;

  return true;
}

