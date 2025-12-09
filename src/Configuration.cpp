
#include <iostream>
#include <fstream>
#include <string.h>
#include "Configuration.h"
#include "String.h"
#include "Utils.h"

using namespace std;

bool ConfigurationSet::ReadFile(const String& strFile)
{
  bool b;
  if( PFMPI::s_master ) {
    b = ReadFileMaster( strFile );
  }
  if( !PFMPI::Bcast(b, 1) ) return false;
  if( !b ) return false;

  int nConf = size();
  if( !PFMPI::Bcast(nConf, 1) ) return false;
  
  resize(nConf);

  for(int i = 0; i < nConf; ++i) {
    if( !(*this)[i].Bcast() ) return false;
  }

  return true;
}

bool Configuration::Bcast()
{
  int nType = m_atomTypes.size();
  int nAtom = m_atoms    .size();

  if( !PFMPI::Bcast(nType, 1) ) return false;
  if( !PFMPI::Bcast(nAtom, 1) ) return false;

  m_atomTypes.resize(nType);
  m_atoms    .resize(nAtom);
  m_force0   .resize(nAtom);

  for(int i = 0; i < nType; ++i) {
    if( !PFMPI::Bcast(m_atomTypes[i], 1) ) return false;
  }

  if( !PFMPI::Bcast(m_atoms[0] , nAtom) ) return false;
  if( !PFMPI::Bcast(m_bUseForce, 1    ) ) return false;
  if( !PFMPI::Bcast(m_dWeight  , 1    ) ) return false;
  if( !PFMPI::Bcast(m_dEnergy0 , 1    ) ) return false;
  if( !PFMPI::Bcast(m_force0[0], nAtom) ) return false;
  if( !PFMPI::Bcast(m_lattice  , 1    ) ) return false;

  return true;
}

bool ConfigurationSet::ModConf(int n, int ntyp,
			       double h00, double h01, double h02,
			       double h10, double h11, double h12,
			       double h20, double h21, double h22,
			       double *rx, double *ry, double *rz, char **asp)
{
  printf("!Modifying config info for ReaxFF...\n");
  clear();
  Configuration conf;
  conf.SetRAtomNum(n);
  conf.m_bUseForce = 1;
  conf.m_atomTypes.push_back("Zr");
  conf.m_atomTypes.push_back("O");
  conf.m_atomTypes.push_back("Ni");
  conf.m_lattice.m_vX.x = h00;  conf.m_lattice.m_vX.y = h10;  conf.m_lattice.m_vX.z = h20;
  conf.m_lattice.m_vY.x = h01;  conf.m_lattice.m_vY.y = h11;  conf.m_lattice.m_vY.z = h21;
  conf.m_lattice.m_vZ.x = h02;  conf.m_lattice.m_vZ.y = h12;  conf.m_lattice.m_vZ.z = h22;
  for (int i = 0; i < n; i++) {
    int ityp = -999;
    for (int j = 0; j < 3; j++ ) {
      if (strcmp(asp[i+1],"Zr")==0) { ityp = 0; }
      if (strcmp(asp[i+1],"O") ==0) { ityp = 1; }
      if (strcmp(asp[i+1],"Ni")==0) { ityp = 2; }
    }
    if (ityp < 0) {
      printf("There is something wrong in ModConf (atom type is not correctly set).\n");
      printf("I stop.\n");
      exit(0);
    }
    conf.m_atoms[i].m_iType = ityp;
    conf.m_atoms[i].m_vPos.x = rx[i+1];
    conf.m_atoms[i].m_vPos.y = ry[i+1];
    conf.m_atoms[i].m_vPos.z = rz[i+1];
    //printf(" %d %e %e %e %s\n",ityp,rx[i+1],ry[i+1],rz[i+1],asp[i+1]);
  }
  push_back(conf);
  
  bool b;
  if( !PFMPI::Bcast(b, 1) ) return false;
  if( !b ) return false;
  int nConf = size();
  if( !PFMPI::Bcast(nConf, 1) ) return false;
  resize(nConf);
  for(int i = 0; i < nConf; ++i) {
    if( !(*this)[i].Bcast() ) return false;
  }
  return true;
}
  
bool ConfigurationSet::ReadFileMaster(const String& strFile)
{
  ifstream ifs( strFile.c_str() );
  if( !ifs ) {
    Utils::Error("Failed to open " + strFile);
    return false;
  }

  cout << "Reading " << strFile << " ... " << flush;

  clear();

  int iLine = 0;
  String line;
  while( getline(ifs, line) ) {
    ++iLine;

    if( line[0] != '#' || line[1] != 'N') continue;
      
    String str = line.substr(2);
    vector<String> tokens = str.Split();
    int nToken = tokens.size();

    if (nToken < 2) {
      Utils::Error("Insufficient number of tokens",
                   iLine, line);
      return false;
    }

    Configuration conf;

    try {
      int nAtom     = tokens[0].ToInt();
      int iUseForce = tokens[1].ToInt();
      conf.SetRAtomNum(nAtom);
      conf.m_bUseForce = (iUseForce != 0);
      printf("READING FILE: NoA = %d\n",nAtom);
    }
    catch(...) {
      Utils::Error("Failed to read number of atoms",
                   iLine, line);
      return false;
    }
    
    bool bHasX = false;
    bool bHasY = false;
    bool bHasZ = false;

    while( getline(ifs, line) ) {
      ++iLine;

      if( line[0] != '#' ) continue;
      if( line.length() < 2 ) continue;

      if( line[1] == 'F' ) break;

      String str = line.substr(3);
      vector<String> tokens = str.Split();
      int nToken = tokens.size();

      String strC = line.substr(1,1).ToUpper();

      try {
        if( strC == "C" ) {
          for(int i = 0; i < nToken; ++i) {
            conf.m_atomTypes.push_back(tokens[i]);
          }
        }
        if( strC == "X" ) {
          if( nToken < 3 ) {
            Utils::Error("Insufficient number of tokens",
                         iLine, line);
            return false;
          }
          Vector3& v = conf.m_lattice.m_vX;
          v.x = tokens[0].ToDouble();
          v.y = tokens[1].ToDouble();
          v.z = tokens[2].ToDouble();
          if( v.y != 0.0 ) {
            Utils::Error(
              "Y component of X lattice vector must be zero",
              iLine, line);
            return false;
          }
          if( v.z != 0.0 ) {
            Utils::Error(
              "Z component of X lattice vector must be zero",
              iLine, line);
            return false;
          }
          bHasX = true;
        }
        else if( strC == "Y" ) {
          if( nToken < 3 ) {
            Utils::Error("Insufficient number of tokens",
                         iLine, line);
            return false;
          }
          Vector3& v = conf.m_lattice.m_vY;
          v.x = tokens[0].ToDouble();
          v.y = tokens[1].ToDouble();
          v.z = tokens[2].ToDouble();
          if( v.z != 0.0 ) {
            Utils::Error(
              "Z component of Y lattice vector must be zero",
              iLine, line);
            return false;
          }
          bHasY = true;
        }
        else if( strC == "Z" ) {
          if( nToken < 3 ) {
            Utils::Error("Insufficient number of tokens",
                         iLine, line);
            return false;
          }
          Vector3& v = conf.m_lattice.m_vZ;
          v.x = tokens[0].ToDouble();
          v.y = tokens[1].ToDouble();
          v.z = tokens[2].ToDouble();
          bHasZ = true;
        }
        else if( strC == "W" ) {
          if( nToken < 1 ) {
            Utils::Error("Insufficient number of tokens",
                         iLine, line);
            return false;
          }
          conf.m_dWeight = tokens[0].ToDouble();
        }
        else if( strC == "E" ) {
          if( nToken < 1 ) {
            Utils::Error("Insufficient number of tokens",
                         iLine, line);
            return false;
          }
          conf.m_dEnergy0 = tokens[0].ToDouble();
        }

      }
      catch(...) {
        Utils::Error("Failed to convert string to number",
                     iLine, line);
        return false;
      }
    }

    if( !(bHasX && bHasY && bHasZ) ) {
      Utils::Error("Lattice vectors are missing");
      return false;
    }

    for (int i = 0; i < conf.GetRAtomNum(); ++i)
    {
      if( !getline(ifs, line) ) {
        Utils::Error("Information for an atom missing",
                     iLine + 1);
        return false;
      }
      ++iLine;

      String str = line;
      vector<String> tokens = str.Split();
      int nToken = tokens.size();

      if( nToken < 7 ) {
        Utils::Error("Insufficiet number of tokens for an atom",
                     iLine, line);
        return false;
      }

      int iType;
      Vector3 v;
      try {
        conf.m_atoms[i].m_iType   = tokens[0].ToInt();
        conf.m_atoms[i].m_vPos.x  = tokens[1].ToDouble();
        conf.m_atoms[i].m_vPos.y  = tokens[2].ToDouble();
        conf.m_atoms[i].m_vPos.z  = tokens[3].ToDouble();
        conf.m_force0[i].x        = tokens[4].ToDouble();
        conf.m_force0[i].y        = tokens[5].ToDouble();
        conf.m_force0[i].z        = tokens[6].ToDouble();
        if( nToken < 8 ) continue;
        conf.m_atoms[i].m_dCharge = tokens[7].ToDouble();
        conf.m_bHasCharge = true;
      }
      catch(...){
        Utils::Error("Failed to convert a string to number",
                     iLine, line);
        return false;        
      }
    }

    push_back(conf);
  }

  cout << "done." << endl;

  return true;
}


