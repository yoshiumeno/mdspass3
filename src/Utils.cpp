
#include "Utils.h"
#include <iostream>

void Utils::Error(const String& msg, int iLine, const String& line)
{
  if( !PFMPI::s_master ) return;

  cerr << endl;
  cerr << "ERROR) " << msg;
  if( -1 < iLine  ) {
    cerr << " on line " << iLine;
  }
  cerr << endl;
  if( line != "" ) {
    cerr << line << endl;
  }
  cerr << endl;
}
