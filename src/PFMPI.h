
#ifndef _PFMPI_H_
#define _PFMPI_H_

#include <vector>

using namespace std;

class PFMPI
{
public:
  static bool Init(int* argc, char*** argv);
  static bool Finalize();
  template <class T> static bool Bcast(T& buffer, int count);
  template <class T> static bool BcastVector(vector<T>& v);
  static bool AllreduceSum(double& d);
  static bool AllgathervVector(vector<double>& v);

public:
  static int  s_rank;
  static int  s_size;
  static bool s_master;
};

#endif // _PFMPI_H_


