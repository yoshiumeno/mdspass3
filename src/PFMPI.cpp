
#ifdef _MPI_
#include <mpi.h>
#endif

#include <numeric>
#include "PFMPI.h"
#include "Configuration.h"
#include "ReaxPotential.h"

int  PFMPI::s_rank;
int  PFMPI::s_size;
bool PFMPI::s_master;

#ifdef _MPI_

bool PFMPI::Init(int* argc, char*** argv)
{
  int ierr = MPI_Init(argc, argv);
  if( ierr != 0 ) return  false;
  MPI_Comm_rank(MPI_COMM_WORLD, &s_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &s_size);
  s_master = (s_rank == 0);
  return true;
}

bool PFMPI::Finalize()
{
  return (MPI_Finalize() == 0);
}

template <class T>
bool PFMPI::Bcast(T& buffer, int count)
{
  count *= sizeof(T);
  int ierr = MPI_Bcast(&buffer, count, MPI_BYTE, 0, 
                       MPI_COMM_WORLD);
  return (ierr == 0);
}

template <>
bool PFMPI::Bcast(String& buffer, int count)
{
  int n = buffer.length();
  if( !Bcast(n, 1) ) return false;
  buffer.resize(n);
  return Bcast(buffer[0], n);
}

template <>
bool PFMPI::Bcast(PotentialParam& buffer, int count)
{
  if( !Bcast(buffer.m_strName, 1) ) return false;
  if( !Bcast(buffer.m_dValue , 1) ) return false;
  if( !Bcast(buffer.m_dMin   , 1) ) return false;
  if( !Bcast(buffer.m_dMax   , 1) ) return false;

  return true;
}

template <class T> 
bool PFMPI::BcastVector(vector<T>& v)
{
  int n = v.size();
  if( !Bcast(n, 1) ) return false;
  v.resize(n);
  for(int i = 0; i < n; ++i) {
    if( !Bcast(v[i], 1) ) return false;
  }
  return true;
}

bool PFMPI::AllreduceSum(double& d)
{
  double dRecv;
  int ierr = MPI_Allreduce(&d, &dRecv, 1, MPI_DOUBLE, MPI_SUM,
                           MPI_COMM_WORLD);
  d = dRecv;
  return (ierr == 0);
}

bool PFMPI::AllgathervVector(vector<double>& v)
{
  int scount = v.size();

  vector<int> rcount(s_size, 0);
  
  int ierr = MPI_Allgather(
    &scount     , 1, MPI_INT, 
    &(rcount[0]), 1, MPI_INT, MPI_COMM_WORLD);

  if( ierr != 0 ) return false;

  vector<int> displs(s_size);
  displs[0] = 0;
  for(int i = 1; i < s_size; ++i) {
    displs[i] = displs[i-1] + rcount[i-1];
  }

  int nSum = accumulate( rcount.begin(), rcount.end(), 0 );

  vector<double> vRecv(nSum, 0.0);
  
  ierr = MPI_Allgatherv(
    &(v[0]), scount, MPI_DOUBLE, &(vRecv[0]), &(rcount[0]),
    &(displs[0]), MPI_DOUBLE, MPI_COMM_WORLD);

  if( ierr != 0 ) return false;

  v = vRecv;
  
  return true;
}

#else

bool PFMPI::Init(int* argc, char*** argv)
{
  s_rank   = 0;
  s_size   = 1;
  s_master = true;
  return true;
}

bool PFMPI::Finalize() 
{
  return true;
}

template <class T> bool PFMPI::Bcast(T& buffer, int count)
{
  return true;
}
template bool PFMPI::Bcast(String& buffer, int count);

template <class T> bool PFMPI::BcastVector(vector<T>& v)
{
  return true;
}

bool PFMPI::AllreduceSum(double& d)
{
  return true;
}

bool PFMPI::AllgathervVector(vector<double>& v)
{
  return true;
}

#endif // _MPI_

template bool PFMPI::Bcast(int&     buffer, int count);
template bool PFMPI::Bcast(double&  buffer, int count);
template bool PFMPI::Bcast(bool&    buffer, int count);
template bool PFMPI::Bcast(RAtom&   buffer, int count);
template bool PFMPI::Bcast(Vector3& buffer, int count);
template bool PFMPI::Bcast(Lattice& buffer, int count);
template bool PFMPI::BcastVector(vector<String>& v);
template bool PFMPI::BcastVector(vector<PotentialParam>& v);

