
#ifndef _POPULATION_H_
#define _POPULATION_H_

#include "ReaxPotential.h"
#include "Configuration.h"
#include "Input.h"
#include "Chromosome.h"

using namespace std;

class Population
{
public:
  Population() {}
  ~Population() {}

  void RandomInit(const ReaxPotential& pot);
  void CheckOpposite(const ReaxPotential& pot, bool bInit);

  const Chromosome& operator[](int i) const
  {
    return m_chroms[i];
  }
  Chromosome& operator[](int i)
  {
    return m_chroms[i];
  }
  int GetChromNum() const { return m_chroms.size(); }
  int GetChromSize() const 
  {
    if( m_chroms.empty() ) {
      return 0;
    }
    else {
      return m_chroms[0].size();
    }
  }

private:
  vector<Chromosome> m_chroms;
};

#endif // _POPULATION_H_

