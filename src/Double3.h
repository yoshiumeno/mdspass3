#ifndef _DOUBLE3_H_
#define _DOUBLE3_H_

class Double3
{
public:
  double& operator[](int i)
  {
    return m_d[i];
  }

  double m_d[3];
};

#endif // _DOUBLE3_H_

