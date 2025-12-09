#include <fstream>
#include <iostream>
#include "myheader.h"


#define BL 1.380662e-23

void vscale()
{
  double temp, factor;
  if (atom.QC==0) {
    temp = atom.Enkin()*2.0/3.0/(double)atom.natom/BL;
  } else {
    temp = atom.Enkin()*2.0/3.0/(double)atom.nrepatom/BL;
  }
  if (temp <= 0) { factor = 0;
  } else if (temp_set <= 0) { factor = 0;
  } else { factor = sqrt(temp_set/temp); }
  for (int i=1; i<=atom.natom; i++) {
    atom.vx[i] = atom.vx[i] * factor;
    atom.vy[i] = atom.vy[i] * factor;
    atom.vz[i] = atom.vz[i] * factor;
  }
}
  
void vscale(double target)
{
  double temp, factor;
  if (atom.QC==0) {
    temp = atom.Enkin()*2.0/3.0/atom.natom/BL;
  } else {
    temp = atom.Enkin()*2.0/3.0/atom.nrepatom/BL;
  }
  if (temp <= 0) { factor = 0;
  } else if (target <= 0) { factor = 0;
  } else { factor = sqrt(target/temp); }
  for (int i=1; i<=atom.natom; i++) {
    atom.vx[i] = atom.vx[i] * factor;
    atom.vy[i] = atom.vy[i] * factor;
    atom.vz[i] = atom.vz[i] * factor;
  }
}
  
void vscale(int i, double target)
{
  double temp, factor;
  if (atom.QC==0) {
    temp = atom.Enkin()*2.0/3.0/atom.natom/BL;
  } else {
    temp = atom.Enkin()*2.0/3.0/atom.nrepatom/BL;
  }
  if (temp <= 0) { factor = 0;
  } else if (target <= 0) { factor = 0;
  } else { factor = sqrt(target/temp); }
  atom.vx[i] = atom.vx[i] * factor;
  atom.vy[i] = atom.vy[i] * factor;
  atom.vz[i] = atom.vz[i] * factor;
}
  
