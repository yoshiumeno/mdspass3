#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "myheader.h"


#define BL 1.380662e-23

void velocity_random(double target)
{
  double temp, factor;
  srand((unsigned int)time(NULL));
  printf("#Velocity is initialized by random distribution (T = %7.2f K)\n",(float)target);
  for (int i=1; i<=atom.natom; i++) {
    double x = 0, y = 0, z = 0;
    for (int k=1; k<=12; k++) { x += (double)rand()/(double)RAND_MAX; }
    for (int k=1; k<=12; k++) { y += (double)rand()/(double)RAND_MAX; }
    for (int k=1; k<=12; k++) { z += (double)rand()/(double)RAND_MAX; }
    atom.vx[i] = sqrt(BL*target/atom.wm[i])*(x-6.0);
    atom.vy[i] = sqrt(BL*target/atom.wm[i])*(y-6.0);
    atom.vz[i] = sqrt(BL*target/atom.wm[i])*(z-6.0);
    if (atom.mfx[i] == true) { atom.vx[i] = 0; }
    if (atom.mfy[i] == true) { atom.vy[i] = 0; }
    if (atom.mfz[i] == true) { atom.vz[i] = 0; }
    if (atom.QC==1) {
      if (atom.repatom[i] == 0) {
	atom.vx[i] = 0; atom.vy[i] = 0; atom.vz[i] = 0; } }
  }
}
  
/*
void vscale(int i, double target)
{
  double temp, factor;
  if (atom.QC==0) {
    temp = atom.Enkin()*2.0/3.0/atom.natom/BL;
  } else {
    temp = atom.Enkin()*2.0/3.0/atom.nrepatom/BL;
  }
  factor = sqrt(target/temp);
  atom.vx[i] = atom.vx[i] * factor;
  atom.vy[i] = atom.vy[i] * factor;
  atom.vz[i] = atom.vz[i] * factor;
}
*/
  
