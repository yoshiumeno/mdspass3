#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#if defined(_WIN32) || defined(_WIN64)
#include <Windows.h>
#endif
#include <GL/gl.h>
#include "myheader.h"

extern int select_atom[10], select_atom_repidx[10];

void get_repnum(int repidx, int &ix, int &iy, int &iz);

void calc_distance(int select_atom[], int select_atom_repidx[])
{
  for (int i=0; i<2; i++) { if (select_atom[i]==0) return; }
  int i0 = select_atom[0], i1 = select_atom[1];
  int ix0, iy0, iz0, ix1, iy1, iz1;
  get_repnum(select_atom_repidx[0],ix0,iy0,iz0);
  get_repnum(select_atom_repidx[1],ix1,iy1,iz1);
  double dist = atom.Dist2(i0, i1, ix1-ix0, iy1-iy0, iz1-iz0);
  dist = sqrt(dist)/ang;
  printf("Distance (Atom %d - %d) = %f (A)\n",i0, i1, dist);
}

void calc_angle(int select_atom[], int select_atom_repidx[])
{
  for (int i=0; i<3; i++) { if (select_atom[i]==0) return; }
  int i0 = select_atom[0], i1 = select_atom[1], i2 = select_atom[2];
  int ix0, iy0, iz0, ix1, iy1, iz1, ix2, iy2, iz2;
  get_repnum(select_atom_repidx[0],ix0,iy0,iz0);
  get_repnum(select_atom_repidx[1],ix1,iy1,iz1);
  get_repnum(select_atom_repidx[2],ix2,iy2,iz2);
  double angle = atom.Angle(i0, i1, i2, ix1-ix0, iy1-iy0, iz1-iz0,
			    ix2-ix0, iy2-ix0, iz2-iz0);
  printf("Angle (Atom %d - %d - %d) = %f (deg)\n",i0, i1, i2, angle);
}

void calc_dihedral(int select_atom[], int select_atom_repidx[])
{
  for (int i=0; i<4; i++) { if (select_atom[i]==0) return; }
  printf("Calculation of dihedral is not implemented yet.\n");
}

void get_repnum(int repidx, int &ix, int &iy, int &iz)
{
  if (repidx >= 4) { ix = 1; } else { ix = 0; }
  if ((repidx % 4 == 2)||(repidx % 4 == 3)) { iy = 1; } else { iy = 0; }
  if (repidx % 2 == 1) { iz = 1; } else { iz = 0; }
}
