/////////////////////
// Functions for QC
/////////////////////
#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "myheader.h"

void mk_helmat(int iel, double rx[], double ry[], double rz[], double mat[][3])
{
  int ix, iy, iz, ix0, iy0, iz0;

  if (atom.elem_v_rep[iel][1] >= 4) { ix0 = 1; } else { ix0 = 0; }
  if ((atom.elem_v_rep[iel][1] % 4 == 2)||(atom.elem_v_rep[iel][1] % 4 == 3)) { iy0 = 1; } else { iy0 = 0; }
  if (atom.elem_v_rep[iel][1] % 2 == 1) { iz0 = 1; } else { iz0 = 0; }

  if (atom.elem_v_rep[iel][2] >= 4) { ix = 1; } else { ix = 0; }
  if ((atom.elem_v_rep[iel][2] % 4 == 2)||(atom.elem_v_rep[iel][2] % 4 == 3)) { iy = 1; } else { iy = 0; }
  if (atom.elem_v_rep[iel][2] % 2 == 1) { iz = 1; } else { iz = 0; }
  ix = ix-ix0; iy = iy-iy0; iz = iz-iz0;
  mat[0][0] = rx[atom.elem_v[iel][2]] - rx[atom.elem_v[iel][1]]
    + cell.hmat[0][0]*ix + cell.hmat[0][1]*iy + cell.hmat[0][2]*iz;
  mat[1][0] = ry[atom.elem_v[iel][2]] - ry[atom.elem_v[iel][1]]
    + cell.hmat[1][0]*ix + cell.hmat[1][1]*iy + cell.hmat[1][2]*iz;
  mat[2][0] = rz[atom.elem_v[iel][2]] - rz[atom.elem_v[iel][1]]
    + cell.hmat[2][0]*ix + cell.hmat[2][1]*iy + cell.hmat[2][2]*iz;

  if (atom.elem_v_rep[iel][3] >= 4) { ix = 1; } else { ix = 0; }
  if ((atom.elem_v_rep[iel][3] % 4 == 2)||(atom.elem_v_rep[iel][3] % 4 == 3)) { iy = 1; } else { iy = 0; }
  if (atom.elem_v_rep[iel][3] % 2 == 1) { iz = 1; } else { iz = 0; }
  ix = ix-ix0; iy = iy-iy0; iz = iz-iz0;
  mat[0][1] = rx[atom.elem_v[iel][3]] - rx[atom.elem_v[iel][1]]
    + cell.hmat[0][0]*ix + cell.hmat[0][1]*iy + cell.hmat[0][2]*iz;
  mat[1][1] = ry[atom.elem_v[iel][3]] - ry[atom.elem_v[iel][1]]
    + cell.hmat[1][0]*ix + cell.hmat[1][1]*iy + cell.hmat[1][2]*iz;
  mat[2][1] = rz[atom.elem_v[iel][3]] - rz[atom.elem_v[iel][1]]
    + cell.hmat[2][0]*ix + cell.hmat[2][1]*iy + cell.hmat[2][2]*iz;

  if (atom.elem_v_rep[iel][4] >= 4) { ix = 1; } else { ix = 0; }
  if ((atom.elem_v_rep[iel][4] % 4 == 2)||(atom.elem_v_rep[iel][4] % 4 == 3)) { iy = 1; } else { iy = 0; }
  if (atom.elem_v_rep[iel][4] % 2 == 1) { iz = 1; } else { iz = 0; }
  ix = ix-ix0; iy = iy-iy0; iz = iz-iz0;
  mat[0][2] = rx[atom.elem_v[iel][4]] - rx[atom.elem_v[iel][1]]
    + cell.hmat[0][0]*ix + cell.hmat[0][1]*iy + cell.hmat[0][2]*iz;
  mat[1][2] = ry[atom.elem_v[iel][4]] - ry[atom.elem_v[iel][1]]
    + cell.hmat[1][0]*ix + cell.hmat[1][1]*iy + cell.hmat[1][2]*iz;
  mat[2][2] = rz[atom.elem_v[iel][4]] - rz[atom.elem_v[iel][1]]
    + cell.hmat[2][0]*ix + cell.hmat[2][1]*iy + cell.hmat[2][2]*iz;
}
double q_x(int i, int iel, double rx[], double ry[], double rz[], double hinelmat[][3])
{
  double x0, y0, z0; int ix0, iy0, iz0;
  x0=rx[atom.elem_v[iel][1]]; y0=ry[atom.elem_v[iel][1]]; z0=rz[atom.elem_v[iel][1]];
  if (atom.elem_v_rep[iel][1] >= 4) { ix0 = 1; } else { ix0 = 0; }
  if ((atom.elem_v_rep[iel][1] % 4 == 2)||(atom.elem_v_rep[iel][1] % 4 == 3)) { iy0 = 1; } else { iy0 = 0; }
  if (atom.elem_v_rep[iel][1] % 2 == 1) { iz0 = 1; } else { iz0 = 0; }
  x0=x0+cell.hmat[0][0]*ix0+cell.hmat[0][1]*iy0+cell.hmat[0][2]*iz0;
  y0=y0+cell.hmat[1][0]*ix0+cell.hmat[1][1]*iy0+cell.hmat[1][2]*iz0;
  z0=z0+cell.hmat[2][0]*ix0+cell.hmat[2][1]*iy0+cell.hmat[2][2]*iz0;
  double xx = rx[i] - x0; double yy = ry[i] - y0; double zz = rz[i] - z0;
  return hinelmat[0][0]*xx + hinelmat[0][1]*yy + hinelmat[0][2]*zz;
}
double q_y(int i, int iel, double rx[], double ry[], double rz[], double hinelmat[][3])
{
  double x0, y0, z0; int ix0, iy0, iz0;
  x0=rx[atom.elem_v[iel][1]]; y0=ry[atom.elem_v[iel][1]]; z0=rz[atom.elem_v[iel][1]];
  if (atom.elem_v_rep[iel][1] >= 4) { ix0 = 1; } else { ix0 = 0; }
  if ((atom.elem_v_rep[iel][1] % 4 == 2)||(atom.elem_v_rep[iel][1] % 4 == 3)) { iy0 = 1; } else { iy0 = 0; }
  if (atom.elem_v_rep[iel][1] % 2 == 1) { iz0 = 1; } else { iz0 = 0; }
  x0=x0+cell.hmat[0][0]*ix0+cell.hmat[0][1]*iy0+cell.hmat[0][2]*iz0;
  y0=y0+cell.hmat[1][0]*ix0+cell.hmat[1][1]*iy0+cell.hmat[1][2]*iz0;
  z0=z0+cell.hmat[2][0]*ix0+cell.hmat[2][1]*iy0+cell.hmat[2][2]*iz0;
  double xx = rx[i] - x0; double yy = ry[i] - y0; double zz = rz[i] - z0;
  return hinelmat[1][0]*xx + hinelmat[1][1]*yy + hinelmat[1][2]*zz;
}
double q_z(int i, int iel, double rx[], double ry[], double rz[], double hinelmat[][3])
{
  double x0, y0, z0; int ix0, iy0, iz0;
  x0=rx[atom.elem_v[iel][1]]; y0=ry[atom.elem_v[iel][1]]; z0=rz[atom.elem_v[iel][1]];
  if (atom.elem_v_rep[iel][1] >= 4) { ix0 = 1; } else { ix0 = 0; }
  if ((atom.elem_v_rep[iel][1] % 4 == 2)||(atom.elem_v_rep[iel][1] % 4 == 3)) { iy0 = 1; } else { iy0 = 0; }
  if (atom.elem_v_rep[iel][1] % 2 == 1) { iz0 = 1; } else { iz0 = 0; }
  x0=x0+cell.hmat[0][0]*ix0+cell.hmat[0][1]*iy0+cell.hmat[0][2]*iz0;
  y0=y0+cell.hmat[1][0]*ix0+cell.hmat[1][1]*iy0+cell.hmat[1][2]*iz0;
  z0=z0+cell.hmat[2][0]*ix0+cell.hmat[2][1]*iy0+cell.hmat[2][2]*iz0;
  double xx = rx[i] - x0; double yy = ry[i] - y0; double zz = rz[i] - z0;
  return hinelmat[2][0]*xx + hinelmat[2][1]*yy + hinelmat[2][2]*zz;
}
double r_x(int i, int iel, double rx[], double ry[], double rz[], double helmat[][3], 
	   double qx, double qy, double qz)
{
  double x0, y0, z0; int ix0, iy0, iz0;
  x0=rx[atom.elem_v[iel][1]]; y0=ry[atom.elem_v[iel][1]]; z0=rz[atom.elem_v[iel][1]];
  if (atom.elem_v_rep[iel][1] >= 4) { ix0 = 1; } else { ix0 = 0; }
  if ((atom.elem_v_rep[iel][1] % 4 == 2)||(atom.elem_v_rep[iel][1] % 4 == 3)) { iy0 = 1; } else { iy0 = 0; }
  if (atom.elem_v_rep[iel][1] % 2 == 1) { iz0 = 1; } else { iz0 = 0; }
  x0=x0+cell.hmat[0][0]*ix0+cell.hmat[0][1]*iy0+cell.hmat[0][2]*iz0;
  y0=y0+cell.hmat[1][0]*ix0+cell.hmat[1][1]*iy0+cell.hmat[1][2]*iz0;
  z0=z0+cell.hmat[2][0]*ix0+cell.hmat[2][1]*iy0+cell.hmat[2][2]*iz0;
  return x0 + helmat[0][0]*qx + helmat[0][1]*qy + helmat[0][2]*qz;
}
double r_y(int i, int iel, double rx[], double ry[], double rz[], double helmat[][3], 
	   double qx, double qy, double qz)
{
  double x0, y0, z0; int ix0, iy0, iz0;
  x0=rx[atom.elem_v[iel][1]]; y0=ry[atom.elem_v[iel][1]]; z0=rz[atom.elem_v[iel][1]];
  if (atom.elem_v_rep[iel][1] >= 4) { ix0 = 1; } else { ix0 = 0; }
  if ((atom.elem_v_rep[iel][1] % 4 == 2)||(atom.elem_v_rep[iel][1] % 4 == 3)) { iy0 = 1; } else { iy0 = 0; }
  if (atom.elem_v_rep[iel][1] % 2 == 1) { iz0 = 1; } else { iz0 = 0; }
  x0=x0+cell.hmat[0][0]*ix0+cell.hmat[0][1]*iy0+cell.hmat[0][2]*iz0;
  y0=y0+cell.hmat[1][0]*ix0+cell.hmat[1][1]*iy0+cell.hmat[1][2]*iz0;
  z0=z0+cell.hmat[2][0]*ix0+cell.hmat[2][1]*iy0+cell.hmat[2][2]*iz0;
  return y0 + helmat[1][0]*qx + helmat[1][1]*qy + helmat[1][2]*qz;
}
double r_z(int i, int iel, double rx[], double ry[], double rz[], double helmat[][3], 
	   double qx, double qy, double qz)
{
  double x0, y0, z0; int ix0, iy0, iz0;
  x0=rx[atom.elem_v[iel][1]]; y0=ry[atom.elem_v[iel][1]]; z0=rz[atom.elem_v[iel][1]];
  if (atom.elem_v_rep[iel][1] >= 4) { ix0 = 1; } else { ix0 = 0; }
  if ((atom.elem_v_rep[iel][1] % 4 == 2)||(atom.elem_v_rep[iel][1] % 4 == 3)) { iy0 = 1; } else { iy0 = 0; }
  if (atom.elem_v_rep[iel][1] % 2 == 1) { iz0 = 1; } else { iz0 = 0; }
  x0=x0+cell.hmat[0][0]*ix0+cell.hmat[0][1]*iy0+cell.hmat[0][2]*iz0;
  y0=y0+cell.hmat[1][0]*ix0+cell.hmat[1][1]*iy0+cell.hmat[1][2]*iz0;
  z0=z0+cell.hmat[2][0]*ix0+cell.hmat[2][1]*iy0+cell.hmat[2][2]*iz0;
  return z0 + helmat[2][0]*qx + helmat[2][1]*qy + helmat[2][2]*qz;
}
