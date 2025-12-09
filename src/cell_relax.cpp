#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <iomanip>
#include "myheader.h"

#if !defined __linux__ && !defined __APPLE__
#define snprintf sprintf_s
#endif

void stretch(double x, double y, double z);
void strain(double strmat[3][3]);
void matcpy(double a[3][3], double b[3][3]);
extern int cellfix_xx, cellfix_yy, cellfix_zz, cellfix_xy, cellfix_yz, cellfix_zx;

void cell_relax_static_reset()
{
  cell.relax_static_initial = true;
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      cell.sgmmat_p[i][j] = cell.sgmmat[i][j];
      cell.hmat_p[i][j]   = cell.hmat[i][j];
    }
  }
  //for (int i=1; i<=atom.natom; i++) {
  //  atom.vx[i] = 0; atom.vy[i] = 0; atom.vz[i] = 0;
  //}
}

void cell_relax_static()
{
  double xstr, ystr, zstr, xystr, yzstr, zxstr;
  double dxstr, dystr, dzstr, dxystr, dyzstr, dzxstr;
  double x, y, z, dx, dy, dz, xy, yz, zx, dxy, dyz, dzx;
  double hmat_s[3][3], strmat[3][3];
  matcpy(cell.hmat,hmat_s);

  x   = cell.sgmmat[0][0]-cell.sgmmat_set[0][0];
  y   = cell.sgmmat[1][1]-cell.sgmmat_set[1][1];
  z   = cell.sgmmat[2][2]-cell.sgmmat_set[2][2];
  xy  = cell.sgmmat[0][1]-cell.sgmmat_set[0][1];
  yz  = cell.sgmmat[1][2]-cell.sgmmat_set[1][2];
  zx  = cell.sgmmat[2][0]-cell.sgmmat_set[2][0];
  dx  = x-(cell.sgmmat_p[0][0]-cell.sgmmat_set[0][0]);
  dy  = y-(cell.sgmmat_p[1][1]-cell.sgmmat_set[1][1]);
  dz  = z-(cell.sgmmat_p[2][2]-cell.sgmmat_set[2][2]);
  dxy = x-(cell.sgmmat_p[0][1]-cell.sgmmat_set[0][1]);
  dyz = x-(cell.sgmmat_p[1][2]-cell.sgmmat_set[1][2]);
  dzx = x-(cell.sgmmat_p[2][0]-cell.sgmmat_set[2][0]);

  if (cell.relax_static_initial) {
    xstr  =  -x*1.0e-6/1000000.;
    ystr  =  -y*1.0e-6/1000000.;
    zstr  =  -z*1.0e-6/1000000.;
    xystr = -xy*1.0e-6/1000000.;
    yzstr = -yz*1.0e-6/1000000.;
    zxstr = -zx*1.0e-6/1000000.;
    cell.relax_static_initial = false;
  } else {
    dxstr  = (cell.hmat[0][0]-cell.hmat_p[0][0])/cell.hmat[0][0];
    dystr  = (cell.hmat[1][1]-cell.hmat_p[1][1])/cell.hmat[1][1];
    dzstr  = (cell.hmat[2][2]-cell.hmat_p[2][2])/cell.hmat[2][2];
    dxystr = (cell.hmat[0][1]-cell.hmat_p[0][1])/cell.hmat[0][0];
    dyzstr = (cell.hmat[1][2]-cell.hmat_p[1][2])/cell.hmat[1][1]; 
    dzxstr = (cell.hmat[2][0]-cell.hmat_p[2][0])/cell.hmat[2][2];
    if (dx  != 0) { xstr  =    -x*dxstr/dx; } else { xstr  =  -x*1.0e-6/1000000.; }
    if (dy  != 0) { ystr  =    -y*dystr/dy; } else { ystr  =  -y*1.0e-6/1000000.; }
    if (dz  != 0) { zstr  =    -z*dzstr/dz; } else { zstr  =  -z*1.0e-6/1000000.; }
    if (dxy != 0) { xystr = -xy*dxystr/dxy; } else { xystr = -xy*1.0e-6/1000000.; }
    if (dyz != 0) { yzstr = -yz*dyzstr/dyz; } else { yzstr = -yz*1.0e-6/1000000.; }
    if (dzx != 0) { zxstr = -zx*dzxstr/dzx; } else { zxstr = -zx*1.0e-6/1000000.; }
  }

  strmat[0][0] =  xstr;
  strmat[1][1] =  ystr;
  strmat[2][2] =  zstr;
  strmat[0][1] = xystr; strmat[1][0] = xystr;
  strmat[1][2] = yzstr; strmat[2][1] = yzstr;
  strmat[2][0] = zxstr; strmat[0][2] = zxstr;

  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      cell.sgmmat_p[i][j] = cell.sgmmat[i][j];
      cell.hmat_p[i][j]   = cell.hmat[i][j];
      if (fabs(strmat[i][j]) > 0.01) { strmat[i][j] /= fabs(strmat[i][j])/0.01;
	cell.relax_static_initial = true;}
    }
  }

  //stretch(1.0+xstr, 1.0+ystr, 1.0+zstr);
  strain(strmat);

  if (cellfix_xx) { cell.hmat[0][0]=hmat_s[0][0]; }
  if (cellfix_yy) { cell.hmat[1][1]=hmat_s[1][1]; }
  if (cellfix_zz) { cell.hmat[2][2]=hmat_s[2][2]; }
  if (cellfix_xy) { cell.hmat[0][1]=hmat_s[0][1]; cell.hmat[1][0]=hmat_s[1][0]; }
  if (cellfix_yz) { cell.hmat[1][2]=hmat_s[1][2]; cell.hmat[2][1]=hmat_s[2][1]; }
  if (cellfix_zx) { cell.hmat[2][0]=hmat_s[2][0]; cell.hmat[0][2]=hmat_s[0][2]; }
}
