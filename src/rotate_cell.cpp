#include <iostream>
#include "myheader.h"

void matmul(double a[3][3], double b[3][3], double c[3][3]);
void inverse(double h[3][3], double hi[3][3]);
void matcpy(double a[3][3], double b[3][3]);

void rotate_cell(float x, float y, float z, int reverse)
{
  double rotx[3][3],roty[3][3],rotz[3][3],hnew[3][3],htmp[3][3];
  double xx,yy,zz,xxx,yyy,zzz;
  x = x/180*M_PI; y = y/180*M_PI; z = z/180*M_PI;
  if (reverse > 0) { x *= -1.0; y *= -1.0; z *= -1.0; }
  rotx[0][0] =     1.; rotx[0][1] =     0.; rotx[0][2] =     0.;
  rotx[1][0] =     0.; rotx[1][1] = cos(x); rotx[1][2] =-sin(x);
  rotx[2][0] =     0.; rotx[2][1] = sin(x); rotx[2][2] = cos(x);
  roty[0][0] = cos(y); roty[0][1] =     0.; roty[0][2] = sin(y);
  roty[1][0] =     0.; roty[1][1] =     1.; roty[1][2] =     0.;
  roty[2][0] =-sin(y); roty[2][1] =     0.; roty[2][2] = cos(y);
  rotz[0][0] = cos(z); rotz[0][1] =-sin(z); rotz[0][2] =     0.;
  rotz[1][0] = sin(z); rotz[1][1] = cos(z); rotz[1][2] =     0.;
  rotz[2][0] =     0.; rotz[2][1] =     0.; rotz[2][2] =     1.;
  inverse(cell.hmat,cell.hinmat);
  matmul(rotx,cell.hmat,hnew);
  matcpy(hnew,htmp);matmul(roty,htmp,hnew);
  matcpy(hnew,htmp);matmul(rotz,htmp,hnew);
  matcpy(hnew,cell.hmat);
  for (int i=1; i<=atom.natom; i++) {
    xx = atom.rx[i]; yy = atom.ry[i]; zz = atom.rz[i];
    xxx = cell.hinmat[0][0]*xx + cell.hinmat[0][1]*yy + cell.hinmat[0][2]*zz;
    yyy = cell.hinmat[1][0]*xx + cell.hinmat[1][1]*yy + cell.hinmat[1][2]*zz;
    zzz = cell.hinmat[2][0]*xx + cell.hinmat[2][1]*yy + cell.hinmat[2][2]*zz;
    atom.rx[i] = cell.hmat[0][0]*xxx + cell.hmat[0][1]*yyy + cell.hmat[0][2]*zzz;
    atom.ry[i] = cell.hmat[1][0]*xxx + cell.hmat[1][1]*yyy + cell.hmat[1][2]*zzz;
    atom.rz[i] = cell.hmat[2][0]*xxx + cell.hmat[2][1]*yyy + cell.hmat[2][2]*zzz;
  }
  inverse(cell.hmat,cell.hinmat);
  cellx=sqrt(cell.hmat[0][0]*cell.hmat[0][0]+cell.hmat[1][0]*cell.hmat[1][0]+cell.hmat[2][0]*cell.hmat[2][0])/ang;
  celly=sqrt(cell.hmat[0][1]*cell.hmat[0][1]+cell.hmat[1][1]*cell.hmat[1][1]+cell.hmat[2][1]*cell.hmat[2][1])/ang;
  cellz=sqrt(cell.hmat[0][2]*cell.hmat[0][2]+cell.hmat[1][2]*cell.hmat[1][2]+cell.hmat[2][2]*cell.hmat[2][2])/ang;
  cell1x=cell.hmat[0][0]/ang;cell1y=cell.hmat[1][0]/ang;cell1z=cell.hmat[2][0]/ang;
  cell2x=cell.hmat[0][1]/ang;cell2y=cell.hmat[1][1]/ang;cell2z=cell.hmat[2][1]/ang;
  cell3x=cell.hmat[0][2]/ang;cell3y=cell.hmat[1][2]/ang;cell3z=cell.hmat[2][2]/ang;
}

