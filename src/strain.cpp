#include <iostream>
#include "myheader.h"

void transpose(double a[3][3]);
void matmul(double a[3][3], double b[3][3], double c[3][3]);

void strain(double strmat[3][3])
{
  double xx, yy, zz;
  double hmat_p[3][3];
  for (int i=0; i<3; i++) { strmat[i][i] += 1.0; }
  for (int i=1; i<=atom.natom; i++) {
    xx = atom.rx[i]; yy = atom.ry[i]; zz = atom.rz[i];
    atom.rx[i] = xx*strmat[0][0]+yy*strmat[1][0]+zz*strmat[2][0];
    atom.ry[i] = xx*strmat[0][1]+yy*strmat[1][1]+zz*strmat[2][1];
    atom.rz[i] = xx*strmat[0][2]+yy*strmat[1][2]+zz*strmat[2][2]; }
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      hmat_p[i][j]=cell.hmat[i][j]; } }
  transpose(strmat);
  matmul(strmat,hmat_p,cell.hmat);
  for (int i=0; i<3; i++) { strmat[i][i] -= 1.0; }
}

