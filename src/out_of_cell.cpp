#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#if defined(_WIN32) || defined(_WIN64)
#include <Windows.h>
#endif
#include <GL/gl.h>
#include "myheader.h"

void inverse(double mat[][3], double imat[][3]);

void out_of_cell()
{
  double hinmat[3][3];
  double a1, a2, a3;
  double xx, yy, zz;
  inverse(cell.hmat, hinmat);
  for (int i=1; i<=atom.natom; i++) {
    xx = hinmat[0][0]*atom.rx[i] + hinmat[0][1]*atom.ry[i] + hinmat[0][2]*atom.rz[i];
    yy = hinmat[1][0]*atom.rx[i] + hinmat[1][1]*atom.ry[i] + hinmat[1][2]*atom.rz[i];
    zz = hinmat[2][0]*atom.rx[i] + hinmat[2][1]*atom.ry[i] + hinmat[2][2]*atom.rz[i];
    a1 = 0.0; a2 = 0.0; a3 = 0.0;
    if (xx > 1.0) { a1 = -1.0; }    if (xx < 0.0) { a1 =  1.0; }
    if (yy > 1.0) { a2 = -1.0; }    if (yy < 0.0) { a2 =  1.0; }
    if (zz > 1.0) { a3 = -1.0; }    if (zz < 0.0) { a3 =  1.0; }
    if (cell.pbcx == 0) { a1 = 0.0; }
    if (cell.pbcy == 0) { a2 = 0.0; }
    if (cell.pbcz == 0) { a3 = 0.0; }
    atom.rx[i]=atom.rx[i]+cell.hmat[0][0]*a1+cell.hmat[0][1]*a2+cell.hmat[0][2]*a3;
    atom.ry[i]=atom.ry[i]+cell.hmat[1][0]*a1+cell.hmat[1][1]*a2+cell.hmat[1][2]*a3;
    atom.rz[i]=atom.rz[i]+cell.hmat[2][0]*a1+cell.hmat[2][1]*a2+cell.hmat[2][2]*a3;
    atom.rx_p[i]=atom.rx_p[i]+cell.hmat[0][0]*a1+cell.hmat[0][1]*a2+cell.hmat[0][2]*a3;
    atom.ry_p[i]=atom.ry_p[i]+cell.hmat[1][0]*a1+cell.hmat[1][1]*a2+cell.hmat[1][2]*a3;
    atom.rz_p[i]=atom.rz_p[i]+cell.hmat[2][0]*a1+cell.hmat[2][1]*a2+cell.hmat[2][2]*a3;
  }
}
