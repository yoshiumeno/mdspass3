#include <iostream>
#include "myheader.h"

void inverse(double mat[][3], double imat[][3]);

void stretch(double x, double y, double z)
{
  double xx,yy,zz,xxx,yyy,zzz;
  inverse(cell.hmat,cell.hinmat);
  cell.hmat[0][0]=cell.hmat[0][0]*x;
  cell.hmat[1][0]=cell.hmat[1][0]*x;
  cell.hmat[2][0]=cell.hmat[2][0]*x;
  cell.hmat[0][1]=cell.hmat[0][1]*y;
  cell.hmat[1][1]=cell.hmat[1][1]*y;
  cell.hmat[2][1]=cell.hmat[2][1]*y;
  cell.hmat[0][2]=cell.hmat[0][2]*z;
  cell.hmat[1][2]=cell.hmat[1][2]*z;
  cell.hmat[2][2]=cell.hmat[2][2]*z;
  for (int i=1; i<=atom.natom; i++) {
    xx=atom.rx[i]; yy=atom.ry[i]; zz=atom.rz[i];
    xxx = cell.hinmat[0][0]*xx + cell.hinmat[0][1]*yy + cell.hinmat[0][2]*zz;
    yyy = cell.hinmat[1][0]*xx + cell.hinmat[1][1]*yy + cell.hinmat[1][2]*zz;
    zzz = cell.hinmat[2][0]*xx + cell.hinmat[2][1]*yy + cell.hinmat[2][2]*zz;
    atom.rx[i] = cell.hmat[0][0]*xxx + cell.hmat[0][1]*yyy + cell.hmat[0][2]*zzz;
    atom.ry[i] = cell.hmat[1][0]*xxx + cell.hmat[1][1]*yyy + cell.hmat[1][2]*zzz;
    atom.rz[i] = cell.hmat[2][0]*xxx + cell.hmat[2][1]*yyy + cell.hmat[2][2]*zzz;
  }
  inverse(cell.hmat,cell.hinmat);
}

void stretch_celladjust(float x, float y, float z)
{
  double strx, stry, strz;
  double len1,len2,len3;
  len1=sqrt(cell.hmat[0][0]*cell.hmat[0][0]+cell.hmat[1][0]*cell.hmat[1][0]+cell.hmat[2][0]*cell.hmat[2][0])/ang;
  len2=sqrt(cell.hmat[0][1]*cell.hmat[0][1]+cell.hmat[1][1]*cell.hmat[1][1]+cell.hmat[2][1]*cell.hmat[2][1])/ang;
  len3=sqrt(cell.hmat[0][2]*cell.hmat[0][2]+cell.hmat[1][2]*cell.hmat[1][2]+cell.hmat[2][2]*cell.hmat[2][2])/ang;
  //strx = x*ang/cell.hmat[0][0];
  //stry = y*ang/cell.hmat[1][1];
  //strz = z*ang/cell.hmat[2][2];
  strx = x/len1;
  stry = y/len2;
  strz = z/len3;
  stretch(strx,stry,strz);
}

void stretch_celladjust(float c1x, float c1y, float c1z,
			float c2x, float c2y, float c2z,
			float c3x, float c3y, float c3z)
{
  double xx,yy,zz,xxx,yyy,zzz;
  inverse(cell.hmat,cell.hinmat);
  cell.hmat[0][0] = (double)c1x*ang;
  cell.hmat[1][0] = (double)c1y*ang;
  cell.hmat[2][0] = (double)c1z*ang;
  cell.hmat[0][1] = (double)c2x*ang;
  cell.hmat[1][1] = (double)c2y*ang;
  cell.hmat[2][1] = (double)c2z*ang;
  cell.hmat[0][2] = (double)c3x*ang;
  cell.hmat[1][2] = (double)c3y*ang;
  cell.hmat[2][2] = (double)c3z*ang;
  if (ifix_atoms == 0) {
    for (int i=1; i<=atom.natom; i++) {
      xx=atom.rx[i]; yy=atom.ry[i]; zz=atom.rz[i];
      xxx = cell.hinmat[0][0]*xx + cell.hinmat[0][1]*yy + cell.hinmat[0][2]*zz;
      yyy = cell.hinmat[1][0]*xx + cell.hinmat[1][1]*yy + cell.hinmat[1][2]*zz;
      zzz = cell.hinmat[2][0]*xx + cell.hinmat[2][1]*yy + cell.hinmat[2][2]*zz;
      atom.rx[i] = cell.hmat[0][0]*xxx + cell.hmat[0][1]*yyy + cell.hmat[0][2]*zzz;
      atom.ry[i] = cell.hmat[1][0]*xxx + cell.hmat[1][1]*yyy + cell.hmat[1][2]*zzz;
      atom.rz[i] = cell.hmat[2][0]*xxx + cell.hmat[2][1]*yyy + cell.hmat[2][2]*zzz;
    }
  }
  inverse(cell.hmat,cell.hinmat);
}
void stretch_celladjust_d(double c1x, double c1y, double c1z,
			double c2x, double c2y, double c2z,
			double c3x, double c3y, double c3z)
{
  double xx,yy,zz,xxx,yyy,zzz;
  inverse(cell.hmat,cell.hinmat);
  cell.hmat[0][0] = (double)c1x*ang;
  cell.hmat[1][0] = (double)c1y*ang;
  cell.hmat[2][0] = (double)c1z*ang;
  cell.hmat[0][1] = (double)c2x*ang;
  cell.hmat[1][1] = (double)c2y*ang;
  cell.hmat[2][1] = (double)c2z*ang;
  cell.hmat[0][2] = (double)c3x*ang;
  cell.hmat[1][2] = (double)c3y*ang;
  cell.hmat[2][2] = (double)c3z*ang;
  if (ifix_atoms == 0) {
    for (int i=1; i<=atom.natom; i++) {
      xx=atom.rx[i]; yy=atom.ry[i]; zz=atom.rz[i];
      xxx = cell.hinmat[0][0]*xx + cell.hinmat[0][1]*yy + cell.hinmat[0][2]*zz;
      yyy = cell.hinmat[1][0]*xx + cell.hinmat[1][1]*yy + cell.hinmat[1][2]*zz;
      zzz = cell.hinmat[2][0]*xx + cell.hinmat[2][1]*yy + cell.hinmat[2][2]*zz;
      atom.rx[i] = cell.hmat[0][0]*xxx + cell.hmat[0][1]*yyy + cell.hmat[0][2]*zzz;
      atom.ry[i] = cell.hmat[1][0]*xxx + cell.hmat[1][1]*yyy + cell.hmat[1][2]*zzz;
      atom.rz[i] = cell.hmat[2][0]*xxx + cell.hmat[2][1]*yyy + cell.hmat[2][2]*zzz;
    }
  }
  inverse(cell.hmat,cell.hinmat);
}

void celladjust(double hmat0[3][3], double hmat[3][3])
{
  double xx,yy,zz, qx,qy,qz;
  double hinmat0[3][3];
  inverse(hmat0,hinmat0);
  for (int i=1; i<=atom.natom; i++) {
    qx = hinmat0[0][0]*atom.rx[i] + hinmat0[0][1]*atom.ry[i] + hinmat0[0][2]*atom.rz[i];
    qy = hinmat0[1][0]*atom.rx[i] + hinmat0[1][1]*atom.ry[i] + hinmat0[1][2]*atom.rz[i];
    qz = hinmat0[2][0]*atom.rx[i] + hinmat0[2][1]*atom.ry[i] + hinmat0[2][2]*atom.rz[i];
    xx = hmat[0][0]*qx + hmat[0][1]*qy + hmat[0][2]*qz;
    yy = hmat[1][0]*qx + hmat[1][1]*qy + hmat[1][2]*qz;
    zz = hmat[2][0]*qx + hmat[2][1]*qy + hmat[2][2]*qz;
    atom.rx[i] = xx;
    atom.ry[i] = yy;
    atom.rz[i] = zz;
  }
}
