#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#if defined(_WIN32) || defined(_WIN64)
#include <Windows.h>
#endif
#include <GL/gl.h>
#include "myheader.h"

double theta(double x, double y);

void rotate_atom(float x, float y, float z)
{
  double center_x = cell.hmat[0][0]/2.0;
  double center_y = cell.hmat[1][1]/2.0;
  double center_z = cell.hmat[2][2]/2.0;
  double xx, yy, zz, rr, t, tt;
  for (int i=1; i<=atom.natom; i++) {
    xx = atom.rx[i] - center_x;
    yy = atom.ry[i] - center_y;
    rr = sqrt(xx*xx+yy*yy);
    t  = theta(xx,yy);
    tt = t + (double)z /180.0 * M_PI;
    atom.rx[i] = center_x + rr * cos(tt);
    atom.ry[i] = center_y + rr * sin(tt);
  }
  for (int i=1; i<=atom.natom; i++) {
    yy = atom.ry[i] - center_y;
    zz = atom.rz[i] - center_z;
    rr = sqrt(yy*yy+zz*zz);
    t  = theta(yy,zz);
    tt = t + (double)x /180.0 * M_PI;
    atom.ry[i] = center_y + rr * cos(tt);
    atom.rz[i] = center_z + rr * sin(tt);
  }
  for (int i=1; i<=atom.natom; i++) {
    zz = atom.rz[i] - center_z;
    xx = atom.rx[i] - center_x;
    rr = sqrt(zz*zz+xx*xx);
    t  = theta(zz,xx);
    tt = t + (double)y /180.0 * M_PI;
    atom.rz[i] = center_z + rr * cos(tt);
    atom.rx[i] = center_x + rr * sin(tt);
  }
}

void rotate_atom(int i, int j, float x, float y, float z)
{
  double center_x = 0.0;
  double center_y = 0.0;
  double center_z = 0.0;
  if ((i>0)&&(i<=atom.natom)) {
    if ((j>0)&&(j<=atom.natom)) {
      if (j>=i) {
	for (int k=i; k<=j; k++) {
	  center_x += atom.rx[k];
	  center_y += atom.ry[k];
	  center_z += atom.rz[k];
	}
	center_x /= double(j-i+1);
	center_y /= double(j-i+1);
	center_z /= double(j-i+1);
      }
    }
  }
  double xx, yy, zz, rr, t, tt;
  if ((i>0)&&(i<=atom.natom)) {
    if ((j>0)&&(j<=atom.natom)) {
      if (j>=i) {
	for (int k=i; k<=j; k++) {
	  xx = atom.rx[k] - center_x;
	  yy = atom.ry[k] - center_y;
	  rr = sqrt(xx*xx+yy*yy);
	  t  = theta(xx,yy);
	  tt = t + (double)z /180.0 * M_PI;
	  atom.rx[k] = center_x + rr * cos(tt);
	  atom.ry[k] = center_y + rr * sin(tt);
	}
	for (int k=i; k<=j; k++) {
	  yy = atom.ry[k] - center_y;
	  zz = atom.rz[k] - center_z;
	  rr = sqrt(yy*yy+zz*zz);
	  t  = theta(yy,zz);
	  tt = t + (double)x /180.0 * M_PI;
	  atom.ry[k] = center_y + rr * cos(tt);
	  atom.rz[k] = center_z + rr * sin(tt);
	}
	for (int k=i; k<=j; k++) {
	  zz = atom.rz[k] - center_z;
	  xx = atom.rx[k] - center_x;
	  rr = sqrt(zz*zz+xx*xx);
	  t  = theta(zz,xx);
	  tt = t + (double)y /180.0 * M_PI;
	  atom.rz[k] = center_z + rr * cos(tt);
	  atom.rx[k] = center_x + rr * sin(tt);
	}
      }
    }
  }
}
