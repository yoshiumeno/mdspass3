#include <string.h>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "myheader.h"

void bookkeep();
void potential();
void loading();
void calculate();
void stretch_celladjust_d(double c1x, double c1y, double c1z,
			  double c2x, double c2y, double c2z,
			  double c3x, double c3y, double c3z);

void calc_center_xy(double &x, double &y);
void calc_cylinder_side_area(double center_x, double center_y, double &rr, double &area);
double theta(double x, double y);

extern float auto_val1, auto_val2;

void autocalc_o()
{
  double rangex = cell.hmat[0][0]/2;
  double rangey = cell.hmat[1][1]/3;
  printf("Auto calculation\n");
  for (int ii=0;ii<=10;ii++) {
  for (int jj=0;jj<=10;jj++) {
  double c1x, c1y, c1z, c2x, c2y, c2z, c3x, c3y, c3z;
  c1x = cell.hmat[0][0]/ang;
  c1y = cell.hmat[1][0]/ang;
  c1z = cell.hmat[2][0]/ang;
  c2x = cell.hmat[0][1]/ang;
  c2y = cell.hmat[1][1]/ang;
  c2z = cell.hmat[2][1]/ang;
  c3x = cell.hmat[0][2]/ang;
  c3y = cell.hmat[1][2]/ang;
  c3z = cell.hmat[2][2]/ang;

  //shift atoms
  double sftx = rangex / 10 * (double)ii;
  double sfty = rangey / 10 * (double)jj;
  for (int i=1; i<=atom.natom; i++) {
    double xx = atom.rz_org[i]/cell.hmat[2][2];
    if (xx > 0.3) {
      atom.rx[i] = atom.rx_org[i] + sftx;
      atom.ry[i] = atom.ry_org[i] + sfty;
    }
  }

  double ytarget = -100;
  double tol = 10.0;
  double xmin = (double)auto_val2 - 0.005;
  double xmax = (double)auto_val2 + 0.005;
  double xcnt;
  double ymin, ymax, ycnt;
  double ep;
  stretch_celladjust_d(c1x,c1y,c1z,c2x,c2y,c2z,c3x,c3y,xmin); 
  calculate(); ymin = cell.sgmmat[2][2]/MPa;
  stretch_celladjust_d(c1x,c1y,c1z,c2x,c2y,c2z,c3x,c3y,xmax);
  calculate(); ymax = cell.sgmmat[2][2]/MPa;
  //printf("Initial: x = %f - %f, y = %f - %f\n",xmin,xmax,ymin,ymax);
START:
  int min_or_max;
  for (int ite=0; ite<20; ite++) {
    xcnt = (xmin+xmax)/2.0;
    stretch_celladjust_d(c1x,c1y,c1z,c2x,c2y,c2z,c3x,c3y,xmin);
    calculate(); ymin = cell.sgmmat[2][2]/MPa;
    if (abs(ymin-ytarget) < tol) { min_or_max=0; break; }
    if (ymin > ytarget) {
      double xx = xmax-xmin;
      xmin -= xx; xmax -= xx; ite--; continue; }
    stretch_celladjust_d(c1x,c1y,c1z,c2x,c2y,c2z,c3x,c3y,xmax);
    calculate(); ymax = cell.sgmmat[2][2]/MPa;
    if (abs(ymax-ytarget) < tol) { min_or_max=1; break; }
    if (ymax < ytarget) {
      double xx = xmax-xmin;
      xmax += xx; xmin += xx; ite--; continue; }
    stretch_celladjust_d(c1x,c1y,c1z,c2x,c2y,c2z,c3x,c3y,xcnt);
    calculate(); ycnt = cell.sgmmat[2][2]/MPa;
    if (ycnt > ytarget) { xmin = xcnt; }
    if (ycnt < ytarget) { xmax = xcnt; }
  }
  //printf("\n");
  if (min_or_max==0) {
  stretch_celladjust_d(c1x,c1y,c1z,c2x,c2y,c2z,c3x,c3y,xmin);
  calculate(); ymin = cell.sgmmat[2][2]/MPa; double zx1 = cell.sgmmat[2][0]/MPa;
  ep = atom.epotsum/eV/(double)atom.natom;
  //printf("%12.8f %12.8f %20.10e %12.8f %12.8f\n",c3x,xmin,ep,ymin,zx1);
  printf("%12.8f %12.8f %12.8f %20.10e %12.8f %12.8f\n",sftx/ang,sfty/ang,xmin,ep,ymin,zx1);
  } else {
  stretch_celladjust_d(c1x,c1y,c1z,c2x,c2y,c2z,c3x,c3y,xmax);
  calculate(); ymax = cell.sgmmat[2][2]/MPa; double zx2 = cell.sgmmat[2][0]/MPa;
  ep = atom.epotsum/eV/(double)atom.natom;
  //printf("%12.8f %12.8f %20.10e %12.8f %12.8f\n",c3x,xmax,ep,ymax,zx2);
  printf("%12.8f %12.8f %12.8f %20.10e %12.8f %12.8f\n",sftx/ang,sfty/ang,xmax,ep,ymax,zx2);
  }

  }printf("\n");}  

  //for (double x=(double)auto_val2-0.0001; x<(double)auto_val2+0.0001; x += 0.00001) {
  //  c3z = x;
  //  stretch_celladjust_d(c1x,c1y,c1z,c2x,c2y,c2z,c3x,c3y,c3z);
  //  calculate();
  //  printf("%f %f %f %f\n",c3z,c3x,cell.sgmmat[2][2]/MPa,cell.sgmmat[2][0]/MPa);
  //}
}

void autocalc_old()
{
  printf("Auto calculation\n");
  double c1x, c1y, c1z, c2x, c2y, c2z, c3x, c3y, c3z;
  c1x = cell.hmat[0][0]/ang;
  c1y = cell.hmat[1][0]/ang;
  c1z = cell.hmat[2][0]/ang;
  c2x = cell.hmat[0][1]/ang;
  c2y = cell.hmat[1][1]/ang;
  c2z = cell.hmat[2][1]/ang;
  //c3x = cell.hmat[0][2]/ang;
  c3x = (double)auto_val1;
  c3y = cell.hmat[1][2]/ang;
  c3z = cell.hmat[2][2]/ang;

  double ytarget = -1000;
  double tol = 10.0;
  double xmin = (double)auto_val2 - 0.005;
  double xmax = (double)auto_val2 + 0.005;
  double xcnt;
  double ymin, ymax, ycnt;
  double ep;
  stretch_celladjust_d(c1x,c1y,c1z,c2x,c2y,c2z,c3x,c3y,xmin); 
  calculate(); ymin = cell.sgmmat[2][2]/MPa;
  stretch_celladjust_d(c1x,c1y,c1z,c2x,c2y,c2z,c3x,c3y,xmax);
  calculate(); ymax = cell.sgmmat[2][2]/MPa;
  printf("Initial: x = %f - %f, y = %f - %f\n",xmin,xmax,ymin,ymax);
START:
  int min_or_max;
  for (int ite=0; ite<20; ite++) {
    xcnt = (xmin+xmax)/2.0;
    stretch_celladjust_d(c1x,c1y,c1z,c2x,c2y,c2z,c3x,c3y,xmin);
    calculate(); ymin = cell.sgmmat[2][2]/MPa;
    if (abs(ymin-ytarget) < tol) { min_or_max=0; break; }
    if (ymin > ytarget) {
      double xx = xmax-xmin;
      xmin -= xx; xmax -= xx; ite--; continue; }
    stretch_celladjust_d(c1x,c1y,c1z,c2x,c2y,c2z,c3x,c3y,xmax);
    calculate(); ymax = cell.sgmmat[2][2]/MPa;
    if (abs(ymax-ytarget) < tol) { min_or_max=1; break; }
    if (ymax < ytarget) {
      double xx = xmax-xmin;
      xmax += xx; xmin += xx; ite--; continue; }
    stretch_celladjust_d(c1x,c1y,c1z,c2x,c2y,c2z,c3x,c3y,xcnt);
    calculate(); ycnt = cell.sgmmat[2][2]/MPa;
    if (ycnt > ytarget) { xmin = xcnt; }
    if (ycnt < ytarget) { xmax = xcnt; }
  }
  printf("\n");
  if (min_or_max==0) {
  stretch_celladjust_d(c1x,c1y,c1z,c2x,c2y,c2z,c3x,c3y,xmin);
  calculate(); ymin = cell.sgmmat[2][2]/MPa; double zx1 = cell.sgmmat[2][0]/MPa;
  ep = atom.epotsum/eV/(double)atom.natom;
  printf("%12.8f %12.8f %20.10e %12.8f %12.8f\n",c3x,xmin,ep,ymin,zx1);
  } else {
  stretch_celladjust_d(c1x,c1y,c1z,c2x,c2y,c2z,c3x,c3y,xmax);
  calculate(); ymax = cell.sgmmat[2][2]/MPa; double zx2 = cell.sgmmat[2][0]/MPa;
  ep = atom.epotsum/eV/(double)atom.natom;
  printf("%12.8f %12.8f %20.10e %12.8f %12.8f\n",c3x,xmax,ep,ymax,zx2);
  }

  //for (double x=(double)auto_val2-0.0001; x<(double)auto_val2+0.0001; x += 0.00001) {
  //  c3z = x;
  //  stretch_celladjust_d(c1x,c1y,c1z,c2x,c2y,c2z,c3x,c3y,c3z);
  //  calculate();
  //  printf("%f %f %f %f\n",c3z,c3x,cell.sgmmat[2][2]/MPa,cell.sgmmat[2][0]/MPa);
  //}
}

void calculate() {
    bookkeep(); potential(); loading();
    // For output in "Status" area
    //cellx=cell.hmat[0][0]/ang;celly=cell.hmat[1][1]/ang;cellz=cell.hmat[2][2]/ang;
    cellx=sqrt(cell.hmat[0][0]*cell.hmat[0][0]+cell.hmat[1][0]*cell.hmat[1][0]+cell.hmat[2][0]*cell.hmat[2][0])/ang;
    celly=sqrt(cell.hmat[0][1]*cell.hmat[0][1]+cell.hmat[1][1]*cell.hmat[1][1]+cell.hmat[2][1]*cell.hmat[2][1])/ang;
    cellz=sqrt(cell.hmat[0][2]*cell.hmat[0][2]+cell.hmat[1][2]*cell.hmat[1][2]+cell.hmat[2][2]*cell.hmat[2][2])/ang;
    cell1x=cell.hmat[0][0]/ang;cell1y=cell.hmat[1][0]/ang;cell1z=cell.hmat[2][0]/ang;
    cell2x=cell.hmat[0][1]/ang;cell2y=cell.hmat[1][1]/ang;cell2z=cell.hmat[2][1]/ang;
    cell3x=cell.hmat[0][2]/ang;cell3y=cell.hmat[1][2]/ang;cell3z=cell.hmat[2][2]/ang;
    f_max=atom.Fmax()/eV*ang;
    epotatom=atom.epotsum/eV/atom.natom;

}

//void autocalc_dwcnt_rotate()
void autocalc()
{
  double outermost_radius, cylinder_side_area;
  double center_x, center_y;
  calc_center_xy(center_x, center_y);
  calc_cylinder_side_area(center_x, center_y, outermost_radius, cylinder_side_area);
  double xx, yy, rr, t, tt;
  
  //int angle = 1;
  for (float angle=0; angle<=18; angle += 0.2) {
    for (int i=1; i<=atom.natom; i++) {
      xx = atom.rx[i] - center_x;
      yy = atom.ry[i] - center_y;
      rr = sqrt(xx*xx+yy*yy);
      if (fabs(rr-outermost_radius)/ang < 0.1) {
	t  = theta(xx,yy);
	tt = t + (double)angle /180.0 * M_PI;
	atom.rx[i] = center_x + rr * cos(tt);
	atom.ry[i] = center_y + rr * sin(tt);
      }
    }
    bookkeep();potential();
    printf("%f %20.10e\n",(float)angle,atom.epotsum/eV/atom.natom);
  }
}
