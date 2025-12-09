// To solve equation of motion
#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "myheader.h"

#define BL 1.380662e-23
void inverse(double mat[3][3], double imat[3][3]);
void mk_helmat(int iel, double rx[], double ry[], double rz[], double mat[][3]);
double q_x(int i, int iel, double rx[], double ry[], double rz[], double hinelmat[][3]);
double q_y(int i, int iel, double rx[], double ry[], double rz[], double hinelmat[][3]);
double q_z(int i, int iel, double rx[], double ry[], double rz[], double hinelmat[][3]);
double r_x(int i, int iel, double rx[], double ry[], double rz[], double helmat[][3], 
	   double qx, double qy, double qz);
double r_y(int i, int iel, double rx[], double ry[], double rz[], double helmat[][3], 
	   double qx, double qy, double qz);
double r_z(int i, int iel, double rx[], double ry[], double rz[], double helmat[][3], 
	   double qx, double qy, double qz);
void recipro(double h[3][3], double r[3][3]);
void transpose(double a[3][3], double b[3][3]);
void matmul(double a[3][3], double b[3][3], double c[3][3]);
void matadd(double a[3][3], double b[3][3], double c[3][3]);
void matsub(double a[3][3], double b[3][3], double c[3][3]);
void resetmat(double a[3][3]);
void matcpy(double a[3][3], double b[3][3]);

extern bool iprdamper;
extern float prdamper_val1, prdamper_val2, prlimit;

void velocity_verlet_a()
{
  double helmat_p[3][3], hinelmat_p[3][3], helmat[3][3];
  double qx, qy, qz, xx, yy, zz;
  for (int i=1; i<=atom.natom; i++)
    {
      atom.rx_p[i] = atom.rx[i]; atom.ry_p[i] = atom.ry[i]; atom.rz_p[i] = atom.rz[i]; 
      if ((atom.QC==0)||(atom.repatom[i] == 1)) {
	atom.rx[i] = atom.rx[i] + dt * atom.vx[i] + (dt*dt/2.0) * atom.fx[i] / atom.wm[i];
	atom.ry[i] = atom.ry[i] + dt * atom.vy[i] + (dt*dt/2.0) * atom.fy[i] / atom.wm[i];
	atom.rz[i] = atom.rz[i] + dt * atom.vz[i] + (dt*dt/2.0) * atom.fz[i] / atom.wm[i];
	atom.vx[i] = atom.vx[i] + dt/2.0 * atom.fx[i] / atom.wm[i];
	atom.vy[i] = atom.vy[i] + dt/2.0 * atom.fy[i] / atom.wm[i];
	atom.vz[i] = atom.vz[i] + dt/2.0 * atom.fz[i] / atom.wm[i];
	if (atom.mfx[i] == true) { atom.rx[i] = atom.rx_p[i]; atom.vx[i] = 0.0; }
	if (atom.mfy[i] == true) { atom.ry[i] = atom.ry_p[i]; atom.vy[i] = 0.0; }
	if (atom.mfz[i] == true) { atom.rz[i] = atom.rz_p[i]; atom.vz[i] = 0.0; }
  if ((ifgrab1 == 1)&&(atom.group[i] == grab1num)) {
    atom.rx[i] = atom.rx_p[i] + grabdxdt1 *(1.0e-10/1.0e-12)* dt;
    atom.ry[i] = atom.ry_p[i] + grabdydt1 *(1.0e-10/1.0e-12)* dt;
    atom.rz[i] = atom.rz_p[i] + grabdzdt1 *(1.0e-10/1.0e-12)* dt;
  }
  if ((ifgrab2 == 1)&&(atom.group[i] == grab2num)) {
    atom.rx[i] = atom.rx_p[i] + grabdxdt2 *(1.0e10/1.0e-12)* dt;
    atom.ry[i] = atom.ry_p[i] + grabdydt2 *(1.0e10/1.0e-12)* dt;
    atom.rz[i] = atom.rz_p[i] + grabdzdt2 *(1.0e10/1.0e-12)* dt;
  }
      }
    }
  if (atom.QC) {
  for (int iel=1; iel<=atom.nelem; iel++)
    {
      mk_helmat(iel, atom.rx_p, atom.ry_p, atom.rz_p, helmat_p);
      mk_helmat(iel, atom.rx, atom.ry, atom.rz, helmat);
      inverse(helmat_p, hinelmat_p);
      for (int i=1; i<=atom.natom; i++)
	{
	  if ((atom.repatom[i] == 0)&&(atom.elem_id[i] == iel)) {
	    //	    mk_qxy(i, iel, atom.rx_p, atom.ry_p, hinelmat_p, qx, qy);
	    qx = q_x(i, iel, atom.rx_p, atom.ry_p, atom.rz_p, hinelmat_p);
	    qy = q_y(i, iel, atom.rx_p, atom.ry_p, atom.rz_p, hinelmat_p);
	    qz = q_z(i, iel, atom.rx_p, atom.ry_p, atom.rz_p, hinelmat_p);
	    //	    mk_rxy(i, iel, atom.rx, atom.ry, helmat, qx, qy);
	    atom.rx[i] = r_x(i, iel, atom.rx, atom.ry, atom.rz, helmat, qx, qy, qz);
	    atom.ry[i] = r_y(i, iel, atom.rx, atom.ry, atom.rz, helmat, qx, qy, qz);
	    atom.rz[i] = r_z(i, iel, atom.rx, atom.ry, atom.rz, helmat, qx, qy, qz);
	    atom.vx[i] = ( atom.rx[i]-atom.rx_p[i] ) / dt;
	    atom.vy[i] = ( atom.ry[i]-atom.ry_p[i] ) / dt;
	    atom.vz[i] = ( atom.rz[i]-atom.rz_p[i] ) / dt;
	  }
	}
    }
  }
}

void velocity_verlet_b()
{
  for (int i=1; i<=atom.natom; i++)
    {
      if ((atom.QC==0)||(atom.repatom[i] == 1)) {
	atom.vx[i] = atom.vx[i] + dt/2.0 * atom.fx[i] / atom.wm[i];
	atom.vy[i] = atom.vy[i] + dt/2.0 * atom.fy[i] / atom.wm[i];
	atom.vz[i] = atom.vz[i] + dt/2.0 * atom.fz[i] / atom.wm[i];
      }
    }
}

void gear_pc_predic()
{
  double helmat_p[3][3], hinelmat_p[3][3], helmat[3][3];
  double qx, qy, qz, xx, yy, zz;
  double c1, c2, c3, c4;
  c1 = dt; c2 = c1 * dt / 2.0; c3 = c2 * dt / 3.0; c4 = c3 * dt / 4.0;
  for (int i=1; i<=atom.natom; i++)
    {
      atom.rx_p[i] = atom.rx[i]; atom.ry_p[i] = atom.ry[i]; atom.rz_p[i] = atom.rz[i]; 
      if ((atom.QC==0)||(atom.repatom[i] == 1)) {
	atom.rx[i] += c1 * atom.vx[i] + c2 * atom.ax[i] + c3 * atom.bx[i] + c4 * atom.cx[i];
	atom.ry[i] += c1 * atom.vy[i] + c2 * atom.ay[i] + c3 * atom.by[i] + c4 * atom.cy[i];
	atom.rz[i] += c1 * atom.vz[i] + c2 * atom.az[i] + c3 * atom.bz[i] + c4 * atom.cz[i];
	atom.vx[i] += c1 * atom.ax[i] + c2 * atom.bx[i] + c3 * atom.cx[i];
	atom.vy[i] += c1 * atom.ay[i] + c2 * atom.by[i] + c3 * atom.cy[i];
	atom.vz[i] += c1 * atom.az[i] + c2 * atom.bz[i] + c3 * atom.cz[i];
	atom.ax[i] += c1 * atom.bx[i] + c2 * atom.cx[i];
	atom.ay[i] += c1 * atom.by[i] + c2 * atom.cy[i];
	atom.az[i] += c1 * atom.bz[i] + c2 * atom.cz[i];
	atom.bx[i] += c1 * atom.cx[i];
	atom.by[i] += c1 * atom.cy[i];
	atom.bz[i] += c1 * atom.cz[i];

	if (atom.mfx[i] == true) { atom.rx[i] = atom.rx_p[i]; atom.vx[i] = 0.0; }
	if (atom.mfy[i] == true) { atom.ry[i] = atom.ry_p[i]; atom.vy[i] = 0.0; }
	if (atom.mfz[i] == true) { atom.rz[i] = atom.rz_p[i]; atom.vz[i] = 0.0; }
  if ((ifgrab1 == 1)&&(i == grab1num)) {
    atom.rx[i] = atom.rx_p[i] + grabdxdt1 * dt;
    atom.ry[i] = atom.ry_p[i] + grabdydt1 * dt;
    atom.rz[i] = atom.rz_p[i] + grabdzdt1 * dt;
  }
  if ((ifgrab1 == 2)&&(i == grab2num)) {
    atom.rx[i] = atom.rx_p[i] + grabdxdt2 * dt;
    atom.ry[i] = atom.ry_p[i] + grabdydt2 * dt;
    atom.rz[i] = atom.rz_p[i] + grabdzdt2 * dt;
  }
      }
    }
  if (atom.QC) {
  for (int iel=1; iel<=atom.nelem; iel++)
    {
      mk_helmat(iel, atom.rx_p, atom.ry_p, atom.rz_p, helmat_p);
      mk_helmat(iel, atom.rx, atom.ry, atom.rz, helmat);
      inverse(helmat_p, hinelmat_p);
      for (int i=1; i<=atom.natom; i++)
	{
	  if ((atom.repatom[i] == 0)&&(atom.elem_id[i] == iel)) {
	    //	    mk_qxy(i, iel, atom.rx_p, atom.ry_p, hinelmat_p, qx, qy);
	    qx = q_x(i, iel, atom.rx_p, atom.ry_p, atom.rz_p, hinelmat_p);
	    qy = q_y(i, iel, atom.rx_p, atom.ry_p, atom.rz_p, hinelmat_p);
	    qz = q_z(i, iel, atom.rx_p, atom.ry_p, atom.rz_p, hinelmat_p);
	    //	    mk_rxy(i, iel, atom.rx, atom.ry, helmat, qx, qy);
	    atom.rx[i] = r_x(i, iel, atom.rx, atom.ry, atom.rz, helmat, qx, qy, qz);
	    atom.ry[i] = r_y(i, iel, atom.rx, atom.ry, atom.rz, helmat, qx, qy, qz);
	    atom.rz[i] = r_z(i, iel, atom.rx, atom.ry, atom.rz, helmat, qx, qy, qz);
	    atom.vx[i] = ( atom.rx[i]-atom.rx_p[i] ) / dt;
	    atom.vy[i] = ( atom.ry[i]-atom.ry_p[i] ) / dt;
	    atom.vz[i] = ( atom.rz[i]-atom.rz_p[i] ) / dt;
	  }
	}
    }
  }
}

void gear_pc_correc()
{
  double helmat_p[3][3], hinelmat_p[3][3], helmat[3][3];
  double qx, qy, qz;
  double gear0, gear1, gear2, gear3, gear4;
  double c1, c2, c3, c4, cr, cv, cb, cc;
  double axi, ayi, azi, corrx, corry, corrz;
  gear0 = 19.0/120.0, gear1 = 3.0/4.0, gear3 = 1.0/2.0, gear4 = 1.0/12.0;
  c1 = dt; c2 = c1 * dt / 2.0; c3 = c2 * dt / 3.0; c4 = c3 * dt / 4.0;
  cr = gear0 * c2; cv = gear1 * c2 / c1; cb = gear3 * c2 / c3; cc = gear4 * c2 / c4;
  for (int i=1; i<=atom.natom; i++)
    {
      if ((atom.QC==0)||(atom.repatom[i] == 1)) {
	axi = atom.fx[i] / atom.wm[i];
	ayi = atom.fy[i] / atom.wm[i];
	azi = atom.fz[i] / atom.wm[i];
	corrx = axi - atom.ax[i];
	corry = ayi - atom.ay[i];
	corrz = azi - atom.az[i];
	atom.rx[i] += cr * corrx;
	atom.ry[i] += cr * corry;
	atom.rz[i] += cr * corrz;
	atom.vx[i] += cv * corrx;
	atom.vy[i] += cv * corry;
	atom.vz[i] += cv * corrz;
	atom.ax[i] = axi;
	atom.ay[i] = ayi;
	atom.az[i] = azi;
	atom.bx[i] += cb * corrx;
	atom.by[i] += cb * corry;
	atom.bz[i] += cb * corrz;
	atom.cx[i] += cc * corrx;
	atom.cy[i] += cc * corry;
	atom.cz[i] += cc * corrz;
      }
    }
  if (atom.QC) {
  for (int iel=1; iel<=atom.nelem; iel++)
    {
      mk_helmat(iel, atom.rx_p, atom.ry_p, atom.rz_p, helmat_p);
      mk_helmat(iel, atom.rx, atom.ry, atom.rz, helmat);
      inverse(helmat_p, hinelmat_p);
      for (int i=1; i<=atom.natom; i++)
	{
	  if ((atom.repatom[i] == 0)&&(atom.elem_id[i] == iel)) {
	    //	    mk_qxy(i, iel, atom.rx_p, atom.ry_p, hinelmat_p, qx, qy);
	    qx = q_x(i, iel, atom.rx_p, atom.ry_p, atom.rz_p, hinelmat_p);
	    qy = q_y(i, iel, atom.rx_p, atom.ry_p, atom.rz_p, hinelmat_p);
	    qz = q_z(i, iel, atom.rx_p, atom.ry_p, atom.rz_p, hinelmat_p);
	    //	    mk_rxy(i, iel, atom.rx, atom.ry, helmat, qx, qy);
	    atom.rx[i] = r_x(i, iel, atom.rx, atom.ry, atom.rz, helmat, qx, qy, qz);
	    atom.ry[i] = r_y(i, iel, atom.rx, atom.ry, atom.rz, helmat, qx, qy, qz);
	    atom.rz[i] = r_z(i, iel, atom.rx, atom.ry, atom.rz, helmat, qx, qy, qz);
	    atom.vx[i] = ( atom.rx[i]-atom.rx_p[i] ) / dt;
	    atom.vy[i] = ( atom.ry[i]-atom.ry_p[i] ) / dt;
	    atom.vz[i] = ( atom.rz[i]-atom.rz_p[i] ) / dt;
	  }
	}
    }
  }
}

void gear_pc_nphpr_predic()
{
  double helmat_p[3][3], hinelmat_p[3][3], helmat[3][3];
  double qx, qy, qz, xx, yy, zz;
  double c1, c2, c3, c4;
  c1 = dt; c2 = c1 * dt / 2.0; c3 = c2 * dt / 3.0; c4 = c3 * dt / 4.0;
  for (int i=1; i<=atom.natom; i++)
    {
      atom.rx_p[i] = atom.rx[i]; atom.ry_p[i] = atom.ry[i]; atom.rz_p[i] = atom.rz[i]; 
      if ((atom.QC==0)||(atom.repatom[i] == 1)) {
	atom.rx[i] += c1 * atom.vx[i] + c2 * atom.ax[i] + c3 * atom.bx[i] + c4 * atom.cx[i];
	atom.ry[i] += c1 * atom.vy[i] + c2 * atom.ay[i] + c3 * atom.by[i] + c4 * atom.cy[i];
	atom.rz[i] += c1 * atom.vz[i] + c2 * atom.az[i] + c3 * atom.bz[i] + c4 * atom.cz[i];
	atom.vx[i] += c1 * atom.ax[i] + c2 * atom.bx[i] + c3 * atom.cx[i];
	atom.vy[i] += c1 * atom.ay[i] + c2 * atom.by[i] + c3 * atom.cy[i];
	atom.vz[i] += c1 * atom.az[i] + c2 * atom.bz[i] + c3 * atom.cz[i];
	atom.ax[i] += c1 * atom.bx[i] + c2 * atom.cx[i];
	atom.ay[i] += c1 * atom.by[i] + c2 * atom.cy[i];
	atom.az[i] += c1 * atom.bz[i] + c2 * atom.cz[i];
	atom.bx[i] += c1 * atom.cx[i];
	atom.by[i] += c1 * atom.cy[i];
	atom.bz[i] += c1 * atom.cz[i];

	if (atom.mfx[i] == true) { atom.rx[i] = atom.rx_p[i]; atom.vx[i] = 0.0; }
	if (atom.mfy[i] == true) { atom.ry[i] = atom.ry_p[i]; atom.vy[i] = 0.0; }
	if (atom.mfz[i] == true) { atom.rz[i] = atom.rz_p[i]; atom.vz[i] = 0.0; }
  if ((ifgrab1 == 1)&&(i == grab1num)) {
    atom.rx[i] = atom.rx_p[i] + grabdxdt1 * dt;
    atom.ry[i] = atom.ry_p[i] + grabdydt1 * dt;
    atom.rz[i] = atom.rz_p[i] + grabdzdt1 * dt;
  }
  if ((ifgrab1 == 2)&&(i == grab2num)) {
    atom.rx[i] = atom.rx_p[i] + grabdxdt2 * dt;
    atom.ry[i] = atom.ry_p[i] + grabdydt2 * dt;
    atom.rz[i] = atom.rz_p[i] + grabdzdt2 * dt;
  }
      }
    }
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      if (!cell.fix[i][j]) { //Jan2015
	cell.hmat[i][j]  += c1 * cell.hvmat[i][j] + c2 * cell.hamat[i][j] + c3 * cell.hbmat[i][j]
	  + c4 * cell.hcmat[i][j];
	cell.hvmat[i][j] += c1 * cell.hamat[i][j] + c2 * cell.hbmat[i][j] + c3 * cell.hcmat[i][j];
	cell.hamat[i][j] += c1 * cell.hbmat[i][j] + c2 * cell.hcmat[i][j];
	if (iprdamper) {
	  cell.hamat[i][j] += cell.hvmat[i][j] * cell.damper_param * (double)prdamper_val1;
	}
	cell.hbmat[i][j] += c1 * cell.hcmat[i][j];
      } else {
	cell.sgmmat_set[i][j] = cell.sgmmat[i][j]; //Jan2015
      }
    }
  }
  if (atom.QC) {
  for (int iel=1; iel<=atom.nelem; iel++)
    {
      mk_helmat(iel, atom.rx_p, atom.ry_p, atom.rz_p, helmat_p);
      mk_helmat(iel, atom.rx, atom.ry, atom.rz, helmat);
      inverse(helmat_p, hinelmat_p);
      for (int i=1; i<=atom.natom; i++)
	{
	  if ((atom.repatom[i] == 0)&&(atom.elem_id[i] == iel)) {
	    //	    mk_qxy(i, iel, atom.rx_p, atom.ry_p, hinelmat_p, qx, qy);
	    qx = q_x(i, iel, atom.rx_p, atom.ry_p, atom.rz_p, hinelmat_p);
	    qy = q_y(i, iel, atom.rx_p, atom.ry_p, atom.rz_p, hinelmat_p);
	    qz = q_z(i, iel, atom.rx_p, atom.ry_p, atom.rz_p, hinelmat_p);
	    //	    mk_rxy(i, iel, atom.rx, atom.ry, helmat, qx, qy);
	    atom.rx[i] = r_x(i, iel, atom.rx, atom.ry, atom.rz, helmat, qx, qy, qz);
	    atom.ry[i] = r_y(i, iel, atom.rx, atom.ry, atom.rz, helmat, qx, qy, qz);
	    atom.rz[i] = r_z(i, iel, atom.rx, atom.ry, atom.rz, helmat, qx, qy, qz);
	    atom.vx[i] = ( atom.rx[i]-atom.rx_p[i] ) / dt;
	    atom.vy[i] = ( atom.ry[i]-atom.ry_p[i] ) / dt;
	    atom.vz[i] = ( atom.rz[i]-atom.rz_p[i] ) / dt;
	  }
	}
    }
  }
}

void gear_pc_nphpr_correc()
{
  double helmat_p[3][3], hinelmat_p[3][3], helmat[3][3];
  double qx, qy, qz, xx, yy, zz;
  double gear0, gear1, gear2, gear3, gear4;
  double c1, c2, c3, c4, cr, cv, cb, cc;
  double axi, ayi, azi, corrx, corry, corrz;
  double gmat[3][3], gvmat[3][3], ginmat[3][3];
  double pt[3][3], tm[3][3], tm2[3][3], tm3[3][3];
  double hf[3][3], rec[3][3], corrm[3][3], tmp[3][3];
  double hmat_t[3][3], hvmat_t[3][3], hinmat_t[3][3]; // transpose
  gear0 = 19.0/90.0, gear1 = 3.0/4.0, gear3 = 1.0/2.0, gear4 = 1.0/12.0;
  c1 = dt; c2 = c1 * dt / 2.0; c3 = c2 * dt / 3.0; c4 = c3 * dt / 4.0;
  cr = gear0 * c2; cv = gear1 * c2 / c1; cb = gear3 * c2 / c3; cc = gear4 * c2 / c4;
  cell.ww = cell.prmass * pow(10,prmass_scale);

  inverse(cell.hmat,cell.hinmat);
  cell.volume = cell.Getvolume();
  recipro(cell.hmat,rec);
  resetmat(pt);
  matmul(cell.hvmat,cell.hinmat,tm);

  double enkin = 0;
  for (int i=1; i<=atom.natom; i++) {
    if ((atom.QC==0)||(atom.repatom[i] == 1)) {
      xx = atom.vx[i] - (tm[0][0]*atom.rx[i]+tm[0][1]*atom.ry[i]+tm[0][2]*atom.rz[i]);
      yy = atom.vy[i] - (tm[1][0]*atom.rx[i]+tm[1][1]*atom.ry[i]+tm[1][2]*atom.rz[i]);
      zz = atom.vz[i] - (tm[2][0]*atom.rx[i]+tm[2][1]*atom.ry[i]+tm[2][2]*atom.rz[i]);
      enkin += 0.5 * atom.wm[i]*(xx*xx+yy*yy+zz*zz); //Enkin must be calculated here
      pt[0][0] += xx * xx * atom.wm[i];
      pt[1][0] += yy * xx * atom.wm[i];
      pt[2][0] += zz * xx * atom.wm[i];
      pt[0][1] += xx * yy * atom.wm[i];
      pt[1][1] += yy * yy * atom.wm[i];
      pt[2][1] += zz * yy * atom.wm[i];
      pt[0][2] += xx * zz * atom.wm[i];
      pt[1][2] += yy * zz * atom.wm[i];
      pt[2][2] += zz * zz * atom.wm[i];
    }
  }
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      pt[i][j] += cell.dmat[i][j];
      pt[i][j] /= cell.volume;
      pt[i][j] += cell.sgmmat_set[i][j];
    }
  }
  matmul(pt,rec,hf);
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      hf[i][j] /= cell.ww;
    }
  }
  transpose(cell.hmat,   hmat_t);
  transpose(cell.hvmat,  hvmat_t);
  transpose(cell.hinmat, hinmat_t);
  matmul(hmat_t,  cell.hmat,  gmat);
  matmul(hvmat_t, cell.hmat,  tmp);
  matmul(hmat_t,  cell.hvmat, gvmat);
  for (int i=0; i<3; i++) { for (int j=0; j<3; j++) {
      gvmat[i][j] += tmp[i][j];  }  }
  matmul(gvmat, cell.hinmat, tm2);
  matcpy(tm2,tmp); matmul(hinmat_t, tmp, tm2);
  matcpy(tm2,tmp); matmul(cell.hinmat, tmp, tm2);
  matcpy(tm2,tmp); matmul(cell.hmat, tmp, tm2);
  matmul(cell.hamat, cell.hinmat, tm3);

  for (int i=1; i<=atom.natom; i++) {
    if ((atom.QC==0)||(atom.repatom[i] == 1)) {
      xx = atom.vx[i] - (tm[0][0]*atom.rx[i]+tm[0][1]*atom.ry[i]+tm[0][2]*atom.rz[i]);
      yy = atom.vy[i] - (tm[1][0]*atom.rx[i]+tm[1][1]*atom.ry[i]+tm[1][2]*atom.rz[i]);
      zz = atom.vz[i] - (tm[2][0]*atom.rx[i]+tm[2][1]*atom.ry[i]+tm[2][2]*atom.rz[i]);
      axi = atom.fx[i] / atom.wm[i] - (tm2[0][0]*xx+tm2[0][1]*yy+tm2[0][2]*zz)
	+ 2*(tm[0][0]*xx+tm[0][1]*yy+tm[0][2]+zz)
	+ (tm3[0][0]*atom.rx[i]+tm3[0][1]*atom.ry[i]+tm3[0][2]*atom.rz[i]);
      ayi = atom.fy[i] / atom.wm[i] - (tm2[1][0]*xx+tm2[1][1]*yy+tm2[1][2]*zz)
	+ 2*(tm[1][0]*xx+tm[1][1]*yy+tm[1][2]+zz)
	+ (tm3[1][0]*atom.rx[i]+tm3[1][1]*atom.ry[i]+tm3[1][2]*atom.rz[i]);
      azi = atom.fz[i] / atom.wm[i] - (tm2[2][0]*xx+tm2[2][1]*yy+tm2[2][2]*zz)
	+ 2*(tm[2][0]*xx+tm[2][1]*yy+tm[2][2]+zz)
	+ (tm3[2][0]*atom.rx[i]+tm3[2][1]*atom.ry[i]+tm3[2][2]*atom.rz[i]);
      corrx = axi - atom.ax[i];
      corry = ayi - atom.ay[i];
      corrz = azi - atom.az[i];
      atom.rx[i] += cr * corrx;
      atom.ry[i] += cr * corry;
      atom.rz[i] += cr * corrz;
      atom.vx[i] += cv * corrx;
      atom.vy[i] += cv * corry;
      atom.vz[i] += cv * corrz;
      atom.ax[i] = axi;
      atom.ay[i] = ayi;
      atom.az[i] = azi;
      atom.bx[i] += cb * corrx;
      atom.by[i] += cb * corry;
      atom.bz[i] += cb * corrz;
      atom.cx[i] += cc * corrx;
      atom.cy[i] += cc * corry;
      atom.cz[i] += cc * corrz;
    }
  }
  matsub(hf, cell.hamat, corrm);
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      if (!cell.fix[i][j]) { //Jan2015
	//if (i==j) printf("%d %e %e %e  \n",i,pt[i][j],rec[i][j],hf[i][j]);
	cell.hmat[i][j]  += cr * corrm[i][j];
	cell.hvmat[i][j] += cv * corrm[i][j];
	cell.hamat[i][j]  = hf[i][j];
	cell.hbmat[i][j] += cb * corrm[i][j];
	cell.hcmat[i][j] += cc * corrm[i][j];
      }
    }
  }
  /*
  if ((iprdamper)&&(prdamper_val2>0)) {
    for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
	cell.hvmat[i][j]=cell.hvmat[i][j]/(1.0+(double)prdamper_val2);
	cell.hamat[i][j]=0;
	cell.hbmat[i][j]=0;
	cell.hcmat[i][j]=0;
      }
    }
  }
  */
  if ((iprdamper)&&(prlimit>0)) {
    for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
	if (abs(cell.hvmat[i][j])>prlimit) {
	  cell.hvmat[i][j] /= 2.0;
	  cell.hamat[i][j]=0;
	  cell.hbmat[i][j]=0;
	  cell.hcmat[i][j]=0;
	}
      }
    }
  }

  cell.volume = cell.Getvolume();
  if (atom.QC) {
  for (int iel=1; iel<=atom.nelem; iel++)
    {
      mk_helmat(iel, atom.rx_p, atom.ry_p, atom.rz_p, helmat_p);
      mk_helmat(iel, atom.rx, atom.ry, atom.rz, helmat);
      inverse(helmat_p, hinelmat_p);
      for (int i=1; i<=atom.natom; i++)
	{
	  if ((atom.repatom[i] == 0)&&(atom.elem_id[i] == iel)) {
	    //	    mk_qxy(i, iel, atom.rx_p, atom.ry_p, hinelmat_p, qx, qy);
	    qx = q_x(i, iel, atom.rx_p, atom.ry_p, atom.rz_p, hinelmat_p);
	    qy = q_y(i, iel, atom.rx_p, atom.ry_p, atom.rz_p, hinelmat_p);
	    qz = q_z(i, iel, atom.rx_p, atom.ry_p, atom.rz_p, hinelmat_p);
	    //	    mk_rxy(i, iel, atom.rx, atom.ry, helmat, qx, qy);
	    atom.rx[i] = r_x(i, iel, atom.rx, atom.ry, atom.rz, helmat, qx, qy, qz);
	    atom.ry[i] = r_y(i, iel, atom.rx, atom.ry, atom.rz, helmat, qx, qy, qz);
	    atom.rz[i] = r_z(i, iel, atom.rx, atom.ry, atom.rz, helmat, qx, qy, qz);
	    atom.vx[i] = ( atom.rx[i]-atom.rx_p[i] ) / dt;
	    atom.vy[i] = ( atom.ry[i]-atom.ry_p[i] ) / dt;
	    atom.vz[i] = ( atom.rz[i]-atom.rz_p[i] ) / dt;
	  }
	}
    }
  }
}

void gear_pc_nvtpr_predic()
{
  double helmat_p[3][3], hinelmat_p[3][3], helmat[3][3];
  double qx, qy, qz, xx, yy, zz;
  double c1, c2, c3, c4;
  c1 = dt; c2 = c1 * dt / 2.0; c3 = c2 * dt / 3.0; c4 = c3 * dt / 4.0;
  for (int i=1; i<=atom.natom; i++)
    {
      atom.rx_p[i] = atom.rx[i]; atom.ry_p[i] = atom.ry[i]; atom.rz_p[i] = atom.rz[i]; 
      if ((atom.QC==0)||(atom.repatom[i] == 1)) {
	atom.rx[i] += c1 * atom.vx[i] + c2 * atom.ax[i] + c3 * atom.bx[i] + c4 * atom.cx[i];
	atom.ry[i] += c1 * atom.vy[i] + c2 * atom.ay[i] + c3 * atom.by[i] + c4 * atom.cy[i];
	atom.rz[i] += c1 * atom.vz[i] + c2 * atom.az[i] + c3 * atom.bz[i] + c4 * atom.cz[i];
	atom.vx[i] += c1 * atom.ax[i] + c2 * atom.bx[i] + c3 * atom.cx[i];
	atom.vy[i] += c1 * atom.ay[i] + c2 * atom.by[i] + c3 * atom.cy[i];
	atom.vz[i] += c1 * atom.az[i] + c2 * atom.bz[i] + c3 * atom.cz[i];
	atom.ax[i] += c1 * atom.bx[i] + c2 * atom.cx[i];
	atom.ay[i] += c1 * atom.by[i] + c2 * atom.cy[i];
	atom.az[i] += c1 * atom.bz[i] + c2 * atom.cz[i];
	atom.bx[i] += c1 * atom.cx[i];
	atom.by[i] += c1 * atom.cy[i];
	atom.bz[i] += c1 * atom.cz[i];

	if (atom.mfx[i] == true) { atom.rx[i] = atom.rx_p[i]; atom.vx[i] = 0.0; }
	if (atom.mfy[i] == true) { atom.ry[i] = atom.ry_p[i]; atom.vy[i] = 0.0; }
	if (atom.mfz[i] == true) { atom.rz[i] = atom.rz_p[i]; atom.vz[i] = 0.0; }
  if ((ifgrab1 == 1)&&(i == grab1num)) {
    atom.rx[i] = atom.rx_p[i] + grabdxdt1 * ang * dt;
    atom.ry[i] = atom.ry_p[i] + grabdydt1 * ang * dt;
    atom.rz[i] = atom.rz_p[i] + grabdzdt1 * ang * dt;
  }
  if ((ifgrab1 == 2)&&(i == grab2num)) {
    atom.rx[i] = atom.rx_p[i] + grabdxdt2 * ang * dt;
    atom.ry[i] = atom.ry_p[i] + grabdydt2 * ang * dt;
    atom.rz[i] = atom.rz_p[i] + grabdzdt2 * ang * dt;
  }
      }
    }
  cell.s  += c1 * cell.sv + c2 * cell.sa + c3 * cell.sb + c4 * cell.sc;
  cell.sv += c1 * cell.sa + c2 * cell.sb + c3 * cell.sc;
  cell.sa += c1 * cell.sb + c2 * cell.sc;
  cell.sb += c1 * cell.sc;
  if (atom.QC) {
  for (int iel=1; iel<=atom.nelem; iel++)
    {
      mk_helmat(iel, atom.rx_p, atom.ry_p, atom.rz_p, helmat_p);
      mk_helmat(iel, atom.rx, atom.ry, atom.rz, helmat);
      inverse(helmat_p, hinelmat_p);
      for (int i=1; i<=atom.natom; i++)
	{
	  if ((atom.repatom[i] == 0)&&(atom.elem_id[i] == iel)) {
	    qx = q_x(i, iel, atom.rx_p, atom.ry_p, atom.rz_p, hinelmat_p);
	    qy = q_y(i, iel, atom.rx_p, atom.ry_p, atom.rz_p, hinelmat_p);
	    qz = q_z(i, iel, atom.rx_p, atom.ry_p, atom.rz_p, hinelmat_p);
	    atom.rx[i] = r_x(i, iel, atom.rx, atom.ry, atom.rz, helmat, qx, qy, qz);
	    atom.ry[i] = r_y(i, iel, atom.rx, atom.ry, atom.rz, helmat, qx, qy, qz);
	    atom.rz[i] = r_z(i, iel, atom.rx, atom.ry, atom.rz, helmat, qx, qy, qz);
	    atom.vx[i] = ( atom.rx[i]-atom.rx_p[i] ) / dt;
	    atom.vy[i] = ( atom.ry[i]-atom.ry_p[i] ) / dt;
	    atom.vz[i] = ( atom.rz[i]-atom.rz_p[i] ) / dt;
	  }
	}
    }
  }
}

void gear_pc_nvtpr_correc()
{
  double helmat_p[3][3], hinelmat_p[3][3], helmat[3][3];
  double qx, qy, qz, xx, yy, zz;
  double gear0, gear1, gear2, gear3, gear4;
  double c1, c2, c3, c4, cr, cv, cb, cc;
  double axi, ayi, azi, corrx, corry, corrz, corr;
  double gmat[3][3], gvmat[3][3], ginmat[3][3];
  double pt[3][3], tm[3][3], tm2[3][3], tm3[3][3];
  double hf[3][3], rec[3][3], corrm[3][3], tmp[3][3];
  double hmat_t[3][3], hvmat_t[3][3], hinmat_t[3][3]; // transpose
  gear0 = 19.0/90.0, gear1 = 3.0/4.0, gear3 = 1.0/2.0, gear4 = 1.0/12.0;
  c1 = dt; c2 = c1 * dt / 2.0; c3 = c2 * dt / 3.0; c4 = c3 * dt / 4.0;
  cr = gear0 * c2; cv = gear1 * c2 / c1; cb = gear3 * c2 / c3; cc = gear4 * c2 / c4;

  cell.wq = cell.nhmass * pow(10,nhmass_scale);

  cell.sf
    = 2.0*atom.Enkin()*cell.s/cell.wq
    - 3.0*(double)atom.natom*BL*(double)temp_set*cell.s/cell.wq
    + cell.sv*cell.sv/cell.s;
  
  for (int i=1; i<=atom.natom; i++) {
    if ((atom.QC==0)||(atom.repatom[i] == 1)) {
      axi = atom.fx[i] / atom.wm[i] - cell.sv * atom.vx[i] / cell.s;
      ayi = atom.fy[i] / atom.wm[i] - cell.sv * atom.vy[i] / cell.s;
      azi = atom.fz[i] / atom.wm[i] - cell.sv * atom.vz[i] / cell.s;
      corrx = axi - atom.ax[i];
      corry = ayi - atom.ay[i];
      corrz = azi - atom.az[i];
      atom.rx[i] += cr * corrx;
      atom.ry[i] += cr * corry;
      atom.rz[i] += cr * corrz;
      atom.vx[i] += cv * corrx;
      atom.vy[i] += cv * corry;
      atom.vz[i] += cv * corrz;
      atom.ax[i] = axi;
      atom.ay[i] = ayi;
      atom.az[i] = azi;
      atom.bx[i] += cb * corrx;
      atom.by[i] += cb * corry;
      atom.bz[i] += cb * corrz;
      atom.cx[i] += cc * corrx;
      atom.cy[i] += cc * corry;
      atom.cz[i] += cc * corrz;
    }
  }
  corr = cell.sf - cell.sa;
  cell.s  += cr * corr;
  cell.sv += cv * corr;
  cell.sa  = cell.sf;
  cell.sb += cb * corr;
  cell.sc += cc * corr;

  if (atom.QC) {
  for (int iel=1; iel<=atom.nelem; iel++)
    {
      mk_helmat(iel, atom.rx_p, atom.ry_p, atom.rz_p, helmat_p);
      mk_helmat(iel, atom.rx, atom.ry, atom.rz, helmat);
      inverse(helmat_p, hinelmat_p);
      for (int i=1; i<=atom.natom; i++)
	{
	  if ((atom.repatom[i] == 0)&&(atom.elem_id[i] == iel)) {
	    qx = q_x(i, iel, atom.rx_p, atom.ry_p, atom.rz_p, hinelmat_p);
	    qy = q_y(i, iel, atom.rx_p, atom.ry_p, atom.rz_p, hinelmat_p);
	    qz = q_z(i, iel, atom.rx_p, atom.ry_p, atom.rz_p, hinelmat_p);
	    atom.rx[i] = r_x(i, iel, atom.rx, atom.ry, atom.rz, helmat, qx, qy, qz);
	    atom.ry[i] = r_y(i, iel, atom.rx, atom.ry, atom.rz, helmat, qx, qy, qz);
	    atom.rz[i] = r_z(i, iel, atom.rx, atom.ry, atom.rz, helmat, qx, qy, qz);
	    atom.vx[i] = ( atom.rx[i]-atom.rx_p[i] ) / dt;
	    atom.vy[i] = ( atom.ry[i]-atom.ry_p[i] ) / dt;
	    atom.vz[i] = ( atom.rz[i]-atom.rz_p[i] ) / dt;
	  }
	}
    }
  }
}
