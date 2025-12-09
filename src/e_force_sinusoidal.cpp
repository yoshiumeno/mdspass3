#include <iostream>
#include<string.h>
#include<stdlib.h>
#include<fstream>
#include "myheader.h"

double vsin(double rmeter, double al, double ro);
double vpsin(double rmeter, double al, double ro);

void e_force_sinusoidal()
{
  double rr, rr2, drx, dry, drz, vp0, v0;
  int j, ix, iy, iz;
  //Parameters  
  double al = 1.0;
  //double ro = 2.83;
  double ro = 2.83;

  if (rcut_f > 0) { rcut = rcut_f; }
  else { rcut = ro * 2.0; }
  rcut = ro * 2.0;
  double rcut2 = rcut * rcut * ang * ang;

  // check force and energy curves
  /*
  double xx=0.0;
  while (xx<ro*2.5) {
    printf("%f %e %e\n",xx,vpsin(xx*ang,al,ro),vsin(xx*ang,al,ro));
    xx += 0.1;
  }
  */
  
  //   Virial term reset
  cell.virx=0.0; cell.viry=0.0; cell.virz=0.0;
  for (int i=0; i<3; i++) { for (int j=0; j<3; j++) { 
      cell.dmat[i][j] = 0.0;
      for (int ii=1; ii<=atom.natom; ii++) { atom.satom[ii][i][j] = 0.0; }  } }

  for (int i=1; i<=atom.natom; i++) {
    atom.fx[i] = 0.0; atom.fy[i] = 0.0; atom.fz[i] = 0.0;  atom.epot[i] = 0.0;
  }
  bool slow = false;
  if (slow) {
  for (int i=1; i<=atom.natom; i++) {
    if ((atom.QC==0)||(atom.repatom[i]==1)) {
      if (book.alistnum[i]>0) {
	for (int k=1; k<=book.alistnum[i]; k++) {
	  j  = book.alist[i][k][0];
	  ix = book.alist[i][k][1]; iy = book.alist[i][k][2]; iz = book.alist[i][k][3];
	  rr2 = atom.Dist2(i,j,ix,iy,iz);
	  //	if ((rr2 < rcut2)&&(rr2>1.0e-30)) {
	  if (rr2 < rcut2) {
	    rr = sqrt(rr2);
	    drx = atom.Dx(i,j,ix,iy,iz); dry = atom.Dy(i,j,ix,iy,iz); drz = atom.Dz(i,j,ix,iy,iz);
	    atom.fx[i] = atom.fx[i]+vpsin(rr,al,ro)/rr*drx;
	    atom.fy[i] = atom.fy[i]+vpsin(rr,al,ro)/rr*dry;
	    atom.fz[i] = atom.fz[i]+vpsin(rr,al,ro)/rr*drz;
	    atom.epot[i]=atom.epot[i]+vsin(rr,al,ro)/2.0;
	    double ad1 = vpsin(rr,al,ro)/rr/2.0; //20191210
	    atom.satom[i][0][0] += ad1*drx*drx;
	    atom.satom[i][0][1] += ad1*dry*drx;
	    atom.satom[i][1][1] += ad1*dry*dry;
	    atom.satom[i][0][2] += ad1*drz*drx;
	    atom.satom[i][1][2] += ad1*drz*dry;
	    atom.satom[i][2][2] += ad1*drz*drz;
	  }
	}
      }
    }
    //    printf("SLOW %d %20.15e %20.15e\n",i,atom.fx[i],rr2);
  }
  //  if (istep==0) {exit(0);}
  } else {

  for (int i=1; i<=atom.natom; i++) {
    if (book.alistnum[i]>0) {
      for (int k=1; k<=book.alistnum[i]; k++) {
	j  = book.alist[i][k][0];
	if (j>=i) {
	  if ((atom.QC==0)||(atom.repatom[i]==1)||(atom.repatom[j]==1)) {
	    ix = book.alist[i][k][1]; iy = book.alist[i][k][2]; iz = book.alist[i][k][3];
	    // rr2 = atom.Dist2(i,j,ix,iy,iz);
	    drx = atom.rx[j]+cell.hmat[0][0]*ix+cell.hmat[0][1]*iy+cell.hmat[0][2]*iz - atom.rx[i];
	    dry = atom.ry[j]+cell.hmat[1][0]*ix+cell.hmat[1][1]*iy+cell.hmat[1][2]*iz - atom.ry[i];
	    drz = atom.rz[j]+cell.hmat[2][0]*ix+cell.hmat[2][1]*iy+cell.hmat[2][2]*iz - atom.rz[i];
	    rr2 = ( drx*drx + dry*dry + drz*drz );
	    //	if ((rr2 < rcut2)&&(rr2>1.0e-30)) {
	    if (rr2 < rcut2) {
	      rr = sqrt(rr2);
	      // drx = atom.Dx(i,j,ix,iy,iz); dry = atom.Dy(i,j,ix,iy,iz); drz = atom.Dz(i,j,ix,iy,iz);
	      vp0 = vpsin(rr,al,ro); v0 = vsin(rr,al,ro); double vp0rr=vp0/rr;
	      if (j==i) { vp0rr /= 2.0; }
	      if ((atom.QC==0)||(atom.repatom[i]==1)) {
		atom.fx[i] = atom.fx[i]+vp0rr*drx;
		atom.fy[i] = atom.fy[i]+vp0rr*dry;
		atom.fz[i] = atom.fz[i]+vp0rr*drz;
		atom.epot[i]=atom.epot[i]+v0/2.0;
	      }
	      if ((atom.QC==0)||(atom.repatom[j]==1)) {
		atom.fx[j] = atom.fx[j]-vp0rr*drx;
		atom.fy[j] = atom.fy[j]-vp0rr*dry;
		atom.fz[j] = atom.fz[j]-vp0rr*drz;
		if (j != i) {
		  atom.epot[j]=atom.epot[j]+v0/2.0;}
	      }
	      if (atom.QC==0) {
		double ad1 = vp0rr/2.0;
		atom.satom[i][0][0] += ad1*drx*drx;
		atom.satom[i][0][1] += ad1*dry*drx;
		atom.satom[i][1][1] += ad1*dry*dry;
		atom.satom[i][0][2] += ad1*drz*drx;
		atom.satom[i][1][2] += ad1*drz*dry;
		atom.satom[i][2][2] += ad1*drz*drz;
		atom.satom[j][0][0] += ad1*drx*drx;
		atom.satom[j][0][1] += ad1*dry*drx;
		atom.satom[j][1][1] += ad1*dry*dry;
		atom.satom[j][0][2] += ad1*drz*drx;
		atom.satom[j][1][2] += ad1*drz*dry;
		atom.satom[j][2][2] += ad1*drz*drz;
	      }

	    }
	  }
	} //j>=i
      }
    }
    //    printf("FAST %d %20.15e %20.15e\n",i,atom.fx[i],rr2);
  }
  //  if (istep==0) {exit(0);}
  }
  // stress
  if (atom.QC==0) {
    for (int i=1; i<=atom.natom; i++) {
      atom.satom[i][1][0] = atom.satom[i][0][1];
      atom.satom[i][2][0] = atom.satom[i][0][2];
      atom.satom[i][2][1] = atom.satom[i][1][2];}
    for (int ii=1; ii<=atom.natom; ii++) {
      for (int i=0; i<3; i++) { for (int j=0; j<3; j++) { 
	  cell.dmat[i][j] = cell.dmat[i][j] - atom.satom[ii][i][j];
	} } }
    cell.virx = cell.dmat[0][0]; cell.viry = cell.dmat[1][1]; cell.virz = cell.dmat[2][2];
    cell.volume = cell.Getvolume();
    for (int ii=1; ii<=atom.natom; ii++) {
      for (int i=0; i<3; i++) { for (int j=0; j<3; j++) {
	  atom.satom[ii][i][j] = atom.satom[ii][i][j] * (double)atom.natom / cell.volume;
	} } }
  }
  atom.epotsum=0.0;
  for (int i=1; i<=atom.natom; i++) {
    atom.epotsum += atom.epot[i]; }
}


double vsin(double rmeter, double al, double ro)
{
  double vval, rang;
  rang = rmeter * 1.0e10;
  //  al = 1.0;  ro = 2.83;
  //  eV = 1.6021892e-19;
  vval = al * ro / M_PI * (-1 + cos(M_PI/ro * rang)) *eV;
  if (rang < ro) vval = vval*5.0 + al*ro/M_PI*2.0*(5.0-1.0)*eV;
  //if (rang > ro * 2.0) vval = 0.0;
  return vval;
}

double vpsin(double rmeter, double al, double ro)
{
  double vpval, rang;
  rang = rmeter * 1.0e10;
  //  al = 1.0;  ro = 2.83;
  //  eV = 1.6021892e-19;
  vpval = -al * sin(M_PI/ro * rang) *eV*1.0e10;
  if (rang < ro ) vpval = vpval*5.0;
  //if (rang > ro * 2.0) vpval = 0;
  return vpval;
}


