#include <iostream>
#include<string.h>
#include<stdlib.h>
#include<fstream>
#include "myheader.h"

double v(double rmeter);
double vp(double rmeter);

/*
void e_force_morse_nobook()
{
  double rr, rr2, drx, dry, drz;

  for (int i=1; i<=atom.natom; i++) {
    atom.fx[i] = 0.0; atom.fy[i] = 0.0; atom.fz[i] = 0.0;  atom.epot[i] = 0.0;
  }
  for (int i=1; i<=atom.natom; i++) {
    for (int j=1; j<=atom.natom; j++) {
      if (j!=i) {
	drx = atom.rx[i] - atom.rx[j]; dry = atom.ry[i] - atom.ry[j]; drz = atom.rz[i] - atom.rz[j];
	rr2 = drx*drx+dry*dry+drz*drz;
	if (rr2 < rcut2) {
	  rr = sqrt(rr2);
	  atom.fx[i] = atom.fx[i]-vp(rr)/rr*drx;
	  atom.fy[i] = atom.fy[i]-vp(rr)/rr*dry;
	  atom.fz[i] = atom.fz[i]-vp(rr)/rr*drz;
	  atom.epot[i]=atom.epot[i]+v(rr)/2.0;
	}
      }
    }
  }
}
*/

 /*
void e_force_morse()
{
  double rr, rr2, drx, dry, drz;
  int j, ix, iy, iz;

  for (int i=1; i<=atom.natom; i++) {
    atom.fx[i] = 0.0; atom.fy[i] = 0.0; atom.fz[i] = 0.0;  atom.epot[i] = 0.0;
  }
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
	    atom.fx[i] = atom.fx[i]+vp(rr)/rr*drx;
	    atom.fy[i] = atom.fy[i]+vp(rr)/rr*dry;
	    atom.fz[i] = atom.fz[i]+vp(rr)/rr*drz;
	    atom.epot[i]=atom.epot[i]+v(rr)/2.0;
	  }
	}
      }
    }
    printf("%d %20.15e %20.15e\n",i,atom.fx[i],rr2);
  }
  exit(0);
}
 */
void e_force_morse()
{
  double rr, rr2, drx, dry, drz, vp0, v0;
  int j, ix, iy, iz;
  if (rcut_f > 0) { rcut = rcut_f; }
  else { rcut = 8.0e0; }
  double rcut2 = rcut * rcut * ang * ang;

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
	    atom.fx[i] = atom.fx[i]+vp(rr)/rr*drx;
	    atom.fy[i] = atom.fy[i]+vp(rr)/rr*dry;
	    atom.fz[i] = atom.fz[i]+vp(rr)/rr*drz;
	    atom.epot[i]=atom.epot[i]+v(rr)/2.0;
	    double ad1 = vp(rr)/rr/2.0; //20191210
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
	      vp0 = vp(rr); v0 = v(rr); double vp0rr=vp0/rr;
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


double v(double rmeter)
{
  double vval, rang, ep, al, ro;
  rang = rmeter * 1.0e10;
  //  ep = 0.2703;  al = 1.1646;  ro = 3.253;
  ep = 0.3429;  al = 1.3588;  ro = 2.866;
  //  eV = 1.6021892e-19;
  vval=ep*(exp(-2.0*al*(rang-ro))-2.0*exp(-al*(rang-ro))) *eV;
  return vval;
}

double vp(double rmeter)
{
  double vpval, rang, ep, al, ro;
  rang = rmeter * 1.0e10;
  //  ep = 0.2703;  al = 1.1646;  ro = 3.253;
  ep = 0.3429;  al = 1.3588;  ro = 2.866;
  //  eV = 1.6021892e-19;
  vpval=-2.0*al*ep*(exp(-2.0*al*(rang-ro))-exp(-al*(rang-ro))) *eV*1.0e10;
  return vpval;
}

