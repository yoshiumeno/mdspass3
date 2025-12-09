#include <iostream>
#include <string.h>
#include <string>
#include <stdlib.h>
#include <fstream>
#include "myheader.h"

#define dsq(a) ((a)*(a))
void tersoff_alloc();
void tersoff_delete();
void tersoff_setparam(const char* arg);
double fcut(double r, double rr, double ss);
double fcutd(double r, double rr, double ss);
double v(double rmeter);
double vp(double rmeter);
void resetmat(double a[3][3]);
int tertyp(int at);

void e_force_tersoff()
{
  tersoff.debug_large = false; // to check large system mode
  double rr, rr2, drx, dry, drz, vp0, v0;
  int j, ix, iy, iz;
  rcut = 4.0e0; double rcut2 = rcut * rcut;

  double pi = M_PI, eps = 1.0e-20, ev2j = 1.6021892e-19, j2ev = 1.0/ev2j;

  if (tersoff.initialize) {
    std::cout<<"Initialize of the Tersoff potential..\n";
    //tersoff_alloc();
    tersoff_delete();
    //printf("bookkeep algo = %d\n",book.algo);
    //if (book.algo != 2) { printf("Error in Tersoff: Bookkeep algo is not 2 (double)\n"); goto OUT; }
    //if (atom.natom > 3000) { printf("ERROR: Please use TersoffNM instead.\n"); goto OUT; }

    tersoff.large = true;
    // 'small' mode turned out to be wrong and unnecessary... (20210713)
    //for (int i=1; i<=atom.natom; i++) {
    //  if (book.alistnum[i]>0) {
    //    int jp, j;
    //    for (int k0=1; k0<=book.alistnum[i]; k0++) {
    //      if (k0>1) jp = j;
    //      j = book.alist[i][k0][0];
    //      if ((k0>1)&&(jp == j)) {
    //        tersoff.large = false;
    //      }
    //    }
    //  }
    //}

    if (tersoff.large)  { printf("## Tersoff: large system mode\n"); }
    if (!tersoff.large) { printf("## Tersoff: small system mode\n"); }
    if (tersoff.nocutoff)  {
      printf("## Tersoff: cut-off removal mode (Shenderova PRB 61 3877)\n");
      printf("   You set 'tersoff_nocutoff yes'. Cut-off radius is removed when step = 1 \n");
      printf("   For this mode, please set very large nbk (e.g. 1000000000)\n");
      printf("   At step = 0, neighbor list is created using frc (in SETDAT) as cut-off radius\n");
      printf("   Due to very large nbk, neighbor list will not be renewed\n");
      printf("   All the atom pairs listed as neighbors at initial will be calculated forever\n");
    }

    tersoff_setparam(atom.potential_arg);
    tersoff.initialize = false;

    std::cout<<"done"<<std::endl;
  }
 OUT:
  int istat = 1;
  double rc = 0;
  for (int i=0; i<tersoff.ntype; i++) { if (rc < tersoff.terss[i]) rc = tersoff.terss[i]; }
  if ((tersoff.nocutoff)&&(istep>0)) rc = 100.0*ang;
  rcut = rc * 1e10;
  double rc2 = rc*rc;
  if ((tersoff.nocutoff)&&(istep>0)) {
    for (int i=0; i<tersoff.nptype; i++) {
      tersoff.terss[i] = rc; tersoff.terrr[i] = rc; } }

  int kns, kn, kns0;
  double rix, riy, riz, rjx, rjy, rjz, fij, dfij, rkx, rky, rkz;
  double dxij, dyij, dzij;
  int typj, typi, typk, pairij, pairik, pairjk;

  //   Virial term reset
  cell.virx=0.0; cell.viry=0.0; cell.virz=0.0;
  for (int i=0; i<3; i++) { for (int j=0; j<3; j++) { 
      cell.dmat[i][j] = 0.0;
      for (int ii=1; ii<=atom.natom; ii++) { atom.satom[ii][i][j] = 0.0; }  } }
  //========================================================
  //   Loop for b and z
  if (tersoff.debug_large||(!tersoff.large)) {
    for (int i=1; i<=atom.natom; i++) {
      for (int j=1; j<=atom.natom; j++) {
	tersoff.zmat[i][j] = 0; tersoff.b[i][j] = 0; } } }
  if (tersoff.debug_large||tersoff.large) {
    for (int i=1; i<=atom.natom; i++) {
      for (int j=1; j<=tersoff.maxnei; j++) {
	tersoff.z_ij[i][j]=0; tersoff.z_ji[i][j]=0; tersoff.b_ij[i][j]=0; tersoff.b_ji[i][j]=0; } } }

  //kns = 0;
  for (int i=1; i<=atom.natom; i++) { // bloop
    typi = tersoff.typen[i];
    if (typi < 0) continue;
    rix = atom.rx[i]; riy = atom.ry[i]; riz = atom.rz[i];
    if (book.alistnum[i]>0) {
      if (book.alistnum[i]>tersoff.maxnei) {
	printf("## Tersoff error: Increase maxnei (now %d) in myclass.cpp\n",tersoff.maxnei);
	printf("##                or set smaller frc in SETDAT\n");
	return;
      }
      for (int k0=1; k0<=book.alistnum[i]; k0++) { // 1200
	int j = book.alist[i][k0][0];
	int ix = book.alist[i][k0][1]; int iy = book.alist[i][k0][2]; int iz = book.alist[i][k0][3];
	typj = tersoff.typen[j];
  if (typj < 0) continue;
	double rij2 = atom.Dist2(i,j,ix,iy,iz);
	if ((rij2 < rc2)&&(rij2 > eps)) {
	  double rij = sqrt(rij2);
	  rjx = atom.rx[i] + atom.Dx(i,j,ix,iy,iz);
	  rjy = atom.ry[i] + atom.Dy(i,j,ix,iy,iz);
	  rjz = atom.rz[i] + atom.Dz(i,j,ix,iy,iz);
	  pairij = tersoff.ptype[typi][typj];
	  fij = fcut(rij,tersoff.terrr[pairij],tersoff.terss[pairij]);
	  if (fij < eps) { continue; }
	  //--------------calc bij---------------------
	  for (int k1=1; k1<=book.alistnum[i]; k1++) {
	    int k = book.alist[i][k1][0];
	    typk = tersoff.typen[k];
      if (typk < 0) continue;
	    ix = book.alist[i][k1][1]; iy = book.alist[i][k1][2]; iz = book.alist[i][k1][3];
	    rkx = atom.rx[i] + atom.Dx(i,k,ix,iy,iz);
	    rky = atom.ry[i] + atom.Dy(i,k,ix,iy,iz);
	    rkz = atom.rz[i] + atom.Dz(i,k,ix,iy,iz);
	    double rik2 = atom.Dist2(i,k,ix,iy,iz);
	    if (rik2 < eps) { continue; }
	    double rik = sqrt(rik2);
	    pairik = tersoff.ptype[typi][typk];
	    double fik = fcut(rik,tersoff.terrr[pairik],tersoff.terss[pairik]);
	    if (fik < eps) { continue; }
	    double rjk2 = (rkx-rjx)*(rkx-rjx)+(rky-rjy)*(rky-rjy)+(rkz-rjz)*(rkz-rjz);
	    if (rjk2 < eps) { continue; }
	    double rjk = sqrt(rjk2);
	    double costh = (rik2+rij2-rjk2)/(2.0*rik*rij);
	    double g;
	    if (tersoff.functype < 2) {
	      g = 1.0+tersoff.terc2[typi]/tersoff.terd2[typi]
		-tersoff.terc2[typi]/
		(tersoff.terd2[typi]+dsq(tersoff.terh[typi]-costh));
	    } else {
	      g = 1.0+tersoff.terc2[pairik]/tersoff.terd2[pairik]
		-tersoff.terc2[pairik]/
		(tersoff.terd2[pairik]+dsq(tersoff.terh[pairik]-costh));
	    }	      
	    double zmatp;
	    if (tersoff.functype == 0) {
	      zmatp = fik*g*exp(tersoff.termum[pairik]*pow(rij-rik,tersoff.term[pairik]));
	    } else { 
	      zmatp = fik*g*tersoff.terw[pairik];
	    }
	    if (tersoff.debug_large||(!tersoff.large))
	      tersoff.zmat[i][j] = tersoff.zmat[i][j] + zmatp;
	    if (tersoff.debug_large||tersoff.large)
	      tersoff.z_ij[i][k0] += zmatp;
	  } // end of k1 loop
	  //if (tersoff.b[i][j] != 0.0) {
	  //  printf("This Tersoff routine does not accept too small cells\n");
	  //  return; }
	  if (tersoff.debug_large||(!tersoff.large)) {
	    tersoff.b[i][j] 
	      = pow(1.0+pow(tersoff.terbeta[typi]*tersoff.zmat[i][j],tersoff.tern[typi])
		    ,-0.5/tersoff.tern[typi]);
	    if (tersoff.functype > 0) { tersoff.b[i][j] *= tersoff.terchi[pairij]; } }
	  if (tersoff.debug_large||tersoff.large) {
	    tersoff.b_ij[i][k0]
	      = pow(1.0+pow(tersoff.terbeta[typi]*tersoff.z_ij[i][k0],tersoff.tern[typi])
		    ,-0.5/tersoff.tern[typi]);
	    if (tersoff.functype > 0) { tersoff.b_ij[i][k0] *= tersoff.terchi[pairij]; } }
	} // 1201
      } // 1200
    }
  } // bloop
  /*
  if (!tersoff.large) {
    printf("zmat---> %e\n",tersoff.zmat[1][2]);
    printf("b   ---> %e\n",tersoff.b[1][2]);
  } else {
    printf("z_ij---> %e\n",tersoff.z_ij[1][2]);
    printf("b_ij---> %e\n",tersoff.b_ij[1][2]);
  }
 */
  // create z_ji and b_ji
  if (tersoff.debug_large||tersoff.large) {
    for (int i=1; i<=atom.natom; i++) {
      if (book.alistnum[i]>0) {
	for (int k0=1; k0<=book.alistnum[i]; k0++) {
	  int j = book.alist[i][k0][0];
	  
	  if (book.alistnum[j]>0) {
	    for (int k1=1; k1<=book.alistnum[j]; k1++) {
	      int k = book.alist[j][k1][0];
	      if (k == i) {
		tersoff.z_ji[i][k0] = tersoff.z_ij[j][k1];
		tersoff.b_ji[i][k0] = tersoff.b_ij[j][k1];
	      }
	    }
	  }
	}
      }
    }
  }

  // Loop
  atom.epotsum = 0.0;
  for (int i=1; i<=atom.natom; i++) {
    atom.fx[i] = 0.0; atom.fy[i] = 0.0; atom.fz[i] = 0.0;  atom.epot[i] = 0.0;
  }
  //kns = 0;
  for (int i=1; i<=atom.natom; i++) { // iloop
    typi = tersoff.typen[i];
    if (typi < 0) continue;
    rix = atom.rx[i]; riy = atom.ry[i]; riz = atom.rz[i];
    if (book.alistnum[i]>0) {
      for (int k0=1; k0<=book.alistnum[i]; k0++) { // 200
	int j = book.alist[i][k0][0];
	int ix = book.alist[i][k0][1]; int iy = book.alist[i][k0][2]; int iz = book.alist[i][k0][3];
	typj = tersoff.typen[j];
  if (typj < 0) continue;
	double rij2 = atom.Dist2(i,j,ix,iy,iz);
	if ((rij2 < rc2)&&(rij2 > eps)) {
	  double rij = sqrt(rij2);
	  rjx = atom.rx[i] + atom.Dx(i,j,ix,iy,iz);
	  rjy = atom.ry[i] + atom.Dy(i,j,ix,iy,iz);
	  rjz = atom.rz[i] + atom.Dz(i,j,ix,iy,iz);
	  pairij = tersoff.ptype[typi][typj];
	  fij = fcut(rij,tersoff.terrr[pairij],tersoff.terss[pairij]);
	  dfij = fcutd(rij,tersoff.terrr[pairij],tersoff.terss[pairij]);
	  if (fij < eps) { continue; }
	  double expmurij=exp(-tersoff.termu[pairij]*rij);
	  double ushtij;
	  if (!tersoff.large) {
	    ushtij=-fij*tersoff.b[i][j]*tersoff.terbb[pairij]*expmurij;
	  } else {
	    ushtij=-fij*tersoff.b_ij[i][k0]*tersoff.terbb[pairij]*expmurij;
	  }
	  double urepij=fij*tersoff.teraa[pairij]*exp(-tersoff.terlambda[pairij]*rij);
	  double uij=urepij+ushtij;
	  atom.epot[i]=atom.epot[i]+uij*0.5;
	  // Atomic stress
	  double ad1 = 0.5*( (-tersoff.terlambda[pairij])*urepij + dfij/fij*urepij )/rij;
	  double adx=rjx-rix; double ady=rjy-riy; double adz=rjz-riz;
	  atom.satom[i][0][0] += ad1*adx*adx;
	  atom.satom[i][0][1] += ad1*ady*adx;
	  atom.satom[i][1][1] += ad1*ady*ady;
	  atom.satom[i][0][2] += ad1*adz*adx;
	  atom.satom[i][1][2] += ad1*adz*ady;
	  atom.satom[i][2][2] += ad1*adz*adz;
	  //
	  double dvdz;
	  if (!tersoff.large) {
	    dvdz=tersoff.terbb[pairij]*fij*expmurij*tersoff.b[i][j]
	      *pow(tersoff.terbeta[typi],tersoff.tern[typi])
	      *pow(tersoff.zmat[i][j],tersoff.tern[typi])/2.0
	      /(1.0+pow(tersoff.terbeta[typi]*tersoff.zmat[i][j],tersoff.tern[typi]))
	      /tersoff.zmat[i][j];
	  } else {
	    dvdz=tersoff.terbb[pairij]*fij*expmurij*tersoff.b_ij[i][k0]
	      *pow(tersoff.terbeta[typi],tersoff.tern[typi])
	      *pow(tersoff.z_ij[i][k0],tersoff.tern[typi])/2.0
	      /(1.0+pow(tersoff.terbeta[typi]*tersoff.z_ij[i][k0],tersoff.tern[typi]))
	      /tersoff.z_ij[i][k0];
	  }
	  if (tersoff.debug_large) {
	    printf("ij %d %d %d  %e %e\n",i,j,k0,tersoff.zmat[i][j],tersoff.z_ij[i][k0]);
	    printf("ji %d %d %d  %e %e\n",i,j,k0,tersoff.zmat[j][i],tersoff.z_ji[i][k0]); }
	  double ffx = 0.0; double ffy = 0.0; double ffz = 0.0;
	  dxij=rix-rjx; dyij=riy-rjy; dzij=riz-rjz; //YU201406xx
	  for (int k1=1; k1<=book.alistnum[i]; k1++) { // kloop (300)
	    int k = book.alist[i][k1][0];
	    typk = tersoff.typen[k];
      if (typk < 0) continue;
	    pairik = tersoff.ptype[typi][typk];
	    ix = book.alist[i][k1][1]; iy = book.alist[i][k1][2]; iz = book.alist[i][k1][3];
	    rkx = atom.rx[i] + atom.Dx(i,k,ix,iy,iz);
	    rky = atom.ry[i] + atom.Dy(i,k,ix,iy,iz);
	    rkz = atom.rz[i] + atom.Dz(i,k,ix,iy,iz);
	    double rik2 = atom.Dist2(i,k,ix,iy,iz);
	    if (rik2 < eps) { continue; }
	    double rik = sqrt(rik2);
	    double fik = fcut(rik,tersoff.terrr[pairik],tersoff.terss[pairik]);
	    if (fik < eps) { continue; }
	    double rjk2 = (rkx-rjx)*(rkx-rjx)+(rky-rjy)*(rky-rjy)+(rkz-rjz)*(rkz-rjz);
	    if (rjk2 < eps) { continue; }
	    double rjk = sqrt(rjk2);
	    double costh = (rik2+rij2-rjk2)/(2.0*rik*rij);
	    double g;
	    if (tersoff.functype < 2) {
	      g = 1.0+tersoff.terc2[typi]/tersoff.terd2[typi]
		-tersoff.terc2[typi]/
		(tersoff.terd2[typi]+dsq(tersoff.terh[typi]-costh));
	    } else {
	      g = 1.0+tersoff.terc2[pairik]/tersoff.terd2[pairik]
		-tersoff.terc2[pairik]/
		(tersoff.terd2[pairik]+dsq(tersoff.terh[pairik]-costh));
	    }
	    double dzdr1, dzdc; // dz/dr_ij, dz/dtheta
	    if (tersoff.functype == 0) {
	      dzdr1=(double)tersoff.term[pairik]*tersoff.termum[pairik]
		*pow(rij-rik,tersoff.term[pairik]-1)
		*fik*exp(tersoff.termum[pairik]*pow(rij-rik,tersoff.term[pairik]))
		*g;
	      dzdc=fik*exp(tersoff.termum[pairik]*pow(rij-rik,tersoff.term[pairik]))
		*(-2.0*tersoff.terc2[typi]*(tersoff.terh[typi]-costh))
		/dsq(tersoff.terd2[typi]+dsq(tersoff.terh[typi]-costh) );
	    } else if (tersoff.functype == 1) {
	      dzdr1=0.0;
	      dzdc=fik*tersoff.terw[pairik]
		*(-2.0*tersoff.terc2[typi]*(tersoff.terh[typi]-costh))
		/dsq(tersoff.terd2[typi]+dsq(tersoff.terh[typi]-costh) );
	    } else {
	      dzdr1=0.0;
	      dzdc=fik*tersoff.terw[pairik]
		*(-2.0*tersoff.terc2[pairik]*(tersoff.terh[pairik]-costh))
		/dsq(tersoff.terd2[pairik]+dsq(tersoff.terh[pairik]-costh) );
	    }
	    double dzdr2; // dz/dr_ik
	    if (tersoff.functype == 0) {
	      dzdr2=fcutd(rik,tersoff.terrr[pairik],tersoff.terss[pairik])
		*exp(tersoff.termum[pairik]*pow(rij-rik,tersoff.term[pairik]))
		*g
		-dzdr1;
		//-(double)tersoff.term[pairik]*tersoff.termum[pairik]
		//*pow(rij-rik,tersoff.term[pairik]-1)
		//*fik*exp(tersoff.termum[pairik]*pow(rij-rik,tersoff.term[pairik]))
		//*g;
	    } else { // dz/dr_ik
	      dzdr2 = fcutd(rik,tersoff.terrr[pairik],tersoff.terss[pairik])
		*tersoff.terw[pairik]*g; //YU
	    }
	    double dxki=rkx-rix; double dyki=rky-riy; double dzki=rkz-riz;
	    double dcdrx=1.0/rik*(dxij/rij+dxki*costh/rik)
	      -1.0/rij*(dxki/rik+dxij*costh/rij);
	    double dcdry=1.0/rik*(dyij/rij+dyki*costh/rik)
	      -1.0/rij*(dyki/rik+dyij*costh/rij);
	    double dcdrz=1.0/rik*(dzij/rij+dzki*costh/rik)
	      -1.0/rij*(dzki/rik+dzij*costh/rij);
	    ffx=ffx-dvdz*(dzdr1*dxij/rij-dzdr2*dxki/rik+dzdc*dcdrx)/2.0;
	    ffy=ffy-dvdz*(dzdr1*dyij/rij-dzdr2*dyki/rik+dzdc*dcdry)/2.0;
	    ffz=ffz-dvdz*(dzdr1*dzij/rij-dzdr2*dzki/rik+dzdc*dcdrz)/2.0;
	    // Atomic stress
	    double ad1 = 0.5*dvdz*dzdr1/rij;
	    double ad2 = 0.5*dvdz*dzdr2/rik;
	    atom.satom[i][0][0] += ad1*dxij*dxij+ad2*dxki*dxki;
	    atom.satom[i][0][1] += ad1*dyij*dxij+ad2*dyki*dxki;
	    atom.satom[i][1][1] += ad1*dyij*dyij+ad2*dyki*dyki;
	    atom.satom[i][0][2] += ad1*dzij*dxij+ad2*dzki*dxki;
	    atom.satom[i][1][2] += ad1*dzij*dyij+ad2*dzki*dyki;
	    atom.satom[i][2][2] += ad1*dzij*dzij+ad2*dzki*dzki;
	    double dxjk=rkx-rjx; double dyjk=rky-rjy; double dzjk=rkz-rjz;
	    double ad0 = 0.5*dvdz*dzdc;
	    ad1 = -1.0/rik/rij;
	    ad2 = 1.0/rik/rij-costh/rik2;
	    double ad3 = 1.0/rik/rij-costh/rij2;
	    atom.satom[i][0][0] += ad0*(ad1*dxjk*dxjk+ad2*dxki*dxki+ad3*dxij*dxij);
	    atom.satom[i][0][1] += ad0*(ad1*dyjk*dxjk+ad2*dyki*dxki+ad3*dyij*dxij);
	    atom.satom[i][1][1] += ad0*(ad1*dyjk*dyjk+ad2*dyki*dyki+ad3*dyij*dyij);
	    atom.satom[i][0][2] += ad0*(ad1*dzjk*dxjk+ad2*dzki*dxki+ad3*dzij*dxij);
	    atom.satom[i][1][2] += ad0*(ad1*dzjk*dyjk+ad2*dzki*dyki+ad3*dzij*dyij);
	    atom.satom[i][2][2] += ad0*(ad1*dzjk*dzjk+ad2*dzki*dzki+ad3*dzij*dzij);
	    //
	  } // kloop (300)
	  double dvdr1, dvdr2;
	  if (!tersoff.large) {
	    dvdr1=tersoff.teraa[pairij]*exp(-tersoff.terlambda[pairij]*rij)
	      *(dfij-tersoff.terlambda[pairij]*fij)
	      -tersoff.b[i][j]*tersoff.terbb[pairij]*expmurij
	      *(dfij-tersoff.termu[pairij]*fij);
	    dvdr2=tersoff.teraa[pairij]*exp(-tersoff.terlambda[pairij]*rij)
	      *(dfij-tersoff.terlambda[pairij]*fij)
	      -tersoff.b[j][i]*tersoff.terbb[pairij]*expmurij
	      *(dfij-tersoff.termu[pairij]*fij);
	  } else {
	    dvdr1=tersoff.teraa[pairij]*exp(-tersoff.terlambda[pairij]*rij)
	      *(dfij-tersoff.terlambda[pairij]*fij)
	      -tersoff.b_ij[i][k0]*tersoff.terbb[pairij]*expmurij
	      *(dfij-tersoff.termu[pairij]*fij);
	    dvdr2=tersoff.teraa[pairij]*exp(-tersoff.terlambda[pairij]*rij)
	      *(dfij-tersoff.terlambda[pairij]*fij)
	      -tersoff.b_ji[i][k0]*tersoff.terbb[pairij]*expmurij
	      *(dfij-tersoff.termu[pairij]*fij);
	  }
	  atom.fx[i]=atom.fx[i]-((dvdr1+dvdr2)*dxij/rij)/2.0;
	  atom.fy[i]=atom.fy[i]-((dvdr1+dvdr2)*dyij/rij)/2.0;
	  atom.fz[i]=atom.fz[i]-((dvdr1+dvdr2)*dzij/rij)/2.0;
	  
	  atom.fx[i]=atom.fx[i]+ffx;
	  atom.fy[i]=atom.fy[i]+ffy;
	  atom.fz[i]=atom.fz[i]+ffz;
	  // Atomic stress
	  if (!tersoff.large) {
	    ad1 = 0.5*(fij*tersoff.termu[pairij]-dfij)*tersoff.terbb[pairij]*
	      (tersoff.b[i][j]+tersoff.b[j][i])/2.0*expmurij/rij;
	  } else {
	    ad1 = 0.5*(fij*tersoff.termu[pairij]-dfij)*tersoff.terbb[pairij]*
	      (tersoff.b_ij[i][k0]+tersoff.b_ji[i][k0])/2.0*expmurij/rij;
	  }
	  adx=rjx-rix; ady=rjy-riy; adz=rjz-riz;
	  atom.satom[i][0][0] += ad1*adx*adx;
	  atom.satom[i][0][1] += ad1*ady*adx;
	  atom.satom[i][1][1] += ad1*ady*ady;
	  atom.satom[i][0][2] += ad1*adz*adx;
	  atom.satom[i][1][2] += ad1*adz*ady;
	  atom.satom[i][2][2] += ad1*adz*adz;
	  //
	} // 201
      } // 200
    }
  } // iloop
  
  // u_ji/r_i
  //kns = 0;
  for (int i=1; i<=atom.natom; i++) { // iloop2
    typi = tersoff.typen[i];
    if (typi < 0) continue;
    rix = atom.rx[i]; riy = atom.ry[i]; riz = atom.rz[i];
    if (book.alistnum[i]>0) {
      for (int k0=1; k0<=book.alistnum[i]; k0++) { // 2200
	int j = book.alist[i][k0][0];
	int ix = book.alist[i][k0][1]; int iy = book.alist[i][k0][2]; int iz = book.alist[i][k0][3];  
	typj = tersoff.typen[j];
  if (typj < 0) continue;
	double rij2 = atom.Dist2(i,j,ix,iy,iz);
	if ((rij2 < rc2)&&(rij2 > eps)) {
	  double rij = sqrt(rij2);
	  rjx = atom.rx[i] + atom.Dx(i,j,ix,iy,iz);
	  rjy = atom.ry[i] + atom.Dy(i,j,ix,iy,iz);
	  rjz = atom.rz[i] + atom.Dz(i,j,ix,iy,iz);
	  pairij = tersoff.ptype[typi][typj];
	  fij = fcut(rij,tersoff.terrr[pairij],tersoff.terss[pairij]);
	  dfij = fcutd(rij,tersoff.terrr[pairij],tersoff.terss[pairij]);
	  if (fij < eps) { continue; }
	  double ffx = 0.0; double ffy = 0.0; double ffz = 0.0;
	  double dvdz0;
	  if (!tersoff.large) {
	    dvdz0=tersoff.terbb[pairij]*fij
	      *exp(-tersoff.termu[pairij]*rij)*tersoff.b[j][i]
	      *pow(tersoff.terbeta[typj],tersoff.tern[typj])
	      *pow(tersoff.zmat[j][i],tersoff.tern[typj])/2.0
	      /(1.0+pow(tersoff.terbeta[typj]*tersoff.zmat[j][i],tersoff.tern[typj]))
	      /tersoff.zmat[j][i];
	  } else {
	    dvdz0=tersoff.terbb[pairij]*fij
	      *exp(-tersoff.termu[pairij]*rij)*tersoff.b_ji[i][k0]
	      *pow(tersoff.terbeta[typj],tersoff.tern[typj])
	      *pow(tersoff.z_ji[i][k0],tersoff.tern[typj])/2.0
	      /(1.0+pow(tersoff.terbeta[typj]*tersoff.z_ji[i][k0],tersoff.tern[typj]))
	      /tersoff.z_ji[i][k0];
	  }
	  for (int k1=1; k1<=book.alistnum[j]; k1++) { // kloop (300)
	    int k = book.alist[j][k1][0];
	    typk = tersoff.typen[k];
      if (typk < 0) continue;
	    pairik = tersoff.ptype[typi][typk];
	    pairjk = tersoff.ptype[typj][typk];
	    ix = book.alist[j][k1][1]; iy = book.alist[j][k1][2]; iz = book.alist[j][k1][3];
	    rkx = rjx + atom.Dx(j,k,ix,iy,iz);
	    rky = rjy + atom.Dy(j,k,ix,iy,iz);
	    rkz = rjz + atom.Dz(j,k,ix,iy,iz);
	    double rik2 = (rkx-rix)*(rkx-rix)+(rky-riy)*(rky-riy)+(rkz-riz)*(rkz-riz);
	    if (rik2 < eps) { continue; }
	    double rik = sqrt(rik2);
	    double rjk2 = (rkx-rjx)*(rkx-rjx)+(rky-rjy)*(rky-rjy)+(rkz-rjz)*(rkz-rjz);
	    if (rjk2 < eps) { continue; }
	    double rjk = sqrt(rjk2);
	    double fjk = fcut(rjk,tersoff.terrr[pairjk],tersoff.terss[pairjk]);
	    if (fjk < eps) { continue; }
	    // u_ij/r_i
	    double costh=(rjk2+rij2-rik2)/(2.0*rjk*rij);
	    double g,g2;
	    if (tersoff.functype < 2) {
	      g = 1.0+tersoff.terc2[typj]/tersoff.terd2[typj]
		-tersoff.terc2[typj]/
		(tersoff.terd2[typj]+dsq(tersoff.terh[typj]-costh) );
	    } else {
	      g = 1.0+tersoff.terc2[pairjk]/tersoff.terd2[pairjk]
		-tersoff.terc2[pairjk]/
		(tersoff.terd2[pairjk]+dsq(tersoff.terh[pairjk]-costh));
	      g2= 1.0+tersoff.terc2[pairij]/tersoff.terd2[pairij]
		-tersoff.terc2[pairij]/
		(tersoff.terd2[pairij]+dsq(tersoff.terh[pairij]-costh));
	    }
	    double dzdr, dzdc;
	    if (tersoff.functype == 0) {
	      dzdr=tersoff.termum[pairjk]*(double)tersoff.term[pairjk]
		*pow(rij-rjk,tersoff.term[pairjk]-1)
		*fjk*exp(tersoff.termum[pairjk]*pow(rij-rjk,tersoff.term[pairjk]))
	 	*g;
	      dzdc=fjk*exp(tersoff.termum[pairjk]*pow(rij-rjk,tersoff.term[pairjk]))
		*(-2.0*tersoff.terc2[typj]*(tersoff.terh[typj]-costh))
		/dsq(tersoff.terd2[typj]+dsq(tersoff.terh[typj]-costh) );
	    } else if (tersoff.functype == 1) {
	      dzdr=0.0;
	      dzdc=fjk*tersoff.terw[pairjk]
		*(-2.0*tersoff.terc2[typj]*(tersoff.terh[typj]-costh))
		/dsq(tersoff.terd2[typj]+dsq(tersoff.terh[typj]-costh) );
	    } else {
	      dzdr=0.0;
	      dzdc=fjk*tersoff.terw[pairjk]
		*(-2.0*tersoff.terc2[pairjk]*(tersoff.terh[pairjk]-costh))
		/dsq(tersoff.terd2[pairjk]+dsq(tersoff.terh[pairjk]-costh) );
	    }
	    double ad1=(rjk2-rik2-rij2)/(2.0*rij2*rij*rjk);
	    double dcdrx=ad1*(rjx-rix)+1.0/(rij*rjk)*(rkx-rix);
	    double dcdry=ad1*(rjy-riy)+1.0/(rij*rjk)*(rky-riy);
	    double dcdrz=ad1*(rjz-riz)+1.0/(rij*rjk)*(rkz-riz);
	    dxij=rix-rjx;dyij=riy-rjy;dzij=riz-rjz;
	    ffx=ffx-dvdz0*(dzdr*dxij/rij+dzdc*dcdrx)/2.0;
	    ffy=ffy-dvdz0*(dzdr*dyij/rij+dzdc*dcdry)/2.0;
	    ffz=ffz-dvdz0*(dzdr*dzij/rij+dzdc*dcdrz)/2.0;
	    //u_jk/r_i
	    double dvdz;
	    if (!tersoff.large) {
	      dvdz=tersoff.terbb[pairjk]*fjk*exp(-tersoff.termu[pairjk]*rjk)*tersoff.b[j][k]
		*pow(tersoff.terbeta[typj],tersoff.tern[typj])
		*pow(tersoff.zmat[j][k],tersoff.tern[typj])/2.0
		/(1.0+pow(tersoff.terbeta[typj]*tersoff.zmat[j][k],tersoff.tern[typj]))
		/tersoff.zmat[j][k];
	    } else {
	      dvdz=tersoff.terbb[pairjk]*fjk*exp(-tersoff.termu[pairjk]*rjk)*tersoff.b_ij[j][k1]
		*pow(tersoff.terbeta[typj],tersoff.tern[typj])
		*pow(tersoff.z_ij[j][k1],tersoff.tern[typj])/2.0
		/(1.0+pow(tersoff.terbeta[typj]*tersoff.z_ij[j][k1],tersoff.tern[typj]))
		/tersoff.z_ij[j][k1];
	    }
	    if (tersoff.functype == 0) { // use [pairij] !
	      dzdr=-(double)tersoff.term[pairij]*tersoff.termum[pairij]
		*pow(rjk-rij,tersoff.term[pairij]-1)
		* fij*g*exp(tersoff.termum[pairij]*pow(rjk-rij,tersoff.term[pairij]))
		+dfij*g*exp(tersoff.termum[pairij]*pow(rjk-rij,tersoff.term[pairij]));
	      dzdc=fij*exp(tersoff.termum[pairij]*pow(rjk-rij,tersoff.term[pairij]))
		*(-2.0*tersoff.terc2[typj]*(tersoff.terh[typj]-costh))
		/dsq(tersoff.terd2[typj]+dsq(tersoff.terh[typj]-costh) );
	    } else if (tersoff.functype == 1) {
	      //dzdr=0.0;
	      dzdr=dfij*tersoff.terw[pairij]*g; //YU201406xx
	      dzdc=fij*tersoff.terw[pairij]
		*(-2.0*tersoff.terc2[typj]*(tersoff.terh[typj]-costh))
		/dsq(tersoff.terd2[typj]+dsq(tersoff.terh[typj]-costh) );
	    } else {
	      dzdr=dfij*tersoff.terw[pairij]*g2; //YU20161008
	      dzdc=fij*tersoff.terw[pairij]
		*(-2.0*tersoff.terc2[pairij]*(tersoff.terh[pairij]-costh))
	      	/dsq(tersoff.terd2[pairij]+dsq(tersoff.terh[pairij]-costh) );
	    }
	    double ad5=(rjk2-rik2-rij2)/2.0/rij2/rij/rjk;
	    double ad6=1.0/rij/rjk;
	    dcdrx=(ad5*(rjx-rix)+ad6*(rkx-rix));
	    dcdry=(ad5*(rjy-riy)+ad6*(rky-riy));
	    dcdrz=(ad5*(rjz-riz)+ad6*(rkz-riz));
	    dxij=rix-rjx;dyij=riy-rjy;dzij=riz-rjz;
	    ffx=ffx-dvdz*(dzdr*dxij/rij+dzdc*dcdrx)/2.0;
	    ffy=ffy-dvdz*(dzdr*dyij/rij+dzdc*dcdry)/2.0;
	    ffz=ffz-dvdz*(dzdr*dzij/rij+dzdc*dcdrz)/2.0;
	  } // 2300
	  atom.fx[i]=atom.fx[i]+ffx;
	  atom.fy[i]=atom.fy[i]+ffy;
	  atom.fz[i]=atom.fz[i]+ffz;
	} // 2201
      } //2200
    }
  } // iloop2

  /*
  for (int i=1; i<=atom.natom; i++) {
    if ((atom.QC==0)||(atom.repatom[i]==1)) {
      if (book.alistnum[i]>0) {
	for (int k=1; k<=book.alistnum[i]; k++) {
	  j  = book.alist[i][k][0];
	  ix = book.alist[i][k][1]; iy = book.alist[i][k][2]; iz = book.alist[i][k][3];
	  rr2 = atom.Dist2(i,j,ix,iy,iz);
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
  }
  */
  //Etot
  atom.epotsum=0.0;
  for (int i=1; i<=atom.natom; i++) {
    atom.epotsum += atom.epot[i]; }
  //Atomic stress
  for (int i=1; i<=atom.natom; i++) {
    atom.satom[i][1][0] = atom.satom[i][0][1];
    atom.satom[i][2][0] = atom.satom[i][0][2];
    atom.satom[i][2][1] = atom.satom[i][1][2]; }
  resetmat(cell.dmat);
  for (int i=1; i<=atom.natom; i++) {
    for (int j=0; j<3; j++) {
      for (int k=0; k<3; k++) {
	cell.dmat[j][k] -= atom.satom[i][j][k]; } } }
  cell.virx = cell.dmat[0][0];
  cell.viry = cell.dmat[1][1];
  cell.virz = cell.dmat[2][2];
  cell.volume = cell.Getvolume();
  for (int i=1; i<=atom.natom; i++) {
    for (int j=0; j<3; j++) {
      for (int k=0; k<3; k++) {
	atom.satom[i][j][k] *= (double)atom.natom / cell.volume; } } }

  if (tersoff.nocutoff) {
    rc = 0;
    for (int i=0; i<tersoff.ntype; i++) {if (rc < tersoff.terss[i]) rc = tersoff.terss[i]; }
    rcut = rc * 1e10; }
}


void tersoff_delete()
{
  if (tersoff.type)  { delete[] tersoff.type;  tersoff.type  = NULL; }
  if (tersoff.typen) { delete[] tersoff.typen; tersoff.typen = NULL; }
  if (tersoff.ptype) {
    for (int i=0; i<tersoff.ntype; i++) { delete[] tersoff.ptype[i]; }
    delete[] tersoff.ptype; tersoff.ptype = NULL; }
  if (tersoff.zmat) {
    //for (int i=0; i<atom.natom+1; i++) { delete[] tersoff.zmat[i]; }
    for (int i=0; i<tersoff.natom+1; i++) { delete[] tersoff.zmat[i]; }
    delete[] tersoff.zmat; tersoff.zmat = NULL; }
  if (tersoff.b) {
    //for (int i=0; i<atom.natom+1; i++) { delete[] tersoff.b[i]; }
    for (int i=0; i<tersoff.natom+1; i++) { delete[] tersoff.b[i]; }
    delete[] tersoff.b; tersoff.b = NULL; }
  if (tersoff.z_ij) {
    for (int i=0; i<tersoff.natom+1; i++) {
      delete[] tersoff.z_ij[i], tersoff.z_ji[i], tersoff.b_ij[i], tersoff.b_ji[i]; }
    delete[] tersoff.z_ij; delete[] tersoff.z_ji; delete[] tersoff.b_ij; delete[] tersoff.b_ji;
    tersoff.z_ij = NULL; tersoff.z_ji = NULL; tersoff.b_ij = NULL; tersoff.b_ji = NULL; }
}

void tersoff_alloc()
{
  tersoff.type  = new int[tersoff.ntype];
  tersoff.typen = new int[atom.natom+1];
  tersoff.ptype = new int*[tersoff.ntype];
  for (int i=0; i<tersoff.ntype; i++) {
    tersoff.ptype[i] = new int[tersoff.ntype]; }
  tersoff.terrr     = new double[tersoff.nptype];
  tersoff.terss     = new double[tersoff.nptype];
  tersoff.teraa     = new double[tersoff.nptype];
  tersoff.terbb     = new double[tersoff.nptype];
  tersoff.terlambda = new double[tersoff.nptype];
  tersoff.termu     = new double[tersoff.nptype];
  tersoff.termum    = new double[tersoff.nptype];
  tersoff.terbeta   = new double[tersoff.ntype];
  tersoff.tern      = new double[tersoff.ntype];
  tersoff.term      = new int[tersoff.nptype];
  if (tersoff.functype < 2) {
    tersoff.terc2     = new double[tersoff.ntype];
    tersoff.terd2     = new double[tersoff.ntype];
    tersoff.terh      = new double[tersoff.ntype];
  } else {
    tersoff.terc2     = new double[tersoff.nptype];
    tersoff.terd2     = new double[tersoff.nptype];
    tersoff.terh      = new double[tersoff.nptype];
  }
  tersoff.terw      = new double[tersoff.nptype];
  tersoff.terchi    = new double[tersoff.nptype];
  if (tersoff.debug_large||(!tersoff.large)) {
    tersoff.zmat  = new double*[atom.natom+1];
    tersoff.b     = new double*[atom.natom+1];
    for (int i=0; i<atom.natom+1; i++) {
      tersoff.zmat[i] = new double[atom.natom+1];
      tersoff.b[i]    = new double[atom.natom+1]; } }
  if (tersoff.debug_large||tersoff.large) {
    tersoff.z_ij  = new double*[atom.natom+1]; tersoff.z_ji  = new double*[atom.natom+1];
    tersoff.b_ij  = new double*[atom.natom+1]; tersoff.b_ji  = new double*[atom.natom+1];
    for (int i=0; i<atom.natom+1; i++) {
      tersoff.z_ij[i] = new double[tersoff.maxnei+1];
      tersoff.z_ji[i] = new double[tersoff.maxnei+1];
      tersoff.b_ij[i] = new double[tersoff.maxnei+1];
      tersoff.b_ji[i] = new double[tersoff.maxnei+1]; } }
  tersoff.natom = atom.natom; // store size of system (for dellocation of memory)
}

int tertyp(int at)
{
  //int i  = atom.anum[at];
  int i = at;
  for (int j=0; j<tersoff.ntype; j++) {
    if (i == tersoff.type[j]) { return j; }
  }
  return -1;
}

void tersoff_setparam(const char* arg)
{
  double pi = M_PI, eps = 1.0e-20, ev2j = 1.6021892e-19, j2ev = 1.0/ev2j;

  //Set # of species and atom type
  if ((strcmp(arg,"Si(B)")  == 0)||
      (strcmp(arg,"Si(C)")  == 0)||
      (strcmp(arg,"Si(B*)") == 0)) { // Si single system (PRB***)
    //    printf(" Si single-atom system: arg = %s  ",arg);
    tersoff.functype = 0;
    tersoff.ntype = 1; tersoff.nptype = 1; tersoff_alloc();
    tersoff.type[0] = 14;
    tersoff.ptype[0][0] = 0;
  } else if (strcmp(arg,"C_Si_Ge") == 0) { // C-Si-Ge (PRB***)
    //    printf(" multi-atom system: arg = %s  ",arg);
    tersoff.functype = 1;
    tersoff.ntype = 3; tersoff.nptype = 6; tersoff_alloc();
    tersoff.type[0] =  6;
    tersoff.type[1] = 14;
    tersoff.type[2] = 32;
    tersoff.ptype[0][0] = 0; // C-C
    tersoff.ptype[1][1] = 1; // Si-Si
    tersoff.ptype[2][2] = 2; // Ge-Ge
    tersoff.ptype[0][1] = 3; tersoff.ptype[1][0] = 3; // C-Si
    tersoff.ptype[1][2] = 4; tersoff.ptype[2][1] = 4; // Si-Ge
    tersoff.ptype[2][0] = 5; tersoff.ptype[0][2] = 5; // Ge-C
  } else if (strcmp(arg,"B_N_C") == 0) { // B-N-C Matsunaga
    //    printf(" multi-atom system: arg = %s  ",arg);
    tersoff.functype = 1;
    tersoff.ntype = 3; tersoff.nptype = 6; tersoff_alloc();
    tersoff.type[0] =  5;
    tersoff.type[1] =  7;
    tersoff.type[2] =  6;
    tersoff.ptype[0][0] = 0; // B-B
    tersoff.ptype[1][1] = 1; // N-N
    tersoff.ptype[2][2] = 2; // C-C
    tersoff.ptype[0][1] = 3; tersoff.ptype[1][0] = 3; // B-N
    tersoff.ptype[1][2] = 4; tersoff.ptype[2][1] = 4; // N-C
    tersoff.ptype[2][0] = 5; tersoff.ptype[0][2] = 5; // C-B
  } else if (strcmp(arg,"SiC_Erhart") == 0) { // Erhart-Albe
    tersoff.functype = 2;
    tersoff.ntype = 2; tersoff. nptype = 3; tersoff_alloc();
    tersoff.type[0] = 14;
    tersoff.type[1] = 6;
    tersoff.ptype[0][0] = 0; // Si-Si
    tersoff.ptype[1][1] = 1; // C-C
    tersoff.ptype[0][1] = 2; tersoff.ptype[1][0] = 2; // Si-C
  }
  

  //Set atom type
  for (int i=1; i<=atom.natom; i++) {
    tersoff.typen[i] = tertyp(atom.anum[i]);
    if (tersoff.typen[i] < 0) {
      printf("##WARNING! Atom %d (%s) is not supported in Tersoff\n",i,atom.asp[i]);
    }
  }

  if (strcmp(arg,"Si(B)") == 0) { // Si(B) (for surface structure)
    printf("Si(B) PRB38-14(1988)p.9902 (for surface)\n");
    tersoff.terrr[0]     = 2.80e-10;
    tersoff.terss[0]     = 3.20e-10;
    tersoff.teraa[0]     = 3.2647e3*ev2j;
    tersoff.terbb[0]     = 9.5373e1*ev2j;
    tersoff.terlambda[0] = 3.2394e10;
    tersoff.termu[0]     = 1.3258e10;
    tersoff.terbeta[0]   = 3.3675e-1;
    tersoff.tern[0]      = 2.2956e1;
    tersoff.term[0]      = 3;
    tersoff.terc2[0]     = dsq(4.8381e0);
    tersoff.terd2[0]     = dsq(2.0417e0);
    tersoff.terh[0]      = 0.0e0;
    tersoff.terw[0]      = 1.0e0;
    tersoff.terchi[0]    = 1.0e0;
  } else if (strcmp(arg,"Si(C)") == 0) { // Si(C) (for elastic property)
    printf("Si(C) PRB38-14(1988)p.9902 (for elastic)\n");
    tersoff.terrr[0]     = 2.70e-10;
    tersoff.terss[0]     = 3.00e-10;
    tersoff.teraa[0]     = 1.8308e3*ev2j;
    tersoff.terbb[0]     = 4.7118e2*ev2j;
    tersoff.terlambda[0] = 2.4799e10;
    tersoff.termu[0]     = 1.7322e10;
    tersoff.terbeta[0]   = 1.0999e-6;
    tersoff.tern[0]      = 7.8734e-1;
    tersoff.term[0]      = 3;
    tersoff.terc2[0]     = dsq(1.0039e5);
    tersoff.terd2[0]     = dsq(1.6218e1);
    tersoff.terh[0]      = -5.9826e-1;
    tersoff.terw[0]      = 1.0e0;
    tersoff.terchi[0]    = 1.0e0;
  } else if (strcmp(arg,"Si(B*)") == 0) { // Si(B*) (for phase transformation)
    printf("Si(B*) PRB 37-12(1988)p.6991 (for phase trans)\n");
    tersoff.terrr[0]     = 2.65e-10;
    tersoff.terss[0]     = 2.85e-10;
    tersoff.teraa[0]     = 3.2647e3*ev2j;
    tersoff.terbb[0]     = 9.5373e1*ev2j;
    tersoff.terlambda[0] = 3.2394e10;
    tersoff.termu[0]     = 1.3258e10;
    tersoff.terbeta[0]   = 3.3675e-1;
    tersoff.tern[0]      = 2.2956e1;
    tersoff.term[0]      = 3;
    tersoff.terc2[0]     = dsq(4.8381e0);
    tersoff.terd2[0]     = dsq(2.0417e0);
    tersoff.terh[0]      = 0.0e0;
    tersoff.terw[0]      = 1.0e0;
    tersoff.terchi[0]    = 1.0e0;
  } else if (strcmp(arg,"C_Si_Ge") == 0) { // C-Si-Ge
    printf("Tersoff multi (C/Si/Ge) PRB 39-8(1989)p.5566\n");
    // for C
    tersoff.terrr[0]     = 1.80e-10;
    tersoff.terss[0]     = 2.10e-10;
    tersoff.teraa[0]     = 1.3936e3*ev2j;
    tersoff.terbb[0]     = 3.467e2*ev2j;
    tersoff.terlambda[0] = 3.4879e10;
    tersoff.termu[0]     = 2.2119e10;
    tersoff.terbeta[0]   = 1.5724e-7;
    tersoff.tern[0]      = 7.2751e-1;
    tersoff.term[0]      = 3;
    tersoff.terc2[0]     = dsq(3.8049e4);
    tersoff.terd2[0]     = dsq(4.384e0);
    tersoff.terh[0]      =-5.7058e-1;
    tersoff.terw[0]      = 1.0e0;
    tersoff.terchi[0]    = 1.0e0;
    // for Si
    tersoff.terrr[1]     = 2.70e-10;
    tersoff.terss[1]     = 3.00e-10;
    tersoff.teraa[1]     = 1.8308e3*ev2j;
    tersoff.terbb[1]     = 4.7118e2*ev2j;
    tersoff.terlambda[1] = 2.4799e10;
    tersoff.termu[1]     = 1.7322e10;
    tersoff.terbeta[1]   = 1.0999e-6;
    tersoff.tern[1]      = 7.8734e-1;
    tersoff.term[1]      = 3;
    tersoff.terc2[1]     = dsq(1.0039e5);
    tersoff.terd2[1]     = dsq(1.6218e1);
    tersoff.terh[1]      =-5.9826e-1;
    tersoff.terw[1]      = 1.0e0;
    tersoff.terchi[1]    = 1.0e0;
    // for Ge
    tersoff.terrr[2]     = 2.80e-10;
    tersoff.terss[2]     = 3.10e-10;
    tersoff.teraa[2]     = 1.769e3*ev2j;
    tersoff.terbb[2]     = 4.1923e2*ev2j;
    tersoff.terlambda[2] = 2.4451e10;
    tersoff.termu[2]     = 1.7047e10;
    tersoff.terbeta[2]   = 9.0166e-7;
    tersoff.tern[2]      = 7.5627e-1;
    tersoff.term[2]      = 3;
    tersoff.terc2[2]     = dsq(1.0643e5);
    tersoff.terd2[2]     = dsq(1.5652e1);
    tersoff.terh[2]      =-4.3884e-1;
    tersoff.terw[2]      = 1.0e0;
    tersoff.terchi[2]    = 1.0e0;
    // Dissimilar atoms
    tersoff.terrr[3] = sqrt(tersoff.terrr[0]*tersoff.terrr[1]);
    tersoff.terrr[4] = sqrt(tersoff.terrr[1]*tersoff.terrr[2]);
    tersoff.terrr[5] = sqrt(tersoff.terrr[2]*tersoff.terrr[0]);
    tersoff.terss[3] = sqrt(tersoff.terss[0]*tersoff.terss[1]);
    tersoff.terss[4] = sqrt(tersoff.terss[1]*tersoff.terss[2]);
    tersoff.terss[5] = sqrt(tersoff.terss[2]*tersoff.terss[0]);
    tersoff.teraa[3] = sqrt(tersoff.teraa[0]*tersoff.teraa[1]);
    tersoff.teraa[4] = sqrt(tersoff.teraa[1]*tersoff.teraa[2]);
    tersoff.teraa[5] = sqrt(tersoff.teraa[2]*tersoff.teraa[0]);
    tersoff.terbb[3] = sqrt(tersoff.terbb[0]*tersoff.terbb[1]);
    tersoff.terbb[4] = sqrt(tersoff.terbb[1]*tersoff.terbb[2]);
    tersoff.terbb[5] = sqrt(tersoff.terbb[2]*tersoff.terbb[0]);
    tersoff.terlambda[3] = (tersoff.terlambda[0]+tersoff.terlambda[1])/2.0;
    tersoff.terlambda[4] = (tersoff.terlambda[1]+tersoff.terlambda[2])/2.0;
    tersoff.terlambda[5] = (tersoff.terlambda[2]+tersoff.terlambda[0])/2.0;
    tersoff.termu[3] = (tersoff.termu[0]+tersoff.termu[1])/2.0;
    tersoff.termu[4] = (tersoff.termu[1]+tersoff.termu[2])/2.0;
    tersoff.termu[5] = (tersoff.termu[2]+tersoff.termu[0])/2.0;
    tersoff.term[3]  = sqrt(tersoff.term[0]*tersoff.term[1]*1.0);
    tersoff.term[4]  = sqrt(tersoff.term[1]*tersoff.term[2]*1.0);
    tersoff.term[5]  = sqrt(tersoff.term[2]*tersoff.term[0]*1.0);
    tersoff.terw[3]  = 1.0;
    tersoff.terw[4]  = 1.0;
    tersoff.terw[5]  = 1.0;
    tersoff.terchi[3]= 0.9776;
    tersoff.terchi[4]= 1.00061;
    tersoff.terchi[5]= 1.0;
  } else if (strcmp(arg,"B_N_C") == 0) { // B-N-C
    printf("Tersoff multi (B/N/C) Matsunaga Jpn.J.Appl.Phys.39 L48-51\n");
    // for B
    tersoff.terrr[0]     = 1.80e-10;
    tersoff.terss[0]     = 2.10e-10;
    tersoff.teraa[0]     = 2.7702e2*ev2j;
    tersoff.terbb[0]     = 1.8349e2*ev2j;
    tersoff.terlambda[0] = 1.9922e10;
    tersoff.termu[0]     = 1.5856e10;
    tersoff.terbeta[0]   = 1.6000e-6;
    tersoff.tern[0]      = 3.9929;
    tersoff.term[0]      = 3;
    tersoff.terc2[0]     = dsq(5.2629e-1);
    tersoff.terd2[0]     = dsq(1.5870e-3);
    tersoff.terh[0]      = 0.5000;
    tersoff.terw[0]      = 1.0e0;
    tersoff.terchi[0]    = 1.0e0;
    // for N
    tersoff.terrr[1]     = 2.00e-10;
    tersoff.terss[1]     = 2.30e-10;
    tersoff.teraa[1]     = 1.1000e4*ev2j;
    tersoff.terbb[1]     = 2.1945e2*ev2j;
    tersoff.terlambda[1] = 5.7708e10;
    tersoff.termu[1]     = 2.5115e10;
    tersoff.terbeta[1]   = 1.0562e-1;
    tersoff.tern[1]      = 12.4498;
    tersoff.term[1]      = 3;
    tersoff.terc2[1]     = dsq(7.9934e4);
    tersoff.terd2[1]     = dsq(1.3432e2);
    tersoff.terh[1]      =-0.9973;
    tersoff.terw[1]      = 1.0e0;
    tersoff.terchi[1]    = 1.0e0;
    // for C
    tersoff.terrr[2]     = 1.80e-10;
    tersoff.terss[2]     = 2.10e-10;
    tersoff.teraa[2]     = 1.3936e3*ev2j;
    tersoff.terbb[2]     = 3.467e2*ev2j;
    tersoff.terlambda[2] = 3.4879e10;
    tersoff.termu[2]     = 2.2119e10;
    tersoff.terbeta[2]   = 1.5724e-7;
    tersoff.tern[2]      = 7.2751e-1;
    tersoff.term[2]      = 3;
    tersoff.terc2[2]     = dsq(3.8049e4);
    tersoff.terd2[2]     = dsq(4.384e0);
    tersoff.terh[2]      =-5.7058e-1;
    tersoff.terw[2]      = 1.0e0;
    tersoff.terchi[2]    = 1.0e0;
    // Dissimilar atoms
    tersoff.terrr[3] = sqrt(tersoff.terrr[0]*tersoff.terrr[1]);
    tersoff.terrr[4] = sqrt(tersoff.terrr[1]*tersoff.terrr[2]);
    tersoff.terrr[5] = sqrt(tersoff.terrr[2]*tersoff.terrr[0]);
    tersoff.terss[3] = sqrt(tersoff.terss[0]*tersoff.terss[1]);
    tersoff.terss[4] = sqrt(tersoff.terss[1]*tersoff.terss[2]);
    tersoff.terss[5] = sqrt(tersoff.terss[2]*tersoff.terss[0]);
    tersoff.teraa[3] = sqrt(tersoff.teraa[0]*tersoff.teraa[1]);
    tersoff.teraa[4] = sqrt(tersoff.teraa[1]*tersoff.teraa[2]);
    tersoff.teraa[5] = sqrt(tersoff.teraa[2]*tersoff.teraa[0]);
    tersoff.terbb[3] = sqrt(tersoff.terbb[0]*tersoff.terbb[1]);
    tersoff.terbb[4] = sqrt(tersoff.terbb[1]*tersoff.terbb[2]);
    tersoff.terbb[5] = sqrt(tersoff.terbb[2]*tersoff.terbb[0]);
    tersoff.terlambda[3] = (tersoff.terlambda[0]+tersoff.terlambda[1])/2.0;
    tersoff.terlambda[4] = (tersoff.terlambda[1]+tersoff.terlambda[2])/2.0;
    tersoff.terlambda[5] = (tersoff.terlambda[2]+tersoff.terlambda[0])/2.0;
    tersoff.termu[3] = (tersoff.termu[0]+tersoff.termu[1])/2.0;
    tersoff.termu[4] = (tersoff.termu[1]+tersoff.termu[2])/2.0;
    tersoff.termu[5] = (tersoff.termu[2]+tersoff.termu[0])/2.0;
    tersoff.term[3]  = sqrt(tersoff.term[0]*tersoff.term[1]*1.0);
    tersoff.term[4]  = sqrt(tersoff.term[1]*tersoff.term[2]*1.0);
    tersoff.term[5]  = sqrt(tersoff.term[2]*tersoff.term[0]*1.0);
    tersoff.terw[3]  = 1.1593;
    tersoff.terw[4]  = 0.9685;
    tersoff.terw[5]  = 1.0025;
    tersoff.terchi[3]= 1.0;
    tersoff.terchi[4]= 0.6381;
    tersoff.terchi[5]= 1.0;
  } else if (strcmp(arg,"SiC_Erhart") == 0) { // Si-C Erhart-Albe
    printf("Tersoff multi (Si/C) Erhart-Albe PRB71-035211\n");
    // for Si
    {
    double D0[3],r0[3],S[3],beta[3],gamma[3];
    D0[0]=3.24;D0[1]=6.00;D0[2]=4.36;
    r0[0]=2.232;r0[1]=1.4276;r0[2]=1.79;
    S[0]=1.842;S[1]=2.167;S[2]=1.847;
    beta[0]=1.4761;beta[1]=2.0099;beta[2]=1.6991;
    gamma[0]=0.114354;gamma[1]=0.11233;gamma[2]=0.011877;
    //    gamma[0]=1;gamma[1]=1;gamma[2]=1;
    tersoff.terrr[0]     = 2.68e-10;
    tersoff.terss[0]     = 2.96e-10;
    tersoff.teraa[0]     =      D0[0]/(S[0]-1)*exp(r0[0]*beta[0]*sqrt(2*S[0])) *ev2j;
    tersoff.terbb[0]     = S[0]*D0[0]/(S[0]-1)*exp(r0[0]*beta[0]*sqrt(2/S[0])) *ev2j;
    tersoff.terlambda[0] = beta[0]*sqrt(2*S[0])*1e10;
    tersoff.termu[0]     = beta[0]*sqrt(2/S[0])*1e10;
    tersoff.terbeta[0]   = 1.0;
    tersoff.tern[0]      = 1.0;
    tersoff.term[0]      = 3;
    tersoff.terc2[0]     = dsq(2.00494);
    tersoff.terd2[0]     = dsq(0.81472);
    tersoff.terh[0]      =-0.259;
    tersoff.terw[0]      = gamma[0];
    tersoff.terchi[0]    = 1.0;
    // for C
    tersoff.terrr[1]     = 1.85e-10;
    tersoff.terss[1]     = 2.15e-10;
    tersoff.teraa[1]     =      D0[1]/(S[1]-1)*exp(r0[1]*beta[1]*sqrt(2*S[1])) *ev2j;
    tersoff.terbb[1]     = S[1]*D0[1]/(S[1]-1)*exp(r0[1]*beta[1]*sqrt(2/S[1])) *ev2j;
    tersoff.terlambda[1] = beta[1]*sqrt(2*S[1])*1e10;
    tersoff.termu[1]     = beta[1]*sqrt(2/S[1])*1e10;
    tersoff.terbeta[1]   = 1.0;
    tersoff.tern[1]      = 1.0;
    tersoff.term[1]      = 3;
    tersoff.terc2[1]     = dsq(181.910);
    tersoff.terd2[1]     = dsq(6.28433);
    tersoff.terh[1]      =-0.5556;
    tersoff.terw[1]      = gamma[1];
    tersoff.terchi[1]    = 1.0;
    // Si-C
    tersoff.terrr[2]     = 2.20;
    tersoff.terss[2]     = 2.60;
    tersoff.teraa[2]     =      D0[2]/(S[2]-1)*exp(r0[2]*beta[2]*sqrt(2*S[2])) *ev2j;
    tersoff.terbb[2]     = S[2]*D0[2]/(S[2]-1)*exp(r0[2]*beta[2]*sqrt(2/S[2])) *ev2j;
    tersoff.terlambda[2] = beta[2]*sqrt(2*S[2])*1e10;
    tersoff.termu[2]     = beta[2]*sqrt(2/S[2])*1e10;
    tersoff.term[2]      = 3;
    tersoff.terc2[2]     = dsq(273987.0);
    tersoff.terd2[2]     = dsq(180.314);
    tersoff.terh[2]      =-0.68;
    tersoff.terw[2]      = gamma[2];
    tersoff.terchi[2]    = 1.0;
    }
    //  } else if (strcmp(arg,"Si/N/B") == 0) { // C-Si-Ge
    //    printf("Si/N/B Matsunaga\n");
  } else {
    printf("Tersoff arg = %s is not supported\n",arg); return;
  }
  //tersoff.terc2=tersoff.terc*tersoff.terc;
  //tersoff.terd2=tersoff.terd*tersoff.terd;
  for (int i=0; i<tersoff.nptype; i++) {
    tersoff.termum[i] = pow(tersoff.termu[i],tersoff.term[i]); } //???
}

double fcut(double r, double rr, double ss)
{
  double pi = M_PI, f;
  if (r < rr) {
    f=1.0;
  } else if (r <= ss) {
    f=0.5+0.5*cos(pi*(r-rr)/(ss-rr));
  } else {
    f=0.0;
  }
  return f;
}
double fcutd(double r, double rr, double ss)
{
  double pi = M_PI, fd;
  if ((r < ss)&&(r > rr)) {
    fd=-pi/2.0*sin(pi*(r-rr)/(ss-rr))/(ss-rr);
  } else {
    fd=0.0;
  }
  return fd;
}
