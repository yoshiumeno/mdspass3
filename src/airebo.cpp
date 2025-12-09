#include <string.h>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include "myheader.h"
#define NCC 10000

void radic(int ki, int kj, double xnt1, double xnt2, double conjug,
	   double &rad, double &drdl, double &drdm, double &drdn);
void tor(double xnt1, double xnt2, double conjug, double &ator,
	 double &drdl, double &drdm, double &drdn);
void bcuint(int kl, int ki, double xx1, double xx2, int nh, int nc,
	    double &ansy, double &ansy1, double &ansy2);

double ss(double t); double ssp(double t);
double dss(double t);double dssp(double t);
double tr(double r, int ki, int kj);
double tc(double r, int ki, int kj);
double tb(double r, int ki, int kj);
double vlj(double r, int ki, int kj);
double dvlj(double r, int ki, int kj);
void get_first_arg(std::string &line, std::string &arg1);
void airebo_lj(int mode);
void remove_carriage_return(std::string& line);

void airebo_lj()
{
  airebo_lj(0);
}

void airebo_lj(int mode)
{
// mode=0: normal, 1: partial
  double pi = 4.0*atan(1.0);
  int ig1 = 0;
  double xk[251][4],xl[251][4],xsik[251],xsjk[251],xsil[251],xsjl[251]
    ,cj[4],ck[4],cl[4],rk[4],rl[4],dt2dik[4],dt2djl[4],dt2dij[4]
    ,dexni[3],dexnj[3],xni[3],xnj[3]
    ,cfuni[NCC+1],cfunj[NCC+1],dcfuni[NCC+1],dcfunj[NCC+1]
    ,cosk[251],cosl[251],sink[251],sinl[251]
    ,dctjk[251],dctij[251],dctik[251],dctil[251],dctji[251],dctjl[251];
  // find number hydrogens and carbons connected to each atom
  for (int i=1; i<=bre.np; i++) {
    int jbegin = bre.nabors[i];
    int jend = bre.nabors[i+1]-1;
    bre.xhc1[i] = 1; bre.xhc2[i] = 1;
    if (jbegin > jend) { continue; }
    for (int j=jbegin; j<=jend; j++) {
      if (bre.lcheck[j]!=1) { continue; } // <<---- wichtig!
      int jn = bre.list[j];
      if (bre.ktype[jn]==1) { bre.xhc1[i] = bre.xhc1[i]+bre.bww[j]; } //C
      if (bre.ktype[jn]==2) { bre.xhc2[i] = bre.xhc2[i]+bre.bww[j]; } //H
    } }
  // sum over bonds between atoms i and j
  for (int i=1; i<=bre.np; i++) {
    int jbegin = bre.nabors[i];
    int jend = bre.nabors[i+1]-1;
    if (jbegin > jend) { continue; }
    int ki = bre.ktype[i];
    for (int j=jbegin; j<=jend; j++) {
      //if (bre.lcheck[j]!=1) { continue; }
      int jn = bre.list[j];
      if (i > jn) { continue; } // changed from if (i >= jn) { continue; } 2015.12.17
      ///
      if (mode > 0) {
	if (atom.lock[i]&&atom.lock[jn]) {
	  //printf("%d - %d skipped\n",i,jn);
	  continue;
	}
      }
      ///
      cj[1] = bre.cor1[j]; cj[2] = bre.cor2[j]; cj[3] = bre.cor3[j];
      double rsqij = cj[1]*cj[1] + cj[2]*cj[2] + cj[3]*cj[3];
      //double sij = bre.rcor[j];
      if (rsqij<0.0001) { continue; }
      //double rsqij = sij*sij;
      double sij = sqrt(rsqij);
      int kj = bre.ktype[jn];
      int kikj = ki+kj;
      double rijm = bre.rb1[ki][kj]; //AIREBO rij_min
      double rijm2 = rijm*rijm;
      double rij=sij;
      sij = rijm; rsqij = rijm2;
      double elj = 0.0 ; double delj = 0.0;
      double wij = ssp(tc(rij,ki,kj));
      double wmax = wij;
      int lbegin = bre.nabors[jn];
      int lend = bre.nabors[jn+1]-1;
      double cij = 1.0;
      if(rij> bre.rb2[ki][kj]*3){
      elj = vlj(rij,ki,kj); delj = dvlj(rij,ki,kj);
      }
      else{
      if(wmax<1){
      if (jbegin != jend ) {
	for (int k=jbegin; k<=jend; k++) {
	  if (k == j) { continue; }
	  if (bre.lcheck[k] != 1) { continue; }
	  int kn = bre.list[k];
	  int kk = bre.ktype[kn];
	  double rik = bre.rcor[k];
	  if(rik>bre.rb2[ki][kk]){continue;}
	  double rkjx = bre.cor1[k]-cj[1];
	  double rkjy = bre.cor2[k]-cj[2];
	  double rkjz = bre.cor3[k]-cj[3];
	  double rkj = sqrt(rkjx*rkjx+rkjy*rkjy+rkjz*rkjz);
	  double ww = ssp(tc(rik,ki,kk))*ssp(tc(rkj,kj,kk));
	  if (ww > wmax) { wmax = ww; }
	  if (wmax==1) { break; }}}}
	  
	  if(wmax<1){
	  if ((jbegin != jend)&&(lbegin != lend)) {
	for (int k=jbegin; k<=jend; k++) {
	  if (k == j) { continue; }
	  if (bre.lcheck[k] != 1) { continue; }
	  int kn = bre.list[k];
	  int kk = bre.ktype[kn];
	  ck[1] = bre.cor1[k]; ck[2] = bre.cor2[k]; ck[3] = bre.cor3[k];
	  double rik = bre.rcor[k];
	  if(rik>bre.rb2[ki][kk]){continue;}
	  for (int l=lbegin; l<=lend; l++) {
	    if (bre.lcheck[l] != 1) { continue; }
	    int ln = bre.list[l];
	    if (ln == i) { continue; }
	    int kl = bre.ktype[ln];
	    cl[1] = bre.cor1[l]; cl[2] = bre.cor2[l]; cl[3] = bre.cor3[l];
	    double rjl = bre.rcor[l];
	    if(rjl>bre.rb2[kj][kl]){continue;}
	    double rklx = cj[1]+cl[1]-ck[1];
	    double rkly = cj[2]+cl[2]-ck[2];
	    double rklz = cj[3]+cl[3]-ck[3];
	    double rkl = sqrt(rklx*rklx+rkly*rkly+rklz*rklz);
	    double ww = ssp(tc(rik,ki,kk))*ssp(tc(rjl,kj,kl))*ssp(tc(rkl,kk,kl));
	    if (ww > wmax) { wmax = ww; }
	    if (wmax==1) { break; }
	  } // loop l
	if (wmax==1) { break; }
	} // loop k
      } // if ((jbegin != jend)&&(lbegin != lend))
      }//if (wmax<1)
      cij = 1 - wmax;
      double vv = vlj(rij,ki,kj); double dvv = dvlj(rij,ki,kj);
      if(cij>0){
      if(rij>1.122462048*bre.siglj[ki][kj]){
      elj = cij*vv;
      delj = cij*dvv;
     }else{
     
      // i side of bond
      int nk = 0;
      double xsij = 0.0;
      double ssumk = 0.0;
      double conk = 0.0;
      xni[1] = bre.xhc1[i]; xni[2] = bre.xhc2[i];
      //xni[kj] = xni[kj]-bre.bww[j]; //AIREBO:modify?
      if (kj==ki) { xni[kj] = xni[kj]-wij; }
      double qi = xni[1]+xni[2]-2.0;
      double sdalik = 0.0;
      if (jbegin != jend ) {
	for (int k=jbegin; k<=jend; k++) {
	  double ali = 0.0;
	  double dali = 0.0;
	  double daldik = 0.0;
	  if (k == j) { continue; }
	  if (bre.lcheck[k] != 1) { continue; }
	  int kn = bre.list[k];
	  int kk = bre.ktype[kn];
	  nk++;
	  double s3 = bre.rcor[k];
	  double costh = (bre.cor1[k]*cj[1]+bre.cor2[k]*cj[2]+bre.cor3[k]*cj[3])/(s3*rij);
	  cosk[nk] = costh;
	  sink[nk] = sqrt(1.0-costh*costh);
	  if (acos(costh) > pi) { sink[nk] = -sink[nk]; }
	  int ig = bre.igc[int(-costh*12)+13];
	  double gangle;
	  if (ki == 1) {
	    if (ig != 4) {
	      gangle = bre.spgc[1][ig]+bre.spgc[2][ig]*costh;
	      for (int jj=3; jj<=6; jj++) {
		gangle = gangle+bre.spgc[jj][ig]*pow(costh,jj-1); }
	    } else {
	      ali = 0.0; dali = 0.0;
	      if (qi < bre.xqm) {
		ali = 1.0;
		if (qi > bre.att) {
		  double dtemp = bre.pq*(qi-bre.att);
		  ali = (1.0+cos(dtemp))/2.0;
		  dali = -bre.pq/2.0*sin(dtemp); } }
	      gangle = bre.spgc[1][ig]+bre.spgc[2][ig]*costh;
	      ig1=ig+1;
	      double gangle1 = bre.spgc[1][ig1]+bre.spgc[2][ig1]*costh;
	      for (int jj=3; jj<=6; jj++) {
		gangle = gangle+bre.spgc[jj][ig]*pow(costh,jj-1);
		gangle1 = gangle1+bre.spgc[jj][ig1]*pow(costh,jj-1); }
	      daldik = dali*(gangle1-gangle);
	      gangle = gangle+ali*(gangle1-gangle);
	    } // if (ig != 4)
	  } else { // if (ki == 1)
	    ig = bre.igh[int(-costh*12.0)+13];
	    gangle = bre.spgh[1][ig]+bre.spgh[2][ig]*costh;
	    for (int jj=3; jj<=6; jj++) {
	      gangle = gangle+bre.spgh[jj][ig]*pow(costh,jj-1); }
	  } // if (ki == 1)
	  //double fc = bre.bww[k]; double dfc = bre.dww[k]; //AIREBO:modify?
	  double fc = ssp(tc(s3,ki,kk));
	  double dfc = dssp(tc(s3,ki,kk))/(bre.rb2[ki][kk]-bre.rb1[ki][kk]);
	  cfuni[nk] = 0.0; dcfuni[nk] = 0.0;
	  if (kk == 1) {
	    double xx = bre.xhc1[kn]+bre.xhc2[kn]-fc-2.0;
	    if (xx <= 3.0) {
	      if (xx <= 2.0) {
		cfuni[nk] = 1.0;
	      } else {
		double px = pi*(xx-2.0);
		cfuni[nk] = (1.0+cos(px))/2.0;
		dcfuni[nk] = -fc*sin(px)*pi/2.0;
	      } } }
	  conk = conk+fc*cfuni[nk];
	  
	  double exx;
	  if (bre.xdb[ki][kj][kk] != 0.0) {
	    exx = bre.reg[ki][kj][kk]*exp(bre.xdb[ki][kj][kk]*(sij-s3));
	  } else {
	    exx = 1.0;
	  }
	  
	  double gs = gangle*exx;
	  ssumk = ssumk+fc*gs;

	} // loop k
      } // if (jbegin != jend )
      
      // j side of bond
      int nl = 0;
      double xsji = 0.0;
      double ssuml = 0.0;
      double conl = 0.0;
            
      xnj[1] = bre.xhc1[jn];
      xnj[2] = bre.xhc2[jn];
      //xnj[ki] = xnj[ki] - bre.bww[j]; //AIREBO:modify?
      if (kj==ki) { xnj[kj] = xnj[kj]-wij; }
      double qj = xnj[1]+xnj[2]-2.0;
      double sdaljl = 0.0;
      if (lbegin != lend) {
	for (int l=lbegin; l<=lend; l++) {
	  double alj = 0.0;
	  double dalj = 0.0;
	  double daldjl = 0.0;
	  int ln = bre.list[l];
	  if (ln == i) { continue; }
	  if (bre.lcheck[l] != 1) { continue; }
	  int kl = bre.ktype[ln];
	  nl++;
	  double s3 = bre.rcor[l];
	  double costh = -(bre.cor1[l]*cj[1]+bre.cor2[l]*cj[2]+bre.cor3[l]*cj[3])/(s3*rij);
	  cosl[nl] = costh;
	  sinl[nl] = sqrt(1.0-costh*costh);
	  if (acos(costh) > pi) { sinl[nl] = -sinl[nl]; }
	  double gangle;
	  if (kj == 1) {
	    int ig = bre.igc[int(-costh*12)+13];
	    if (ig != 4) {
	      gangle = bre.spgc[1][ig]+bre.spgc[2][ig]*costh;
	      for (int jj=3; jj<=6; jj++) {
		gangle = gangle+bre.spgc[jj][ig]*pow(costh,jj-1); }
	    } else {
	      alj = 0.0; dalj = 0.0;
	      if (qj < bre.xqm) {
		alj = 1.0;
		if (qj > bre.att) {
		  double dtemp = bre.pq*(qj-bre.att);
		  alj = (1.0+cos(dtemp))/2.0;
		  dalj = -bre.pq/2.0*sin(dtemp); } }
	      gangle = bre.spgc[1][ig]+bre.spgc[2][ig]*costh;
	      ig1=ig+1;
	      double gangle1 = bre.spgc[1][ig1]+bre.spgc[2][ig1]*costh;
	      for (int jj=3; jj<=6; jj++) {
		gangle = gangle+bre.spgc[jj][ig]*pow(costh,jj-1);
		gangle1 = gangle1+bre.spgc[jj][ig1]*pow(costh,jj-1); }
	      daldjl = dalj*(gangle1-gangle);
	      gangle = gangle+alj*(gangle1-gangle);
	    } // if (ig != 4) (2)
	  } else { // if (ki == 1) (2)
	    int ig = bre.igh[int(-costh*12.0)+13];
	    gangle = bre.spgh[1][ig]+bre.spgh[2][ig]*costh;
	    for (int jj=3; jj<=6; jj++) {
	      gangle = gangle+bre.spgh[jj][ig]*pow(costh,jj-1); }
	  } // if (ki == 1) (2)
	  //double fc = bre.bww[l]; double dfc = bre.dww[l]; //AIREBO:modify?
	  double fc = ssp(tc(s3,kl,kj));
	  double dfc = dssp(tc(s3,kl,kj))/(bre.rb2[kl][kj]-bre.rb1[kl][kj]);
	  cfunj[nl] = 0.0; dcfunj[nl] = 0.0;
	  if (kl == 1) {
	    double xx = bre.xhc1[ln]+bre.xhc2[ln]-fc-2.0;
	    if (xx <= 3.0) {
	      if (xx <= 2.0) {
		cfunj[nl] = 1.0;
	      } else {
		double px = pi*(xx-2.0);
		cfunj[nl] = (1.0+cos(px))/2.0;
		dcfunj[nl] = -fc*sin(px)*pi/2.0;
	      } } }
	  conl = conl+fc*cfunj[nl];

	  double exx;
	  if (bre.xdb[kj][ki][kl] != 0.0) {
	    exx = bre.reg[kj][ki][kl]*exp(bre.xdb[kj][ki][kl]*(sij-s3));
	  } else {
	    exx = 1.0;
	  }

	  double gs = gangle*exx;
	  ssuml = ssuml+fc*gs;

	} // loop l
      } // if (lbegin != lend)


      double exnij = 0.0;
      dexni[1] = 0.0; dexni[2] = 0.0;
      if (ki == 1) {
	int nh = int(xni[2]+1.0e-12);
	int nc = int(xni[1]+1.0e-12);
	if ((fabs((float)nh-xni[2]) > 1.0e-8)||(fabs((float)nc-xni[1])) > 1.0e-8) {
	  bcuint(ki,kj,xni[2],xni[1],nh,nc,exnij,dexni[2],dexni[1]);
	} else {
	  exnij = bre.xh[kj][nh][nc];
	  dexni[2] = bre.xh1[kj][nh][nc];
	  dexni[1] = bre.xh2[kj][nh][nc];
	}
      }
      
      double exnji = 0.0;
      dexnj[1] = 0.0; dexnj[2] = 0.0;
      if (kj == 1) {
	int nh = int(xnj[2]+1.0e-12);
	int nc = int(xnj[1]+1.0e-12);
	if ((fabs((float)nh-xnj[2]) > 1.0e-8)||(fabs((float)nc-xnj[1])) > 1.0e-8) {
	  bcuint(kj,ki,xnj[2],xnj[1],nh,nc,exnji,dexnj[2],dexnj[1]);
	} else {
	  exnji = bre.xh[ki][nh][nc];
	  dexnj[2] = bre.xh1[ki][nh][nc];
	  dexnj[1] = bre.xh2[ki][nh][nc];
	}
      }

      double dij, bij, dji, bji, dbdzi, dbdzj, vatt, dradi, dradj, drdc, conjug,
	xnt1, xnt2, rad, btot = 0;
      	dij = (1.0+exnij+ssumk);
	bij = pow(dij,-0.5);
	dji = (1.0+exnji+ssuml);
	bji = pow(dji,-0.5);
	dbdzi = -0.5*bij/dij;
	dbdzj = -0.5*bji/dji;
	vatt = bre.exx1[j];
	dradi = 0.0;
	dradj = 0.0;
	drdc = 0.0;
	conjug = 1.0+conk*conk+conl*conl;
	xnt1 = xni[1]+xni[2]-1.0;
	xnt2 = xnj[1]+xnj[2]-1.0;
	rad = 0.0; // ???
	radic(ki,kj,xnt1,xnt2,conjug,rad,dradi,dradj,drdc);
	btot = bji+bij+rad;
      

      // dihedral terms
      if (kikj == bre.ndihed) {
	//dbtori = 0.0; dbtorj = 0.0; dbtorc = 0.0;
	double datori = 0.0; double datorj = 0.0; double datorc = 0.0; // ???
	double btor = 0.0;
	double ator = 0.0; // ???
	tor(xnt1,xnt2,conjug,ator,datori,datorj,datorc);
	if ((jbegin != jend)&&(lbegin != lend)) {
	  nk = 0;
	  for (int k=jbegin; k<=jend; k++) {
	    if (k == j) { continue; } // go to 220
	    if (bre.lcheck[k] != 1) { continue; } // go to 220
	    nk++;
	    if (fabs(sink[nk]) < 1.0e-1) { continue; } // go to 220
	    int kn = bre.list[k];
	    ck[1] = bre.cor1[k]; ck[2] = bre.cor2[k]; ck[3] = bre.cor3[k];
	    double rck = bre.rcor[k];
	    double fck;
	    if (bre.ktype[kn] == 2) {
	      fck = 1.0;
	      if (rck >= 1.60) { continue; } // go to 220
	      if (rck >= 1.30) {
		double dtemp = bre.pidt*(rck-1.30);
		fck = (1.0+cos(dtemp))/2.0; }
	    } else {
	      fck = bre.bww[k];
	    }
	    nl = 0;
	    for (int l=lbegin; l<=lend; l++) {
	      int ln = bre.list[l];
	      if (ln == i) { continue; } // go to 210
	      if (bre.lcheck[l] != 1) { continue; } // go to 210
	      nl++;
	      if (fabs(sinl[nl]) < 1.0e-1) { continue; } // go to 210
	      double sinl2 = sinl[nl]*sinl[nl];
	      cl[1] = bre.cor1[l]; cl[2] = bre.cor2[l]; cl[3] = bre.cor3[l];
	      double rcl = bre.rcor[l];
	      double fcl;
	      if (bre.ktype[ln] == 2) {
		fcl = 1.0;
		if (rcl >= 1.60) { continue; } // go to 210
		if (rcl >= 1.30) {
		  double dtemp = bre.pidt*(rcl-1.30);
		  fcl = (1.0+cos(dtemp))/2.0;
		}
	      } else {
		fcl = bre.bww[l];
	      }
	      double crkx=ck[2]*cj[3]-cj[2]*ck[3];
	      double crlx=cj[2]*cl[3]-cl[2]*cj[3];
	      double crky=ck[3]*cj[1]-cj[3]*ck[1];
	      double crly=cj[3]*cl[1]-cl[3]*cj[1];
	      double crkz=ck[1]*cj[2]-cj[1]*ck[2];
	      double crlz=cj[1]*cl[2]-cl[1]*cj[2];
	      double tik=crkx*crkx+crky*crky+crkz*crkz;
	      double tjl=crlx*crlx+crly*crly+crlz*crlz;
	      double t1=sqrt(tik*tjl);
	      double t2=crkx*crlx+crky*crly+crkz*crlz;
	      double cw=-t2/t1;
	      double bt=(1.0-cw*cw);
	      btor=btor+bt*fck*fcl;

	      //printf("HOGE1\n");
	    } // loop l (210)
	  } // loop k (220)
	} // if ((jbegin != jend)&&(lbegin != lend)) (230)
	btot=btot+btor*ator;//if(i==1){printf("%f %f\n",btor,ator);}
      } // if (kikj == ndihed) (231)
      // end dihedral forces
      // dh-loop 2 for wmax

    btot/=2.0;
    if(btot<bre.bmin[ki][kj]){
      elj = cij*vv;
      delj = cij*dvv;
      }else{
      elj = (ss(tr(rij,ki,kj))*ss(tb(btot,ki,kj))+1-ss(tr(rij,ki,kj)))
	*cij*vv;
      delj = (ss(tr(rij,ki,kj))*ss(tb(btot,ki,kj))+1-ss(tr(rij,ki,kj)))
	*cij*dvv
      	+dss(tr(rij,ki,kj))/(0.122462048*bre.siglj[ki][kj])
	*(ss(tb(btot,ki,kj))-1)*cij*vv;
	
	}}}}
      if(i==jn){elj/=2;delj/=2;} // added 2015.12.17
      atom.epot[i]  += elj/2;
      atom.epot[jn] += elj/2;
      double rp; double repx,repy,repz;
      rp = -delj/rij;
      repx = rp*cj[1];
      bre.rnp1[i]=bre.rnp1[i]+repx; bre.rnp1[jn]=bre.rnp1[jn]-repx;
      repy = rp*cj[2];
      bre.rnp2[i]=bre.rnp2[i]+repy; bre.rnp2[jn]=bre.rnp2[jn]-repy;
      repz = rp*cj[3];
      bre.rnp3[i]=bre.rnp3[i]+repz; bre.rnp3[jn]=bre.rnp3[jn]-repz;
      // stress
      repx /= 2; repy /= 2; repz /= 2;
      atom.satom[i][0][0] -= repx * cj[1];
      atom.satom[i][0][1] -= repx * cj[2];
      atom.satom[i][0][2] -= repx * cj[3];
      atom.satom[i][1][1] -= repy * cj[2];
      atom.satom[i][1][2] -= repy * cj[3];
      atom.satom[i][2][2] -= repz * cj[3];
      atom.satom[jn][0][0] -= repx * cj[1];
      atom.satom[jn][0][1] -= repx * cj[2];
      atom.satom[jn][0][2] -= repx * cj[3];
      atom.satom[jn][1][1] -= repy * cj[2];
      atom.satom[jn][1][2] -= repy * cj[3];
      atom.satom[jn][2][2] -= repz * cj[3];
      
      
  	//E torsion
	if(wij>0){double dwij = dssp(tc(rij,ki,kj))/(bre.rb2[ki][kj]-bre.rb1[ki][kj]);
	if ((jbegin != jend)&&(lbegin != lend)) {
	  
	  for (int k=jbegin; k<=jend; k++) {
	    if (k == j) { continue; }
	    if (bre.lcheck[k] != 1) { continue; } 
	    int kn = bre.list[k];
	    int kk = bre.ktype[kn];
	    ck[1] = bre.cor1[k]; ck[2] = bre.cor2[k]; ck[3] = bre.cor3[k];
	    double rck = sqrt(ck[1]*ck[1]+ck[2]*ck[2]+ck[3]*ck[3]);
	    
	    double fck = ssp(tc(rck,ki,kk));
	    double dfck = dssp(tc(rck,ki,kk))/(bre.rb2[ki][kk]-bre.rb1[ki][kk]);
	     if(fck==0){continue;}
	   
	    for (int l=lbegin; l<=lend; l++) {
	      int ln = bre.list[l];
	      if (ln == i) { continue; }
	      if (bre.lcheck[l] != 1) { continue; } 
	      int kl = bre.ktype[ln];
	      cl[1] = bre.cor1[l]; cl[2] = bre.cor2[l]; cl[3] = bre.cor3[l];
	      double rcl = sqrt(cl[1]*cl[1]+cl[2]*cl[2]+cl[3]*cl[3]);
	      double fcl= ssp(tc(rcl,kj,kl));
	      double dfcl = dssp(tc(rcl,kj,kl))/(bre.rb2[kj][kl]-bre.rb1[kj][kl]);
	       if(fcl==0){continue;}
	     
	      double crkx=ck[2]*cj[3]-cj[2]*ck[3];
	      double crlx=cj[2]*cl[3]-cl[2]*cj[3];
	      double crky=ck[3]*cj[1]-cj[3]*ck[1];
	      double crly=cj[3]*cl[1]-cl[3]*cj[1];
	      double crkz=ck[1]*cj[2]-cj[1]*ck[2];
	      double crlz=cj[1]*cl[2]-cl[1]*cj[2];
	      double tik=crkx*crkx+crky*crky+crkz*crkz;
	      double tjl=crlx*crlx+crly*crly+crlz*crlz;
	      if(tik*tjl==0){continue;}
	      double t1=sqrt(tik*tjl);
	      double t2=crkx*crlx+crky*crly+crkz*crlz;
	      double cw=-t2/t1;
	      double cw2=(cw+1.00000)/2.00000;
	      double bt=256*cw2*cw2*cw2*cw2*cw2/405-0.10;
	      double dbt=256*cw2*cw2*cw2*cw2/81;
	      double rijik=cj[1]*ck[1]+cj[2]*ck[2]+cj[3]*ck[3];
	      double rijjl=cj[1]*cl[1]+cj[2]*cl[2]+cj[3]*cl[3];
	      double rikjl=cl[1]*ck[1]+cl[2]*ck[2]+cl[3]*ck[3];
	      double dcwdij[4],dcwdik[4],dcwdjl[4],detordij[4],detordik[4],detordjl[4];
	      double etor=bre.epsts[kk][kl]*fck*fcl*wij*bt;
	      for(int n=1;n<4;n++){
	      dcwdij[n]=cw*((ck[n]*rijjl+cl[n]*rijik-2*cj[n]*rikjl)/t2-(cj[n]*rck*rck-ck[n]*rijik)/tik-(cj[n]*rcl*rcl-cl[n]*rijjl)/tjl);
	      dcwdik[n]=cw*((cj[n]*rijjl-cl[n]*rij*rij)/t2-(ck[n]*rij*rij-cj[n]*rijik)/tik);
	      dcwdjl[n]=cw*((cj[n]*rijik-ck[n]*rij*rij)/t2-(cl[n]*rij*rij-cj[n]*rijjl)/tjl);
	      detordij[n]=bre.epsts[kk][kl]*fck*fcl*(dwij*cj[n]/rij*bt+wij*dbt*0.5*dcwdij[n]);
	      detordik[n]=bre.epsts[kk][kl]*wij*fcl*(dfck*ck[n]/rck*bt+fck*dbt*0.5*dcwdik[n]);
	      detordjl[n]=bre.epsts[kk][kl]*wij*fck*(dfcl*cl[n]/rcl*bt+fcl*dbt*0.5*dcwdjl[n]);
	       }
	      if(i==jn){elj/=2;delj/=2;} // added 2015.12.17
	     atom.epot[i]  += etor/2;
	     atom.epot[jn] += etor/2;
	     
	     	     
	     repx=-detordij[1];repy=-detordij[2];repz=-detordij[3];
	     bre.rnp1[i]=bre.rnp1[i]+repx; bre.rnp1[jn]=bre.rnp1[jn]-repx;
            bre.rnp2[i]=bre.rnp2[i]+repy; bre.rnp2[jn]=bre.rnp2[jn]-repy;
            bre.rnp3[i]=bre.rnp3[i]+repz; bre.rnp3[jn]=bre.rnp3[jn]-repz;
            
            repx /= 2; repy /= 2; repz /= 2;
            atom.satom[i][0][0] -= repx * cj[1];
            atom.satom[i][0][1] -= repx * cj[2];
            atom.satom[i][0][2] -= repx * cj[3];
            atom.satom[i][1][1] -= repy * cj[2];
            atom.satom[i][1][2] -= repy * cj[3];
            atom.satom[i][2][2] -= repz * cj[3];
            atom.satom[jn][0][0] -= repx * cj[1];
            atom.satom[jn][0][1] -= repx * cj[2];
            atom.satom[jn][0][2] -= repx * cj[3];
            atom.satom[jn][1][1] -= repy * cj[2];
            atom.satom[jn][1][2] -= repy * cj[3];
            atom.satom[jn][2][2] -= repz * cj[3];
            
            repx=-detordik[1];repy=-detordik[2];repz=-detordik[3];
	     bre.rnp1[i]=bre.rnp1[i]+repx; bre.rnp1[kn]=bre.rnp1[kn]-repx;
            bre.rnp2[i]=bre.rnp2[i]+repy; bre.rnp2[kn]=bre.rnp2[kn]-repy;
            bre.rnp3[i]=bre.rnp3[i]+repz; bre.rnp3[kn]=bre.rnp3[kn]-repz;
            
            repx /= 2; repy /= 2; repz /= 2;
            atom.satom[i][0][0] -= repx * ck[1];
            atom.satom[i][0][1] -= repx * ck[2];
            atom.satom[i][0][2] -= repx * ck[3];
            atom.satom[i][1][1] -= repy * ck[2];
            atom.satom[i][1][2] -= repy * ck[3];
            atom.satom[i][2][2] -= repz * ck[3];
            atom.satom[kn][0][0] -= repx * ck[1];
            atom.satom[kn][0][1] -= repx * ck[2];
            atom.satom[kn][0][2] -= repx * ck[3];
            atom.satom[kn][1][1] -= repy * ck[2];
            atom.satom[kn][1][2] -= repy * ck[3];
            atom.satom[kn][2][2] -= repz * ck[3];
            
            repx=-detordjl[1];repy=-detordjl[2];repz=-detordjl[3];
	     bre.rnp1[jn]=bre.rnp1[jn]+repx; bre.rnp1[ln]=bre.rnp1[ln]-repx;
            bre.rnp2[jn]=bre.rnp2[jn]+repy; bre.rnp2[ln]=bre.rnp2[ln]-repy;
            bre.rnp3[jn]=bre.rnp3[jn]+repz; bre.rnp3[ln]=bre.rnp3[ln]-repz;
            
            repx /= 2; repy /= 2; repz /= 2;
            atom.satom[jn][0][0] -= repx * cl[1];
            atom.satom[jn][0][1] -= repx * cl[2];
            atom.satom[jn][0][2] -= repx * cl[3];
            atom.satom[jn][1][1] -= repy * cl[2];
            atom.satom[jn][1][2] -= repy * cl[3];
            atom.satom[jn][2][2] -= repz * cl[3];
            atom.satom[ln][0][0] -= repx * cl[1];
            atom.satom[ln][0][1] -= repx * cl[2];
            atom.satom[ln][0][2] -= repx * cl[3];
            atom.satom[ln][1][1] -= repy * cl[2];
            atom.satom[ln][1][2] -= repy * cl[3];
            atom.satom[ln][2][2] -= repz * cl[3];
              } // loop l (210)
	  } // loop k (220)
	} // if ((jbegin != jend)&&(lbegin != lend)) (230)
	}//end E torsion
	
    } // loop j
  } // loop i
} // end of pibond()

double vlj(double r, int ki, int kj)
{
  double x, y, y3, y6, y12;
  y = bre.siglj[ki][kj]/r;
  y3 = y * y * y; y6 = y3 * y3; y12 = y6 * y6;
  x = 4 * bre.epslj[ki][kj] * (y12 - y6);
  return x;
}
double dvlj(double r, int ki, int kj)
{
  double x, y, y3, y6, y12;
  y = bre.siglj[ki][kj]/r;
  y3 = y * y * y; y6 = y3 * y3; y12 = y6 * y6;
  x = 4 * bre.epslj[ki][kj] * (-12 * y12 + 6 * y6) / r;
  return x;
}
double ss(double t)
{
  double x;
  if (t < 0) { x = 1; }
  else if (t > 1) { x = 0; }
  else { x = 1.0-t*t*(3-2*t); }
  return x;
}
double dss(double t)
{
  double x;
  if (t < 0) { x = 0; }
  else if (t > 1) { x = 0; }
  else { x = 6 * t * (t - 1); }
  return x;
}
double ssp(double t)
{
  double x;
  if (t < 0) { x = 1; }
  else if (t > 1) { x = 0; }
  else { x = (1 + cos(M_PI * t)) / 2; }
  return x;
}
double dssp(double t)
{
  double x;
  if (t < 0) { x = 0; }
  else if (t > 1) { x = 0; }
  else { x = -M_PI * sin(M_PI * t) / 2; }
  return x;
}
double tr(double r, int ki, int kj)
{
  double x = (r-bre.siglj[ki][kj])/(0.122462048*bre.siglj[ki][kj]);
  return x;
}
double tc(double r, int ki, int kj)
{
  double x = (r-bre.rb1[ki][kj])/(bre.rb2[ki][kj]-bre.rb1[ki][kj]);
  return x;
}
double tb(double r, int ki, int kj)
{
  double x = (r-bre.bmin[ki][kj])/(bre.bmax[ki][kj]-bre.bmin[ki][kj]);
  return x;
}

void airebo_lj_read()
{
  int mesh_lj = 0;
  FILE *fp;
  char fname[80] = "aaa";
  std::string line, arg1, arg2, arg3;
  char larg1[40], larg2[40], larg3[40]; double xx, yy; int ix, iy;
  bool ifmeshdata = false;
  std::ifstream fin(fname); fin.close(); fin.clear(); // "clear" needed for Win
#if defined __linux__ || defined __APPLE__
  strcpy(fname,"pot/AIREBO.lj");
#else
  strcpy(fname,"pot\\AIREBO.lj");
#endif
  // If mesh data exists, it is used.
  fin.open(fname);
  if (fin) {
    printf("%s exists\n", fname);
    ifmeshdata = true; ix = -1;
    while (getline(fin,line)) { remove_carriage_return(line);
      if (ix == -1) {
	get_first_arg(line,arg1); get_first_arg(line, arg2);
	strcpy (larg1, arg1.c_str()); strcpy (larg2, arg2.c_str());
	bre.ljmin = atof(larg1); bre.ljmax = atof(larg2);
      } else {
	get_first_arg(line,arg1); get_first_arg(line, arg2); get_first_arg(line, arg3);
	strcpy (larg1, arg1.c_str()); strcpy (larg2, arg2.c_str()); strcpy (larg3, arg3.c_str());
	bre.lj[ix] = atof(larg2); bre.ljd[ix] = atof(larg3);
      }
      ix++; }
    bre.ljmesh = ix;
  }
  fin.close(); fin.clear();
  //If fine-mesh data do not exist, they are created.
  if (!ifmeshdata) {
    printf("%s does not exists\n", fname);
    std::ofstream fout;
    fout.open(fname);
    fout.precision(10);
    fout << bre.ljmin << " " << bre.ljmax << std::endl;
    for (int i=0; i<bre.ljmesh; i++) {
      double xx = bre.ljmin + (bre.ljmax-bre.ljmin)/(double)bre.ljmesh * (double)i;
      double yy = vlj(xx, 1, 1);
      double zz = dvlj(xx, 1, 1);
      bre.lj[i] = yy; bre.ljd[i] = zz;
      fout << xx << " " << yy << " " << zz << std::endl;
    }
    fout.close(); fout.clear();
  }

  
}
