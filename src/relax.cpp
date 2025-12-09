#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "myheader.h"
#define DIS_LIMIT 0.1

void vscale();
void vscale(double target);
void vscale(int i, double target);
void intpl_cg (double x1, double y1, double dy1,
	       double x2, double y2, double dy2, double &x3);
void relax_fixatom_set();
void relax_damper(float x);

extern float fp_alph_ini, fp_ffinc, fp_ffdec, fp_ffalph, fp_ffdtmax;
extern int fp_nfmin;
extern int istep;
extern int relax_accel, relax_accel_interval;
extern float relax_accel_threshold;

void relax_damper(float x)
{
  for (int i=1; i<=atom.natom; i++) {
    atom.vx[i] /= (double)x + 1.0;
    atom.vy[i] /= (double)x + 1.0;
    atom.vz[i] /= (double)x + 1.0;
  }
}

void relax_fixatom_set()
{
  if (istep % relax_accel_interval == 0) {
    for (int i=1; i<=atom.natom; i++) {
      atom.lock[i] = false;
    }
  } else if (istep % relax_accel_interval == 1) {
    double fmx = atom.Fmax(); fmx *= relax_accel_threshold; fmx = fmx * fmx;
    int i0=0; int i1=0;
    for (int i=1; i<=atom.natom; i++) {
      double xx = atom.fx[i]*atom.fx[i] + atom.fy[i]*atom.fy[i] + atom.fz[i]*atom.fz[i];
      if (xx < fmx) {
	atom.lock[i] = true;
	i0++;
      } else {
	atom.lock[i] = false;
	i1++;
      }
    }
    //printf("accel %f %\n",(float)i0/(float)(i0+i1)*100);
  }
}

void relax_gloc()
{
  if (relax_accel == 1) {
    relax_fixatom_set();
  }
    double dp = 0.0;
    if (atom.QC==0) {
      for (int i=1; i<=atom.natom; i++) {
	//dp = dp + atom.vx[i]*atom.fx[i] + atom.vy[i]*atom.fy[i] + atom.vz[i]*atom.fz[i];
	if (!atom.mfx[i]) { dp += atom.vx[i]*atom.fx[i]; }
	if (!atom.mfy[i]) { dp += atom.vy[i]*atom.fy[i]; }
	if (!atom.mfz[i]) { dp += atom.vz[i]*atom.fz[i]; }
      }
      if (dp < 0.0) { vscale(0.0); }
    } else {
      for (int i=1; i<=atom.natom; i++) {
	if (atom.repatom[i]) {
	  //dp = dp + atom.vx[i]*atom.fx[i] + atom.vy[i]*atom.fy[i] + atom.vz[i]*atom.fz[i];
	  if (!atom.mfx[i]) { dp += atom.vx[i]*atom.fx[i]; }
	  if (!atom.mfy[i]) { dp += atom.vy[i]*atom.fy[i]; }
	  if (!atom.mfz[i]) { dp += atom.vz[i]*atom.fz[i]; }
	}
	// if (dp < 0.0) { vscale(i, 0.0); }
      }
      if (dp < 0.0) { vscale(0.0); }
    }
}

void relax_fire_reset()
{
  fire.sp = 0.0; fire.fnorm = 0.0; fire.vnorm = 0.0;
  fire.fire_alph = fp_alph_ini;
  fire.fire_alph_ini = fp_alph_ini;
  fire.nfmin = fp_nfmin;
  fire.ffinc = fp_ffinc; fire.ffdec = fp_ffdec; fire.ffalph = fp_ffalph;
  fire.ifire = 0;
  fire.ffdtmax = fp_ffdtmax*1.0e-15;
}

void relax_fire()
{
  /*
    Atom relaxation algorithm by Fast Inertial Relaxation Engine (FIRE)
    E.Bitzek, P.Koshinen, F.Gahler, M.Moseler and P.Gumbsch
    Physical Review Letters, 97, 170201 (2006) 
  */
  // Parameters for FIRE
  // Initial values are written in "myclass.h"
  //if (istep == 0) { fire.ffdtmax=dt*10.0; } // ROUGHLY ESTIMATED VALUE 
  //if (fire.ffdtmax < 1.0e-20) { fire.ffdtmax=dt*10.0; }

  //-----F1 : P=F*v 
  fire.ifire++;
  fire.sp=0.0;
  fire.fnorm=0.0;
  fire.vnorm=0.0;
  for (int i=1; i<=atom.natom; i++) {
    fire.sp = fire.sp + atom.fx[i]*atom.vx[i]+atom.fy[i]*atom.vy[i]+atom.fz[i]*atom.vz[i];
    fire.fnorm += atom.fx[i]*atom.fx[i]+atom.fy[i]*atom.fy[i]+atom.fz[i]*atom.fz[i];
    fire.vnorm += atom.vx[i]*atom.vx[i]+atom.vy[i]*atom.vy[i]+atom.vz[i]*atom.vz[i];
  }
  fire.fnorm=sqrt(fire.fnorm);
  fire.vnorm=sqrt(fire.vnorm);

  //-----F2 : v = (1-a)*v + a*F*|v|
  for (int i=1; i<=atom.natom; i++) {
    atom.vx[i]=(1.0-fire.fire_alph)*atom.vx[i]
      +fire.fire_alph*atom.fx[i]/fire.fnorm*fire.vnorm;
    atom.vy[i]=(1.0-fire.fire_alph)*atom.vy[i]
      +fire.fire_alph*atom.fy[i]/fire.fnorm*fire.vnorm;
    atom.vz[i]=(1.0-fire.fire_alph)*atom.vz[i]
      +fire.fire_alph*atom.fz[i]/fire.fnorm*fire.vnorm;
  }

  //-----F3 : P > 0
  if ((fire.sp > 0.0)&&(fire.ifire >= fire.nfmin)) {
    //dt = min(dt*ffinc,ffdtmax);
    dt = dt*fire.ffinc; if (dt > fire.ffdtmax) { dt = fire.ffdtmax; }
    fire.fire_alph = fire.fire_alph * fire.ffalph;
    //-----F4 : P <= 0
  }  else if (fire.sp <= 0.0) {
    fire.ifire = 0;
    dt = dt*fire.ffdec;
    fire.fire_alph = fire.fire_alph_ini;
    vscale(0.0);
  }
  // To avoid burst...
  double xx = atom.Enkin()*2.0/3.0/(double)atom.natom/1.380662e-23;
  if (xx > 500) {
    fire.ifire = 0;
    dt = dt*fire.ffdec;
    fire.fire_alph = fire.fire_alph_ini;
    vscale(0.0);
  }
  //      print *,'FIRE:', dt, ifire, fire_alph, ffdtmax
}

void relax_cg_reset()
{
  cg.lamtr = cg.lam_init;
  cg.step = 0;
  for (int i=1; i<=atom.natom; i++) {
    cg.rstx[i]=atom.rx[i]; cg.rsty[i]=atom.ry[i]; cg.rstz[i]=atom.rz[i];
    for (int j=0; j<2; j++) {
      cg.sgx[i][j]=0; cg.sgy[i][j]=0; cg.sgz[i][j]=0;
      cg.shx[i][j]=0; cg.shy[i][j]=0; cg.shz[i][j]=0;
    }
  }
}

void relax_cg()
{
  double gamcg, gamcg1, gamcg2;
  double lamcrr;
  cg.step++;

  // Initial step (PREVIOUS step data)
  if (cg.step==1) {
    cg.pot0 = atom.epotsum;
    cg.grad0 = 0.0;
    for (int i=1; i<=atom.natom; i++) {
      cg.grad0 -= (atom.fx[i]*atom.fx[i]+atom.fy[i]*atom.fy[i]+atom.fz[i]*atom.fz[i]); }
    cg.normsh = sqrt(fabs(cg.grad0));
    cg.grad0 /= cg.normsh;
    for (int i=1; i<=atom.natom; i++) {
      cg.rstx[i] = atom.rx[i];
      cg.rsty[i] = atom.ry[i];
      cg.rstz[i] = atom.rz[i]; }
  } else if (cg.step % 2 == 0) { // Trial step (PREVIOUS step data)
    cg.pot1  = atom.epotsum;
    cg.grad1 = 0.0;
    for (int i=1; i<=atom.natom; i++) {
      cg.grad1 -= (atom.fx[i]*cg.shx[i][0]+atom.fy[i]*cg.shy[i][0]+atom.fz[i]*cg.shz[i][0]); }
    cg.grad1 /= cg.normsh;
  } else if (cg.step % 2 == 1) { // Corrector step (PREVIOUS step data)
    cg.pot0  = atom.epotsum;
    cg.grad0 = 0.0;
    for (int i=1; i<=atom.natom; i++) {
      cg.grad0 -= (atom.fx[i]*cg.shx[i][0]+atom.fy[i]*cg.shy[i][0]+atom.fz[i]*cg.shz[i][0]); }
    cg.grad0 /= cg.normsh;
    for (int i=1; i<=atom.natom; i++) {
      cg.rstx[i] = atom.rx[i];
      cg.rsty[i] = atom.ry[i];
      cg.rstz[i] = atom.rz[i]; }
  }

  //-----Search vector shxyz
  if (cg.step % 2 == 1) {  // trial step
    for (int i=1; i<=atom.natom; i++) {
      cg.sgx[i][0] = atom.fx[i];
      cg.sgy[i][0] = atom.fy[i];
      cg.sgz[i][0] = atom.fz[i];
      if (atom.mfx[i] == true) { cg.sgx[i][0] = 0.0; }
      if (atom.mfy[i] == true) { cg.sgy[i][0] = 0.0; }
      if (atom.mfz[i] == true) { cg.sgz[i][0] = 0.0; }
    }
    // set GAMMA
    if (cg.step == 1) {
      gamcg = 0.0;
    } else {
      gamcg1 = 0.0;
      gamcg2 = 0.0;
      for (int i=1; i<=atom.natom; i++) {
	gamcg1 += (cg.sgx[i][0]-cg.sgx[i][1])*cg.sgx[i][0]
	  +(cg.sgy[i][0]-cg.sgy[i][1])*cg.sgy[i][0]
	  +(cg.sgz[i][0]-cg.sgz[i][1])*cg.sgz[i][0];
	gamcg2 += cg.sgx[i][1]*cg.sgx[i][1]+cg.sgy[i][1]*cg.sgy[i][1]+cg.sgz[i][1]*cg.sgz[i][1];
      }
      gamcg = gamcg1/gamcg2;
      // Steepest Descent method : set gamcg=0
      //gamcg=0.0;
    }
    // set search vector SH
    for (int i=1; i<=atom.natom; i++) {
      cg.shx[i][0] = cg.sgx[i][0] + gamcg * cg.shx[i][1];
      cg.shy[i][0] = cg.sgy[i][0] + gamcg * cg.shy[i][1];
      cg.shz[i][0] = cg.sgz[i][0] + gamcg * cg.shz[i][1];
    }
    cg.normsh = 0.0;
    for (int i=1; i<=atom.natom; i++) {
      cg.normsh +=
	cg.shx[i][0]*cg.shx[i][0]+cg.shy[i][0]*cg.shy[i][0]+cg.shz[i][0]*cg.shz[i][0]; }
    cg.normsh = sqrt(cg.normsh);
    for (int i=1; i<=atom.natom; i++) {
      cg.sgx[i][1]=cg.sgx[i][0];
      cg.sgy[i][1]=cg.sgy[i][0];
      cg.sgz[i][1]=cg.sgz[i][0];
      cg.shx[i][1]=cg.shx[i][0];
      cg.shy[i][1]=cg.shy[i][0];
      cg.shz[i][1]=cg.shz[i][0];
    }
  }

  // Check if displacement is small
  double dis = 0; double xx;
  for (int i=1; i<=atom.natom; i++) {
    xx = cg.shx[i][0]*cg.shx[i][0]+cg.shy[i][0]*cg.shy[i][0]+cg.shz[i][0]*cg.shz[i][0];
    if (xx > dis) { dis = xx; }
  }
  dis = sqrt(dis); dis *= cg.lamtr;
  if (dis/ang > DIS_LIMIT) {
    printf("CG over limit %d  %e %e\n",istep,dis/ang,cg.lamtr);
    //for (int i=1; i<=atom.natom; i++) {
      //cg.shx[i][0] /= 2.0; cg.shy[i][0] /= 2.0; cg.shz[i][0] /= 2.0;
    //}
    //cg.lamtr = cg.lam_init;
    //relax_cg_reset();
    //return;
  }

  //-----Line minimization
  // Trial step
  if (cg.step % 2 == 1) {
    for (int i=1; i<=atom.natom; i++) {
      atom.rx[i] = cg.rstx[i]+cg.shx[i][0]*cg.lamtr;
      atom.ry[i] = cg.rsty[i]+cg.shy[i][0]*cg.lamtr;
      atom.rz[i] = cg.rstz[i]+cg.shz[i][0]*cg.lamtr;
    }
  } else if (cg.step % 2 == 0) { // Corrector step
    intpl_cg(0.0,cg.pot0,cg.grad0,cg.lamtr*cg.normsh,cg.pot1,cg.grad1,lamcrr);
    lamcrr /= cg.normsh;
    cg.lamtr = lamcrr;
    //printf("%d %e\n",cg.step,lamcrr);
    // RESET LAMDA_TRIAL
    if (cg.lamtr <= 1.0e-8) { cg.lamtr = cg.lam_init; }
    for (int i=1; i<=atom.natom; i++) {
      atom.rx[i] = cg.rstx[i]+cg.shx[i][0]*lamcrr;
      atom.ry[i] = cg.rsty[i]+cg.shy[i][0]*lamcrr;
      atom.rz[i] = cg.rstz[i]+cg.shz[i][0]*lamcrr;
    }
  }
}

//-----------------------------------------------------------------------
void intpl_cg (double x1, double y1, double dy1,
	       double x2, double y2, double dy2, double &x3)
{
  double a, b, c, d, xx, ddy, root;
      
  //     ** Only 'x1=0.0d0' is acceptable **

  //     Cubic funciton interpolation (Default)
  a = -(2.0*y2-dy2*x2-dy1*x2-2.0*y1)/(x2*x2*x2);
  b = (3.0*y2-3.0*y1-2.0*dy1*x2-dy2*x2)/(x2*x2);
  c = dy1;
  d = y1;
  root = b*b-3.0*a*c;
  if (root > 0.0) {
    xx = (-b+sqrt(root))/(3.0*a);
    ddy = 6.0*a*xx+2.0*b;
    if (ddy > 0.0) {
      x3 = xx;
    } else {
      x3 = (-b-sqrt(root))/(3.0*a);
    }
  } else { //Quadratic function (if cubic interpolation failed)
    a = (y2-y1-dy1*x2)/(x2*x2);
    b = dy1;
    x3 = -b/(2.0*a);
    if (a <= 0.0) {
      a = -(y2-y1-dy2*x2)/(x2*x2);
      b = (2.0*y2-2.0*y1-dy2*x2)/x2;
      x3 = -b/(2.0*a);
      if (a <= 0.0) { x3 = x2; }
    }
  }
}
