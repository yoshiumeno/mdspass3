#include <iostream>
#include<string.h>
#include<stdlib.h>
#include<fstream>
#include <sstream>
#include "myheader.h"

void adp_alloc();
void adp_param_tin();
void adp_param_NdFeB();

void mishin     (double rr, double **param, double &val, double &grad, int p);
void mishin_sc  (double rr, double **param, double &val, double &grad, int p);
void csw2       (double rr, double *param,  double &val, double &grad);
void csw2_sc    (double rr, double **param, double &val, double &grad, int p);
void poly_5     (double rr, double **param, double &val, double &grad, int p);
void exp_plus   (double rr, double *param,  double &val, double &grad);
void exp_plus_sc(double rr, double **param, double &val, double &grad, int p);
void remove_carriage_return(std::string& line);

// Kubo 20140124 ---------------
/*
void psi(double rr, double rc, double h, double &val, double &grad);
void eopp_exp   (double rr, double **param, double &val, double &grad, int p);
void eopp_exp_sc(double rr, double **param, double &val, double &grad, int p);
void meopp      (double rr, double **param, double &val, double &grad, int p);
void meopp_sc   (double rr, double **param, double &val, double &grad, int p);
*/
// End -------------------------

double dsquare(double d);
void resetmat(double a[3][3]);
int atom_number(char* at);
int adptyp(char* at);
int adp_read_number(const char* fname);
int adp_read_param(const char* fname);
void get_first_arg(std::string &line, std::string &arg1);
int remove_head_spaces(std::string &line);
int count_arg_number(std::string line);
int remove_after_sharp(std::string &line);
int adp_find_param(const char* fname, const char* tag, std::string &line);
int adp_find_param(std::istream& fin, const char* tag, std::string &line);

void e_force_adp()
{
  double rr, rr2, drx, dry, drz, vp0, v0;
  int j, ix, iy, iz;
  //double rcut = 8.0e0; double rcut2 = rcut * rcut;
  int    self;
  int    h, k, l, typ1, typ2, uf, us, stresses, ipair;
  double value, grad, value_tail, grad_tail, grad_i, grad_j, p_sr_tail;
  double value_el, grad_el, ggrad_el;

  double phi_val, phi_grad, u_val, u_grad, w_val, w_grad, rho_val, rho_grad;
  double eb_val, eb_grad, rho_val_j, rho_grad_j;

  // Initialization and setup
  if (adp.initialize) {
    printf("ADP initialize...\n");
    printf(" ADP potential file = %s\n",adp.fname);
    adp_read_number(adp.fname);
    printf(" # of species = %d, # of pairs = %d\n",adp.ntype, adp.nptype);
    //adp.ntype = 3; adp.nptype = 6; //NdFeB
    //adp.ntype = 1; adp.nptype = 1; //Sn

// Kubo 20140226 ---------------
    adp.Phi = NULL;
    adp.Rho = NULL;
    adp.F   = NULL;
    adp.U   = NULL;
    adp.W   = NULL;
// End -------------------------

    adp_alloc();
    //adp_param_tin();
    //adp_param_NdFeB();
    if (adp_read_param(adp.fname) == 0) {
      printf("ADP error (unsupported atom)\n"); mdmotion=0; return; }
    adp.initialize = false;
    printf("ADP initialize done.\n");
  }
  rcut = adp.cut; rcut2 = rcut * rcut;

  // First loop (reset variables)
  for (int i=1; i<=atom.natom; i++) {
    atom.fx[i] = 0.0; atom.fy[i] = 0.0; atom.fz[i] = 0.0;  atom.epot[i] = 0.0;
    adp.rho[i] = 0.0; adp.gradF[i] = 0.0;
    adp.mu_x[i] = 0.0; adp.mu_y[i] = 0.0; adp.mu_z[i] = 0.0;
    adp.lambda_xx[i] = 0.0; adp.lambda_yy[i] = 0.0; adp.lambda_zz[i] = 0.0;
    adp.lambda_xy[i] = 0.0; adp.lambda_yz[i] = 0.0; adp.lambda_zx[i] = 0.0;
    for (int j=0; j<3; j++) { for (int k=0; k<3; k++) {
	atom.satom[i][j][k] = 0.0; } }
  }// FIRST LOOP END

  // Second loop (pair term and atomic density)
  for (int i=1; i<=atom.natom; i++) {
    //typ1 = adptyp(atom.asp[i]);
    typ1 = adp.typen[i];
    if (book.alistnum[i]>0) { for (int k=1; k<=book.alistnum[i]; k++) {
	j  = book.alist[i][k][0];
	if (j<i) continue; //OK???
	ix = book.alist[i][k][1]; iy = book.alist[i][k][2]; iz = book.alist[i][k][3];
	rr2 = atom.Dist2(i,j,ix,iy,iz)/ang/ang; //rr = sqrt(rr2);
	//typ2 = adptyp(atom.asp[j]);
	typ2 = adp.typen[j];
	ipair = adp.ptype[typ1][typ2];
	drx = atom.Dx(i,j,ix,iy,iz); dry = atom.Dy(i,j,ix,iy,iz); drz = atom.Dz(i,j,ix,iy,iz);
	drx /= ang; dry /= ang; drz /= ang;
	// pair potential
	if (rr2 < rcut2) { rr = sqrt(rr2);
	  self = 0; if (i==j) { self = 1; }

// Kubo 20140225 ---------------
//	  mishin_sc(rr, adp.param_pair, phi_val, phi_grad, ipair);
    adp.Phi[ipair]->Calc(rr,phi_val,phi_grad);
// End -------------------------

	  if (self) { phi_val *= 0.5; phi_grad *= 0.5; }
	  phi_grad /= rr;
	  atom.epot[i] += phi_val/2.0; atom.epot[j] += phi_val/2.0;
	  atom.fx[i] += drx*phi_grad; atom.fy[i] += dry*phi_grad; atom.fz[i] += drz*phi_grad;
	  atom.fx[j] -= drx*phi_grad; atom.fy[j] -= dry*phi_grad; atom.fz[j] -= drz*phi_grad;
	  atom.satom[i][0][0] += drx*drx*phi_grad/2.0;
	  atom.satom[i][0][1] += dry*drx*phi_grad/2.0;
	  atom.satom[i][0][2] += drz*drx*phi_grad/2.0;
	  atom.satom[i][1][1] += dry*dry*phi_grad/2.0;
	  atom.satom[i][1][2] += drz*dry*phi_grad/2.0;
	  atom.satom[i][2][2] += drz*drz*phi_grad/2.0;
	  atom.satom[j][0][0] += drx*drx*phi_grad/2.0;
	  atom.satom[j][0][1] += dry*drx*phi_grad/2.0;
	  atom.satom[j][0][2] += drz*drx*phi_grad/2.0;
	  atom.satom[j][1][1] += dry*dry*phi_grad/2.0;
	  atom.satom[j][1][2] += drz*dry*phi_grad/2.0;
	  atom.satom[j][2][2] += drz*drz*phi_grad/2.0;

	  // dipole distortion

// Kubo 20140225 ---------------
//	  exp_plus_sc(rr, adp.param_dp, u_val, u_grad, ipair);
    adp.U[ipair]->Calc(rr,u_val,u_grad);
// End -------------------------

	  if (self) { u_val *= 0.5; u_grad *= 0.5; }
	  //u_val /= rr; ###
	  adp.mu_x[i] += u_val*drx; adp.mu_x[j] -= u_val*drx;
	  adp.mu_y[i] += u_val*dry; adp.mu_y[j] -= u_val*dry;
	  adp.mu_z[i] += u_val*drz; adp.mu_z[j] -= u_val*drz;

	  // quadrupole distortion

// Kubo 20140225 ---------------
//	  exp_plus_sc(rr, adp.param_qp, w_val, w_grad, ipair);
    adp.W[ipair]->Calc(rr,w_val,w_grad);
// End -------------------------

	  if (self) { w_val *= 0.5; w_grad *= 0.5; }
	  //w_val /= rr2; ###
	  adp.lambda_xx[i] += w_val*drx*drx; adp.lambda_xx[j] += w_val*drx*drx;
	  adp.lambda_yy[i] += w_val*dry*dry; adp.lambda_yy[j] += w_val*dry*dry;
	  adp.lambda_zz[i] += w_val*drz*drz; adp.lambda_zz[j] += w_val*drz*drz;
	  adp.lambda_xy[i] += w_val*drx*dry; adp.lambda_xy[j] += w_val*drx*dry;
	  adp.lambda_yz[i] += w_val*dry*drz; adp.lambda_yz[j] += w_val*dry*drz;
	  adp.lambda_zx[i] += w_val*drz*drx; adp.lambda_zx[j] += w_val*drz*drx;
	  // atomic density (???)
	  if (typ1 == typ2) {

// Kubo 20140225 ---------------
//	    csw2_sc(rr, adp.param_den, rho_val, rho_grad, typ2);
      adp.Rho[typ2]->Calc(rr,rho_val,rho_grad);
// End -------------------------

	    if (self) { rho_val *= 0.5; rho_grad *= 0.5; }
	    adp.rho[i] += rho_val; adp.rho[j] += rho_val;
	  } else {

// Kubo 20140225 ---------------
//	    csw2_sc(rr, adp.param_den, rho_val, rho_grad, typ2);
//	    csw2_sc(rr, adp.param_den, rho_val_j, rho_grad_j, typ1);
      adp.Rho[typ2]->Calc(rr,rho_val,rho_grad);
      adp.Rho[typ1]->Calc(rr,rho_val_j,rho_grad_j);
// End -------------------------

	    adp.rho[i] += rho_val; adp.rho[j] += rho_val_j;
	  }
	} //endif rr2<rcut2
      } } } // Second loop end
  for (int i=1; i<=atom.natom; i++) { // Second loop supplement
    typ1 = adp.typen[i];

    // embedding energy

// Kubo 20140225 ---------------
//    poly_5(adp.rho[i], adp.param_eb, eb_val, eb_grad, typ1);
    adp.F[typ1]->Calc(adp.rho[i],eb_val,eb_grad);
//    std::cout<<"(^Q^)~~ "<< eb_val <<std::endl; exit(0);
// End -------------------------

    atom.epot[i] += eb_val;
    adp.gradF[i] += eb_grad;
    // ADP energy
    double tmp = 0.0;
    tmp += dsquare(adp.mu_x[i]);
    tmp += dsquare(adp.mu_y[i]);
    tmp += dsquare(adp.mu_z[i]);
    adp.nu[i] = adp.lambda_xx[i] + adp.lambda_yy[i] + adp.lambda_zz[i];
    double tr = adp.nu[i] / 3.0;
    tmp += dsquare(adp.lambda_xx[i] - tr);
    tmp += dsquare(adp.lambda_yy[i] - tr);
    tmp += dsquare(adp.lambda_zz[i] - tr);
    tmp += dsquare(adp.lambda_xy[i])*2.0;
    tmp += dsquare(adp.lambda_yz[i])*2.0;
    tmp += dsquare(adp.lambda_zx[i])*2.0;
    tmp /= 2.0;
    atom.epot[i] += tmp;
  } // Second loop supplement end

//  for(int i=1;i<=atom.natom;i++){std::cout<<adp.rho[i]<<std::endl;} exit(0); // Kubo

  // Third loop (APD force)
  for (int i=1; i<=atom.natom; i++) {
    //typ1 = adptyp(atom.asp[i]);
    typ1 = adp.typen[i];
    if (book.alistnum[i]>0) { for (int k=1; k<=book.alistnum[i]; k++) {
	j  = book.alist[i][k][0];
	if (j<i) continue; //OK???
	ix = book.alist[i][k][1]; iy = book.alist[i][k][2]; iz = book.alist[i][k][3];
	rr2 = atom.Dist2(i,j,ix,iy,iz)/ang/ang; //rr = sqrt(rr2);
	//typ2 = adptyp(atom.asp[j]);
	typ2 = adp.typen[j];
	ipair = adp.ptype[typ1][typ2];
	drx = atom.Dx(i,j,ix,iy,iz); dry = atom.Dy(i,j,ix,iy,iz); drz = atom.Dz(i,j,ix,iy,iz);
	drx /= ang; dry /= ang; drz /= ang;
	if (rr2 < rcut2) { rr = sqrt(rr2);
	  self = 0; if (i==j) { self = 1; }

	  // rho (???)

// Kubo 20140227 ---------------
//	  csw2_sc(rr, adp.param_den, rho_val, rho_grad, typ2);
    adp.Rho[typ2]->Calc(rr,rho_val,rho_grad);
// End -------------------------

	  if (typ1 == typ2) {
	    rho_grad_j = rho_grad;
	  } else {
// Kubo 20140227 ---------------
//	    csw2_sc(rr, adp.param_den, rho_val, rho_grad_j, typ1);
    adp.Rho[typ1]->Calc(rr,rho_val,rho_grad_j);
// End -------------------------
	  }
	  double tmp =
	    rho_grad * adp.gradF[i] + rho_grad_j * adp.gradF[j];
	  if (self) { tmp /= 2.0; }
	  tmp /= rr;
	  atom.fx[i] += drx*tmp; atom.fx[j] -= drx*tmp;
	  atom.fy[i] += dry*tmp; atom.fy[j] -= dry*tmp;
	  atom.fz[i] += drz*tmp; atom.fz[j] -= drz*tmp;
	  atom.satom[i][0][0] += drx*drx*tmp/2.0;
	  atom.satom[i][0][1] += dry*drx*tmp/2.0;
	  atom.satom[i][0][2] += drz*drx*tmp/2.0;
	  atom.satom[i][1][1] += dry*dry*tmp/2.0;
	  atom.satom[i][1][2] += drz*dry*tmp/2.0;
	  atom.satom[i][2][2] += drz*drz*tmp/2.0;
	  atom.satom[j][0][0] += drx*drx*tmp/2.0;
	  atom.satom[j][0][1] += dry*drx*tmp/2.0;
	  atom.satom[j][0][2] += drz*drx*tmp/2.0;
	  atom.satom[j][1][1] += dry*dry*tmp/2.0;
	  atom.satom[j][1][2] += drz*dry*tmp/2.0;
	  atom.satom[j][2][2] += drz*drz*tmp/2.0;
	  // dipole
	  double utmpx = adp.mu_x[i] - adp.mu_x[j];
	  double utmpy = adp.mu_y[i] - adp.mu_y[j];
	  double utmpz = adp.mu_z[i] - adp.mu_z[j];
	  //if (self) { utmpx /= 2.0; utmpy /= 2.0; utmpz /= 2.0; } // Bug in potfit? 
	  // (for small cell it does harm)

// Kubo 20140225 ---------------
//	  exp_plus_sc(rr, adp.param_dp, u_val, u_grad, ipair); // calculated again..
    adp.U[ipair]->Calc(rr,u_val,u_grad);
// End -------------------------

	  if (self) { u_val *= 0.5; u_grad *= 0.5; }
	  double utmp = (utmpx*drx + utmpy*dry + utmpz*drz) * u_grad;
	  utmp /= rr;
	  double tmpx = utmpx * u_val + utmp * drx;
	  double tmpy = utmpy * u_val + utmp * dry;
	  double tmpz = utmpz * u_val + utmp * drz;
	  atom.fx[i] += tmpx; atom.fx[j] -= tmpx;
	  atom.fy[i] += tmpy; atom.fy[j] -= tmpy;
	  atom.fz[i] += tmpz; atom.fz[j] -= tmpz;
	  atom.satom[i][0][0] += drx*tmpx/2.0;
	  atom.satom[i][0][1] += dry*tmpx/2.0;
	  atom.satom[i][0][2] += drz*tmpx/2.0;
	  atom.satom[i][1][1] += dry*tmpy/2.0;
	  atom.satom[i][1][2] += drz*tmpy/2.0;
	  atom.satom[i][2][2] += drz*tmpz/2.0;
	  atom.satom[j][0][0] += drx*tmpx/2.0;
	  atom.satom[j][0][1] += dry*tmpx/2.0;
	  atom.satom[j][0][2] += drz*tmpx/2.0;
	  atom.satom[j][1][1] += dry*tmpy/2.0;
	  atom.satom[j][1][2] += drz*tmpy/2.0;
	  atom.satom[j][2][2] += drz*tmpz/2.0;
	  // quadrupole
	  double wxx = adp.lambda_xx[i] + adp.lambda_xx[j];
	  double wyy = adp.lambda_yy[i] + adp.lambda_yy[j];
	  double wzz = adp.lambda_zz[i] + adp.lambda_zz[j];
	  double wxy = adp.lambda_xy[i] + adp.lambda_xy[j];
	  double wyz = adp.lambda_yz[i] + adp.lambda_yz[j];
	  double wzx = adp.lambda_zx[i] + adp.lambda_zx[j];
	  //if (self) { // Bug in potfit? (for small cell it does harm)
	  //  wxx /= 2.0; wyy /= 2.0; wzz /= 2.0;
	  //  wxy /= 2.0; wyz /= 2.0; wzx /= 2.0; }
	  double vx = (wxx*drx + wxy*dry + wzx*drz);
	  double vy = (wxy*drx + wyy*dry + wyz*drz);
	  double vz = (wzx*drx + wyz*dry + wzz*drz);
	  double nu = (adp.nu[i] + adp.nu[j])/3.0;

// Kubo 20140225 ---------------
//	  exp_plus_sc(rr, adp.param_qp, w_val, w_grad, ipair); // calculated again..
    adp.W[ipair]->Calc(rr,w_val,w_grad);
// End -------------------------

	  if (self) { w_val /= 2.0; w_grad /= 2.0; }
	  double f1 = w_val * 2.0;
	  double f2 = ((vx*drx+vy*dry+vz*drz)-nu*rr2)*w_grad-nu*f1*rr;
	  tmpx = f1 * vx + f2 * drx / rr;
	  tmpy = f1 * vy + f2 * dry / rr;
	  tmpz = f1 * vz + f2 * drz / rr;
	  atom.fx[i] += tmpx; atom.fx[j] -= tmpx;
	  atom.fy[i] += tmpy; atom.fy[j] -= tmpy;
	  atom.fz[i] += tmpz; atom.fz[j] -= tmpz;
	  atom.satom[i][0][0] += drx*tmpx/2.0;
	  atom.satom[i][0][1] += dry*tmpx/2.0;
	  atom.satom[i][0][2] += drz*tmpx/2.0;
	  atom.satom[i][1][1] += dry*tmpy/2.0;
	  atom.satom[i][1][2] += drz*tmpy/2.0;
	  atom.satom[i][2][2] += drz*tmpz/2.0;
	  atom.satom[j][0][0] += drx*tmpx/2.0;
	  atom.satom[j][0][1] += dry*tmpx/2.0;
	  atom.satom[j][0][2] += drz*tmpx/2.0;
	  atom.satom[j][1][1] += dry*tmpy/2.0;
	  atom.satom[j][1][2] += drz*tmpy/2.0;
	  atom.satom[j][2][2] += drz*tmpz/2.0;
	  

	} // rr2<rcut2
      } } } // Third loop end

 OUT:
  atom.epotsum=0.0;
  for (int i=1; i<=atom.natom; i++) {
    atom.epot[i] *= eV;
    atom.epotsum += atom.epot[i];
    atom.fx[i] *= eV/ang; atom.fy[i] *= eV/ang; atom.fz[i] *= eV/ang;
    //printf("%10d %20.12e %20.12e %20.12e\n",i-1,atom.fx[i]/eV*ang,atom.fy[i]/eV*ang,atom.fz[i]/eV*ang);
  }
    
  //Atomic stress
  for (int i=1; i<=atom.natom; i++) {
    atom.satom[i][1][0] = atom.satom[i][0][1];
    atom.satom[i][2][0] = atom.satom[i][0][2];
    atom.satom[i][2][1] = atom.satom[i][1][2]; }
  resetmat(cell.dmat);
  for (int i=1; i<=atom.natom; i++) {
    for (int j=0; j<3; j++) {
      for (int k=0; k<3; k++) {
	atom.satom[i][j][k] *= eV;
	cell.dmat[j][k] -= atom.satom[i][j][k]; } } }
  cell.virx = cell.dmat[0][0];
  cell.viry = cell.dmat[1][1];
  cell.virz = cell.dmat[2][2];
  cell.volume = cell.Getvolume();
  for (int i=1; i<=atom.natom; i++) {
    for (int j=0; j<3; j++) {
      for (int k=0; k<3; k++) {
	atom.satom[i][j][k] *= (double)atom.natom / cell.volume; } } }

}


// Potential Functions =========================================

//void mishin(double rr, double *param, double &val, double &grad)
void mishin(double rr, double **param, double &val, double &grad, int ipair)
{
  double z  = rr - param[ipair][3];
  double e  = exp(-param[ipair][5] * z);
  double p  = pow(z,param[ipair][4]);
  double pp = param[ipair][4] * p / z;
  double ep = -param[ipair][5] * e;
  val  = param[ipair][0] * p * e * (1.0 + param[ipair][1] * e) + param[ipair][2];
  grad = param[ipair][0] * ( pp  *  e * (1.0 + param[ipair][1] *  e)
			     + p * ep * (1.0 + param[ipair][1] *  e)
			     + p *  e * (      param[ipair][1] * ep) );
}
void mishin_sc(double rr, double **param, double &val, double &grad, int ipair)
{
  double x = (rr - adp.cut)/param[ipair][6];
  double x2 = x*x; double x3 = x2*x; double x4 = x2*x2;
  double sc = x4 / (1.0 + x4); double scp = 4.0*x3/(1.0+x4)/(1.0+x4)/param[ipair][6];
  double z  = rr - param[ipair][3];
  //double e  = exp(-param[ipair][5] * z); //<-- should be this
  double e  = exp(-param[ipair][5] * rr); //<-- definition in potfit
  double p  = pow(z,param[ipair][4]);
  double pp = param[ipair][4] * p / z;
  double ep = -param[ipair][5] * e;
  double tmp  = param[ipair][0] * p * e * (1.0 + param[ipair][1] * e) + param[ipair][2];
  val  = tmp * sc;
  grad = param[ipair][0] * ( pp  *  e * (1.0 + param[ipair][1] *  e)
			     + p * ep * (1.0 + param[ipair][1] *  e)
			     + p *  e * (      param[ipair][1] * ep) ) * sc + tmp * scp;
}
void csw2(double rr, double *param, double &val, double &grad)
{
  double p  = pow(rr, param[3]);
  val  = ( 1.0 + param[0] * cos(param[1] * rr + param[2]) ) / p;
  grad = ( -param[0] * param[1] * sin(param[1] * rr + param[2]) ) / p
    - param[3] * val / rr;
}
void csw2_sc(double rr, double **param, double &val, double &grad, int ityp)
{
  double x = (rr - adp.cut)/param[ityp][4];
  double x2 = x*x; double x3 = x2*x; double x4 = x2*x2;
  double sc = x4 / (1.0 + x4); double scp = 4.0*x3/(1.0+x4)/(1.0+x4)/param[ityp][4];
  double p  = pow(rr, param[ityp][3]);
  double val0 = ( 1.0 + param[ityp][0] * cos(param[ityp][1] * rr + param[ityp][2]) ) / p;
  //val  = ( 1.0 + param[ityp][0] * cos(param[ityp][1] * rr + param[ityp][2]) ) / p * sc;
  val = val0 * sc;
  grad = (( -param[ityp][0] * param[ityp][1] * sin(param[ityp][1] * rr + param[ityp][2]) ) / p
	  - param[ityp][3] * val0 / rr) * sc + val0 * scp;
}
void poly_5(double rr, double **param, double &val, double &grad, int type)
{
  double  r1 = rr - 1.0;
  double dr1 = r1 * r1;
  val  = param[type][0] + 0.5 * param[type][1] * dr1
    + param[type][2]*r1*dr1 + param[type][3]*dr1*dr1 + param[type][4]*r1*dr1*dr1;
  grad = param[type][1] * r1
    + 3.0*param[type][2]*dr1 + 4.0*param[type][3]*r1*dr1 + 5.0*param[type][4]*dr1*dr1;
}
void exp_plus(double rr, double *param, double &val, double &grad)
{
  double tmp = param[0] * exp(-param[1]*rr);
  val  = tmp + param[2];
  grad = -param[1] * tmp;
}
void exp_plus_sc(double rr, double **param, double &val, double &grad, int ipair)
{
  double x = (rr - adp.cut)/param[ipair][3];
  double x2 = x*x; double x3 = x2*x; double x4 = x2*x2;
  double sc = x4 / (1.0 + x4); double scp = 4.0*x3/(1.0+x4)/(1.0+x4)/param[ipair][3];
  double tmp = param[ipair][0] * exp(-param[ipair][1]*rr);
  val  = (tmp + param[ipair][2]) * sc;
  grad = (-param[ipair][1] * tmp) * sc + (tmp + param[ipair][2]) * scp;
}


// Potential Functions Above ===================================







void adp_param_tin()
{
  adp.cut = 7.0;
  adp.type[1] = 50;
  adp.ptype[1][1] = 1;
  for (int i=1; i<=atom.natom; i++) {
    adp.typen[i] = adptyp(atom.asp[i]);
  }
  adp.param_pair[1][0] = 0.00691025;
  adp.param_pair[1][1] = 29.92791121;
  adp.param_pair[1][2] = -0.19277842;
  adp.param_pair[1][3] = 0.00001165;
  adp.param_pair[1][4] = 0.98767802;
  adp.param_pair[1][5] = 0.15430821;
  adp.param_pair[1][6] = 0.5;
  adp.param_den[1][0] = -0.75198545;
  adp.param_den[1][1] = 0.69706211;
  adp.param_den[1][2] = -2.93829476;
  adp.param_den[1][3] = 1.36634379;
  adp.param_den[1][4] = 2.0;
  adp.param_eb[1][0] = -4.41231862;
  adp.param_eb[1][1] = 6.16520171;
  adp.param_eb[1][2] = -4.09634234;
  adp.param_eb[1][3] = -2.60663587;
  adp.param_eb[1][4] = 6.89527313;
  adp.param_dp[1][0] = 102.41749562;
  adp.param_dp[1][1] = 2.30668629;
  adp.param_dp[1][2] = -0.01064225;
  adp.param_dp[1][3] = 1.53956319;
  adp.param_qp[1][0] = -103.46596987;
  adp.param_qp[1][1] = 2.56364479;
  adp.param_qp[1][2] = 0.00128067;
  adp.param_qp[1][3] = 0.76080274;
}

void adp_param_NdFeB()
{
  adp.cut = 7.0;
  adp.type[1] = 60;
  adp.type[2] = 26;
  adp.type[3] =  5;
  adp.ptype[1][1] = 1;
  adp.ptype[1][2] = 2;
  adp.ptype[1][3] = 3;
  adp.ptype[2][2] = 4;
  adp.ptype[2][3] = 5;
  adp.ptype[3][3] = 6;
  adp.ptype[2][1] = adp.ptype[1][2];
  adp.ptype[3][1] = adp.ptype[1][3];
  adp.ptype[3][2] = adp.ptype[2][3];
  for (int i=1; i<=atom.natom; i++) {
    adp.typen[i] = adptyp(atom.asp[i]);
  }
  // ADP for Nd-Fe-B
  // Development_Code: NdFeB_001-35b
  int i;
  // Phi ---------------------------------------------------------
  // Nd-Nd
  i = 1;
  adp.param_pair[i][0] =      221.10457017;
  adp.param_pair[i][1] =    16016.40327088;
  adp.param_pair[i][2] =       -0.01618936;
  adp.param_pair[i][3] =        0.41538704;
  adp.param_pair[i][4] =        7.41680115;
  adp.param_pair[i][5] =        4.05315715;
  adp.param_pair[i][6] =        0.50001857;
  // Nd-Fe
  i = 2;
  adp.param_pair[i][0] =    18020.32373548;
  adp.param_pair[i][1] =  3999954.54384480;
  adp.param_pair[i][2] =       -0.02212760;
  adp.param_pair[i][3] =        0.11236013;
  adp.param_pair[i][4] =       12.70425999;
  adp.param_pair[i][5] =        8.53320160;
  adp.param_pair[i][6] =        1.99995711;
  // Nd-B
  i = 3;
  adp.param_pair[i][0] =        0.01287612;
  adp.param_pair[i][1] =  2997699.14805589;
  adp.param_pair[i][2] =        0.09412598;
  adp.param_pair[i][3] =        0.00005888;
  adp.param_pair[i][4] =        8.74290757;
  adp.param_pair[i][5] =        3.93289203;
  adp.param_pair[i][6] =        1.58859291;
  // Fe-Fe
  i = 4;
  adp.param_pair[i][0] =        0.10730396;
  adp.param_pair[i][1] =      322.69235837;
  adp.param_pair[i][2] =       -1.26453922;
  adp.param_pair[i][3] =        0.95067629;
  adp.param_pair[i][4] =        6.60868963;
  adp.param_pair[i][5] =        1.42792364;
  adp.param_pair[i][6] =        1.99999904;
  // Fe-B
  i = 5;
  adp.param_pair[i][0] =        0.04211440;
  adp.param_pair[i][1] = 78453654.41735350;
  adp.param_pair[i][2] =       -0.30711725;
  adp.param_pair[i][3] =        0.00001649;
  adp.param_pair[i][4] =       10.81116364;
  adp.param_pair[i][5] =        6.32414689;
  adp.param_pair[i][6] =        1.99999951;
  // B-B
  i = 6;
  adp.param_pair[i][0] =        0.97097443;
  adp.param_pair[i][1] =   298600.58981664;
  adp.param_pair[i][2] =       -0.66438945;
  adp.param_pair[i][3] =        0.55623000;
  adp.param_pair[i][4] =       14.70598170;
  adp.param_pair[i][5] =        4.71612278;
  adp.param_pair[i][6] =        1.99993115;
  // rho ---------------------------------------------------------
  // Nd
  i = 1;
  adp.param_den[i][0] =  -0.24731054;
  adp.param_den[i][1] =   1.82649822;
  adp.param_den[i][2] =   3.14129199;
  adp.param_den[i][3] =   2.82482455;
  adp.param_den[i][4] =   1.99973870;
  // Fe
  i = 2;
  adp.param_den[i][0] =  -0.16358490;
  adp.param_den[i][1] =   3.41038264;
  adp.param_den[i][2] =  -1.29258358;
  adp.param_den[i][3] =   2.73216487;
  adp.param_den[i][4] =   1.99999706;
  // B
  i = 3;
  adp.param_den[i][0] =  -0.38585879;
  adp.param_den[i][1] =   3.13645134;
  adp.param_den[i][2] =   0.09229026;
  adp.param_den[i][3] =   2.85385686;
  adp.param_den[i][4] =   1.99944137;
  // F -----------------------------------------------------------
  // Nd
  i = 1;
  adp.param_eb[i][0] =  -12.05769815;
  adp.param_eb[i][1] =   18.61562351;
  adp.param_eb[i][2] =   -7.47690514;
  adp.param_eb[i][3] =   -7.77562852;
  adp.param_eb[i][4] =    0.00280314;
  // Fe
  i = 2;
  adp.param_eb[i][0] =   -2.58650536;
  adp.param_eb[i][1] =   11.93719370;
  adp.param_eb[i][2] =    1.95157877;
  adp.param_eb[i][3] =   -8.22475289;
  adp.param_eb[i][4] =    4.43693905;
  // B
  i = 3;
  adp.param_eb[i][0] =   -4.36894265;
  adp.param_eb[i][1] =    2.57853337;
  adp.param_eb[i][2] =   -0.59532339;
  adp.param_eb[i][3] =    0.31918245;
  adp.param_eb[i][4] =    0.00048675;
  // u -----------------------------------------------------------
  // Nd-Nd
  i = 1;
  adp.param_dp[i][0] = -1.52620830;
  adp.param_dp[i][1] =  0.01460996;
  adp.param_dp[i][2] =  1.43732218;
  adp.param_dp[i][3] =  1.67240374;
  // Nd-Fe
  i = 2;
  adp.param_dp[i][0] = -0.74178636;
  adp.param_dp[i][1] =  0.02627717;
  adp.param_dp[i][2] =  0.63990280;
  adp.param_dp[i][3] =  0.71822670;
  // Nd-B
  i = 3;
  adp.param_dp[i][0] = -2.77207819;
  adp.param_dp[i][1] = -0.01428758;
  adp.param_dp[i][2] =  3.00687447;
  adp.param_dp[i][3] =  1.99990662;
  // Fe-Fe
  i = 4;
  adp.param_dp[i][0] = 48.31013255;
  adp.param_dp[i][1] =  2.59579481;
  adp.param_dp[i][2] =  0.00174346;
  adp.param_dp[i][3] =  0.50000398;
  // Fe-B
  i = 5;
  adp.param_dp[i][0] = -6.05491786;
  adp.param_dp[i][1] =  0.00424405;
  adp.param_dp[i][2] =  5.92955891;
  adp.param_dp[i][3] =  1.99999157;
  // B-B
  i = 6;
  adp.param_dp[i][0] =  4.10939731;
  adp.param_dp[i][1] =  0.03873231;
  adp.param_dp[i][2] = -3.59639067;
  adp.param_dp[i][3] =  1.99998753;
  // w -----------------------------------------------------------
  // Nd-Nd
  i = 1;
  adp.param_qp[i][0] = -989.03111616;
  adp.param_qp[i][1] =    5.07115169;
  adp.param_qp[i][2] =   -0.00626730;
  adp.param_qp[i][3] =    1.35019925;
  // Nd-Fe
  i = 2;
  adp.param_qp[i][0] =   -6.85737561;
  adp.param_qp[i][1] =    1.77017779;
  adp.param_qp[i][2] =   -0.00326452;
  adp.param_qp[i][3] =    1.11944983;
  // Nd-B
  i = 3;
  adp.param_qp[i][0] =   -0.33707649;
  adp.param_qp[i][1] =    0.00067309;
  adp.param_qp[i][2] =    0.33845345;
  adp.param_qp[i][3] =    1.99985367;
  // Fe-Fe
  i = 4;
  adp.param_qp[i][0] =   -8.80329569;
  adp.param_qp[i][1] =    1.90125390;
  adp.param_qp[i][2] =    0.00008390;
  adp.param_qp[i][3] =    0.50010750;
  // Fe-B
  i = 5;
  adp.param_qp[i][0] =    0.27724858;
  adp.param_qp[i][1] =    0.01682020;
  adp.param_qp[i][2] =   -0.25292302;
  adp.param_qp[i][3] =    0.80554165;
  // B-B
  i = 6;
  adp.param_qp[i][0] =  290.40620566;
  adp.param_qp[i][1] =    4.29054527;
  adp.param_qp[i][2] =    0.01740259;
  adp.param_qp[i][3] =    0.68518535;
}

void adp_alloc()
{
  if (adp.rho)        { delete[] adp.rho;       adp.rho       = NULL; }
  if (adp.mu_x)       { delete[] adp.mu_x;      adp.mu_x      = NULL; }
  if (adp.mu_y)       { delete[] adp.mu_y;      adp.mu_y      = NULL; }
  if (adp.mu_z)       { delete[] adp.mu_z;      adp.mu_z      = NULL; }
  if (adp.nu)         { delete[] adp.nu;        adp.nu        = NULL; }
  if (adp.gradF)      { delete[] adp.gradF;     adp.gradF     = NULL; }
  if (adp.lambda_xx)  { delete[] adp.lambda_xx; adp.lambda_xx = NULL; }
  if (adp.lambda_yy)  { delete[] adp.lambda_yy; adp.lambda_yy = NULL; }
  if (adp.lambda_zz)  { delete[] adp.lambda_zz; adp.lambda_zz = NULL; }
  if (adp.lambda_xy)  { delete[] adp.lambda_xy; adp.lambda_xy = NULL; }
  if (adp.lambda_yz)  { delete[] adp.lambda_yz; adp.lambda_yz = NULL; }
  if (adp.lambda_zx)  { delete[] adp.lambda_zx; adp.lambda_zx = NULL; }
  if (adp.param_pair) {
    for (int i=0; i<adp.nptype+1; i++) { delete[] adp.param_pair[i]; }
    delete[] adp.param_pair; adp.param_pair= NULL; }
  if (adp.param_den)  { 
    for (int i=0; i<adp.ntype+1; i++) { delete[] adp.param_den[i]; }
    delete[] adp.param_den; adp.param_den = NULL; }
  if (adp.param_eb)   { 
    for (int i=0; i<adp.ntype+1; i++) { delete[] adp.param_eb[i]; }
    delete[] adp.param_eb;  adp.param_eb  = NULL; }
  if (adp.param_dp)   {
    for (int i=0; i<adp.nptype+1; i++) { delete[] adp.param_dp[i]; }
    delete[] adp.param_dp;  adp.param_dp  = NULL; }
  if (adp.param_qp)   {
    for (int i=0; i<adp.nptype+1; i++) { delete[] adp.param_qp[i]; }
    delete[] adp.param_qp;  adp.param_qp  = NULL; }
  if (adp.type)     { delete[] adp.type;     adp.type     = NULL; }
  if (adp.typen)    { delete[] adp.typen;    adp.typen    = NULL; }
  if (adp.ptype) {
    for (int i=0; i<adp.ntype+1; i++) { delete[] adp.ptype[i]; }
    delete[] adp.ptype; adp.ptype = NULL; }

// Kubo 20140224 -------------//
  for(int i=0; i<adp.nptype+1; i++){
    if(adp.Phi!=NULL){
    if(adp.Phi[i]!=NULL){
      delete adp.Phi[i];
      adp.Phi[i] = NULL;
    }
    }
    if(adp.U!=NULL){
    if(adp.U[i]!=NULL){
      delete adp.U[i];
      adp.U[i] = NULL;
    }
    }
    if(adp.W!=NULL){
    if(adp.W[i]!=NULL){
      delete adp.W[i];
      adp.W[i] = NULL;
    }
    }
  }
  for(int i=0; i<adp.ntype+1; i++){
    if(adp.Rho!=NULL){
    if(adp.Rho[i]!=NULL){
      delete adp.Rho[i];
      adp.Rho[i] = NULL;
    }
    }
    if(adp.F!=NULL){
    if(adp.F[i]!=NULL){
      delete adp.F[i];
      adp.F[i] = NULL;
    }
    }
  }
  if(adp.Phi!=NULL){delete [] adp.Phi; adp.Phi = NULL;}
  if(adp.Rho!=NULL){delete [] adp.Rho; adp.Rho = NULL;}
  if(adp.F  !=NULL){delete [] adp.F;   adp.F   = NULL;}
  if(adp.U  !=NULL){delete [] adp.U;   adp.U   = NULL;}
  if(adp.W  !=NULL){delete [] adp.W;   adp.W   = NULL;}

// End -----------------------//

  adp.rho        = new double[atom.natom+1];
  adp.mu_x       = new double[atom.natom+1];
  adp.mu_y       = new double[atom.natom+1];
  adp.mu_z       = new double[atom.natom+1];
  adp.nu         = new double[atom.natom+1];
  adp.gradF      = new double[atom.natom+1];
  adp.lambda_xx  = new double[atom.natom+1];
  adp.lambda_yy  = new double[atom.natom+1];
  adp.lambda_zz  = new double[atom.natom+1];
  adp.lambda_xy  = new double[atom.natom+1];
  adp.lambda_yz  = new double[atom.natom+1];
  adp.lambda_zx  = new double[atom.natom+1];

// Kubo 20140225 ---------------
  adp.param_pair = new double*[adp.nptype+1];
//  for (int i=0; i<adp.nptype+1; i++) { adp.param_pair[i] = new double[7]; }
  adp.param_den  = new double*[adp.ntype+1];
//  for (int i=0; i<adp.ntype+1; i++) { adp.param_den[i] = new double[5]; }
  adp.param_eb   = new double*[adp.ntype+1];
//  for (int i=0; i<adp.ntype+1; i++) { adp.param_eb[i] = new double[5]; }
  adp.param_dp   = new double*[adp.nptype+1];
//  for (int i=0; i<adp.nptype+1; i++) { adp.param_dp[i] = new double[4]; }
  adp.param_qp   = new double*[adp.nptype+1];
//  for (int i=0; i<adp.nptype+1; i++) { adp.param_qp[i] = new double[4]; }
// End -------------------------

  adp.type     = new int[adp.ntype+1];
  adp.ptype    = new int*[adp.ntype+1];
  for (int i=0; i<adp.ntype+1; i++) {
    adp.ptype[i] = new int[adp.ntype+1]; }
  adp.typen    = new int[atom.natom+1];

// Kubo 20140224 -------------//
  adp.Phi = new PairFunction*[adp.nptype+1];
  adp.Rho = new PairFunction*[adp.ntype+1];
  adp.F   = new PairFunction*[adp.ntype+1];
  adp.U   = new PairFunction*[adp.nptype+1];
  adp.W   = new PairFunction*[adp.nptype+1];

/* Later Create Instances
  for(int i=0; i<adp.nptype+1; i++){
    adp.Phi[i] = new PairFunction();
    adp.U[i]   = new PairFunction();
    adp.W[i]   = new PairFunction();
  }
  for(int i=0; i<adp.ntype+1; i++){
    adp.Rho[i] = new PairFunction();
    adp.F[i]   = new PairFunction();
  }
*/

// End -----------------------//
}

int adptyp(char* at)
{
  int i = atom_number(at);
  /*
  if (i == 40) { return 1;
  } else if (i ==  8) { return 2;
  } else if (i == 39) { return 3;
  } else { return 0; }
  */
  for (int j=1; j<=adp.ntype; j++) {
    if (i == adp.type[j]) { return j; }
  }
  return 0;
}
int adptyp(int at)
{
  int i = atom.anum[at];
  /*
  if (i == 40) { return 1;
  } else if (i ==  8) { return 2;
  } else if (i == 39) { return 3;
  } else { return 0; }
  */
  for (int j=1; j<=adp.ntype; j++) {
    if (i == adp.type[j]) { return j; }
  }
  return 0;
}

int adp_read_number(const char* fname)
{
  std::ifstream fin( fname );
  std::string line, arg1;
  while (getline(fin, line)) { remove_carriage_return(line);
    if (line.length()==0) { continue; }
    remove_head_spaces(line);
    if (line.at(0)!='#') { break; }
  }
  remove_after_sharp(line);
  int narg = count_arg_number(line); adp.ntype = narg;
  while (getline(fin, line)) { remove_carriage_return(line);
    if (line.length()==0) { continue; }
    remove_head_spaces(line);
    if (line.at(0)!='#') { break; }
  }
  int np = narg * (narg -1) / 2 + narg;
  remove_after_sharp(line);
  narg = count_arg_number(line);
  if (narg < np) {
    printf("########################################\n");
    printf("ADP: Potential file error\n (incorrect flag of pair interaction)\n");
    printf(" narg = %d, np = %d\n",narg,np);
    printf("########################################\n");
    fin.close(); fin.clear();
    return 0;
  }
  adp.nptype = 0;
  for (int i=0; i<np; i++) {
    get_first_arg(line,arg1);
    if (strcmp(arg1.c_str(), "1") == 0) { adp.nptype++;
    } else if (strcmp(arg1.c_str(), "Y") == 0) { adp.nptype++;
    } else if (strcmp(arg1.c_str(), "y") == 0) { adp.nptype++;
    } else if (strcmp(arg1.c_str(), "T") == 0) { adp.nptype++;
    } else if (strcmp(arg1.c_str(), "t") == 0) { adp.nptype++;
    }
  }
  fin.close(); fin.clear();
  return 1;
}
int adp_read_param(const char* fname)
{
  std::ifstream fin( fname );
  std::string line;
  std::string arg1;
  char larg1[40];
  // Atom species
  while (getline(fin, line)) { remove_carriage_return(line);
    if (line.length()==0) { continue; }
    remove_head_spaces(line);
    if (line.at(0)!='#') { remove_after_sharp(line); break; }
  }
  for (int i=1; i<=adp.ntype; i++) {
    get_first_arg(line,arg1);
    char sp[3] = "aa"; strcpy(sp,arg1.c_str());
    adp.type[i] = atom_number(sp); //printf("%d %d\n",i,adp.type[i]);
  }
  // Pair-interaction flags
  while (getline(fin, line)) { remove_carriage_return(line);
    if (line.length()==0) { continue; }
    remove_head_spaces(line);
    if (line.at(0)!='#') { remove_after_sharp(line); break; }
  }
  int narg = adp.ntype; 
  int pt = 1;
  for (int i=1; i<=narg; i++) {
    for (int j=i; j<=narg; j++) {
      get_first_arg(line,arg1);
      if (strcmp(arg1.c_str(), "1") == 0) { adp.ptype[i][j] = pt;
      } else if (strcmp(arg1.c_str(), "Y") == 0) { adp.ptype[i][j] = pt;
      } else if (strcmp(arg1.c_str(), "y") == 0) { adp.ptype[i][j] = pt;
      } else if (strcmp(arg1.c_str(), "T") == 0) { adp.ptype[i][j] = pt;
      } else if (strcmp(arg1.c_str(), "t") == 0) { adp.ptype[i][j] = pt;
      } else { adp.ptype[i][j] = 0; }
      //printf("%d %d %d\n",i,j,pt);
      if (i != j) { adp.ptype[j][i] = adp.ptype[i][j]; }
      if (adp.ptype[i][j] > 0) { pt++; }
    }
  }
  pt--;
  // read and set parameters
  double xx;
  // cutoff
  adp_find_param(fin,"r_cut",line);
  get_first_arg(line, arg1); strcpy (larg1, arg1.c_str()); xx = atof(larg1);
  adp.cut = xx; //printf("rcut %f\n",xx);

//std::cout<<"(^Q^) Find Me!"<<std::endl; exit(0);

  // pair ------------------------------------------------------
  adp_find_param(fin,"pair",line);
  get_first_arg(line, arg1);

/* Original Edition
  if (strcmp(arg1.c_str(), "mishin_sc")!=0) { printf("Error in ADP param\n"); return 0; }
  for (int j=0; j<7; j++) {
    getline (fin, line); //std::cout<<"L: "<<line<<std::endl;
    for (int i=1; i<=pt; i++) {
      get_first_arg(line, arg1);
      strcpy (larg1, arg1.c_str());
      xx = atof(larg1);
      adp.param_pair[i][j] = xx; //printf("pair %d %d %f\n",i,j,xx);
    }
  }
//*/

// Kubo 20140224 -------------//
  {
  using namespace std;

  PairFunction *ghost = new PairFunction();
  int iNumParam = ghost->GetNumOfParam(arg1);
  string functype = arg1;


  for(int i=0; i<=adp.nptype; i++){adp.param_pair[i] = new double[iNumParam];}

  for(int j=0; j<iNumParam; j++){
    getline(fin, line); remove_carriage_return(line);
    for(int i=1; i<=adp.nptype; i++){
      get_first_arg(line, arg1);
      strcpy(larg1, arg1.c_str());
      xx = atof(larg1);
      adp.param_pair[i][j] = xx;
    }
  }

  adp.Phi[0] = new PairFunction();
  for(int i=1; i<=adp.nptype; i++){
//    adp.Phi[i] = new PairFunction(functype,adp.param_pair[i],adp.cut);
    stringstream ss; ss << "PHI_" << i << "(" << functype << ")";
    adp.Phi[i] = new PairFunction(functype,adp.param_pair[i],adp.cut,ss.str());
  }

  delete ghost; ghost = NULL;

//  adp.Phi[1]->Plot(2,7,1000,"test.out"); exit(0);
  }
// End -----------------------//


  // den -------------------------------------------------------
  adp_find_param(fin,"den",line);
  get_first_arg(line, arg1);

/* Original Edition
  if (strcmp(arg1.c_str(), "csw2_sc")!=0) { printf("Error in ADP param\n"); return 0; }
  for (int j=0; j<5; j++) {
    getline (fin, line); //std::cout<<"L2: "<<line<<std::endl;
    for (int i=1; i<=adp.ntype; i++) {
      get_first_arg(line, arg1);
      strcpy (larg1, arg1.c_str());
      xx = atof(larg1);
      adp.param_den[i][j] = xx; //printf("den %d %d %f\n",i,j,xx);
    }
  }
//*/

// Kubo 20140225 -------------//
  {
  using namespace std;

  PairFunction *ghost = new PairFunction();
  int iNumParam = ghost->GetNumOfParam(arg1);
  string functype = arg1;

  for(int i=0; i<=adp.ntype; i++){adp.param_den[i] = new double[iNumParam];}

  for(int j=0; j<iNumParam; j++){
    getline(fin, line); remove_carriage_return(line);
    for(int i=1; i<=adp.ntype; i++){
      get_first_arg(line, arg1);
      strcpy(larg1, arg1.c_str());
      xx = atof(larg1);
      adp.param_den[i][j] = xx;
    }
  }

  adp.Rho[0] = new PairFunction();
  for(int i=1; i<=adp.ntype; i++){
//    adp.Rho[i] = new PairFunction(functype,adp.param_den[i],adp.cut);
    stringstream ss; ss << "RHO_" << i << "(" << functype << ")";
    adp.Rho[i] = new PairFunction(functype,adp.param_den[i],adp.cut,ss.str());
  }

  delete ghost; ghost = NULL;

//  adp.Rho[1]->Plot(2,7,1000,"test.out"); exit(0);
  }
// End 20140225 --------------//


  //embed ------------------------------------------------------
  adp_find_param(fin,"embed",line);
  get_first_arg(line, arg1);

/* Original Edition
  if(strcmp(arg1.c_str(), "poly_5")!=0){ printf("Error in ADP param\n"); return 0; }
  for(int j=0; j<5; j++){
    getline(fin, line);
    for(int i=1; i<=adp.ntype; i++) {
      get_first_arg(line, arg1);
      strcpy(larg1, arg1.c_str());
      xx = atof(larg1);
      adp.param_eb[i][j] = xx; //printf("eb %d %d %f\n",i,j,xx);
    }
  }
//*/

// Kubo 20140225 -------------//
  {
  using namespace std;

  PairFunction *ghost = new PairFunction();
  int iNumParam = ghost->GetNumOfParam(arg1);
  string functype = arg1;

  for(int i=0; i<=adp.ntype; i++){adp.param_eb[i] = new double[iNumParam];}

  for(int j=0; j<iNumParam; j++){
    getline(fin, line); remove_carriage_return(line);
    for(int i=1; i<=adp.ntype; i++){
      get_first_arg(line, arg1);
      strcpy(larg1, arg1.c_str());
      xx = atof(larg1);
      adp.param_eb[i][j] = xx;
    }
  }

  adp.F[0] = new PairFunction();
  for(int i=1; i<=adp.ntype; i++){
//    adp.F[i] = new PairFunction(functype,adp.param_eb[i],0.0);
    stringstream ss; ss << "F_" << i << "(" << functype << ")";
    adp.F[i] = new PairFunction(functype,adp.param_eb[i],0.0,ss.str());
  }

  delete ghost; ghost = NULL;

//  adp.F[1]->Plot(0,5,1000,"test.out"); exit(0);
  }
// End 20140225 --------------//


  //dipole -----------------------------------------------------
  adp_find_param(fin,"dipole",line);
  get_first_arg(line, arg1);

/* Original Edition
  if (strcmp(arg1.c_str(), "exp_plus_sc")!=0) { printf("Error in ADP param\n"); return 0; }
  for (int j=0; j<4; j++) {
    getline (fin, line);
    for (int i=1; i<=pt; i++) {
      get_first_arg(line, arg1);
      strcpy (larg1, arg1.c_str());
      xx = atof(larg1);
      adp.param_dp[i][j] = xx; //printf("dp %d %d %f\n",i,j,xx);
    }
  }
//*/

// Kubo 20140225 -------------//
  {
  using namespace std;

  PairFunction *ghost = new PairFunction();
  int iNumParam = ghost->GetNumOfParam(arg1);
  string functype = arg1;

  for(int i=0; i<=adp.nptype; i++){adp.param_dp[i] = new double[iNumParam];}

  for(int j=0; j<iNumParam; j++){
    getline(fin, line); remove_carriage_return(line);
    for(int i=1; i<=adp.nptype; i++){
      get_first_arg(line, arg1);
      strcpy(larg1, arg1.c_str());
      xx = atof(larg1);
      adp.param_dp[i][j] = xx;
    }
  }

  adp.U[0] = new PairFunction();
  for(int i=1; i<=adp.nptype; i++){
//    adp.U[i] = new PairFunction(functype,adp.param_dp[i],adp.cut);
    stringstream ss; ss << "U_" << i << "(" << functype << ")";
    adp.U[i] = new PairFunction(functype,adp.param_dp[i],adp.cut,ss.str());
  }

  delete ghost; ghost = NULL;

//  adp.U[1]->Plot(2,7,1000,"test.out"); exit(0);
  }
// End 20140225 --------------//


  //quadrupole -------------------------------------------------
  adp_find_param(fin,"quad",line);
  get_first_arg(line, arg1);

/* Original Edition
  if (strcmp(arg1.c_str(), "exp_plus_sc")!=0) { printf("Error in ADP param\n"); return 0; }
  for (int j=0; j<4; j++) {
    getline (fin, line);
    for (int i=1; i<=pt; i++) {
      get_first_arg(line, arg1);
      strcpy (larg1, arg1.c_str());
      xx = atof(larg1);
      adp.param_qp[i][j] = xx; //printf("qp %d %d %f\n",i,j,xx);
    }
  }
//*/

// Kubo 20140225 -------------//
  {
  using namespace std;

  PairFunction *ghost = new PairFunction();
  int iNumParam = ghost->GetNumOfParam(arg1);
  string functype = arg1;

  for(int i=0; i<=adp.nptype; i++){adp.param_qp[i] = new double[iNumParam];}

  for(int j=0; j<iNumParam; j++){
    getline(fin, line); remove_carriage_return(line);
    for(int i=1; i<=adp.nptype; i++){
      get_first_arg(line, arg1);
      strcpy(larg1, arg1.c_str());
      xx = atof(larg1);
      adp.param_qp[i][j] = xx;
    }
  }

  adp.W[0] = new PairFunction();
  for(int i=1; i<=adp.nptype; i++){
//    adp.W[i] = new PairFunction(functype,adp.param_qp[i],adp.cut);
    stringstream ss; ss << "W_" << i << "(" << functype << ")";
    adp.W[i] = new PairFunction(functype,adp.param_qp[i],adp.cut,ss.str());
  }

  delete ghost; ghost = NULL;

//  adp.W[1]->Plot(2,7,1000,"test.out"); exit(0);
  }
// End 20140225 --------------//


// Kubo 20140226 -------------//
/* Plot All Functions
{ using namespace std;
  for(int i=1;i<=adp.nptype;i++){
    stringstream ss;
    ss << "PHI_" << i << ".out";
    adp.Phi[i]->Plot(2,7,1000,ss.str());
  }
  for(int i=1;i<=adp.ntype;i++){
    stringstream ss;
    ss << "RHO_" << i << ".out";
    adp.Rho[i]->Plot(2,7,1000,ss.str());
  }
  for(int i=1;i<=adp.ntype;i++){
    stringstream ss;
    ss << "F_" << i << ".out";
    adp.F[i]->Plot(0,5,1000,ss.str());
  }
  for(int i=1;i<=adp.nptype;i++){
    stringstream ss;
    ss << "U_" << i << ".out";
    adp.U[i]->Plot(2,7,1000,ss.str());
  }
  for(int i=1;i<=adp.nptype;i++){
    stringstream ss;
    ss << "W_" << i << ".out";
    adp.W[i]->Plot(2,7,1000,ss.str());
  }
  exit(0);
}
//*/
// End -----------------------//






  for (int i=1; i<=atom.natom; i++) {
    adp.typen[i] = adptyp(atom.asp[i]);
    if (adp.typen[i] == 0) { return 0; }
  }
  return 1;
}

int adp_find_param(const char* fname, const char* tag, std::string &line)
{
  std::ifstream fin(fname);
  std::string arg1;
  int flg = 0;
  for (int i=1; i<=100; i++) {
    getline (fin, line); remove_carriage_return(line);
    get_first_arg(line, arg1);
    if (strcmp(arg1.c_str(), tag)==0) {
      flg = 1; break;
    }
  }
  fin.close(); fin.clear();
  return flg;
}
int adp_find_param(std::istream& fin, const char* tag, std::string &line)
{
  //  std::ifstream fin(fname);
  std::string arg1;
  int flg = 0;
  for (int i=1; i<=100; i++) {
    getline (fin, line); remove_carriage_return(line);
    get_first_arg(line, arg1);
    if (strcmp(arg1.c_str(), tag)==0) {
      flg = 1; break;
    }
  }
  //fin.close(); fin.clear();
  return flg;
}

