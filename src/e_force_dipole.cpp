#include <iostream>
#include<string>
#include<stdlib.h>
#include<fstream>
#include "myheader.h"

#if !defined __linux__ && !defined __APPLE__
#include <math.h>
double erf(double x);
double erfc(double x);
#endif

#ifndef signbit
#define signbit(x) ((x) < 0.0)
#endif

#define MB_SMOOTH_RANGE 0.1

void dipole_delete();
void dipole_alloc();
int dipole_param();
void bmh(double r, double &e, double &f, double a, double sig, double c);
void buckingham(double r, double &e, double &f, double a, double sig, double c);
void bucksft(double r, double &e, double &f, double a, double sig, double c,
	     double rc, double vrc, double vprc);
void elstat_shift(double r, double dp_kappa, double &value_tail, double &grad_tail,
		  double &ggrad_tail);
void elstat_value(double r, double dp_kappa, double &ftail, double &gtail,
		  double &ggtail);
double shortrange_value(double r, double a, double b, double c);
void shortrange_term(double r, double b, double c, double &srval_tail,
		     double &srgrad_tail);
int atom_number(char* at);
int dptyp(char* at);
//int dptyp(int at);
int pairtyp(char* ati, char* atj);
int pairtyp(int ati, int atj);
double dsquare(double d);
void resetmat(double a[3][3]);
int dipole_read_number(const char* fname);
int dipole_read_param(const char* fname);
void get_first_arg(std::string &line, std::string &arg1);
int remove_head_spaces(std::string &line);
int count_arg_number(std::string line);
int remove_after_sharp(std::string &line);
int dipole_find_param(const char* fname, const char* tag, std::string &line);
void remove_carriage_return(std::string& line);

void e_force_dipole(PotentialMode mode)
{
  double rr, rr2, drx, dry, drz, vp0, v0;
  int j, ix, iy, iz;
  //double rcut = 8.0e0; double rcut2 = rcut * rcut;
  int    self;
  int    h, k, l, typ1, typ2, uf, us, stresses, ipair;
  double value, grad, value_tail, grad_tail, grad_i, grad_j, p_sr_tail;
  double value_el, grad_el, ggrad_el;
  //dipole.dp_tol = 1.0e-12;
  //dipole.shift = true;

  if (dipole.initialize) {
    std::cout<<"Initialize of the DIPOLE potential.."<<std::endl;
    dipole_delete();
    printf(" Dipole potential file = %s\n",dipole.fname);
    dipole_read_number(dipole.fname);
    printf(" # of species = %d, # of pairs = %d\n",dipole.ntype, dipole.nptype);
    //dipole.ntype = 3; dipole.nptype = 3;
    dipole_alloc();
    if (dipole_read_param(dipole.fname) == 0 && mode != COMBINEDPOTMODE) {
      printf("Dipole error (unsupported atom)\n"); mdmotion=0; return; }
    //if (dipole_param() == 0) { printf("Dipole error (unsupported atom)\n"); mdmotion=0; return; }
    dipole.initialize = false;
    std::cout<<"DIOPLE potential initialization done."<<std::endl;
  }

  if (mode != COMBINEDPOTMODE)
    rcut = 0;

  for (int i=1; i<=dipole.nptype; i++) {
    if (rcut < dipole.r_cut[i]) { rcut = dipole.r_cut[i]; } }

  if (rcut < dipole.dp_cut) rcut = dipole.dp_cut;
  rcut *=ang;
  
  rcut2 = rcut * rcut;

  /* General rcut is used, OK if only one potential is used, othewise be careful */
  /* make sure that distance for nieghbor criterion is large enough*/

  // Reset dipoles and fields: LOOP ZERO
  for (int i=1; i<=atom.natom; i++) {
    //dipole.E_oldx[i] = 0; dipole.E_oldy[i] = 0; dipole.E_oldz[i] = 0;
    dipole.E_totx[i] = 0; dipole.E_toty[i] = 0; dipole.E_totz[i] = 0;
    dipole.E_indx[i] = 0; dipole.E_indy[i] = 0; dipole.E_indz[i] = 0;
    dipole.p_indx[i] = 0; dipole.p_indy[i] = 0; dipole.p_indz[i] = 0;
    dipole.E_statx[i] = 0; dipole.E_staty[i] = 0; dipole.E_statz[i] = 0;
    dipole.p_srx[i]   = 0; dipole.p_sry[i]   = 0; dipole.p_srz[i]   = 0; }

  // FIRST LOOP (set zero atomic forces is it is not a combined mode, otherwise it was done before)
  if (mode != COMBINEDPOTMODE) {
    for (int i=1; i<=atom.natom; i++) {
      atom.fx[i] = 0.0; atom.fy[i] = 0.0; atom.fz[i] = 0.0;  atom.epot[i] = 0.0;
      for (int j=0; j<3; j++) { for (int k=0; k<3; k++) {
	  atom.satom[i][j][k] = 0.0; } }
    }// FIRST LOOP END
  }

  // SECOND LOOP: calculate short-range and monopole forces,
  // calculate static field- and dipole-contributions
  for (int i=1; i<=atom.natom; i++) {
    //typ1 = dptyp(atom.asp[i]);
    //typ1 = dptyp(i);
    typ1 = dipole.typen[i];

    //if this atomic species should not be accounted in the dipole model, then skip it
    if (typ1 <= 0)
    {   
    	continue;
    }

    if (book.alistnum[i]>0) { for (int k=1; k<=book.alistnum[i]; k++) {
	j  = book.alist[i][k][0];
	if (j<i) continue;
	ix = book.alist[i][k][1]; iy = book.alist[i][k][2]; iz = book.alist[i][k][3];
	rr2 = atom.Dist2(i,j,ix,iy,iz)/ang/ang; rr = sqrt(rr2);
	//typ2 = dptyp(atom.asp[j]);
	//typ2 = dptyp(j);
	typ2 = dipole.typen[j];
        //if there atomic species should not be accounted in the dipole model, then skip it
    	if (typ2 <= 0) continue;
	//ipair = pairtyp(atom.asp[i],atom.asp[j]); // type of pair interaction
	//ipair = pairtyp(i,j); // type of pair interaction
	ipair = dipole.ptype[typ1][typ2];
	elstat_shift(rr, dipole.dp_kappa, value_el, grad_el, ggrad_el); //tail-functions
	drx = atom.Dx(i,j,ix,iy,iz); dry = atom.Dy(i,j,ix,iy,iz); drz = atom.Dz(i,j,ix,iy,iz);
	drx /= ang; dry /= ang; drz /= ang;
	//if ((rr < rcut)&&(ipair!=0)) { // calculate short-range forces
	if ((rr < dipole.r_cut[ipair])&&(ipair!=0)) { // calculate short-range forces
	  self = 0; if (i==j) { self = 1; }
	  // For consistency with POTFIT, buckingham should be replaced with bucksft.
	  // For bucksft, buck_vprc and buck_vrc must be given in initialization.
	  // (However, bucksft seems weird...)
	  if (!dipole.shift) {
	    buckingham(rr, value, grad, dipole.buck_a[ipair], dipole.buck_s[ipair]
		       ,dipole.buck_c[ipair]); //Buckingham type
	  } else {
	    bucksft(rr, value, grad, dipole.buck_a[ipair], dipole.buck_s[ipair]
		    ,dipole.buck_c[ipair]
		    ,dipole.r_cut[ipair],dipole.buck_vrc[ipair],dipole.buck_vprc[ipair]);
	  }
	  if (self) { value *= 0.5; grad *= 0.5; }
	  grad /= rr;
	  atom.epot[i] += value/2.0; atom.epot[j] += value/2.0; // -3 contribution. OK!!
	  atom.fx[i] += drx*grad; atom.fy[i] += dry*grad; atom.fz[i] += drz*grad;
	  atom.fx[j] -= drx*grad; atom.fy[j] -= dry*grad; atom.fz[j] -= drz*grad;
	  // stress
	  atom.satom[i][0][0] += drx*drx*grad/2.0;
	  atom.satom[i][0][1] += dry*drx*grad/2.0;
	  atom.satom[i][0][2] += drz*drx*grad/2.0;
	  atom.satom[i][1][1] += dry*dry*grad/2.0;
	  atom.satom[i][1][2] += drz*dry*grad/2.0;
	  atom.satom[i][2][2] += drz*drz*grad/2.0;
	  atom.satom[j][0][0] += drx*drx*grad/2.0;
	  atom.satom[j][0][1] += dry*drx*grad/2.0;
	  atom.satom[j][0][2] += drz*drx*grad/2.0;
	  atom.satom[j][1][1] += dry*dry*grad/2.0;
	  atom.satom[j][1][2] += drz*dry*grad/2.0;
	  atom.satom[j][2][2] += drz*drz*grad/2.0;
	}
	if (rr < dipole.dp_cut) { // calculate monopole forces
	  self = 0; if (i==j) { self = 1; }
	  value_tail = value_el;
	  grad_tail  = grad_el;
	  grad_i = dipole.charge[typ2] * grad_tail;
	  if (typ1 == typ2) { grad_j = grad_i;
	  } else { grad_j = dipole.charge[typ1] * grad_tail; }
	  value = dipole.charge[typ1] * dipole.charge[typ2] * value_tail;
	  grad  = dipole.charge[typ1] * grad_i;
	  if (self) { grad_i *= 0.5; grad_j *= 0.5; value *= 0.5; grad *= 0.5; }
	  atom.epot[i] += value/2.0; atom.epot[j] += value/2.0; // -18 contribution DUBIOUS!!
	  atom.fx[i] += drx*grad; atom.fy[i] += dry*grad; atom.fz[i] += drz*grad;
	  atom.fx[j] -= drx*grad; atom.fy[j] -= dry*grad; atom.fz[j] -= drz*grad;
	  // stress
	  atom.satom[i][0][0] += drx*drx*grad/2.0;
	  atom.satom[i][0][1] += dry*drx*grad/2.0;
	  atom.satom[i][0][2] += drz*drx*grad/2.0;
	  atom.satom[i][1][1] += dry*dry*grad/2.0;
	  atom.satom[i][1][2] += drz*dry*grad/2.0;
	  atom.satom[i][2][2] += drz*drz*grad/2.0;
	  atom.satom[j][0][0] += drx*drx*grad/2.0;
	  atom.satom[j][0][1] += dry*drx*grad/2.0;
	  atom.satom[j][0][2] += drz*drx*grad/2.0;
	  atom.satom[j][1][1] += dry*dry*grad/2.0;
	  atom.satom[j][1][2] += drz*dry*grad/2.0;
	  atom.satom[j][2][2] += drz*drz*grad/2.0;
	  // calculate static field-contributions
	  dipole.E_statx[i] += drx*grad_i; dipole.E_statx[j] -= drx*grad_j;
	  dipole.E_staty[i] += dry*grad_i; dipole.E_staty[j] -= dry*grad_j;
	  dipole.E_statz[i] += drz*grad_i; dipole.E_statz[j] -= drz*grad_j;
	  // calculate short-range dipoles
	  if (ipair!=0) {
	  p_sr_tail = grad_tail*rr
	    *shortrange_value(rr, dipole.dp_alpha[typ1], dipole.dp_b[ipair],
			      dipole.dp_c[ipair]);
	  dipole.p_srx[i] += dipole.charge[typ2]*drx/rr*p_sr_tail;
	  dipole.p_sry[i] += dipole.charge[typ2]*dry/rr*p_sr_tail;
	  dipole.p_srz[i] += dipole.charge[typ2]*drz/rr*p_sr_tail;
	  if (!self) {
	    p_sr_tail = grad_tail*rr
	      *shortrange_value(rr, dipole.dp_alpha[typ2], dipole.dp_b[ipair],
				dipole.dp_c[ipair]);
	    dipole.p_srx[j] -= dipole.charge[typ1]*drx/rr*p_sr_tail;
	    dipole.p_sry[j] -= dipole.charge[typ1]*dry/rr*p_sr_tail;
	    dipole.p_srz[j] -= dipole.charge[typ1]*drz/rr*p_sr_tail; }
	  }
	}
      } }
  } // SECOND LOOP END

  // THIRD LOOP: calculate whole dipole moment for every atom
  double rp, dp_sum;
  int dp_converged = 0, dp_it = 0;
  double max_diff = 10;
  //goto OUT;
  while (dp_converged == 0) {
    dp_sum = 0;
    for (int i=1; i<=atom.natom; i++) {
      //typ1 = dptyp(atom.asp[i]);
      //typ1 = dptyp(i);
      typ1 = dipole.typen[i];
      //if this atomic species should not be accounted in the dipole model, then skip it
      if (typ1 <= 0) continue;
      if (dipole.dp_alpha[typ1]!=0) {
      if (dp_it) {
	dipole.E_totx[i] = (1-dipole.dp_mix)*dipole.E_indx[i]
	  + dipole.dp_mix*dipole.E_oldx[i] + dipole.E_statx[i];
	dipole.E_toty[i] = (1-dipole.dp_mix)*dipole.E_indy[i]
	  + dipole.dp_mix*dipole.E_oldy[i] + dipole.E_staty[i];
	dipole.E_totz[i] = (1-dipole.dp_mix)*dipole.E_indz[i]
	  + dipole.dp_mix*dipole.E_oldz[i] + dipole.E_statz[i];
      } else {
	dipole.E_totx[i] = dipole.E_indx[i] + dipole.E_statx[i];
	dipole.E_toty[i] = dipole.E_indy[i] + dipole.E_staty[i];
	dipole.E_totz[i] = dipole.E_indz[i] + dipole.E_statz[i];
      } // if dp_it
      dipole.p_indx[i] = dipole.dp_alpha[typ1] * dipole.E_totx[i] + dipole.p_srx[i];
      dipole.p_indy[i] = dipole.dp_alpha[typ1] * dipole.E_toty[i] + dipole.p_sry[i];
      dipole.p_indz[i] = dipole.dp_alpha[typ1] * dipole.E_totz[i] + dipole.p_srz[i];
      dipole.E_oldx[i] = dipole.E_indx[i]; dipole.E_indx[i] = 0.;
      dipole.E_oldy[i] = dipole.E_indy[i]; dipole.E_indy[i] = 0.;
      dipole.E_oldz[i] = dipole.E_indz[i]; dipole.E_indz[i] = 0.;
      }
    } // loop of i (atoms)

    for (int i=1; i<=atom.natom; i++) {
      //typ1 = dptyp(atom.asp[i]);
      //typ1 = dptyp(i);
      typ1 = dipole.typen[i];
      if (typ1 <= 0) continue;
      if (book.alistnum[i]>0) { for (int k=1; k<=book.alistnum[i]; k++) {
	  j  = book.alist[i][k][0];
	  if (j<i) continue;
	  ix = book.alist[i][k][1]; iy = book.alist[i][k][2]; iz = book.alist[i][k][3];
	  rr2 = atom.Dist2(i,j,ix,iy,iz)/ang/ang; rr = sqrt(rr2);
	  //typ2 = dptyp(atom.asp[j]);
	  //typ2 = dptyp(j);
	  typ2 = dipole.typen[j];
          if (typ2 == 0) continue;
	  //ipair = pairtyp(atom.asp[i],atom.asp[j]); // type of pair interaction
	  //ipair = pairtyp(i,j); // type of pair interaction
	  ipair = dipole.ptype[typ1][typ2];
	  elstat_shift(rr, dipole.dp_kappa, value_el, grad_el, ggrad_el); //tail-functions
	  if (rr < dipole.dp_cut) {
	    if ((dipole.dp_alpha[typ1]!=0)&&(dipole.dp_alpha[typ2]!=0)) {
	    self = 0; if (i==j) { self = 1; }
	    drx = atom.Dx(i,j,ix,iy,iz); dry = atom.Dy(i,j,ix,iy,iz); drz = atom.Dz(i,j,ix,iy,iz);
	    drx /= ang; dry /= ang; drz /= ang;
	    drx /=rr; dry /=rr; drz /=rr;
	    rp = dipole.p_indx[j]*drx + dipole.p_indy[j]*dry + dipole.p_indz[j]*drz;
	    // sign below is fixed YU(AK) 2013.10.23
	    dipole.E_indx[i] -= grad_el*(3*rp*drx-dipole.p_indx[j]);
	    dipole.E_indy[i] -= grad_el*(3*rp*dry-dipole.p_indy[j]);
	    dipole.E_indz[i] -= grad_el*(3*rp*drz-dipole.p_indz[j]);
	    if (!self) {
	      rp = dipole.p_indx[i]*drx + dipole.p_indy[i]*dry + dipole.p_indz[i]*drz;
	      // sign below is fixed YU(AK) 2013.10.23
	      dipole.E_indx[j] -= grad_el*(3*rp*drx-dipole.p_indx[i]);
	      dipole.E_indy[j] -= grad_el*(3*rp*dry-dipole.p_indy[i]);
	      dipole.E_indz[j] -= grad_el*(3*rp*drz-dipole.p_indz[i]);
	    }
	    }}
	} } // j
    } // loop of i (atoms)
    for (int i=1; i<=atom.natom; i++) {
      //typ1 = dptyp(atom.asp[i]);
      //typ1 = dptyp(i);
      typ1 = dipole.typen[i];
      if (typ1 <= 0) continue;
      dp_sum += dsquare(dipole.dp_alpha[typ1]*(dipole.E_oldx[i]-dipole.E_indx[i]));
      dp_sum += dsquare(dipole.dp_alpha[typ1]*(dipole.E_oldy[i]-dipole.E_indy[i]));
      dp_sum += dsquare(dipole.dp_alpha[typ1]*(dipole.E_oldz[i]-dipole.E_indz[i]));
    } // loop of i (atoms)

    //dp_sum /= 3*(double)atom.natom; // old version without combined mode
    dp_sum /= 3*(double)dipole.natoms; // new version with combined mode
    dp_sum = sqrt(dp_sum);

    if (dp_it) {
      if ((dp_sum > max_diff) || (dp_it > 50)) {
	dp_converged = 1;
	printf("dp (failure) dp_sum=%e  dp_it=%d\n",dp_sum,dp_it);
	for (int i=1; i<=atom.natom; i++) {
	  //typ1 = dptyp(atom.asp[i]);
	  //typ1 = dptyp(i);
	  typ1 = dipole.typen[i];
          if (typ1 <= 0) continue;
	  if (dipole.dp_alpha[typ1]!=0) {
	    dipole.p_indx[i] = dipole.dp_alpha[typ1]*dipole.E_statx[i]+dipole.p_srx[i];
	    dipole.p_indy[i] = dipole.dp_alpha[typ1]*dipole.E_staty[i]+dipole.p_sry[i];
	    dipole.p_indz[i] = dipole.dp_alpha[typ1]*dipole.E_statz[i]+dipole.p_srz[i];
	    dipole.E_indx[i] = dipole.E_statx[i];
	    dipole.E_indy[i] = dipole.E_staty[i];
	    dipole.E_indz[i] = dipole.E_statz[i]; }
	}
      }
    } // if dp_it
    if (dp_sum < dipole.dp_tol) {
      dp_converged = 1;
      //printf("dp (success) dp_it=%d\n",dp_it);
    }
    dp_it++;
  } // while dp_converged (THIRD LOOP END)
  //printf("%d %e\n",dp_it,dp_sum);

  // FOURTH LOOP: calculate monopole-dipole and dipole-dipole forces
  double  rp_i, rp_j, pp_ij, tmp_1, tmp_2, tmpx, tmpy, tmpz;
  double  grad_1, grad_2, srval, srgrad, srval_tail, srgrad_tail, value_sum,
    grad_sum;
  for (int i=1; i<=atom.natom; i++) {
    //typ1 = dptyp(atom.asp[i]);
    //typ1 = dptyp(i);
    typ1 = dipole.typen[i];
    if (typ1 <= 0) continue;
    if (book.alistnum[i]>0) { for (int k=1; k<=book.alistnum[i]; k++) {
	j  = book.alist[i][k][0];
	if (j<i) continue;
	ix = book.alist[i][k][1]; iy = book.alist[i][k][2]; iz = book.alist[i][k][3];
	rr2 = atom.Dist2(i,j,ix,iy,iz)/ang/ang; rr = sqrt(rr2);
	//typ2 = dptyp(atom.asp[j]);
	//typ2 = dptyp(j);
	typ2 = dipole.typen[j];
        if (typ2 <= 0) continue;
	//ipair = pairtyp(atom.asp[i],atom.asp[j]); // type of pair interaction
	//ipair = pairtyp(i,j); // type of pair interaction
	ipair = dipole.ptype[typ1][typ2];
	elstat_shift(rr, dipole.dp_kappa, value_el, grad_el, ggrad_el); //tail-functions
	if ((rr < dipole.dp_cut)&&((dipole.dp_alpha[typ1]!=0)||(dipole.dp_alpha[typ2]!=0))) {
	  self = 0; if (i==j) { self = 1; }
	  value_tail = -grad_el;
	  //grad_tail  = -ggrad_el;
	  grad_tail  = -ggrad_el*rr; //fixed by YU 2013.09.25
	  if (ipair!=0) {
	    shortrange_term(rr, dipole.dp_b[ipair], dipole.dp_c[ipair], srval_tail, srgrad_tail);
	    srval = value_tail * srval_tail;
	    srgrad = value_tail * srgrad_tail + grad_tail * srval_tail;
	  }
	  if (self) { value_tail *= 0.5; grad_tail *= 0.5; }
	  // monopole-dipole contributions //if (charge[typ1] && dp_alpha[typ2]) 
	  if (ipair!=0) {
	    value_sum = value_tail + srval;
	    grad_sum  = grad_tail + srgrad;
	  } else {
	    value_sum = value_tail;
	    grad_sum  = grad_tail;
	  }
	  drx = atom.Dx(i,j,ix,iy,iz); dry = atom.Dy(i,j,ix,iy,iz); drz = atom.Dz(i,j,ix,iy,iz);
	  drx /= ang; dry /= ang; drz /= ang;
	  drx /=rr; dry /=rr; drz /=rr;
	  rp_j = dipole.p_indx[j]*drx + dipole.p_indy[j]*dry + dipole.p_indz[j]*drz;
	  value  = dipole.charge[typ1] * rp_j * rr * value_sum;
	  //grad_1 = dipole.charge[typ1] * rp_j * rr2 * grad_sum ;
	  grad_1 = dipole.charge[typ1] * rp_j * rr * grad_sum ; //YU 2013.10.03
	  grad_2 = dipole.charge[typ1] * value_sum;
	  atom.epot[i] -= value/2.0; atom.epot[j] -= value/2.0;  // DUBIOUS!!!
	  tmpx = drx * grad_1 + dipole.p_indx[j] * grad_2;
	  tmpy = dry * grad_1 + dipole.p_indy[j] * grad_2;
	  tmpz = drz * grad_1 + dipole.p_indz[j] * grad_2;
	  atom.fx[i] -= tmpx; atom.fy[i] -= tmpy; atom.fz[i] -= tmpz;
	  atom.fx[j] += tmpx; atom.fy[j] += tmpy; atom.fz[j] += tmpz;
	  // stress
	  tmpx *= rr; tmpy *= rr; tmpz *= rr; // YU 2013.10.18
	  atom.satom[i][0][0] -= tmpx*drx/2.0;
	  atom.satom[i][0][1] -= tmpx*dry/2.0;
	  atom.satom[i][0][2] -= tmpx*drz/2.0;
	  atom.satom[i][1][1] -= tmpy*dry/2.0;
	  atom.satom[i][1][2] -= tmpy*drz/2.0;
	  atom.satom[i][2][2] -= tmpz*drz/2.0;
	  atom.satom[j][0][0] -= tmpx*drx/2.0;
	  atom.satom[j][0][1] -= tmpx*dry/2.0;
	  atom.satom[j][0][2] -= tmpx*drz/2.0;
	  atom.satom[j][1][1] -= tmpy*dry/2.0;
	  atom.satom[j][1][2] -= tmpy*drz/2.0;
	  atom.satom[j][2][2] -= tmpz*drz/2.0;
	  // dipole-monopole contributions //if (dp_alpha[typ1] && charge[typ2])
	  if (ipair!=0) {
	    value_sum = value_tail + srval;
	    grad_sum  = grad_tail + srgrad;
	  } else {
	    value_sum = value_tail;
	    grad_sum  = grad_tail;
	  }
	  rp_i = dipole.p_indx[i]*drx + dipole.p_indy[i]*dry + dipole.p_indz[i]*drz;
	  value  = dipole.charge[typ2] * rp_i * rr * value_sum;
	  //grad_1 = dipole.charge[typ2] * rp_i * rr2 * grad_sum ;
	  grad_1 = dipole.charge[typ2] * rp_i * rr * grad_sum ; //YU 2013.10.03
	  grad_2 = dipole.charge[typ2] * value_sum;
	  atom.epot[i] += value/2.0;  atom.epot[j] += value/2.0;
	  tmpx = drx * grad_1 + dipole.p_indx[i] * grad_2;
	  tmpy = dry * grad_1 + dipole.p_indy[i] * grad_2;
	  tmpz = drz * grad_1 + dipole.p_indz[i] * grad_2;
	  atom.fx[i] += tmpx; atom.fy[i] += tmpy; atom.fz[i] += tmpz;
	  atom.fx[j] -= tmpx; atom.fy[j] -= tmpy; atom.fz[j] -= tmpz;
	  // stress
	  tmpx *= rr; tmpy *= rr; tmpz *= rr; // YU 2013.10.18
	  atom.satom[i][0][0] += tmpx*drx/2.0;
	  atom.satom[i][0][1] += tmpx*dry/2.0;
	  atom.satom[i][0][2] += tmpx*drz/2.0;
	  atom.satom[i][1][1] += tmpy*dry/2.0;
	  atom.satom[i][1][2] += tmpy*drz/2.0;
	  atom.satom[i][2][2] += tmpz*drz/2.0;
	  atom.satom[j][0][0] += tmpx*drx/2.0;
	  atom.satom[j][0][1] += tmpx*dry/2.0;
	  atom.satom[j][0][2] += tmpx*drz/2.0;
	  atom.satom[j][1][1] += tmpy*dry/2.0;
	  atom.satom[j][1][2] += tmpy*drz/2.0;
	  atom.satom[j][2][2] += tmpz*drz/2.0;
	  // dipole-dipole contributions //if (dp_alpha[typ1] && dp_alpha[typ2])
	  if ((dipole.dp_alpha[typ1]!=0)&&(dipole.dp_alpha[typ2]!=0)) {
	  pp_ij = dipole.p_indx[i]*dipole.p_indx[j]
	    + dipole.p_indy[i]*dipole.p_indy[j] + dipole.p_indz[i]*dipole.p_indz[j];
	  tmp_1 = 3 * rp_i * rp_j;
	  tmp_2 = 3 * value_tail / rr2;
	  value = -(tmp_1 - pp_ij) * value_tail;
	  grad_1 = (tmp_1 - pp_ij) * grad_tail;
	  grad_2 = 2 * rp_i * rp_j;
	  atom.epot[i] += value/2.0; atom.epot[j] += value/2.0;
	  //tmpx = grad_1 * rr * drx - tmp_2 * rr *
	  tmpx = grad_1 * drx - tmp_2 * rr * // YU
	    (grad_2 * drx - rp_i * dipole.p_indx[j] - rp_j * dipole.p_indx[i]);
	  //tmpy = grad_1 * rr * dry - tmp_2 * rr *
	  tmpy = grad_1 * dry - tmp_2 * rr * // YU
	    (grad_2 * dry - rp_i * dipole.p_indy[j] - rp_j * dipole.p_indy[i]);
	  //tmpz = grad_1 * rr * drz - tmp_2 * rr *
	  tmpz = grad_1 * drz - tmp_2 * rr * // YU
	    (grad_2 * drz - rp_i * dipole.p_indz[j] - rp_j * dipole.p_indz[i]);
	  atom.fx[i] -= tmpx; atom.fx[j] += tmpx;
	  atom.fy[i] -= tmpy; atom.fy[j] += tmpy;
	  atom.fz[i] -= tmpz; atom.fz[j] += tmpz;
	  // stress
	  tmpx *= rr; tmpy *= rr; tmpz *= rr; // YU 2013.10.18
	  atom.satom[i][0][0] -= tmpx*drx/2.0;
	  atom.satom[i][0][1] -= tmpx*dry/2.0;
	  atom.satom[i][0][2] -= tmpx*drz/2.0;
	  atom.satom[i][1][1] -= tmpy*dry/2.0;
	  atom.satom[i][1][2] -= tmpy*drz/2.0;
	  atom.satom[i][2][2] -= tmpz*drz/2.0;
	  atom.satom[j][0][0] -= tmpx*drx/2.0;
	  atom.satom[j][0][1] -= tmpx*dry/2.0;
	  atom.satom[j][0][2] -= tmpx*drz/2.0;
	  atom.satom[j][1][1] -= tmpy*dry/2.0;
	  atom.satom[j][1][2] -= tmpy*drz/2.0;
	  atom.satom[j][2][2] -= tmpz*drz/2.0;
	  }
	}
      } }
  } // FOURTH LOOP END

  // FIFTH LOOP: self energy contributions and sum-up force contributions
  double qq, pp;
  for (int i=1; i<=atom.natom; i++) {
    //typ1 = dptyp(atom.asp[i]);
    //typ1 = dptyp(i);
    typ1 = dipole.typen[i];
    if (typ1 <= 0) continue;
    // self energy contributions
    qq = dipole.charge[typ1]*dipole.charge[typ1];
    value = dipole.dp_eps * dipole.dp_kappa * qq / sqrt(M_PI);
    atom.epot[i] -= value;
    if (dipole.dp_alpha[typ1] != 0.0) {
      pp = dipole.p_indx[i]*dipole.p_indx[i]
	+ dipole.p_indy[i]*dipole.p_indy[i] + dipole.p_indz[i]*dipole.p_indz[i];
      //printf("%d %e %e\n",i,pp,dipole.p_indx[i]);
      value = pp / (2 * dipole.dp_alpha[typ1]);
      atom.epot[i] += value;
    }
  } // FIFTH LOOP END

  //atom.epotsum=0.0;
  //for (int i=1; i<=atom.natom; i++) {
  //  atom.epotsum += atom.epot[i]*eV;
  //}
  //printf("%e   %e %e %e\n",atom.epotsum/atom.natom,atom.fx[1],atom.fy[1],atom.fz[1]);

 OUT:
  atom.epotsum=0.0;
  for (int i=1; i<=atom.natom; i++) {
    typ1 = dipole.typen[i];
    if (typ1 <= 0) continue;
    atom.epot[i] *= eV;
    atom.epotsum += atom.epot[i];
    atom.fx[i] *= eV/ang; atom.fy[i] *= eV/ang; atom.fz[i] *= eV/ang;
    //printf("%10d %20.12e %20.12e %20.12e\n",i-1,atom.fx[i]/eV*ang,atom.fy[i]/eV*ang,atom.fz[i]/eV*ang);
  }
  //int i=1;
  //printf("%10d %20.12e %20.12e %20.12e %20.12e\n",i-1,atom.epotsum/eV/atom.natom,atom.fx[i]/eV*ang,atom.fy[i]/eV*ang,atom.fz[i]/eV*ang);
    
  //Atomic stress
  for (int i=1; i<=atom.natom; i++) {
    typ1 = dipole.typen[i];
    if (typ1 <= 0) continue;
    atom.satom[i][1][0] = atom.satom[i][0][1];
    atom.satom[i][2][0] = atom.satom[i][0][2];
    atom.satom[i][2][1] = atom.satom[i][1][2]; }
  resetmat(cell.dmat);
  for (int i=1; i<=atom.natom; i++) {
    typ1 = dipole.typen[i];
    if (typ1 == 0) continue;
    for (int j=0; j<3; j++) {
      for (int k=0; k<3; k++) {
	atom.satom[i][j][k] *= eV;
	cell.dmat[j][k] -= atom.satom[i][j][k]; } } }
  cell.virx = cell.dmat[0][0];
  cell.viry = cell.dmat[1][1];
  cell.virz = cell.dmat[2][2];
  cell.volume = cell.Getvolume();
  for (int i=1; i<=atom.natom; i++) {
    typ1 = dipole.typen[i];
    if (typ1 <= 0) continue;
    for (int j=0; j<3; j++) {
      for (int k=0; k<3; k++) {
	//atom.satom[i][j][k] *= (double)atom.natom / cell.volume; } } } // old version without combined mode
	atom.satom[i][j][k] *= (double)dipole.natoms / cell.volume; } } } // new version with combined mode

} // end of e_force_dipole

void dipole_delete()
{
  if (dipole.r_cut)    { delete[] dipole.r_cut;    dipole.r_cut    = NULL; }
  if (dipole.dp_alpha) { delete[] dipole.dp_alpha; dipole.dp_alpha = NULL; }
  if (dipole.dp_b)     { delete[] dipole.dp_b;     dipole.dp_b     = NULL; }
  if (dipole.dp_c)     { delete[] dipole.dp_c;     dipole.dp_c     = NULL; }
  if (dipole.E_statx)  { delete[] dipole.E_statx;  dipole.E_statx  = NULL; }
  if (dipole.E_staty)  { delete[] dipole.E_staty;  dipole.E_staty  = NULL; }
  if (dipole.E_statz)  { delete[] dipole.E_statz;  dipole.E_statz  = NULL; }
  if (dipole.E_indx)   { delete[] dipole.E_indx;   dipole.E_indx   = NULL; }
  if (dipole.E_indy)   { delete[] dipole.E_indy;   dipole.E_indy   = NULL; }
  if (dipole.E_indz)   { delete[] dipole.E_indz;   dipole.E_indz   = NULL; }
  if (dipole.E_oldx)   { delete[] dipole.E_oldx;   dipole.E_oldx   = NULL; }
  if (dipole.E_oldy)   { delete[] dipole.E_oldy;   dipole.E_oldy   = NULL; }
  if (dipole.E_oldz)   { delete[] dipole.E_oldz;   dipole.E_oldz   = NULL; }
  if (dipole.E_totx)   { delete[] dipole.E_totx;   dipole.E_totx   = NULL; }
  if (dipole.E_toty)   { delete[] dipole.E_toty;   dipole.E_toty   = NULL; }
  if (dipole.E_totz)   { delete[] dipole.E_totz;   dipole.E_totz   = NULL; }
  if (dipole.p_srx)    { delete[] dipole.p_srx;    dipole.p_srx    = NULL; }
  if (dipole.p_sry)    { delete[] dipole.p_sry;    dipole.p_sry    = NULL; }
  if (dipole.p_srz)    { delete[] dipole.p_srz;    dipole.p_srz    = NULL; }
  if (dipole.p_indx)   { delete[] dipole.p_indx;   dipole.p_indx   = NULL; }
  if (dipole.p_indy)   { delete[] dipole.p_indy;   dipole.p_indy   = NULL; }
  if (dipole.p_indz)   { delete[] dipole.p_indz;   dipole.p_indz   = NULL; }
  if (dipole.charge)   { delete[] dipole.charge;   dipole.charge   = NULL; }
  if (dipole.buck_a)   { delete[] dipole.buck_a;   dipole.buck_a   = NULL; }
  if (dipole.buck_s)   { delete[] dipole.buck_s;   dipole.buck_s   = NULL; }
  if (dipole.buck_c)   { delete[] dipole.buck_c;   dipole.buck_c   = NULL; }
  if (dipole.buck_vrc) { delete[] dipole.buck_vrc; dipole.buck_vrc = NULL; }
  if (dipole.buck_vprc){ delete[] dipole.buck_vprc;dipole.buck_vprc= NULL; }
  if (dipole.type)     { delete[] dipole.type;     dipole.type     = NULL; }
  if (dipole.typen)    { delete[] dipole.typen;    dipole.typen    = NULL; }
  if (dipole.ptype) {
    for (int i=0; i<dipole.ntype+1; i++) { delete[] dipole.ptype[i]; }
    delete[] dipole.ptype; dipole.ptype = NULL; }
}
void dipole_alloc()
{
  dipole.r_cut    = new double[dipole.nptype+1];
  dipole.dp_alpha = new double[dipole.ntype+1];
  dipole.dp_b     = new double[dipole.nptype+1];
  dipole.dp_c     = new double[dipole.nptype+1];
  dipole.E_statx  = new double[atom.natom+1];
  dipole.E_staty  = new double[atom.natom+1];
  dipole.E_statz  = new double[atom.natom+1];
  dipole.E_indx   = new double[atom.natom+1];
  dipole.E_indy   = new double[atom.natom+1];
  dipole.E_indz   = new double[atom.natom+1];
  dipole.E_oldx   = new double[atom.natom+1];
  dipole.E_oldy   = new double[atom.natom+1];
  dipole.E_oldz   = new double[atom.natom+1];
  dipole.E_totx   = new double[atom.natom+1];
  dipole.E_toty   = new double[atom.natom+1];
  dipole.E_totz   = new double[atom.natom+1];
  dipole.p_srx    = new double[atom.natom+1];
  dipole.p_sry    = new double[atom.natom+1];
  dipole.p_srz    = new double[atom.natom+1];
  dipole.p_indx   = new double[atom.natom+1];
  dipole.p_indy   = new double[atom.natom+1];
  dipole.p_indz   = new double[atom.natom+1];
  dipole.charge   = new double[dipole.ntype+1];
  dipole.buck_a   = new double[dipole.nptype+1];
  dipole.buck_s   = new double[dipole.nptype+1];
  dipole.buck_c   = new double[dipole.nptype+1];
  dipole.buck_vrc = new double[dipole.nptype+1];
  dipole.buck_vprc= new double[dipole.nptype+1];
  dipole.type     = new int[dipole.ntype+1];
  dipole.ptype    = new int*[dipole.ntype+1];
  for (int i=0; i<dipole.ntype+1; i++) {
    dipole.ptype[i] = new int[dipole.ntype+1]; }
  dipole.typen    = new int[atom.natom+1];
}

int dptyp(char* at)
{
  int i = atom_number(at);
  for (int j=1; j<=dipole.ntype; j++) {
    if (i == dipole.type[j]) { return j; }
  }
  return -1;
}
int dptyp(int at)
{
  //int i = atom.anum[at];
  int i = at;
  for (int j=1; j<=dipole.ntype; j++) {
    if (i == dipole.type[j]) { return j; }
  }
  return -1;
}
int pairtyp(char* ati, char* atj)
{
  int i = dptyp(ati);
  int j = dptyp(atj);
  if (((i==1)&&(j==2))||((i==2)&&(j==1))) { return 1;
  } else if ((i==2)&&(j==2)) { return 2;
  } else if (((i==2)&&(j==3))||((i==3)&&(j==2))) { return 3;
  } else { return 0; }
}
int pairtyp(int ati, int atj)
{
  int i = dptyp(ati);
  int j = dptyp(atj);
  if (((i==1)&&(j==2))||((i==2)&&(j==1))) { return 1;
  } else if ((i==2)&&(j==2)) { return 2;
  } else if (((i==2)&&(j==3))||((i==3)&&(j==2))) { return 3;
  } else { return 0; }
}
int dipole_param()
{ //Zr=1, O=2, Y=3,  Zr-O=1, O-O=2, Y-O=3
  dipole.type[1] = 40;
  dipole.type[2] =  8;
  dipole.type[3] = 39;
  dipole.ptype[1][1] = 0;
  dipole.ptype[1][2] = 1;
  dipole.ptype[1][3] = 0;
  dipole.ptype[2][2] = 2;
  dipole.ptype[2][3] = 3;
  dipole.ptype[3][3] = 0;
  dipole.ptype[2][1] = dipole.ptype[1][2];
  dipole.ptype[3][1] = dipole.ptype[1][3];
  dipole.ptype[3][2] = dipole.ptype[2][3];
  dipole.dp_kappa = 0.10;
  dipole.charge[1] = 2.116022;
  dipole.charge[2] =-1.055614;
  dipole.charge[3] = 1.511511;
  dipole.dp_alpha[1] = 0.;
  dipole.dp_alpha[2] = 0.010019;
  dipole.dp_alpha[3] = 0.;
  dipole.dp_b[1] = 13.439888;
  dipole.dp_b[2] = 59.408009;
  dipole.dp_b[3] = 68.516735;
  dipole.dp_c[1] =  -0.170177;
  dipole.dp_c[2] = -28.550774;
  dipole.dp_c[3] =   0.00000;
  dipole.buck_a[1] = 7299.50008961;
  dipole.buck_a[2] = 8175.01587504;
  dipole.buck_a[3] = 6968.55253502;
  dipole.buck_s[1] = 0.22914290;
  dipole.buck_s[2] = 0.28736600;
  dipole.buck_s[3] = 0.23183401;
  dipole.buck_c[1] =  87383.66707018;
  dipole.buck_c[2] = 850244.09929026;
  dipole.buck_c[3] =  4700.69794845;
  dipole.r_cut[1] = 10.0;
  dipole.r_cut[2] = 10.0;
  dipole.r_cut[3] = 10.0;

  dipole.dp_cut = 14.0;
  if (dipole.shift) {
    double value = 0.0, grad = 0.0;
    for (int i=1; i<=3; i++) {
      buckingham(dipole.dp_cut, value, grad, dipole.buck_a[i], dipole.buck_s[i],
		 dipole.buck_c[i]);
      dipole.buck_vrc[i]  = value;
      dipole.buck_vprc[i] = grad;
      //printf("Buckingham shift %d %e %e\n",i,dipole.buck_vrc[i],dipole.buck_vprc[i]);
    }
  }
  for (int i=1; i<=atom.natom; i++) {
    dipole.typen[i] = dptyp(atom.asp[i]);
    if (dipole.typen[i] == 0) { return 0; }
  }
  return 1;
}

void bmh(double r, double &e, double &f, double a, double sig, double c)
{
  //  double rang = r;
  double exprs = exp((a-r)/sig);
  e = c * sig * exprs;
  f = -c * exprs;
}

void buckingham(double r, double &e, double &f, double a, double sig, double c)
{
  //  double rang = r;
  double x = sig/r; double x2 = x * x; double x3 = x2 * x;
  double x6 = x3 * x3; double x5 = x3 * x2;
  double exprs = exp(-r/sig);
  e = a * exprs - c * x6;
  f = -1.0/sig * a * exprs + c * 6.0 * x6 / r;
}
void bucksft(double r, double &e, double &f, double a, double sig, double c,
	     double rc, double vrc, double vprc)
{
  //  double rang = r;
  double x = sig/r; double x2 = x * x; double x3 = x2 * x;
  double x6 = x3 * x3; double x5 = x3 * x2;
  double exprs = exp(-r/sig);
  e = a * exprs - c * x6;
  f = -1.0/sig * a * exprs + c * 6.0 * x6 / r;
  // The following two lines are weird, but identical with potfit. (YU 2014.1.9)
  // We regard it as a new function (e.g. Modified Buckingam)
  // e += -vrc - r*(r-rc)*vprc;
  // f += (rc-2*r)*vprc;
  // We decided to go with correct shift below instead of the above lines (YU 2014.6.2)
  e += -vrc - (r-rc)*vprc;
  f += -vprc;

  return;

  // 'artificial' smoothing: Needed for MD simulation to avoid 'jumps' in stress and force
  if (rc-r < 0) { return; }
  if (rc-r < MB_SMOOTH_RANGE) {
    double xx = (r - rc + MB_SMOOTH_RANGE) / MB_SMOOTH_RANGE * M_PI / 2;
    f = f*cos(xx) - M_PI/2/MB_SMOOTH_RANGE*e*sin(xx);
    e *= cos(xx);
  }
  // The following should be 'correct shift', but unused. (YU 2014.1.9)
  // e += -vrc - (r-rc)*vprc;
  // f += -vprc;
}

double shortrange_value(double r, double a, double b, double c)
{
  static double x[5];

  x[0] = b * r;
  x[1] = x[0] * x[0];
  x[2] = x[1] * x[0];
  x[3] = x[1] * x[1];
  x[4] = 1 + x[0] + x[1] / 2 + x[2] / 6 + x[3] / 24;

  return a * c * x[4] * exp(-x[0]) / dipole.dp_eps;
}
void shortrange_term(double r, double b, double c, double &srval_tail,
  double &srgrad_tail)
{
  static double x[6];

  x[0] = b * r;
  x[1] = x[0] * x[0];
  x[2] = x[1] * x[0];
  x[3] = x[1] * x[1];
  x[4] = 1 + x[0] + x[1] / 2 + x[2] / 6 + x[3] / 24;
  x[5] = exp(-x[0]);

  srval_tail = c * x[4] * x[5] / dipole.dp_eps;
  //srgrad_tail = -c * b * x[3] * x[5] / (24 * dipole.dp_eps * r); // YU 2013.09.18
  srgrad_tail = -c * b * x[3] * x[5] / (24 * dipole.dp_eps);
}

void elstat_shift(double r, double dp_kappa, double &value_tail, double &grad_tail,
  double &ggrad_tail)
{
  static double ftail, gtail, ggtail, ftail_cut, gtail_cut, ggtail_cut;
  static double x[3];

  x[0] = r * r;
  x[1] = dipole.dp_cut * dipole.dp_cut;
  x[2] = x[0] - x[1];

  elstat_value(r, dp_kappa, ftail, gtail, ggtail);
  elstat_value(dipole.dp_cut, dp_kappa, ftail_cut, gtail_cut, ggtail_cut);

  value_tail = ftail - ftail_cut - x[2] * gtail_cut / 2;
  grad_tail = gtail - gtail_cut;
  ggrad_tail = 0.;

  value_tail -= x[2] * x[2] * ggtail_cut / 8;
  grad_tail -= x[2] * ggtail_cut / 2;
  ggrad_tail = ggtail - ggtail_cut;

  /* YU: why below lines do not work?
  ggrad_tail = ggtail - ggtail_cut;
  grad_tail  = gtail - gtail_cut - (r-dipole.dp_cut) * ggtail_cut;
  value_tail = ftail - ftail_cut - (r-dipole.dp_cut) * gtail_cut - x[2] * gtail_cut / 2;
  */
}

void elstat_value(double r, double dp_kappa, double &ftail, double &gtail, double &ggtail)
{
  static double x[4];

  x[0] = r * r;
  x[1] = dp_kappa * dp_kappa;
  x[2] = 2 * dipole.dp_eps * dp_kappa / sqrt(M_PI);
  x[3] = exp(-x[0] * x[1]);


  ftail = dipole.dp_eps * erfc(dp_kappa * r) / r;
  gtail = -(ftail + x[2] * x[3]) / x[0];
  ggtail = (2 * x[1] * x[2] * x[3] - gtail * 3) / x[0];

  /*
  ggtail = (ftail/x[0] + x[2]*x[3]/x[0] +x[1]*x[2]*x[3]/2.0/r)/r*2.0;
  */
}
double dsquare(double d)
{
  return d * d;
}

#if !defined __linux__ && !defined __APPLE__
static const double rel_error= 1E-12;
double erf(double x)
{
    static const double two_sqrtpi=  1.128379167095512574;        // 2/sqrt(pi)
    if (fabs(x) > 2.2) {
        return 1.0 - erfc(x);        //use continued fraction when fabs(x) > 2.2
    }
    double sum= x, term= x, xsqr= x*x;
    int j= 1;
    do {
        term*= xsqr/j;
        sum-= term/(2*j+1);
        ++j;
        term*= xsqr/j;
        sum+= term/(2*j+1);
        ++j;
    } while (fabs(term)/sum > rel_error);
    return two_sqrtpi*sum;
}
double erfc(double x)
{
    static const double one_sqrtpi=  0.564189583547756287;        // 1/sqrt(pi)
    if (fabs(x) < 2.2) {
        return 1.0 - erf(x);        //use series when fabs(x) < 2.2
    }
    if (signbit(x)) {               //continued fraction only valid for x>0
        return 2.0 - erfc(-x);
    }
    double a=1, b=x;                //last two convergent numerators
    double c=x, d=x*x+0.5;          //last two convergent denominators
    double q1,q2;                   //last two convergents (a/c and b/d)
    double n= 1.0, t;
    do {
        t= a*n+b*x;
        a= b;
         b= t;
        t= c*n+d*x;
        c= d;
        d= t;
        n+= 0.5;
        q1= q2;
        q2= b/d;
      } while (fabs(q1-q2)/q2 > rel_error);
    return one_sqrtpi*exp(-x*x)*q2;
}
#endif

int dipole_read_number(const char* fname)
{
  std::ifstream fin( fname );
  std::string line, arg1;
  while (getline(fin, line)) { remove_carriage_return(line);
    if (line.length()==0) { continue; }
    remove_head_spaces(line);
    if (line.at(0)!='#') { break; }
  }
  remove_after_sharp(line);
  int narg = count_arg_number(line); dipole.ntype = narg;
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
    printf("Dipole: Potential file error\n (incorrect flag of pair interaction)\n");
    printf(" narg = %d, np = %d\n",narg,np);
    printf("########################################\n");
    fin.close(); fin.clear();
    return 0;
  }
  dipole.nptype = 0;
  for (int i=0; i<np; i++) {
    get_first_arg(line,arg1);
    if (strcmp(arg1.c_str(), "1") == 0) { dipole.nptype++;
    } else if (strcmp(arg1.c_str(), "Y") == 0) { dipole.nptype++;
    } else if (strcmp(arg1.c_str(), "y") == 0) { dipole.nptype++;
    } else if (strcmp(arg1.c_str(), "T") == 0) { dipole.nptype++;
    } else if (strcmp(arg1.c_str(), "t") == 0) { dipole.nptype++;
    }
  }
  fin.close(); fin.clear();
  return 1;
}
int dipole_read_param(const char* fname)
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
  for (int i=1; i<=dipole.ntype; i++) {
    get_first_arg(line,arg1);
    char sp[3] = "aa"; strcpy(sp,arg1.c_str());
    dipole.type[i] = atom_number(sp); //printf("%d %d\n",i,dipole.type[i]);
  }
  // Pair-interaction flags
  while (getline(fin, line)) { remove_carriage_return(line);
    if (line.length()==0) { continue; }
    remove_head_spaces(line);
    if (line.at(0)!='#') { remove_after_sharp(line); break; }
  }
  int narg = dipole.ntype; 
  int pt = 1;
  for (int i=1; i<=narg; i++) {
    for (int j=i; j<=narg; j++) {
      get_first_arg(line,arg1);
      if (strcmp(arg1.c_str(), "1") == 0) { dipole.ptype[i][j] = pt;
      } else if (strcmp(arg1.c_str(), "Y") == 0) { dipole.ptype[i][j] = pt;
      } else if (strcmp(arg1.c_str(), "y") == 0) { dipole.ptype[i][j] = pt;
      } else if (strcmp(arg1.c_str(), "T") == 0) { dipole.ptype[i][j] = pt;
      } else if (strcmp(arg1.c_str(), "t") == 0) { dipole.ptype[i][j] = pt;
      } else { dipole.ptype[i][j] = 0; }
      //printf("%d %d %d\n",i,j,pt);
      if (i != j) { dipole.ptype[j][i] = dipole.ptype[i][j]; }
      if (dipole.ptype[i][j] > 0) { pt++; }
    }
  }
  pt--;
  // read and set parameters
  double xx;
  // charge
  dipole_find_param(fname,"charge",line);
  for (int i=1; i<=dipole.ntype; i++) {
    get_first_arg(line, arg1);
    strcpy (larg1, arg1.c_str());
    xx = atof(larg1);
    dipole.charge[i] = xx; //printf("charge %d %f\n",i,xx);
  }
  //buck_a, buck_sigma, buck_c
  dipole_find_param(fname,"buck_a",line);
  for (int i=1; i<=pt; i++) {
    get_first_arg(line, arg1);
    strcpy (larg1, arg1.c_str());
    xx = atof(larg1);
    dipole.buck_a[i] = xx; //printf("buck_a %d %f\n",i,xx);
  }
  dipole_find_param(fname,"buck_sigma",line);
  for (int i=1; i<=pt; i++) {
    get_first_arg(line, arg1);
    strcpy (larg1, arg1.c_str());
    xx = atof(larg1);
    dipole.buck_s[i] = xx; //printf("buck_s %d %f\n",i,xx);
  }
  dipole_find_param(fname,"buck_c",line);
  for (int i=1; i<=pt; i++) {
    get_first_arg(line, arg1);
    strcpy (larg1, arg1.c_str());
    xx = atof(larg1);
    dipole.buck_c[i] = xx; //printf("buck_c %d %f\n",i,xx);
  }
  // ew_rcut, ew_kappa, r_cut
  dipole_find_param(fname,"ew_rcut",line);
  get_first_arg(line, arg1);
  strcpy (larg1, arg1.c_str());
  xx = atof(larg1);
  dipole.dp_cut = xx; //printf("ew_cut %f\n",xx);
  //dipole.ew_rcut = xx; //printf("ew_cut %f\n",xx);
  dipole_find_param(fname,"ew_kappa",line);
  get_first_arg(line, arg1);
  strcpy (larg1, arg1.c_str());
  xx = atof(larg1);
  dipole.dp_kappa = xx; //printf("dp_kappa %f\n",xx);
  dipole_find_param(fname,"r_cut",line);
  for (int i=1; i<=pt; i++) {
    get_first_arg(line, arg1);
    strcpy (larg1, arg1.c_str());
    xx = atof(larg1);
    dipole.r_cut[i] = xx; //printf("r_cut %d %f\n",i,xx);
  }
  // dp_alpha, dp_b, dp_c
  dipole_find_param(fname,"dp_alpha",line);
  for (int i=1; i<=dipole.ntype; i++) {
    get_first_arg(line, arg1);
    strcpy (larg1, arg1.c_str());
    xx = atof(larg1);
    dipole.dp_alpha[i] = xx; //printf("dp_alpha %d %f\n",i,xx);
  }
  dipole_find_param(fname,"dp_b",line);
  for (int i=1; i<=pt; i++) {
    get_first_arg(line, arg1);
    strcpy (larg1, arg1.c_str());
    xx = atof(larg1);
    dipole.dp_b[i] = xx; //printf("dp_b %d %f\n",i,xx);
  }
  dipole_find_param(fname,"dp_c",line);
  for (int i=1; i<=pt; i++) {
    get_first_arg(line, arg1);
    strcpy (larg1, arg1.c_str());
    xx = atof(larg1);
    dipole.dp_c[i] = xx; //printf("dp_c %d %f\n",i,xx);
  }
  //if (dipole_find_param(fname,"dp_cut",line)) {
  //  get_first_arg(line, arg1);
  //  strcpy (larg1, arg1.c_str());
  //  xx = atof(larg1);
  //  dipole.dp_cut = xx; //printf("dp_cut %f\n",xx);
  //}
  if (dipole.shift) {
    double value = 0.0, grad = 0.0;
    for (int i=1; i<=pt; i++) {
      buckingham(dipole.dp_cut, value, grad, dipole.buck_a[i], dipole.buck_s[i],
		 dipole.buck_c[i]);
      dipole.buck_vrc[i]  = value;
      dipole.buck_vprc[i] = grad;
      //printf("Buckingham shift %d %e %e\n",i,dipole.buck_vrc[i],dipole.buck_vprc[i]);
    }
  }

  //check if all species are mentioned in the potential with dipole interactions
  //and calculate number of atoms for dipole interaction
  dipole.natoms = 0;

  for (int i=1; i<=atom.natom; i++) {
    dipole.typen[i] = dptyp(atom.asp[i]);
    //if (dipole.typen[i] == 0) { return 0; }
    //if (dipole.typen[i] != 0) { dipole.natoms++; }
    if (dipole.typen[i] > 0) { dipole.natoms++; } // YU20140612
  }

  // return zero if there are some atoms that are not handled in dipole potential
  if (dipole.natoms != atom.natom) return 0;

  return 1;
}

int dipole_find_param(const char* fname, const char* tag, std::string &line)
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
