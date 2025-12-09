#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "myheader.h"

void cnt_wall_setnormalvectors(); void cnt_wall_centering();
extern int mode_cnt_corrugation, cnt_load_algo;
extern float cnt_ring_radius, cnt_ring_fmax, cnt_ring_sharpness;
extern float cnt_pressure, cnt_pressure_ftot, cnt_pressure_gpa, cnt_pressure_gpa2;
extern float outermost_radius_f, cylinder_side_area_f, cylinder_side_area0_f;
extern float load_aux_ftot;
void loading();
void loading_old();
double repulsion(double r, double fmax, double sharpness);
void calc_center_xy(double &x, double &y);
void calc_cylinder_side_area(double center_x, double center_y, double &rr, double &area);
void calc_cylinder_side_area_arb(double &area);
void calc_cylinder_side_area_arb2(double &area);
void calc_cylinder_side_area_arb_slice(double &area);
double dsquare(double d);
double theta(double x, double y);
void sort(int n, double *val, int *nod);
void potential();

void loading()
{
  for (int i=1; i<=atom.natom; i++) { // Reset F by loading
    atom.fx_l[i] = 0; atom.fy_l[i] = 0; atom.fz_l[i] = 0;
  }

  if (!mode_cnt_corrugation) { 
    if (ifpush1 == 1) {
      for (int i = 1; i <= atom.natom; i++) {
        if (atom.group[i] == push1num) {
          atom.fx_l[i] += pushfx1 * 1.0e-9;
          atom.fy_l[i] += pushfy1 * 1.0e-9;
          atom.fz_l[i] += pushfz1 * 1.0e-9;
        }
      }
    }
    if (ifpush2 == 1) {
      for (int i = 1; i <= atom.natom; i++) {
        if (atom.group[i] == push2num) {
          atom.fx_l[i] += pushfx2 * 1.0e-9;
          atom.fy_l[i] += pushfy2 * 1.0e-9;
          atom.fz_l[i] += pushfz2 * 1.0e-9;
        }
      }
    }
    for (int i=1; i<=atom.natom; i++) { // Add F by loading
	    atom.fx[i] += atom.fx_l[i]; atom.fy[i] += atom.fy_l[i]; atom.fz[i] += atom.fz_l[i];
    }
  }

  if (mode_cnt_corrugation) { 
    double center_x, center_y, xx, yy, rr, ff, rr0, x0, y0;
    /*
    center_x = cell.hmat[0][0]/2.0;
    center_y = cell.hmat[1][1]/2.0;
    rr = 0;
    for (int i=1; i<=atom.natom; i++) {
      xx = atom.rx[i] - center_x;
      yy = atom.ry[i] - center_y;
      rr += sqrt(xx*xx+yy*yy);
    }
    rr /= (double)atom.natom;
    double cylinder_side_area = 2*M_PI*rr*cell.hmat[2][2]; //(in m^2)
    */
    double outermost_radius, cylinder_side_area;
    double cylinder_side_area0;
    calc_center_xy(center_x, center_y);
    // Side area calculation for cylinder shape
    //calc_cylinder_side_area(center_x, center_y, outermost_radius, cylinder_side_area);
    //printf("area %e\n",cylinder_side_area);
    calc_cylinder_side_area_arb2(cylinder_side_area); // Side area calculation for arbitrary shape
    calc_cylinder_side_area_arb(cylinder_side_area0); // Side area calculation (ver.1)
    //printf("area %e\n",cylinder_side_area);
    //outermost_radius_f = (float)outermost_radius/ang;
    cylinder_side_area_f = (float)cylinder_side_area /ang/ang;
    cylinder_side_area0_f= (float)cylinder_side_area0/ang/ang;
    if (cnt_load_algo == 0) { // Hydrostatic pressure (Wall) mode
      double factor0 = 1.0e-10*cnt_pressure;
      double area0 = 0.8731e-20;
      cnt_wall_setnormalvectors();
      //cnt_pressure_ftot = 0.0;
      for (int i=1; i<=atom.nwall; i++) {
	int i0 = atom.wall_v[i][0]; int i1 = atom.wall_v[i][1]; int i2 = atom.wall_v[i][2];
	i0 = abs(i0); i1 = abs(i1); i2 = abs(i2);
	double factor = factor0*atom.wall_area[i]/area0;
	//factor = factor0; // <== If commented out, identical with old version
	double fxx = atom.wall_nvec[i][0]*factor;
	double fyy = atom.wall_nvec[i][1]*factor;
	double fzz = atom.wall_nvec[i][2]*factor;
	//cnt_pressure_ftot = cnt_pressure_ftot + sqrt(fxx*fxx+fyy*fyy+fzz*fzz); // Erroneous
	atom.fx_l[i0] += fxx; atom.fy_l[i0] += fyy; atom.fz_l[i0] += fzz;
	atom.fx_l[i1] += fxx; atom.fy_l[i1] += fyy; atom.fz_l[i1] += fzz;
	atom.fx_l[i2] += fxx; atom.fy_l[i2] += fyy; atom.fz_l[i2] += fzz;
      }
      //cnt_pressure_gpa = cnt_pressure_ftot/cylinder_side_area * 1e-9; // pressure (in GPa)
      //cnt_pressure_ftot *= 1e9; // total force (in nN)
      //printf("CNT pressure: Force per atom = %f (nN)\n",cnt_pressure_ftot);
    } else if (cnt_load_algo == 1) { // Ring-force mode
      double fmax, sharpness;
      //cnt_pressure_ftot = 0.0;
      for (int i=1; i<=atom.natom; i++) {
	center_x = cell.hmat[0][0]/2.0;
	center_y = cell.hmat[1][1]/2.0;
	xx = atom.rx[i] - center_x;
	yy = atom.ry[i] - center_y;
	rr = sqrt(xx*xx+yy*yy);
	rr0 = (double)cnt_ring_radius - rr*1.0e10;
	fmax = (double)cnt_ring_fmax;
	sharpness = (double)cnt_ring_sharpness;
	ff = repulsion(rr0, fmax, sharpness)*eV/ang;
	//cnt_pressure_ftot = cnt_pressure_ftot + sqrt(ff*ff);
	atom.fx_l[i] += -ff * xx / rr;
	atom.fy_l[i] += -ff * yy / rr;
      }
      //cnt_pressure_gpa = cnt_pressure_ftot/cylinder_side_area * 1e-9; // pressure (in GPa)
      //cnt_pressure_ftot *= 1e9; // total force (in nN)
    } else if (cnt_load_algo == 2) { // Ring-rigid mode
      cnt_pressure_ftot = 0.0;
      for (int i=1; i<=atom.natom; i++) {
	atom.fx_l[i] = atom.fx[i];
	atom.fy_l[i] = atom.fy[i];
	center_x = cell.hmat[0][0]/2.0;
	center_y = cell.hmat[1][1]/2.0;
	xx = atom.rx[i] - center_x;
	yy = atom.ry[i] - center_y;
	rr = sqrt(xx*xx+yy*yy);
	rr0 = (double)cnt_ring_radius - rr*1.0e10;
	if (rr0 < 0) {
	  x0 = atom.rx[i]; y0 = atom.ry[i];
	  atom.rx[i] = center_x + xx / rr * (double)cnt_ring_radius * ang;
	  atom.ry[i] = center_y + yy / rr * (double)cnt_ring_radius * ang;
	  //atom.vx[i] += (atom.rx[i] - x0) / dt;
	  //atom.vy[i] += (atom.ry[i] - y0) / dt;
	  atom.vx[i] = 0; atom.vy[i] = 0;
	}
      }
      potential();
      for (int i=1; i<=atom.natom; i++) {
	xx = atom.rx[i] - center_x;
	yy = atom.ry[i] - center_y;
	rr = sqrt(xx*xx+yy*yy);
	rr0 = (double)cnt_ring_radius - rr*1.0e10;
	if (rr0 <= 1.0e-5) {
	  x0 = atom.fx_l[i]; y0 = atom.fy_l[i];
	  atom.fx_l[i] = -atom.fx[i]; xx = atom.fx_l[i];
	  atom.fy_l[i] = -atom.fy[i]; yy = atom.fy_l[i];
	  //atom.fx[i] -= x0; atom.fy[i] -= y0;
	  cnt_pressure_ftot += sqrt(xx*xx+yy*yy);
	}
      }
      cnt_pressure_gpa = cnt_pressure_ftot/cylinder_side_area * 1e-9; // pressure (in GPa)
      cnt_pressure_gpa2= cnt_pressure_ftot/cylinder_side_area0* 1e-9; // pressure (in GPa)
      cnt_pressure_ftot *= 1e9; // total force (in nN)
    }
    if (cnt_load_algo != 2) {
      // Bug fix 20140413
      cnt_pressure_ftot = 0.0;
      for (int i=1; i<=atom.natom; i++) {
        double fxx = atom.fx_l[i]; double fyy = atom.fy_l[i]; double fzz = atom.fz_l[i]; 
        cnt_pressure_ftot += sqrt(fxx*fxx+fyy*fyy+fzz*fzz);
      }
      cnt_pressure_gpa = cnt_pressure_ftot/cylinder_side_area * 1e-9; // pressure (in GPa)
      cnt_pressure_gpa2= cnt_pressure_ftot/cylinder_side_area0* 1e-9; // pressure (in GPa)
      cnt_pressure_ftot *= 1e9; // total force (in nN)
      //

      for (int i=1; i<=atom.natom; i++) { // Add F by loading
	atom.fx[i] += atom.fx_l[i]; atom.fy[i] += atom.fy_l[i]; atom.fz[i] += atom.fz_l[i];
      }
    }
  } // end of mode_cnt_corrugation

 // auxiliary constraint
  if (yzplane_punch) {
    yzplane_punch_ftot = 0.0;
    double center_x = cell.hmat[0][0]/2.0;
    double center_y = cell.hmat[1][1]/2.0;
    for (int i=1; i<=atom.natom; i++) {
      atom.fx_l[i] = atom.fx[i];
      double xx = atom.rx[i] - center_x;
      double xx0 = (double)yzplane_punch_d/2 - fabs(xx/ang);
      if (xx0 < 0) {
	double x0 = atom.rx[i];
	atom.rx[i] = center_x + (double)yzplane_punch_d/2 * ang * xx / fabs(xx);
	//atom.vx[i] += (atom.rx[i] - x0) / dt;
	atom.vx[i] = 0;
      }
    }
    potential();
    for (int i=1; i<=atom.natom; i++) {
      double xx = atom.rx[i] - center_x;
      double xx0 = (double)yzplane_punch_d/2 - fabs(xx/ang);
      if (xx0 <= 1.0e-5) {
	double x0 = atom.fx_l[i];
	atom.fx_l[i] = -atom.fx[i]; xx = atom.fx_l[i];
	yzplane_punch_ftot += sqrt(xx*xx);
      }
    }
    yzplane_punch_ftot *= 1e9; // total force (in nN)
  }    
  // end of auxiliary constraint
}

double repulsion(double r, double fmax, double sharpness)
{
  double f;
  //double fmax = 10.0;
  //double sharpness = 2.0;
  if (r < 0) {
    f = fmax; 
  } else {
    f = fmax*exp(-r*sharpness);
  }
  return f;
}

void loading_old()
{
 if (mode_cnt_corrugation) {
   double factor = 1.0e-10*cnt_pressure;
   cnt_wall_setnormalvectors();
   cnt_pressure_ftot = 0.0;
   for (int i=1; i<=atom.nwall; i++) {
     int i0 = atom.wall_v[i][0]; int i1 = atom.wall_v[i][1]; int i2 = atom.wall_v[i][2];
     double fxx = atom.wall_nvec[i][0]*factor;
     double fyy = atom.wall_nvec[i][1]*factor;
     double fzz = atom.wall_nvec[i][2]*factor;
     cnt_pressure_ftot = cnt_pressure_ftot + fxx*fxx+fyy*fyy+fzz*fzz;
     atom.fx[i0] = atom.fx[i0] + fxx; atom.fy[i0] = atom.fy[i0] + fyy; atom.fz[i0] = atom.fz[i0] + fzz;
     atom.fx[i1] = atom.fx[i1] + fxx; atom.fy[i1] = atom.fy[i1] + fyy; atom.fz[i1] = atom.fz[i1] + fzz;
     atom.fx[i2] = atom.fx[i2] + fxx; atom.fy[i2] = atom.fy[i2] + fyy; atom.fz[i2] = atom.fz[i2] + fzz;
   }
   cnt_pressure_ftot = sqrt(cnt_pressure_ftot)/(double)atom.natom*1e9; // force per atom (in nN)
   //printf("CNT pressure: Force per atom = %f (nN)\n",cnt_pressure_ftot);
 }
}

void calc_center_xy(double &x, double &y)
{
  x = 0; y = 0;
  for (int i=1; i<=atom.natom; i++) {
    x += atom.rx[i]; y += atom.ry[i]; }
  x /= (double)atom.natom; y /= (double)atom.natom;
}

// Side area calculation assuming the cross section is a right circle
void calc_cylinder_side_area(double center_x, double center_y, double &rr, double &area)
{
  double d = 1.0*ang;
  double xx, yy, rr2, rr2_max;
  rr2_max = 0;
  for (int i=1; i<=atom.natom; i++) {
    xx = atom.rx[i] - center_x;
    yy = atom.ry[i] - center_y;
    rr2 = xx*xx+yy*yy;
    if (rr2_max < rr2 ) { rr2_max = rr2; }
  }
  rr = 0;
  int cnt = 0;
  double rr_thres = sqrt(rr2_max) - d;
  double rr_thres2 = rr_thres * rr_thres;
  for (int i=1; i<=atom.natom; i++) {
    xx = atom.rx[i] - center_x;
    yy = atom.ry[i] - center_y;
    rr2 = xx*xx+yy*yy;
    if (rr2 > rr_thres2) {
      rr += sqrt(rr2); cnt++;
    }
  }
  rr /= (double)cnt; // Calculated radius of outermost CNT (in m)
  //printf(" %d %f \n", cnt, rr/ang);
  area = 2.0*M_PI*rr*cell.hmat[2][2]; //(in m^2)
}

// Side area calculation for an arbitrary cross section shape (ver.1)
// Atoms are projected to x-y plane, outermost points are
// then sieved out, crowded points are erased and the
// shape of the outermost points are found in the form of
// r(theta). If the function becomes two-valued, this does not
// work. Therefore, this is not suitable for the case where
// the CNT cross-section is in a dumbbell shape.
void calc_cylinder_side_area_arb(double &area)
{
  double xx, yy;
  double *rad, *the, *dum; // radius, theta, dummy
  int *nod; // order
  // nod[] (1d array) is used to store order of values
  int num_elem = atom.natom; // # of mesh for r(theta)
  rad = new double[atom.natom]; // r
  the = new double[atom.natom]; // theta
  dum = new double[atom.natom]; // dummy
  nod = new int[atom.natom]; // order
  // calculate center (x,y)
  double center_x, center_y;
  calc_center_xy(center_x, center_y);
  // store r(theta)
  for (int i=1; i<=atom.natom; i++) {
    xx = atom.rx[i] - center_x;
    yy = atom.ry[i] - center_y;
    rad[i-1] = sqrt(xx*xx + yy*yy);
    the[i-1] = theta(xx, yy)/M_PI*180.0;
  }
  // sort r(theta) in order of theta
  sort(num_elem, the, nod);
  for (int i=0; i<num_elem; i++) { dum[i] = rad[i]; }
  for (int i=0; i<num_elem; i++) { rad[i] = dum[nod[i]]; }
  //for (int i=0; i<num_elem; i++) { printf("%d %e %e\n",i,the[i],rad[i]); }
  // sieve out outermost wall
  double threshold = 0.2 * ang;
  for (int i=1; i<num_elem; i++) {
    double drad = rad[i] - rad[i-1];
    if (drad >  threshold) {
      for (int j=i; j<num_elem; j++) {
	rad[j-i] = rad[j]; the[j-i] = the[j];
      }
      num_elem -= i; i=0;
    } else if (drad < -threshold) {
      for (int j=i; j<num_elem-1; j++) {
	rad[j] = rad[j+1]; the[j] = the[j+1];
      }
      num_elem--; i--;
    }
  }
  //for (int i=0; i<num_elem; i++) { printf("%d %e %e\n",i,the[i],rad[i]); }
  // remove crowded points
  threshold = 0.2;
  for (int i=1; i<num_elem; i++) {
    double dthe = the[i] - the[i-1];
    if (fabs(dthe) < threshold) {
      for (int j=i; j<num_elem-1; j++) {
	rad[j] = rad[j+1]; the[j] = the[j+1];
      }
      num_elem--; i--;
    }
  }
  //for (int i=0; i<num_elem; i++) { printf("%d %e %e\n",i,the[i],rad[i]); }
  // calculate path length
  double path = 0.0;
  for (int i=0; i<num_elem; i++) {
    double dthe;
    if (i<num_elem-1) {
      dthe = the[i+1]-the[i];
      if (dthe < 0.0) { dthe += 360.0; }
      dthe *= M_PI/180.0;
      path += (rad[i+1]+rad[i])*dthe/2.0;
    } else {
      dthe = the[0]-the[i];
      if (dthe < 0.0) { dthe += 360.0; }
      dthe *= M_PI/180.0;
      path += (rad[0]+rad[i])*dthe/2.0;
    }
  }
  area = path*cell.hmat[2][2]; //(in m^2)

  // (for debugging) output line shape
  FILE *fp = fopen("cnt_cross_section.d","w");
  for (int i=0; i<num_elem; i++) {
    fprintf(fp,"%f %f\n",the[i],rad[i]*1.0e10);
  }
  fclose(fp);

  // delete array
  delete[] rad; delete[] the; delete[] dum; delete[] nod;
}

// Side area calculation for an arbitrary cross section shape (ver.2)
// Starting from one atom, a line is traced by searching neighbor points
// (ver.1 does not work for the case where the cross section is in a dumbbell shape)
// When CNT contains defect(s), this may not work because a tiny loop around the
// defect can be formed.
void calc_cylinder_side_area_arb2(double &area)
{
  double dsmin = 0.2e-10, dsmax = 2.0e-10;
  double dsminsq = dsmin*dsmin;
  double dsmaxsq = dsmax*dsmax;
  bool debug = true;
  double *crx, *cry; // coordinates of points for cross sectional shape
  crx = new double[atom.natom+1];
  cry = new double[atom.natom+1];
  bool *flg; // flag for search
  flg = new bool[atom.natom+1];
  double rr, dss; // for calculation
  for (int i=0; i<=atom.natom; i++) { flg[i] = true; }

  // calculate center (x,y)
  double center_x, center_y;
  calc_center_xy(center_x, center_y);
  // remove crowded points (flg = false means out of search list)
  for (int i=1; i<=atom.natom-1; i++) {
    if (flg[i]) {
      for (int j=i+1; j<=atom.natom; j++) {
	if (flg[j]) { // if i & j are too close, j is removed from search list
	  dss = dsquare(atom.rx[j]-atom.rx[i])+dsquare(atom.ry[j]-atom.ry[i]);
	  if (dss<dsminsq) { flg[j]=false; }
	}
      }
    }
  }
  // count number of points
  int npts = 0;
  for (int i=1; i<=atom.natom; i++) { if (flg[i]) npts++; }
  // find farthest atom (which should be in outermost wall)
  int ifar = 0; double rrfar = 0;
  for (int i=1; i<=atom.natom; i++) {
    if (flg[i]) {
      rr = dsquare(atom.rx[i]-center_x)+dsquare(atom.ry[i]-center_y);
      if (rr > rrfar) { rrfar = rr; ifar = i; }
    }
  }
  // (for debugging) output line shape
  if (debug) {
    FILE *fp = fopen("cnt_cross_section_shape.d","w");
    for (int i=1; i<atom.natom; i++) {
      if (flg[i]) fprintf(fp,"%f %f\n",atom.rx[i]*1.0e10,atom.ry[i]*1.0e10);
    }
    fclose(fp);
  }
  
  // find points for outermost wall (starting from ifar)
  flg[ifar] = false;
  int itag = ifar; int npts_out = 1; int inn; double xnn, ynn;
  crx[1]=atom.rx[ifar]; cry[1]=atom.ry[ifar];
  for (int k=1; k<=npts; k++) { // k is dummy
    double rrmin=1.0e10;
    for (int i=1; i<=atom.natom; i++) {
      if (flg[i]) {
	dss = dsquare(atom.rx[i]-atom.rx[itag])+dsquare(atom.ry[i]-atom.ry[itag]);
	if (dss<rrmin) { rrmin = dss; inn = i; xnn = atom.rx[i]; ynn = atom.ry[i]; }
      }
    }
    if (rrmin<dsmaxsq) {
      flg[inn] = false; npts_out++; crx[npts_out] = xnn; cry[npts_out] = ynn; itag = inn;
    } else {
      break;
    }
  }
  crx[0] = crx[npts_out];cry[0] = cry[npts_out];
  
  // (for debugging) output line shape
  if (debug) {
    FILE *fpp = fopen("cnt_cross_section_shape_outermost_sort.d","w");
    for (int i=0; i<=npts_out; i++) {
      fprintf(fpp,"%f %f\n",crx[i]*1.0e10,cry[i]*1.0e10);
    }
    fclose(fpp);
  }

  // calculate length of line (and side area)
  double path = 0.0;
  for (int i=1; i<=npts_out; i++) {
    path += sqrt(dsquare(crx[i]-crx[i-1])+dsquare(cry[i]-cry[i-1]));
  }
  area = path*cell.hmat[2][2]; //(in m^2)

  delete[] crx; delete[] cry; delete[] flg;
}

// Side area calculation for an arbitrary cross section shape
// Model is sliced normal to z axis and each slice is evaluated.
// Traced lines for slices are written in 'cnt_cross_section_shape_outermost_sort.d'
// (Use splot in gnuplot to check the traced lines)
void calc_cylinder_side_area_arb_slice(double &area)
{
  double xx, yy;
  double *rad, *the, *dum; // radius, theta, dummy
  int *nod; // order
  // nod[] (1d array) is used to store order of values
  int num_elem = atom.natom; // # of mesh for r(theta)
  rad = new double[atom.natom]; // r
  the = new double[atom.natom]; // theta
  dum = new double[atom.natom]; // dummy
  nod = new int[atom.natom]; // order
  // calculate center (x,y)
  double center_x, center_y;
  calc_center_xy(center_x, center_y);
  // store r(theta)
  for (int i=1; i<=atom.natom; i++) {
    xx = atom.rx[i] - center_x;
    yy = atom.ry[i] - center_y;
    rad[i-1] = sqrt(xx*xx + yy*yy);
    the[i-1] = theta(xx, yy)/M_PI*180.0;
  }
  // sort r(theta) in order of theta
  sort(num_elem, the, nod);
  for (int i=0; i<num_elem; i++) { dum[i] = rad[i]; }
  for (int i=0; i<num_elem; i++) { rad[i] = dum[nod[i]]; }
  //for (int i=0; i<num_elem; i++) { printf("%d %e %e\n",i,the[i],rad[i]); }
  // sieve out outermost wall
  double threshold = 0.2 * ang;
  for (int i=1; i<num_elem; i++) {
    double drad = rad[i] - rad[i-1];
    if (drad >  threshold) {
      for (int j=i; j<num_elem; j++) {
	rad[j-i] = rad[j]; the[j-i] = the[j];
      }
      num_elem -= i; i=0;
    } else if (drad < -threshold) {
      for (int j=i; j<num_elem-1; j++) {
	rad[j] = rad[j+1]; the[j] = the[j+1];
      }
      num_elem--; i--;
    }
  }
  //for (int i=0; i<num_elem; i++) { printf("%d %e %e\n",i,the[i],rad[i]); }
  // remove crowded points
  threshold = 0.2;
  for (int i=1; i<num_elem; i++) {
    double dthe = the[i] - the[i-1];
    if (fabs(dthe) < threshold) {
      for (int j=i; j<num_elem-1; j++) {
	rad[j] = rad[j+1]; the[j] = the[j+1];
      }
      num_elem--; i--;
    }
  }
  //for (int i=0; i<num_elem; i++) { printf("%d %e %e\n",i,the[i],rad[i]); }
  // calculate path length
  double path = 0.0;
  for (int i=0; i<num_elem; i++) {
    double dthe;
    if (i<num_elem-1) {
      dthe = the[i+1]-the[i];
      if (dthe < 0.0) { dthe += 360.0; }
      dthe *= M_PI/180.0;
      path += (rad[i+1]+rad[i])*dthe/2.0;
    } else {
      dthe = the[0]-the[i];
      if (dthe < 0.0) { dthe += 360.0; }
      dthe *= M_PI/180.0;
      path += (rad[0]+rad[i])*dthe/2.0;
    }
  }
  area = path*cell.hmat[2][2]; //(in m^2)

  // (for debugging) output line shape
  FILE *fp = fopen("cnt_cross_section.d","w");
  for (int i=0; i<num_elem; i++) {
    fprintf(fp,"%f %f\n",the[i],rad[i]*1.0e10);
  }
  fclose(fp);

  // delete array
  delete[] rad; delete[] the; delete[] dum; delete[] nod;
}
