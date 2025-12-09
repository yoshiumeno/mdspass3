#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <iomanip>
#include "myheader.h"

#if !defined __linux__ && !defined __APPLE__
#define snprintf sprintf_s
#endif

//Definition for curve fitting
//#define N 3             //N-1 th order (N>1)
#define S 11             //number of data


extern "C"
  long int dgeev_(const char *jobvl, const char *jobvr, long int *n, double *a, 
	     long int *lda, double *wr, double *wi, double *vl, 
	     long int *ldvl, double *vr, long int *ldvr, double *work, 
	     long int *lwork, long int *info);
extern "C"
  long int sgeev_(char *jobvl, char *jobvr, long int *n, double *a, 
	     long int *lda, double *wr, double *wi, double *vl, 
	     long int *ldvl, double *vr, long int *ldvr, double *work, 
	     long int *lwork, long int *info);

extern float cntshell_dmu, cntshell_dnu, cntshell_eps;
extern int cntshell_n;
extern float cnt_pressure;

void calc_center_xy(double &x, double &y);
double theta(double x, double y);
double dsquare(double d);
void cntshell_deform(double dmu, double dnu, int n);
void cntshell_deform_2(double dmu1, double dnu1, int n1, double dmu2, double dnu2, int n2);
void cntshell_calc_de(double dmu, double dnu, int n, double &ext_work, double &de_pot, double &de_tot,
		   bool nomove);
void cntshell_calc_de_2(double dmu1, double dnu1, int n1, double dmu2, double dnu2, int n2, double &ext_work, double &de_pot, double &de_tot,
		   bool nomove);
void bookkeep(); void potential(); void loading();
void gauss_n5(double a[5][6],double xx[5]);
void gauss_n3(double a[3][4],double xx[3]);
void gauss_n2(double a[2][3],double xx[2]);
void curvefit_n5(double x[S],double y[S],double xx[5]);
void curvefit_n3(double x[S],double y[S],double xx[3]);
void curvefit_n2(double x[S],double y[S],double xx[2]);

void cntshell_calc_1_old()
{
  FILE *fp;
  fp = fopen("shell.d","w");
  printf("### CNT 'shell theory' calculation  ###\n");
  double ext_work, de_pot, de_tot;
  double eps = cntshell_eps;
  int rng = 10;
  bookkeep(); potential(); loading();
  printf("# i  j   dW   dE_p    dE_p - dW \n");
  for (int i=-rng; i<=rng; i++) {
    for (int j=-rng; j<=rng; j++) {
      cntshell_calc_de(eps*(double)i, eps*(double)j, cntshell_n, ext_work, de_pot, de_tot, true);
      printf("%d %d   %e  %e    %e\n",i,j,ext_work,de_pot,de_tot);
      fprintf(fp,"%d %d   %e  %e    %e\n",i,j,ext_work,de_pot,de_tot);
    }
    printf("\n");
    fprintf(fp,"\n");
  }
  fclose(fp); 
}
void cntshell_calc_1()
{
  double mat[2][2], evec[2][2];
  printf("### CNT 'shell theory' calculation (single n) ###\n");
  double ext_work, de_pot, de_tot;
  double eps = cntshell_eps; double eps2 = eps*eps;
  double epsd = eps;
  double de, de1, de2;
  bookkeep(); potential(); loading();

  double x[S],y[S],coe_n3[3],coe_n2[2],coe_n5[5];
  int n0=-S/2; int n1=n0+S-1;

  if (true) {
    FILE *fp;
    fp = fopen("esurf_check.d","w");
    for (int i=n0; i<=n1; i++) {
      for (int j=n0; j<=n1; j++) {
	double x0 = double(i)*epsd/double(n1-n0);
	double y0 = double(j)*eps/double(n1-n0);
	cntshell_calc_de(x0,y0,2, ext_work, de_pot, de_tot, true);
	fprintf(fp,"%e %e %e \n",x0,y0,de_tot);
      }
      fprintf(fp,"\n");
    }
    fclose(fp);
  }

  for (int i=n0; i<=n1; i++) {
    double xx = double(i)*eps/double(n1-n0);
    cntshell_calc_de(xx,0,2, ext_work, de_pot, de_tot, true);
    x[i-n0]=xx;y[i-n0]=de_tot;
  }
  curvefit_n3(x,y,coe_n3);// printf("%e\n",coe_n3[2]*2);
  //curvefit_n5(x,y,coe_n5);// printf("%e\n",coe_n3[2]*2);
  //printf("n5 %f   n3 %f\n",coe_n5[2]*2,coe_n3[2]*2);
  mat[0][0] = coe_n3[2]*2;

  for (int i=n0; i<=n1; i++) {
    double xx = double(i)*eps/double(n1-n0);
    cntshell_calc_de(0,xx,2, ext_work, de_pot, de_tot, true);
    x[i-n0]=xx;y[i-n0]=de_tot;
  }
  curvefit_n3(x,y,coe_n3);// printf("%e\n",coe_n3[2]*2);
  mat[1][1] = coe_n3[2]*2;

  double xp[S],yp[S];
  for (int ip=n0; ip<=n1; ip++) {
    double xxp = double(ip)*eps/double(n1-n0);

    for (int i=n0; i<=n1; i++) {
      double xx = double(i)*eps/double(n1-n0);
      cntshell_calc_de(xxp,xx,2, ext_work, de_pot, de_tot, true);
      x[i-n0]=xx;y[i-n0]=de_tot;
    }
    curvefit_n2(x,y,coe_n2);// printf("%e\n",coe_n3[2]*2);
    xp[ip-n0]=xxp;yp[ip-n0]=coe_n2[1];
  /*
  mat[0][1] = coe_n2[1];
  for (int i=n0; i<=n1; i++) {
    double xx = double(i)*eps/double(n1-n0);
    cntshell_calc_de(-eps,xx,2, ext_work, de_pot, de_tot, true);
    x[i-n0]=xx;y[i-n0]=de_tot;
  }
  curvefit_n2(x,y,coe_n2);// printf("%e\n",coe_n3[2]*2);
  */
  }
  curvefit_n2(xp,yp,coe_n2);
  mat[0][1] = coe_n2[1];
  mat[1][0] = mat[0][1];
  printf(" %f %f %f\n",mat[0][0],mat[1][1],mat[0][1]);
  /*
  de = 0.0;
  cntshell_calc_de( eps,0,2, ext_work, de_pot, de_tot, true); de += de_tot;
  cntshell_calc_de(-eps,0,2, ext_work, de_pot, de_tot, true); de += de_tot;
  cntshell_calc_de(   0,0,2, ext_work, de_pot, de_tot, true); de -= de_tot*2; de /= eps2;
  mat[0][0] = de;
  de = 0.0;
  cntshell_calc_de(0, eps,2, ext_work, de_pot, de_tot, true); de += de_tot;
  cntshell_calc_de(0,-eps,2, ext_work, de_pot, de_tot, true); de += de_tot;
  cntshell_calc_de(0,   0,2, ext_work, de_pot, de_tot, true); de -= de_tot*2; de /= eps2;
  mat[1][1] = de;

  de = 0.0; de1 = 0.0; de2 = 0.0;
  cntshell_calc_de( eps, eps,2, ext_work, de_pot, de_tot, true); de1 += de_tot;
  cntshell_calc_de( eps,-eps,2, ext_work, de_pot, de_tot, true); de1 -= de_tot; de1 /= eps*2;
  cntshell_calc_de(-eps, eps,2, ext_work, de_pot, de_tot, true); de2 += de_tot;
  cntshell_calc_de(-eps,-eps,2, ext_work, de_pot, de_tot, true); de2 -= de_tot; de2 /= eps*2;
  mat[0][1] = (de1 - de2) / (eps*2); mat[1][0] = mat[0][1];
  printf(" %f %f %f\n",mat[0][0],mat[1][1],mat[0][1]);
  */

  long int nn, nn4, info;
  nn=2; nn4=nn*4;
  double *ar = (double *)calloc(nn*nn, sizeof(double)); 
  double *e  = (double *)calloc(nn, sizeof(double)); 
  double *ei = (double *)calloc(nn, sizeof(double)); 
  double *evi = (double *)calloc(nn*nn, sizeof(double)); 
  double *evr = (double *)calloc(nn*nn, sizeof(double)); 
  double *work = (double *)calloc(nn4, sizeof(double)); 

  printf("Elements:\n");
  for (int i=0; i<nn; i++) {
    printf("%d %f,  ",i,mat[i][i]);
  }
  printf("%d %f,  ",2,mat[0][1]);
  printf("\n");

  int ii=0;
  for (int i=0; i<nn; i++) {
    for (int j=0; j<nn; j++) {
      ar[ii]=mat[i][j]; ii++;
      //printf("%d %d %e\n",i,j,mat[i][j]);
    }
    //printf("\n");
  }
  dgeev_("N","V",&nn,ar,&nn,e,ei,evi,&nn,evr,&nn,work,&nn4,&info);
  ii=0;
  printf("CNT pressure value: %f\n",cnt_pressure);
  for (int i=0; i<nn; i++) {
    printf("eigenvalue: %d %f",i,e[i]);
    printf(" (");
    for (int j=0; j<nn; j++) {
      evec[i][j] = evr[ii];
      printf("  %f ", evr[ii]);
      ii++;
    }
    printf(")\n");
  }

  // Verification of eigenvalue solution
  if (false) {
    double xx;
    for (int n=0; n<nn; n++) { printf("%d -th mode\n",n);
      for (int i=0; i<nn; i++) {
	xx = 0;
	for (int k=0; k<nn; k++) { xx += mat[i][k]*evec[n][k]; }
	printf("%d %f --  %f\n",i,xx,e[n]*evec[n][i]);
      }
      printf("\n",n);
    }
  }

  // Calculate curve along eigenmodes
  if (true) {
    FILE *fp0;
    fp0 = fopen("ecurv_check.d","w");
    double xx, yy;
    for (int i=-50; i<=50; i++) {
      double fac = (double)i*eps/50 * 5;
      double x0 = evec[0][0]*fac; double y0 = evec[0][1]*fac;
      cntshell_calc_de(x0,y0,2, ext_work, de_pot, de_tot, true);
      fprintf(fp0,"%e  %e ",fac,de_tot);
      x0 = evec[1][0]*fac; y0 = evec[1][1]*fac;
      cntshell_calc_de(x0,y0,2, ext_work, de_pot, de_tot, true);
      fprintf(fp0," %e \n",de_tot);
    }
    fclose(fp0);
  }


}
void cntshell_calc_2()
{
  double mat[4][4], evec[4][4];
  printf("### CNT 'shell theory' calculation (double n) ###\n");
  double ext_work, de_pot, de_tot;
  double eps = cntshell_eps; double eps2 = eps*eps;
  double de, de1, de2;
  bookkeep(); potential(); loading();

  de = 0.0;
  cntshell_calc_de_2( eps,0,2, 0,0,4, ext_work, de_pot, de_tot, true); de += de_tot;
  cntshell_calc_de_2(-eps,0,2, 0,0,4, ext_work, de_pot, de_tot, true); de += de_tot;
  cntshell_calc_de_2(   0,0,2, 0,0,4, ext_work, de_pot, de_tot, true); de -= de_tot*2; de /= eps2;
  mat[0][0] = de;
  de = 0.0;
  cntshell_calc_de_2(0, eps,2, 0,0,4, ext_work, de_pot, de_tot, true); de += de_tot;
  cntshell_calc_de_2(0,-eps,2, 0,0,4, ext_work, de_pot, de_tot, true); de += de_tot;
  cntshell_calc_de_2(0,   0,2, 0,0,4, ext_work, de_pot, de_tot, true); de -= de_tot*2; de /= eps2;
  mat[1][1] = de;
  de = 0.0;
  cntshell_calc_de_2(0,0,2,  eps,0,4, ext_work, de_pot, de_tot, true); de += de_tot;
  cntshell_calc_de_2(0,0,2, -eps,0,4, ext_work, de_pot, de_tot, true); de += de_tot;
  cntshell_calc_de_2(0,0,2,    0,0,4, ext_work, de_pot, de_tot, true); de -= de_tot*2; de /= eps2;
  mat[2][2] = de;
  de = 0.0;
  cntshell_calc_de_2(0,0,2, 0, eps,4, ext_work, de_pot, de_tot, true); de += de_tot;
  cntshell_calc_de_2(0,0,2, 0,-eps,4, ext_work, de_pot, de_tot, true); de += de_tot;
  cntshell_calc_de_2(0,0,2, 0,   0,4, ext_work, de_pot, de_tot, true); de -= de_tot*2; de /= eps2;
  mat[3][3] = de;

  de = 0.0; de1 = 0.0; de2 = 0.0;
  cntshell_calc_de_2( eps, eps,2, 0,0,4, ext_work, de_pot, de_tot, true); de1 += de_tot;
  cntshell_calc_de_2( eps,-eps,2, 0,0,4, ext_work, de_pot, de_tot, true); de1 -= de_tot; de1 /= eps*2;
  cntshell_calc_de_2(-eps, eps,2, 0,0,4, ext_work, de_pot, de_tot, true); de2 += de_tot;
  cntshell_calc_de_2(-eps,-eps,2, 0,0,4, ext_work, de_pot, de_tot, true); de2 -= de_tot; de2 /= eps*2;
  mat[0][1] = (de1 - de2) / (eps*2); mat[1][0] = mat[0][1];
  de = 0.0; de1 = 0.0; de2 = 0.0;
  cntshell_calc_de_2( eps,0,2, eps,0,4, ext_work, de_pot, de_tot, true); de1 += de_tot;
  cntshell_calc_de_2( eps,0,2,-eps,0,4, ext_work, de_pot, de_tot, true); de1 -= de_tot; de1 /= eps*2;
  cntshell_calc_de_2(-eps,0,2, eps,0,4, ext_work, de_pot, de_tot, true); de2 += de_tot;
  cntshell_calc_de_2(-eps,0,2,-eps,0,4, ext_work, de_pot, de_tot, true); de2 -= de_tot; de2 /= eps*2;
  mat[0][2] = (de1 - de2) / (eps*2); mat[2][0] = mat[0][2];
  de = 0.0; de1 = 0.0; de2 = 0.0;
  cntshell_calc_de_2( eps,0,2, 0, eps,4, ext_work, de_pot, de_tot, true); de1 += de_tot;
  cntshell_calc_de_2( eps,0,2, 0,-eps,4, ext_work, de_pot, de_tot, true); de1 -= de_tot; de1 /= eps*2;
  cntshell_calc_de_2(-eps,0,2, 0, eps,4, ext_work, de_pot, de_tot, true); de2 += de_tot;
  cntshell_calc_de_2(-eps,0,2, 0,-eps,4, ext_work, de_pot, de_tot, true); de2 -= de_tot; de2 /= eps*2;
  mat[0][3] = (de1 - de2) / (eps*2); mat[3][0] = mat[0][3];

  de = 0.0; de1 = 0.0; de2 = 0.0;
  cntshell_calc_de_2(0, eps,2, eps,0,4, ext_work, de_pot, de_tot, true); de1 += de_tot;
  cntshell_calc_de_2(0, eps,2,-eps,0,4, ext_work, de_pot, de_tot, true); de1 -= de_tot; de1 /= eps*2;
  cntshell_calc_de_2(0,-eps,2, eps,0,4, ext_work, de_pot, de_tot, true); de2 += de_tot;
  cntshell_calc_de_2(0,-eps,2,-eps,0,4, ext_work, de_pot, de_tot, true); de2 -= de_tot; de2 /= eps*2;
  mat[1][2] = (de1 - de2) / (eps*2); mat[2][1] = mat[1][2];
  de = 0.0; de1 = 0.0; de2 = 0.0;
  cntshell_calc_de_2(0, eps,2, 0, eps,4, ext_work, de_pot, de_tot, true); de1 += de_tot;
  cntshell_calc_de_2(0, eps,2, 0,-eps,4, ext_work, de_pot, de_tot, true); de1 -= de_tot; de1 /= eps*2;
  cntshell_calc_de_2(0,-eps,2, 0, eps,4, ext_work, de_pot, de_tot, true); de2 += de_tot;
  cntshell_calc_de_2(0,-eps,2, 0,-eps,4, ext_work, de_pot, de_tot, true); de2 -= de_tot; de2 /= eps*2;
  mat[1][3] = (de1 - de2) / (eps*2); mat[3][1] = mat[1][3];

  de = 0.0; de1 = 0.0; de2 = 0.0;
  cntshell_calc_de_2(0,0,2, eps, eps,4, ext_work, de_pot, de_tot, true); de1 += de_tot;
  cntshell_calc_de_2(0,0,2, eps,-eps,4, ext_work, de_pot, de_tot, true); de1 -= de_tot; de1 /= eps*2;
  cntshell_calc_de_2(0,0,2,-eps, eps,4, ext_work, de_pot, de_tot, true); de2 += de_tot;
  cntshell_calc_de_2(0,0,2,-eps,-eps,4, ext_work, de_pot, de_tot, true); de2 -= de_tot; de2 /= eps*2;
  mat[2][3] = (de1 - de2) / (eps*2); mat[3][2] = mat[2][3];

  long int nn, nn4, info;
  nn=4; nn4=nn*4;
  double *ar = (double *)calloc(nn*nn, sizeof(double)); 
  double *e  = (double *)calloc(nn, sizeof(double)); 
  double *ei = (double *)calloc(nn, sizeof(double)); 
  double *evi = (double *)calloc(nn*nn, sizeof(double)); 
  double *evr = (double *)calloc(nn*nn, sizeof(double)); 
  double *work = (double *)calloc(nn4, sizeof(double)); 

  printf("Diagonal elements:\n");
  for (int i=0; i<nn; i++) {
    printf("%d %f,  ",i,mat[i][i]);
  }
  printf("\n");

  int ii=0;
  for (int i=0; i<nn; i++) {
    for (int j=0; j<nn; j++) {
      ar[ii]=mat[i][j]; ii++;
      //printf("%d %d %e\n",i,j,mat[i][j]);
    }
    //printf("\n");
  }
  dgeev_("N","V",&nn,ar,&nn,e,ei,evi,&nn,evr,&nn,work,&nn4,&info);
  ii=0;
  for (int i=0; i<nn; i++) {
    printf("eigenvalue: %d %f\n",i,e[i]);
    printf(" (");
    for (int j=0; j<nn; j++) {
      evec[i][j] = evr[ii];
      printf("  %f ", evr[ii]);
      ii++;
    }
    printf(")\n");
  }

  // Verification of eigenvalue solution
  if (false) {
    double xx;
    for (int n=0; n<nn; n++) { printf("%d -th mode\n",n);
      for (int i=0; i<nn; i++) {
	xx = 0;
	for (int k=0; k<nn; k++) { xx += mat[i][k]*evec[n][k]; }
	printf("%d %f --  %f\n",i,xx,e[n]*evec[n][i]);
      }
      printf("\n",n);
    }
  }

}

void cntshell_calc_test(double dmu, double dnu, int n)
{
  double ext_work, de_pot, de_tot;
  cntshell_calc_de(dmu, dnu, n, ext_work, de_pot, de_tot, false);
  printf("dW (work by load) = %e   dE_pot = %f  dE_tot (dE_pot-dW) %f\n",
  	 ext_work,de_pot,de_tot);
}

void cntshell_calc_de(double dmu, double dnu, int n, double &ext_work, double &de_pot, double &de_tot,
		   bool nomove)
{
  double dx, dy, dz;
  for (int i=1; i<=atom.natom; i++) {
    atom.rx_p[i] = atom.rx[i]; atom.ry_p[i] = atom.ry[i]; atom.rz_p[i] = atom.rz[i];
  }
  bookkeep(); potential(); loading();epotatom=atom.epotsum/eV/atom.natom;
  double epot_bfr = atom.epotsum/eV/atom.natom;
  ext_work = 0;
  cntshell_deform(dmu, dnu, n);
  for (int i=1; i<=atom.natom; i++) {
    dx = atom.rx[i] - atom.rx_p[i];
    dy = atom.ry[i] - atom.ry_p[i];
    dz = atom.rz[i] - atom.rz_p[i];
    //ext_work += (atom.fx_l[i] * dx + atom.fy_l[i] * dy + atom.fz_l[i] * dz)/eV/ang;
    ext_work += (atom.fx_l[i] * dx + atom.fy_l[i] * dy + atom.fz_l[i] * dz)/eV;
  }
  ext_work /= (double)atom.natom;
  bookkeep(); potential(); loading();epotatom=atom.epotsum/eV/atom.natom;
  double epot_aft = atom.epotsum/eV/atom.natom;
  de_pot = epot_aft-epot_bfr;
  de_tot = de_pot - ext_work;
  //printf("E_bfr = %f  E_aft = %f  E_ext = %e   dE_pot = %e  dE_tot = %e\n",
  //	 epot_bfr,epot_aft,ext_work,de_pot,de_tot);
  if (nomove) {
    for (int i=1; i<=atom.natom; i++) {
      atom.rx[i] = atom.rx_p[i]; atom.ry[i] = atom.ry_p[i]; atom.rz[i] = atom.rz_p[i];
    }
    bookkeep(); potential(); loading();epotatom=atom.epotsum/eV/atom.natom;
  }
}

void cntshell_calc_de_2(double dmu1, double dnu1, int n1, double dmu2, double dnu2, int n2, double &ext_work, double &de_pot, double &de_tot,
		   bool nomove)
{
  double dx, dy, dz;
  for (int i=1; i<=atom.natom; i++) {
    atom.rx_p[i] = atom.rx[i]; atom.ry_p[i] = atom.ry[i]; atom.rz_p[i] = atom.rz[i];
  }
  bookkeep(); potential(); loading();epotatom=atom.epotsum/eV/atom.natom;
  double epot_bfr = atom.epotsum/eV/atom.natom;
  ext_work = 0;
  cntshell_deform_2(dmu1, dnu1, n1, dmu2, dnu2, n2);
  for (int i=1; i<=atom.natom; i++) {
    dx = atom.rx[i] - atom.rx_p[i];
    dy = atom.ry[i] - atom.ry_p[i];
    dz = atom.rz[i] - atom.rz_p[i];
    //ext_work += (atom.fx_l[i] * dx + atom.fy_l[i] * dy + atom.fz_l[i] * dz)/eV/ang;
    ext_work += (atom.fx_l[i] * dx + atom.fy_l[i] * dy + atom.fz_l[i] * dz)/eV;
  }
  ext_work /= (double)atom.natom;
  bookkeep(); potential(); loading();epotatom=atom.epotsum/eV/atom.natom;
  double epot_aft = atom.epotsum/eV/atom.natom;
  de_pot = epot_aft-epot_bfr;
  de_tot = de_pot - ext_work;
  //printf("E_bfr = %f  E_aft = %f  E_ext = %e   dE_pot = %f  dE_tot = %f\n",
  //	 epot_bfr,epot_aft,ext_work,de_pot,de_tot);
  if (nomove) {
    for (int i=1; i<=atom.natom; i++) {
      atom.rx[i] = atom.rx_p[i]; atom.ry[i] = atom.ry_p[i]; atom.rz[i] = atom.rz_p[i];
    }
    bookkeep(); potential(); loading();epotatom=atom.epotsum/eV/atom.natom;
  }
}

void cntshell_deform(double dmu, double dnu, int n)
{
  //dmu *= 1.13*1.13;
  double dx[2]; dx[0] = (double)dmu*ang; dx[1] = (double)dnu*ang;
  double center_x, center_y;
  calc_center_xy(center_x, center_y);
  for (int i=1; i<=atom.natom; i++) {
    double x = atom.rx[i]; double y = atom.ry[i];
    double r = sqrt(dsquare(x-center_x)+dsquare(y-center_y));
    double t = theta(x-center_x,y-center_y);
    double du = dx[0] * cos(double(n)*t);
    double dv = dx[1] * sin(double(n)*t);
    atom.rx[i] = center_x + (r + du) * cos(t) - dv * sin(t);
    atom.ry[i] = center_y + (r + du) * sin(t) + dv * cos(t);
  }

}

void cntshell_deform_2(double dmu1, double dnu1, int n1, double dmu2, double dnu2, int n2)
{
  double dx1[2]; dx1[0] = (double)dmu1*ang; dx1[1] = (double)dnu1*ang;
  double dx2[2]; dx2[0] = (double)dmu2*ang; dx2[1] = (double)dnu2*ang;
  double center_x, center_y;
  calc_center_xy(center_x, center_y);
  for (int i=1; i<=atom.natom; i++) {
    double x = atom.rx[i]; double y = atom.ry[i];
    double r = sqrt(dsquare(x-center_x)+dsquare(y-center_y));
    double t = theta(x-center_x,y-center_y);
    double du = dx1[0] * cos(double(n1)*t);
    double dv = dx1[1] * sin(double(n1)*t);
    du += dx2[0] * cos(double(n2)*t);
    dv += dx2[1] * sin(double(n2)*t);
    atom.rx[i] = center_x + (r + du) * cos(t) - dv * sin(t);
    atom.ry[i] = center_y + (r + du) * sin(t) + dv * cos(t);
  }

}





void curvefit_n2(double x[S],double y[S],double xx[3])
{
  //  #define n 3
        int i,j,k;
        double X,Y;
	//        double A[n][n+1],xx[N];
        double A[2][2+1];//,xx[N];
	int n = 2;

        //FILE *output1;
        //FILE *output2;
        //output1=fopen("output1.data","w");
        //output2=fopen("output2.data","w");

        for(i=0;i<n;i++) {
                for(j=0;j<n+1;j++) {
                        A[i][j]=0.0;
                }
        }

        for(i=0;i<n;i++) {
                for(j=0;j<n;j++) {
                        for(k=0;k<S;k++) {
                                A[i][j]+=pow(x[k],i+j);
                        }
                }
        }
        for(i=0;i<n;i++) {
                for(k=0;k<S;k++) {
                        A[i][n]+=pow(x[k],i)*y[k];
                }
        }

        gauss_n2(A,xx);
	//printf("%e %e %e\n",xx[0],xx[1],xx[2]);

}

void curvefit_n3(double x[S],double y[S],double xx[3])
{
  //  #define n 3
        int i,j,k;
        double X,Y;
	//        double A[n][n+1],xx[N];
        double A[3][3+1];//,xx[N];
	int n = 3;

        //FILE *output1;
        //FILE *output2;
        //output1=fopen("output1.data","w");
        //output2=fopen("output2.data","w");

        for(i=0;i<n;i++) {
                for(j=0;j<n+1;j++) {
                        A[i][j]=0.0;
                }
        }

        for(i=0;i<n;i++) {
                for(j=0;j<n;j++) {
                        for(k=0;k<S;k++) {
                                A[i][j]+=pow(x[k],i+j);
                        }
                }
        }
        for(i=0;i<n;i++) {
                for(k=0;k<S;k++) {
                        A[i][n]+=pow(x[k],i)*y[k];
                }
        }

        gauss_n3(A,xx);
	//printf("%e %e %e\n",xx[0],xx[1],xx[2]);

}

void curvefit_n5(double x[S],double y[S],double xx[5])
{
  //  #define n 3
        int i,j,k;
        double X,Y;
	//        double A[n][n+1],xx[N];
        double A[5][5+1];//,xx[N];
	int n = 5;

        //FILE *output1;
        //FILE *output2;
        //output1=fopen("output1.data","w");
        //output2=fopen("output2.data","w");

        for(i=0;i<n;i++) {
                for(j=0;j<n+1;j++) {
                        A[i][j]=0.0;
                }
        }

        for(i=0;i<n;i++) {
                for(j=0;j<n;j++) {
                        for(k=0;k<S;k++) {
                                A[i][j]+=pow(x[k],i+j);
                        }
                }
        }
        for(i=0;i<n;i++) {
                for(k=0;k<S;k++) {
                        A[i][n]+=pow(x[k],i)*y[k];
                }
        }

        gauss_n5(A,xx);
	//printf("%e %e %e\n",xx[0],xx[1],xx[2]);

}

void gauss_n2(double a[2][3],double xx[3])
{
        int i,j,k,l,pivot;
        double x[2];
        double p,q,m,b[1][3];

        for(i=0;i<2;i++) {
                m=0;
                pivot=i;

                for(l=i;l<2;l++) {
                        if(fabs(a[l][i])>m) {   //i列の中で一番値が大きい行を選ぶ
                                m=fabs(a[l][i]);
                                pivot=l;
                        }
                }

                if(pivot!=i) {                          //pivotがiと違えば、行の入れ替え
                        for(j=0;j<2+1;j++) {
                                b[0][j]=a[i][j];        
                                a[i][j]=a[pivot][j];
                                a[pivot][j]=b[0][j];
                        }
                }
        }

        for(k=0;k<2;k++) {
                p=a[k][k];              //対角要素を保存
                a[k][k]=1;              //対角要素は１になることがわかっているから

                for(j=k+1;j<2+1;j++) {
                        a[k][j]/=p;
                }

                for(i=k+1;i<2;i++) {
                        q=a[i][k];

                        for(j=k+1;j<2+1;j++) {
                                a[i][j]-=q*a[k][j];
                        }
                a[i][k]=0;              //０となることがわかっているところ
                }
        }

//解の計算
        for(i=2-1;i>=0;i--) {
                x[i]=a[i][2];
                for(j=2-1;j>i;j--) {
                        x[i]-=a[i][j]*x[j];
                }
        }

        for(i=0;i<2;i++) {
                xx[i]=x[i];
        }

}
void gauss_n3(double a[3][4],double xx[3])
{
        int i,j,k,l,pivot;
        double x[3];
        double p,q,m,b[1][3+1];

        for(i=0;i<3;i++) {
                m=0;
                pivot=i;

                for(l=i;l<3;l++) {
                        if(fabs(a[l][i])>m) {   //i列の中で一番値が大きい行を選ぶ
                                m=fabs(a[l][i]);
                                pivot=l;
                        }
                }

                if(pivot!=i) {                          //pivotがiと違えば、行の入れ替え
                        for(j=0;j<3+1;j++) {
                                b[0][j]=a[i][j];        
                                a[i][j]=a[pivot][j];
                                a[pivot][j]=b[0][j];
                        }
                }
        }

        for(k=0;k<3;k++) {
                p=a[k][k];              //対角要素を保存
                a[k][k]=1;              //対角要素は１になることがわかっているから

                for(j=k+1;j<3+1;j++) {
                        a[k][j]/=p;
                }

                for(i=k+1;i<3;i++) {
                        q=a[i][k];

                        for(j=k+1;j<3+1;j++) {
                                a[i][j]-=q*a[k][j];
                        }
                a[i][k]=0;              //０となることがわかっているところ
                }
        }

//解の計算
        for(i=3-1;i>=0;i--) {
                x[i]=a[i][3];
                for(j=3-1;j>i;j--) {
                        x[i]-=a[i][j]*x[j];
                }
        }

        for(i=0;i<3;i++) {
                xx[i]=x[i];
        }

}
void gauss_n5(double a[5][6],double xx[5])
{
        int i,j,k,l,pivot;
        double x[5];
        double p,q,m,b[1][5+1];

        for(i=0;i<5;i++) {
                m=0;
                pivot=i;

                for(l=i;l<5;l++) {
                        if(fabs(a[l][i])>m) {   //i列の中で一番値が大きい行を選ぶ
                                m=fabs(a[l][i]);
                                pivot=l;
                        }
                }

                if(pivot!=i) {                          //pivotがiと違えば、行の入れ替え
                        for(j=0;j<5+1;j++) {
                                b[0][j]=a[i][j];        
                                a[i][j]=a[pivot][j];
                                a[pivot][j]=b[0][j];
                        }
                }
        }

        for(k=0;k<5;k++) {
                p=a[k][k];              //対角要素を保存
                a[k][k]=1;              //対角要素は１になることがわかっているから

                for(j=k+1;j<5+1;j++) {
                        a[k][j]/=p;
                }

                for(i=k+1;i<5;i++) {
                        q=a[i][k];

                        for(j=k+1;j<5+1;j++) {
                                a[i][j]-=q*a[k][j];
                        }
                a[i][k]=0;              //０となることがわかっているところ
                }
        }

//解の計算
        for(i=5-1;i>=0;i--) {
                x[i]=a[i][5];
                for(j=5-1;j>i;j--) {
                        x[i]-=a[i][j]*x[j];
                }
        }

        for(i=0;i<5;i++) {
                xx[i]=x[i];
        }

}
