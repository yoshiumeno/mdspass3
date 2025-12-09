#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include "myheader.h"

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

void potential();
void bookkeep();
void sort(int n, double *val, int *nod);
void deform_element(int iel);
double q_x(int i, int iel, double rx[], double ry[], double rz[], double hinelmat[][3]);
double q_y(int i, int iel, double rx[], double ry[], double rz[], double hinelmat[][3]);
double q_z(int i, int iel, double rx[], double ry[], double rz[], double hinelmat[][3]);
void mk_helmat(int iel, double rx[], double ry[], double rz[], double mat[][3]);
void inverse(double mat[][3], double imat[][3]);
double r_x(int i, int iel, double rx[], double ry[], double rz[], double helmat[][3], 
	   double qx, double qy, double qz);
double r_y(int i, int iel, double rx[], double ry[], double rz[], double helmat[][3], 
	   double qx, double qy, double qz);
double r_z(int i, int iel, double rx[], double ry[], double rz[], double helmat[][3], 
	   double qx, double qy, double qz);
void write_evec_xsf(const char* fname, int mod);


void instability_QC()
{
  std::cout<<"Instability analysis.."<<std::endl;
  long int nn, nn4, info;
  double dis = 0.05e-10; // Infinitesimal displacement
  //  double dis = 0.0005e-10; // Infinitesimal displacement

  nn=atom.nrepatom*3-3; nn4=nn*4;
  //  nn=903; nn4=nn*4;
  double *ar = (double *)calloc(nn*nn, sizeof(double)); 
  double *e  = (double *)calloc(nn, sizeof(double)); 
  double *ei = (double *)calloc(nn, sizeof(double)); 
  double *evi = (double *)calloc(nn*nn, sizeof(double)); 
  double *evr = (double *)calloc(nn*nn, sizeof(double)); 
  double *work = (double *)calloc(nn4, sizeof(double)); 
  int *nod = (int *)calloc(nn, sizeof(int));
  int *ia = (int *)calloc(atom.nrepatom+1, sizeof(int));

  // Allocate hessian matrix
  double **hessian;
  hessian = (double**)malloc(sizeof(double*)*(atom.nrepatom*3+1));
  for (int i=1; i<=atom.nrepatom*3; i++) {
    hessian[i] = (double *)malloc(sizeof(double)*(atom.nrepatom*3+1));
  }

  // Store atom position
  for (int i=1; i<=atom.natom; i++) {
    atom.rx_p[i] = atom.rx[i]; atom.ry_p[i] = atom.ry[i]; atom.rz_p[i] = atom.rz[i];
  }

  // Numbering repatoms
  int inum=0;
  for (int i=1; i<=atom.natom; i++) {
    if (atom.repatom[i]) { inum++; ia[inum]=i; }
  }

  bookkeep();
  std::cout<<"Hessian calc (1/3).."<<std::endl;
  for (int i=1; i<=atom.nrepatom; i++) {
    int ii = ia[i];
    atom.rx[ii]=atom.rx[ii]+dis;
    if (atom.elem_id[ii]) { for (int iel=1; iel<=atom.nelem; iel++) { deform_element(iel); } }
    potential();
    for (int j=1; j<=atom.nrepatom; j++) {
      int jj = ia[j];
      hessian[(i-1)*3+1][(j-1)*3+1] = atom.fx[jj];
      hessian[(i-1)*3+1][(j-1)*3+2] = atom.fy[jj];
      hessian[(i-1)*3+1][(j-1)*3+3] = atom.fz[jj];
    }
    atom.rx[ii]=atom.rx[ii]-dis*2;
    if (atom.elem_id[ii]) { for (int iel=1; iel<=atom.nelem; iel++) { deform_element(iel); } }
    potential();
    for (int j=1; j<=atom.nrepatom; j++) {
      int jj = ia[j];
      hessian[(i-1)*3+1][(j-1)*3+1] = (atom.fx[jj]-hessian[(i-1)*3+1][(j-1)*3+1])/dis/2;
      hessian[(i-1)*3+1][(j-1)*3+2] = (atom.fy[jj]-hessian[(i-1)*3+1][(j-1)*3+2])/dis/2;
      hessian[(i-1)*3+1][(j-1)*3+3] = (atom.fz[jj]-hessian[(i-1)*3+1][(j-1)*3+3])/dis/2;
    }
    atom.rx[ii]=atom.rx[ii]+dis;
    if (atom.elem_id[ii]) { for (int iel=1; iel<=atom.nelem; iel++) { deform_element(iel); } }
  }
  std::cout<<"Hessian calc (2/3).."<<std::endl;
  for (int i=1; i<=atom.nrepatom; i++) {
    int ii = ia[i];
    atom.ry[ii]=atom.ry[ii]+dis;
    if (atom.elem_id[ii]) { for (int iel=1; iel<=atom.nelem; iel++) { deform_element(iel); } }
    potential();
    for (int j=1; j<=atom.nrepatom; j++) {
      int jj = ia[j];
      hessian[(i-1)*3+2][(j-1)*3+1] = atom.fx[jj];
      hessian[(i-1)*3+2][(j-1)*3+2] = atom.fy[jj];
      hessian[(i-1)*3+2][(j-1)*3+3] = atom.fz[jj];
    }
    atom.ry[ii]=atom.ry[ii]-dis*2;
    if (atom.elem_id[ii]) { for (int iel=1; iel<=atom.nelem; iel++) { deform_element(iel); } }
    potential();
    for (int j=1; j<=atom.nrepatom; j++) {
      int jj = ia[j];
      hessian[(i-1)*3+2][(j-1)*3+1] = (atom.fx[jj]-hessian[(i-1)*3+2][(j-1)*3+1])/dis/2;
      hessian[(i-1)*3+2][(j-1)*3+2] = (atom.fy[jj]-hessian[(i-1)*3+2][(j-1)*3+2])/dis/2;
      hessian[(i-1)*3+2][(j-1)*3+3] = (atom.fz[jj]-hessian[(i-1)*3+2][(j-1)*3+3])/dis/2;
    }
    atom.ry[ii]=atom.ry[ii]+dis;
    if (atom.elem_id[ii]) { for (int iel=1; iel<=atom.nelem; iel++) { deform_element(iel); } }
  }
  std::cout<<"Hessian calc (3/3).."<<std::endl;
  for (int i=1; i<=atom.nrepatom; i++) {
    int ii = ia[i];
    atom.rz[ii]=atom.rz[ii]+dis;
    if (atom.elem_id[ii]) { for (int iel=1; iel<=atom.nelem; iel++) { deform_element(iel); } }
    potential();
    for (int j=1; j<=atom.nrepatom; j++) {
      int jj = ia[j];
      hessian[(i-1)*3+3][(j-1)*3+1] = atom.fx[jj];
      hessian[(i-1)*3+3][(j-1)*3+2] = atom.fy[jj];
      hessian[(i-1)*3+3][(j-1)*3+3] = atom.fz[jj];
    }
    atom.rz[ii]=atom.rz[ii]-dis*2;
    if (atom.elem_id[ii]) { for (int iel=1; iel<=atom.nelem; iel++) { deform_element(iel); } }
    potential();
    for (int j=1; j<=atom.nrepatom; j++) {
      int jj = ia[j];
      hessian[(i-1)*3+3][(j-1)*3+1] = (atom.fx[jj]-hessian[(i-1)*3+3][(j-1)*3+1])/dis/2;
      hessian[(i-1)*3+3][(j-1)*3+2] = (atom.fy[jj]-hessian[(i-1)*3+3][(j-1)*3+2])/dis/2;
      hessian[(i-1)*3+3][(j-1)*3+3] = (atom.fz[jj]-hessian[(i-1)*3+3][(j-1)*3+3])/dis/2;
    }
    atom.rz[ii]=atom.rz[ii]+dis;
    if (atom.elem_id[ii]) { for (int iel=1; iel<=atom.nelem; iel++) { deform_element(iel); } }
  }

  // Copy hessian to ar.
  std::cout<<atom.instcenter<<" -th atom is fixed."<<std::endl;
  int iii=0;
  for (int i=1; i<=atom.nrepatom; i++) {
    int ii=ia[i];
    if (ii!=atom.instcenter) {
      for (int j=1; j<=atom.nrepatom; j++) {
	int jj=ia[j];
	if (jj!=atom.instcenter) {
	  iii=0;
	  if (jj<atom.instcenter) { iii=iii+(j-1)*3; } else { iii=iii+(j-2)*3; }
	  if (ii<atom.instcenter) { iii=iii+(i-1)*(atom.nrepatom-1)*9; } 
	  else { iii=iii+(i-2)*(atom.nrepatom-1)*9; }
	  ar[iii]=hessian[(i-1)*3+1][(j-1)*3+1]; iii++;
	  ar[iii]=hessian[(i-1)*3+1][(j-1)*3+2]; iii++;
	  ar[iii]=hessian[(i-1)*3+1][(j-1)*3+3]; iii++; iii=iii+nn-3;
	  ar[iii]=hessian[(i-1)*3+2][(j-1)*3+1]; iii++;
	  ar[iii]=hessian[(i-1)*3+2][(j-1)*3+2]; iii++;
	  ar[iii]=hessian[(i-1)*3+2][(j-1)*3+3]; iii++; iii=iii+nn-3;
	  ar[iii]=hessian[(i-1)*3+3][(j-1)*3+1]; iii++;
	  ar[iii]=hessian[(i-1)*3+3][(j-1)*3+2]; iii++;
	  ar[iii]=hessian[(i-1)*3+3][(j-1)*3+3]; iii++; iii=iii+nn-3;
	}
      }
    }
  }

  std::cout<<"Solving eigenvalue problem."<<std::endl;
#undef FLOAT
#undef CUDA
#ifdef CUDA  
  sgeev_("N","V",&nn,ar,&nn,e,ei,evi,&nn,evr,&nn,work,&nn4,&info);
#else
#ifdef FLOAT
  sgeev_("N","V",&nn,ar,&nn,e,ei,evi,&nn,evr,&nn,work,&nn4,&info);
#else
  dgeev_("N","V",&nn,ar,&nn,e,ei,evi,&nn,evr,&nn,work,&nn4,&info);
#endif
#endif
 
  sort(nn, e, nod);

  FILE *outputfile = fopen("eminall.d","w");
  for (int i=0; i<nn; i++) { fprintf (outputfile, "%d %20.15e\n",i,e[i]); }
  fclose(outputfile);
  std::cout<<"Eigenvalues written in eminall.d"<<std::endl;


  for (int mod=0; mod<20; mod++) {
    int imod = nod[mod];
    for (int i=1; i<=atom.natom; i++) {
      atom.evecx[i][mod]=0; atom.evecy[i][mod]=0; atom.evecz[i][mod]=0;
    }
    for (int i=1; i<=atom.nrepatom; i++) {
      int ii=ia[i];
      if (ii!=atom.instcenter) {
	int iii=0; if (ii<atom.instcenter) { iii=(i-1)*3; } else { iii=(i-2)*3; }
	iii=iii+nn*imod;
	atom.evecx[ii][mod]=evr[iii]; atom.evecy[ii][mod]=evr[iii+1]; atom.evecz[ii][mod]=evr[iii+2];
      }
    }
    if (mod== 0) { write_evec_xsf("atominst01.xsf", 0); }
    if (mod== 1) { write_evec_xsf("atominst02.xsf", 1); }
    if (mod== 2) { write_evec_xsf("atominst03.xsf", 2); }
    if (mod== 3) { write_evec_xsf("atominst04.xsf", 3); }
    if (mod== 4) { write_evec_xsf("atominst05.xsf", 4); }
    if (mod== 5) { write_evec_xsf("atominst06.xsf", 5); }
    if (mod== 6) { write_evec_xsf("atominst07.xsf", 6); }
    if (mod== 7) { write_evec_xsf("atominst08.xsf", 7); }
    if (mod== 8) { write_evec_xsf("atominst09.xsf", 8); }
    if (mod== 9) { write_evec_xsf("atominst10.xsf", 9); }
    if (mod==10) { write_evec_xsf("atominst11.xsf",10); }
    if (mod==11) { write_evec_xsf("atominst12.xsf",11); }
    if (mod==12) { write_evec_xsf("atominst13.xsf",12); }
    if (mod==13) { write_evec_xsf("atominst14.xsf",13); }
    if (mod==14) { write_evec_xsf("atominst15.xsf",14); }
    if (mod==15) { write_evec_xsf("atominst16.xsf",15); }
    if (mod==16) { write_evec_xsf("atominst17.xsf",16); }
    if (mod==17) { write_evec_xsf("atominst18.xsf",17); }
    if (mod==18) { write_evec_xsf("atominst19.xsf",18); }
    if (mod==19) { write_evec_xsf("atominst20.xsf",19); }
  }

  std::cout<<"Instability analysis done."<<std::endl;
  potential();
  free(ar); free(e); free(ei); free(evi), free(evr), free(work);
  for (int i=1; i<=atom.nrepatom*3; i++) { free(hessian[i]); }; free(hessian);
}

void deform_element(int iel)
{
  double qx, qy, qz, helmat_p[3][3], helmat[3][3], hinelmat_p[3][3];
  mk_helmat(iel, atom.rx_p, atom.ry_p, atom.rz_p, helmat_p);
  mk_helmat(iel, atom.rx, atom.ry, atom.rz, helmat);
  inverse(helmat_p, hinelmat_p);
  for (int i=1; i<=atom.natom; i++) {
    if ((atom.repatom[i] == 0)&&(atom.elem_id[i] == iel)) { // if the atom is inside iel-th element
      qx = q_x(i, iel, atom.rx_p, atom.ry_p, atom.rz_p, hinelmat_p);
      qy = q_y(i, iel, atom.rx_p, atom.ry_p, atom.rz_p, hinelmat_p);
      qz = q_z(i, iel, atom.rx_p, atom.ry_p, atom.rz_p, hinelmat_p);
      atom.rx[i] = r_x(i, iel, atom.rx, atom.ry, atom.rz, helmat, qx, qy, qz);
      atom.ry[i] = r_y(i, iel, atom.rx, atom.ry, atom.rz, helmat, qx, qy, qz);
      atom.rz[i] = r_z(i, iel, atom.rx, atom.ry, atom.rz, helmat, qx, qy, qz);
    }
  }
}
