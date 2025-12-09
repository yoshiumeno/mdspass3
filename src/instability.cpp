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

extern int hessian_read, hessian_write;
void potential();
void bookkeep();
void write_evec_xsf(const char* fname, int mod);
int atom_number(char* at);

void sort(int n, double *val, int *nod)
{
  double tmp;
  int itmp;

  for (int i=0; i<n; i++) { nod[i] = i; }

  for (int i=0; i<n-1; i++) {
    for (int j=i+1; j<n; j++) {
      if (val[i]>val[j]) {
	tmp = val[i]; val[i] = val[j]; val[j] = tmp;
	itmp = nod[i]; nod[i] = nod[j]; nod[j] = itmp;
      }
    }
  }
}

void instability()
{
  std::cout<<"Instability analysis.."<<std::endl;
  long int nn, nn4, info;
  double dis = 0.05e-10; // Infinitesimal displacement
  //  double dis = 0.0005e-10; // Infinitesimal displacement
  double *rx_store = (double *)calloc(atom.natom+1, sizeof(double));
  double *ry_store = (double *)calloc(atom.natom+1, sizeof(double));
  double *rz_store = (double *)calloc(atom.natom+1, sizeof(double));

  nn=atom.natom*3-3; nn4=nn*4;
  //  nn=903; nn4=nn*4;
  double *ar = (double *)calloc(nn*nn, sizeof(double)); 
  //  double *aro = (double *)calloc(nn*nn, sizeof(double)); 
  double *e  = (double *)calloc(nn, sizeof(double)); 
  double *ei = (double *)calloc(nn, sizeof(double)); 
  double *evi = (double *)calloc(nn*nn, sizeof(double)); 
  double *evr = (double *)calloc(nn*nn, sizeof(double)); 
  double *work = (double *)calloc(nn4, sizeof(double)); 

  int *nod = (int *)calloc(nn, sizeof(int));

  // Allocate hessian matrix
  double **hessian;
  hessian = (double**)malloc(sizeof(double*)*(atom.natom*3+1));
  for (int i=1; i<=atom.natom*3; i++) {
    hessian[i] = (double *)malloc(sizeof(double)*(atom.natom*3+1));
  }
    
  for (int i=1; i<=atom.natom; i++) {
    rx_store[i]=atom.rx[i]; ry_store[i]=atom.ry[i]; rz_store[i]=atom.rz[i]; 
  }

  if (!hessian_read) {

  bookkeep();
  //  atom.rx[1]=atom.rx[1]+dis;
  //  e_force_morse();
  //  printf("%20.15e\n",atom.rx[1]); printf("%20.15e\n",atom.fx[1]); exit(0);
  std::cout<<"Hessian calc (1/3).."<<std::endl;
  for (int i=1; i<=atom.natom; i++) {
    atom.rx[i]=atom.rx[i]+dis;
    potential();
    for (int j=1; j<=atom.natom; j++) {
      hessian[(i-1)*3+1][(j-1)*3+1] = atom.fx[j];
      hessian[(i-1)*3+1][(j-1)*3+2] = atom.fy[j];
      hessian[(i-1)*3+1][(j-1)*3+3] = atom.fz[j];
    }
    atom.rx[i]=atom.rx[i]-dis*2;
    potential();
    for (int j=1; j<=atom.natom; j++) {
      hessian[(i-1)*3+1][(j-1)*3+1] = (atom.fx[j]-hessian[(i-1)*3+1][(j-1)*3+1])/dis/2;
      hessian[(i-1)*3+1][(j-1)*3+2] = (atom.fy[j]-hessian[(i-1)*3+1][(j-1)*3+2])/dis/2;
      hessian[(i-1)*3+1][(j-1)*3+3] = (atom.fz[j]-hessian[(i-1)*3+1][(j-1)*3+3])/dis/2;
    }
    atom.rx[i]=atom.rx[i]+dis;
    //    printf("%20.15e\n",hessian[1][1]);exit(0);
  }
  std::cout<<"Hessian calc (2/3).."<<std::endl;
  for (int i=1; i<=atom.natom; i++) {
    atom.ry[i]=atom.ry[i]+dis;
    potential();
    for (int j=1; j<=atom.natom; j++) {
      hessian[(i-1)*3+2][(j-1)*3+1] = atom.fx[j];
      hessian[(i-1)*3+2][(j-1)*3+2] = atom.fy[j];
      hessian[(i-1)*3+2][(j-1)*3+3] = atom.fz[j];
    }
    atom.ry[i]=atom.ry[i]-dis*2;
    potential();
    for (int j=1; j<=atom.natom; j++) {
      hessian[(i-1)*3+2][(j-1)*3+1] = (atom.fx[j]-hessian[(i-1)*3+2][(j-1)*3+1])/dis/2;
      hessian[(i-1)*3+2][(j-1)*3+2] = (atom.fy[j]-hessian[(i-1)*3+2][(j-1)*3+2])/dis/2;
      hessian[(i-1)*3+2][(j-1)*3+3] = (atom.fz[j]-hessian[(i-1)*3+2][(j-1)*3+3])/dis/2;
    }
    atom.ry[i]=atom.ry[i]+dis;
  }
  std::cout<<"Hessian calc (3/3).."<<std::endl;
  for (int i=1; i<=atom.natom; i++) {
    atom.rz[i]=atom.rz[i]+dis;
    potential();
    for (int j=1; j<=atom.natom; j++) {
      hessian[(i-1)*3+3][(j-1)*3+1] = atom.fx[j];
      hessian[(i-1)*3+3][(j-1)*3+2] = atom.fy[j];
      hessian[(i-1)*3+3][(j-1)*3+3] = atom.fz[j];
    }
    atom.rz[i]=atom.rz[i]-dis*2;
    potential();
    for (int j=1; j<=atom.natom; j++) {
      hessian[(i-1)*3+3][(j-1)*3+1] = (atom.fx[j]-hessian[(i-1)*3+3][(j-1)*3+1])/dis/2;
      hessian[(i-1)*3+3][(j-1)*3+2] = (atom.fy[j]-hessian[(i-1)*3+3][(j-1)*3+2])/dis/2;
      hessian[(i-1)*3+3][(j-1)*3+3] = (atom.fz[j]-hessian[(i-1)*3+3][(j-1)*3+3])/dis/2;
    }
    atom.rz[i]=atom.rz[i]+dis;
  }

  if (hessian_write) {
    std::ofstream ofs("hessian.d");
    for (int i=1; i<=atom.natom*3; i++) {
      for (int j=1; j<=atom.natom*3; j++) {
	ofs << hessian[i][j] << std::endl;
    }
    }
  }

  } else { // if hessian_read

  std::ifstream ifs( "hessian.d" );
  for (int i=1; i<=atom.natom*3; i++) {
    for (int j=1; j<=atom.natom*3; j++) {
      ifs >> hessian[i][j];
    }
  }

  } // endif hessian_read

  // Copy hessian to ar.
  std::cout<<atom.instcenter<<" -th atom is fixed."<<std::endl;
  int ii=0;
  for (int i=1; i<=atom.natom; i++) {
    if (i!=atom.instcenter) {
      for (int j=1; j<=atom.natom; j++) {
	if (j!=atom.instcenter) {
	  ii=0;
	  if (j<atom.instcenter) { ii=ii+(j-1)*3; } else { ii=ii+(j-2)*3; }
	  if (i<atom.instcenter) { ii=ii+(i-1)*(atom.natom-1)*9; } else { ii=ii+(i-2)*(atom.natom-1)*9; }
	  ar[ii]=hessian[(i-1)*3+1][(j-1)*3+1]; ii++;
	  ar[ii]=hessian[(i-1)*3+1][(j-1)*3+2]; ii++;
	  ar[ii]=hessian[(i-1)*3+1][(j-1)*3+3]; ii++; ii=ii+nn-3;
	  ar[ii]=hessian[(i-1)*3+2][(j-1)*3+1]; ii++;
	  ar[ii]=hessian[(i-1)*3+2][(j-1)*3+2]; ii++;
	  ar[ii]=hessian[(i-1)*3+2][(j-1)*3+3]; ii++; ii=ii+nn-3;
	  ar[ii]=hessian[(i-1)*3+3][(j-1)*3+1]; ii++;
	  ar[ii]=hessian[(i-1)*3+3][(j-1)*3+2]; ii++;
	  ar[ii]=hessian[(i-1)*3+3][(j-1)*3+3]; ii++; ii=ii+nn-3;
	}
      }
    }
  }


  //for (int i=0; i<nn*nn; i++) {
  //  printf("%d %20.15e\n",i,ar[i]);
  //}

  /*
  std::ifstream ifs( "ar.d" );
  for (int i=0; i<nn*nn; i++) {
    ifs >> ar[i];
  }
  */
  /*
  for (int i=0; i<nn*nn; i++) {
    printf("%d %20.15e\n",i,ar[i]);
  }
  */

  /*
  std::ofstream ofs("ar.d");
  for (int i=0; i<nn*nn; i++) {
    ofs << ar[i] << std::endl;
  }
  */

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

  //  for (int i=0; i<nn; i++) { printf("%d %20.15e\n",i,e[i]); }
  FILE *outputfile = fopen("eminall.d","w");
  for (int i=0; i<nn; i++) { fprintf (outputfile, "%d %20.15e\n",i,e[i]); }
  fclose(outputfile);
  std::cout<<"Eigenvalues written in eminall.d"<<std::endl;


  for (int mod=0; mod<20; mod++) {
    int imod = nod[mod];
    for (int i=1; i<=atom.natom; i++) {
      if (i!=atom.instcenter) {
	int ii=0; if (i<atom.instcenter) { ii=(i-1)*3; } else { ii=(i-2)*3; }
	ii=ii+nn*imod;
	atom.evecx[i][mod]=evr[ii]; atom.evecy[i][mod]=evr[ii+1]; atom.evecz[i][mod]=evr[ii+2];
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

  //for (int i=0; i<nn*nn; i++) {
  //  printf("%d %20.15e\n",i,ar[i]);
  //}
  //  printf("%20.15f\n",e[nod[0]]);printf("%20.15f\n",e[nod[2]]);printf("%20.15f\n",e[nod[2]]);

  std::cout<<"Instability analysis done."<<std::endl;
  potential();
  free(rx_store); free(ry_store); free(rz_store);
  free(ar); free(e); free(ei); free(evi), free(evr), free(work);
  for (int i=1; i<=atom.natom*3; i++) { free(hessian[i]); }; free(hessian);
}

void write_evec_xsf(const char* fname, int mod)
{
  FILE *fp = fopen(fname,"w");
  fprintf(fp,"CRYSTAL\n");
  fprintf(fp,"PRIMVEC\n");
  fprintf(fp,"%10.5f %10.5f %10.5f\n",
	  cell.hmat[0][0]/ang, cell.hmat[1][0]/ang, cell.hmat[2][0]/ang);
  fprintf(fp,"%10.5f %10.5f %10.5f\n",
	  cell.hmat[0][1]/ang, cell.hmat[1][1]/ang, cell.hmat[2][1]/ang);
  fprintf(fp,"%10.5f %10.5f %10.5f\n",
	  cell.hmat[0][2]/ang, cell.hmat[1][2]/ang, cell.hmat[2][2]/ang);
  fprintf(fp,"PRIMCOORD\n");
  fprintf(fp,"%d 1\n",atom.natom);
  for (int i=1; i<=atom.natom; i++) {
    fprintf(fp,"%d %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
	    atom_number(atom.asp[i]),atom.rx[i]/ang,atom.ry[i]/ang,atom.rz[i]/ang,
	    atom.evecx[i][mod],atom.evecy[i][mod],atom.evecz[i][mod]);
  }
  fclose(fp);
}
