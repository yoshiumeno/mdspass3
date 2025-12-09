#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include "myheader.h"

#if !defined __linux__ && !defined __APPLE__
#define snprintf    _snprintf
#ifdef _DEBUG
#pragma comment(lib, "libpng16d.lib")
#else
#pragma comment(lib, "libpng16.lib")
#endif
#endif

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
extern "C"
void  dsyevr_(const char *JOBZp, const char *RANGEp, const char *UPLOp, int *Np,
	      double *A, int *LDAp, double *VLp, double *VUp,
	      int *ILp, int *IUp, double *ABSTOLp, int *Mp,
	      double *W, double *Z, int *LDZp, int *ISUPPZ,
	      double *WORK, int *LWORKp, int *IWORK, int *LIWORKp,
	      int *INFOp);

extern int hessian_read, hessian_write;
void potential();
void bookkeep();
void write_evec_xsf(const char* fname, int mod);
int atom_number(char* at);

void sort(int n, double *val, int *nod);
void write_evec_xsf(const char* fname, int mod);
void matcpy(double a[3][3], double b[3][3]);
void strain(double strmat[3][3]);
void calcstresstensor();

void instability_atom_noremove()
{
  std::cout<<"Instability analysis (atom only).."<<std::endl;
  int nn, info;
  int n = atom.natom;
  int iedg = 0;
  int icnt = atom.instcenter;

  //Wolfram
  // variables which are necessary for dsyevr call
  // il/iu: lower and higher index of eigenvalues to be calculated
  // health warning: we assume that il=1, else the allocation of
  // issupz must be changed!
  int il, iu;
  // output: array which stores eigenvectors
  double *evr;
  // output: supporting info for evr (cf, dsyevr docu)
  int **issupz;
  // output: how many ev. have been returned?
  int mret;
  // integer work array
  int *iwrk;
  // safe minimum for eigenvalue calculation
  double safmin;
  il=1;
  iu=MAXMODE;
  //safmin =DLAMCH( 'Safe minimum' )
  safmin = 1.17549e-38;

  //degree of freedom: atom
  //nn=n*3-3;
  nn=n*3;
  int nofpbc = 0;
  if (cell.pbcx) nofpbc++; if (cell.pbcy) nofpbc++; if (cell.pbcz) nofpbc++;
  /*
  if (nofpbc==1) {
    if (iedg==icnt) { nn=n*3-3+3; }
    if (iedg!=icnt) { nn=n*3-4+3; }
  }
  if (nofpbc==2) { nn=n*3-3+3; }
  if (nofpbc==3) { nn=n*3-3+6; }
  if (nofpbc==0) {
    printf("instability_atomcell.F: NO PBC is currently not supported.\n");
    return;
  }
  */

  // modified by Wolfram:
  // new allocations due to dsyevr call:
  //evr = new double*[nn];
  //for (int i=0; i<nn; i++) { evr[i] = new double[iu]; }
  evr = new double[nn*iu];
  //issupz = new int*[1];
  //for (int i=0; i<1; i++) { issupz[i] = new int[2*iu]; }
  iwrk = new int[nn*10];
  // end of new allocations


  double *ar   = new double[nn*nn];
  double *e    = new double[nn];
  double *ei   = new double[nn];
  double *work = new double[nn*61];

  double *rx_store = new double[atom.natom+1];
  double *ry_store = new double[atom.natom+1];
  double *rz_store = new double[atom.natom+1];
  int *nindx = new int[nn+1];

  //int *nod = (int *)calloc(nn, sizeof(int));

  // Allocate hessian matrix
  double **hessian;
  hessian = new double*[n*3];
  for (int i=0; i<n*3; i++) { hessian[i] = new double[n*3]; }
  /*
  double hessian_cell[6][6];
  double **hessian_ac;
  hessian_ac = new double*[n*3];
  for (int i=0; i<n*3; i++) { hessian_ac[i] = new double[6]; }
  */

  double hmat_store[3][3];
  for (int i=0; i<3; i++) { for (int j=0; j<3; j++) { hmat_store[i][j]=cell.hmat[i][j]; }}
  for (int i=1; i<=atom.natom; i++) {
    rx_store[i]=atom.rx[i]; ry_store[i]=atom.ry[i]; rz_store[i]=atom.rz[i]; 
  }

  double str = 0.01; // Infinitesimal strain
  double dis = 0.05*ang; // Infinitesimal displacement

  bookkeep();

  /*
  std::cout<<"Hessian (for cell) .."<<std::endl;
  double strmat[3][3];
  for (int i=0; i<6; i++) {
    matcpy(hmat_store, cell.hmat);
    for (int k=1; k<=atom.natom; k++) {
      atom.rx[k]=rx_store[k]; atom.ry[k]=ry_store[k]; atom.rz[k]=rz_store[k]; }
    double exx=0.0; double eyy=0.0; double ezz=0.0;
    double eyz=0.0; double ezx=0.0; double exy=0.0;
    if (i==0) { exx=str; } if (i==1) { eyy=str; } if (i==2) { ezz=str; }
    if (i==3) { eyz=str; } if (i==4) { ezx=str; } if (i==5) { exy=str; }
    strmat[0][0] = exx;   strmat[1][1] = eyy;   strmat[2][2] = ezz;
    strmat[1][2] = eyz/2; strmat[0][2] = ezx/2; strmat[0][1] = exy/2;
    strmat[2][1] = eyz/2; strmat[2][0] = ezx/2; strmat[1][0] = exy/2;
    strain(strmat);
    potential();
    calcstresstensor();
    for (int j=0; j<6; j++) {
      if (j==0) { hessian_cell[i][j] = cell.sgmmat[0][0]; }
      if (j==1) { hessian_cell[i][j] = cell.sgmmat[1][1]; }
      if (j==2) { hessian_cell[i][j] = cell.sgmmat[2][2]; }
      if (j==3) { hessian_cell[i][j] = cell.sgmmat[1][2]*2; }
      if (j==4) { hessian_cell[i][j] = cell.sgmmat[2][0]*2; }
      if (j==5) { hessian_cell[i][j] = cell.sgmmat[0][1]*2; }
    }
    for (int j=0; j<n; j++) {
      hessian_ac[j*3+0][i] = -atom.fx[j];
      hessian_ac[j*3+1][i] = -atom.fy[j];
      hessian_ac[j*3+2][i] = -atom.fz[j];
    }
    matcpy(hmat_store, cell.hmat);
    for (int k=1; k<=atom.natom; k++) {
      atom.rx[k]=rx_store[k]; atom.ry[k]=ry_store[k]; atom.rz[k]=rz_store[k]; }
    for (int k=0; k<3; k++) { for (int j=0; j<3; j++) { strmat[k][j] *= -1; } }
    strain(strmat);
    potential();
    calcstresstensor();
    for (int j=0; j<6; j++) {
      if (j==0) { hessian_cell[i][j] -= cell.sgmmat[0][0]; }
      if (j==1) { hessian_cell[i][j] -= cell.sgmmat[1][1]; }
      if (j==2) { hessian_cell[i][j] -= cell.sgmmat[2][2]; }
      if (j==3) { hessian_cell[i][j] -= cell.sgmmat[1][2]*2; }
      if (j==4) { hessian_cell[i][j] -= cell.sgmmat[2][0]*2; }
      if (j==5) { hessian_cell[i][j] -= cell.sgmmat[0][1]*2; }
    }
    for (int j=0; j<n; j++) {
      hessian_ac[j*3+0][i] += atom.fx[j];
      hessian_ac[j*3+1][i] += atom.fy[j];
      hessian_ac[j*3+2][i] += atom.fz[j];
    }
  }
  for (int i=0; i<3; i++) { for (int j=0; j<3; j++) { hessian_cell[i][j] /= str*2; } }
  */

  if (!hessian_read) {

  std::cout<<"Hessian calc (1/3).."<<std::endl;
  for (int i=1; i<=atom.natom; i++) {
    atom.rx[i]=atom.rx[i]+dis;
    potential();
    for (int j=1; j<=atom.natom; j++) {
      hessian[(i-1)*3+0][(j-1)*3+0] = atom.fx[j];
      hessian[(i-1)*3+0][(j-1)*3+1] = atom.fy[j];
      hessian[(i-1)*3+0][(j-1)*3+2] = atom.fz[j];
    }
    atom.rx[i]=atom.rx[i]-dis*2;
    potential();
    for (int j=1; j<=atom.natom; j++) {
      hessian[(i-1)*3+0][(j-1)*3+0] = (atom.fx[j]-hessian[(i-1)*3+0][(j-1)*3+0])/dis/2;
      hessian[(i-1)*3+0][(j-1)*3+1] = (atom.fy[j]-hessian[(i-1)*3+0][(j-1)*3+1])/dis/2;
      hessian[(i-1)*3+0][(j-1)*3+2] = (atom.fz[j]-hessian[(i-1)*3+0][(j-1)*3+2])/dis/2;
    }
    atom.rx[i]=atom.rx[i]+dis;
  }
  std::cout<<"Hessian calc (2/3).."<<std::endl;
  for (int i=1; i<=atom.natom; i++) {
    atom.ry[i]=atom.ry[i]+dis;
    potential();
    for (int j=1; j<=atom.natom; j++) {
      hessian[(i-1)*3+1][(j-1)*3+0] = atom.fx[j];
      hessian[(i-1)*3+1][(j-1)*3+1] = atom.fy[j];
      hessian[(i-1)*3+1][(j-1)*3+2] = atom.fz[j];
    }
    atom.ry[i]=atom.ry[i]-dis*2;
    potential();
    for (int j=1; j<=atom.natom; j++) {
      hessian[(i-1)*3+1][(j-1)*3+0] = (atom.fx[j]-hessian[(i-1)*3+1][(j-1)*3+0])/dis/2;
      hessian[(i-1)*3+1][(j-1)*3+1] = (atom.fy[j]-hessian[(i-1)*3+1][(j-1)*3+1])/dis/2;
      hessian[(i-1)*3+1][(j-1)*3+2] = (atom.fz[j]-hessian[(i-1)*3+1][(j-1)*3+2])/dis/2;
    }
    atom.ry[i]=atom.ry[i]+dis;
  }
  std::cout<<"Hessian calc (3/3).."<<std::endl;
  for (int i=1; i<=atom.natom; i++) {
    atom.rz[i]=atom.rz[i]+dis;
    potential();
    for (int j=1; j<=atom.natom; j++) {
      hessian[(i-1)*3+2][(j-1)*3+0] = atom.fx[j];
      hessian[(i-1)*3+2][(j-1)*3+1] = atom.fy[j];
      hessian[(i-1)*3+2][(j-1)*3+2] = atom.fz[j];
    }
    atom.rz[i]=atom.rz[i]-dis*2;
    potential();
    for (int j=1; j<=atom.natom; j++) {
      hessian[(i-1)*3+2][(j-1)*3+0] = (atom.fx[j]-hessian[(i-1)*3+2][(j-1)*3+0])/dis/2;
      hessian[(i-1)*3+2][(j-1)*3+1] = (atom.fy[j]-hessian[(i-1)*3+2][(j-1)*3+1])/dis/2;
      hessian[(i-1)*3+2][(j-1)*3+2] = (atom.fz[j]-hessian[(i-1)*3+2][(j-1)*3+2])/dis/2;
    }
    atom.rz[i]=atom.rz[i]+dis;
  }

  if (hessian_write) {
    std::ofstream ofs("hessian.d");
    for (int i=0; i<atom.natom*3; i++) {
      for (int j=0; j<atom.natom*3; j++) {
	ofs << hessian[i][j] << std::endl;
      }
    }
  }

  } else { // if hessian_read

  std::ifstream ifs( "hessian.d" );
  for (int i=0; i<atom.natom*3; i++) {
    for (int j=0; j<atom.natom*3; j++) {
      ifs >> hessian[i][j];
    }
  }

  } // endif hessian_read

  /*
  double emax = -1.0e100;
  for (int i=0; i<3*n; i++) { for (int j=0; j<3*n; j++) {
      if (emax < fabs(hessian[i][j])) { emax = fabs(hessian[i][j]); }
    } }
  for (int i=0; i<3*n; i++) { for (int j=0; j<3*n; j++) {
      hessian[i][j] /= emax;
    } }
  */
  /*
  emax = -1.0e100;
  for (int i=0; i<6; i++) { for (int j=0; j<6; j++) {
      if (emax < fabs(hessian_cell[i][j])) { emax = fabs(hessian_cell[i][j]); }
    } }
  for (int i=0; i<6; i++) { for (int j=0; j<6; j++) {
      hessian_cell[i][j] /= emax;
    } }
  emax = -1.0e100;
  for (int i=0; i<3*n; i++) { for (int j=0; j<6; j++) {
      if (emax < fabs(hessian_ac[i][j])) { emax = fabs(hessian_ac[i][j]); }
    } }
  for (int i=0; i<3*n; i++) { for (int j=0; j<6; j++) {
      hessian_ac[i][j] /= emax;
    } }
  */

  //search one atom at edge and one near cell center
  if (iedg == 0) {
    double xmin=1.0e10; double ymin=1.0e10; double zmin=1.0e10; 
    for (int i=1; i<=n; i++) {
      if (atom.rx[i] < xmin) { xmin=atom.rx[i]; iedg=i; }
      if (atom.ry[i] < ymin) { ymin=atom.ry[i]; iedg=i; }
      if (atom.rz[i] < zmin) { zmin=atom.rz[i]; iedg=i; }
    }
    double xcnt = (cell.hmat[0][0]+cell.hmat[1][0]+cell.hmat[2][0])/2.0;
    double ycnt = (cell.hmat[0][1]+cell.hmat[1][1]+cell.hmat[2][1])/2.0;
    double zcnt = (cell.hmat[0][2]+cell.hmat[1][2]+cell.hmat[2][2])/2.0;
    double rng = 1.0e-11;
  POINT:
    icnt = 0;
    for (int i=1; i<=n; i++) {
      if ((fabs(atom.rx[i]-xcnt)<rng)&&(fabs(atom.ry[i]-ycnt)<rng)
	  &&(fabs(atom.rz[i]-zcnt)<rng)) {
	icnt=i; break;
      }
    }
    if (icnt == 0) {rng *= 2; goto POINT; }
  }
  atom.instcenter=icnt;
  printf("Edge atom = %d (%f %f %f [ang])\n",
	 iedg,(atom.rx[iedg]/ang),(atom.ry[iedg]/ang),(atom.rz[iedg]/ang));
  printf("Cent atom = %d (%f %f %f [ang])\n",
	 icnt,(atom.rx[icnt]/ang),(atom.ry[icnt]/ang),(atom.rz[icnt]/ang));

  // --- Matrix set (removing some lines to meet DoF) ---
  /*
  int nrmv, irmv[7];
  if (nofpbc==3) {      // --> atom icnt is fixed
    nrmv=3; irmv[1]=3*(icnt-1)+1; irmv[2]=3*(icnt-1)+2; irmv[3]=3*(icnt-1)+3;
  } else if (nofpbc==2) { // --> atom icnt is fixed 
    nrmv=3; irmv[1]=3*(icnt-1)+1; irmv[2]=3*(icnt-1)+2; irmv[3]=3*(icnt-1)+3;
  } else if (nofpbc==1) { // --> atom icnt is fixed, iedg is constrained
    printf("instability atomcell: one-direction PBC is currently not supported.\n");
    if (iedg!=icnt) {
      nrmv=4; irmv[1]=3*(icnt-1)+1; irmv[2]=3*(icnt-1)+2; irmv[3]=3*(icnt-1)+3;
      if (cell.pbcx) { irmv[4]=3*(iedg-1); }
      if (cell.pbcy) { irmv[4]=3*(iedg-1)+1; }
      if (cell.pbcz) { irmv[4]=3*(iedg-1)+2; }
    } else {
      nrmv=3; irmv[1]=3*(icnt-1)+1; irmv[2]=3*(icnt-1)+2; irmv[3]=3*(icnt-1)+3;
    }
  }
  */

  // Copy hessian to ar.
  int indx, jndx, dummy;
  indx = 0;
  for (int i=1; i<=n*3; i++) {
    //for (int ii=1; ii<=nrmv; ii++) {
      //if (i==irmv[ii]) { goto P990; }
    //}
    indx++;
    nindx[indx]=i;

    jndx = 0;
    for (int j=1; j<=n*3; j++) {
      //for (int ii=1; ii<=nrmv; ii++) {
	//if (j==irmv[ii]) { goto P991; }
      //}
      jndx++;
      
      //ar[indx-1][jndx-1]=hessian[i-1][j-1];
      ar[indx-1 + (jndx-1)*nn]=hessian[i-1][j-1];
      
    P991:
      dummy=0;
    }
    
  P990:
    dummy=0;
  }

  /*

  int nhes = n*3-nrmv;
  for (int i=0; i<nhes; i++) {
    if (nofpbc==3) {
      for (int j=0; j<6; j++) {
	ar[i + (nhes+j)*nn]=hessian_ac[i][j];
	ar[nhes+j + (i)*nn]=hessian_ac[i][j];
      }
    } else if (nofpbc==2) {
      if (!cell.pbcz) {
	ar[i + (nhes+0)*nn]=hessian_ac[i][0]; ar[nhes+0 + i*nn]=hessian_ac[i][0];
	ar[i + (nhes+1)*nn]=hessian_ac[i][1]; ar[nhes+1 + i*nn]=hessian_ac[i][1];
	ar[i + (nhes+2)*nn]=hessian_ac[i][5]; ar[nhes+2 + i*nn]=hessian_ac[i][5];
      } else if (!cell.pbcx) {
	ar[i + (nhes+0)*nn]=hessian_ac[i][1]; ar[nhes+0 + i*nn]=hessian_ac[i][1];
	ar[i + (nhes+1)*nn]=hessian_ac[i][2]; ar[nhes+1 + i*nn]=hessian_ac[i][2];
	ar[i + (nhes+2)*nn]=hessian_ac[i][3]; ar[nhes+2 + i*nn]=hessian_ac[i][3];
      } else if (!cell.pbcy) {
	ar[i + (nhes+0)*nn]=hessian_ac[i][0]; ar[nhes+0 + i*nn]=hessian_ac[i][0];
	ar[i + (nhes+1)*nn]=hessian_ac[i][2]; ar[nhes+1 + i*nn]=hessian_ac[i][2];
	ar[i + (nhes+2)*nn]=hessian_ac[i][4]; ar[nhes+2 + i*nn]=hessian_ac[i][4];
      }
    }
  }
  if (nofpbc==3) {
    for (int i=0; i<6; i++) {
      for (int j=0; j<6; j++) {
	ar[nhes+i + (nhes+j)*nn]=hessian_cell[i][j];
      }
    }
  } else if (nofpbc==2) {
    if (!cell.pbcz) {
      ar[nhes+0 + (nhes+0)*nn]=hessian_cell[0][0];
      ar[nhes+0 + (nhes+1)*nn]=hessian_cell[0][1];
      ar[nhes+0 + (nhes+2)*nn]=hessian_cell[0][5];
      ar[nhes+1 + (nhes+0)*nn]=hessian_cell[1][0];
      ar[nhes+1 + (nhes+1)*nn]=hessian_cell[1][1];
      ar[nhes+1 + (nhes+2)*nn]=hessian_cell[1][5];
      ar[nhes+2 + (nhes+0)*nn]=hessian_cell[5][0];
      ar[nhes+2 + (nhes+1)*nn]=hessian_cell[5][1];
      ar[nhes+2 + (nhes+2)*nn]=hessian_cell[5][5];
    } else if (!cell.pbcx) {
      ar[nhes+0 + (nhes+0)*nn]=hessian_cell[1][1];
      ar[nhes+0 + (nhes+1)*nn]=hessian_cell[1][2];
      ar[nhes+0 + (nhes+2)*nn]=hessian_cell[1][3];
      ar[nhes+1 + (nhes+0)*nn]=hessian_cell[2][1];
      ar[nhes+1 + (nhes+1)*nn]=hessian_cell[2][2];
      ar[nhes+1 + (nhes+2)*nn]=hessian_cell[2][3];
      ar[nhes+2 + (nhes+0)*nn]=hessian_cell[3][1];
      ar[nhes+2 + (nhes+1)*nn]=hessian_cell[3][2];
      ar[nhes+2 + (nhes+2)*nn]=hessian_cell[3][3];
    } else if (!cell.pbcy) {
      ar[nhes+0 + (nhes+0)*nn]=hessian_cell[0][0];
      ar[nhes+0 + (nhes+1)*nn]=hessian_cell[0][2];
      ar[nhes+0 + (nhes+2)*nn]=hessian_cell[0][4];
      ar[nhes+1 + (nhes+0)*nn]=hessian_cell[2][0];
      ar[nhes+1 + (nhes+1)*nn]=hessian_cell[2][2];
      ar[nhes+1 + (nhes+2)*nn]=hessian_cell[2][4];
      ar[nhes+2 + (nhes+0)*nn]=hessian_cell[4][0];
      ar[nhes+2 + (nhes+1)*nn]=hessian_cell[4][2];
      ar[nhes+2 + (nhes+2)*nn]=hessian_cell[4][4];
    }
  }

  */

  int noc = nn;

  printf("Solving eigenvalue problem for atom+cell matrix..\n");
  printf("Matrix size nn: %d\n",nn);
  printf("error tolerance: %e\n",safmin);

  std::cout<<"Solving eigenvalue problem."<<std::endl;

  //sgeev_("N","V",&nn,ar,&nn,e,ei,evi,&nn,evr,&nn,work,&nn4,&info);
  int isuppz, nn61=nn*61, nn10=nn*10;
  double vlp=-1.0, vup=1.0;
  dsyevr_("V","I","U",&nn,ar,&nn,&vlp,&vup,&il,&iu,&safmin,&mret,e,evr,&nn,&isuppz,work,&nn61,iwrk,&nn10,&info);
 
  int *nod = new int[iu];
  for (int i=0; i<iu; i++) { nod[i]=i; }
  sort(iu, e, nod);

  //  for (int i=0; i<nn; i++) { printf("%d %20.15e\n",i,e[i]); }
  FILE *outputfile = fopen("eminall.d","w");
  for (int i=0; i<iu; i++) { fprintf (outputfile, "%d %20.15e\n",i,e[i]); atom.eigval[i] = e[i]; }
  fclose(outputfile);
  std::cout<<"Eigenvalues written in eminall.d"<<std::endl;

  for (int mod=0; mod<MAXMODE; mod++) {
    int imod = nod[mod];
    for (int i=1; i<=atom.natom; i++) {
      //if (i!=atom.instcenter) {
        int ii=0; ii=(i-1)*3; //if (i<atom.instcenter) { ii=(i-1)*3; } else { ii=(i-2)*3; }
	ii=ii+nn*imod;
	atom.evecx[i][mod]=evr[ii]; atom.evecy[i][mod]=evr[ii+1]; atom.evecz[i][mod]=evr[ii+2];
	//}
    }
    char filepath[80] = "atominst";
    char numc[3];
    snprintf(numc, sizeof(numc), "%02d", mod+1);
    strcat(filepath,numc); strcat(filepath,".xsf");
    write_evec_xsf(filepath, mod);
  }

  //for (int i=0; i<nn*nn; i++) {
  //  printf("%d %20.15e\n",i,ar[i]);
  //}
  //  printf("%20.15f\n",e[nod[0]]);printf("%20.15f\n",e[nod[2]]);printf("%20.15f\n",e[nod[2]]);

  std::cout<<"Instability analysis done."<<std::endl;
  potential();

  delete[] evr; delete[] iwrk; delete[] ar; delete[] e; delete[] ei; delete[] work;
  delete[] rx_store; delete[] ry_store; delete[] rz_store, nindx;
  for (int i=0; i<n*3; i++) { delete[] hessian[i]; }; delete[] hessian;
  /*
  for (int i=0; i<n*3; i++) { delete[] hessian_ac[i]; }; delete[] hessian_ac;
  */
  delete[] nod;
}

