#include <iostream>
#include <stdio.h>
#include <fstream>
#include "myheader.h"

void potential();
void inverse(double a[3][3], double b[3][3]);
void matmul(double a[3][3], double b[3][3], double c[3][3]);
void multiply_cell(int ix, int iy, int iz, int mode);
void divide_cell(int ix, int iy, int iz, int mode);
void forceconst(double **fcmat);
void kgen_fcc(double *kx, double *ky, double *kz, int phonon_knum, int &size);
void kgen_bcc(double *kx, double *ky, double *kz, int phonon_knum, int &size);
void kgen_hcp(double *kx, double *ky, double *kz, int phonon_knum, int &size);
void kgen_sc(double *kx, double *ky, double *kz, int phonon_knum, int &size);
void kgen_chain(double *kx, double *ky, double *kz, int phonon_knum, int &size);
void kgen_fcc_dense(double *kx, double *ky, double *kz, double *fkp, int phonon_knum, int &size);
void phonon_calc(int mode);

struct _dcomplex { double re, im; };
typedef struct _dcomplex dcomplex;
extern "C"
void zgeev_( const char* jobvl, const char* jobvr, int* n, dcomplex* a,
	     int* lda, dcomplex* w, dcomplex* vl, int* ldvl, dcomplex* vr, int* ldvr,
	     dcomplex* work, int* lwork, double* rwork, int* info);

extern int phonon_rep, phonon_knum, phonon_kp;
extern float dos_gauss_width;
extern int dos_kmesh;

double gauss(double x, double a, double c)
{
  return a * exp(-x*x/2.0/(c*c));
}

void phonon_calc()
{
  phonon_calc(0);
}


void phonon_calc(int mode)
{ // mode: 0 = dispersion curve, 1 = dos
  //multiply cell size
  int n0 = atom.natom*3;
  int ix = phonon_rep, iy = phonon_rep, iz = phonon_rep;
  multiply_cell(ix, iy, iz, 1);
  //allocate array
  double **fcmat;
  dcomplex *dm;
  int nn = atom.natom*3;
  fcmat = new double*[nn];
  dm = new dcomplex[n0*n0];
  for (int i=0; i<nn; i++) {
    fcmat[i] = new double[nn];
  }
  int nn2 = n0*2, info = 1, ldvl = 1, ldvr = 1;
  dcomplex *work; work = new dcomplex[nn2];
  double *rwork; rwork = new double[nn2];
  dcomplex *w;
  dcomplex *vl, *vr;
  w = new dcomplex[n0];
  vl = new dcomplex[n0*ldvl]; vr = new dcomplex[nn*ldvr];

  // for DoS
  double a, c, range;
  int mesh;
  double *dos;

  forceconst(fcmat);
  std::ofstream ofs;
  if (mode == 0) { // E-k mode
    ofs.open("phonon.d");
    printf("Phonon dispersion data is written in phonon.d\n");
  } else { // Dos mode
    ofs.open("dos.d");
    printf("Density of states data is written in dos.d\n");

    a = 1.0, c = dos_gauss_width; // Amplitude and width of gauss function
    range = c * 5.0;
    mesh = 10000;
    dos = new double[mesh];
    for (int i=0; i<mesh; i++) { dos[i] = 0; }
  }

  double *kx, *ky, *kz, *fkp, *eval;

  int size = 0;
  if (phonon_kp == 0) {
    if (mode == 0) {
      kgen_fcc(kx, ky, kz, phonon_knum, size);
    } else {
      kgen_fcc_dense(kx, ky, kz, fkp, phonon_knum, size);
    }
  } else if (phonon_kp == 1) {
    kgen_bcc(kx, ky, kz, phonon_knum, size);
  } else if (phonon_kp == 2) {
    kgen_sc(kx, ky, kz, phonon_knum, size);
  } else if (phonon_kp == 3) {
    kgen_hcp(kx, ky, kz, phonon_knum, size);
  } else if (phonon_kp == 4) {
    kgen_chain(kx, ky, kz, phonon_knum, size);
  }
  printf("# of k-points = %d\n",size);
  kx = new double[size];
  ky = new double[size];
  kz = new double[size];
  if (mode == 1) {
    fkp = new double[size];
    eval = new double[size*n0];
  }
  if (phonon_kp == 0) {
    if (mode == 0) {
      kgen_fcc(kx, ky, kz, phonon_knum, size);
    } else {
      kgen_fcc_dense(kx, ky, kz, fkp, phonon_knum, size);
    }
  } else if (phonon_kp == 1) {
    kgen_bcc(kx, ky, kz, phonon_knum, size);
  } else if (phonon_kp == 2) {
    kgen_sc(kx, ky, kz, phonon_knum, size);
  } else if (phonon_kp == 3) {
    kgen_hcp(kx, ky, kz, phonon_knum, size);
  } else if (phonon_kp == 4) {
    kgen_chain(kx, ky, kz, phonon_knum, size);
  }

  int kk=0;
  for (int ik=0; ik<size; ik++) {

    for (int i=0; i<n0*n0; i++) {
      dm[i].re = 0; dm[i].im = 0; }
    
    for (int i=0; i<n0; i++) {
      int ia = i/3 + 1;
      for (int j=0; j<nn; j++) {
	int ja = j/3 + 1;
	double dx = atom.rx[ia] - atom.rx[ja];
	double dy = atom.ry[ia] - atom.ry[ja];
	double dz = atom.rz[ia] - atom.rz[ja];
	double xx = - (kx[ik]*dx + ky[ik]*dy + kz[ik]*dz);
	double ms = sqrt(atom.wm[ia]*atom.wm[ja]);
	double xxr = fcmat[i][j] * cos(xx) / ms;
	double xxi = fcmat[i][j] * sin(xx) / ms;
	int ii = i+(j%n0)*n0;
	dm[ii].re += xxr; dm[ii].im += xxi;
      }
    }
    zgeev_("N", "N", &n0, dm, &n0, w, vl, &ldvl, vr, &ldvr, work, &nn2, rwork, &info);
    for (int i=0; i<n0; i++) {
      if (w[i].re < 0) { w[i].re = 0; }
      else { w[i].re = sqrt(w[i].re); }
    }

    //double scale = 1e-27;
    double scale = 1e-13;
    if (mode == 0) { // for E-k mode, eigenvalues are simply written to file
      for (int i=0; i<n0; i++) {
	ofs << ik << " " << w[i].re * scale << std::endl;
      }
      ofs << std::endl;
    } else { // for DoS mode, eigenvalues are stored in eval[]
      for (int i=0; i<n0; i++) {
	eval[kk] = w[i].re * scale;
	kk++;
      }
    }
  }

  // Integral for DoS 
  if (mode == 1) {
    double emax, emin;
    for (int i=0; i<kk; i++) {
      if (emax < eval[i]) { emax = eval[i]; }
      if (emin > eval[i]) { emin = eval[i]; }
    }
    emin -= range; emax += range;
    printf("E_min, E_max: %f , %f\n",(float)emin, (float)emax);
    double dx = (emax-emin)/(double)(mesh-1);
    int i0, i1;
    int i=0;
    for (int ik=0; ik<size; ik++) {
      for (int ie=0; ie<n0; ie++) {
	i0 = int(((eval[i]-range)-emin)/dx)+1;
	i1 = int(((eval[i]+range)-emin)/dx)+1;
	for (int j=i0; j<=i1; j++) {
	  double y = -range + dx * (double)(j-i0);
	  dos[j] += gauss(y,a*fkp[ik],c);
	}
	i++;
      }
    }
    for (int i=0; i<mesh; i++) {
      double x = emin + (double)i * dx;
      ofs << x << " " << dos[i] << std::endl;
    }
  }

  //deallocate array
  for (int i=0; i<nn; i++) {
    delete[] fcmat[i];
  }
  delete[] fcmat; delete[] dm; delete[] work; delete[] rwork;
  delete[] w; delete[] vl; delete[] vr;
  delete[] kx; delete[] ky; delete[] kz;
  if (mode == 1) { delete[] fkp; delete[] dos; delete[] eval; }
  ofs.close();

  //Divide cell
  divide_cell(ix, iy, iz, 1);

}

void kgen_fcc(double *kx, double *ky, double *kz, int knum, int &size)
{
  double len[5] = {1,1,1,0.5,1}; // ratio of k-point density
  double spx[10],spy[10],spz[10];
  double a;
  if (cell.hmat[1][0]/cell.hmat[0][0] > 0.5) {
    a = cell.hmat[0][0]/(double)(phonon_rep*2-1)*2;
  } else {
    a = cell.hmat[0][0]/(double)(phonon_rep*2-1);
  }
  const char *sp[] = {"W", "L", "L", "Gamma", "Gamma", "X", "X", "U", "K", "Gamma"};
  printf("lattice constant a = %f (A)\n",a/ang);
  spx[0] = 1.00; spy[0] = 0.00; spz[0] = 0.50; // W
  spx[1] = 0.50; spy[1] = 0.50; spz[1] = 0.50; // L
  spx[2] = 0.50; spy[2] = 0.50; spz[2] = 0.50; // L
  spx[3] = 0.00; spy[3] = 0.00; spz[3] = 0.00; // Gamma
  spx[4] = 0.00; spy[4] = 0.00; spz[4] = 0.00; // Gamma
  spx[5] = 1.00; spy[5] = 0.00; spz[5] = 0.00; // X
  spx[6] = 1.00; spy[6] = 0.00; spz[6] = 0.00; // X
  spx[7] = 1.00; spy[7] = 0.25; spz[7] = 0.25; // U
  spx[8] = 0.75; spy[8] = 0.00; spz[8] = 0.75; // K
  spx[9] = 0.00; spy[9] = 0.00; spz[9] = 0.00; // Gamma
  bool count = true;
  if (size > 0) count = false;
  size=0;
  for (int i=0; i<10; i+=2) {
    int knumr = knum*len[i/2];
    for (int k=0; k<knumr; k++) {
      if (!count) {
	double x = (double)k/(double)knumr;
	kx[size] = 2*M_PI/a* (spx[i] + (spx[i+1]-spx[i])*x );
	ky[size] = 2*M_PI/a* (spy[i] + (spy[i+1]-spy[i])*x );
	kz[size] = 2*M_PI/a* (spz[i] + (spz[i+1]-spz[i])*x );
	if (k==0) {
	  printf("%4d - %4d: %s - %s\n",size,size+knumr,sp[i],sp[i+1]);
	}
	//printf("%d %d %f %f %f\n",size,i,kx[size]*a/2/M_PI,ky[size]*a/2/M_PI,kz[size]*a/2/M_PI);
      }
      size++;
    }
  }
  if (!count) {
    kx[size]=spx[9]; ky[size]=spy[9]; kz[size]=spz[9];
  }
  size++;
}

void kgen_bcc(double *kx, double *ky, double *kz, int knum, int &size)
{
  double len[5] = {1,1,0.5,1,1}; // ratio of k-point density
  double spx[10],spy[10],spz[10];
  double a;
  if (fabs(cell.hmat[0][1]/cell.hmat[0][0]) > 0.5) {
    a = cell.hmat[0][0]/(double)(phonon_rep*2-1)*2;
  } else {
    a = cell.hmat[0][0]/(double)(phonon_rep*2-1);
  }
  const char *sp[] = {"Gamma", "H", "H", "N", "N", "P", "P", "Gamma", "Gamma", "N"};
  printf("lattice constant a = %e\n",a);
  spx[0] = 0.00; spy[0] = 0.00; spz[0] = 0.00; // Gamma
  spx[1] = 1.00; spy[1] = 0.00; spz[1] = 0.00; // H
  spx[2] = 1.00; spy[2] = 0.00; spz[2] = 0.00; // H
  spx[3] = 0.50; spy[3] = 0.00; spz[3] = 0.50; // N
  spx[4] = 0.50; spy[4] = 0.00; spz[4] = 0.50; // N
  spx[5] = 0.50; spy[5] = 0.50; spz[5] = 0.50; // P
  spx[6] = 0.50; spy[6] = 0.50; spz[6] = 0.50; // P
  spx[7] = 0.00; spy[7] = 0.00; spz[7] = 0.00; // Gamma
  spx[8] = 0.00; spy[8] = 0.00; spz[8] = 0.00; // Gamma
  spx[9] = 0.50; spy[9] = 0.00; spz[9] = 0.50; // N
  bool count = true;
  if (size > 0) count = false;
  size=0;
  for (int i=0; i<10; i+=2) {
    int knumr = knum*len[i/2];
    for (int k=0; k<knumr; k++) {
      if (!count) {
	double x = (double)k/(double)knumr;
	kx[size] = 2*M_PI/a* (spx[i] + (spx[i+1]-spx[i])*x );
	ky[size] = 2*M_PI/a* (spy[i] + (spy[i+1]-spy[i])*x );
	kz[size] = 2*M_PI/a* (spz[i] + (spz[i+1]-spz[i])*x );
	if (k==0) {
	  printf("%4d - %4d: %s - %s\n",size,size+knumr,sp[i],sp[i+1]);
	}
	//printf("%d %d %f %f %f\n",size,i,kx[size]*a/2/M_PI,ky[size]*a/2/M_PI,kz[size]*a/2/M_PI);
      }
      size++;
    }
  }
  if (!count) {
    kx[size]=spx[9]; ky[size]=spy[9]; kz[size]=spz[9];
  }
  size++;
}

void kgen_sc(double *kx, double *ky, double *kz, int knum, int &size)
{
  double len[4] = {1,1,1,1}; // ratio of k-point density
  double spx[8],spy[8],spz[8];
  double a = cell.hmat[0][0]/(double)(phonon_rep*2-1)*2;
  const char *sp[] = {"Gamma", "X", "X", "M", "M", "R", "R", "Gamma"};
  printf("lattice constant a = %e\n",a);
  spx[0] = 0.00; spy[0] = 0.00; spz[0] = 0.00; // Gamma
  spx[1] = 1.00; spy[1] = 0.00; spz[1] = 0.00; // X
  spx[2] = 1.00; spy[2] = 0.00; spz[2] = 0.00; // X
  spx[3] = 1.00; spy[3] = 1.00; spz[3] = 0.00; // M
  spx[4] = 1.00; spy[4] = 1.00; spz[4] = 0.00; // M
  spx[5] = 1.00; spy[5] = 1.00; spz[5] = 1.00; // R
  spx[6] = 1.00; spy[6] = 1.00; spz[6] = 1.00; // R
  spx[7] = 0.00; spy[7] = 0.00; spz[7] = 0.00; // Gamma
  bool count = true;
  if (size > 0) count = false;
  size=0;
  for (int i=0; i<10; i+=2) {
    int knumr = knum*len[i/2];
    for (int k=0; k<knumr; k++) {
      if (!count) {
	double x = (double)k/(double)knumr;
	kx[size] = 2*M_PI/a* (spx[i] + (spx[i+1]-spx[i])*x );
	ky[size] = 2*M_PI/a* (spy[i] + (spy[i+1]-spy[i])*x );
	kz[size] = 2*M_PI/a* (spz[i] + (spz[i+1]-spz[i])*x );
	if (k==0) {
	  printf("%4d - %4d: %s - %s\n",size,size+knumr,sp[i],sp[i+1]);
	}
	//printf("%d %d %f %f %f\n",size,i,kx[size]*a/2/M_PI,ky[size]*a/2/M_PI,kz[size]*a/2/M_PI);
      }
      size++;
    }
  }
  if (!count) {
    kx[size]=spx[7]; ky[size]=spy[7]; kz[size]=spz[7];
  }
  size++;
}

void kgen_hcp(double *kx, double *ky, double *kz, int knum, int &size)
{
  double len[7] = {1,0.5,1,0.5,1,0.5,1}; // ratio of k-point density
  double spx[14],spy[14],spz[14];
  double a = cell.hmat[0][0]/(double)(phonon_rep*2-1);
  double c = cell.hmat[2][2]/(double)(phonon_rep*2-1);
  const char *sp[] = {"Gamma", "M", "M", "K", "K", "Gamma", "Gamma", "A", "A", "L", "L", "H", "H", "A"};
  printf("lattice constant a = %e\n",a);
  spx[ 0] = 0.00; spy[ 0] = 0.00; spz[ 0] = 0.00; // Gamma
  spx[ 1] = 1.00; spy[ 1] = 0.00; spz[ 1] = 0.00; // M
  spx[ 2] = 1.00; spy[ 2] = 0.00; spz[ 2] = 0.00; // M
  spx[ 3] = 1.00; spy[ 3] = 0.50; spz[ 3] = 0.00; // K
  spx[ 4] = 1.00; spy[ 4] = 0.50; spz[ 4] = 0.00; // K
  spx[ 5] = 0.00; spy[ 5] = 0.00; spz[ 5] = 0.00; // Gamma
  spx[ 6] = 0.00; spy[ 6] = 0.00; spz[ 6] = 0.00; // Gamma
  spx[ 7] = 0.00; spy[ 7] = 0.00; spz[ 7] = 1.00; // A
  spx[ 8] = 0.00; spy[ 8] = 0.00; spz[ 8] = 1.00; // A
  spx[ 9] = 1.00; spy[ 9] = 0.00; spz[ 9] = 1.00; // L
  spx[10] = 1.00; spy[10] = 0.00; spz[10] = 1.00; // L
  spx[11] = 1.00; spy[11] = 0.50; spz[11] = 1.00; // H
  spx[12] = 1.00; spy[12] = 0.50; spz[12] = 1.00; // H
  spx[13] = 0.00; spy[13] = 0.00; spz[13] = 1.00; // A
  bool count = true;
  if (size > 0) count = false;
  size=0;
  for (int i=0; i<14; i+=2) {
    int knumr = knum*len[i/2];
    for (int k=0; k<knumr; k++) {
      if (!count) {
	double x = (double)k/(double)knumr;
	kx[size] = 2*M_PI/a* (spx[i] + (spx[i+1]-spx[i])*x );
	ky[size] = 2*M_PI/a* (spy[i] + (spy[i+1]-spy[i])*x );
	kz[size] = 2*M_PI/c* (spz[i] + (spz[i+1]-spz[i])*x );
	if (k==0) {
	  printf("%4d - %4d: %s - %s\n",size,size+knumr,sp[i],sp[i+1]);
	}
	//printf("%d %d %f %f %f\n",size,i,kx[size]*a/2/M_PI,ky[size]*a/2/M_PI,kz[size]*a/2/M_PI);
      }
      size++;
    }
  }
  if (!count) {
    kx[size]=spx[13]; ky[size]=spy[13]; kz[size]=spz[13];
  }
  size++;
}

void kgen_chain(double *kx, double *ky, double *kz, int knum, int &size)
{
  double len[1] = {1}; // ratio of k-point density
  double spx[2],spy[2],spz[2];
  double a;
  a = cell.hmat[0][0]/(double)(phonon_rep*2-1);
  printf("lattice constant a = %e\n",a);
  spx[0] = 0.00; spy[0] = 0.00; spz[0] = 0.00; // Gamma
  spx[1] = 1.00; spy[1] = 0.00; spz[1] = 0.00; // X
  bool count = true;
  if (size > 0) count = false;
  size=0;
  for (int i=0; i<2; i+=2) {
    int knumr = knum*len[i/2];
    for (int k=0; k<knumr; k++) {
      if (!count) {
	double x = (double)k/(double)knumr;
	kx[size] = M_PI/a* (spx[i] + (spx[i+1]-spx[i])*x );
	ky[size] = M_PI/a* (spy[i] + (spy[i+1]-spy[i])*x );
	kz[size] = M_PI/a* (spz[i] + (spz[i+1]-spz[i])*x );
	//printf("%d %d %f %f %f\n",size,i,kx[size]*a/2/M_PI,ky[size]*a/2/M_PI,kz[size]*a/2/M_PI);
      }
      size++;
    }
  }
  if (!count) {
    kx[size]=spx[1]; ky[size]=spy[1]; kz[size]=spz[1];
  }
  size++;
}

void kgen_fcc_dense(double *kx, double *ky, double *kz, double *fkp, int knum, int &size)
{ // LINEAR DENSE MESH IN B.Z. OF FCC LATTICE
  bool count = true;
  if (size > 0) count = false;
  size=0;

  double a;
  if (cell.hmat[1][0]/cell.hmat[0][0] > 0.5) {
    a = cell.hmat[0][0]/(double)(phonon_rep*2-1)*2;
  } else {
    a = cell.hmat[0][0]/(double)(phonon_rep*2-1);
  }
  double ppp = 2.0 * M_PI / a;
  int mesh = dos_kmesh;
  for (int ix=0; ix<mesh; ix++) {
    for (int iy=0; iy<mesh; iy++) {
      for (int iz=0; iz<mesh; iz++) {
	if (!count) {
	  kx[size] = ppp / (double)mesh * (double)ix;
	  ky[size] = ppp / (double)mesh * (double)iy;
	  kz[size] = ppp / (double)mesh * (double)iz;
	  fkp[size] = 1.0;
	}
	size++;
      }
    }
  }

  if (!count) {
    double at = 0.0;
    for (int i=0; i<size; i++) { at += fkp[i]; }
    for (int i=0; i<size; i++) { fkp[i] /= at; }
  }

}

void kgen_fcc_dense89(double *kx, double *ky, double *kz, double *fkp, int knum, int &size)
{ // LINEAR DENSE MESH IN B.Z. OF FCC LATTICE (89 points)
  if (size != 89) {
    size = 89;
    return;
  }
  double a;
  if (cell.hmat[1][0]/cell.hmat[0][0] > 0.5) {
    a = cell.hmat[0][0]/(double)(phonon_rep*2-1)*2;
  } else {
    a = cell.hmat[0][0]/(double)(phonon_rep*2-1);
  }
  double ppp = 2.0 * M_PI / a; double p8 = ppp/8.0;
  for (int i= 1; i<=39; i++) { kz[i-1] =         0.0; }
  for (int i=40; i<=66; i++) { kz[i-1] =     ppp/8.0; }
  for (int i=67; i<=82; i++) { kz[i-1] =     ppp/4.0; }
  for (int i=83; i<=88; i++) { kz[i-1] = 3.0*ppp/8.0; }
  for (int i=89; i<=89; i++) { kz[i-1] =     ppp/2.0; }
  for (int i= 1; i<= 9; i++) { ky[i-1] =         0.0; }
  for (int i=10; i<=17; i++) { ky[i-1] =     ppp/8.0; }
  for (int i=18; i<=24; i++) { ky[i-1] =     ppp/4.0; }
  for (int i=25; i<=30; i++) { ky[i-1] = 3.0*ppp/8.0; }
  for (int i=31; i<=35; i++) { ky[i-1] =     ppp/2.0; }
  for (int i=36; i<=38; i++) { ky[i-1] = 5.0*ppp/8.0; }
  for (int i=39; i<=39; i++) { ky[i-1] = 3.0*ppp/4.0; }
  for (int i=40; i<=47; i++) { ky[i-1] =     ppp/8.0; }
  for (int i=48; i<=54; i++) { ky[i-1] =     ppp/4.0; }
  for (int i=55; i<=60; i++) { ky[i-1] = 3.0*ppp/8.0; }
  for (int i=61; i<=64; i++) { ky[i-1] =     ppp/2.0; }
  for (int i=65; i<=66; i++) { ky[i-1] = 5.0*ppp/8.0; }
  for (int i=67; i<=73; i++) { ky[i-1] =     ppp/4.0; }
  for (int i=74; i<=78; i++) { ky[i-1] = 3.0*ppp/8.0; }
  for (int i=79; i<=81; i++) { ky[i-1] =     ppp/2.0; }
  for (int i=82; i<=82; i++) { ky[i-1] = 5.0*ppp/8.0; }
  for (int i=83; i<=86; i++) { ky[i-1] = 3.0*ppp/8.0; }
  for (int i=87; i<=89; i++) { ky[i-1] =     ppp/2.0; }
  for (int i= 1; i<= 9; i++) { kx[i-1] = p8*(double)(i- 1); }
  for (int i=10; i<=17; i++) { kx[i-1] = p8*(double)(i- 9); }
  for (int i=18; i<=24; i++) { kx[i-1] = p8*(double)(i-16); }
  for (int i=25; i<=30; i++) { kx[i-1] = p8*(double)(i-22); }
  for (int i=31; i<=35; i++) { kx[i-1] = p8*(double)(i-27); }
  for (int i=36; i<=38; i++) { kx[i-1] = p8*(double)(i-31); }
  for (int i=39; i<=39; i++) { kx[i-1] = p8*6.0;            }
  for (int i=40; i<=47; i++) { kx[i-1] = p8*(double)(i-39); }
  for (int i=48; i<=54; i++) { kx[i-1] = p8*(double)(i-46); }
  for (int i=55; i<=60; i++) { kx[i-1] = p8*(double)(i-52); }
  for (int i=61; i<=64; i++) { kx[i-1] = p8*(double)(i-57); }
  for (int i=65; i<=66; i++) { kx[i-1] = p8*(double)(i-60); }
  for (int i=67; i<=73; i++) { kx[i-1] = p8*(double)(i-65); }
  for (int i=74; i<=78; i++) { kx[i-1] = p8*(double)(i-71); }
  for (int i=79; i<=81; i++) { kx[i-1] = p8*(double)(i-75); }
  for (int i=82; i<=82; i++) { kx[i-1] = p8*5.0;            }
  for (int i=83; i<=86; i++) { kx[i-1] = p8*(double)(i-80); }
  for (int i=87; i<=88; i++) { kx[i-1] = p8*(double)(i-83); }
  for (int i=89; i<=89; i++) { kx[i-1] = p8*4.0;            }

  fkp[ 1-1] = 1.0/48.0; // Gamma
  for (int i= 2; i<= 8; i++) { fkp[i-1] = 1.0/8.0; } // Gamma - X
  fkp[ 9-1] = 1.0/16.0; // X
  fkp[10-1] = 1.0/ 4.0;
  for (int i=11; i<=16; i++) { fkp[i-1] = 1.0/2.0; }
  fkp[17-1] = 1.0/ 4.0;
  fkp[18-1] = 1.0/ 4.0;
  for (int i=19; i<=23; i++) { fkp[i-1] = 1.0/2.0; }
  for (int i=24; i<=25; i++) { fkp[i-1] = 1.0/4.0; }
  for (int i=26; i<=29; i++) { fkp[i-1] = 1.0/2.0; }
  for (int i=30; i<=31; i++) { fkp[i-1] = 1.0/4.0; }
  for (int i=32; i<=34; i++) { fkp[i-1] = 1.0/2.0; }
  fkp[35-1] = 1.0/ 8.0;
  fkp[36-1] = 1.0/ 4.0;
  fkp[37-1] = 1.0/ 2.0;
  fkp[38-1] = 1.0/ 6.0;
  fkp[39-1] = 1.0/12.0;
  fkp[40-1] = 1.0/ 6.0;
  for (int i=41; i<=46; i++) { fkp[i-1] = 1.0/2.0; }
  fkp[47-1] = 1.0/ 4.0;
  fkp[48-1] = 1.0/ 2.0;
  for (int i=49; i<=53; i++) { fkp[i-1] = 1.0    ; }
  for (int i=54; i<=55; i++) { fkp[i-1] = 1.0/2.0; }
  for (int i=56; i<=59; i++) { fkp[i-1] = 1.0    ; }
  fkp[60-1] = 1.0/ 3.0;
  fkp[61-1] = 1.0/ 2.0;
  for (int i=62; i<=63; i++) { fkp[i-1] = 1.0    ; }
  for (int i=64; i<=66; i++) { fkp[i-1] = 1.0/2.0; }
  fkp[67-1] = 1.0/ 6.0;
  for (int i=68; i<=88; i++) { fkp[i-1] = 1.0/2.0; }
  fkp[73-1] = 1.0/ 6.0;
  for (int i=75; i<=77; i++) { fkp[i-1] = 1.0;     }
  fkp[80-1] = 1.0;
  fkp[82-1] = 1.0/ 4.0;
  fkp[83-1] = 1.0/ 6.0;
  fkp[86-1] = 1.0/ 4.0;
  fkp[89-1] = 1.0/12.0;

  double at = 0.0;
  for (int i=1; i<=89; i++) { at += fkp[i-1]; }
  for (int i=1; i<=89; i++) { fkp[i-1] /= at; }

}

void kgen_bcc_dense(double *kx, double *ky, double *kz, int knum, int &size)
{
  /*  
      LINEAR DENSE MESH IN B.Z. OF BCC LATTICE (55 points)
  */

  /*
      PPP=2.0D0*PI/A
      DO I=1,25                                                      
      bzk(I,3)=0.0D0                                                       
      enddo
      DO I=26,41                                                     
      bzk(I,3)=PPP/8.0D0                                                     
      enddo
      DO I=42,50                                                     
      bzk(I,3)=PPP/4.0D0                                                     
      enddo
      DO I=51,54                                                     
      bzk(I,3)=3.0D0*PPP/8.0D0                                               
      enddo
      bzk(55,3)=PPP/2.0D0                                                    
C
      DO I=1,9                                                       
      bzk(I,2)=0.0D0                                                       
      enddo
      DO I=10,16                                                     
      bzk(I,2)=PPP/8.0D0                                                     
      enddo
      DO I=17,21                                                     
      bzk(I,2)=PPP/4.0D0                                                     
      enddo
      DO I=22,24                                                     
      bzk(I,2)=PPP*3.0D0/8.0D0                                               
      enddo
      DO I=25,25                                                     
      bzk(I,2)=PPP/2.0D0                                                     
      enddo
C
      DO I=26,32                                                     
      bzk(I,2)=PPP/8.0D0                                                     
      enddo
      DO I=33,37                                                     
      bzk(I,2)=PPP/4.0D0                                                     
      enddo
      DO I=38,40                                                     
      bzk(I,2)=PPP*3.0D0/8.0D0                                               
      enddo
      DO I=41,41                                                     
      bzk(I,2)=PPP/2.0D0                                                     
      enddo
C
      DO I=42,46                                                     
      bzk(I,2)=PPP/4.0D0                                                     
      enddo
      DO I=47,49                                                     
      bzk(I,2)=PPP*3.0D0/8.0D0                                               
      enddo
      DO I=50,50                                                     
      bzk(I,2)=PPP/2.0D0                                                     
      enddo
C
      DO I=51,53                                                     
      bzk(I,2)=PPP*3.0D0/8.0D0                                               
      enddo
      DO I=54,55                                                     
      bzk(I,2)=PPP/2.0D0                                                     
      enddo
C
      P8=PPP/8.0D0                                                        
      DO I=1,9                                                      
      bzk(I,1)=P8*(I-1)                                                    
      enddo
      DO I=10,16                                                    
      bzk(I,1)=P8*(I-9)                                                    
      enddo
      DO I=17,21                                                    
      bzk(I,1)=P8*(I-15)                                                   
      enddo
      DO I=22,24                                                    
      bzk(I,1)=P8*(I-19)                                                   
      enddo
      DO I=25,25                                                    
      bzk(I,1)=P8*(I-21)                                                   
      enddo
C
      DO I=26,41                                                    
      bzk(I,1)=bzk(I-16,1)                                                    
      enddo
C
      DO I=42,50                                                    
      bzk(I,1)=bzk(I-25,1)                                                    
      enddo
      DO I=51,54                                                    
      bzk(I,1)=bzk(I-29,1)                                                    
      enddo
      DO I=55,55                                                    
      bzk(I,1)=bzk(25,1)                                                      
      enddo
C                                                                       
C                                                                       
      fkp(1)=1.0D0/48.0D0                                                 
      DO I=2,8                                                      
      fkp(I)=1.0D0/8.0D0                                                  
      enddo
      fkp(9)=1.0D0/48.0D0                                                 
      fkp(10)=1.0D0/4.0D0                                                 
      fkp(16)=fkp(10)                                                       
      DO I=11,15                                                    
      fkp(I)=1.0D0/2.0D0                                                  
      enddo
      fkp(17)=fkp(10)                                                       
      DO I=18,20                                                    
      fkp(I)=1.0D0/2.0D0                                                  
      enddo
      fkp(21)=1.0D0/4.0D0                                                 
      fkp(22)=fkp(21)                                                       
      fkp(23)=fkp(20)                                                       
      fkp(24)=fkp(22)                                                       
      fkp(25)=1.0D0/8.0D0                                                 
C
      fkp(26)=1.0D0/6.0D0                                                 
      DO I=27,31                                                    
      fkp(I)=1.0D0/2.0D0                                                  
      enddo
      fkp(32)=1.0D0/6.0D0                                                 
      fkp(33)=1.0D0/2.0D0                                                 
      DO I=34,36                                                    
      fkp(I)=1.0D0                                                        
      enddo
      fkp(37)=1.0D0/2.0D0                                                 
      fkp(38)=1.0D0/2.0D0                                                 
      fkp(40)=1.0D0/2.0D0                                                 
      fkp(41)=1.0D0/4.0D0                                                 
      fkp(39)=1.0D0                                                       
C
      fkp(48)=1.0D0                                                       
      fkp(42)=1.0D0/6.0D0                                                 
      fkp(46)=fkp(42)                                                       
      fkp(47)=1.0D0/2.0D0                                                 
      fkp(49)=fkp(47)                                                       
      fkp(50)=1.0D0/4.0D0                                                 
      DO I=43,45                                                    
      fkp(I)=1.0D0/2.0D0                                                  
      enddo
C
      fkp(51)=1.0D0/6.0D0                                                 
      fkp(52)=1.0D0/2.0D0                                                 
      fkp(53)=1.0D0/6.0D0                                                 
      fkp(54)=1.0D0/4.0D0                                                 
C
      fkp(55)=1.0D0/24.0D0                                                
C                                                                       
C                                                                       
      AT=0.0D0                                                          
      DO I=1,55                                                     
      AT=AT+fkp(I)                                                        
      enddo
      DO I=1,55                                                     
      fkp(I)=fkp(I)/AT                                                     
      enddo
c      WRITE(6,610) AT                                                   
c  610 FORMAT(1H /1H ,'BCCKP  AT=',D24.16)                               
c      RETURN                                                            
c      END                                                               
  */
}
