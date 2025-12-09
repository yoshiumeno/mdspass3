//
// Brenner potential
// Modified for stress calculation and for a non-orthogonal cell (Nov. 2012)
// Stress calculation is tested for cnt and diamond
//
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "myheader.h"

#define NCC 10000
#if defined __linux__ || defined __APPLE__
#include <unistd.h>
#else
#include <direct.h>
#endif

void get_first_arg(std::string &line, std::string &arg1);
void pibond(); void airebo_lj(); void airebo_lj_read();
void airebo_lj(int mode);
void radic(int ki, int kj, double xnt1, double xnt2, double conjug,
	   double &rad, double &drdl, double &drdm, double &drdn);
void tor(double xnt1, double xnt2, double conjug, double &ator,
	 double &drdl, double &drdm, double &drdn);
void bcuint(int kl, int ki, double xx1, double xx2, int nh, int nc,
	    double &ansy, double &ansy1, double &ansy2);
double anint(double x);
void bond_set();
void e_force_brenner(int mode);
extern char cwdname[80];
extern int relax_accel;

void brealloc()
{
  if (bre.npmax<atom.natom) {
    //std::cout<<"Error: npmax (in Bre) is smaller than natom"<<std::endl; exit(0); }
    printf("\n");
    printf("  Force-Brenner: npmax (in Bre, %d) is smaller than natom (=%d)\n",bre.npmax,atom.natom);
    printf("  so it is now increased to natom.\n");
    bre.npmax = atom.natom;
    bre.npm1 = bre.npmax + 1;
  }
  //1st order
  if (bre.xmass) { delete[] bre.xmass; bre.xmass = NULL; }
  if (bre.noa) { delete[] bre.noa; bre.noa = NULL; }
  if (bre.bww) { delete[] bre.bww; bre.bww = NULL; }
  if (bre.dww) { delete[] bre.dww; bre.dww = NULL; }
  if (bre.rcor) { delete[] bre.rcor; bre.rcor = NULL; }
  if (bre.cor1) { delete[] bre.cor1; bre.cor1 = NULL; }
  if (bre.cor2) { delete[] bre.cor2; bre.cor2 = NULL; }
  if (bre.cor3) { delete[] bre.cor3; bre.cor3 = NULL; }
  if (bre.rep1) { delete[] bre.rep1; bre.rep1 = NULL; }
  if (bre.rep2) { delete[] bre.rep2; bre.rep2 = NULL; }
  if (bre.rep3) { delete[] bre.rep3; bre.rep3 = NULL; }
  if (bre.list) { delete[] bre.list; bre.list = NULL; }
  if (bre.nabors) { delete[] bre.nabors; bre.nabors = NULL; }
  if (bre.lcheck) { delete[] bre.lcheck; bre.lcheck = NULL; }
  if (bre.exx1) { delete[] bre.exx1; bre.exx1 = NULL; }
  if (bre.dexx1) { delete[] bre.dexx1; bre.dexx1 = NULL; }
  if (bre.exx2) { delete[] bre.exx2; bre.exx2 = NULL; }
  if (bre.ivct2b) { delete[] bre.ivct2b; bre.ivct2b = NULL; }
  if (bre.jvct2b) { delete[] bre.jvct2b; bre.jvct2b = NULL; }
  if (bre.xhc1) { delete[] bre.xhc1; bre.xhc1 = NULL; }
  if (bre.xhc2) { delete[] bre.xhc2; bre.xhc2 = NULL; }
  if (bre.ktype) { delete[] bre.ktype; bre.ktype = NULL; }
  if (bre.r01) { delete[] bre.r01; bre.r01 = NULL; }
  if (bre.r02) { delete[] bre.r02; bre.r02 = NULL; }
  if (bre.r03) { delete[] bre.r03; bre.r03 = NULL; }
  if (bre.rnp1) { delete[] bre.rnp1; bre.rnp1 = NULL; }
  if (bre.rnp2) { delete[] bre.rnp2; bre.rnp2 = NULL; }
  if (bre.rnp3) { delete[] bre.rnp3; bre.rnp3 = NULL; }
  if (bre.rpp1) { delete[] bre.rpp1; bre.rpp1 = NULL; }
  if (bre.rpp2) { delete[] bre.rpp2; bre.rpp2 = NULL; }
  if (bre.rpp3) { delete[] bre.rpp3; bre.rpp3 = NULL; }
  if (bre.sig) {
    for (int i=0; i<bre.ntypes+1; i++) { delete[] bre.sig[i]; }
    delete[] bre.sig; bre.sig = NULL; }
  if (bre.eps) {
    for (int i=0; i<bre.ntypes+1; i++) { delete[] bre.eps[i]; }
    delete[] bre.eps; bre.eps = NULL; }
  if (bre.atable) {
    for (int i=0; i<5; i++) {
      for (int j=0; j<5; j++) { delete[] bre.atable[i][j]; }
      delete[] bre.atable[i]; } delete[] bre.atable; bre.atable = NULL; }
  if (bre.datable) {
    for (int i=0; i<5; i++) {
      for (int j=0; j<5; j++) { delete[] bre.datable[i][j]; }
      delete[] bre.datable[i]; } delete[] bre.datable; bre.datable = NULL; }
  if (bre.rtable) {
    for (int i=0; i<5; i++) {
      for (int j=0; j<5; j++) { delete[] bre.rtable[i][j]; }
      delete[] bre.rtable[i]; } delete[] bre.rtable; bre.rtable = NULL; }
  if (bre.drtable) {
    for (int i=0; i<5; i++) {
      for (int j=0; j<5; j++) { delete[] bre.drtable[i][j]; }
      delete[] bre.drtable[i]; } delete[] bre.drtable; bre.drtable = NULL; }
  if (bre.tabfc) {
    for (int i=0; i<5; i++) {
      for (int j=0; j<5; j++) { delete[] bre.tabfc[i][j]; }
      delete[] bre.tabfc[i]; } delete[] bre.tabfc; bre.tabfc = NULL; }
  if (bre.tabdfc) {
    for (int i=0; i<5; i++) {
      for (int j=0; j<5; j++) { delete[] bre.tabdfc[i][j]; }
      delete[] bre.tabdfc[i]; } delete[] bre.tabdfc; bre.tabdfc = NULL; }

  bre.xmass = new double[bre.ntypes+1];
  bre.noa = new int[bre.ntypes+1];
  bre.bww = new double[bre.nlmax+1];
  bre.dww = new double[bre.nlmax+1];
  bre.rcor = new double[bre.nlmax+1];
  bre.cor1 = new double[bre.nlmax+1];
  bre.cor2 = new double[bre.nlmax+1];
  bre.cor3 = new double[bre.nlmax+1];
  bre.rep1 = new int[bre.nlmax+1];
  bre.rep2 = new int[bre.nlmax+1];
  bre.rep3 = new int[bre.nlmax+1];
  bre.list = new int[bre.nlmax+1];
  bre.nabors = new int[bre.npm1+1];
  bre.lcheck = new int[bre.nlmax+1];
  bre.exx1 = new double[bre.nlmax+1];
  bre.dexx1 = new double[bre.nlmax+1];
  bre.exx2 = new double[bre.nlmax+1];
  bre.ivct2b = new int[bre.nlmax+1];
  bre.jvct2b = new int[bre.nlmax+1];
  bre.xhc1 = new double[bre.npmax+1];
  bre.xhc2 = new double[bre.npmax+1];
  bre.ktype = new int[bre.npmax+1];
  bre.r01 = new double[bre.npmax+1];
  bre.r02 = new double[bre.npmax+1];
  bre.r03 = new double[bre.npmax+1];
  bre.rnp1 = new double[bre.npmax+1];
  bre.rnp2 = new double[bre.npmax+1];
  bre.rnp3 = new double[bre.npmax+1];
  bre.rpp1 = new double[bre.nlmax+1];
  bre.rpp2 = new double[bre.nlmax+1];
  bre.rpp3 = new double[bre.nlmax+1];
  //2nd order
  bre.sig = new double*[bre.ntypes+1];
  bre.eps = new double*[bre.ntypes+1];
  for (int i=0; i<bre.ntypes+1; i++) {
    bre.sig[i] = new double[bre.ntypes+1];
    bre.eps[i] = new double[bre.ntypes+1]; }
  //3rd order
  bre.atable = new double**[5];
  bre.datable = new double**[5];
  bre.rtable = new double**[5];
  bre.drtable = new double**[5];
  bre.tabfc = new double**[5];
  bre.tabdfc = new double**[5];
  for (int i=0; i<5; i++) {
    bre.atable[i] = new double*[5];
    bre.datable[i] = new double*[5];
    bre.rtable[i] = new double*[5];
    bre.drtable[i] = new double*[5];
    bre.tabfc[i] = new double*[5];
    bre.tabdfc[i] = new double*[5];
    for (int j=0; j<5; j++) {
      bre.atable[i][j] = new double[bre.ntab+1];
      bre.datable[i][j] = new double[bre.ntab+1];
      bre.rtable[i][j] = new double[bre.ntab+1];
      bre.drtable[i][j] = new double[bre.ntab+1];
      bre.tabfc[i][j] = new double[bre.ntab+1];
      bre.tabdfc[i][j] = new double[bre.ntab+1];
    } }
}

void bresetin(int mode)
{
  for (int i=1; i<=bre.ntypes; i++) {
    bre.xmass[i] = 0.0;
    bre.noa[i] = 0;
    for (int j=1; j<=bre.ntypes; j++) {
      bre.sig[i][j] = 0.0; bre.eps[i][j] = 0.0;
    } }
  bre.xmass[1] = 12.0; bre.xmass[2] = 1.0; bre.xmass[3] = 28.0; bre.xmass[4] = 72.0;
  bre.kt[1] =2; bre.kt[6] =1; bre.kt[14] =3; bre.kt[32] =4;
  bre.kt2[1]=6; bre.kt2[2]=1; bre.kt2[3]=14; bre.kt2[4]=32;
  for (int i= 1; i<=16; i++) { bre.igc[i]=4; }
  for (int i=17; i<=18; i++) { bre.igc[i]=3; }
  for (int i=19; i<=20; i++) { bre.igc[i]=2; }
  for (int i=21; i<=25; i++) { bre.igc[i]=1; }
  for (int i= 1; i<=18; i++) { bre.igh[i]=3; }
  for (int i=19; i<=22; i++) { bre.igh[i]=2; }
  for (int i=23; i<=25; i++) { bre.igh[i]=1; }
  for (int i=1; i<=2; i++) { 
    for (int j=1; j<=10; j++) { 
      for (int k=1; k<=10; k++) {
	bre.xh[i][j][k]=0.0; bre.xh1[i][j][k]=0.0; bre.xh2[i][j][k]=0.0; } } }
  bre.att = 3.20;
  bre.xqm = 3.70;
  bre.spgc[1][1] = 0.2817216000000e+00;
  bre.spgc[2][1] = 0.1062912000000e+01;
  bre.spgc[3][1] = 0.2136736000000e+01;
  bre.spgc[4][1] = 0.2533952000000e+01;
  bre.spgc[5][1] = 0.1554736000000e+01;
  bre.spgc[6][1] = 0.3863296000000e+00;
  bre.spgc[1][2] = 0.2817216000000e+00;
  bre.spgc[2][2] = 0.1062912000000e+01;
  bre.spgc[3][2] = 0.2136736000000e+01;
  bre.spgc[4][2] = 0.2533952000000e+01;
  bre.spgc[5][2] = 0.1554736000000e+01;
  bre.spgc[6][2] = 0.3863296000000e+00;
  bre.spgc[1][3] = 0.6900668660000e+00;
  bre.spgc[2][3] = 0.5460691360000e+01;
  bre.spgc[3][3] = 0.2301345680000e+02;
  bre.spgc[4][3] = 0.5491519344000e+02;
  bre.spgc[5][3] = 0.6862037040000e+02;
  bre.spgc[6][3] = 0.3470897779200e+02;
  if (mode==0){
  bre.spgc[1][4] = 0.3754490870000e+00;
  bre.spgc[2][4] = 0.1407252749388e+01;
  bre.spgc[3][4] = 0.2255103926323e+01;
  bre.spgc[4][4] = 0.2028902219952e+01;
  bre.spgc[5][4] = 0.1426981217906e+01;
  bre.spgc[6][4] = 0.5063107994308e+00;
  bre.spgc[1][5] = 0.2718558000000e+00;
  bre.spgc[2][5] = 0.4892727456293e+00;
  bre.spgc[3][5] =-0.4328199017473e+00;
  bre.spgc[4][5] =-0.5616795197048e+00;
  bre.spgc[5][5] = 0.1270874966906e+01;
  bre.spgc[6][5] =-0.3750409108350e-01; //brenner
  } else {
  bre.spgc[1][5] = 0.3754490870000e+00;
  bre.spgc[2][5] = 0.1407252749388e+01;
  bre.spgc[3][5] = 0.2255103926323e+01;
  bre.spgc[4][5] = 0.2028902219952e+01;
  bre.spgc[5][5] = 0.1426981217906e+01;
  bre.spgc[6][5] = 0.5063107994308e+00;
  bre.spgc[1][4] = 0.2718558000000e+00;
  bre.spgc[2][4] = 0.4892727456293e+00;
  bre.spgc[3][4] =-0.4328199017473e+00;
  bre.spgc[4][4] =-0.5616795197048e+00;
  bre.spgc[5][4] = 0.1270874966906e+01;
  bre.spgc[6][4] =-0.3750409108350e-01; //stuart
  }
  bre.spgh[1][1] =  270.467795364007301;
  bre.spgh[2][1] = 1549.701314596994564;
  bre.spgh[3][1] = 3781.927258631323866;
  bre.spgh[4][1] = 4582.337619544424228;
  bre.spgh[5][1] = 2721.538161662818368;
  bre.spgh[6][1] =  630.658598136730774;
  bre.spgh[1][2] =   16.956325544514659;
  bre.spgh[2][2] =  -21.059084522755980;
  bre.spgh[3][2] = -102.394184748124742;
  bre.spgh[4][2] = -210.527926707779059;
  bre.spgh[5][2] = -229.759473570467513;
  bre.spgh[6][2] =  -94.968528666251945;
  bre.spgh[1][3] =   19.065031149937783;
  bre.spgh[2][3] =    2.017732531534021;
  bre.spgh[3][3] =   -2.566444502991983;
  bre.spgh[4][3] =    3.291353893907436;
  bre.spgh[5][3] =   -2.653536801884563;
  bre.spgh[6][3] =    0.837650930130006;
}
void breparam(int mode)
{
  double pi = 4.0*atan(1.0);
  //FILE *fp14 = fopen("pot/inter3d_iv_new.d","r");
  //FILE *fp15 = fopen("pot/inter2d_iv.d","r");
  //FILE *fp16 = fopen("pot/inter3dtors.d","r");
  //FILE *fp17 = fopen("pot/inter3d_h.d","r");
  //FILE *fp18 = fopen("pot/inter3d_ch.d","r");
#if defined __linux__ || defined __APPLE__
  chdir(cwdname);
  std::ifstream f14("pot/inter3d_iv_new.d");
  std::ifstream f15("pot/inter2d_iv.d");
  printf("  MODE (0:Brenner, 1:AIREBO) = %d\n",mode);
  if (mode==1) {
  f15.close();
  f15.open("pot/inter2d_iv_stuart.d");
  //std::ifstream f15("pot/inter2d_iv_stuart.d");
  }
  std::ifstream f16("pot/inter3dtors.d");
  std::ifstream f17("pot/inter3d_h.d");
  std::ifstream f18("pot/inter3d_ch.d");
#else
  std::ifstream f14("pot\\inter3d_iv_new.d");
  std::ifstream f15("pot\\inter2d_iv.d");
  if (mode==1) {
  f15.close();
  f15.open("pot\\inter2d_iv_stuart.d");
  //std::ifstream f15("pot\\inter2d_iv_stuart.d");
  }
  std::ifstream f16("pot\\inter3dtors.d");
  std::ifstream f17("pot\\inter3d_h.d");
  std::ifstream f18("pot\\inter3d_ch.d");
#endif
  std::string line, arg1, arg2;
  
  for (int i=1; i<=2; i++) {
    for (int l=1; l<=10; l++) {
      for (int m=1; m<=10; m++) {
	for (int j=1; j<=16; j++) {
	  bre.clm[i][l][m][j] = 0.0; } } } }
  int i2d, i3d, itd, ii,jj,kk,ll,mm,nn;
  double xhh;
  f15>>i2d;
  while(1) {
    f15>>ii; f15>>jj; f15>>kk; f15>>xhh;
    if (ii<=0) { break; }
    bre.xh[ii][jj][kk] = xhh; }
  bre.xh1[2][3][1]=(bre.xh[2][4][1]-bre.xh[2][2][1])/2.0;
  bre.xh1[2][2][2]=(bre.xh[2][3][2]-bre.xh[2][1][2])/2.0;
  bre.xh2[2][2][2]=(bre.xh[2][2][3]-bre.xh[2][2][1])/2.0;
  bre.xh2[2][1][3]=(bre.xh[2][1][4]-bre.xh[2][1][2])/2.0;
  while (!f15.eof()) {
    f15>>ii; f15>>ll; f15>>mm;
    for (int j=1; j<=16; j++) { f15>>bre.clm[ii][ll][mm][j]; }
  }
  int ic = 0;
  for (int i=1; i<=4; i++) {
    for (int j=1; j<=4; j++) {
      ic++;
      bre.in2[ic][1] = i-1;
      bre.in2[ic][2] = j-1; } }
  ic = 0;
  for (int i=1; i<=4; i++) {
    for (int j=1; j<=4; j++) {
      for (int k=1; k<=4; k++) {
	ic++;
	bre.in3[ic][1] = i-1;
	bre.in3[ic][2] = j-1;
	bre.in3[ic][3] = k-1; } } }
  f14>>i3d;
  while (!f14.eof()) {
    f14>>ll; f14>>mm; f14>>nn;
    for (int i=1; i<=64; i++) { f14>>bre.clmn[1][ll][mm][nn][i]; }
  }
  f17>>i3d;
  while (!f17.eof()) {
    f17>>ll; f17>>mm; f17>>nn;
    for (int i=1; i<=64; i++) { f17>>bre.clmn[3][ll][mm][nn][i]; }
  }
  f18>>i3d;
  while (!f18.eof()) {
    f18>>ll; f18>>mm; f18>>nn;
    for (int i=1; i<=64; i++) { f18>>bre.clmn[2][ll][mm][nn][i]; }
  }
  for (int l=1; l<=10; l++) {
    for (int m=1; m<=10; m++) {
      for (int i=1; i<=64; i++) {
	bre.clmn[1][l][m][10][i] = bre.clmn[1][l][m][9][i];
	bre.clmn[2][l][m][10][i] = bre.clmn[2][l][m][9][i];
	for (int n=6; n<=10; n++) {
	  bre.clmn[3][l][m][n][i] = bre.clmn[3][l][m][5][i]; } } } }
  f16>>itd;
  while (!f16.eof()) {
    f16>>ll; f16>>mm; f16>>nn;
    for (int i=1; i<=64; i++) { f16>>bre.tlmn[ll][mm][nn][i]; }
  }
  for (int l=1; l<=10; l++) {
    for (int m=1; m<=10; m++) {
      for (int n=4; n<=10; n++) {
	for (int i=1; i<=64; i++) {
	  bre.tlmn[l][m][n][i] = bre.tlmn[l][m][3][i]; } } } }
  if ((itd!=i2d)||(itd!=i3d)) {
    std::cout<<"Brenner error: Incompatible potential types!!!"<<std::endl;
    f14.close(); f15.close(); f16.close(); f17.close(); f18.close();
    f14.clear(); f15.clear(); f16.clear(); f17.clear(); f18.clear();
    exit(0);
  }
  bre.pq = pi/(bre.xqm-bre.att);
  for (int i=1; i<=4; i++) {
    for (int j=1; j<=4; j++) {
      bre.ad[i][j] = 0.0; bre.axl[i][j] = 0.0; bre.bd[i][j] = 0.0;
      bre.bxl[i][j] = 0.0; bre.cd[i][j] = 0.0; bre.cxl[i][j] = 0.0;
      bre.dd[i][j] = 0.0; bre.dxl[i][j] = 0.0; bre.ed[i][j] = 0.0;
      bre.rb1[i][j] = 0.0; bre.rb2[i][j] = 0.0; bre.rmax[i][j] = 0.0;
      bre.pid[i][j] = 0.0; bre.chi[i][j] = 0.0;
      for (int k=1; k<=4; k++) {
	bre.xdb[i][j][k] = 0.0;
      } } }
  bre.ndihed = 2;
  //Carbon
  bre.ad[1][1]=12388.79197798375;
  bre.axl[1][1]=4.720452312717397;
  bre.bd[1][1]=17.56740646508968;
  bre.bxl[1][1]=1.433213249951261;
  bre.cd[1][1]=30.71493208065162;
  bre.cxl[1][1]=1.382691250599169;
  bre.dd[1][1]=10953.54416216992;
  bre.dxl[1][1]=4.746539060659529;
  bre.ed[1][1]=0.3134602960832605;
  bre.rb1[1][1]=1.7;
  bre.rb2[1][1]=2.0;
  bre.rmax[1][1]=bre.rb2[1][1];
  bre.pid[1][1]=pi/(bre.rb2[1][1]-bre.rb1[1][1]);
  bre.epslj[1][1]=0.00284; //AIREBO, epsilon LJ
  bre.siglj[1][1]=3.40; //AIREBO, sigma LJ
  bre.epsts[1][1]=0.3079; //AIREBO, epsilon torsional
  bre.bmin[1][1]=0.77; bre.bmax[1][1]=0.81; //AIREBO, b-switch
  //Hydrogen
  bre.ad[2][2]=29.6325931;
  bre.axl[2][2]=1.715892169856421;
  bre.bd[2][2]=0.0;
  bre.bxl[2][2]=1.0;
  bre.cd[2][2]=0.0;
  bre.cxl[2][2]=1.0;
  bre.dd[2][2]=32.81735574722296;
  bre.dxl[2][2]=3.536298648376465;
  bre.ed[2][2]=0.3704714870452888;
  bre.rb1[2][2]=1.10;
  bre.rb2[2][2]=1.70;
  bre.rmax[2][2]=bre.rb2[2][2];
  bre.pid[2][2]=pi/(bre.rb2[2][2]-bre.rb1[2][2]);
  bre.epslj[2][2]=0.00150; //AIREBO, epsilon LJ
  bre.siglj[2][2]=2.65; //AIREBO, sigma LJ
  bre.epsts[2][2]=0.1250; //AIREBO, epsilon torsional
  bre.bmin[2][2]=0.32; bre.bmax[2][2]=0.42; //AIREBO, b-switch
  //Carbon-Hydrogen
  bre.ad[2][1]=32.35518665873256;
  bre.axl[2][1]=1.434458059249837;
  bre.dd[2][1]=149.9409872288120;
  bre.dxl[2][1]= 4.102549828548784;
  bre.ed[2][1]=0.3407757282257080;
  bre.bd[2][1]=0.0;
  bre.bxl[2][1]=1.0;
  bre.cd[2][1]=0.0;
  bre.cxl[2][1]=1.0;
  bre.ad[1][2]=bre.ad[2][1];
  bre.axl[1][2]=bre.axl[2][1];
  bre.bd[1][2]=bre.bd[2][1];
  bre.bxl[1][2]=bre.bxl[2][1];
  bre.cd[1][2]=bre.cd[2][1];
  bre.cxl[1][2]=bre.cxl[2][1];
  bre.dd[1][2]=bre.dd[2][1];
  bre.dxl[1][2]=bre.dxl[2][1];
  bre.ed[1][2]=bre.ed[2][1];
  bre.rb1[2][1]=1.3;
  bre.rb2[2][1]=1.8;
  bre.rmax[2][1]=bre.rb2[2][1];
  bre.pid[2][1]=pi/(bre.rb2[2][1]-bre.rb1[2][1]);
  bre.pidt=pi/0.30;
  bre.rb1[1][2]=bre.rb1[2][1];
  bre.rb2[1][2]=bre.rb2[2][1];
  bre.rmax[1][2]=bre.rb2[1][2];
  bre.pid[1][2]=pi/(bre.rb2[1][2]-bre.rb1[1][2]);
  bre.epslj[1][2]=0.002063976744; //AIREBO, epsilon LJ
  bre.siglj[1][2]=3.025; //AIREBO, sigma LJ
  bre.epsts[1][2]=0.1787; //AIREBO, epsilon torsional
  bre.epslj[2][1]=bre.epslj[1][2]; //AIREBO
  bre.siglj[2][1]=bre.siglj[1][2]; //AIREBO
  bre.epsts[2][1]=bre.epsts[1][2]; //AIREBO
  bre.bmin[1][2]=0.75; bre.bmax[1][2]=0.90; //AIREBO, b-switch
  bre.bmin[2][1]=0.75; bre.bmax[2][1]=0.90; //AIREBO, b-switch
  for (int i=1; i<=2; i++) {
    for (int j=1; j<=2; j++) {
      for (int k=1; k<=2; k++) {
	bre.xdb[i][j][k]=0.0;
	bre.reg[i][j][k]=1.0; } } }
  bre.xxdb=4.0;
  bre.rhh=0.7415886997;
  bre.rch=1.09;
  bre.xdb[2][2][2]=4.0;
  bre.xdb[2][1][2]=4.0;
  bre.xdb[2][2][1]=4.0;
  bre.xdb[2][1][1]=4.0;
  bre.xdb[1][2][1]=0.0;
  bre.xdb[1][2][2]=0.0;
  bre.reg[2][1][2]=exp(bre.xdb[2][1][2]*(bre.rhh-bre.rch));
  bre.reg[2][2][1]=exp(bre.xdb[2][2][1]*(bre.rch-bre.rhh));
  
  bre.sigma=1.0;
  bre.epsi=11605.0;
  bre.rll=0.0; // ???

  for (int i=1; i<=4; i++) {
    for (int j=1; j<=4; j++) {
      double xx = bre.rmax[i][j]+bre.rll;
      bre.rlist[i][j]=xx*xx;
      bre.rmax[i][j]=bre.rmax[i][j]*bre.rmax[i][j]; } }
  
    f14.close(); f15.close(); f16.close(); f17.close(); f18.close();
    f14.clear(); f15.clear(); f16.clear(); f17.clear(); f18.clear();
}
void bremtable()
{
  double va, dva, vb, dvb, vc, dvc, vv, dvv;
  double fc, dfc, dtemp, ff1, df1, ff2, df2, dvm;
  for (int ki=1; ki<=4; ki++) {
    for (int kj=1; kj<=4; kj++) {
      bre.ddtab[ki][kj] = bre.rb2[ki][kj]/(double)(bre.ntab-2);
      bre.ddtab[kj][ki] = bre.ddtab[ki][kj];
      bre.rc = 0.0;
      for (int i=2; i<=bre.ntab-1; i++) {
	if (bre.ddtab[ki][kj]!=0.0) {
	  bre.rc = bre.rc + bre.ddtab[ki][kj];
	  double rsq = bre.rc*bre.rc;
	  fc = 0.0;
	  dfc = 0.0;
	  if (bre.rc<bre.rb2[ki][kj]) {
	    dtemp = bre.pid[ki][kj]*(bre.rc-bre.rb1[ki][kj]);
	    fc = (1.0+cos(dtemp))/2.0;
	    dfc = -bre.pid[ki][kj]/2.0*sin(dtemp); }
	  if (bre.rc<=bre.rb1[ki][kj]) {
	    fc = 1.0; dfc = 0.0; }
	  bre.tabfc[ki][kj][i] = fc;
	  bre.tabfc[kj][ki][i] = bre.tabfc[ki][kj][i];
	  bre.tabdfc[ki][kj][i] = dfc;
	  bre.tabdfc[kj][ki][i] = bre.tabdfc[ki][kj][i];
	  va = bre.ad[ki][kj]*exp(-bre.axl[ki][kj]*bre.rc);
	  dva = -bre.axl[ki][kj]*va;
	  vb = bre.bd[ki][kj]*exp(-bre.bxl[ki][kj]*bre.rc);
	  dvb = -bre.bxl[ki][kj]*vb;
	  vc = bre.cd[ki][kj]*exp(-bre.cxl[ki][kj]*bre.rc);
	  dvc = -bre.cxl[ki][kj]*vc;
	  vv = (va+vb+vc)/2.0;
	  dvv = (dva+dvb+dvc)/2.0;
	  bre.atable[ki][kj][i] = fc*vv;
	  bre.atable[kj][ki][i] = bre.atable[ki][kj][i];
	  bre.datable[ki][kj][i] = (fc*dvv+dfc*vv)/bre.rc;
	  bre.datable[kj][ki][i] = bre.datable[ki][kj][i];
	  ff1 = bre.dd[ki][kj]*exp(-bre.dxl[ki][kj]*bre.rc);
	  df1 = -bre.dxl[ki][kj]*ff1;
	  ff2 = (1.0+bre.ed[ki][kj]/bre.rc);
	  df2 = -bre.ed[ki][kj]/rsq;
	  vv = ff1*ff2;
	  dvm = (df1*ff2 + ff1*df2);
	  bre.rtable[ki][kj][i] = vv*fc;
	  bre.rtable[kj][ki][i] = bre.rtable[ki][kj][i];
	  bre.drtable[ki][kj][i] = -(fc*dvm+dfc*vv)/bre.rc;
	  bre.drtable[kj][ki][i] = bre.drtable[ki][kj][i];
	} else {
	  bre.tabfc[ki][kj][i]=0.0;
	  bre.tabfc[kj][ki][i] = bre.tabfc[ki][kj][i];
	  bre.tabdfc[ki][kj][i] = 0.0;
	  bre.tabdfc[kj][ki][i] = bre.tabdfc[ki][kj][i];
	  bre.atable[ki][kj][i] = 0.0;
	  bre.atable[kj][ki][i] = bre.atable[ki][kj][i];
	  bre.datable[ki][kj][i] = 0.0;
	  bre.datable[kj][ki][i] = bre.datable[ki][kj][i];
	  bre.rtable[ki][kj][i] = 0.0;
	  bre.rtable[kj][ki][i] = bre.rtable[ki][kj][i];
	  bre.drtable[ki][kj][i] = 0.0;
	  bre.drtable[kj][ki][i] = bre.drtable[ki][kj][i];
	}
      }
      bre.atable[ki][kj][1] = bre.atable[ki][kj][2];
      bre.atable[kj][ki][1] = bre.atable[ki][kj][1];
      bre.datable[ki][kj][1] = bre.datable[ki][kj][2];
      bre.datable[kj][ki][1] = bre.datable[ki][kj][1];
      bre.rtable[ki][kj][1] = bre.rtable[ki][kj][2];
      bre.rtable[kj][ki][1] = bre.rtable[ki][kj][1];
      bre.drtable[ki][kj][1] = bre.drtable[ki][kj][2];
      bre.drtable[kj][ki][1] = bre.drtable[ki][kj][1];
      bre.tabfc[ki][kj][1] = 0.0;
      bre.tabfc[kj][ki][1] = 0.0;
      bre.tabdfc[ki][kj][1] = 0.0;
      bre.tabdfc[kj][ki][1] = 0.0;

      bre.atable[ki][kj][bre.ntab] = 0.0;
      bre.atable[kj][ki][bre.ntab] = bre.atable[ki][kj][bre.ntab];
      bre.datable[ki][kj][bre.ntab] = 0.0;
      bre.datable[kj][ki][bre.ntab] = bre.datable[ki][kj][bre.ntab];
      bre.rtable[ki][kj][bre.ntab] = 0.0;
      bre.rtable[kj][ki][bre.ntab] = bre.rtable[ki][kj][bre.ntab];
      bre.drtable[ki][kj][bre.ntab] = 0.0;
      bre.drtable[kj][ki][bre.ntab] = bre.drtable[ki][kj][bre.ntab];
      bre.tabfc[ki][kj][bre.ntab] = 0.0;
      bre.tabfc[kj][ki][bre.ntab] = bre.tabfc[ki][kj][bre.ntab];
      bre.tabdfc[ki][kj][bre.ntab] = 0.0;
      bre.tabdfc[kj][ki][bre.ntab] = bre.tabdfc[ki][kj][bre.ntab];
    }
  }
}
void e_force_brenner()
{
  e_force_brenner(0);
}

void e_force_brenner(int mode)
{
  //
  // mode: 0 = Brenner (REBO), 1 = AIREBO (2nd generation Brenner)
  //
  //double rr, rr2, drx, dry, drz, vp0, v0;
  //int j, ix, iy, iz;
  double pi = 4.0*atan(1.0);
  double bolz = 1.380662, epsi=11604.5, avo=6.02205;
  double econv=(1.0/(avo*bolz*epsi))*1.0e7;
  double ri1, ri2, ri3, rr1, rr2, rr3;
  double cube1, cube2, cube3;
  //int kend;
  //bre.rc = 6.0e0; bre.rsq = bre.rc * bre.rc;
  bool force_bookkeep = false;

 TOP:
  
  if (bre.initialize) {
    std::cout<<"Initialize of the Brenner potential..";
    brealloc();
    bresetin(mode);
    breparam(mode);
    bremtable();
    //if (mode>0) { 
    //  bre.lj  = new double[bre.ljmesh];
    //  bre.ljd = new double[bre.ljmesh];
    //  airebo_lj_read(); }
    //breparam();
    for (int i=1; i<=atom.natom; i++) {
      //bre.ktype[i] = bre.kt[atom.getAtomNumber(atom.asp[i])];
      bre.ktype[i] = bre.kt[atom.anum[i]]; // H:ktype=2, C:ktype=1
      bre.noa[bre.ktype[i]] = bre.noa[bre.ktype[i]]+1;
    }
    if (bre.noa[3]!=0) { std::cout<<"\n### W A R N I N G ### Brenner-Si does not work"<<std::endl; }
    if (bre.noa[4]!=0) { std::cout<<"\n### W A R N I N G ### Brenner-Ge does not work"<<std::endl; }
    bre.initialize = false;
    force_bookkeep = true;
    std::cout<<" done."<<std::endl;
  } // end of bre.initialize
  
  for (int i=1; i<=atom.natom; i++) {
    //if ((atom.getAtomNumber(atom.asp[i])!=6)&&(atom.getAtomNumber(atom.asp[i])!=1)) {
    if ((atom.anum[i]!=6)&&(atom.anum[i]!=1)) {
      mdmotion = 0; printf("Brenner does not work for this atom\n"); break; } }

  bre.np=atom.natom;
  cube1=cell.hmat[0][0]/ang;
  cube2=cell.hmat[1][1]/ang;
  cube3=cell.hmat[2][2]/ang;
  //
  int nrpx = 0; int nrpy = 0; int nrpz = 0;
  if (cell.pbcx) { nrpx=(int)(book.frc*2/cell.hmat[0][0]+1); }
  if (cell.pbcy) { nrpy=(int)(book.frc*2/cell.hmat[1][1]+1); }
  if (cell.pbcz) { nrpz=(int)(book.frc*2/cell.hmat[2][2]+1); }
  //printf("%e %e  %e\n",book.frc*2, cell.hmat[0][0], book.frc*2/cell.hmat[0][0]);
  //printf("# of replica: %d x %d x %d\n",nrpx,nrpy,nrpz);
  //
  for (int i=1; i<=atom.natom; i++) {
    bre.r01[i] = atom.rx[i]/ang;
    bre.r02[i] = atom.ry[i]/ang;
    bre.r03[i] = atom.rz[i]/ang; }
  // Reset energy, force and others
  for (int i=1; i<=bre.np; i++) {
    atom.epot[i] = 0.0;
    atom.fx[i] = 0.0; atom.fy[i] = 0.0; atom.fz[i] = 0.0;  atom.epot[i] = 0.0;
    bre.rnp1[i] = 0.0; bre.rnp2[i] = 0.0; bre.rnp3[i] = 0.0; }
  // Reset stress
  cell.virx = 0.0; cell.viry = 0.0; cell.virz = 0.0;
  for (int i=0; i<3; i++) { for (int j=0; j<3; j++) { 
      cell.dmat[i][j] = 0.0;
      for (int ii=1; ii<=atom.natom; ii++) { atom.satom[ii][i][j] = 0.0; }  } }
  //if ((istep%book.nbk==0)||(istep==0)) { // Bookkeep
  if ((istep%book.nbk==0)||(istep==0)||force_bookkeep) { // Bookkeep
    if (force_bookkeep) { force_bookkeep = false; }
    //std::cout<<"Bookkeep for the Brenner potential..";
    int k=0;
    //printf("======= frc = %e\n",book.frc);
    for (int i=1; i<=bre.np; i++) {
      bre.nabors[i] = k+1;
      ri1 = bre.r01[i]; ri2 = bre.r02[i]; ri3 = bre.r03[i];
      int ki = bre.ktype[i];
      if (ki>=5) { continue; }
      for (int j=1; j<=bre.np; j++) {
	// if (i==j) { continue; } // removed 2015.12.17
	int kj = bre.ktype[j];
	if (kj>=5) { continue; }
	bre.rlis = bre.rlist[ki][kj];
	bre.rlis = (book.frc/ang)*(book.frc/ang);
	double rsq = 0.0;
	rr1 = ri1-bre.r01[j]; rr2 = ri2-bre.r02[j]; rr3 = ri3-bre.r03[j];
	//if (cell.pbcx) { rr1 = rr1-cube1*anint(rr1/cube1); }
	//if (cell.pbcy) { rr2 = rr2-cube2*anint(rr2/cube2); }
	//if (cell.pbcz) { rr3 = rr3-cube3*anint(rr3/cube3); }
	double xx=0,yy=0,zz=0;
	int xx0=0,yy0=0,zz0=0;
	int xx1=0,yy1=0,zz1=0;
	if (cell.pbcx) { //xx = -(double)anint(rr1/cube1);
	  xx0 = (-book.frc/ang-rr1)/cube1; xx1 = (book.frc/ang-rr1)/cube1;
	}
	if (cell.pbcy) { //yy = -(double)anint(rr2/cube2);
	  yy0 = (-book.frc/ang-rr2)/cube2; yy1 = (book.frc/ang-rr2)/cube2;
	}
	if (cell.pbcz) { //zz = -(double)anint(rr3/cube3);
	  zz0 = (-book.frc/ang-rr3)/cube3; zz1 = (book.frc/ang-rr3)/cube3;
	}
	//printf("%d %d   %d %d   %d %d\n",xx0,xx1,yy0,yy1,zz0,zz1);

	/*
	rr1 += xx*cell.hmat[0][0]/ang + yy*cell.hmat[0][1]/ang + zz*cell.hmat[0][2]/ang;
	rr2 += xx*cell.hmat[1][0]/ang + yy*cell.hmat[1][1]/ang + zz*cell.hmat[1][2]/ang;
	rr3 += xx*cell.hmat[2][0]/ang + yy*cell.hmat[2][1]/ang + zz*cell.hmat[2][2]/ang;
	rsq = rsq + rr1*rr1; //if (rsq > bre.rlis) { continue; }
	rsq = rsq + rr2*rr2; //if (rsq > bre.rlis) { continue; }
	rsq = rsq + rr3*rr3; //if (rsq > bre.rlis) { continue; }
	if (rsq > bre.rlis) { continue; }
	k++;
	bre.list[k] = j; bre.ivct2b[k] = i; bre.jvct2b[k] = j;
	*/

	for (int ix=xx0; ix<=xx1; ix++) { double xx=(double)ix;
	for (int iy=yy0; iy<=yy1; iy++) { double yy=(double)iy;
	for (int iz=zz0; iz<=zz1; iz++) { double zz=(double)iz;
	  rr1 = ri1-bre.r01[j]; rr2 = ri2-bre.r02[j]; rr3 = ri3-bre.r03[j];
	  rr1 += xx*cell.hmat[0][0]/ang + yy*cell.hmat[0][1]/ang + zz*cell.hmat[0][2]/ang;
	  rr2 += xx*cell.hmat[1][0]/ang + yy*cell.hmat[1][1]/ang + zz*cell.hmat[1][2]/ang;
	  rr3 += xx*cell.hmat[2][0]/ang + yy*cell.hmat[2][1]/ang + zz*cell.hmat[2][2]/ang;
	  rsq = rr1*rr1 + rr2*rr2 + rr3*rr3;
	  if (rsq<0.0001) { continue; } // added 2015.12.17
	  if (rsq > bre.rlis) { continue; }
	  k++;
	  if (k<=bre.nlmax) {
	    bre.list[k] = j; bre.ivct2b[k] = i; bre.jvct2b[k] = j;
	    bre.rep1[k] = ix; bre.rep2[k] = iy; bre.rep3[k] = iz;
	  }
	} } }
	
      } // end of loop j
    } // end of loop i
    bre.nabors[bre.np+1] = k+1;
    bre.kend = k;
    //printf("=== Brenner: kend=%d, nlmax=%d\n",bre.kend,bre.nlmax);
    //std::cout<<bre.kend<<std::endl;
    if (bre.kend> bre.nlmax) {
      //std::cout<<"Error: kend exceeds nlmax"<<std::endl; exit(0); }
      printf("Brenner: kend (%d) exceeds nlmax (%d).\n",bre.kend,bre.nlmax);
      printf("  Now kend is increased to nlmax and arrays are reallocated.\n");
      bre.nlmax = bre.kend; bre.initialize=true; goto TOP;}
    //std::cout<<" done."<<std::endl;
    bond_set();
  } // Bookkeep end

  for (int k=1; k<=bre.kend; k++) {
    int i=bre.ivct2b[k]; int j=bre.jvct2b[k];
    int ki=bre.ktype[i]; int kj=bre.ktype[j];
    bre.lcheck[k]=0;
    double rsq=0.0;

    rr1 = bre.r01[i]-bre.r01[j]; //if (cell.pbcx) { rr1 = rr1-cube1*anint(rr1/cube1); }
    rr2 = bre.r02[i]-bre.r02[j]; //if (cell.pbcy) { rr2 = rr2-cube2*anint(rr2/cube2); }
    rr3 = bre.r03[i]-bre.r03[j]; //if (cell.pbcz) { rr3 = rr3-cube3*anint(rr3/cube3); }
    double xx=0,yy=0,zz=0;
    /*
    if (cell.pbcx) { xx = -(double)anint(rr1/cube1); }
    if (cell.pbcy) { yy = -(double)anint(rr2/cube2); }
    if (cell.pbcz) { zz = -(double)anint(rr3/cube3); }
    */
    xx = (double)bre.rep1[k]; yy = (double)bre.rep2[k]; zz = (double)bre.rep3[k];
    rr1 += xx*cell.hmat[0][0]/ang + yy*cell.hmat[0][1]/ang + zz*cell.hmat[0][2]/ang;
    rr2 += xx*cell.hmat[1][0]/ang + yy*cell.hmat[1][1]/ang + zz*cell.hmat[1][2]/ang;
    rr3 += xx*cell.hmat[2][0]/ang + yy*cell.hmat[2][1]/ang + zz*cell.hmat[2][2]/ang;
    rsq = rsq+rr1*rr1+rr2*rr2+rr3*rr3;
    bre.cor1[k] = rr1; bre.cor2[k] = rr2; bre.cor3[k] = rr3; 

    if (rsq > bre.rmax[ki][kj]) { continue; }
    if ((kj<=2)&&(ki<=2)) { bre.lcheck[k]=1; }
    if ((kj>=3)&&(ki>=3)) { bre.lcheck[k]=2; }
    bre.rc = sqrt(rsq);
    double rt = bre.rc/bre.ddtab[ki][kj];
    int it = (std::min)(int(rt)+1, bre.ntab-1);    
    bre.rcor[k] = bre.rc;

    // Modified by M. Satoh, Aug. 2014
    /*double xrt = rt-(double)it+1.0;
    bre.bww[k] = bre.tabfc[ki][kj][it]
      + (bre.tabfc[ki][kj][it+1]-bre.tabfc[ki][kj][it])*xrt;
    bre.dww[k] = bre.tabdfc[ki][kj][it]
      + (bre.tabdfc[ki][kj][it+1]-bre.tabdfc[ki][kj][it])*xrt;
    bre.exx1[k] = bre.atable[ki][kj][it]
      + (bre.atable[ki][kj][it+1]-bre.atable[ki][kj][it])*xrt;
    bre.dexx1[k] = bre.datable[ki][kj][it]
      + (bre.datable[ki][kj][it+1]-bre.datable[ki][kj][it])*xrt;
    if (i>=j) { continue; }
    double vv =bre.rtable[ki][kj][it] + (bre.rtable[ki][kj][it+1]-bre.rtable[ki][kj][it])*xrt;
    double rp =bre.drtable[ki][kj][it] + (bre.drtable[ki][kj][it+1]-bre.drtable[ki][kj][it])*xrt;*/
    double wij = 1.0, dwij = 0.0,vv=0.0,rp=0.0;
    bre.bww[k]=0.0;
    bre.dww[k]=0.0; 
    bre.exx1[k] = 0.0;
    bre.dexx1[k] = 0.0;
    if (bre.rc>bre.rb2[ki][kj]) {wij=0.0;}
    if (bre.rc>bre.rb1[ki][kj]) {
    double dtemp = bre.pid[ki][kj]*(bre.rc-bre.rb1[ki][kj]);    
    wij = (1.0+cos(dtemp))/2.0;
    dwij = -bre.pid[ki][kj]/2.0*sin(dtemp);}
    if(wij>0){
    double va = bre.ad[ki][kj]*exp(-bre.axl[ki][kj]*bre.rc);
    double dva = -bre.axl[ki][kj]*va;
    double vb = bre.bd[ki][kj]*exp(-bre.bxl[ki][kj]*bre.rc);
    double dvb = -bre.bxl[ki][kj]*vb;
    double vc = bre.cd[ki][kj]*exp(-bre.cxl[ki][kj]*bre.rc);
    double dvc = -bre.cxl[ki][kj]*vc;
    double vd = (va+vb+vc)/2.0;
    double dvd = (dva+dvb+dvc)/2.0;
    bre.bww[k]=wij;
    bre.dww[k]=dwij; 
    bre.exx1[k] = wij*vd;
    bre.dexx1[k] = (wij*dvd+dwij*vd)/bre.rc;}
    if (i>=j) { continue; }
    if(wij>0){
    double f1 = bre.dd[ki][kj]*exp(-bre.dxl[ki][kj]*bre.rc);
    double df1 = -bre.dxl[ki][kj]*f1;
    double f2 = (1.0+bre.ed[ki][kj]/bre.rc);
    double df2 = -bre.ed[ki][kj]/rsq;
    vv = wij*f1*f2;
    rp = -(dwij*f1*f2+wij*(f1*df2+df1*f2))/bre.rc;}
    // End of modification

    bre.tote = bre.tote + vv;
    atom.epot[i] = atom.epot[i] + vv/2.0;
    atom.epot[j] = atom.epot[j] + vv/2.0;
    bre.rpp1[k] = rp * rr1; bre.rpp2[k] = rp * rr2; bre.rpp3[k] = rp * rr3;
    // stress
    rp /= 2;
    atom.satom[i][0][0] -= rp * rr1 * rr1;
    atom.satom[i][0][1] -= rp * rr1 * rr2;
    atom.satom[i][0][2] -= rp * rr1 * rr3;
    atom.satom[i][1][1] -= rp * rr2 * rr2;
    atom.satom[i][1][2] -= rp * rr2 * rr3;
    atom.satom[i][2][2] -= rp * rr3 * rr3;
    atom.satom[j][0][0] -= rp * rr1 * rr1;
    atom.satom[j][0][1] -= rp * rr1 * rr2;
    atom.satom[j][0][2] -= rp * rr1 * rr3;
    atom.satom[j][1][1] -= rp * rr2 * rr2;
    atom.satom[j][1][2] -= rp * rr2 * rr3;
    atom.satom[j][2][2] -= rp * rr3 * rr3;
  } // End of loop k (320)
  for (int k=1; k<=bre.kend; k++) {
    if (bre.lcheck[k]==0) { continue; }
    int i = bre.ivct2b[k];
    int j = bre.jvct2b[k];
    if (i>=j) { continue; }
    bre.rnp1[i] = bre.rnp1[i]+bre.rpp1[k];
    bre.rnp1[j] = bre.rnp1[j]-bre.rpp1[k];
    bre.rnp2[i] = bre.rnp2[i]+bre.rpp2[k];
    bre.rnp2[j] = bre.rnp2[j]-bre.rpp2[k];
    bre.rnp3[i] = bre.rnp3[i]+bre.rpp3[k];
    bre.rnp3[j] = bre.rnp3[j]-bre.rpp3[k];
  } // End of loop k (321)
  if ((bre.noa[1]+bre.noa[2])!=0) { pibond(); }
  if ((bre.noa[3]+bre.noa[4])!=0) {
    std::cout<<"Brenner-Si/Ge contains serious bug!"<<std::endl; }
  if (mode > 0) {
    if (relax_accel==1) {
      airebo_lj(1);
    } else {
      airebo_lj();
    }
  }
  for (int i=1; i<=bre.np; i++) {
    atom.fx[i] = bre.rnp1[i]*eV/ang;
    atom.fy[i] = bre.rnp2[i]*eV/ang;
    atom.fz[i] = bre.rnp3[i]*eV/ang;
    atom.epot[i] = atom.epot[i]*eV;
    ///
    if (relax_accel==1) {
      if (atom.lock[i]) { atom.epot[i] = atom.epot_p[i]; }
    }
    ///
    for (int ii=0; ii<3; ii++) { for (int jj=ii; jj<3; jj++) {
	atom.satom[i][ii][jj] *= eV/ang*ang;
	if (ii!=jj) { atom.satom[i][jj][ii] = atom.satom[i][ii][jj]; } } }
  }
  // Virial
  for (int ii=1; ii<=atom.natom; ii++) {
    for (int i=0; i<3; i++) { for (int j=0; j<3; j++) { 
	cell.dmat[i][j] = cell.dmat[i][j] - atom.satom[ii][i][j];
      } } }
  cell.virx = cell.dmat[0][0]; cell.viry = cell.dmat[1][1]; cell.virz = cell.dmat[2][2];
  cell.volume = cell.Getvolume();
  for (int ii=1; ii<=atom.natom; ii++) {
    for (int i=0; i<3; i++) { for (int j=0; j<3; j++) {
	atom.satom[ii][i][j] = atom.satom[ii][i][j] * (double)atom.natom / cell.volume;
      } } }
  atom.epotsum=0.0;
  for (int i=1; i<=atom.natom; i++) {
    atom.epotsum += atom.epot[i]; }
  // rcut
  rcut = 2.0;
  //if (mode > 0) { rcut = 3.82; } // AIREBO
  if (mode > 0) { rcut = 10.00; } // AIREBO
} // end of e_force_brenner

void pibond()
{
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
      if (bre.lcheck[j]!=1) { continue; }
      int jn = bre.list[j];
      if (bre.ktype[jn]==1) { bre.xhc1[i] = bre.xhc1[i]+bre.bww[j]; } //C
      if (bre.ktype[jn]==2) { bre.xhc2[i] = bre.xhc2[i]+bre.bww[j]; } //H
      if (bre.xhc1[i]>12) { printf("??????\n"); }
    } }
  // sum over bonds between atoms i and j
  for (int i=1; i<=bre.np; i++) {
    int jbegin = bre.nabors[i];
    int jend = bre.nabors[i+1]-1;
    if (jbegin > jend) { continue; }
    int ki = bre.ktype[i];
    for (int j=jbegin; j<=jend; j++) {
      if (bre.lcheck[j]!=1) { continue; }
      int jn = bre.list[j];
      if (i >= jn) { continue; }
      cj[1] = bre.cor1[j]; cj[2] = bre.cor2[j]; cj[3] = bre.cor3[j];
      double sij = bre.rcor[j];
      double rsqij = sij*sij;
      int kj = bre.ktype[jn];
      int kikj = ki+kj;

      // i side of bond
      int nk = 0;
      double xsij = 0.0;
      double ssumk = 0.0;
      double conk = 0.0;
      xni[1] = bre.xhc1[i]; xni[2] = bre.xhc2[i];
      xni[kj] = xni[kj]-bre.bww[j];
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
	  double rsq3 = s3*s3;
	  double rsq2 = 0.0;
	  xk[nk][1] = bre.cor1[k]-cj[1];
	  xk[nk][2] = bre.cor2[k]-cj[2];
	  xk[nk][3] = bre.cor3[k]-cj[3];
	  rsq2 = rsq2+xk[nk][1]*xk[nk][1];
	  rsq2 = rsq2+xk[nk][2]*xk[nk][2];
	  rsq2 = rsq2+xk[nk][3]*xk[nk][3];
	  double ss = 2.0*sij*s3;
	  double rr = rsqij-rsq3;
	  double costh = (rsqij+rsq3-rsq2)/ss;
	  if (costh > 1.0) { costh = 1.0; }
	  if (costh <-1.0) { costh =-1.0; }
	  cosk[nk] = costh;
	  sink[nk] = sqrt(1.0-costh*costh);
	  if (acos(costh) > pi) { sink[nk] = -sink[nk]; }
	  int ig = bre.igc[int(-costh*12)+13];
	  double gangle; double dgdthet;
	  if (ki == 1) {
	    if (ig != 4) {
	      gangle = bre.spgc[1][ig]+bre.spgc[2][ig]*costh;
	      dgdthet = bre.spgc[2][ig];
	      for (int jj=3; jj<=6; jj++) {
		gangle = gangle+bre.spgc[jj][ig]*pow(costh,jj-1);
		dgdthet = dgdthet+bre.spgc[jj][ig]*(double)(jj-1)*pow(costh,jj-2); }
	    } else {
	      ali = 0.0; dali = 0.0;
	      if (qi < bre.xqm) {
		ali = 1.0;
		if (qi > bre.att) {
		  double dtemp = bre.pq*(qi-bre.att);
		  ali = (1.0+cos(dtemp))/2.0;
		  dali = -bre.pq/2.0*sin(dtemp); } }
	      gangle = bre.spgc[1][ig]+bre.spgc[2][ig]*costh;
	      dgdthet = bre.spgc[2][ig];
	      ig1=ig+1;
	      double gangle1 = bre.spgc[1][ig1]+bre.spgc[2][ig1]*costh;
	      double dgdthet1 = bre.spgc[2][ig1];
	      for (int jj=3; jj<=6; jj++) {
		gangle = gangle+bre.spgc[jj][ig]*pow(costh,jj-1);
		dgdthet = dgdthet+bre.spgc[jj][ig]*(double)(jj-1)*pow(costh,jj-2);
		gangle1 = gangle1+bre.spgc[jj][ig1]*pow(costh,jj-1);
		dgdthet1 = dgdthet1+bre.spgc[jj][ig1]*(double)(jj-1)*pow(costh,jj-2); }
	      daldik = dali*(gangle1-gangle);
	      gangle = gangle+ali*(gangle1-gangle);
	      dgdthet = dgdthet+ali*(dgdthet1-dgdthet);
	    } // if (ig != 4)
	  } else { // if (ki == 1)
	    ig = bre.igh[int(-costh*12.0)+13];
	    gangle = bre.spgh[1][ig]+bre.spgh[2][ig]*costh;
	    dgdthet = bre.spgh[2][ig];
	    for (int jj=3; jj<=6; jj++) {
	      gangle = gangle+bre.spgh[jj][ig]*pow(costh,jj-1);
	      dgdthet = dgdthet+bre.spgh[jj][ig]*(double)(jj-1)*pow(costh,jj-2); }
	  } // if (ki == 1)
	  double fc = bre.bww[k]; double dfc = bre.dww[k];
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
	  
	  double dctdjk = -2.0/ss;
	  double dctdij = (rr+rsq2)/(ss*rsqij);
	  double dctdik = (-rr+rsq2)/(ss*rsq3);
	  dctjk[nk] = dctdjk;
	  dctij[nk] = dctdij;
	  dctik[nk] = dctdik;
	  double gs = gangle*exx;
	  ssumk = ssumk+fc*gs;
	  double xtemp = fc*exx*dgdthet;
	  double gfx = gs*fc*bre.xdb[ki][kj][kk];
	  xsij = xsij+xtemp*dctdij+gfx/sij;
	  xsik[nk] = (gs*dfc-gfx)/s3+xtemp*dctdik;
	  sdalik = sdalik+exx*fc*daldik;
	  xsjk[nk] = xtemp*dctdjk;
	} // loop k
      } // if (jbegin != jend )
      
      // j side of bond
      int nl = 0;
      double xsji = 0.0;
      double ssuml = 0.0;
      double conl = 0.0;
      int lbegin = bre.nabors[jn];
      int lend = bre.nabors[jn+1]-1;
      xnj[1] = bre.xhc1[jn];
      xnj[2] = bre.xhc2[jn];
      xnj[ki] = xnj[ki] - bre.bww[j];
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
	  double rsq3 = s3*s3;
	  double rsq2 = 0.0;
	  xl[nl][1] = bre.cor1[l]+cj[1];
	  xl[nl][2] = bre.cor2[l]+cj[2];
	  xl[nl][3] = bre.cor3[l]+cj[3];
	  rsq2 = rsq2+xl[nl][1]*xl[nl][1];
	  rsq2 = rsq2+xl[nl][2]*xl[nl][2];
	  rsq2 = rsq2+xl[nl][3]*xl[nl][3];
	  double ss = 2.0*sij*s3;
	  double rr = rsqij-rsq3;
	  double costh = (rsqij+rsq3-rsq2)/ss;
	  if (costh > 1.0) { costh = 1.0; }
	  if (costh <-1.0) { costh =-1.0; }
	  cosl[nl] = costh;
	  sinl[nl] = sqrt(1.0-costh*costh);
	  if (acos(costh) > pi) { sinl[nl] = -sinl[nl]; }
	  double gangle; double dgdthet;
	  if (kj == 1) {
	    int ig = bre.igc[int(-costh*12)+13];
	    if (ig != 4) {
	      gangle = bre.spgc[1][ig]+bre.spgc[2][ig]*costh;
	      dgdthet = bre.spgc[2][ig];
	      for (int jj=3; jj<=6; jj++) {
		gangle = gangle+bre.spgc[jj][ig]*pow(costh,jj-1);
		dgdthet = dgdthet+bre.spgc[jj][ig]*(double)(jj-1)*pow(costh,jj-2); }
	    } else {
	      alj = 0.0; dalj = 0.0;
	      if (qj < bre.xqm) {
		alj = 1.0; // 'double' removed 2014.11.06
		if (qj > bre.att) {
		  double dtemp = bre.pq*(qj-bre.att);
		  alj = (1.0+cos(dtemp))/2.0;
		  dalj = -bre.pq/2.0*sin(dtemp); } }
	      gangle = bre.spgc[1][ig]+bre.spgc[2][ig]*costh;
	      dgdthet = bre.spgc[2][ig];
	      ig1=ig+1;
	      double gangle1 = bre.spgc[1][ig1]+bre.spgc[2][ig1]*costh;
	      double dgdthet1 = bre.spgc[2][ig1];
	      for (int jj=3; jj<=6; jj++) {
		gangle = gangle+bre.spgc[jj][ig]*pow(costh,jj-1);
		dgdthet = dgdthet+bre.spgc[jj][ig]*(double)(jj-1)*pow(costh,jj-2);
		// ig -> ig1 in the following two lines 2014.11.06
		gangle1 = gangle1+bre.spgc[jj][ig1]*pow(costh,jj-1);
		dgdthet1 = dgdthet1+bre.spgc[jj][ig1]*(double)(jj-1)*pow(costh,jj-2); }
	      daldjl = dalj*(gangle1-gangle);
	      gangle = gangle+alj*(gangle1-gangle);
	      dgdthet = dgdthet+alj*(dgdthet1-dgdthet);
	    } // if (ig != 4) (2)
	  } else { // if (ki == 1) (2)
	    int ig = bre.igh[int(-costh*12.0)+13];
	    gangle = bre.spgh[1][ig]+bre.spgh[2][ig]*costh;
	    dgdthet = bre.spgh[2][ig];
	    for (int jj=3; jj<=6; jj++) {
	      gangle = gangle+bre.spgh[jj][ig]*pow(costh,jj-1);
	      dgdthet = dgdthet+bre.spgh[jj][ig]*(double)(jj-1)*pow(costh,jj-2); }
	  } // if (ki == 1) (2)
	  double fc = bre.bww[l]; double dfc = bre.dww[l];
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

	  double dctdil = -2.0/ss;
	  double dctdji = (rr+rsq2)/(ss*rsqij);
	  double dctdjl = (-rr+rsq2)/(ss*rsq3);
	  dctil[nl] = dctdil;
	  dctji[nl] = dctdji;
	  dctjl[nl] = dctdjl;
	  double gs = gangle*exx;
	  ssuml = ssuml+fc*gs;
	  double xtemp = fc*exx*dgdthet;
	  double gfx = gs*fc*bre.xdb[kj][ki][kl];
	  xsji = xsji+xtemp*dctdji+gfx/sij;
	  xsjl[nl] = (gs*dfc-gfx)/s3+xtemp*dctdjl;
	  sdaljl = sdaljl+exx*fc*daldjl;
	  xsil[nl] = xtemp*dctdil;
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

      double dij = (1.0+exnij+ssumk);
      double bij = pow(dij,-0.5);
      double dji = (1.0+exnji+ssuml);
      double bji = pow(dji,-0.5);
      double dbdzi = -0.5*bij/dij;
      double dbdzj = -0.5*bji/dji;
      double vatt = bre.exx1[j];
      double dradi = 0.0;
      double dradj = 0.0;
      double drdc = 0.0;
      double conjug = 1.0+conk*conk+conl*conl;
      double xnt1 = xni[1]+xni[2]-1.0;
      double xnt2 = xnj[1]+xnj[2]-1.0;
      double rad = 0.0; // ???
      radic(ki,kj,xnt1,xnt2,conjug,rad,dradi,dradj,drdc);
      double btot = bji+bij+rad;

      //goto OUT;
      // dihedral terms
      if (kikj == bre.ndihed) {
	//dbtori = 0.0; dbtorj = 0.0; dbtorc = 0.0;
	double datori = 0.0; double datorj = 0.0; double datorc = 0.0; // ???
	double btor = 0.0;
	double ator = 0.0; // ???
	tor(xnt1,xnt2,conjug,ator,datori,datorj,datorc);
	if (fabs(ator) <= 1.0e-08) { goto OUT; }
	
	if ((jbegin != jend)&&(lbegin != lend)) {
	  nk = 0;
	  for (int k=jbegin; k<=jend; k++) {
	    if (k == j) { continue; } // go to 220
	    if (bre.lcheck[k] != 1) { continue; } // go to 220
	    nk++;
	    if (fabs(sink[nk]) < 1.0e-1) { continue; } // go to 220
	    double sink2 = sink[nk]*sink[nk];
	    int kn = bre.list[k];
	    ck[1] = bre.cor1[k]; ck[2] = bre.cor2[k]; ck[3] = bre.cor3[k];
	    double rck = bre.rcor[k];
	    double fck; double dfck;
	    if (bre.ktype[kn] == 2) {
	      fck = 1.0;
	      dfck = 0.0;
	      if (rck >= 1.60) { continue; } // go to 220
	      if (rck >= 1.30) {
		double dtemp = bre.pidt*(rck-1.30);
		fck = (1.0+cos(dtemp))/2.0;
		dfck = -bre.pidt/2.0*sin(dtemp); }
	    } else {
	      fck = bre.bww[k];
	      dfck = bre.dww[k];
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
	      double fcl; double dfcl;
	      if (bre.ktype[ln] == 2) {
		fcl = 1.0;
		dfcl = 0.0;
		if (rcl >= 1.60) { continue; } // go to 210
		if (rcl >= 1.30) {
		  double dtemp = bre.pidt*(rcl-1.30);
		  fcl = (1.0+cos(dtemp))/2.0;
		  dfcl = -bre.pidt/2.0*sin(dtemp);
		}
	      } else {
		fcl = bre.bww[l];
		dfcl = bre.dww[l];
	      }
	      
	      double t1=rck*rcl*sij*sij*sink[nk]*sinl[nl];
	      double dt1dik=1.0/rck/rck-dctik[nk]/sink2*cosk[nk];
	      double dt1djk=-dctjk[nk]/sink2*cosk[nk];
	      double dt1djl=1.0/rcl/rcl-dctjl[nl]/sinl2*cosl[nl];
	      double dt1dil=-dctil[nl]/sinl2*cosl[nl];
	      double dt1dij=2.0/sij/sij
		-dctij[nk]/sink2*cosk[nk]-dctji[nl]/sinl2*cosl[nl];
	      double crkx=ck[2]*cj[3]-cj[2]*ck[3];
	      double crlx=cj[2]*cl[3]-cl[2]*cj[3];
	      double crky=ck[3]*cj[1]-cj[3]*ck[1];
	      double crly=cj[3]*cl[1]-cl[3]*cj[1];
	      double crkz=ck[1]*cj[2]-cj[1]*ck[2];
	      double crlz=cj[1]*cl[2]-cl[1]*cj[2];
	      double t2=crkx*crlx+crky*crly+crkz*crlz;
	      double cw=t2/t1;
	      double bt=(1.0-cw*cw);
	      btor=btor+bt*fck*fcl;
	      dt2dik[1]=-cj[3]*crly+cj[2]*crlz;
	      dt2dik[2]=-cj[1]*crlz+cj[3]*crlx;
	      dt2dik[3]=-cj[2]*crlx+cj[1]*crly;
	      dt2djl[1]=-cj[2]*crkz+cj[3]*crky;
	      dt2djl[2]=-cj[3]*crkx+cj[1]*crkz;
	      dt2djl[3]=-cj[1]*crky+cj[2]*crkx;
	      dt2dij[1]=ck[3]*crly-cl[3]*crky-ck[2]*crlz+cl[2]*crkz;
	      dt2dij[2]=ck[1]*crlz-cl[1]*crkz-ck[3]*crlx+cl[3]*crkx;
	      dt2dij[3]=ck[2]*crlx-cl[2]*crkx-ck[1]*crly+cl[1]*crky;
	      double aa=-vatt*2.0*cw/t1*ator*fcl*fck;
	      double aaa1=vatt*bt*ator;
	      double at2=aa*t2;
	      double rp1=-dt1dij*at2;
	      double rp2=-dt1dik*at2+aaa1*fcl*dfck/rck;
	      double rp3=-dt1djl*at2+aaa1*fck*dfcl/rcl;
	      double rp4=-dt1djk*at2;
	      double rp5=-dt1dil*at2;
	      double rep; double repx,repy,repz;
	      rep=rp1*cj[1]+aa*dt2dij[1];//
	      bre.rnp1[i]=bre.rnp1[i]+rep; bre.rnp1[jn]=bre.rnp1[jn]-rep;
	      rep=rp1*cj[2]+aa*dt2dij[2];
	      bre.rnp2[i]=bre.rnp2[i]+rep; bre.rnp2[jn]=bre.rnp2[jn]-rep;
	      rep=rp1*cj[3]+aa*dt2dij[3];
	      bre.rnp3[i]=bre.rnp3[i]+rep; bre.rnp3[jn]=bre.rnp3[jn]-rep;
	      // stress
	      atom.satom[i][0][0] -= (rp1*cj[1]*cj[1]+aa*dt2dij[1]*cj[1])/2;
	      atom.satom[i][0][1] -= (rp1*cj[1]*cj[2]+aa*dt2dij[1]*cj[2])/2;
	      atom.satom[i][0][2] -= (rp1*cj[1]*cj[3]+aa*dt2dij[1]*cj[3])/2;
	      atom.satom[i][1][1] -= (rp1*cj[2]*cj[2]+aa*dt2dij[2]*cj[2])/2;
	      atom.satom[i][1][2] -= (rp1*cj[2]*cj[3]+aa*dt2dij[2]*cj[3])/2;
	      atom.satom[i][2][2] -= (rp1*cj[3]*cj[3]+aa*dt2dij[3]*cj[3])/2;
	      atom.satom[jn][0][0] -= (rp1*cj[1]*cj[1]+aa*dt2dij[1]*cj[1])/2;
	      atom.satom[jn][0][1] -= (rp1*cj[1]*cj[2]+aa*dt2dij[1]*cj[2])/2;
	      atom.satom[jn][0][2] -= (rp1*cj[1]*cj[3]+aa*dt2dij[1]*cj[3])/2;
	      atom.satom[jn][1][1] -= (rp1*cj[2]*cj[2]+aa*dt2dij[2]*cj[2])/2;
	      atom.satom[jn][1][2] -= (rp1*cj[2]*cj[3]+aa*dt2dij[2]*cj[3])/2;
	      atom.satom[jn][2][2] -= (rp1*cj[3]*cj[3]+aa*dt2dij[3]*cj[3])/2;
	      rep=rp2*ck[1] + aa*dt2dik[1];//
	      bre.rnp1[i]=bre.rnp1[i]+rep; bre.rnp1[kn]=bre.rnp1[kn]-rep;
	      rep=rp2*ck[2] + aa*dt2dik[2];
	      bre.rnp2[i]=bre.rnp2[i]+rep; bre.rnp2[kn]=bre.rnp2[kn]-rep;
	      rep=rp2*ck[3] + aa*dt2dik[3];
	      bre.rnp3[i]=bre.rnp3[i]+rep; bre.rnp3[kn]=bre.rnp3[kn]-rep;
	      // stress
	      atom.satom[i][0][0] -= (rp2*ck[1]*ck[1]+aa*dt2dik[1]*ck[1])/2;
	      atom.satom[i][0][1] -= (rp2*ck[1]*ck[2]+aa*dt2dik[1]*ck[2])/2;
	      atom.satom[i][0][2] -= (rp2*ck[1]*ck[3]+aa*dt2dik[1]*ck[3])/2;
	      atom.satom[i][1][1] -= (rp2*ck[2]*ck[2]+aa*dt2dik[2]*ck[2])/2;
	      atom.satom[i][1][2] -= (rp2*ck[2]*ck[3]+aa*dt2dik[2]*ck[3])/2;
	      atom.satom[i][2][2] -= (rp2*ck[3]*ck[3]+aa*dt2dik[3]*ck[3])/2;
	      atom.satom[kn][0][0] -= (rp2*ck[1]*ck[1]+aa*dt2dik[1]*ck[1])/2;
	      atom.satom[kn][0][1] -= (rp2*ck[1]*ck[2]+aa*dt2dik[1]*ck[2])/2;
	      atom.satom[kn][0][2] -= (rp2*ck[1]*ck[3]+aa*dt2dik[1]*ck[3])/2;
	      atom.satom[kn][1][1] -= (rp2*ck[2]*ck[2]+aa*dt2dik[2]*ck[2])/2;
	      atom.satom[kn][1][2] -= (rp2*ck[2]*ck[3]+aa*dt2dik[2]*ck[3])/2;
	      atom.satom[kn][2][2] -= (rp2*ck[3]*ck[3]+aa*dt2dik[3]*ck[3])/2;
	      rep=rp3*cl[1] + aa*dt2djl[1];//
	      bre.rnp1[jn]=bre.rnp1[jn]+rep; bre.rnp1[ln]=bre.rnp1[ln]-rep;
	      rep=rp3*cl[2] + aa*dt2djl[2];
	      bre.rnp2[jn]=bre.rnp2[jn]+rep; bre.rnp2[ln]=bre.rnp2[ln]-rep;
	      rep=rp3*cl[3] + aa*dt2djl[3];
	      bre.rnp3[jn]=bre.rnp3[jn]+rep; bre.rnp3[ln]=bre.rnp3[ln]-rep;
	      // stress
	      atom.satom[jn][0][0] -= (rp3*cl[1]*cl[1]+aa*dt2djl[1]*cl[1])/2;
	      atom.satom[jn][0][1] -= (rp3*cl[1]*cl[2]+aa*dt2djl[1]*cl[2])/2;
	      atom.satom[jn][0][2] -= (rp3*cl[1]*cl[3]+aa*dt2djl[1]*cl[3])/2;
	      atom.satom[jn][1][1] -= (rp3*cl[2]*cl[2]+aa*dt2djl[2]*cl[2])/2;
	      atom.satom[jn][1][2] -= (rp3*cl[2]*cl[3]+aa*dt2djl[2]*cl[3])/2;
	      atom.satom[jn][2][2] -= (rp3*cl[3]*cl[3]+aa*dt2djl[3]*cl[3])/2;
	      atom.satom[ln][0][0] -= (rp3*cl[1]*cl[1]+aa*dt2djl[1]*cl[1])/2;
	      atom.satom[ln][0][1] -= (rp3*cl[1]*cl[2]+aa*dt2djl[1]*cl[2])/2;
	      atom.satom[ln][0][2] -= (rp3*cl[1]*cl[3]+aa*dt2djl[1]*cl[3])/2;
	      atom.satom[ln][1][1] -= (rp3*cl[2]*cl[2]+aa*dt2djl[2]*cl[2])/2;
	      atom.satom[ln][1][2] -= (rp3*cl[2]*cl[3]+aa*dt2djl[2]*cl[3])/2;
	      atom.satom[ln][2][2] -= (rp3*cl[3]*cl[3]+aa*dt2djl[3]*cl[3])/2;
	      rep=rp4*xk[nk][1];//
	      bre.rnp1[jn]=bre.rnp1[jn]+rep; bre.rnp1[kn]=bre.rnp1[kn]-rep;
	      rep=rp4*xk[nk][2];
	      bre.rnp2[jn]=bre.rnp2[jn]+rep; bre.rnp2[kn]=bre.rnp2[kn]-rep;
	      rep=rp4*xk[nk][3];
	      bre.rnp3[jn]=bre.rnp3[jn]+rep; bre.rnp3[kn]=bre.rnp3[kn]-rep;
	      // stress
	      atom.satom[jn][0][0] -= (rp4*xk[nk][1]*xk[nk][1])/2;
	      atom.satom[jn][0][1] -= (rp4*xk[nk][1]*xk[nk][2])/2;
	      atom.satom[jn][0][2] -= (rp4*xk[nk][1]*xk[nk][3])/2;
	      atom.satom[jn][1][1] -= (rp4*xk[nk][2]*xk[nk][2])/2;
	      atom.satom[jn][1][2] -= (rp4*xk[nk][2]*xk[nk][3])/2;
	      atom.satom[jn][2][2] -= (rp4*xk[nk][3]*xk[nk][3])/2;
	      atom.satom[kn][0][0] -= (rp4*xk[nk][1]*xk[nk][1])/2;
	      atom.satom[kn][0][1] -= (rp4*xk[nk][1]*xk[nk][2])/2;
	      atom.satom[kn][0][2] -= (rp4*xk[nk][1]*xk[nk][3])/2;
	      atom.satom[kn][1][1] -= (rp4*xk[nk][2]*xk[nk][2])/2;
	      atom.satom[kn][1][2] -= (rp4*xk[nk][2]*xk[nk][3])/2;
	      atom.satom[kn][2][2] -= (rp4*xk[nk][3]*xk[nk][3])/2;
	      rep=rp5*xl[nl][1];//
	      bre.rnp1[i]=bre.rnp1[i]+rep; bre.rnp1[ln]=bre.rnp1[ln]-rep;
	      rep=rp5*xl[nl][2];
	      bre.rnp2[i]=bre.rnp2[i]+rep; bre.rnp2[ln]=bre.rnp2[ln]-rep;
	      rep=rp5*xl[nl][3];
	      bre.rnp3[i]=bre.rnp3[i]+rep; bre.rnp3[ln]=bre.rnp3[ln]-rep;
	      // stress
	      atom.satom[i][0][0] -= (rp5*xl[nl][1]*xl[nl][1])/2;
	      atom.satom[i][0][1] -= (rp5*xl[nl][1]*xl[nl][2])/2;
	      atom.satom[i][0][2] -= (rp5*xl[nl][1]*xl[nl][3])/2;
	      atom.satom[i][1][1] -= (rp5*xl[nl][2]*xl[nl][2])/2;
	      atom.satom[i][1][2] -= (rp5*xl[nl][2]*xl[nl][3])/2;
	      atom.satom[i][2][2] -= (rp5*xl[nl][3]*xl[nl][3])/2;
	      atom.satom[ln][0][0] -= (rp5*xl[nl][1]*xl[nl][1])/2;
	      atom.satom[ln][0][1] -= (rp5*xl[nl][1]*xl[nl][2])/2;
	      atom.satom[ln][0][2] -= (rp5*xl[nl][1]*xl[nl][3])/2;
	      atom.satom[ln][1][1] -= (rp5*xl[nl][2]*xl[nl][2])/2;
	      atom.satom[ln][1][2] -= (rp5*xl[nl][2]*xl[nl][3])/2;
	      atom.satom[ln][2][2] -= (rp5*xl[nl][3]*xl[nl][3])/2;
	      //printf("HOGE1\n");
	    } // loop l (210)
	  } // loop k (220)
	} // if ((jbegin != jend)&&(lbegin != lend)) (230)
	btot=btot+btor*ator;
	dradi=dradi+datori*btor;
	dradj=dradj+datorj*btor;
	drdc=drdc+datorc*btor;
      } // if (kikj == ndihed) (231)
    OUT:
      // end dihedral forces
      
      bre.tote = bre.tote-btot*vatt;
      atom.epot[i] = atom.epot[i]-btot*vatt/2.0;
      atom.epot[jn] = atom.epot[jn]-btot*vatt/2.0;
      double vdbdi=vatt*dbdzi;
      double vdbdj=vatt*dbdzj;
      double vdrdc=vatt*drdc;
      double vdrdi=vatt*dradi;
      double vdrdj=vatt*dradj;
      double rp = vdbdi*xsij + vdbdj*xsji + btot*bre.dexx1[j];
      double rep; double repx,repy,repz;
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
      //printf("-HOGE2\n");
      // add many-body forces
      // i side of bond
      if (jbegin != jend) {
	nk = 0;
	for (int k=jbegin; k<=jend; k++) {
	  if (k == j) { continue; }
	  if (bre.lcheck[k] != 1) { continue; }
	  int kn = bre.list[k];
	  int kk = bre.ktype[kn];
	  double dwr = bre.dww[k]/bre.rcor[k];
	  nk++;
	  // first neighbors
	  double rp1 = vdbdi*(xsik[nk]+dwr*dexni[kk])
	    +dwr*(vdrdi+vdrdc*cfuni[nk])+vdbdi*dwr*sdalik;
	  double rp2=vdbdi*xsjk[nk];
	  double rep; double repx,repy,repz;
	  repx = rp1*bre.cor1[k];
	  bre.rnp1[i]=bre.rnp1[i]+repx; bre.rnp1[kn]=bre.rnp1[kn]-repx;
	  repy = rp1*bre.cor2[k];
	  bre.rnp2[i]=bre.rnp2[i]+repy; bre.rnp2[kn]=bre.rnp2[kn]-repy;
	  repz = rp1*bre.cor3[k];
	  bre.rnp3[i]=bre.rnp3[i]+repz; bre.rnp3[kn]=bre.rnp3[kn]-repz;
	  // stress
	  repx /= 2; repy /= 2; repz /= 2;
	  atom.satom[i][0][0] -= repx * bre.cor1[k];
	  atom.satom[i][0][1] -= repx * bre.cor2[k];
	  atom.satom[i][0][2] -= repx * bre.cor3[k];
	  atom.satom[i][1][1] -= repy * bre.cor2[k];
	  atom.satom[i][1][2] -= repy * bre.cor3[k];
	  atom.satom[i][2][2] -= repz * bre.cor3[k];
	  atom.satom[kn][0][0] -= repx * bre.cor1[k];
	  atom.satom[kn][0][1] -= repx * bre.cor2[k];
	  atom.satom[kn][0][2] -= repx * bre.cor3[k];
	  atom.satom[kn][1][1] -= repy * bre.cor2[k];
	  atom.satom[kn][1][2] -= repy * bre.cor3[k];
	  atom.satom[kn][2][2] -= repz * bre.cor3[k];
	  repx = rp2*xk[nk][1];
	  bre.rnp1[jn]=bre.rnp1[jn]+repx; bre.rnp1[kn]=bre.rnp1[kn]-repx;
	  repy = rp2*xk[nk][2];
	  bre.rnp2[jn]=bre.rnp2[jn]+repy; bre.rnp2[kn]=bre.rnp2[kn]-repy;
	  repz = rp2*xk[nk][3];
	  bre.rnp3[jn]=bre.rnp3[jn]+repz; bre.rnp3[kn]=bre.rnp3[kn]-repz;
	  // stress
	  repx /= 2; repy /= 2; repz /= 2;
	  atom.satom[jn][0][0] -= repx * xk[nk][1];
	  atom.satom[jn][0][1] -= repx * xk[nk][2];
	  atom.satom[jn][0][2] -= repx * xk[nk][3];
	  atom.satom[jn][1][1] -= repy * xk[nk][2];
	  atom.satom[jn][1][2] -= repy * xk[nk][3];
	  atom.satom[jn][2][2] -= repz * xk[nk][3];
	  atom.satom[kn][0][0] -= repx * xk[nk][1];
	  atom.satom[kn][0][1] -= repx * xk[nk][2];
	  atom.satom[kn][0][2] -= repx * xk[nk][3];
	  atom.satom[kn][1][1] -= repy * xk[nk][2];
	  atom.satom[kn][1][2] -= repy * xk[nk][3];
	  atom.satom[kn][2][2] -= repz * xk[nk][3];
	  //printf("*****HOGE3-1\n");
	  // second neighbors via radic
	  double ddr=vdrdc*dcfuni[nk]*2.0*conk;
	  if (ddr == 0.0) { continue; }
	  int mbegin = bre.nabors[kn];
	  int mend = bre.nabors[kn+1]-1;
	  if (mbegin >= mend ) { continue; }
	  for (int m=mbegin; m<=mend; m++) {
	    if (bre.lcheck[m] != 1) { continue; }
	    int mn = bre.list[m];
	    if (mn == kn) { continue; }
	    double rp = ddr*bre.dww[m]/bre.rcor[m];
	    repx = rp*bre.cor1[m];
	    bre.rnp1[kn]=bre.rnp1[kn]+repx; bre.rnp1[mn]=bre.rnp1[mn]-repx;
	    repy = rp*bre.cor2[m];
	    bre.rnp2[kn]=bre.rnp2[kn]+repy; bre.rnp2[mn]=bre.rnp2[mn]-repy;
	    repz = rp*bre.cor3[m];
	    bre.rnp3[kn]=bre.rnp3[kn]+repz; bre.rnp3[mn]=bre.rnp3[mn]-repz;
	    // stress
	    repx /= 2; repy /= 2; repz /= 2;
	    atom.satom[kn][0][0] -= repx * bre.cor1[m];
	    atom.satom[kn][0][1] -= repx * bre.cor2[m];
	    atom.satom[kn][0][2] -= repx * bre.cor3[m];
	    atom.satom[kn][1][1] -= repy * bre.cor2[m];
	    atom.satom[kn][1][2] -= repy * bre.cor3[m];
	    atom.satom[kn][2][2] -= repz * bre.cor3[m];
	    atom.satom[mn][0][0] -= repx * bre.cor1[m];
	    atom.satom[mn][0][1] -= repx * bre.cor2[m];
	    atom.satom[mn][0][2] -= repx * bre.cor3[m];
	    atom.satom[mn][1][1] -= repy * bre.cor2[m];
	    atom.satom[mn][1][2] -= repy * bre.cor3[m];
	    atom.satom[mn][2][2] -= repz * bre.cor3[m];
	    //printf("+++HOGE3\n");
	  } // loop m (17)
	} // loop k (22)
      } // if (jbegin != jend) (23)

      // j side of bond
      if (lbegin != lend) {
	int nl = 0;
	for (int l=lbegin; l<=lend; l++) {
	  int ln = bre.list[l];
	  if (ln == i) { continue; }
	  if (bre.lcheck[l] != 1) { continue; }
	  int kl = bre.ktype[ln];
	  double dwr = bre.dww[l]/bre.rcor[l];
	  nl++;
	  // first neighbors
	  double rp1 = vdbdj*(xsjl[nl]+dwr*dexnj[kl])
	    +dwr*(vdrdj+vdrdc*cfunj[nl])+vdbdj*dwr*sdaljl;
	  double rp2 = vdbdj*xsil[nl];
	  double rep; double repx,repy,repz;
	  repx = rp1*bre.cor1[l];
	  bre.rnp1[jn]=bre.rnp1[jn]+repx; bre.rnp1[ln]=bre.rnp1[ln]-repx;
	  repy = rp1*bre.cor2[l];
	  bre.rnp2[jn]=bre.rnp2[jn]+repy; bre.rnp2[ln]=bre.rnp2[ln]-repy;
	  repz = rp1*bre.cor3[l];
	  bre.rnp3[jn]=bre.rnp3[jn]+repz; bre.rnp3[ln]=bre.rnp3[ln]-repz;
	  // stress
	  repx /= 2; repy /= 2; repz /= 2;
	  atom.satom[jn][0][0] -= repx * bre.cor1[l];
	  atom.satom[jn][0][1] -= repx * bre.cor2[l];
	  atom.satom[jn][0][2] -= repx * bre.cor3[l];
	  atom.satom[jn][1][1] -= repy * bre.cor2[l];
	  atom.satom[jn][1][2] -= repy * bre.cor3[l];
	  atom.satom[jn][2][2] -= repz * bre.cor3[l];
	  atom.satom[ln][0][0] -= repx * bre.cor1[l];
	  atom.satom[ln][0][1] -= repx * bre.cor2[l];
	  atom.satom[ln][0][2] -= repx * bre.cor3[l];
	  atom.satom[ln][1][1] -= repy * bre.cor2[l];
	  atom.satom[ln][1][2] -= repy * bre.cor3[l];
	  atom.satom[ln][2][2] -= repz * bre.cor3[l];
	  repx = rp2*xl[nl][1];
	  bre.rnp1[i]=bre.rnp1[i]+repx; bre.rnp1[ln]=bre.rnp1[ln]-repx;
	  repy = rp2*xl[nl][2];
	  bre.rnp2[i]=bre.rnp2[i]+repy; bre.rnp2[ln]=bre.rnp2[ln]-repy;
	  repz = rp2*xl[nl][3];
	  bre.rnp3[i]=bre.rnp3[i]+repz; bre.rnp3[ln]=bre.rnp3[ln]-repz;
	  // stress
	  repx /= 2; repy /= 2; repz /= 2;
	  atom.satom[i][0][0] -= repx * xl[nl][1];
	  atom.satom[i][0][1] -= repx * xl[nl][2];
	  atom.satom[i][0][2] -= repx * xl[nl][3];
	  atom.satom[i][1][1] -= repy * xl[nl][2];
	  atom.satom[i][1][2] -= repy * xl[nl][3];
	  atom.satom[i][2][2] -= repz * xl[nl][3];
	  atom.satom[ln][0][0] -= repx * xl[nl][1];
	  atom.satom[ln][0][1] -= repx * xl[nl][2];
	  atom.satom[ln][0][2] -= repx * xl[nl][3];
	  atom.satom[ln][1][1] -= repy * xl[nl][2];
	  atom.satom[ln][1][2] -= repy * xl[nl][3];
	  atom.satom[ln][2][2] -= repz * xl[nl][3];
	  //printf("===---==HOGE4-1\n");
	  // second neighbors via radic
	  double ddr = vdrdc*dcfunj[nl]*2.0*conl;
	  if (ddr == 0.0) { continue; }
	  int nbegin = bre.nabors[ln];
	  int nend = bre.nabors[ln+1]-1;
	  if (nbegin >= nend) { continue; }
	  for (int n=nbegin; n<=nend; n++) {
	    if (bre.lcheck[n] != 1) { continue; }
	    int nn = bre.list[n];
	    if (nn == ln) { continue; }
	    double rp = ddr*bre.dww[n]/bre.rcor[n];
	    repx = rp*bre.cor1[n];
	    bre.rnp1[ln]=bre.rnp1[ln]+repx; bre.rnp1[nn]=bre.rnp1[nn]-repx;
	    repy = rp*bre.cor2[n];
	    bre.rnp2[ln]=bre.rnp2[ln]+repy; bre.rnp2[nn]=bre.rnp2[nn]-repy;
	    repz = rp*bre.cor3[n];
	    bre.rnp3[ln]=bre.rnp3[ln]+repz; bre.rnp3[nn]=bre.rnp3[nn]-repz;
	    // stress
	    repx /= 2; repy /= 2; repz /= 2;
	    atom.satom[ln][0][0] -= repx * bre.cor1[n];
	    atom.satom[ln][0][1] -= repx * bre.cor2[n];
	    atom.satom[ln][0][2] -= repx * bre.cor3[n];
	    atom.satom[ln][1][1] -= repy * bre.cor2[n];
	    atom.satom[ln][1][2] -= repy * bre.cor3[n];
	    atom.satom[ln][2][2] -= repz * bre.cor3[n];
	    atom.satom[nn][0][0] -= repx * bre.cor1[n];
	    atom.satom[nn][0][1] -= repx * bre.cor2[n];
	    atom.satom[nn][0][2] -= repx * bre.cor3[n];
	    atom.satom[nn][1][1] -= repy * bre.cor2[n];
	    atom.satom[nn][1][2] -= repy * bre.cor3[n];
	    atom.satom[nn][2][2] -= repz * bre.cor3[n];
	    //printf("=========HOGE4\n");
	  } // loop n (18)

	} // loop l (12)
      } // if (lbegin != lend) (30)

    } // loop j
  } // loop i
} // end of pibond()

void radic(int ki, int kj, double xnt1, double xnt2, double conjug,
	   double &rad, double &drdl, double &drdm, double &drdn)
{
  int l = int(xnt1);
  int m = int(xnt2);
  int n = int(conjug);
  rad = 0.0; drdl = 0.0; drdm = 0.0; drdn = 0.0;
  int kikj = ki+kj-1;
  if (l >= 4) { l = 4; xnt1 = 4.0; }
  if (m >= 4) { m = 4; xnt2 = 4.0; }
  if (n >= 9) { n = 9; conjug = 9.0; }
  for (int j=1; j<=64; j++) {
    double x = bre.clmn[kikj][l][m][n][j]*
      pow(xnt1,bre.in3[j][1])*pow(xnt2,bre.in3[j][2])*pow(conjug,bre.in3[j][3]);
    rad=rad+x;
    drdl=drdl+x*bre.in3[j][1]/xnt1;
    drdm=drdm+x*bre.in3[j][2]/xnt2;
    drdn=drdn+x*bre.in3[j][3]/conjug;
  }
}
void tor(double xnt1, double xnt2, double conjug, double &ator,
	 double &drdl, double &drdm, double &drdn)
{
  ator = 0.0; drdl = 0.0; drdm = 0.0; drdn = 0.0;
  if ((xnt1 >= 4.0)||(xnt2 >= 4.0)) { return; }
  int l = int(xnt1);
  int m = int(xnt2);
  int n = int(conjug);
  for (int j=1; j<=64; j++) {
    double x=bre.tlmn[l][m][n][j]*
      pow(xnt1,bre.in3[j][1])*pow(xnt2,bre.in3[j][2])*pow(conjug,bre.in3[j][3]);
    ator=ator+x;
    drdl=drdl+x*bre.in3[j][1]/xnt1;
    drdm=drdm+x*bre.in3[j][2]/xnt2;
    drdn=drdn+x*bre.in3[j][3]/conjug;
  }
}
void bcuint(int kl, int ki, double xx1, double xx2, int nh, int nc,
	    double &ansy, double &ansy1, double &ansy2)
{
  ansy = 0.0; ansy1 = 0.0; ansy2 = 0.0;
  if ((ki == 0)||(nh == 0)||(nc ==0)) {
    //    std::cout<<"bcuint error"<<std::endl; exit(0); }
    std::cout<<"bcuint error"<<std::endl; return; }
  for (int j=1; j<=16; j++) {
    double x=bre.clm[ki][nh][nc][j]*
      pow(xx1,bre.in2[j][1])*pow(xx2,bre.in2[j][2]);
    ansy=ansy+x;
    ansy1=ansy1+x*bre.in2[j][1]/xx1;
    ansy2=ansy2+x*bre.in2[j][2]/xx2;
  }
}
double anint(double x)
{
  double y;
  if (x==0.0) { return 0; }
  if (x>0.0) { y=x; } else { y=-x; }
  int i=int(y);
  if (y-(double)i>=0.5) { i++; }
  if (x<0.0) { i=-i; }
  return (double)i;
}
