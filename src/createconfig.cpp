#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#if defined(_WIN32) || defined(_WIN64)
#include <Windows.h>
#endif
#include <GL/gl.h>
#include <complex>
#include "myheader.h"

#define ZERO 1e-5

using namespace std;

extern GLuint objects;
//extern GLfloat *color[NOBJECTS];
extern GLfloat **color, **color0;
extern int *iatom, *repidx;
extern GLfloat red[], yellow[];
extern float ex;
extern char config_atom[], config_atom2[];
extern int milX1, milX2, milX3;
extern int milY1, milY2, milY3;
extern int milZ1, milZ2, milZ3;


void configmake_cub(int type, int mode);
void configmake_wurtzite(int mode);
void configmake_2dtri(int x, int y);
void configmake_graphene();
void cellmultiply(int ncell);
void configmake_cnt();
int ncount_cnt();
int lcm(int n1, int n2);
void set_atom_color();
void set_atom_weight();
void potential();
void bookkeep();
void pot_initialize_all();
void cnt_wall_discard();
int atom_number(char* at);
void deallocate_arrays();
void allocate_arrays();
bool ifzero(double x);
void inverse(double mat[][3], double imat[][3]);

void createconfig()
{
  //If config_atom starts with lower case, it's capitalized here.
  if ((int)config_atom[0]  > 0x60) { config_atom[0]  -=0x20; }
  if ((int)config_atom2[0] > 0x60) { config_atom2[0] -=0x20; }
  
  deallocate_arrays();

  //if (glIsList(objects+1)) {glDeleteLists(objects, atom.natom*3+1); printf("glDeleteLists\n"); }
  //if (glIsList(objects+1)) {glDeleteLists(objects, atom.natom*3+1); }
  if (glIsList(objects)) {glDeleteLists(objects, atom.natom*3); }

  if (config_type==0) { // fcc
    configmake_cub(0,0);
    atom.natom *= irepx*irepy*irepz;
  } else if (config_type==1) { // bcc
    configmake_cub(1,0);
    atom.natom *= irepx*irepy*irepz;
  } else if (config_type==2) { // dia
    configmake_cub(2,0);
    atom.natom *= irepx*irepy*irepz;
  } else if (config_type==3) { // Wurtzite
    configmake_wurtzite(0);
    atom.natom *= irepx*irepy*irepz;
  } else if (config_type==4) { // CNT
    atom.natom = ncount_cnt();
  } else if (config_type==5) { // 2-D triangular lattice
    atom.natom = 2*irepx*irepy;
  } else if (config_type==6) { // Graphene
    atom.natom = 4*irepx*irepy*irepz;
  }

  allocate_arrays();

  objects = glGenLists(atom.natom*3); //printf("objects = %d\n",objects);
  for (int i=1; i<=atom.natom*3; i++) { memcpy(color[i],yellow,sizeof(GLfloat)*4); }

  printf("======= CREATION OF NEW CONFIG ======\n");
  if (config_type==0) {
    printf("FCC %d x %d x %d x, N = %d\n",irepx,irepy,irepz,atom.natom);
    configmake_cub(0,1);
  } else if (config_type==1) {
    printf("BCC %d x %d x %d x, N = %d\n",irepx,irepy,irepz,atom.natom);
    configmake_cub(1,1);
  } else if (config_type==2) {
    printf("Diamond %d x %d x %d x, N = %d\n",irepx,irepy,irepz,atom.natom);
    configmake_cub(2,1);
  } else if (config_type==3) {
    printf("Wurtzite %d x %d x %d x, N = %d\n",irepx,irepy,irepz,atom.natom);
    configmake_wurtzite(1);
  } else if (config_type==4) {
    printf("Nanotube (%d, %d) x %d, N = %d\n",icntm,icntn,irepz,atom.natom);
    configmake_cnt(); cnt_wall_discard();
  } else if (config_type==5) {
    printf("2-D triangle lattice %d x %d, N = %d\n",irepx,irepy,atom.natom);
    irepz = 1;
    configmake_2dtri(irepx, irepy);
  } else if (config_type==6) {
    printf("Graphene-type honeycomb %d x %d x %d, N = %d\n",irepx,irepy,irepz,atom.natom);
    configmake_graphene();
  }

  if (strcmp(config_atom,"M")==0) {
    printf("Multi atom system\n");
  } else {
    printf("Atom type: %s\n",config_atom);
    //for (int i=1; i<=atom.natom; i++) {
    //  strcpy(atom.asp[i], config_atom);
    //}
  }
  printf("Cell dimension (ang):      x          y          z\n");
  printf("                  1 : %8.3f   %8.3f   %8.3f\n",
	 cell.hmat[0][0]/ang,cell.hmat[1][0]/ang,cell.hmat[2][0]/ang);
  printf("                  2 : %8.3f   %8.3f   %8.3f\n",
	 cell.hmat[0][1]/ang,cell.hmat[1][1]/ang,cell.hmat[2][1]/ang);
  printf("                  3 : %8.3f   %8.3f   %8.3f\n",
	 cell.hmat[0][2]/ang,cell.hmat[1][2]/ang,cell.hmat[2][2]/ang);
  printf("=====================================\n");
  cellx=cell.hmat[0][0]/ang;celly=cell.hmat[1][1]/ang;cellz=cell.hmat[2][2]/ang;
  cellx=sqrt(cell.hmat[0][0]*cell.hmat[0][0]+cell.hmat[1][0]*cell.hmat[1][0]+cell.hmat[2][0]*cell.hmat[2][0])/ang;
  celly=sqrt(cell.hmat[0][1]*cell.hmat[0][1]+cell.hmat[1][1]*cell.hmat[1][1]+cell.hmat[2][1]*cell.hmat[2][1])/ang;
  cellz=sqrt(cell.hmat[0][2]*cell.hmat[0][2]+cell.hmat[1][2]*cell.hmat[1][2]+cell.hmat[2][2]*cell.hmat[2][2])/ang;
  cell1x=cell.hmat[0][0]/ang;cell1y=cell.hmat[1][0]/ang;cell1z=cell.hmat[2][0]/ang;
  cell2x=cell.hmat[0][1]/ang;cell2y=cell.hmat[1][1]/ang;cell2z=cell.hmat[2][1]/ang;
  cell3x=cell.hmat[0][2]/ang;cell3y=cell.hmat[1][2]/ang;cell3z=cell.hmat[2][2]/ang;

  // Atom color
  set_atom_color();
  // Atom weight
  set_atom_weight();
  // Atom number
  for (int i=1; i<=atom.natom; i++) {
    atom.anum[i] = atom_number(atom.asp[i]);
  }
  for (int i=1; i<=atom.natom; i++) {
    atom.rx_org[i]=atom.rx[i]; atom.ry_org[i]=atom.ry[i]; atom.rz_org[i]=atom.rz[i];
    atom.vx[i] = 0.0; atom.vy[i] = 0.0; atom.vz[i] = 0.0;
    atom.ax[i] = 0.0; atom.ay[i] = 0.0; atom.az[i] = 0.0;
    atom.bx[i] = 0.0; atom.by[i] = 0.0; atom.bz[i] = 0.0;
    atom.cx[i] = 0.0; atom.cy[i] = 0.0; atom.cz[i] = 0.0;
  }
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      cell.hmat_org[i][j]=cell.hmat[i][j];
      cell.hvmat[i][j] = 0; cell.hamat[i][j] = 0; cell.hbmat[i][j] = 0;
      cell.hcmat[i][j] = 0; cell.sgmmat_set[i][j] = 0;
    }
  }
  // Let all atoms move
  for (int i=1; i<=atom.natom; i++) {
    atom.mfx[i]=false;atom.mfy[i]=false;atom.mfz[i]=false;
  }
  istep = 0;
  ex = 0.0;
  pot_initialize_all();
  bookkeep();
  potential();
}

void configmake_cub(int type, int mode)
{
  //type: 0 = fcc, 1 = bcc, 2 = diamond
  //mode: 0 = count, 1 = actual set
  double ao=alat;//*1.0e-10;
  cell.alat = (double)alat;
  double rrx[9], rry[9], rrz[9];
  int atom_type[9];
  int ncell;
  if (type == 0) {
    ncell = 4;
    rrx[1]=   0.0;rry[1]=   0.0;rrz[1]=   0.0;
    rrx[2]=   0.0;rry[2]=ao/2.0;rrz[2]=ao/2.0;
    rrx[3]=ao/2.0;rry[3]=   0.0;rrz[3]=ao/2.0;
    rrx[4]=ao/2.0;rry[4]=ao/2.0;rrz[4]=   0.0;
    for (int i=1; i<=ncell; i++) { atom_type[1] = 1; }
  } else if (type == 1) {
    ncell = 2;
    rrx[1]=   0.0;rry[1]=   0.0;rrz[1]=   0.0;
    rrx[2]=ao/2.0;rry[2]=ao/2.0;rrz[2]=ao/2.0;
    for (int i=1; i<=ncell; i++) { atom_type[i] = 1; }
  } else if (type == 2) {
    ncell = 8;
    rrx[1]=   0.0;rry[1]=   0.0;rrz[1]=   0.0;
    rrx[2]=   0.0;rry[2]=ao/2.0;rrz[2]=ao/2.0;
    rrx[3]=ao/2.0;rry[3]=   0.0;rrz[3]=ao/2.0;
    rrx[4]=ao/2.0;rry[4]=ao/2.0;rrz[4]=   0.0;
    for (int i=5; i<=8; i++) {
      rrx[i]=rrx[i-4]+ao/4.0;
      rry[i]=rry[i-4]+ao/4.0;
      rrz[i]=rrz[i-4]+ao/4.0;
    }
    for (int i=1; i<=4; i++) { atom_type[i] = 1; }
    for (int i=5; i<=8; i++) { atom_type[i] = 2; }
  }

  double millerX[3], millerY[3], millerZ[3];
  //millerX[0]= 1; millerX[1]= 0; millerX[2]= 2;
  //millerY[0]= 0; millerY[1]= 1; millerY[2]= 0;
  //millerZ[0]= 2; millerZ[1]= 0; millerZ[2]=-1;
  if (milX1 != 0) { millerX[0] = (double)milX1; } else { millerX[0] = 0.; }
  if (milX2 != 0) { millerX[1] = (double)milX2; } else { millerX[1] = 0.; }
  if (milX3 != 0) { millerX[2] = (double)milX3; } else { millerX[2] = 0.; }
  if (milY1 != 0) { millerY[0] = (double)milY1; } else { millerY[0] = 0.; }
  if (milY2 != 0) { millerY[1] = (double)milY2; } else { millerY[1] = 0.; }
  if (milY3 != 0) { millerY[2] = (double)milY3; } else { millerY[2] = 0.; }
  if (milZ1 != 0) { millerZ[0] = (double)milZ1; } else { millerZ[0] = 0.; }
  if (milZ2 != 0) { millerZ[1] = (double)milZ2; } else { millerZ[1] = 0.; }
  if (milZ3 != 0) { millerZ[2] = (double)milZ3; } else { millerZ[2] = 0.; }

  /*
  if (milX1*milY1+milX2*milY2+milX3*milY3 != 0) {
    printf("### WARNING: X and Y are not orthogonal\n"); }
  if (milZ1*milY1+milZ2*milY2+milZ3*milY3 != 0) {
    printf("### WARNING: Y and Z are not orthogonal\n"); }
  if (milX1*milZ1+milX2*milZ2+milX3*milZ3 != 0) {
    printf("### WARNING: Z and X are not orthogonal\n"); }
  */

  double xx, yy, zz, xxx, yyy, zzz;
  int nx, ny, nz;
  for (int i=1; i<=100; i++) {
    xx = (double)millerX[0]*(double)i ; xxx = xx-floor(xx);
    if (ifzero(xxx)) {
      yy = (double)millerX[1]*(double)i ; yyy = yy-floor(yy);
      if (ifzero(yyy)) {
	zz = (double)millerX[2]*(double)i ; zzz = zz-floor(zz);
	if (ifzero(zzz)) {
	  cell.hmat[0][0] = (double)millerX[0]*(double)i*ao;
	  cell.hmat[1][0] = (double)millerX[1]*(double)i*ao;
	  cell.hmat[2][0] = (double)millerX[2]*(double)i*ao;
	  nx = i;
	  goto OUT1;
	}
      }
    }
  }
 OUT1:
  for (int i=1; i<=100; i++) {
    xx = (double)millerY[0]*(double)i ; xxx = xx-floor(xx);
    if (ifzero(xxx)) {
      yy = (double)millerY[1]*(double)i ; yyy = yy-floor(yy);
      if (ifzero(yyy)) {
	zz = (double)millerY[2]*(double)i ; zzz = zz-floor(zz);
	if (ifzero(zzz)) {
	  cell.hmat[0][1] = (double)millerY[0]*(double)i*ao;
	  cell.hmat[1][1] = (double)millerY[1]*(double)i*ao;
	  cell.hmat[2][1] = (double)millerY[2]*(double)i*ao;
	  ny = i;
	  goto OUT2;
	}
      }
    }
  }
 OUT2:
  for (int i=1; i<=100; i++) {
    xx = (double)millerZ[0]*(double)i ; xxx = xx-floor(xx);
    if (ifzero(xxx)) {
      yy = (double)millerZ[1]*(double)i ; yyy = yy-floor(yy);
      if (ifzero(yyy)) {
	zz = (double)millerZ[2]*(double)i ; zzz = zz-floor(zz);
	if (ifzero(zzz)) {
	  cell.hmat[0][2] = (double)millerZ[0]*(double)i*ao;
	  cell.hmat[1][2] = (double)millerZ[1]*(double)i*ao;
	  cell.hmat[2][2] = (double)millerZ[2]*(double)i*ao;
	  nz = i;
	  goto OUT3;
	}
      }
    }
  }
 OUT3:
  if (cell.Getvolume() < 0) { // check if right-hand system
    cell.hmat[0][2] *= -1.0; cell.hmat[1][2] *= -1.0; cell.hmat[2][2] *= -1.0;
  }
  inverse(cell.hmat,cell.hinmat);
  int icnt = 0;
  for (int ix=-10; ix<=10; ix++) {
    for (int iy=-10; iy<=10; iy++) {
      for (int iz=-10; iz<=10; iz++) {
	for (int i=1; i<=ncell; i++) {
	  xx = rrx[i]+ao*(double)ix;
	  yy = rry[i]+ao*(double)iy;
	  zz = rrz[i]+ao*(double)iz;
	  xxx = cell.hinmat[0][0]*xx + cell.hinmat[0][1]*yy + cell.hinmat[0][2]*zz;
	  yyy = cell.hinmat[1][0]*xx + cell.hinmat[1][1]*yy + cell.hinmat[1][2]*zz;
	  zzz = cell.hinmat[2][0]*xx + cell.hinmat[2][1]*yy + cell.hinmat[2][2]*zz;
	  if ((xxx>-ZERO)&&(xxx<1.0-ZERO)&&(yyy>-ZERO)&&(yyy<1.0-ZERO)&&(zzz>-ZERO)&&(zzz<1.0-ZERO)) {
	    icnt++;
	    if (mode == 1) {
	      atom.rx[icnt] = xx*ang; atom.ry[icnt] = yy*ang; atom.rz[icnt] = zz*ang;
	      if (atom_type[i] == 1) {
		strcpy(atom.asp[icnt], config_atom);
	      } else {
		strcpy(atom.asp[icnt], config_atom2);
	      }
	    }
	  }
	}
      }
    }
  }
  if (mode == 0) {
    atom.natom = icnt;
  }
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      cell.hmat[i][j] *= ang;
    } }
  inverse(cell.hmat,cell.hinmat);

  cell.Setlen();
  xx = cell.len[0];  yy = cell.len[1];  zz = cell.len[2];
  //cell.hmat[0][0]=xx;cell.hmat[1][0]=0.;cell.hmat[2][0]=0.;
  //cell.hmat[0][1]=0.;cell.hmat[1][1]=yy;cell.hmat[2][1]=0.;
  //cell.hmat[0][2]=0.;cell.hmat[1][2]=0.;cell.hmat[2][2]=zz;

  if (mode == 1) {
    for (int i=1; i<=icnt; i++) {
      xx=atom.rx[i]; yy=atom.ry[i]; zz=atom.rz[i];
      xxx = cell.hinmat[0][0]*xx + cell.hinmat[0][1]*yy + cell.hinmat[0][2]*zz;
      yyy = cell.hinmat[1][0]*xx + cell.hinmat[1][1]*yy + cell.hinmat[1][2]*zz;
      zzz = cell.hinmat[2][0]*xx + cell.hinmat[2][1]*yy + cell.hinmat[2][2]*zz;
      atom.rx[i] = cell.hmat[0][0]*xxx + cell.hmat[0][1]*yyy + cell.hmat[0][2]*zzz;
      atom.ry[i] = cell.hmat[1][0]*xxx + cell.hmat[1][1]*yyy + cell.hmat[1][2]*zzz;
      atom.rz[i] = cell.hmat[2][0]*xxx + cell.hmat[2][1]*yyy + cell.hmat[2][2]*zzz;
    }
  }

  inverse(cell.hmat,cell.hinmat);
  if (mode == 1) {
    cellmultiply(icnt);
  }
}

void configmake_wurtzite(int mode)
{
  //mode: 0 = count, 1 = actual set
  double ao=alat;//*1.0e-10;
  double co=clat;
  cell.alat = (double)alat;
  cell.clat = (double)clat;
  double rrx[5], rry[5], rrz[5];
  double xx, yy, zz;
  int atom_type[5];
  int ncell;
  ncell = 4;
  if (mode == 0) {
    atom.natom = ncell;
    return;
  }
  cell.hmat[0][0] =  ao;
  cell.hmat[1][0] = 0.0;
  cell.hmat[2][0] = 0.0;
  cell.hmat[0][1] = ao*0.5;
  cell.hmat[1][1] = ao*sqrt(3.0)/2.0;
  cell.hmat[2][1] = 0.0;
  cell.hmat[0][2] = 0.0;
  cell.hmat[1][2] = 0.0;
  cell.hmat[2][2] =  co;
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      cell.hmat[i][j] *= ang;
    } }
  rrx[1]=   0.0; rry[1]=             0.0; rrz[1]=   0.000;
  rrx[2]=ao*0.5; rry[2]=ao/2.0/sqrt(3.0); rrz[2]=co*0.500;
  rrx[3]=ao*0.5; rry[3]=ao/2.0/sqrt(3.0); rrz[3]=co*0.125;
  rrx[4]=   0.0; rry[4]=             0.0; rrz[4]=co*0.625;
  for (int i=1; i<=2; i++) { atom_type[i] = 1; }
  for (int i=3; i<=4; i++) { atom_type[i] = 2; }

  for (int i=1; i<=ncell; i++) {
    xx = rrx[i]; yy = rry[i]; zz = rrz[i];
    atom.rx[i] = xx*ang; atom.ry[i] = yy*ang; atom.rz[i] = zz*ang;
    if (atom_type[i] == 1) {
      strcpy(atom.asp[i], config_atom);
    } else {
      strcpy(atom.asp[i], config_atom2);
    }
  }

  inverse(cell.hmat,cell.hinmat);
  if (mode == 1) {
    cellmultiply(ncell);
  }
}

void configmake_2dtri(int x, int y)
{
  double ao=alat*1.0e-10;
  cell.alat = (double)alat;
  atom.rx[1]=   0.0;atom.ry[1]=             0.0;atom.rz[1]=ao*2.5;
  atom.rx[2]=ao/2.0;atom.ry[2]=ao/2.0*sqrt(3.0);atom.rz[2]=ao*2.5;
  strcpy(atom.asp[1],config_atom);
  strcpy(atom.asp[2],config_atom);
  int ncell=2;
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) { cell.hmat[i][j]=0.0; }  }
  cell.hmat[0][0] = ao;
  cell.hmat[1][1] = ao*sqrt(3.0);
  cell.hmat[2][2] = ao*5.0;
  cellmultiply(ncell);
  for (int i=1; i<=atom.natom; i++) {
    atom.mfx[i]=false;atom.mfy[i]=false;atom.mfz[i]=true;
  }
}

void configmake_graphene()
{
  double ao=alat*1.0e-10;
  cell.alat = (double)alat;
  double az=ao*1.36;
  atom.rx[1]=             0.0;atom.ry[1]=    0.0;atom.rz[1]=0.5*az;
  atom.rx[2]=    ao/sqrt(3.0);atom.ry[2]=    0.0;atom.rz[2]=0.5*az;
  atom.rx[3]=  sqrt(3.0)/2*ao;atom.ry[3]= 0.5*ao;atom.rz[3]=0.5*az;
  atom.rx[4]=sqrt(3.0)*5/6*ao;atom.ry[4]= 0.5*ao;atom.rz[4]=0.5*az;
  strcpy(atom.asp[1],config_atom);
  strcpy(atom.asp[2],config_atom2);
  strcpy(atom.asp[3],config_atom);
  strcpy(atom.asp[4],config_atom2);
  int ncell=4;
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) { cell.hmat[i][j]=0.0; }  }
  cell.hmat[0][0] = ao*sqrt(3.0);
  cell.hmat[1][1] = ao;
  cell.hmat[2][2] = az;
  cellmultiply(ncell);
  inverse(cell.hmat,cell.hinmat);
}

void cellmultiply(int ncell)
{
  int iii=ncell;
  int jjj=iii;
  if (irepx > 1) {
    for (int k=2; k<=irepx; k++) {
      for (int j=1; j<=jjj; j++) {
	iii++;
	atom.rx[iii]=atom.rx[iii-jjj]+cell.hmat[0][0];
	atom.ry[iii]=atom.ry[iii-jjj]+cell.hmat[1][0];
	atom.rz[iii]=atom.rz[iii-jjj]+cell.hmat[2][0];
	strcpy(atom.asp[iii],atom.asp[iii-jjj]);
      }
    }
    jjj=iii;
  }
  if (irepy > 1) {
    for (int k=2; k<=irepy; k++) {
      for (int j=1; j<=jjj; j++) {
	iii++;
	atom.rx[iii]=atom.rx[iii-jjj]+cell.hmat[0][1];
	atom.ry[iii]=atom.ry[iii-jjj]+cell.hmat[1][1];
	atom.rz[iii]=atom.rz[iii-jjj]+cell.hmat[2][1];
	strcpy(atom.asp[iii],atom.asp[iii-jjj]);
      }
    }
    jjj=iii;
  }
  if (irepz > 1) {
    for (int k=2; k<=irepz; k++) {
      for (int j=1; j<=jjj; j++) {
	iii++;
	atom.rx[iii]=atom.rx[iii-jjj]+cell.hmat[0][2];
	atom.ry[iii]=atom.ry[iii-jjj]+cell.hmat[1][2];
	atom.rz[iii]=atom.rz[iii-jjj]+cell.hmat[2][2];
	strcpy(atom.asp[iii],atom.asp[iii-jjj]);
      }
    }
  }
  for (int i=0; i<3; i++) {
    cell.hmat[i][0]=cell.hmat[i][0]*(double)irepx;
    cell.hmat[i][1]=cell.hmat[i][1]*(double)irepy;
    cell.hmat[i][2]=cell.hmat[i][2]*(double)irepz;
  }
}

void configmake_cnt()
{
  int na=4;
  double x[5],y[5],z[5];
  double a1[4],a2[4],a3[4];
  int atom_type[5];
  complex<double> cr, crcnj, ctube, ctubecnj, ctrans;

  int natom = 1;
  int iatom = 0;
  double ao = 1.42e-10;
  int nwall = 0;

  int iwall = 1;
  int n1 = icntm;
  int n2 = icntn;
  int nlayer = irepz;

  double pi = 4.0*atan(1.0);

  if (n1+n2 < 2) { return; }

  a1[1]=ao*sqrt(3.0);
  a1[2]=0.0;
  a1[3]=0.0;
  a2[1]=0.0;
  a2[2]=ao*3.0;
  a2[3]=0.0;
  a3[1]=0.0;
  a3[2]=0.0;
  a3[3]=6.696/1.42*ao;
  
  x[1]=0.0;
  y[1]=0.0;
  z[1]=0.0;
  x[2]=0.0;
  y[2]=ao;
  z[2]=0.0;
  x[3]=ao*sqrt(3.0)/2.0;
  y[3]=-ao/2.0;
  z[3]=0.0;
  x[4]=ao*sqrt(3.0)/2.0;
  y[4]=ao*3.0/2.0;
  z[4]=0.0;
  atom_type[1]=1;
  atom_type[2]=2;
  atom_type[3]=2;
  atom_type[4]=1;

  double at11=a1[1];
  double at12=0.0;
  double at21=a1[1]/2.0;
  double at22=-ao*3.0/2.0;

  double atube1=at11*(double)n1+at21*(double)n2;
  double atube2=at12*(double)n1+at22*(double)n2;
  double rtubex=sqrt(atube1*atube1+atube2*atube2);
  
  double tubecos=atube1/rtubex;
  double tubesin=atube2/rtubex;
  ctube=complex<double>(tubecos,tubesin);
  ctubecnj=complex<double>(tubecos,-tubesin);
  for (int i=1; i<=na; i++) {
    cr=complex<double>(x[i],y[i])*ctubecnj;
    x[i]=cr.real();
    y[i]=cr.imag();
  }

  int nn1=2*n1+n2;
  int nn2=2*n2+n1;
  //  call lcm(nn1,nn2,nc)
  int nc = lcm(nn1,nn2);
  
  double as1=at11*(double)(nc/nn1)-at21*(double)(nc/nn2);
  double as2=at12*(double)(nc/nn1)-at22*(double)(nc/nn2);
  double rtubey=sqrt(as1*as1+as2*as2);
  //  printf("%e %e\n",as1,as2);

  ctrans=complex<double>(a1[1],a1[2]);
  ctrans=ctrans*ctubecnj;
  a1[1]=ctrans.real();
  a1[2]=ctrans.imag();
  ctrans=complex<double>(a2[1],a2[2]);
  ctrans=ctrans/ctube;
  a2[1]=ctrans.real();
  a2[2]=ctrans.imag();

  for (int i=1; i<=na; i++) {
    for (int ix=-200; ix<=200; ix++) {
      for (int iy=-200; iy<=200; iy++) {
	double xx=x[i]+(double)ix*a1[1]+(double)iy*a2[1];
	double yy=y[i]+(double)ix*a1[2]+(double)iy*a2[2];
	if ((xx >= -1.0e-12)&&(xx < rtubex-1.0e-12)) {
	  if ((yy >= -1.0e-12)&&(yy < rtubey-1.0e-12)) {
	    iatom++;
	    if (iatom > atom.natom) { printf("Warning in configmake_cnt\n"); }
	    atom.rx[iatom]=xx;
	    atom.ry[iatom]=yy;
	    atom.rz[iatom]=0.0;
	    if (atom_type[i] == 1) {
	      strcpy(atom.asp[iatom],config_atom);
	    } else {
	      strcpy(atom.asp[iatom],config_atom2);
	    }
	  }
	}
      }
    }
  }

  for (int i=natom; i<=iatom; i++) {
    atom.rz[i]=atom.ry[i];
    double xx=atom.rx[i];
    atom.rx[i]=rtubex/(2.0*pi)*cos(2.0*pi*xx/rtubex+2.0*pi/96.0);
    atom.ry[i]=rtubex/(2.0*pi)*sin(2.0*pi*xx/rtubex+2.0*pi/96.0);
  }
  natom=iatom+1;

  int katom=iatom;
  for (int i=1; i<=nlayer-1; i++) {
    for (int j=1; j<=katom; j++) {
      iatom++;
      if (iatom > atom.natom) { printf("Warning in configmake_cnt (2)\n"); }
      atom.rx[iatom]=atom.rx[j];
      atom.ry[iatom]=atom.ry[j];
      atom.rz[iatom]=atom.rz[j]+rtubey*(double)i;
      strcpy(atom.asp[iatom],atom.asp[j]);
    }
  }

  // Rotation
  double rad = sqrt(atom.rx[1]*atom.rx[1]+atom.ry[1]*atom.ry[1]);
  for (int i=1; i<=iatom; i++) {
    double theta = 180.0/pi*atan(atom.ry[i]/atom.rx[i]);
    if ((atom.ry[i] > 0.0)&&(theta < 0.0)) { theta = theta + 180.0; }
    if ((atom.ry[i] < 0.0)&&(theta > 0.0)) { theta = theta + 180.0; }
    if ((atom.ry[i] < 0.0)&&(theta < 0.0)) { theta = theta + 360.0; }
    if ((atom.ry[i] == 0.0)&&(atom.rx[i]>0.0)) { theta = 0; }
    if ((atom.ry[i] == 0.0)&&(atom.rx[i]<0.0)) { theta = 180; }
    if ((atom.rx[i] == 0.0)&&(atom.ry[i]>0.0)) { theta = 90; }
    if ((atom.rx[i] == 0.0)&&(atom.ry[i]<0.0)) { theta = 270; }
    atom.rx[i]=rad*cos((theta+rotz)/180.0*pi);
    atom.ry[i]=rad*sin((theta+rotz)/180.0*pi);
  }
  for (int i=0; i<3; i++) { for (int j=0; j<3; j++) { cell.hmat[i][j]=0.0; } }
  cell.hmat[0][0]=cscnt*1e-10;
  cell.hmat[1][1]=cscnt*1e-10;
  cell.hmat[2][2]=rtubey*(double)nlayer;
  // Shift
  for (int i=1; i<=iatom; i++) {
    atom.rx[i]=atom.rx[i]+cell.hmat[0][0]/2.0;
    atom.ry[i]=atom.ry[i]+cell.hmat[1][1]/2.0;
    atom.rz[i]=atom.rz[i]+cell.hmat[2][2]*shiftz/100.0;
    if (atom.rz[i] > cell.hmat[2][2]) { atom.rz[i] = atom.rz[i] - cell.hmat[2][2]; }
  }
  printf("CNT diameter (ang) = %f\n",rad*2/ang);
}

int lcm(int n1, int n2)
{
  int m1, m2, mm1, nc;
  int dif = n1-n2;
  if (dif < 0) {
    m1=n1;
    m2=n2;
  } else {
    m1=n2;
    m2=n1;
  }
  mm1=m1;
  nc=1;

  int ifg=1;
  while(ifg == 1)
    {
      for (int i=2; i<=mm1; i++) {
	if ((m1%i == 0)&&(m2%i == 0)) {
	  m1=m1/i;
	  m2=m2/i;
	  nc=nc*i;
	  ifg = 1;
	  break;
	} else {
	  ifg = 0;
	}
      }
    }
  nc=nc*m1*m2;
  return nc;
}

int ncount_cnt()
{
  int na=4;
  double x[5],y[5],z[5];
  double a1[4],a2[4],a3[4];
  complex<double> cr, crcnj, ctube, ctubecnj, ctrans;

  int natom = 1;
  int iatom = 0;
  double ao = 1.42e-10;
  int nwall = 0;

  int iwall = 1;
  int n1 = icntm;
  int n2 = icntn;
  int nlayer = irepz;

  double pi = 4.0*atan(1.0);

  if (n1+n2 < 2) { return 1; }

  a1[1]=ao*sqrt(3.0);
  a1[2]=0.0;
  a1[3]=0.0;
  a2[1]=0.0;
  a2[2]=ao*3.0;
  a2[3]=0.0;
  a3[1]=0.0;
  a3[2]=0.0;
  a3[3]=6.696/1.42*ao;
  
  x[1]=0.0;
  y[1]=0.0;
  z[1]=0.0;
  x[2]=0.0;
  y[2]=ao;
  z[2]=0.0;
  x[3]=ao*sqrt(3.0)/2.0;
  y[3]=-ao/2.0;
  z[3]=0.0;
  x[4]=ao*sqrt(3.0)/2.0;
  y[4]=ao*3.0/2.0;
  z[4]=0.0;

  double at11=a1[1];
  double at12=0.0;
  double at21=a1[1]/2.0;
  double at22=-ao*3.0/2.0;

  double atube1=at11*(double)n1+at21*(double)n2;
  double atube2=at12*(double)n1+at22*(double)n2;
  double rtubex=sqrt(atube1*atube1+atube2*atube2);
  
  double tubecos=atube1/rtubex;
  double tubesin=atube2/rtubex;
  ctube=complex<double>(tubecos,tubesin);
  ctubecnj=complex<double>(tubecos,-tubesin);
  for (int i=1; i<=na; i++) {
    cr=complex<double>(x[i],y[i])*ctubecnj;
    x[i]=cr.real();
    y[i]=cr.imag();
  }

  int nn1=2*n1+n2;
  int nn2=2*n2+n1;
  int nc = lcm(nn1,nn2);
  
  double as1=at11*(double)(nc/nn1)-at21*(double)(nc/nn2);
  double as2=at12*(double)(nc/nn1)-at22*(double)(nc/nn2);
  double rtubey=sqrt(as1*as1+as2*as2);

  ctrans=complex<double>(a1[1],a1[2]);
  ctrans=ctrans*ctubecnj;
  a1[1]=ctrans.real();
  a1[2]=ctrans.imag();
  ctrans=complex<double>(a2[1],a2[2]);
  ctrans=ctrans/ctube;
  a2[1]=ctrans.real();
  a2[2]=ctrans.imag();

  for (int i=1; i<=na; i++) {
    for (int ix=-200; ix<=200; ix++) {
      for (int iy=-200; iy<=200; iy++) {
	double xx=x[i]+(double)ix*a1[1]+(double)iy*a2[1];
	double yy=y[i]+(double)ix*a1[2]+(double)iy*a2[2];
	if ((xx >= -1.0e-12)&&(xx < rtubex-1.0e-12)) {
	  if ((yy >= -1.0e-12)&&(yy < rtubey-1.0e-12)) {
	    iatom++;
	  }
	}
      }
    }
  }
  natom=iatom+1;

  int katom=iatom;
  for (int i=1; i<=nlayer-1; i++) {
    for (int j=1; j<=katom; j++) {
      iatom++;
    }
  }
  return iatom;
}
 
bool ifzero(double x)
{
  bool b;
  if ((x < ZERO)&&(x > -ZERO)) {
    b = true;
  } else {
    b = false;
  }
  return b;
}

