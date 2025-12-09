#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#if defined(_WIN32) || defined(_WIN64)
#include <Windows.h>
#endif
#include <GL/gl.h>
#include "myheader.h"

void set_atom_weight();
void potential();
void bookkeep();
void pot_initialize_all();
void set_atom_color();
int atom_number(char* at);
void deallocate_arrays();
void allocate_arrays();
void multiply_cell(int ix, int iy, int iz, int mode);

extern GLuint objects;
extern GLfloat **color;
extern GLfloat yellow[];
extern float ex;

void multiply_cell(int ix, int iy, int iz)
{
  multiply_cell(ix, iy, iz, 0);
}

void multiply_cell(int ix, int iy, int iz, int mode)
{
  // mode = 0: normal (1 ~ ix) 1: double (-ix ~ +ix)
  //           n->n*ix*iy*iz      n->n*(2*ix-1)*(2*iy-1)*(2*iz-1)
  printf("Cell multiplication: %d x %d x %d\n",ix,iy,iz);
  if (ix < 1) { ix = 1; }
  if (iy < 1) { iy = 1; }
  if (iz < 1) { iz = 1; }
  // arrays for storing data
  double *rxs, *rys, *rzs;
  char **asps; int *anums;
  bool *mfxs, *mfys, *mfzs;
  rxs = new double[atom.natom+1];
  rys = new double[atom.natom+1];
  rzs = new double[atom.natom+1];
  asps = new char*[atom.natom+1];
  for (int i=0; i<=atom.natom; i++) { asps[i] = new char[3]; }
  anums = new int[atom.natom+1];
  mfxs = new bool[atom.natom+1];
  mfys = new bool[atom.natom+1];
  mfzs = new bool[atom.natom+1];
  // Store coordinates and others
  for (int i=1; i<=atom.natom; i++) {
    rxs[i] = atom.rx[i]; rys[i] = atom.ry[i]; rzs[i] = atom.rz[i];
    strcpy(asps[i], atom.asp[i]);
    anums[i] = atom.anum[i];
    mfxs[i] = atom.mfx[i]; mfys[i] = atom.mfy[i]; mfzs[i] = atom.mfz[i];
  }
  if (glIsList(objects+1)) {glDeleteLists(objects, atom.natom*3+1); } //printf("glDeleteLists\n"); }
  deallocate_arrays();
  int natom_org = atom.natom;
  if (mode == 0) {
    atom.natom *= ix*iy*iz;
  } else {
    atom.natom *= (ix*2-1)*(iy*2-1)*(iz*2-1);
  }
  allocate_arrays();
  // copy original coordinates and others
  for (int i=1; i<=natom_org; i++) {
    atom.rx[i] = rxs[i]; atom.ry[i] = rys[i]; atom.rz[i] = rzs[i];
    strcpy(atom.asp[i], asps[i]);
    atom.anum[i] = anums[i];
    atom.mfx[i] = mfxs[i]; atom.mfy[i] = mfys[i]; atom.mfz[i] = mfzs[i];
  }
  // multiply atoms
  int natom_base = natom_org;
  if (ix > 1) {
    if (mode == 0) {
      for (int k=1; k<ix; k++) {
	for (int i=1; i<=natom_base; i++) {
	  atom.rx[i+natom_base*k] = atom.rx[i] + cell.hmat[0][0]*(double)k;
	  atom.ry[i+natom_base*k] = atom.ry[i] + cell.hmat[1][0]*(double)k;
	  atom.rz[i+natom_base*k] = atom.rz[i] + cell.hmat[2][0]*(double)k;
	  strcpy(atom.asp[i+natom_base*k], atom.asp[i]);
	  atom.mfx[i+natom_base*k] = atom.mfx[i];
	  atom.mfy[i+natom_base*k] = atom.mfy[i];
	  atom.mfz[i+natom_base*k] = atom.mfz[i];
	}
      }
    } else {
      int kk=0;
      for (int k=-ix+1; k<ix; k++) {
	if (k!=0) {
	  kk++;
	  for (int i=1; i<=natom_base; i++) {
	    atom.rx[i+natom_base*kk] = atom.rx[i] + cell.hmat[0][0]*(double)k;
	    atom.ry[i+natom_base*kk] = atom.ry[i] + cell.hmat[1][0]*(double)k;
	    atom.rz[i+natom_base*kk] = atom.rz[i] + cell.hmat[2][0]*(double)k;
	    strcpy(atom.asp[i+natom_base*kk], atom.asp[i]);
	    atom.mfx[i+natom_base*kk] = atom.mfx[i];
	    atom.mfy[i+natom_base*kk] = atom.mfy[i];
	    atom.mfz[i+natom_base*kk] = atom.mfz[i];
	  }
	}
      }
    }
  }
  if (mode == 0) {
    natom_base *= ix;
  } else {
    natom_base *= (ix*2-1);
  }
  if (iy > 1) {
    if (mode == 0) {
      for (int k=1; k<iy; k++) {
	for (int i=1; i<=natom_base; i++) {
	  atom.rx[i+natom_base*k] = atom.rx[i] + cell.hmat[0][1]*(double)k;
	  atom.ry[i+natom_base*k] = atom.ry[i] + cell.hmat[1][1]*(double)k;
	  atom.rz[i+natom_base*k] = atom.rz[i] + cell.hmat[2][1]*(double)k;
	  strcpy(atom.asp[i+natom_base*k], atom.asp[i]);
	  atom.mfx[i+natom_base*k] = atom.mfx[i];
	  atom.mfy[i+natom_base*k] = atom.mfy[i];
	  atom.mfz[i+natom_base*k] = atom.mfz[i];
	}
      }
    } else {
      int kk=0;
      for (int k=-iy+1; k<iy; k++) {
	if (k!=0) {
	  kk++;
	  for (int i=1; i<=natom_base; i++) {
	    atom.rx[i+natom_base*kk] = atom.rx[i] + cell.hmat[0][1]*(double)k;
	    atom.ry[i+natom_base*kk] = atom.ry[i] + cell.hmat[1][1]*(double)k;
	    atom.rz[i+natom_base*kk] = atom.rz[i] + cell.hmat[2][1]*(double)k;
	    strcpy(atom.asp[i+natom_base*kk], atom.asp[i]);
	    atom.mfx[i+natom_base*kk] = atom.mfx[i];
	    atom.mfy[i+natom_base*kk] = atom.mfy[i];
	    atom.mfz[i+natom_base*kk] = atom.mfz[i];
	  }
	}
      }
    }
  }
  if (mode == 0) {
    natom_base *= iy;
  } else {
    natom_base *= (iy*2-1);
  }
  if (iz > 1) {
    if (mode == 0) {
      for (int k=1; k<iz; k++) {
	for (int i=1; i<=natom_base; i++) {
	  atom.rx[i+natom_base*k] = atom.rx[i] + cell.hmat[0][2]*(double)k;
	  atom.ry[i+natom_base*k] = atom.ry[i] + cell.hmat[1][2]*(double)k;
	  atom.rz[i+natom_base*k] = atom.rz[i] + cell.hmat[2][2]*(double)k;
	  strcpy(atom.asp[i+natom_base*k], atom.asp[i]);
	  atom.mfx[i+natom_base*k] = atom.mfx[i];
	  atom.mfy[i+natom_base*k] = atom.mfy[i];
	  atom.mfz[i+natom_base*k] = atom.mfz[i];
	}
      }
    } else {
      int kk=0;
      for (int k=-iz+1; k<iz; k++) {
	if (k!=0) {
	  kk++;
	  for (int i=1; i<=natom_base; i++) {
	    atom.rx[i+natom_base*kk] = atom.rx[i] + cell.hmat[0][2]*(double)k;
	    atom.ry[i+natom_base*kk] = atom.ry[i] + cell.hmat[1][2]*(double)k;
	    atom.rz[i+natom_base*kk] = atom.rz[i] + cell.hmat[2][2]*(double)k;
	    strcpy(atom.asp[i+natom_base*kk], atom.asp[i]);
	    atom.mfx[i+natom_base*kk] = atom.mfx[i];
	    atom.mfy[i+natom_base*kk] = atom.mfy[i];
	    atom.mfz[i+natom_base*kk] = atom.mfz[i];
	  }
	}
      }
    }
  }
  // cell multiplication
  for (int i=0; i<3; i++) {
    if (mode == 0) {
      cell.hmat[i][0] *= (double)ix;
      cell.hmat[i][1] *= (double)iy;
      cell.hmat[i][2] *= (double)iz;
    } else {
      cell.hmat[i][0] *= (double)(ix*2-1);
      cell.hmat[i][1] *= (double)(iy*2-1);
      cell.hmat[i][2] *= (double)(iz*2-1);
    }
  }
  // deallocate arrays for stroing data
  delete[] rxs; delete[] rys; delete[] rzs;
  for (int i=0; i<=natom_org; i++) { delete[] asps[i]; }
  delete[] asps;
  delete[] anums;
  delete[] mfxs; delete[] mfys; delete[] mfzs;

  // general initialization...
  objects = glGenLists(atom.natom*3); //printf("objects = %d\n",objects);
  for (int i=1; i<=atom.natom*3; i++) { memcpy(color[i],yellow,sizeof(GLfloat)*4); }
  set_atom_color();
  set_atom_weight();
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
  istep = 0;
  ex = 0.0;
  pot_initialize_all();
  bookkeep();
  potential();
  cellx=sqrt(cell.hmat[0][0]*cell.hmat[0][0]+cell.hmat[1][0]*cell.hmat[1][0]+cell.hmat[2][0]*cell.hmat[2][0])/ang;
  celly=sqrt(cell.hmat[0][1]*cell.hmat[0][1]+cell.hmat[1][1]*cell.hmat[1][1]+cell.hmat[2][1]*cell.hmat[2][1])/ang;
  cellz=sqrt(cell.hmat[0][2]*cell.hmat[0][2]+cell.hmat[1][2]*cell.hmat[1][2]+cell.hmat[2][2]*cell.hmat[2][2])/ang;
  cell1x=cell.hmat[0][0]/ang;cell1y=cell.hmat[1][0]/ang;cell1z=cell.hmat[2][0]/ang;
  cell2x=cell.hmat[0][1]/ang;cell2y=cell.hmat[1][1]/ang;cell2z=cell.hmat[2][1]/ang;
  cell3x=cell.hmat[0][2]/ang;cell3y=cell.hmat[1][2]/ang;cell3z=cell.hmat[2][2]/ang;
  f_max=atom.Fmax()/eV*ang;
  epotatom=atom.epotsum/eV/atom.natom;
  
  printf("Number of atoms: %d --> %d\n",natom_org,atom.natom);
}
