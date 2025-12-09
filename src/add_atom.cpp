#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#if defined(_WIN32) || defined(_WIN64)
#include <Windows.h>
#endif
#include <GL/gl.h>
#include "myheader.h"

void potential();
void bookkeep();
void pot_initialize_all();
void set_atom_weight();
void set_atom_color();
int atom_number(char* at);
void deallocate_arrays();
void allocate_arrays();
extern GLuint objects;
extern GLfloat **color;
extern GLfloat yellow[];

void add_atom(const char *aasp, float arx, float ary, float arz)
{
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
  //if (glIsList(objects+1)) {glDeleteLists(objects, atom.natom*3+1); } //printf("glDeleteLists\n"); }
  if (glIsList(objects)) {glDeleteLists(objects, atom.natom*3); } //printf("glDeleteLists\n"); }
  deallocate_arrays();
  int natom_org = atom.natom;
  atom.natom++;
  allocate_arrays();
  // copy original coordinates and others
  for (int i=1; i<=atom.natom; i++) {
    if (i<atom.natom) {
      atom.rx[i] = rxs[i]; atom.ry[i] = rys[i]; atom.rz[i] = rzs[i];
      strcpy(atom.asp[i], asps[i]);
      atom.anum[i] = anums[i];
      atom.mfx[i] = mfxs[i]; atom.mfy[i] = mfys[i]; atom.mfz[i] = mfzs[i];
    } else {
      atom.rx[i] = (double)arx*ang; atom.ry[i] = (double)ary*ang; atom.rz[i] = (double)arz*ang;
      strcpy(atom.asp[i], aasp);
      atom.anum[i] = atom_number(atom.asp[i]);
      atom.mfx[i] = false; atom.mfy[i] = false; atom.mfz[i] = false;
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
  pot_initialize_all();
  bookkeep();
  potential();
}


