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

void slice(float r1min, float r1max, float r2min, float r2max, float r3min, float r3max)
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

  double xx,yy,zz,xxx,yyy,zzz;
  int icnt=0;
  float d=0.01;
  for (int i=1; i<=atom.natom; i++) {
    xx = atom.rx[i];
    yy = atom.ry[i];
    zz = atom.rz[i];
    xxx = cell.hinmat[0][0]*xx + cell.hinmat[0][1]*yy + cell.hinmat[0][2]*zz;
    yyy = cell.hinmat[1][0]*xx + cell.hinmat[1][1]*yy + cell.hinmat[1][2]*zz;
    zzz = cell.hinmat[2][0]*xx + cell.hinmat[2][1]*yy + cell.hinmat[2][2]*zz;
    if ((xxx>=r1min-d)&&(xxx<=r1max+d)&&(yyy>=r2min-d)&&(yyy<=r2max+d)&&(zzz>=r3min-d)&&(zzz<=r3max+d)) {
      icnt++;
      rxs[icnt] = atom.rx[i]; rys[icnt] = atom.ry[i]; rzs[icnt] = atom.rz[i];
      strcpy(asps[icnt], atom.asp[i]);
      anums[icnt] = atom.anum[i];
      mfxs[icnt] = atom.mfx[i]; mfys[icnt] = atom.mfy[i]; mfzs[icnt] = atom.mfz[i];
    }
  }
  if (glIsList(objects+1)) {glDeleteLists(objects, atom.natom*3+1); } //printf("glDeleteLists\n"); }
  deallocate_arrays();
  int natom_org = atom.natom;
  atom.natom = icnt;
  allocate_arrays();
  for (int i=1; i<=atom.natom; i++) {
    atom.rx[i] = rxs[i]; atom.ry[i] = rys[i]; atom.rz[i] = rzs[i];
    strcpy(atom.asp[i], asps[i]);
    atom.anum[i] = anums[i];
    atom.mfx[i] = mfxs[i]; atom.mfy[i] = mfys[i]; atom.mfz[i] = mfzs[i];
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

