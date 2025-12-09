#include <iostream>
#if defined(_WIN32) || defined(_WIN64)
#define NOMINMAX
#include <Windows.h>
#endif
#include <GL/gl.h>
#include "myheader.h"

void potential();
void bookkeep();
void pot_initialize_all();
void inverse(double mat[][3], double imat[][3]);
void rigid_shift(float z, float dx, float dy, float dz);
void rigid_shift_cfg(float z, float dx, float dy, float dz);
int rigid_shift_relax(float z, int maxstep);
void md();
extern int itolfor, relax_accel, relax_algo;
extern float tolfor;
extern int rigid_shift_relax_maxstep;

void rigid_shift_xy(float z, float dx, float dy, float dz, int ix, int iy)
{
  FILE *gsffile = fopen("gsf.d","w");
  for (int i=1; i<=ix; i++) {
    for (int j=1; j<=iy; j++) {
      rigid_shift(z,dx,0,0);
      int ii = rigid_shift_relax(z,rigid_shift_relax_maxstep);
      fprintf(gsffile, "%f %f %20.10e\n",
	      cell.hmat[0][2]/ang, cell.hmat[1][2]/ang, epotatom);
      fflush(gsffile);
      printf("GSF calc: X: %d Y: %d  Relaxation took %d steps.\n",i,j,ii);
    }
    for (int j=1; j<=iy; j++) {
      rigid_shift_cfg(z,-dx,0,0);
    }
    rigid_shift_cfg(z,0,dy,0);
    fprintf(gsffile, "\n"); fflush(gsffile);
  }
  for (int i=1; i<=ix; i++) {
    rigid_shift_cfg(z,0,-dy,0);
  }
  fclose(gsffile);

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
}

void rigid_shift(float z, float dx, float dy, float dz)
{
  rigid_shift_cfg(z, dx, dy, dz);
  /*
  double qz;
  cell.hmat[0][2] += (double)dx*ang;
  cell.hmat[1][2] += (double)dy*ang;
  cell.hmat[2][2] += (double)dz*ang;
  inverse(cell.hmat,cell.hinmat);
  
  for (int i=1; i<=atom.natom; i++) {
    qz = cell.hinmat[2][0]*atom.rx[i]
      +  cell.hinmat[2][1]*atom.ry[i]
      +  cell.hinmat[2][2]*atom.rz[i];
    if (qz > z) {
      atom.rx[i] += (double)dx*ang;
      atom.ry[i] += (double)dy*ang;
      atom.rz[i] += (double)dz*ang;
    }
  }
  */

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
  //printf("%f %f %20.10e\n", cell.hmat[0][2]/ang, cell.hmat[1][2]/ang, epotatom);
}

void rigid_shift_cfg(float z, float dx, float dy, float dz)
{
  double qz;
  cell.hmat[0][2] += (double)dx*ang;
  cell.hmat[1][2] += (double)dy*ang;
  cell.hmat[2][2] += (double)dz*ang;
  inverse(cell.hmat,cell.hinmat);
  
  for (int i=1; i<=atom.natom; i++) {
    qz = cell.hinmat[2][0]*atom.rx[i]
      +  cell.hinmat[2][1]*atom.ry[i]
      +  cell.hinmat[2][2]*atom.rz[i];
    if (qz > z) {
      atom.rx[i] += (double)dx*ang;
      atom.ry[i] += (double)dy*ang;
      atom.rz[i] += (double)dz*ang;
    }
  }
}

int rigid_shift_relax(float z, int maxstep)
{
  bool *mfx0, *mfy0, *mfz0; // to store mf{x,y,z}
  double qz;
  mfx0 = new bool[atom.natom+1];
  mfy0 = new bool[atom.natom+1];
  mfz0 = new bool[atom.natom+1];
  for (int i=0; i<=atom.natom; i++) {
    mfx0[i] = atom.mfx[i]; mfy0[i] = atom.mfy[i]; mfz0[i] = atom.mfz[i];
  }
  for (int i=0; i<=atom.natom; i++) {
    qz = cell.hinmat[2][0]*atom.rx[i]
      +  cell.hinmat[2][1]*atom.ry[i]
      +  cell.hinmat[2][2]*atom.rz[i];
    if (qz > z) {
      atom.mfx[i] = true; atom.mfy[i] = true; atom.mfz[i] = false;
    } else {
      atom.mfx[i] = true; atom.mfy[i] = true; atom.mfz[i] = true;
    }
  }
  ensemble = 3;
  relax_accel = 0;
  relax_algo = 0;
  mdmotion = 1;
  int icount = 0;
  //pot_initialize_all();
  bookkeep();
  potential();
  double maxfz;
  while ((mdmotion == 1)&&(icount < maxstep)) {
    maxfz = 0.0;
    md(); istep++;
    for (int i=1; i<=atom.natom; i++) {
      if (!atom.mfz[i]) {
      	maxfz = std::max(maxfz,atom.fz[i]*atom.fz[i]);
      }
    }
    if (sqrt(maxfz)/eV*ang < 0.01) {
      mdmotion = 0;
    }
    icount++;
  }
  maxfz = sqrt(maxfz);
  printf("rigid_shift_relax residual fz: %f[eV/ang]\n",maxfz/eV*ang);
  mdmotion = 0;
  // restore setting
  for (int i=0; i<=atom.natom; i++) {
    atom.mfx[i] = mfx0[i]; atom.mfy[i] = mfy0[i]; atom.mfz[i] = mfz0[i];
  }
  delete[] mfx0; delete[] mfy0; delete[] mfz0;
  return icount;
}
