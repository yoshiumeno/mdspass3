#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#if defined(_WIN32) || defined(_WIN64)
#include <Windows.h>
#endif
#include <GL/gl.h>
#include "myheader.h"

extern GLfloat **color, **color0;
extern int *iatom, *repidx;

void deallocate_arrays()
{
 // Deallocate arrays
  if (color)  { for (int i=0; i<=atom.natom*3; i++) { delete[] color[i];  } delete[] color;  color = NULL;  }
  if (color0) { for (int i=0; i<=atom.natom*3; i++) { delete[] color0[i]; } delete[] color0; color0 = NULL; }
  if (iatom)  { delete[]  iatom;  iatom = NULL; }
  if (repidx) { delete[] repidx; repidx = NULL; }
  if (atom.asp) { for (int i=0; i<=atom.natom; i++) { delete[] atom.asp[i]; } delete[] atom.asp; atom.asp = NULL; }
  if (atom.anum) { delete[] atom.anum; atom.anum = NULL; }
  if (atom.wm) { delete[] atom.wm; atom.wm = NULL; }
  if (atom.visible) { delete[] atom.visible; atom.visible = NULL; }
  if (atom.rx) { delete[] atom.rx; atom.rx = NULL; }
  if (atom.ry) { delete[] atom.ry; atom.ry = NULL; }
  if (atom.rz) { delete[] atom.rz; atom.rz = NULL; }
  if (atom.rx_float) { delete[] atom.rx_float; atom.rx_float = NULL; }
  if (atom.ry_float) { delete[] atom.ry_float; atom.ry_float = NULL; }
  if (atom.rz_float) { delete[] atom.rz_float; atom.rz_float = NULL; }
  if (atom.rx_org) { delete[] atom.rx_org; atom.rx_org = NULL; }
  if (atom.ry_org) { delete[] atom.ry_org; atom.ry_org = NULL; }
  if (atom.rz_org) { delete[] atom.rz_org; atom.rz_org = NULL; }
  if (atom.rx_p) { delete[] atom.rx_p; atom.rx_p = NULL; }
  if (atom.ry_p) { delete[] atom.ry_p; atom.ry_p = NULL; }
  if (atom.rz_p) { delete[] atom.rz_p; atom.rz_p = NULL; }
  if (atom.qx) { delete[] atom.qx; atom.qx = NULL; }
  if (atom.qy) { delete[] atom.qy; atom.qy = NULL; }
  if (atom.qz) { delete[] atom.qz; atom.qz = NULL; }
  if (atom.vx) { delete[] atom.vx; atom.vx = NULL; }
  if (atom.vy) { delete[] atom.vy; atom.vy = NULL; }
  if (atom.vz) { delete[] atom.vz; atom.vz = NULL; }
  if (atom.fx) { delete[] atom.fx; atom.fx = NULL; }
  if (atom.fy) { delete[] atom.fy; atom.fy = NULL; }
  if (atom.fz) { delete[] atom.fz; atom.fz = NULL; }
  if (atom.ax) { delete[] atom.ax; atom.ax = NULL; }
  if (atom.ay) { delete[] atom.ay; atom.ay = NULL; }
  if (atom.az) { delete[] atom.az; atom.az = NULL; }
  if (atom.bx) { delete[] atom.bx; atom.bx = NULL; }
  if (atom.by) { delete[] atom.by; atom.by = NULL; }
  if (atom.bz) { delete[] atom.bz; atom.bz = NULL; }
  if (atom.cx) { delete[] atom.cx; atom.cx = NULL; }
  if (atom.cy) { delete[] atom.cy; atom.cy = NULL; }
  if (atom.cz) { delete[] atom.cz; atom.cz = NULL; }
  if (atom.fx_l) { delete[] atom.fx_l; atom.fx_l = NULL; }
  if (atom.fy_l) { delete[] atom.fy_l; atom.fy_l = NULL; }
  if (atom.fz_l) { delete[] atom.fz_l; atom.fz_l = NULL; }
  if (atom.group) { delete[] atom.group; atom.group = NULL; }
  if (atom.fx_float) { delete[] atom.fx_float; atom.fx_float = NULL; }
  if (atom.fy_float) { delete[] atom.fy_float; atom.fy_float = NULL; }
  if (atom.fz_float) { delete[] atom.fz_float; atom.fz_float = NULL; }
  if (atom.mfx) { delete[] atom.mfx; atom.mfx = NULL; }
  if (atom.mfy) { delete[] atom.mfy; atom.mfy = NULL; }
  if (atom.mfz) { delete[] atom.mfz; atom.mfz = NULL; }
  if (atom.lock) { delete[] atom.lock; atom.lock = NULL; }
  if (atom.epot)       { delete[] atom.epot;       atom.epot = NULL;       }
  if (atom.epot_p)     { delete[] atom.epot_p;     atom.epot_p = NULL;     }
  if (atom.epot_float) { delete[] atom.epot_float; atom.epot_float = NULL; }
  if (atom.repatom) { delete[] atom.repatom; atom.repatom = NULL; }
  if (atom.elem_id) { delete[] atom.elem_id; atom.elem_id = NULL; }
  if (atom.evecx) { for (int i=0; i<=atom.natom; i++) { delete[] atom.evecx[i]; } delete[] atom.evecx; atom.evecx = NULL; }
  if (atom.evecy) { for (int i=0; i<=atom.natom; i++) { delete[] atom.evecy[i]; } delete[] atom.evecy; atom.evecx = NULL; }
  if (atom.evecz) { for (int i=0; i<=atom.natom; i++) { delete[] atom.evecz[i]; } delete[] atom.evecz; atom.evecx = NULL; }
  if (atom.satom) { for (int i=0; i<=atom.natom; i++) {
      for (int j=0; j<3; j++) { delete[] atom.satom[i][j]; } delete[] atom.satom[i]; }
    delete[] atom.satom; atom.satom = NULL; }
  if (atom.nneighbor) { delete[] atom.nneighbor; atom.nneighbor = NULL;
    for (int i=0; i<=atom.natom; i++) { delete[] atom.neighbor[i]; }
    delete[] atom.neighbor; atom.neighbor = NULL; }
  delete[] book.alistnum;
  if (cg.sgx) { for (int i=0; i<=atom.natom; i++) { delete[] cg.sgx[i]; }
    delete[] cg.sgx; cg.sgx = NULL; }
  if (cg.sgy) { for (int i=0; i<=atom.natom; i++) { delete[] cg.sgy[i]; }
    delete[] cg.sgy; cg.sgy = NULL; }
  if (cg.sgz) { for (int i=0; i<=atom.natom; i++) { delete[] cg.sgz[i]; }
    delete[] cg.sgz; cg.sgz = NULL; }
  if (cg.shx) { for (int i=0; i<=atom.natom; i++) { delete[] cg.shx[i]; }
    delete[] cg.shx; cg.shx = NULL; }
  if (cg.shy) { for (int i=0; i<=atom.natom; i++) { delete[] cg.shy[i]; }
    delete[] cg.shy; cg.shy = NULL; }
  if (cg.shz) { for (int i=0; i<=atom.natom; i++) { delete[] cg.shz[i]; }
    delete[] cg.shz; cg.shz = NULL; }
  if (cg.rstx) { delete[] cg.rstx; cg.rstx = NULL; }
  if (cg.rsty) { delete[] cg.rsty; cg.rsty = NULL; }
  if (cg.rstz) { delete[] cg.rstz; cg.rstz = NULL; }
 // Deallocate arrays end
}

void allocate_arrays()
{
  // Allocate arrays
  color  = new GLfloat*[atom.natom*3+1]; for (int i=0; i<=atom.natom*3; i++) { color[i]  = new GLfloat[4]; }
  color0 = new GLfloat*[atom.natom*3+1]; for (int i=0; i<=atom.natom*3; i++) { color0[i] = new GLfloat[4]; }
  iatom  = new int[atom.natom*2+1];
  repidx = new int[atom.natom*2+1];
  atom.asp = new char*[atom.natom+1]; for (int i=0; i<=atom.natom; i++) { atom.asp[i] = new char[3]; }
  atom.anum = new int[atom.natom+1];
  atom.wm = new double[atom.natom+1];
  atom.visible = new bool[atom.natom+1]; for (int i=0; i<=atom.natom; i++) { atom.visible[i] = true; }
  atom.rx = new double[atom.natom+1];
  atom.ry = new double[atom.natom+1];
  atom.rz = new double[atom.natom+1];
  atom.qx = new double[atom.natom+1];
  atom.qy = new double[atom.natom+1];
  atom.qz = new double[atom.natom+1];
  atom.ax = new double[atom.natom+1];
  atom.ay = new double[atom.natom+1];
  atom.az = new double[atom.natom+1];
  atom.bx = new double[atom.natom+1];
  atom.by = new double[atom.natom+1];
  atom.bz = new double[atom.natom+1];
  atom.cx = new double[atom.natom+1];
  atom.cy = new double[atom.natom+1];
  atom.cz = new double[atom.natom+1];
  atom.rx_float = new FLOAT[atom.natom+1];
  atom.ry_float = new FLOAT[atom.natom+1];
  atom.rz_float = new FLOAT[atom.natom+1];
  atom.rx_org = new double[atom.natom+1];
  atom.ry_org = new double[atom.natom+1];
  atom.rz_org = new double[atom.natom+1];
  atom.rx_p = new double[atom.natom+1];
  atom.ry_p = new double[atom.natom+1];
  atom.rz_p = new double[atom.natom+1];
  atom.vx = new double[atom.natom+1];
  atom.vy = new double[atom.natom+1];
  atom.vz = new double[atom.natom+1];
  atom.fx = new double[atom.natom+1];
  atom.fy = new double[atom.natom+1];
  atom.fz = new double[atom.natom+1];
  atom.fx_l = new double[atom.natom+1];
  atom.fy_l = new double[atom.natom+1];
  atom.fz_l = new double[atom.natom+1];
  atom.group = new int[atom.natom+1];
  atom.fx_float = new FLOAT[atom.natom+1];
  atom.fy_float = new FLOAT[atom.natom+1];
  atom.fz_float = new FLOAT[atom.natom+1];
  atom.mfx  = new bool[atom.natom+1];
  atom.mfy  = new bool[atom.natom+1];
  atom.mfz  = new bool[atom.natom+1];
  atom.lock = new bool[atom.natom+1];
  atom.epot       = new double[atom.natom+1];
  atom.epot_p     = new double[atom.natom+1];
  atom.epot_float = new FLOAT[atom.natom+1];
  atom.repatom = new int[atom.natom+1];
  atom.elem_id = new int[atom.natom+1];
  atom.evecx = new double*[atom.natom+1]; for (int i=0; i<=atom.natom; i++) { atom.evecx[i] = new double[MAXMODE]; }
  atom.evecy = new double*[atom.natom+1]; for (int i=0; i<=atom.natom; i++) { atom.evecy[i] = new double[MAXMODE]; }
  atom.evecz = new double*[atom.natom+1]; for (int i=0; i<=atom.natom; i++) { atom.evecz[i] = new double[MAXMODE]; }
  atom.satom = new double**[atom.natom+1];
  for (int i=0; i<=atom.natom; i++) { atom.satom[i]    = new double*[3];
    for (int j=0; j<3; j++)         { atom.satom[i][j] = new double[3]; } }
  atom.nneighbor = new int[atom.natom+1]; atom.neighbor = new int*[atom.natom+1];
  for (int i=0; i<=atom.natom; i++) { atom.neighbor[i] = new int[MAXNEIGHBOR]; }
  book.alistnum = new int[atom.natom+1]; book.alloc = true;
  
  cg.sgx = new double*[atom.natom+1];
  cg.sgy = new double*[atom.natom+1];
  cg.sgz = new double*[atom.natom+1];
  cg.shx = new double*[atom.natom+1];
  cg.shy = new double*[atom.natom+1];
  cg.shz = new double*[atom.natom+1];
  for (int i=0; i<=atom.natom; i++) {
    cg.sgx[i] = new double[2];
    cg.sgy[i] = new double[2];
    cg.sgz[i] = new double[2];
    cg.shx[i] = new double[2];
    cg.shy[i] = new double[2];
    cg.shz[i] = new double[2];
  }
  cg.rstx = new double[atom.natom+1];
  cg.rsty = new double[atom.natom+1];
  cg.rstz = new double[atom.natom+1];
  // Allocate arrays end
}
