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
void out_of_cell();
void bookkeep();

void forceconst(double **fcmat)
{
  //double del = 1.0e-2;
  double del = 0.05*ang;
  //double **fcmat;
  //fcmat = new double*[atom.natom*3];
  //for (int i=0; i<atom.natom*3; i++) {
  //  fcmat[i] = new double[atom.natom*3];
  //}
  for (int i=0; i<atom.natom*3; i++) {
    for (int j=0; j<atom.natom*3; j++) {
      fcmat[i][j] = 0;
    }
  }
  for (int ia=0; ia<atom.natom; ia++) {
    //x
    atom.rx[ia+1] += del;
    potential();
    for (int ib=0; ib<atom.natom; ib++) {
      fcmat[ia*3+0][ib*3+0] -= atom.fx[ib+1]/del/2;
      fcmat[ia*3+0][ib*3+1] -= atom.fy[ib+1]/del/2;
      fcmat[ia*3+0][ib*3+2] -= atom.fz[ib+1]/del/2;
    }
    atom.rx[ia+1] -= del*2;
    potential();
    for (int ib=0; ib<atom.natom; ib++) {
      fcmat[ia*3+0][ib*3+0] += atom.fx[ib+1]/del/2;
      fcmat[ia*3+0][ib*3+1] += atom.fy[ib+1]/del/2;
      fcmat[ia*3+0][ib*3+2] += atom.fz[ib+1]/del/2;
    }
    atom.rx[ia+1] += del;
    //y
    atom.ry[ia+1] += del;
    potential();
    for (int ib=0; ib<atom.natom; ib++) {
      fcmat[ia*3+1][ib*3+0] -= atom.fx[ib+1]/del/2;
      fcmat[ia*3+1][ib*3+1] -= atom.fy[ib+1]/del/2;
      fcmat[ia*3+1][ib*3+2] -= atom.fz[ib+1]/del/2;
    }
    atom.ry[ia+1] -= del*2;
    potential();
    for (int ib=0; ib<atom.natom; ib++) {
      fcmat[ia*3+1][ib*3+0] += atom.fx[ib+1]/del/2;
      fcmat[ia*3+1][ib*3+1] += atom.fy[ib+1]/del/2;
      fcmat[ia*3+1][ib*3+2] += atom.fz[ib+1]/del/2;
    }
    atom.ry[ia+1] += del;
    //z
    atom.rz[ia+1] += del;
    potential();
    for (int ib=0; ib<atom.natom; ib++) {
      fcmat[ia*3+2][ib*3+0] -= atom.fx[ib+1]/del/2;
      fcmat[ia*3+2][ib*3+1] -= atom.fy[ib+1]/del/2;
      fcmat[ia*3+2][ib*3+2] -= atom.fz[ib+1]/del/2;
    }
    atom.rz[ia+1] -= del*2;
    potential();
    for (int ib=0; ib<atom.natom; ib++) {
      fcmat[ia*3+2][ib*3+0] += atom.fx[ib+1]/del/2;
      fcmat[ia*3+2][ib*3+1] += atom.fy[ib+1]/del/2;
      fcmat[ia*3+2][ib*3+2] += atom.fz[ib+1]/del/2;
    }
    atom.rz[ia+1] += del;
  }

  //for (int i=0; i<atom.natom; i++) {
  //  delete[] fcmat[i];
  //}
  //delete[] fcmat;
}

