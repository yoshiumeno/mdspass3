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

void move_atom_shift(float frx, float fry, float frz)
{
  for (int i=1; i<=atom.natom; i++) {
    atom.rx[i] += frx*ang;
    atom.ry[i] += fry*ang;
    atom.rz[i] += frz*ang;
    atom.rx_org[i]=atom.rx[i]; atom.ry_org[i]=atom.ry[i]; atom.rz_org[i]=atom.rz[i];
  }
  pot_initialize_all();
  bookkeep();
  potential();
}

void move_atom(int i, float frx, float fry, float frz)
{
  if ((i>0)&&(i<=atom.natom)) {
    atom.rx[i] = frx*ang;
    atom.ry[i] = fry*ang;
    atom.rz[i] = frz*ang;
    atom.rx_org[i]=atom.rx[i]; atom.ry_org[i]=atom.ry[i]; atom.rz_org[i]=atom.rz[i];
  }
  pot_initialize_all();
  bookkeep();
  potential();
}

void move_atom_shift(int i, float frx, float fry, float frz)
{
  if ((i>0)&&(i<=atom.natom)) {
    atom.rx[i] += frx*ang;
    atom.ry[i] += fry*ang;
    atom.rz[i] += frz*ang;
    atom.rx_org[i]=atom.rx[i]; atom.ry_org[i]=atom.ry[i]; atom.rz_org[i]=atom.rz[i];
  }
  pot_initialize_all();
  bookkeep();
  potential();
}

void move_atom_shift(int i, int j, float frx, float fry, float frz)
{
  if ((i>0)&&(i<=atom.natom)) {
    if ((j>0)&&(j<=atom.natom)) {
      if (j>=i) {
	for (int k=i; k<=j; k++) {
	  atom.rx[k] += frx*ang;
	  atom.ry[k] += fry*ang;
	  atom.rz[k] += frz*ang;
	  atom.rx_org[k]=atom.rx[k];
	  atom.ry_org[k]=atom.ry[k];
	  atom.rz_org[k]=atom.rz[k];
	}
      }
    }
  }
  pot_initialize_all();
  bookkeep();
  potential();
}

void move_atom(int i, const char *aasp, float frx, float fry, float frz)
{
  if ((i>0)&&(i<=atom.natom)) {
    atom.rx[i] = frx*ang;
    atom.ry[i] = fry*ang;
    atom.rz[i] = frz*ang;
    atom.rx_org[i]=atom.rx[i]; atom.ry_org[i]=atom.ry[i]; atom.rz_org[i]=atom.rz[i];
    strcpy(atom.asp[i], aasp);
    atom.anum[i] = atom_number(atom.asp[i]);
  }
  set_atom_color();
  set_atom_weight();
  pot_initialize_all();
  bookkeep();
  potential();
}


