#include <stdio.h>
#include <math.h>
#include "myheader.h"

void show_e_force_all()
{
  printf("Energy and Force output in misc.d\n");
  fprintf(miscfile,"Position (A)\n");
  for (int i=1; i<=atom.natom; i++) {
    fprintf(miscfile,"%4d [%2s] %25.15e %25.15e %25.15e\n", 
	   i,atom.asp[i],atom.rx[i]/ang,atom.ry[i]/ang,atom.rz[i]/ang); }
  fprintf(miscfile,"\n");
  fprintf(miscfile,"Force (eV/A)\n");
  for (int i=1; i<=atom.natom; i++) {
    fprintf(miscfile,"%4d [%2s] %25.15e %25.15e %25.15e\n",
	   i,atom.asp[i],atom.fx[i]/eV*ang,atom.fy[i]/eV*ang,atom.fz[i]/eV*ang); }
  fprintf(miscfile,"\n");
  fprintf(miscfile,"Energy (eV)\n");
  for (int i=1; i<=atom.natom; i++) {
    fprintf(miscfile,"%4d [%2s] %25.15e\n",
	   i,atom.asp[i],atom.epot[i]/eV); }
  fprintf(miscfile,"\n");
  fflush(miscfile);
}
