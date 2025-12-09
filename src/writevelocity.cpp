#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include "myheader.h"

#define dsq(a) ((a)*(a))

void inverse(double mat[][3], double imat[][3]);
int count_num_species(int nesp[50], char **aesp);

void writevelocity(const char* fname)
{
  printf ("Write velocity file = %s\n", fname);
  /*
  std::ofstream foutconfig( "config.end", std::ios::out );
  foutconfig << "ATOM_TYPE" << std::endl;
  foutconfig << "POTENTIAL" << std::endl;
  foutconfig << atom.natom << std::endl;
  */
  FILE *fp;
  double xx, yy, zz;
  double hinmat[3][3];
  char species[3];
  //  fp = fopen("config.end", "w");
  fp = fopen(fname, "w");
  strcpy(species,atom.asp[1]);
  if (atom.natom > 1) {
    for (int i=2; i<=atom.natom; i++) {
      if(strcmp(atom.asp[i], species)==1) {strcpy(species,"M");} } }
  fprintf(fp," %s\n",species);
  fprintf(fp," %s %s\n",atom.potential_func,atom.potential_arg);
  fprintf(fp," %d\n",atom.natom);
  fprintf(fp," %e\n",cell.alat);
  xx=cell.hmat[0][0]/cell.alat/1.0e-10; yy=cell.hmat[1][0]/cell.alat/1.0e-10; zz=cell.hmat[2][0]/cell.alat/1.0e-10; 
  fprintf(fp," %20.15e %20.15e %20.15e\n",xx,yy,zz);
  xx=cell.hmat[0][1]/cell.alat/1.0e-10; yy=cell.hmat[1][1]/cell.alat/1.0e-10; zz=cell.hmat[2][1]/cell.alat/1.0e-10; 
  fprintf(fp," %20.15e %20.15e %20.15e\n",xx,yy,zz);
  xx=cell.hmat[0][2]/cell.alat/1.0e-10; yy=cell.hmat[1][2]/cell.alat/1.0e-10; zz=cell.hmat[2][2]/cell.alat/1.0e-10; 
  fprintf(fp," %20.15e %20.15e %20.15e\n",xx,yy,zz);
  fprintf(fp," %s\n","fra");
  inverse(cell.hmat, hinmat);
  for (int i=1; i<=atom.natom; i++) {
    xx = hinmat[0][0]*atom.rx[i] + hinmat[0][1]*atom.ry[i] + hinmat[0][2]*atom.rz[i];
    yy = hinmat[1][0]*atom.rx[i] + hinmat[1][1]*atom.ry[i] + hinmat[1][2]*atom.rz[i];
    zz = hinmat[2][0]*atom.rx[i] + hinmat[2][1]*atom.ry[i] + hinmat[2][2]*atom.rz[i];
    fprintf(fp," %20.15e %20.15e %20.15e %s\n",xx,yy,zz,atom.asp[i]);
  }

  fclose(fp);
}

