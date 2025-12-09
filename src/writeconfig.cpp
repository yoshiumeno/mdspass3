#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include "myheader.h"

#define dsq(a) ((a)*(a))

void inverse(double mat[][3], double imat[][3]);
int count_num_species(int nesp[50], char **aesp);

//char tf[2] = "T";

void FT(bool val, char *tf) {
  if (val == true)  { strcpy(tf,"F"); }
  if (val == false) { strcpy(tf,"T"); }
  //return tf;
}

void writeconfig(const char* fname)
{
  printf ("Write config file = %s\n", fname);
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
  char mfx[2], mfy[2], mfz[2];
  //  fp = fopen("config.end", "w");
  fp = fopen(fname, "w");
  strcpy(species,atom.asp[1]);
  int fixmode = 0; int groupmode = 0;
  for (int i=1; i<=atom.natom; i++) {
    if (atom.mfx[i] == true) { fixmode = 1; }
    if (atom.mfy[i] == true) { fixmode = 1; }
    if (atom.mfz[i] == true) { fixmode = 1; }
    if (atom.group[i] != 0) { groupmode = 1; }
  }
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
    if ((fixmode==0)&&(groupmode==0)) {
      fprintf(fp," %20.15e %20.15e %20.15e %s\n",xx,yy,zz,atom.asp[i]);
    } else if ((fixmode==1)&&(groupmode==0)) {
      FT(atom.mfx[i],mfx); FT(atom.mfy[i],mfy); FT(atom.mfz[i],mfz);
      fprintf(fp," %20.15e %20.15e %20.15e %s %s %s %s\n",xx,yy,zz,
      mfx,mfy,mfz,atom.asp[i]);
    } else if ((fixmode==0)&&(groupmode==1)) {
      fprintf(fp," %20.15e %20.15e %20.15e %s %d\n",xx,yy,zz,
      atom.asp[i],atom.group[i]);
    } else if ((fixmode==1)&&(groupmode==1)) {
      FT(atom.mfx[i],mfx); FT(atom.mfy[i],mfy); FT(atom.mfz[i],mfz);
      fprintf(fp," %20.15e %20.15e %20.15e %s %s %s %s %d\n",xx,yy,zz,
      mfx,mfy,mfz,atom.asp[i],atom.group[i]);
    }
  }

  fclose(fp);
}

void writeconfig_abs(const char* fname)
{
  printf ("Write config file = %s\n", fname);
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
  printf("%s %s\n",atom.potential_func,atom.potential_arg);
  fprintf(fp," %d\n",atom.natom);
  fprintf(fp," %e\n",cell.alat);
  xx=cell.hmat[0][0]/cell.alat/1.0e-10; yy=cell.hmat[1][0]/cell.alat/1.0e-10; zz=cell.hmat[2][0]/cell.alat/1.0e-10; 
  fprintf(fp," %20.15e %20.15e %20.15e\n",xx,yy,zz);
  xx=cell.hmat[0][1]/cell.alat/1.0e-10; yy=cell.hmat[1][1]/cell.alat/1.0e-10; zz=cell.hmat[2][1]/cell.alat/1.0e-10; 
  fprintf(fp," %20.15e %20.15e %20.15e\n",xx,yy,zz);
  xx=cell.hmat[0][2]/cell.alat/1.0e-10; yy=cell.hmat[1][2]/cell.alat/1.0e-10; zz=cell.hmat[2][2]/cell.alat/1.0e-10; 
  fprintf(fp," %20.15e %20.15e %20.15e\n",xx,yy,zz);
  fprintf(fp," %s\n","abs");
  inverse(cell.hmat, hinmat);
  for (int i=1; i<=atom.natom; i++) {
    fprintf(fp," %20.15e %20.15e %20.15e %s\n",atom.rx[i]/ang,atom.ry[i]/ang,atom.rz[i]/ang,atom.asp[i]);
  }

  fclose(fp);
}

void writeposcar(const char* fname)
{
  printf ("Write POSCAR file = %s\n", fname);
  FILE *fp;
  double xx, yy, zz;
  double hinmat[3][3];
  char species[3];
  int nesp[50];
  char **aesp;
  aesp = new char*[atom.natom+1];
  for (int i=0; i<=atom.natom; i++) { aesp[i] = new char[3]; }

  fp = fopen(fname, "w");
  int nsp = count_num_species(nesp,aesp);
  for (int i=1; i<=nsp; i++) {
    fprintf(fp, " %s", aesp[i]);
  }
  fprintf(fp, "\n");
  fprintf(fp, "1.000\n");
  for (int j=0; j<3; j++) {
    double xx = cell.hmat[0][j]/ang;
    double yy = cell.hmat[1][j]/ang;
    double zz = cell.hmat[2][j]/ang;
    fprintf(fp, " %20.15e  %20.15e  %20.15e \n",xx,yy,zz);
  }
  for (int i=1; i<=nsp; i++) {
    fprintf(fp, " %d", nesp[i]);
  }
  fprintf(fp, "\n");
  fprintf(fp, " direct\n");
  inverse(cell.hmat, hinmat);
  for (int j=1; j<=nsp; j++) {
    for (int i=1; i<=atom.natom; i++) {
      if (strcmp(atom.asp[i],aesp[j])==0) {
	double xx = hinmat[0][0]*atom.rx[i]+hinmat[0][1]*atom.ry[i]+hinmat[0][2]*atom.rz[i];
	double yy = hinmat[1][0]*atom.rx[i]+hinmat[1][1]*atom.ry[i]+hinmat[1][2]*atom.rz[i];
	double zz = hinmat[2][0]*atom.rx[i]+hinmat[2][1]*atom.ry[i]+hinmat[2][2]*atom.rz[i];
	fprintf(fp, " %20.15e  %20.15e  %20.15e \n",xx,yy,zz);
      }
    }
  }
  fclose(fp);
  for (int i=0; i<=atom.natom; i++) { delete[] aesp[i]; }
  delete[] aesp;
}

void writecfgs(const char* fname) // Standard cfg
{
  printf ("Write standard cfg file = %s\n", fname);
  double hinmat[3][3];
  FILE* fp;
  fp = fopen(fname, "w");
  fprintf(fp, "Number of particles = %d\n", atom.natom);
  fprintf(fp, "A = %f Angstrom\n", cell.alat);
  fprintf(fp, "H0(1,1) = %f\n", cell.hmat[0][0]/ang/cell.alat);
  fprintf(fp, "H0(1,2) = %f\n", cell.hmat[1][0]/ang/cell.alat);
  fprintf(fp, "H0(1,3) = %f\n", cell.hmat[2][0]/ang/cell.alat);
  fprintf(fp, "H0(2,1) = %f\n", cell.hmat[0][1]/ang/cell.alat);
  fprintf(fp, "H0(2,2) = %f\n", cell.hmat[1][1]/ang/cell.alat);
  fprintf(fp, "H0(2,3) = %f\n", cell.hmat[2][1]/ang/cell.alat);
  fprintf(fp, "H0(3,1) = %f\n", cell.hmat[0][2]/ang/cell.alat);
  fprintf(fp, "H0(3,2) = %f\n", cell.hmat[1][2]/ang/cell.alat);
  fprintf(fp, "H0(3,3) = %f\n", cell.hmat[2][2]/ang/cell.alat);
  
  inverse(cell.hmat,hinmat);
  double xx, yy, zz;
  double am = 1.6726485e-27;
  for (int i=1; i<=atom.natom; i++) {
    xx = hinmat[0][0]*atom.rx[i] + hinmat[0][1]*atom.ry[i] + hinmat[0][2]*atom.rz[i];
    yy = hinmat[1][0]*atom.rx[i] + hinmat[1][1]*atom.ry[i] + hinmat[1][2]*atom.rz[i];
    zz = hinmat[2][0]*atom.rx[i] + hinmat[2][1]*atom.ry[i] + hinmat[2][2]*atom.rz[i];
    fprintf(fp, "%f %s %f %f %f 0. 0. 0\n",atom.wm[i]/am,atom.asp[i],xx,yy,zz);
  }
  fclose(fp);
}
void writecfge(const char* fname) // Extended cfg
{
  printf ("Write extended cfg file = %s\n", fname);
  double hinmat[3][3];
  FILE* fp;
  fp = fopen(fname, "w");
  fprintf(fp, "Number of particles = %d\n", atom.natom);
  fprintf(fp, "A = %f Angstrom\n", cell.alat);
  fprintf(fp, "H0(1,1) = %f\n", cell.hmat[0][0]/ang/cell.alat);
  fprintf(fp, "H0(1,2) = %f\n", cell.hmat[1][0]/ang/cell.alat);
  fprintf(fp, "H0(1,3) = %f\n", cell.hmat[2][0]/ang/cell.alat);
  fprintf(fp, "H0(2,1) = %f\n", cell.hmat[0][1]/ang/cell.alat);
  fprintf(fp, "H0(2,2) = %f\n", cell.hmat[1][1]/ang/cell.alat);
  fprintf(fp, "H0(2,3) = %f\n", cell.hmat[2][1]/ang/cell.alat);
  fprintf(fp, "H0(3,1) = %f\n", cell.hmat[0][2]/ang/cell.alat);
  fprintf(fp, "H0(3,2) = %f\n", cell.hmat[1][2]/ang/cell.alat);
  fprintf(fp, "H0(3,3) = %f\n", cell.hmat[2][2]/ang/cell.alat);
  fprintf(fp, ".NO_VELOCITY.\n");
  fprintf(fp, "entry_count = 7\n");
  fprintf(fp, "auxiliary[0] = kine [eV]\n");
  fprintf(fp, "auxiliary[1] = pote [eV]\n");
  fprintf(fp, "auxiliary[2] = hydro [MPa]\n");
  fprintf(fp, "auxiliary[3] = mises [MPa]\n");
  
  inverse(cell.hmat,hinmat);
  double xx, yy, zz;
  double am = 1.6726485e-27;
  double ekin, shydro, smises;
  for (int i=1; i<=atom.natom; i++) {
    xx = hinmat[0][0]*atom.rx[i] + hinmat[0][1]*atom.ry[i] + hinmat[0][2]*atom.rz[i];
    yy = hinmat[1][0]*atom.rx[i] + hinmat[1][1]*atom.ry[i] + hinmat[1][2]*atom.rz[i];
    zz = hinmat[2][0]*atom.rx[i] + hinmat[2][1]*atom.ry[i] + hinmat[2][2]*atom.rz[i];
    fprintf(fp, "%f\n %s\n",atom.wm[i]/am,atom.asp[i]);
    ekin = atom.wm[i]/2.0
      *(atom.vx[i]*atom.vx[i]+atom.vy[i]*atom.vy[i]+atom.vz[i]*atom.vz[i]);
    shydro = atom.satom[i][0][0] + atom.satom[i][1][1] + atom.satom[i][2][2];
    smises
      = dsq(atom.satom[i][0][0]-atom.satom[i][1][1])
      + dsq(atom.satom[i][1][1]-atom.satom[i][2][2])
      + dsq(atom.satom[i][2][2]-atom.satom[i][0][0])
      + 6.0 * dsq(atom.satom[i][0][1])
      + 6.0 * dsq(atom.satom[i][1][2])
      + 6.0 * dsq(atom.satom[i][2][0]);
    smises = sqrt(smises / 2.0);
    fprintf(fp, "%f %f %f  %f %f  %f %f\n"
	    ,xx,yy,zz,ekin/eV,atom.epot[i]/eV,shydro/MPa,smises/MPa);
  }
  fclose(fp);
}

int count_num_species(int nesp[50], char **aesp)
{
  int isnewspecies, atom_number;
  int nsp = 1;
  //  int nesp[50];

  strcpy(aesp[1],atom.asp[1]);
  for (int i=2; i<=atom.natom; i++) {
    // is this new species of atom?
    isnewspecies=1;
    for (int j=1; j<=i-1; j++) {
      if (strcmp(atom.asp[j],atom.asp[i])==0) {
	isnewspecies=0;
	break;
      }
    }
    if (isnewspecies==1) {
      nsp++;
      strcpy(aesp[nsp],atom.asp[i]);
    }
  }
        
  for (int i=1; i<=nsp; i++) {
    nesp[i]=0;
    for (int j=1; j<=atom.natom; j++) {
      if (strcmp(atom.asp[j],aesp[i])==0) { nesp[i]++; }
    }
  }
  return nsp;
}

