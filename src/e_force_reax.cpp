#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <fstream>
#include "myheader.h"
#include "ReaxFF.h"

void resetmat(double a[3][3]);

void e_force_reax()
{
  if (myreax.initialize) {
    if (ReaxFF::GetNumberOfAtoms()!=atom.natom) {
      myreax.myreaxpot.ReadFile(myreax.potfname);
      printf("=========================\n");
      printf("Resize ReaxFF conf needed\n");
      printf("=========================\n");
      //myreax.myreaxconf.ReadFile("pot/ReaxFF.config2");
      myreax.myreaxconf.ModConf(atom.natom,3, // <-- 3 should be modified later
				cell.hmat[0][0],cell.hmat[0][1],cell.hmat[0][2],
				cell.hmat[1][0],cell.hmat[1][1],cell.hmat[1][2],
				cell.hmat[2][0],cell.hmat[2][1],cell.hmat[2][2],
				atom.rx,atom.ry,atom.rz,atom.asp);
      ReaxFF::Initialize(myreax.myreaxpot,myreax.myreaxconf,myreax.myreaxinput);
      printf("ReaxFF initialized (system has %d atoms)\n",ReaxFF::GetNumberOfAtoms());
    }
    
  }
  
  double ha[3][3];
  for (int i=0; i<3; i++) { for (int j=0; j<3; j++) { ha[i][j] = cell.hmat[i][j]/ang; } }
  ReaxFF::SetLattice(ha[0][0],ha[0][1],ha[0][2],
		     ha[1][0],ha[1][1],ha[1][2],
		     ha[2][0],ha[2][1],ha[2][2]);
  if (ReaxFF::GetNumberOfAtoms()!=atom.natom) {
    printf("#########################\n");
    printf("Resize ReaxFF conf needed\n");
    printf("#########################\n");
    exit(0);
  }
  for (int i=1; i<=atom.natom; i++) {
    ReaxFF::SetRAtomPos(i,atom.rx[i]/ang,atom.ry[i]/ang,atom.rz[i]/ang);
  }
  //ReaxFF::CalcAll();
  //printf("Etot = %e\n",atom.epotsum);
  //ReaxFF::GetNumberOfAtoms();
  ReaxFF::CalcAll(atom.epotsum);
  for (int i=1; i<=atom.natom; i++) {
    atom.epot[i] = atom.epotsum / (double)atom.natom;
    ReaxFF::GetForce(i,atom.fx[i],atom.fy[i],atom.fz[i]);
    atom.fx[i] *= eV/ang; atom.fy[i] *= eV/ang; atom.fz[i] *= eV/ang;
    for(int iDim=0;iDim<3;iDim++){ for(int jDim=0;jDim<3;jDim++){
	atom.satom[i][iDim][jDim] = 0.0; } }
    ReaxFF::GetStress(i,atom.satom[i][0][0],atom.satom[i][1][1],atom.satom[i][2][2]);
    ReaxFF::GetDStress(i,atom.satom[i][0][1],atom.satom[i][1][2],atom.satom[i][2][0]);
    ReaxFF::GetDStress(i,atom.satom[i][1][0],atom.satom[i][2][1],atom.satom[i][0][2]);
  }
  resetmat(cell.dmat);
  for(int iDim=0;iDim<3;iDim++){
    for(int jDim=0;jDim<3;jDim++){
      cell.dmat[iDim][jDim] = 0.0;
      for(int i=1; i<=atom.natom; i++){
	atom.satom[i][iDim][jDim] *= eV;
	cell.dmat[iDim][jDim] -= atom.satom[i][iDim][jDim];
	atom.satom[i][iDim][jDim] *= (double)atom.natom / cell.volume;
      } // ILOOP
    } // for: jDim
  } // for: iDim
    //ReaxFF::ShowForce(0);
}
