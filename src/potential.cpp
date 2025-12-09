#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#if defined __linux__ || defined __APPLE__
#include <unistd.h>
#else
#include <direct.h>
#endif
#include "myheader.h"
#include "ReaxFF.h"
//#include "Input.h"
//#include "ReaxPotential.h"
//#include "Configuration.h"

//extern vector<ReaxFF> s_reax;
//extern Input input;
//extern ReaxPotential reaxpot;
//extern ConfigurationSet confSet;


void resetmat(double a[3][3]);
void e_force_morse();
void e_force_geam_am();
void e_force_brenner();
void e_force_brenner(int mode);
void e_force_tersoff();
void e_force_mishin(PotentialMode mode);
void e_force_dipole(PotentialMode mode);
void e_force_combined();
void e_force_adp();
//void e_force_vashishta(); // Kubo 20140619 X0.1
void e_force_sw();
void e_force_reax();
void e_force_sinusoidal();
void e_force_shellmodel();
void calcstresstensor();
void bookkeep();
void get_first_arg(std::string &line, std::string &arg1);
int count_arg_number(std::string line);
void set_potfile();
void set_potfile(const char* fname);
void pot_initialize_all();

extern std::string potfile;
extern char cwdname[80];

void potential_set (int control)
{
  printf("Potential: No. %d, %s\n",ipottype,atom.potstring_list[ipottype]);
  //strcpy(atom.potential_func, potstring_list[ipottype]);
  //std::string line = std::string(potstring_list[ipottype]);
  //std::string arg1, arg2; int narg = count_arg_number(line);
  std::string arg1 = std::string(atom.potstring_list[ipottype]);
  if ((atom.potarg_number[ipottype]==0)&&(!atom.potarg_readable[ipottype])) {
    strcpy(atom.potential_func, arg1.c_str());
    strcpy(atom.potential_arg, "");
  } else if (atom.potarg_readable[ipottype]) {
    strcpy(atom.potential_func, arg1.c_str());
    strcpy(atom.potential_arg, "");
  } else {
    std::string arg2 = std::string(atom.potargstring_list[ipottype][ipotarg]);
    strcpy(atom.potential_func, arg1.c_str());
    strcpy(atom.potential_arg,  arg2.c_str()); }
  book.algo = 1;
  //if (strcmp(atom.potential_func,"Tersoff")==0) { book.algo = 2; }
  if (strcmp(atom.potential_func,"Brenner")==0) { book.algo = 3; }
  if (strcmp(atom.potential_func,"AIREBO" )==0) { book.algo = 3; }
  pot_initialize_all();
}

void potential()
{
 CALC:
  if (strcmp(atom.potential_func, "Morse") == 0) {
    e_force_morse();
  } else if (strcmp(atom.potential_func, "EAM") == 0) {
    if (strcmp(atom.potential_arg, "GEAM") == 0) {
      e_force_geam_am();
    } else if (strcmp(atom.potential_arg, "Mishin") == 0) {
      e_force_mishin(NORMALPOTMODE);
    }
  } else if (strcmp(atom.potential_func, "Brenner") == 0) {
    e_force_brenner();
  } else if (strcmp(atom.potential_func, "AIREBO") == 0) {
    e_force_brenner(1);
  } else if (strcmp(atom.potential_func, "Tersoff") == 0) {
    e_force_tersoff();
    //} else if ((strcmp(atom.potential_func, "EAM") == 0)&&
    //	     (strcmp(atom.potential_arg, "Mishin") == 0)) {
    //    e_force_mishin(NORMALPOTMODE);
  } else if (strcmp(atom.potential_func, "Dipole") == 0) {
    e_force_dipole(NORMALPOTMODE);
  } else if (strcmp(atom.potential_func, "ADP") == 0) {
    e_force_adp();
  } else if (strcmp(atom.potential_func, "NiYSZ") == 0) {
    e_force_combined();
  } else if (strcmp(atom.potential_func, "AIREBO_BN") == 0) {
    e_force_brenner(1);
  } else if (strcmp(atom.potential_func, "SW") == 0) {
    e_force_sw();
    /*
    if (strcmp(atom.potential_arg, "SiC_Vashishta") == 0) {
      e_force_vashishta();
    } else if (strcmp(atom.potential_arg, "GaN_Bere") == 0) {
      e_force_vashishta();
    } else if (strcmp(atom.potential_arg, "GaN_Aichoune") == 0) {
      e_force_vashishta();
    } 
    */
  } else if (strcmp(atom.potential_func, "ReaxFF") == 0) {
    e_force_reax();
    /*
    //printf("## ReaxFF called ##\n");
    //ReaxFF::Initialize(reaxpot,confSet,input);
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
    */
  } else if (strcmp(atom.potential_func, "Sinusoidal") == 0) {
    e_force_sinusoidal();
  } else if (strcmp(atom.potential_func, "ShellModel") == 0) {
    e_force_shellmodel();
  } else {
    printf("##### Invalid potential for this system #####\n");
    mdmotion = 0;
    return;
  }
  
  // Set rcut_f
  if (log10(rcut)<-8) { rcut_f = rcut*1e10; } // Some routines use rcut in [m]
  else { rcut_f = rcut; }

  // Is frc enough?
  //printf("%f %f\n",rcut_f, book.frc/ang);
  if (!tersoff.nocutoff) {
  if (rcut_f > book.frc / ang) {
    book.frc = ((double)rcut_f + frcmar) * ang;
    frc_f = book.frc / ang;
    bookkeep();
    goto CALC;
  }
  }
  // For output in "Stresss" area
  calcstresstensor();
  strs_xx = cell.sgmmat[0][0]*1e-6;
  strs_yy = cell.sgmmat[1][1]*1e-6;
  strs_zz = cell.sgmmat[2][2]*1e-6;
  strs_xy = cell.sgmmat[0][1]*1e-6;
  strs_yz = cell.sgmmat[1][2]*1e-6;
  strs_zx = cell.sgmmat[2][0]*1e-6;
}

void set_potfile()
{
  if (strcmp(atom.potential_func, "Dipole") == 0) {
    printf("Dipole parameter file is set to %s\n",potfile.c_str());
    strcpy(dipole.fname,potfile.c_str());
    dipole.initialize = true;
  } else if (strcmp(atom.potential_func, "ADP") == 0) {
    printf("ADP parameter file is set to %s\n",potfile.c_str());
    strcpy(adp.fname,potfile.c_str());
    adp.initialize = true;
  } else if (strcmp(atom.potential_func, "Buckingham") == 0) {
    printf("Buckingham parameter file is set to %s\n",potfile.c_str());
    if (pairPot == NULL) pairPot = new Buckingham(4, 1, COMBINEDPOTMODE);
    strcpy(pairPot->fname,potfile.c_str());
    pairPot->initialized = false;
  } else if (strcmp(atom.potential_func, "ReaxFF") == 0) {
    printf("ReaxFF parameter file is set to %s\n",potfile.c_str());
  }

}
void set_potfile(const char* fname)
{
  char fname0[120] = "aaa";
  //getcwd(cwdname,80); // Abolished Oct2015
  strcpy(fname0,cwdname);
#if defined __linux__ || defined __APPLE__
  strcat(fname0,"/pot/");
#else
  strcat(fname0,"\\pot\\");
#endif
  strcat(fname0,fname);
  if (strcmp(atom.potential_func, "Dipole") == 0) {
    printf("Dipole parameter file is set to %s\n",fname);
    strcpy(dipole.fname,fname0);
    dipole.initialize = true;
  } else if (strcmp(atom.potential_func, "ADP") == 0) {
    printf("ADP parameter file is set to %s\n",fname);
    strcpy(adp.fname,fname0);
    adp.initialize = true;
  } else if (strcmp(atom.potential_func, "Pair") == 0) {
    printf("Pair parameter file is set to %s\n",fname);
    if (pairPot == NULL) pairPot = new Buckingham(4, 1, COMBINEDPOTMODE);
    strcpy(pairPot->fname,fname0);
    pairPot->initialized = false;
  }
}

void pot_initialize_all()
{
  bre.initialize = true;
  tersoff.initialize = true;
  eammis.initialize = true;
  geam.initialize = true;
  adp.initialize = true;
  dipole.initialize = true;
  sw.initialize = true;
  myreax.initialize = true;
  //ReaxFF::Initialize(reaxpot,reaxconf,input);
  /*
  printf("#############ReaxFF Initialize in potential.cpp\n");
  ReaxFF::Initialize(myreax.myreaxpot,myreax.myreaxconf,myreax.myreaxinput);
  if (ReaxFF::GetNumberOfAtoms()!=atom.natom) {
      printf("=========================\n");
      printf("Resize ReaxFF conf needed\n");
      printf("=========================\n");
      //myreax.myreaxconf.ReadFile("pot/ReaxFF.config2");
      myreax.myreaxconf.ModConf(atom.natom,3,
				cell.hmat[0][0],cell.hmat[0][1],cell.hmat[0][2],
				cell.hmat[1][0],cell.hmat[1][1],cell.hmat[1][2],
				cell.hmat[2][0],cell.hmat[2][1],cell.hmat[2][2],
				atom.rx,atom.ry,atom.rz,atom.asp
				);
      ReaxFF::Initialize(myreax.myreaxpot,myreax.myreaxconf,myreax.myreaxinput);
      printf("%d\n",ReaxFF::GetNumberOfAtoms());
  }
  */
  if (pairPot == NULL) pairPot = new Buckingham(4, 1, COMBINEDPOTMODE);
}
