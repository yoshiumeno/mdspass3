#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<string>
#include<iostream>
#include<sstream>
#include<fstream>
#include"myheader.h"

#define ILOOP int i=1;i<=atom.natom;i++
#define JLOOP int k0=1;k0<=book.alistnum[i];k0++
#define KLOOP int k1=1;k1<=book.alistnum[i];k1++

using namespace std;

//// Funcs in Other Files.
void resetmat(double a[3][3]);
double innerproduct(double a[3], double b[3]);

//// Funcs in This File.
double Norm(double *v);
double Norm2(double *v);
int delta(int i,int j);
//void SetParam(void);
void SetParam(const char* arg);
void Initialize(void);
void DerivCOS(double *rik, double *rij, double *dcos_ri);
void DerivCOS(double *rik, double *rij, double *dcos_ri, double *dcos_rj, double *dcos_rk);
int GetNeighbor(int *bk, int *iTrans);
void ShowVector(double *v);


////////////////////////////////////////////////////////////////X
////                                                        ////X
////               Stillinger-Weber   Model                 ////X
////      implemented by Kubo originally for Vashishta      ////X
////       modified by Umeno for SW                         ////X
////                                                        ////X
////////////////////////////////////////////////////////////////X

void e_force_sw(void){ //// Principal Function //////////
  double eps = 1.0e-10; // [A]

  Initialize();
  //if(sw.initialize){SetParam();}
  if(sw.initialize){SetParam(atom.potential_arg);}

  for(ILOOP){
    if(book.alistnum[i]<=0){continue;}

    double rveci[3]; // Position of i-th Atom [A]
    rveci[0] = atom.rx[i]/ang;
    rveci[1] = atom.ry[i]/ang;
    rveci[2] = atom.rz[i]/ang;
    int iSpec = sw.AtomNum2Index[atom.anum[i]];

    for(JLOOP){
      int j;                 // Atom ID for Bookkeeping
      int iTrans[3];         // Cell Translation
      int ijPair;            // Pair ID (e.g., Si-Si = 1, Si-C = 2, C-C = 3, ...)
      double rvecj[3];       // Position of Replica of j-th Atom [A]
      double rvecij[3];      // rj-ri [A]
      double nvecij[3];      // rij/|rij| [-]
      double rij;            // |rij| [A]

      j = GetNeighbor(book.alist[i][k0],iTrans);
      int jSpec = sw.AtomNum2Index[atom.anum[j]];
      ijPair = sw.PairType[iSpec][jSpec];
      if(ijPair==0){continue;} // Skip Non-Reactive Pair.

      rvecij[0] = atom.Dx(i,j,iTrans[0],iTrans[1],iTrans[2])/ang; // rij := rj - ri.
      rvecij[1] = atom.Dy(i,j,iTrans[0],iTrans[1],iTrans[2])/ang;
      rvecij[2] = atom.Dz(i,j,iTrans[0],iTrans[1],iTrans[2])/ang;
      rij = Norm(rvecij);
      if(rij<eps){continue;} // Skip Too Close (or Identical) Pair.

      for(int iDim=0;iDim<3;iDim++){
        rvecj[iDim]  = rveci[iDim] + rvecij[iDim];
        nvecij[iDim] = rvecij[iDim]/rij;
      } // for: iDim

//////// Pair Term Vij -----------------------------------------
      double Vij, dVij; // Value and its Derivative.
      sw.Phi[ijPair]->Calc(rij,Vij,dVij);

//Vij=0; dVij=0;

      // Energy [eV]
      atom.epot[i] += Vij/2.0;

      // Force [eV/A]
      atom.fx[i] += nvecij[0]*dVij;
      atom.fy[i] += nvecij[1]*dVij;
      atom.fz[i] += nvecij[2]*dVij;

      // Stress [eV/A3]
      for(int iDim=0;iDim<3;iDim++){
      for(int jDim=0;jDim<3;jDim++){
        atom.satom[i][iDim][jDim] += rvecij[iDim]*rvecij[jDim]*dVij/rij/2.0;
      } // for: jDim
      } // for: iDim

//// Trio Term Vjik --------------------------------------------
//*
      for(KLOOP){
        int k = GetNeighbor(book.alist[i][k1],iTrans);
        int kSpec = sw.AtomNum2Index[atom.anum[k]];
        int ikPair  = sw.PairType[iSpec][kSpec];
        int jikTrio = sw.TrioType[jSpec][iSpec][kSpec];
        if( ikPair==0){continue;} // Skip Non-Reactive Pair.
        if(jikTrio==0){continue;} // Skip Non-Reactive Trio.

        if(rij>sw.R[jikTrio]->GetCutoffRadius()){continue;}

        double rveck[3];       // Position of Replica of k-th Atom [A]
        double rvecik[3];      // rk-ri [A]
        double rik;            // |rik| [A]
        double nvecik[3];      // rik/|rik| [-]
        double rvecjk[3];      // rk-rj [A]
        double rjk;            // |rjk| [A]

        rvecik[0] = atom.Dx(i,k,iTrans[0],iTrans[1],iTrans[2])/ang; // rik := rk - ri.
        rvecik[1] = atom.Dy(i,k,iTrans[0],iTrans[1],iTrans[2])/ang;
        rvecik[2] = atom.Dz(i,k,iTrans[0],iTrans[1],iTrans[2])/ang;
        rik = Norm(rvecik);
        if(rik<eps){continue;} // Skip Too Close (or Identical) Pair.
        if(rik>sw.R[jikTrio]->GetCutoffRadius()){continue;}

        for(int iDim=0;iDim<3;iDim++){
          rveck[iDim]  = rveci[iDim] + rvecik[iDim];
          nvecik[iDim] = rvecik[iDim]/rik;
          rvecjk[iDim] = rveck[iDim] - rvecj[iDim];
        } // for: iDim
        rjk = Norm(rvecjk);
        if(rjk<eps){continue;} // Skip Too Close (= Identical) Pair.

        double Rij,  dRij;  // Radial  Function in V3.
        double Rik,  dRik;  // Radial  Function in V3.
        double Pjik, dPjik = 0; // Angular Function in V3.
        double costh = innerproduct(nvecij,nvecik);
        sw.R[jikTrio]->Calc(rij,Rij,dRij);
        sw.R[jikTrio]->Calc(rik,Rik,dRik);
        sw.P[jikTrio]->Calc(costh,Pjik,dPjik);
        double Vjik = Rij*Rik*Pjik;

        double dVjik_dri[3], dVjik_drj[3], dVjik_drk[3];
        double  dcos_dri[3],  dcos_drj[3],  dcos_drk[3];
        DerivCOS(rvecij,rvecik,dcos_dri,dcos_drj,dcos_drk);
        for(int iDim=0;iDim<3;iDim++){
          dVjik_dri[iDim] = - dRij* Rik* Pjik*nvecij[iDim]
                            -  Rij*dRik* Pjik*nvecik[iDim]
                            +  Rij* Rik*dPjik*dcos_dri[iDim];
          dVjik_drj[iDim] = + dRij* Rik* Pjik*nvecij[iDim]
                            +  Rij* Rik*dPjik*dcos_drj[iDim];
          dVjik_drk[iDim] = +  Rij*dRik* Pjik*nvecik[iDim]
                            +  Rij* Rik*dPjik*dcos_drk[iDim];
        } // for: iDim

        // Energy [eV]
        atom.epot[i] += Vjik/2.0;

        // Force [eV/A]
        atom.fx[i] -= dVjik_dri[0]/2.0;
        atom.fy[i] -= dVjik_dri[1]/2.0;
        atom.fz[i] -= dVjik_dri[2]/2.0;

        atom.fx[j] -= dVjik_drj[0]/2.0;
        atom.fy[j] -= dVjik_drj[1]/2.0;
        atom.fz[j] -= dVjik_drj[2]/2.0;

        atom.fx[k] -= dVjik_drk[0]/2.0;
        atom.fy[k] -= dVjik_drk[1]/2.0;
        atom.fz[k] -= dVjik_drk[2]/2.0;

        // Stress [eV/A3]
        for(int iDim=0;iDim<3;iDim++){
        for(int jDim=0;jDim<3;jDim++){
          //atom.satom[i][iDim][jDim] += rveci[iDim]*dVjik_dri[jDim]/2.0;
          //atom.satom[j][iDim][jDim] += rvecj[iDim]*dVjik_drj[jDim]/2.0;
          //atom.satom[k][iDim][jDim] += rveck[iDim]*dVjik_drk[jDim]/2.0;
	  atom.satom[i][iDim][jDim] += 
	    + dRij* Rik* Pjik*nvecij[iDim] * nvecij[jDim] * rij / 2.0
	    +  Rij*dRik* Pjik*nvecik[iDim] * nvecik[jDim] * rik / 2.0
	    +  Rij* Rik*dPjik / 2.0
	    *((-rvecjk[iDim]*rvecjk[jDim]
	       +rvecik[iDim]*rvecik[jDim]
	       +rvecij[iDim]*rvecij[jDim])/rik/rij
	      -costh*( nvecik[iDim]*nvecik[jDim]
		      +nvecij[iDim]*nvecij[jDim]));
        } // for: jDim
        } // for: iDim

      } // KLOOP
//*/

    } // JLOOP

  } // ILOOP


//// Ending: Cenvert Units eV->J, A->m
  for(ILOOP){
    atom.epot[i] *= eV;
    atom.epotsum += atom.epot[i];
    atom.fx[i] *= eV/ang;
    atom.fy[i] *= eV/ang;
    atom.fz[i] *= eV/ang;
  } // ILOOP

  resetmat(cell.dmat);
  for(int iDim=0;iDim<3;iDim++){
  for(int jDim=0;jDim<3;jDim++){
    cell.dmat[iDim][jDim] = 0.0;
    for(ILOOP){
    	atom.satom[i][iDim][jDim] *= eV;
    	cell.dmat[iDim][jDim] -= atom.satom[i][iDim][jDim];
     atom.satom[i][iDim][jDim] *= (double)atom.natom / cell.volume;
    } // ILOOP
  } // for: jDim
  } // for: iDim

}





//// UTILITIES /////////////////////////////////////////////////

int GetNeighbor( //=============================================
  int *bk,    // IN:  Bookkeeping Information
  int *iTrans // OUT: Translation Vector in unit of Cell.
){
  for(int iDim=0;iDim<3;iDim++){iTrans[iDim] = bk[iDim+1];}
  return(bk[0]);
}

void DerivCOS( //===============================================
  double *rvecij,  // IN: Vector rij := rj - ri.
  double *rvecik,  // IN: Vector rik := rk - ri.
  double *dcos_dri // OUT: Vector dcos(theta)/dri.
){
  double rij, rik;
  double nij[3], nik[3];

  rij = Norm(rvecij);
  rik = Norm(rvecik);

  for(int iDim=0;iDim<3;iDim++){
    nij[iDim] = rvecij[iDim]/rij;
    nik[iDim] = rvecik[iDim]/rik;
  }

  double costh = innerproduct(nij,nik);
  for(int iDim=0;iDim<3;iDim++){
    double term[2];
    term[0] = nik[iDim] - costh*nij[iDim];
    term[1] = nij[iDim] - costh*nik[iDim];
    dcos_dri[iDim] = - term[0]/rij - term[1]/rik;
  }
}

void DerivCOS( //===============================================
  double *rvecij, // IN: Vector rij := rj - ri.
  double *rvecik, // IN: Vector rik := rk - ri.
  double *dcos_dri, // OUT: Vector dcos(theta)/dri.
  double *dcos_drj, // OUT: Vector dcos(theta)/drj.
  double *dcos_drk  // OUT: Vector dcos(theta)/drk.
){

  double rij, rik;
  double direcij[3], direcik[3];

  rij = Norm(rvecij);
  rik = Norm(rvecik);

  for(int iDim=0;iDim<3;iDim++){
    direcij[iDim] = rvecij[iDim]/rij;
    direcik[iDim] = rvecik[iDim]/rik;
  }

  for(int iDim=0;iDim<3;iDim++){
    dcos_dri[iDim] = 0.0;
    dcos_drj[iDim] = 0.0;
    dcos_drk[iDim] = 0.0;
    for(int jDim=0;jDim<3;jDim++){
      double tensorij = delta(iDim,jDim) - direcij[iDim]*direcij[jDim];
      double tensorik = delta(iDim,jDim) - direcik[iDim]*direcik[jDim];

      dcos_drj[iDim] += tensorij*direcik[jDim]/rij;
      dcos_drk[iDim] += tensorik*direcij[jDim]/rik;
    }
    dcos_dri[iDim] = - dcos_drj[iDim] - dcos_drk[iDim];
  }
}

double Norm(double *vec){ // |vec|
  double n2 = Norm2(vec);
  return(sqrt(n2));
}

double Norm2(double *vec){
  return(vec[0]*vec[0]
        +vec[1]*vec[1]
        +vec[2]*vec[2]);
}

int delta(int i, int j){ // Kronecker's delta
  if(i==j){return(1);}
  else    {return(0);}
}






//////// Initialization ////////////////////////////////////////
void Initialize(void){
  atom.epotsum=0.0;
  for(ILOOP){
    atom.epot[i] = 0.0;
    atom.fx[i]   = 0.0;
    atom.fy[i]   = 0.0;
    atom.fz[i]   = 0.0;
    for(int iDim=0;iDim<3;iDim++){
    for(int jDim=0;jDim<3;jDim++){
      atom.satom[i][iDim][jDim] = 0.0;
    }
    }
  }
}

//void SetParam(void){
void SetParam(const char* arg){
  sw.DeleteArray();

  if ((strcmp(arg,"SiC_Vashishta") == 0)||
      (strcmp(arg,"GaN_Bere") == 0)||
      (strcmp(arg,"GaN_Aichoune") == 0)) {
  sw.nSpec = 2;
  sw.nPair = 3;
  sw.nTrio = 2;
  }
  if (strcmp(arg,"Si_Stillinger") == 0) {
    sw.nSpec = 1;
    sw.nPair = 1;
    sw.nTrio = 1;
  }
  sw.AllocateArray();
  sw.initialize=false;

  if (strcmp(arg,"SiC_Vashishta") == 0) {
  sw.AtomNum2Index[14] = 1; // Si = 1
  sw.AtomNum2Index[ 6] = 2; // C  = 2
  sw.PairType[1][1] = 1; // Si-Si
  sw.PairType[1][2] = 2; // Si-C
  sw.PairType[2][1] = 2; // C-Si
  sw.PairType[2][2] = 3; // C-C
  // Notice! Order is j-i-k.
  sw.TrioType[1][2][1] = 1; //  Si-C-Si
  sw.TrioType[2][1][2] = 2; //  C-Si-C
  } else if ((strcmp(arg,"GaN_Bere") == 0)||
	     (strcmp(arg,"GaN_Aichoune") == 0)) {
  sw.AtomNum2Index[31] = 1; // Ga = 1
  sw.AtomNum2Index[ 7] = 2; // N  = 2
  sw.PairType[1][1] = 1; // Ga-Ga
  sw.PairType[1][2] = 2; // Ga-N
  sw.PairType[2][1] = 2; // N-Ga
  sw.PairType[2][2] = 3; // N-N
  // Notice! Order is j-i-k.
  sw.TrioType[1][2][1] = 1; //  Ga-N-Ga
  sw.TrioType[2][1][2] = 2; //  N-Ga-N
  } else if (strcmp(arg,"Si_Stillinger") == 0) {
  sw.AtomNum2Index[14] = 1; // Si = 1
  sw.PairType[1][1] = 1; // Si-Si
  // Notice! Order is j-i-k.
  sw.TrioType[1][1][1] = 1; //  Si-Si-Si
  }

  if (strcmp(arg,"SiC_Vashishta") == 0) { // SiC_Vashishta
  printf("Stillinger-Weber potential: Vashishta for SiC\n");
  //// General Parameters
  double zSi    = +1.201; // [e]
  double zC     = -zSi;   // [e]
  double diel   = 14.40;  // [AeV/e2]
  double lambda = 5.0;    // [A]
  double xi     = 3.0;    // [A]
  double rCut   = 7.35;   // [A]
  rcut = rCut;
  //// Pair Function =============================================
  sw.Phi[0] = new PairFunction();
  { //// Si-Si -----------------
  double param[7];
  string functype = "vashpair_shift";
  param[0] = 23.67291;     // H_ij
  param[1] = zSi*zSi*diel; // Z_i x Z_j
  param[2] = 2.1636*diel;  // D_ij
  param[3] = 0.00000;      // W_ij
  param[4] = 7.00000;      // eta_ij
  param[5] = lambda;       // lambda
  param[6] = xi;           // xi
  sw.Phi[1] = new PairFunction(functype,param,rCut);
  }
  { //// Si-C ------------------
  double param[7];
  string functype = "vashpair_shift";
  param[0] = 447.09026;   // H_ij
  param[1] = zSi*zC*diel; // Z_i x Z_j
  param[2] = 1.0818*diel; // D_ij
  param[3] = 61.46940;    // W_ij
  param[4] = 9.00000;     // eta_ij
  param[5] = lambda;      // lambda
  param[6] = xi;          // xi
  sw.Phi[2] = new PairFunction(functype,param,rCut);
  }
  { //// C-C -------------------
  double param[7];
  string functype = "vashpair_shift";
  param[0] = 471.74538;  // H_ij
  param[1] = zC*zC*diel; // Z_i x Z_j
  param[2] = 0.0*diel;   // D_ij
  param[3] = 0.00000;    // W_ij
  param[4] = 7.00000;    // eta_ij
  param[5] = lambda;     // lambda
  param[6] = xi;         // xi
  sw.Phi[3] = new PairFunction(functype,param,rCut);
  }
  //// Trio Function =============================================
  sw.R[0] = new PairFunction();
  { //// Si-C-Si ---------------
  double param[2];
  string functype = "vashtrior";
  param[0] = 1.00; // gamma_jik
  param[1] = 2.90; // r0_jik
  sw.R[1] = new PairFunction(functype,param,param[1]);
  }
  { //// C-Si-C ----------------
  double param[2];
  string functype = "vashtrior";
  param[0] = 1.00; // gamma_jik
  param[1] = 2.90; // r0_jik
  sw.R[2] = new PairFunction(functype,param,param[1]);
  }
  sw.P[0] = new PairFunction();
  { //// Si-C-Si ---------------
  double param[3];
  string functype = "vashtriop";
  param[0] =  9.003;   // B_jik
  param[1] =  5.0;     // C_jik
  param[2] = -1.0/3.0; // cos(theta0_jik)
  sw.P[1] = new PairFunction(functype,param,2.0);
  }
  { //// C-Si-C ----------------
  double param[3];
  string functype = "vashtriop";
  param[0] =  9.003;   // B_jik
  param[1] =  5.0;     // C_jik
  param[2] = -1.0/3.0; // cos(theta0_jik)
  sw.P[2] = new PairFunction(functype,param,2.0);
  }

  } else if ((strcmp(arg,"GaN_Bere") == 0)||
	     (strcmp(arg,"GaN_Aichoune") == 0)) { // GaN_Bere or GaN_Aichoune
  //// General Parameters
  //double rCut   = 7.35;   // [A]
  double rCut   = 3.35;   // [A]
  rcut = rCut;
  //// Pair Function =============================================
  sw.Phi[0] = new PairFunction();
  { //// Ga-Ga -----------------
  double param[7];
  string functype = "swpair";
  if (strcmp(arg,"GaN_Bere") == 0) {
  printf("Stillinger-Weber potential: Bere for SiC\n");
  param[0] = 1.200; // epsilon
  param[1] = 7.917; // A
  param[2] = 0.720; // B
  param[3] = 1.600; // a
  param[4] = 2.100; // sigma
  } else {
  printf("Stillinger-Weber potential: Aichoune for SiC\n");
  param[0] = 0.655; // epsilon
  param[1] = 7.917; // A
  param[2] = 0.720; // B
  param[3] = 1.800; // a
  param[4] = 2.038; // sigma
  }
  //sw.Phi[1] = new PairFunction(functype,param,rCut);
  sw.Phi[1] = new PairFunction(functype,param,param[3]*param[4]);
  }
  { //// Ga-N ------------------
  double param[7];
  string functype = "swpair";
  param[0] = 2.170; // epsilon
  param[1] = 7.917; // A
  param[2] = 0.720; // B
  param[3] = 1.800; // a
  param[4] = 1.695; // sigma
  //sw.Phi[2] = new PairFunction(functype,param,rCut);
  sw.Phi[2] = new PairFunction(functype,param,param[3]*param[4]);
  }
  { //// N-N -------------------
  double param[7];
  string functype = "swpair";
  if (strcmp(arg,"GaN_Bere") == 0) {
  param[0] = 1.200; // epsilon
  param[1] = 7.917; // A
  param[2] = 0.720; // B
  param[3] = 1.800; // a
  param[4] = 1.300; // sigma
  } else {
  param[0] = 0.655; // epsilon
  param[1] = 7.917; // A
  param[2] = 0.720; // B
  param[3] = 1.800; // a
  param[4] = 1.302; // sigma
  }
  //sw.Phi[3] = new PairFunction(functype,param,rCut);
  sw.Phi[3] = new PairFunction(functype,param,param[3]*param[4]);
  }
  //// Trio Function =============================================
  sw.R[0] = new PairFunction();
  { //// Ga-N-Ga ---------------
  double param[2];
  string functype = "swtrior";
  param[0] = 1.000; // gamma
  param[1] = 1.800; // a
  param[2] = 1.695; // sigma
  sw.R[1] = new PairFunction(functype,param,param[1]*param[2]);
  }
  { //// N-Ga-N ----------------
  double param[2];
  string functype = "swtrior";
  param[0] = 1.000; // gamma
  param[1] = 1.800; // a
  param[2] = 1.695; // sigma
  sw.R[2] = new PairFunction(functype,param,param[1]*param[2]);
  }
  sw.P[0] = new PairFunction();
  { //// Ga-N-Ga ---------------
  double param[3];
  string functype = "swtriop";
  param[0] = 2.17; // epsilon
  param[1] = 32.5; // lambda
  sw.P[1] = new PairFunction(functype,param,2.0);
  }
  { //// N-Ga-N ----------------
  double param[3];
  string functype = "swtriop";
  param[0] = 2.17; // epsilon
  param[1] = 32.5; // lambda
  sw.P[2] = new PairFunction(functype,param,2.0);
  }
  }

  else if (strcmp(arg,"Si_Stillinger") == 0) { // Si_Stillinger
  printf("Stillinger-Weber potential: Original for Si\n");
  //// General Parameters
  double rCut   = 2.0951 * 1.8;   // 3.77118 [A]
  rcut = rCut;
  //// Pair Function =============================================
  sw.Phi[0] = new PairFunction();
  { //// Si-Si -----------------
  double param[7];
  string functype = "swpair";
  param[0] = 2.1685; // epsilon [eV] (=50kcal/mol)
  param[1] = 7.049556; // A
  param[2] = 0.602224; // B
  param[3] = 1.800;  // a
  param[4] = 2.0951; // sigma [A]
  sw.Phi[1] = new PairFunction(functype,param,param[3]*param[4]);
  }
  //// Trio Function =============================================
  sw.R[0] = new PairFunction();
  { //// Si-Si-Si ---------------
  double param[2];
  string functype = "swtrior";
  param[0] = 1.200; // gamma
  param[1] = 1.800; // a
  param[2] = 2.0951; // sigma
  sw.R[1] = new PairFunction(functype,param,param[1]*param[2]);
  }
  sw.P[0] = new PairFunction();
  { //// Si-Si-Si ---------------
  double param[3];
  string functype = "swtriop";
  param[0] = 2.1685; // epsilon
  param[1] = 21.0; // lambda
  sw.P[1] = new PairFunction(functype,param,2.0);
  }

  }

/*
  sw.Phi[1]->Plot(1,8,1000,"Phi_Si-Si");
  sw.Phi[2]->Plot(1,8,1000,"Phi_Si-C");
  sw.Phi[3]->Plot(1,8,1000,"Phi_C-C");

  sw.R[1]->Plot(1,3,1000,"R_Si-C-Si");
  sw.R[2]->Plot(1,3,1000,"R_C-Si-C");

  sw.P[1]->Plot(-1,1,1000,"P_Si-C-Si");
  sw.P[2]->Plot(-1,1,1000,"P_C-Si-C");
*/

}


// For Debugging
void ShowVector(double *v){
    cout<<"("<<v[0]<<","<<v[1]<<","<<v[2]<<")"<<endl;
}


