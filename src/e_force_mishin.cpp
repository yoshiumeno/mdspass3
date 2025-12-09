#include <string>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include "myheader.h"

double density(double);
double density_d(double);

double embedded(double);
double embedded_d(double);

double eam_pair(double);
double eam_pair_d(double);

void eammis_alloc();
void eammis_read(const char*);
void get_first_arg(std::string &line, std::string &arg1);
int remove_head_spaces(std::string &line);
double splint(double xx[], double yy[], int n, double xi);
void sortdi1(double a[], int n, int jun[]);
void subspl(double x[], double y[], int n, double h[], double sp[]);
int atom_number(char* at);
void remove_carriage_return(std::string& line);

void e_force_mishin(PotentialMode mode)
{


  double rr, rr2, drx, dry, drz, vp0, v0, rcut_eam;
  int j, ix, iy, iz;
  double ki, kif; 

  if (eammis.initialize) {
    std::cout<<"Initialize of the EAM (Mishin) potential.."<<std::endl;
    eammis_alloc();
    if (mode == NORMALPOTMODE) {
      strcpy(eammis.species, atom.asp[1]);
    } else {
      strcpy(eammis.species, "Ni");
    }
    eammis_read(eammis.species);
    eammis.initialize = false;
    std::cout<<"Done with EAM potential."<<std::endl;
  }
  
  // set cutoff radius (rcut and frc) only if mode is not combined (otherwise it could be set to a higher value by other potentials) or in combined mode when it is smaller then needed for eam
  if (mode == NORMALPOTMODE || (mode == COMBINEDPOTMODE && rcut < eammis.pairmax*ang))
    rcut = eammis.pairmax*ang;
  
  /* set cut-off values used only inside this function (only for eam) */
  rcut_eam = eammis.pairmax*ang;
  double rcut_eam2 = rcut_eam * rcut_eam;
  
  if (book.frc < rcut) {
    book.frc = rcut * 1.15; printf("### WARNING: frc is corrected (EAM Mishin)\n");}
  
  //  double rcut2 = rcut * rcut;
  
  // Initialize charge densities at atomic sites
  for (int i=1; i<=atom.natom; i++) {
    if (strcmp(atom.asp[i], eammis.species) != 0) continue;
    eammis.rhob[i] = 0.0;
  }
  
  /* reset all values in one potential mode */
  if (mode == NORMALPOTMODE)
    {
      // Initialize stress-related variables
      cell.virx = 0.0; cell.viry = 0.0; cell.virz = 0.0;
      for (int i=0; i<3; i++) { for (int j=0; j<3; j++) { 
	  cell.dmat[i][j] = 0.0;
	  for (int ii=1; ii<=atom.natom; ii++) {
	    if (eammis.atom_number == atom.anum[ii]) atom.satom[ii][i][j] = 0.0; }  } }
      // Initialize Force and E_pot
      for (int i=1; i<=atom.natom; i++) { if (eammis.atom_number == atom.anum[i]) {
	  atom.fx[i] = 0.0; atom.fy[i] = 0.0; atom.fz[i] = 0.0;  atom.epot[i] = 0.0; } }
      
      // total potential energy
      atom.epotsum=0.0;
    }

  // First loop: Chg dns calc
  for (int i=1; i<=atom.natom; i++) {
    if (eammis.atom_number != atom.anum[i]) continue;
    if (book.alistnum[i]>0) {
      for (int k=1; k<=book.alistnum[i]; k++) {
	j = book.alist[i][k][0];
	
	if (eammis.atom_number != atom.anum[j]) continue;
	
	if (j>=i) {
	  ix = book.alist[i][k][1]; iy = book.alist[i][k][2]; iz = book.alist[i][k][3];
	  rr2 = atom.Dist2(i,j,ix,iy,iz);
	  if (rr2 < rcut_eam2) {
	    rr = sqrt(rr2);
	    double rho_ij = density(rr);
	    //    rhob[i] = rhob[i] + geamden(rr/ang, geamsp(atom.asp[j]));
	    eammis.rhob[i]=eammis.rhob[i]+rho_ij;
	    if (j!=i) {
	      //      rhob[j] = rhob[j] + geamden(rr/ang, geamsp(atom.asp[i])); }
	      eammis.rhob[j]=eammis.rhob[j]+rho_ij;}
	  }
	} //j>=i
      }
    }
  }

  for (int i=1; i<=atom.natom; i++) {
    
    if (eammis.atom_number != atom.anum[i]) continue;
    
    //printf("i: %d, rho: %10.5f, embedd: %10.5f\n", i, eammis.rhob[i], embedded(eammis.rhob[i])/eV);
    atom.epot[i]=embedded(eammis.rhob[i]);
  }
  
  // Second loop: Force calc
  for (int i=1; i<=atom.natom; i++) {

    if (eammis.atom_number != atom.anum[i]) continue;
    
    if (book.alistnum[i]>0) {
      for (int k=1; k<=book.alistnum[i]; k++) {
	j  = book.alist[i][k][0];
	
	if (eammis.atom_number != atom.anum[j]) continue;
	
	if (j>=i) {
	  if ((atom.QC==0)||(atom.repatom[i]==1)||(atom.repatom[j]==1)) {
	    ix = book.alist[i][k][1]; iy = book.alist[i][k][2]; iz = book.alist[i][k][3];
	    rr2 = atom.Dist2(i,j,ix,iy,iz);
	    //	if ((rr2 < rcut2)&&(rr2>1.0e-30)) {
	    if (rr2 < rcut_eam2) {
	      rr = sqrt(rr2);
	      drx=-atom.Dx(i,j,ix,iy,iz); dry=-atom.Dy(i,j,ix,iy,iz); drz=-atom.Dz(i,j,ix,iy,iz);

	      double vij=eam_pair(rr);
	      double dv=eam_pair_d(rr);

	      double df_i=embedded_d(eammis.rhob[i]);
              double df_j=embedded_d(eammis.rhob[j]);
	 
	     double dr=density_d(rr);
	      double f=-(df_i*dr+df_j*dr+dv)*eV/ang;
	      double ki=f/rr;
	      if (j==i) ki=ki/2.0;
	      double kif=ki/2.0;

	      if ((atom.QC==0)||(atom.repatom[i]==1)) {
		atom.fx[i] = atom.fx[i] + ki*drx;
		atom.fy[i] = atom.fy[i] + ki*dry;
		atom.fz[i] = atom.fz[i] + ki*drz;
		atom.epot[i]+=vij/2.0;
	      }
	      if ((atom.QC==0)||(atom.repatom[j]==1)) {
		atom.fx[j] = atom.fx[j] - ki*drx;
		atom.fy[j] = atom.fy[j] - ki*dry;
		atom.fz[j] = atom.fz[j] - ki*drz;
		if (j != i) {
		  atom.epot[j]+=vij/2.0;}
	      }
	      /* This causes a different result from non-QC.
	      if ((atom.QC)&&(atom.repatom[i]==0)) {
		int iel = atom.elem_id[i];
		for (int iv=1; iv<=4; iv++) {
		  atom.fx[ atom.elem_v[iel][iv] ] = atom.fx[ atom.elem_v[iel][iv] ] + ki*drx/4;
		  atom.fy[ atom.elem_v[iel][iv] ] = atom.fy[ atom.elem_v[iel][iv] ] + ki*dry/4;
		  atom.fz[ atom.elem_v[iel][iv] ] = atom.fz[ atom.elem_v[iel][iv] ] + ki*drz/4;
		}
	      }
	      if ((atom.QC)&&(atom.repatom[j]==0)) {
		int iel = atom.elem_id[j];
		for (int iv=1; iv<=4; iv++) {
		  atom.fx[ atom.elem_v[iel][iv] ] = atom.fx[ atom.elem_v[iel][iv] ] - ki*drx/4;
		  atom.fy[ atom.elem_v[iel][iv] ] = atom.fy[ atom.elem_v[iel][iv] ] - ki*dry/4;
		  atom.fz[ atom.elem_v[iel][iv] ] = atom.fz[ atom.elem_v[iel][iv] ] - ki*drz/4;
		}
	      }
	      */
	      if (atom.QC==0) {
		atom.satom[i][0][0] -= kif*drx*drx;
		atom.satom[i][0][1] -= kif*dry*drx;
		atom.satom[i][1][1] -= kif*dry*dry;
		atom.satom[i][0][2] -= kif*drz*drx;
		atom.satom[i][1][2] -= kif*drz*dry;
		atom.satom[i][2][2] -= kif*drz*drz;
		atom.satom[j][0][0] -= kif*drx*drx;
		atom.satom[j][0][1] -= kif*dry*drx;
		atom.satom[j][1][1] -= kif*dry*dry;
		atom.satom[j][0][2] -= kif*drz*drx;
		atom.satom[j][1][2] -= kif*drz*dry;
		atom.satom[j][2][2] -= kif*drz*drz;
	      }
	    }
	  }
	} //j>=i
      }
    }
    //    printf("FAST %d %20.15e %20.15e\n",i,atom.fx[i],rr2);
  }
  //  if (istep==0) {exit(0);}
  
  
  if (atom.QC==0) {
    for (int i=1; i<=atom.natom; i++) {
      if (eammis.atom_number != atom.anum[i]) continue;
      
      atom.satom[i][1][0] = atom.satom[i][0][1];
      atom.satom[i][2][0] = atom.satom[i][0][2];
      atom.satom[i][2][1] = atom.satom[i][1][2];}
    for (int ii=1; ii<=atom.natom; ii++) {
      
      if (eammis.atom_number != atom.anum[ii]) continue;
      
      for (int i=0; i<3; i++) { for (int j=0; j<3; j++) { 
	  cell.dmat[i][j] = cell.dmat[i][j] - atom.satom[ii][i][j];
	} } }
    cell.virx = cell.dmat[0][0]; cell.viry = cell.dmat[1][1]; cell.virz = cell.dmat[2][2];
    cell.volume = cell.Getvolume();
    for (int ii=1; ii<=atom.natom; ii++) {
      
      if (eammis.atom_number != atom.anum[ii]) continue;
      
      for (int i=0; i<3; i++) { for (int j=0; j<3; j++) {
	  atom.satom[ii][i][j] = atom.satom[ii][i][j] * (double)eammis.natoms / cell.volume;
	} } }
  }
  
  // calculate potential energy (in combined mode in will be calculated after all potentials)
  if (mode != COMBINEDPOTMODE)
    {
      atom.epotsum=0.0;
      
      for (int i=1; i<=atom.natom; i++) {
	if (strcmp(atom.asp[i], eammis.species) != 0) continue;
	
	atom.epotsum += atom.epot[i]; }
    }
} // e_force_mishin

void eammis_alloc()
{
  if (eammis.rhob)   { delete[] eammis.rhob;  eammis.rhob   = NULL; }
  if (eammis.pair)   { delete[] eammis.pair;  eammis.pair   = NULL; }
  if (eammis.paird)  { delete[] eammis.paird; eammis.paird  = NULL; }
  if (eammis.den)    { delete[] eammis.den;   eammis.den    = NULL; }
  if (eammis.dend)   { delete[] eammis.dend;  eammis.dend   = NULL; }
  if (eammis.embed)  { delete[] eammis.embed; eammis.embed  = NULL; }
  if (eammis.embedd) { delete[] eammis.embedd;eammis.embedd = NULL; }

  if (eammis.aden)   { delete[] eammis.aden;  eammis.aden   = NULL; }
  if (eammis.apair)  { delete[] eammis.apair; eammis.apair  = NULL; }
  if (eammis.aembed) { delete[] eammis.aembed;eammis.aembed = NULL; }
 
  if (eammis.bden)   { delete[] eammis.bden;  eammis.bden   = NULL; }
  if (eammis.bpair)  { delete[] eammis.bpair; eammis.bpair  = NULL; }
  if (eammis.bembed) { delete[] eammis.bembed;eammis.bembed = NULL; }

  if (eammis.cden)   { delete[] eammis.cden;  eammis.cden   = NULL; }
  if (eammis.cpair)  { delete[] eammis.cpair; eammis.cpair  = NULL; }
  if (eammis.cembed) { delete[] eammis.cembed;eammis.cembed = NULL; }

  if (eammis.dden)   { delete[] eammis.dden;  eammis.dden   = NULL; }
  if (eammis.dpair)  { delete[] eammis.dpair; eammis.dpair  = NULL; }
  if (eammis.dembed) { delete[] eammis.dembed;eammis.dembed = NULL; }

  eammis.rhob   = new double[atom.natom+1];
  eammis.pair   = new double[eammis.mesh+1];
  eammis.paird  = new double[eammis.mesh+1];
  eammis.den    = new double[eammis.mesh+1];
  eammis.dend   = new double[eammis.mesh+1];
  eammis.embed  = new double[eammis.mesh+1];
  eammis.embedd = new double[eammis.mesh+1];

  eammis.aden   = new double[eammis.mesh];
  eammis.apair  = new double[eammis.mesh];
  eammis.aembed = new double[eammis.mesh];

  eammis.bden   = new double[eammis.mesh];
  eammis.bpair  = new double[eammis.mesh];
  eammis.bembed = new double[eammis.mesh];

  eammis.cden   = new double[eammis.mesh];
  eammis.cpair  = new double[eammis.mesh];
  eammis.cembed = new double[eammis.mesh];

  eammis.dden   = new double[eammis.mesh];
  eammis.dpair  = new double[eammis.mesh];
  eammis.dembed = new double[eammis.mesh];
}

void eammis_read(const char* atom_type)
{
  int mesh_pair = 0, mesh_den = 0, mesh_embed = 0;
  FILE *fp;
  char fname_p[80] = "aaa", fnamf_p[80] = "aaa";
  char fname_d[80] = "aaa", fnamf_d[80] = "aaa";
  char fname_e[80] = "aaa", fnamf_e[80] = "aaa";
  std::string line, arg1, arg2;
  char larg1[40], larg2[40]; double xx, yy; int ix, iy;
  bool fine_pair = false, fine_den = false, fine_embed = false;

  std::ifstream fin(fnamf_p); fin.close(); fin.clear(); // "clear" needed for Win

#if defined __linux__ || defined __APPLE__
    strcpy(fname_p,"pot/EAMPARAM.");
    strcpy(fname_d,"pot/EAMPARAM.");
    strcpy(fname_e,"pot/EAMPARAM.");
    strcpy(fnamf_p,"pot/EAMPARAM.");
    strcpy(fnamf_d,"pot/EAMPARAM.");
    strcpy(fnamf_e,"pot/EAMPARAM.");
#else
    strcpy(fname_p,"pot\\EAMPARAM.");
    strcpy(fname_d,"pot\\EAMPARAM.");
    strcpy(fname_e,"pot\\EAMPARAM.");
    strcpy(fnamf_p,"pot\\EAMPARAM.");
    strcpy(fnamf_d,"pot\\EAMPARAM.");
    strcpy(fnamf_e,"pot\\EAMPARAM.");
#endif
  if ((strcmp(atom_type,"Cu")==0)) {
    strcat(fname_p,"PRB63-224106.cu.pair");
    strcat(fname_d,"PRB63-224106.cu.den");
    strcat(fname_e,"PRB63-224106.cu.embed");
    strcat(fnamf_p,"PRB63-224106.cu.pair.fine");
    strcat(fnamf_d,"PRB63-224106.cu.den.fine");
    strcat(fnamf_e,"PRB63-224106.cu.embed.fine");
  } else if ((strcmp(atom_type,"Ni")==0)) {
    strcat(fname_p,"PRB59-3393.ni.pair");
    strcat(fname_d,"PRB59-3393.ni.den");
    strcat(fname_e,"PRB59-3393.ni.embed");
    strcat(fnamf_p,"PRB59-3393.ni.pair.fine");
    strcat(fnamf_d,"PRB59-3393.ni.den.fine");
    strcat(fnamf_e,"PRB59-3393.ni.embed.fine");
  } else if ((strcmp(atom_type,"Al")==0)) {
    strcat(fname_p,"PRB59-3393.al.pair");
    strcat(fname_d,"PRB59-3393.al.den");
    strcat(fname_e,"PRB59-3393.al.embed");
    strcat(fnamf_p,"PRB59-3393.al.pair.fine");
    strcat(fnamf_d,"PRB59-3393.al.den.fine");
    strcat(fnamf_e,"PRB59-3393.al.embed.fine");
  } else if ((strcmp(atom_type,"Au")==0)) {
    strcat(fname_p,"JCP123-204719.au.pair");
    strcat(fname_d,"JCP123-204719.au.den");
    strcat(fname_e,"JCP123-204719.au.embed");
    strcat(fnamf_p,"JCP123-204719.au.pair.fine");
    strcat(fnamf_d,"JCP123-204719.au.den.fine");
    strcat(fnamf_e,"JCP123-204719.au.embed.fine");
  } else { printf("EAM Mishin is only for Cu, Ni, Al and Au\n");
    printf("%s is not supported\n",atom_type); return;
  }

  printf("EAM Mishin for %s has been initialized\n", atom_type);

  // If fine-mesh data exist, they are used.
  fin.open(fnamf_p);
  if (fin) { fine_pair = true; ix = 0;
    while (getline(fin,line)) { remove_carriage_return(line);
      get_first_arg(line,arg1); get_first_arg(line, arg2);
      strcpy (larg1, arg1.c_str()); strcpy (larg2, arg2.c_str());
      if (ix == 0) { eammis.pairmin = atof(larg1); eammis.pairmax = atof(larg2);
      } else { eammis.pair[ix] = atof(larg1); eammis.paird[ix] = atof(larg2); } 
      ix++; } }
  fin.close(); fin.clear();
  fin.open(fnamf_d);
  if (fin) { fine_den = true; ix = 0;
    while (getline(fin,line)) { remove_carriage_return(line);
      get_first_arg(line,arg1); get_first_arg(line, arg2);
      strcpy (larg1, arg1.c_str()); strcpy (larg2, arg2.c_str());
      if (ix == 0) { eammis.denmin = atof(larg1); eammis.denmax = atof(larg2);
      } else { eammis.den[ix] = atof(larg1); eammis.dend[ix] = atof(larg2); } 
      ix++; } }
  fin.close(); fin.clear();
  fin.open(fnamf_e);
  if (fin) { fine_embed = true; ix = 0;
    while (getline(fin,line)) { remove_carriage_return(line);
      get_first_arg(line,arg1); get_first_arg(line, arg2);
      strcpy (larg1, arg1.c_str()); strcpy (larg2, arg2.c_str());
      if (ix == 0) { eammis.embedmin = atof(larg1); eammis.embedmax = atof(larg2);
      } else { eammis.embed[ix] = atof(larg1); eammis.embedd[ix] = atof(larg2); } 
      ix++; } }
  fin.close(); fin.clear();

  //If fine-mesh data do not exist, they are created.
  if (!fine_pair) {
  fin.open(fname_p);
  if (!fin) { printf("### FILE [%s] does not exist ###\n",fname_p); return; }
  while (getline(fin,line)) { remove_carriage_return(line);
    remove_head_spaces(line);
    if (line.at(0)!='#') {  mesh_pair++; } }
  eammis.pairx = new double[mesh_pair];
  eammis.pairy = new double[mesh_pair];
  //eammis.pairc = new double*[mesh_pair];
  //for (int i=0; i<=mesh_pair; i++) { eammis.pairc[i] = new double[4]; }
  fin.clear(); fin.seekg(0); ix = 0;
  while (getline(fin,line)) { remove_carriage_return(line);
    remove_head_spaces(line);
    if (line.at(0)!='#') {
      get_first_arg(line,arg1); get_first_arg(line, arg2);
      strcpy (larg1, arg1.c_str()); strcpy (larg2, arg2.c_str());
      if (ix==0) { eammis.pairmin=atof(larg1); }
      eammis.pairx[ix] = atof(larg1); eammis.pairy[ix] = atof(larg2); ix++; }
  }
  eammis.pairmax=atof(larg1); fin.close(); fin.clear();
  }

  if (!fine_den) {
  fin.open(fname_d);
  if (!fin) { printf("### FILE [%s] does not exist ###\n",fname_d); return; }
  while (getline(fin,line)) { remove_carriage_return(line);
    remove_head_spaces(line);
    if (line.at(0)!='#') {  mesh_den++; } }
  eammis.denx = new double[mesh_den];
  eammis.deny = new double[mesh_den];
  //eammis.denc = new double*[mesh_den];
  //for (int i=0; i<=mesh_den; i++) { eammis.denc[i] = new double[4]; }
  fin.clear(); fin.seekg(0); ix = 0;
  while (getline(fin,line)) { remove_carriage_return(line);
    remove_head_spaces(line);
    if (line.at(0)!='#') {
      get_first_arg(line,arg1); get_first_arg(line, arg2);
      strcpy (larg1, arg1.c_str()); strcpy (larg2, arg2.c_str());
      if (ix==0) { eammis.denmin=atof(larg1); }
      eammis.denx[ix] = atof(larg1); eammis.deny[ix] = atof(larg2); ix++; }
  }
  eammis.denmax=atof(larg1); fin.close(); fin.clear();
  }

  if (!fine_embed) {
  fin.open(fname_e);
  if (!fin) { printf("### FILE [%s] does not exist ###\n",fname_e); return; }
  while (getline(fin,line)) { remove_carriage_return(line);
    remove_head_spaces(line);
    if (line.at(0)!='#') {  mesh_embed++; } }
  eammis.embedx = new double[mesh_embed];
  eammis.embedy = new double[mesh_embed];
  //eammis.embedc = new double*[mesh_embed];
  //for (int i=0; i<=mesh_embed; i++) { eammis.embedc[i] = new double[4]; }
  fin.clear(); fin.seekg(0); ix = 0;
  while (getline(fin,line)) { remove_carriage_return(line);
    remove_head_spaces(line);
    if (line.at(0)!='#') {
      get_first_arg(line,arg1); get_first_arg(line, arg2);
      strcpy (larg1, arg1.c_str()); strcpy (larg2, arg2.c_str());
      if (ix==0) { eammis.embedmin=atof(larg1); }
      eammis.embedx[ix] = atof(larg1); eammis.embedy[ix] = atof(larg2); ix++; }
  }
  eammis.embedmax=atof(larg1); fin.close(); fin.clear();
  }

  double d, x;
  if (!fine_pair) {
  d=(eammis.pairmax-eammis.pairmin)/(double)(eammis.mesh-1);
  for (int i=0; i<=eammis.mesh-1; i++) {
    x=eammis.pairmin+d*(double)i;
    eammis.pair[i]=splint(eammis.pairx, eammis.pairy, mesh_pair, x); }
  for (int i=0; i<=eammis.mesh-1; i++) {
    if (i==0) {eammis.paird[i]=(eammis.pair[i+1]-eammis.pair[i])/d;
    } else if (i==eammis.mesh-1) {eammis.paird[i]=(eammis.pair[i]-eammis.pair[i-1])/d;
    } else { eammis.paird[i]=(eammis.pair[i+1]-eammis.pair[i-1])/d/2; } }
  delete[] eammis.pairx; delete[] eammis.pairy;
  fp = fopen(fnamf_p, "w");
  fprintf(fp,"%20.15e %20.15e\n",eammis.pairmin,eammis.pairmax);
  for (int i=0; i<=eammis.mesh-1; i++) {
    fprintf(fp,"%20.15e %20.15e\n",eammis.pair[i],eammis.paird[i]); }
  fclose(fp);
  }
  if (!fine_den) {
  d=(eammis.denmax-eammis.denmin)/(double)(eammis.mesh-1);
  for (int i=0; i<=eammis.mesh-1; i++) {
    x=eammis.denmin+d*(double)i;
    eammis.den[i]=splint(eammis.denx, eammis.deny, mesh_den, x); }
  for (int i=0; i<=eammis.mesh-1; i++) {
    if (i==0) {eammis.dend[i]=(eammis.den[i+1]-eammis.den[i])/d;
    } else if (i==eammis.mesh-1) {eammis.dend[i]=(eammis.den[i]-eammis.den[i-1])/d;
    } else { eammis.dend[i]=(eammis.den[i+1]-eammis.den[i-1])/d/2; } }
  delete[] eammis.denx; delete[] eammis.deny;
  fp = fopen(fnamf_d, "w");
  fprintf(fp,"%20.15e %20.15e\n",eammis.denmin,eammis.denmax);
  for (int i=0; i<=eammis.mesh-1; i++) {
    fprintf(fp,"%20.15e %20.15e\n",eammis.den[i],eammis.dend[i]); }
  fclose(fp);
  }
  if (!fine_embed) {
  d=(eammis.embedmax-eammis.embedmin)/(double)(eammis.mesh-1);
  for (int i=0; i<=eammis.mesh-1; i++) {
    x=eammis.embedmin+d*(double)i;
    eammis.embed[i]=splint(eammis.embedx, eammis.embedy, mesh_embed, x); }
  for (int i=0; i<=eammis.mesh-1; i++) {
    if (i==0) {eammis.embedd[i]=(eammis.embed[i+1]-eammis.embed[i])/d;
    } else if (i==eammis.mesh-1) {eammis.embedd[i]=(eammis.embed[i]-eammis.embed[i-1])/d;
    } else { eammis.embedd[i]=(eammis.embed[i+1]-eammis.embed[i-1])/d/2; } }
  delete[] eammis.embedx; delete[] eammis.embedy;
  fp = fopen(fnamf_e, "w");
  fprintf(fp,"%20.15e %20.15e\n",eammis.embedmin,eammis.embedmax);
  for (int i=0; i<=eammis.mesh-1; i++) {
    fprintf(fp,"%20.15e %20.15e\n",eammis.embed[i],eammis.embedd[i]); }
  fclose(fp);
  }
  //std::ofstream fout("test.d");
  //for (int i=0; i<=eammis.mesh-1; i++) {
  //  fout<<i<<" "<<eammis.embed[i]<<" "<<eammis.embedd[i]<<std::endl; }
  //fout.close();

  /*set atom number */
  eammis.atom_number = atom_number(eammis.species);

  /* calculate to how many atoms eam potential should be applied */
  eammis.natoms = 0;
  for (int i=1; i<=atom.natom; i++) 
     if (eammis.atom_number == atom.anum[i])
      	eammis.natoms++;

 double dden =(eammis.denmax-eammis.denmin)/(double)(eammis.mesh-1);
 double dembed =(eammis.embedmax-eammis.embedmin)/(double)(eammis.mesh-1);
 double dpair =(eammis.pairmax-eammis.pairmin)/(double)(eammis.mesh-1);

 // create cubic interpolation
 for (int i=1; i<=eammis.mesh - 1; i++)
 {
	eammis.aden[i] = eammis.den[i] - eammis.den[i + 1] + eammis.dend[i + 1]*dden;
	eammis.bden[i] = -2 * (eammis.den[i] - eammis.den[i + 1]) - (eammis.dend[i] + eammis.dend[i+1])*dden;
	eammis.cden[i] = eammis.dend[i]*dden;
	eammis.dden[i] = eammis.den[i];

	eammis.apair[i] = eammis.pair[i] - eammis.pair[i + 1] + eammis.paird[i + 1]*dpair;
	eammis.bpair[i] = -2 * (eammis.pair[i] - eammis.pair[i + 1]) - (eammis.paird[i] + eammis.paird[i+1])*dpair;
	eammis.cpair[i] = eammis.paird[i]*dpair;
	eammis.dpair[i] = eammis.pair[i];

	eammis.aembed[i] = eammis.embed[i] - eammis.embed[i + 1] + eammis.embedd[i + 1]*dembed;
	eammis.bembed[i] = -2 * (eammis.embed[i] - eammis.embed[i + 1]) - (eammis.embedd[i] + eammis.embedd[i+1])*dembed;
	eammis.cembed[i] = eammis.embedd[i]*dembed;
	eammis.dembed[i] = eammis.embed[i];
 }
}

double density(double x)
{
        if (x > eammis.denmax*ang) { x = eammis.denmax*ang; }
        if (x < eammis.denmin*ang) { x = eammis.denmin*ang; }
	double xi= (double)(eammis.mesh - 1)*(x/ang - eammis.denmin)/(eammis.denmax - eammis.denmin) + 1.0;
	int ii=(int)xi;

	double xx = xi-(double)ii;

	return eammis.aden[ii] * xx * xx * xx + eammis.bden[ii] * xx * xx + eammis.cden[ii] * xx + eammis.dden[ii];
}

double density_d(double x)
{
        if (x > eammis.denmax*ang) { x = eammis.denmax*ang; }
        if (x < eammis.denmin*ang) { x = eammis.denmin*ang; }
	double xi= (double)(eammis.mesh - 1)*(x/ang - eammis.denmin)/(eammis.denmax - eammis.denmin) + 1.0;
	int ii=(int)xi;

	double xx = xi-(double)ii;	

	return (eammis.mesh - 1) * (3*eammis.aden[ii] * xx * xx + 2 * eammis.bden[ii] * xx + eammis.cden[ii]) / (eammis.denmax - eammis.denmin);
}

double embedded(double x)
{
        if (x > eammis.embedmax) { x = eammis.embedmax; }
        if (x < eammis.embedmin) { x = eammis.embedmin; }
	double xi= (double)(eammis.mesh - 1)*(x - eammis.embedmin)/(eammis.embedmax - eammis.embedmin) + 1.0;
	int ii=(int)xi;

	double xx = xi-(double)ii;	

	
	return eV*(eammis.aembed[ii] * xx * xx * xx + eammis.bembed[ii] * xx * xx + eammis.cembed[ii] * xx + eammis.dembed[ii]);
}

double embedded_d(double x)
{
        if (x > eammis.embedmax) { x = eammis.embedmax; }
        if (x < eammis.embedmin) { x = eammis.embedmin; }
	double xi= (double)(eammis.mesh - 1)*(x - eammis.embedmin)/(eammis.embedmax - eammis.embedmin)+1.0;
	int ii=(int)xi;

	double xx = xi-(double)ii;	

	return  (eammis.mesh - 1) * (3*eammis.aembed[ii] * xx * xx + 2*eammis.bembed[ii] * xx + eammis.cembed[ii])/(eammis.embedmax - eammis.embedmin);
}


double eam_pair(double x)
{
        if (x > eammis.pairmax*ang) { x = eammis.pairmax*ang; }
        if (x < eammis.pairmin*ang) { x = eammis.pairmin*ang; }
	double xi= (double)(eammis.mesh - 1)*(x/ang - eammis.pairmin)/(eammis.pairmax - eammis.pairmin) + 1.0;
	int ii=(int)xi;

	double xx = xi-(double)ii;	
	return eV*(eammis.apair[ii] * xx * xx * xx + eammis.bpair[ii] * xx * xx + eammis.cpair[ii] * xx + eammis.dpair[ii]);
}

double eam_pair_d(double x)
{
        if (x > eammis.pairmax*ang) { x = eammis.pairmax*ang; }
        if (x < eammis.pairmin*ang) { x = eammis.pairmin*ang; }
	double xi= (double)(eammis.mesh - 1)*(x/ang - eammis.pairmin)/(eammis.pairmax - eammis.pairmin) + 1.0;
	int ii=(int)xi;

	double xx = xi-(double)ii;	

	return (eammis.mesh - 1) * (3*eammis.apair[ii] * xx * xx + 2*eammis.bpair[ii] * xx + eammis.cpair[ii])/(eammis.pairmax - eammis.pairmin);
}


double splint(double xx[], double yy[], int n, double xi)
{
	int flag, i, *jun, *k;
	double dxp, dxm, *h, hi, hi2, *p, *q, *r, *s, sm, si, *sp, *x, *y, yi;

	//if (xi==xx[0]) return yy[0];
	//if (xi==xx[n-1]) return yy[n-1];
	if ((xi-xx[0]>-1.0e-8)&&(xi-xx[0]<1.0e-8)) return yy[0];
	if ((xi-xx[n-1]>-1.0e-8)&&(xi-xx[n-1]<1.0e-8)) return yy[n-1];

	x = (double *)malloc(n * sizeof(double));
	if(x == NULL)
	{
		fprintf(stderr, "Error : Out of memory  in splint()\n");
		return 0.;
	}
	y = (double *)malloc(n * sizeof(double));
	if(y == NULL)
	{
		fprintf(stderr, "Error : Out of memory  in splint()\n");
		free((char *)x);
		return 0.;
	}
	for(i = 0, p = xx, flag = 0; i < n - 1; i++, p++)
	{
		if(*p >= *(p + 1))
		{
			flag = 1;
			break;
		}
	}
	if(flag)
	{
		jun = (int *)malloc(n * sizeof(int));
		if(jun == NULL)
		{
			fprintf(stderr, "Error : Out of memory  in splint()\n");
			free((char *)x);
			free((char *)y);
			return 0.;
		}
		sortdi1(xx, n, jun);
		for(i = 0, p = x, q = y, k = jun; i < n; i++)
		{
			*p++ = *(xx + *k);
			*q++ = *(yy + *k++);
		}
		free((char *)jun);
	}
	else
	{
		for(i = 0, p = x, q = y, r = xx, s = yy; i < n; i++)
		{
			*p++ = *r++;
			*q++ = *s++;
		}
	}

	if(n < 2 || xi < *x || xi > *(x + n - 1))
	{
		fprintf(stderr, "Error : Illegal parameter  in splint()\n");
		free((char *)x);
		free((char *)y);
		return 0.;
	}
	for(i = 0, p = x; i < n; i++)
	{
		if(*p++ == xi)
		{
			free((char *)x);
			free((char *)y);
			return *(y + i);
		}
	}
	h = (double *)malloc(n * sizeof(double));
	if(h == NULL)
	{
		fprintf(stderr, "Error : Out of memory  in splint()\n");
		free((char *)x);
		free((char *)y);
		return 0.;
	}
	sp = (double *)malloc(n * sizeof(double));
	if(sp == NULL)
	{
		fprintf(stderr, "Error : Out of memory  in splint()\n");
		free((char *)x);
		free((char *)y);
		free((char *)h);
		return 0.;
	}

	subspl(x, y, n - 1, h, sp);
	for(i = 1, p = x + 1; i <= n; i++, p++)
	{
		if(*(p - 1) <= xi && xi < *p)
		{
			sm = *(sp + i - 1);
			si = *(sp + i);
			hi = *(h + i);
			hi2 = hi * hi;
			dxp = *p - xi;
			dxm = xi - *(p - 1);
			yi = (  sm * dxp * dxp * dxp + si * dxm * dxm * dxm
				  + (6. * *(y + i - 1) - hi2 * sm) * dxp
				  + (6. * *(y + i)     - hi2 * si) * dxm) / hi / 6.;
			free((char *)x);
			free((char *)y);
			free((char *)h);
			free((char *)sp);
			return yi;
		}
	}
	return 0;
}

void sortdi1(double a[], int n, int jun[])
{
	int i, j, k, l, *p, r, s, w, st1[32], st2[32];
	double x;

	if(n <= 1)
	{
		fprintf(stderr, "Error : n <= 1  in sortdi1()\n");
		return;
	}
	for(i = 0, p = jun; i < n; i++)	*p++ = i;
	s = 0;
	st1[0] = 0;
	st2[0] = n - 1;
	do
	{
		l = st1[s];
		r = st2[s];
		s--;
		if(r - l < 11)
		{
			i = l;
			while(i < r)
			{
				k = i++;
				j = i;
				while(*(a + *(jun + k)) > *(a + *(jun + j)))
				{
					w = *(jun + k);
					*(jun + k) = *(jun + j);
					*(jun + j) = w;
					if(j <= l + 1)	break;
					j = k--;
				}
			}
		}
		else
		{
			for(;;)
			{
				i = l;
				j = r;
				x = *(a + *(jun + (int)((l + r) / 2)));
				for(;;)
				{
					while(x > *(a + *(jun + i)))	i++;
					while(*(a + *(jun + j)) > x)	j--;
					if(i > j)	break;
					w = *(jun + i);
					*(jun + i++) = *(jun + j);
					*(jun + j--) = w;
					if(i > j)	break;
				}
				if(j - l < r - i)
				{
					if(i < r)
					{
						s++;
						st1[s] = i;
						st2[s] = r;
					}
					r = j;
				}
				else
				{
					if(l < j)
					{
						s++;
						st1[s] = l;
						st2[s] = j;
					}
					l = i;
				}
				if(l >= r)	break;
			}
		}
	} while(s >= 0);
	return;
}

void subspl(double x[], double y[], int n, double h[], double sp[])
{
	int i;
	double *dc, *du, g, hip1;

	dc = (double *)calloc(n + 1, sizeof(double));
	du = (double *)calloc(n + 1, sizeof(double));
	if(dc == NULL || du == NULL)
	{
		fprintf(stderr, "Error : out of memory in subspl().\n");
		return;
	}
	for(i = 1; i <= n; i++)	h[i] = x[i] - x[i - 1];
	for(i = 1; i < n; i++)
	{
		hip1 = x[i + 1] - x[i];
		sp[i] = 6.* ((y[i + 1] - y[i]) / h[i] / hip1 - (y[i] - y[i - 1]) / h[i] / h[i]);
		du[i] = h[i + 1] / h[i];
		dc[i] = 2. * (1. + du[i]);
	}
	dc[n] += 1.;
	dc[n - 1] += (h[n] / h[n - 1]);
	du[1] /= dc[1];
	sp[1] /= dc[1];
	for(i = 2; i <= n; i++)
	{
		g = 1. / (dc[i] - du[i - 1]);
		du[i] *= g;
		sp[i] = (sp[i] - sp[i - 1]) * g;
	}
	for(i = n - 1; i >= 1; i--)	sp[i] -= sp[i + 1] * du[i];
	sp[0] = sp[1];
	sp[n] = sp[n - 1];
	free(du);
	free(dc);
}

