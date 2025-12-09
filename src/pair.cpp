#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>
#include "myheader.h"

int atom_number(char* at);
void remove_carriage_return(std::string& line);

/* class constructor */
Pair::Pair(int ntype, double rcut, PotentialMode mode)
{
  this->ntype = ntype;
  this->rcut = rcut;
  this->mode = mode;
  initialized = false;
  fname[0] = NULL;
}

Pair::~Pair()
{
  delete[] types;
  delete[] params;
  delete[] pairnums;
  delete[] atomtypes;
}

/* read potential parameters from a POTFIT-format file*/
int Pair::ReadParams()
{
  printf("Starting pair potential initialization...\n");

  if (fname[0] == NULL) {
  if (type == BUCKINGHAMPOTTYPE)
    sprintf(fname, "pot/Buckingham.pot");
  else if (type == MORSEPOTTYPE)
    sprintf(fname, "pot/Morse.pot");
  else
    {
      printf("Pair potential type %d is not known.\n", type);
      mdmotion = 0;
      return 0;
    }
  }

  printf("Pair potential type %d (%s)\n", type, type == BUCKINGHAMPOTTYPE ? "Buckingham" : type == MORSEPOTTYPE ? "MORSE" : "!!!WRONG POTENTIAL TYPE!!!");
  printf("Pair potential parameters are being read from %s\n", fname);
  types = new int[ntype];
  
  //theoretical number of pair interactions that must present in potential file
  npairs = ntype * (ntype -1) / 2 + ntype;
  
  std::ifstream fin( this->fname );
  std::string line, arg1, arg2, arg3, arg4;
  char asp[3];
  
  // read the first lines with number of interactions and etc.
  getline(fin, line); remove_carriage_return(line);
  
  if (line.at(0)!='#' || line.at(1)!='F') {
    Output::ShowError("Incorrect first line in potential file, it sould start from #F ."); return 0;
  }
  else {
    // theoretically thehe must be npairs pair interactions in the file, so further we check it
    StringUtility::GetFirstArg(line, arg1);	StringUtility::GetFirstArg(line, arg1); StringUtility::GetFirstArg(line, arg1);
    
    if (atoi(arg1.c_str()) != npairs) {Output::ShowError("Wrong number of pairs in potential file");return 0; }
  }
  
  // read second line
  getline(fin, line); remove_carriage_return(line);
  if (line.at(0)!='#' || line.at(1)!='T')
    Output::ShowError("Incorrect second line in potential file, it sould start from #T .");
  
  // read third line and remember atomic species
  getline(fin, line); remove_carriage_return(line);
  if (line.at(0)!='#' || line.at(1)!='C') 
    Output::ShowError("Incorrect third line in potential file, it sould start from #C .");
  else {
    // remember atomic species
    for (int i = 0; i <= ntype; i++) {
      StringUtility::GetFirstArg(line, arg1);	
      if (i > 0) {strcpy(asp, arg1.c_str()); types[i-1] = atom_number(asp);}
    }
    
    atomtypes = new int[atom.natom + 1];
    
    for (int i = 1; i <= atom.natom; i++) {
      bool isTypeFound = false;
      
      for (int j = 0; j < ntype; j++) {
	if (types[j] == atom_number(atom.asp[i])) {
	  atomtypes[i] = j;
	  isTypeFound = true;
	  break;
	}
      }
      //std::cout << i << " " << atom.asp[i] << " " << atomtypes[i] <<std::endl;
      if (!isTypeFound)
	{
	  Output::ShowError("Atom type not found.");
	  return 0;
	}
    }
  }
  
  //read  lines from forth to sixth without any processing (only format check)
  getline(fin, line); remove_carriage_return(line);
  if (line.at(0)!='#' || line.at(1)!='#') 
    Output::ShowError("Incorrect forth line in potential file, it sould start from ## .");
  
  getline(fin, line); remove_carriage_return(line);
  if (line.at(0)!='#' || line.at(1)!='I') 
    Output::ShowError("Incorrect fifth line in potential file, it sould start from #I .");
  
  getline(fin, line); remove_carriage_return(line);
  if (line.at(0)!='#' || line.at(1)!='E') 
    Output::ShowError("Incorrect sixth line in potential file, it sould start from #E .");
  
  // allocate memory for potential parameters
  params = new double*[npairs];
  for (int i = 0; i< npairs; i++) {
    params[i] = new double[nparams];
  }
  
  //numerate atomic pairs
  pairnums = new int*[ntype];
  for (int i = 0; i< ntype; i++) {
    pairnums[i] =  new int[ntype];
  }	
  
  int counter = 0;	
  // read parameters of all pairwise interactions
  for (int i = 0; i < ntype; i++)
    for (int j = i; j < ntype; j++) {
      pairnums[i][j] = counter;
      pairnums[j][i] = pairnums[i][j];
      
      //skip blank lines
      while (getline(fin, line)) {  if (line.length()!=0) break; }
      
      //read cut-off distance
      getline(fin, line); remove_carriage_return(line);
      
      StringUtility::GetFirstArg(line, arg1);	StringUtility::GetFirstArg(line, arg2);
      
      if (strcmp(arg1.c_str(),"cutoff")!=0) { Output::ShowError("Incorrect line in potential file, cutoff radius not found"); return 0; }
      
      rcut = (double)atof(arg2.c_str());
      //rcut2 = rcut * rcut;
      
      //read line with minimum distance (no processing)
      getline(fin, line); remove_carriage_return(line);
      
      // !!! IMPORTANT PAIR INTERACTIONS MUST BE ORDERED LIKE THE FOLLOWING: for Zr O Y Ni the proper order is Zr-Zr Zr-O Zr-Y Zr-Ni O-O O-Y O-Ni Y-Y Y-Ni Ni-Ni !!!
      for (int k = 0; k < nparams; k++) {
	getline(fin, line); remove_carriage_return(line);
	StringUtility::GetFirstArg(line, arg1);	StringUtility::GetFirstArg(line, arg2);StringUtility::GetFirstArg(line, arg3);	StringUtility::GetFirstArg(line, arg4);			
	
	params[counter][k] = atof(arg2.c_str());
      }
      
      counter++;				
    }
  
  /*for (int i = 0; i < ntype; i++)
    {
    for (int j = 0; j < ntype; j++) 
    {
    std::cout << types[i]<<"-"<<types[j] <<":"<< params[i][j][0] << "	";
    }
    std::cout << "\n";
    }*/
  
  
  //std::cout << ntype << " " << counter << std::endl;
  fin.close(); fin.clear();
  
  initialized = true;
  printf("Pair potential initialization done.\n");
  return 1;
}

void Pair::Calculate()
{
  // if potential parameters were not read from potential file, then do it
  if (!initialized)
    ReadParams();
  
  // if potential is still not initialized, then something went wrong during reading potential file
  if (!initialized)
    {
      mdmotion = 0;
      return;
    }
  
  // reset some values in not combined mode, otherwise they are acuumulated
  if (mode != COMBINEDPOTMODE)
    Reset();
  

  int j, ix, iy, iz, pairnumber;
  double drx, dry, drz, rr, rr2, v0, vp0;
  
  /*loop over atoms*/
  for (int i=1; i<=atom.natom; i++)
    {
      if (book.alistnum[i]>0)
	{
	  /* loop over its neighbors */
	  for (int k=1; k<=book.alistnum[i]; k++)
	    {
	      j = book.alist[i][k][0];
	      
	      if (j>=i)
		{
		  if ((atom.QC==0)||(atom.repatom[i]==1)||(atom.repatom[j]==1))
		    {
		      ix = book.alist[i][k][1]; iy = book.alist[i][k][2]; iz = book.alist[i][k][3];
		      
		      drx = atom.rx[j]+cell.hmat[0][0]*ix+cell.hmat[0][1]*iy+cell.hmat[0][2]*iz - atom.rx[i];
		      dry = atom.ry[j]+cell.hmat[1][0]*ix+cell.hmat[1][1]*iy+cell.hmat[1][2]*iz - atom.ry[i];
		      drz = atom.rz[j]+cell.hmat[2][0]*ix+cell.hmat[2][1]*iy+cell.hmat[2][2]*iz - atom.rz[i];
		      rr2 = ( drx*drx + dry*dry + drz*drz );
		      
		      if (rr2 < this->rcut * this->rcut) 
			{
			  rr = sqrt(rr2);
			  pairnumber = pairnums[atomtypes[i]][atomtypes[j]];
			  //std::cout << i << " " << j << " " << atomtypes[i] << " " << atomtypes[j] << " " << pairnums[atomtypes[i]][atomtypes[j]] << std::endl;
			  //if (pairnums[atomtypes[i]][atomtypes[j]] > 20)
			  //	return; 
			  
			  vp0 = pairForce(rr, pairnumber); v0 = pairEnergy(rr, pairnumber); 
			  double vp0rr=vp0/rr;
			  
			  
			  if (j==i) { vp0rr /= 2.0; }
			  if ((atom.QC==0)||(atom.repatom[i]==1))
			    {
			      atom.fx[i]+=vp0rr*drx;
			      atom.fy[i]+=vp0rr*dry;
			      atom.fz[i]+=vp0rr*drz;
			      atom.epot[i]+=v0/2.0;
			      
			      //if (v0 > 0)
			      //std::cout << "Pair energy: " << v0 << "\n";
			      
			    }
			  if ((atom.QC==0)||(atom.repatom[j]==1))
			    {
			      atom.fx[j]-=vp0rr*drx;
			      atom.fy[j]-=vp0rr*dry;
			      atom.fz[j]-=vp0rr*drz;
			      if (j != i) {atom.epot[j]+=v0/2.0;}
			    }
			  if (atom.QC==0)
			    {
			      double ad1 = vp0rr/2.0;
			      atom.satom[i][0][0] += ad1*drx*drx;
			      atom.satom[i][0][1] += ad1*dry*drx;
			      atom.satom[i][1][1] += ad1*dry*dry;
			      atom.satom[i][0][2] += ad1*drz*drx;
			      atom.satom[i][1][2] += ad1*drz*dry;
			      atom.satom[i][2][2] += ad1*drz*drz;
			      atom.satom[j][0][0] += ad1*drx*drx;
			      atom.satom[j][0][1] += ad1*dry*drx;
			      atom.satom[j][1][1] += ad1*dry*dry;
			      atom.satom[j][0][2] += ad1*drz*drx;
			      atom.satom[j][1][2] += ad1*drz*dry;
			      atom.satom[j][2][2] += ad1*drz*drz;
			    } 
			  
			}
		    }
		}
	    }
	}
    }
  
  if (mode != COMBINEDPOTMODE)
    {
      atom.epotsum = 0;
      for (int i=1; i<=atom.natom; i++) { atom.epotsum += atom.epot[i]; }
    }
}

void Pair::Reset()
{
  //   Virial term reset
  cell.virx=0.0; cell.viry=0.0; cell.virz=0.0;
  
  for (int i=0; i<3; i++) { for (int j=0; j<3; j++) { 
      cell.dmat[i][j] = 0.0;
      for (int ii=1; ii<=atom.natom; ii++) { atom.satom[ii][i][j] = 0.0; }  }
  }
  
  for (int i=1; i<=atom.natom; i++) { atom.fx[i] = 0.0; atom.fy[i] = 0.0; atom.fz[i] = 0.0;  atom.epot[i] = 0.0; }
}

/* MORSE POTENTIAL */
Morse::Morse(int ntype, double rcut, PotentialMode mode):Pair(ntype, rcut, mode)
{
  nparams = 7;
  type = MORSEPOTTYPE;
}

/* Implement virtual function for pairwise energy */
/* first parameter - distance between atoms, second - number of interatomic pair */
double Morse::pairEnergy(double rr, int pn)
{
  double rang = rr * 1.0e10;
  return eV * (params[pn][0]*MorseExp(params[pn][1], rang, params[pn][2]) + params[pn][3]*MorseExp(params[pn][4], rang, params[pn][5]) + params[pn][6]); 
}

/* Implement virtual function for pairwise force */
double Morse::pairForce(double rr, int pn)
{
  double rang = rr * 1.0e10;
  return eV * 1.0e10 * (params[pn][0]*MorseExpDeriv(params[pn][1], rang, params[pn][2]) + params[pn][3]*MorseExpDeriv(params[pn][4], rang, params[pn][5]) + params[pn][6]); 	
}

double Morse::MorseExp(double a, double r, double r0)
{
  double s = a*(r-r0);
  return exp(-2*s)-2*exp(-s);
}

double Morse::MorseExpDeriv(double a, double r, double r0)
{
  double s = a*(r-r0);
  return 2*a*(exp(-s)-exp(-2*s));
}

/* Buckingham POTENTIAL */
Buckingham::Buckingham(int ntype, double rcut, PotentialMode mode):Pair(ntype, rcut, mode)
{
  nparams = 3;
  type = BUCKINGHAMPOTTYPE;
}

/* Implement virtual function for pairwise energy */
/* first parameter - distance between atoms, second - number of interatomic pair */
double Buckingham::pairEnergy(double rr, int pn)
{
  double r = (rr / params[pn][1]) / ang;	
  //if ((params[pn][0]*exp(-r) - params[pn][2]/(r*r*r*r*r*r)) > 0.01)
  //	std::cout << "OK" << pn << " " << params[pn][0] << " " << params[pn][1] << " " << params[pn][3] << "\n";
  
  return eV * (params[pn][0]*exp(-r) - params[pn][2]/(r*r*r*r*r*r)); 
}

/* Implement virtual function for pairwise force */
double Buckingham::pairForce(double rr, int pn)
{
  double r = (rr / params[pn][1]) / ang;
  
  return eV* 1.0e10 *((-params[pn][0]/params[pn][1]) * exp(-r) + 6 * (params[pn][2]/params[pn][1])/(r*r*r*r*r*r*r)); 	
}


