#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include "myheader.h"

double geamden(double r, int i);
double geamdend(double r, int i);
double geampair(double r, int i, int j);
double geampaird(double r, int i, int j);
double geamembed(double r, int i);
double geamembedd(double r, int i);
int geamsp(char *species);
int geamsp(int i);

void e_force_geam_am()
{
  double rr, rr2, drx, dry, drz, vp0, v0;
  int j, ix, iy, iz;
  /*
  double gre[16] = {2.556162,2.891814,2.885034,2.488746,2.750897,2.771916,2.886166,3.499723
		    ,2.481987,2.728100,2.860082,2.740840,3.196291,2.505979,2.933872,3.199978};
  double gfe[16] = {1.554485,1.106232,1.529021,2.007018,1.595417,2.336509,1.392302,0.647872 
		    ,1.885957,2.723710,3.086341,3.487340,0.544323,1.975299,1.863200,2.230909};
  double grhoe[16] = {22.150141,15.539255,21.319637,27.984706,22.770550,34.108882,20.226537,8.906840 
		      ,20.041463,29.354065,33.787168,37.234847,7.132600,27.206789,25.565138,30.879991};
  double grhos[16] = {22.150141,15.539255,21.319637,27.984706,22.770550,34.108882,20.226537,8.906840
		      ,20.041463,29.354065,33.787168,37.234847,7.132600,27.206789,25.565138,30.879991};
  double galp[16] = {7.669911,7.944536,8.086176,8.029633,7.605017,7.079952,6.942419,8.468412 
		     ,9.818270,8.393531,8.489528,8.900114,10.228708,8.679625,8.775431,8.559190};
  double gbet[16] = {4.090619,4.237086,4.312627,4.282471,4.056009,3.775974,3.702623,4.516486
		     ,5.236411,4.476550,4.527748,4.746728,5.455311,4.629134,4.680230,4.564902};
  double gaa[16] = {0.327584,0.266074,0.230728,0.439664,0.385412,0.449644,0.251519,0.134878
		    ,0.392811,0.708787,0.611679,0.882435,0.137518,0.421378,0.373601,0.424667};
  double gbb[16] = {0.468735,0.386272,0.336695,0.632771,0.545121,0.593713,0.313394,0.203093
		    ,0.646243,1.120373,1.032101,1.394592,0.225930,0.640107,0.570968,0.640054};
  double gkai[16] = {0.431307,0.425351,0.420755,0.413436,0.425578,0.413484,0.395132,0.425877 
		     ,0.170306,0.137640,0.176977,0.139209,0.5,0.5,0.5,0.5};
  double glam[16] = {0.862140,0.850703,0.841511,0.826873,0.851156,0.826967,0.790264,0.851753
		     ,0.340613,0.275280,0.353954,0.278417,1.0,1.0,1.0,1.0};
  double gffn0[16] = {-2.176490,-1.729619,-2.930281,-2.693996,-2.320473,-4.099542,-2.806783,-1.419644
		      ,-2.534992,-3.692913,-5.103845,-4.946281,-0.896473,-2.541799,-3.203773,-4.485793};
  double gffn1[16] = {-0.140035,-0.221025,-0.554034,-0.066073,-0.421263,-0.754764,-0.276173,-0.228622
		      ,-0.059605,-0.178812,-0.405524,-0.148818,-0.044291,-0.219415,-0.198262,-0.293129};
  double gffn2[16] = {0.285621,0.541558,1.489437,0.170482,0.966525,1.766503,0.893409,0.630069
		      ,0.193065,0.380450,1.112997,0.365057,0.162232,0.733381,0.683779,0.990148};
  double gffn3[16] = {-1.750834,-0.967036,-0.886809,-2.457442,-0.932685,1.578274,-1.637201,-0.560952
		      ,-2.282322,-3.133650,-3.585325,-4.432406,-0.689950,-1.589003,-2.321732,-3.202516};
  double gff0[16] = {-2.19,-1.75,-2.98,-2.70,-2.36,-4.17,-2.83,-1.44
		     ,-2.54,-3.71,-5.14,-4.96,-0.90,-2.56,-3.22,-4.51};
  double gff1[16] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double gff2[16] = {0.702991,0.983967,2.283863,0.282257,1.966273,3.474733,0.929508,0.921049
		     ,0.200269,0.875874,1.640098,0.661935,0.122838,0.705845,0.608587,0.928602};
  double gff3[16] = {0.683705, 0.520904, 0.494127, 0.102879, 1.396717, 2.288323,-0.682320, 0.108847
		     ,-0.148770,0.776222,0.221375,0.348147,-0.226010,-0.687140,-0.750710,-0.981870};
  double geta[16] = {0.921150,1.149461,1.286960,0.509860,1.399758,1.393490,0.779208,1.172361
		     ,0.391750,0.790879,0.848843,-0.582714,0.431425,0.694608,0.558572,0.597133};
  double gffe[16] = {-2.191675,-1.751274,-2.981365,-2.700493,-2.362609,-4.174332,-2.829437,-1.440494
		     ,-2.539945,-3.712093,-5.141526,-4.961306,-0.899702,-2.559307,-3.219176,-4.509025};
  */
  double ki, kif;
   if (geam.initialize) {
    if (geam.rhob) { delete[] geam.rhob; geam.rhob = NULL; }
    if (geam.sp)   { delete[] geam.sp;   geam.sp   = NULL; }
   geam.rhob = new double[atom.natom+1];
   geam.sp   = new int[atom.natom+1];
    for (int i=1; i<=atom.natom; i++) {
      geam.sp[i] = geamsp(i);
    }
    geam.initialize = false;
  }

  // set cutoff radius (rcut and frc)
  rcut = 3.499723*2*1.0e-10;
  if (book.frc < rcut) {
    book.frc = rcut * 1.15; printf("### WARNING: frc is corrected (GEAM)\n");}
  double rcut2 = rcut * rcut;
  // Initialize rhob
  for (int i=1; i<=atom.natom; i++) { geam.rhob[i] = 0.0; }
  // Initialize stress-related variables
  cell.virx = 0.0; cell.viry = 0.0; cell.virz = 0.0;
  for (int i=0; i<3; i++) { for (int j=0; j<3; j++) { 
      cell.dmat[i][j] = 0.0;
      for (int ii=1; ii<=atom.natom; ii++) { atom.satom[ii][i][j] = 0.0; }  } }
  // Initialize Force and E_pot
  for (int i=1; i<=atom.natom; i++) {
    atom.fx[i] = 0.0; atom.fy[i] = 0.0; atom.fz[i] = 0.0;  atom.epot[i] = 0.0;
  }

  // First loop: Chg dns calc
  for (int i=1; i<=atom.natom; i++) {
    if (book.alistnum[i]>0) {
      for (int k=1; k<=book.alistnum[i]; k++) {
	j  = book.alist[i][k][0];
	if (j>=i) {
	  ix = book.alist[i][k][1]; iy = book.alist[i][k][2]; iz = book.alist[i][k][3];
	  rr2 = atom.Dist2(i,j,ix,iy,iz);
	  if (rr2 < rcut2) {
	    rr = sqrt(rr2);
	    //rhob[i] = rhob[i] + geamden(rr/ang, geamsp(atom.asp[j]));
	    geam.rhob[i] = geam.rhob[i] + geamden(rr/ang, geam.sp[j]);
	    if (j!=i) {
	      //rhob[j] = rhob[j] + geamden(rr/ang, geamsp(atom.asp[i])); }
	      geam.rhob[j] = geam.rhob[j] + geamden(rr/ang, geam.sp[i]); }
	  }
	} //j>=i
      }
    }
  }
  for (int i=1; i<=atom.natom; i++) {
    //atom.epot[i] = geamembed(rhob[i], geamsp(atom.asp[i])) * eV; }
    atom.epot[i] = geamembed(geam.rhob[i], geam.sp[i]) * eV; }
  
  // Second loop: Force calc
  for (int i=1; i<=atom.natom; i++) {
    if (book.alistnum[i]>0) {
      for (int k=1; k<=book.alistnum[i]; k++) {
	j  = book.alist[i][k][0];
	if (j>=i) {
	  if ((atom.QC==0)||(atom.repatom[i]==1)||(atom.repatom[j]==1)) {
	    ix = book.alist[i][k][1]; iy = book.alist[i][k][2]; iz = book.alist[i][k][3];
	    rr2 = atom.Dist2(i,j,ix,iy,iz);
	    //	if ((rr2 < rcut2)&&(rr2>1.0e-30)) {
	    if (rr2 < rcut2) {
	      rr = sqrt(rr2);
	      drx=-atom.Dx(i,j,ix,iy,iz); dry=-atom.Dy(i,j,ix,iy,iz); drz=-atom.Dz(i,j,ix,iy,iz);
	      //double vij = geampair(rr/ang, geamsp(atom.asp[i]), geamsp(atom.asp[j])) * eV;
	      //double dv = geampaird(rr/ang, geamsp(atom.asp[i]), geamsp(atom.asp[j]));
	      //double df_i = geamembedd(rhob[i],geamsp(atom.asp[i]));
	      //double df_j = geamembedd(rhob[j],geamsp(atom.asp[j]));
	      //double dr_i = geamdend(rr/ang, geamsp(atom.asp[i]));
	      //double dr_j = geamdend(rr/ang, geamsp(atom.asp[j]));
	      double vij = geampair(rr/ang, geam.sp[i], geam.sp[j]) * eV;
	      double dv = geampaird(rr/ang, geam.sp[i], geam.sp[j]);
	      double df_i = geamembedd(geam.rhob[i],geam.sp[i]);
	      double df_j = geamembedd(geam.rhob[j],geam.sp[j]);
	      double dr_i = geamdend(rr/ang, geam.sp[i]);
	      double dr_j = geamdend(rr/ang, geam.sp[j]);
	      double f = -(df_i*dr_j+df_j*dr_i+dv)*eV/ang;
	      double ki=f/rr;
	      if (j==i) ki=ki/2.0;
	      double kif=ki/2.0;

	      if ((atom.QC==0)||(atom.repatom[i]==1)) {
		atom.fx[i] = atom.fx[i] + ki*drx;
		atom.fy[i] = atom.fy[i] + ki*dry;
		atom.fz[i] = atom.fz[i] + ki*drz;
		atom.epot[i]=atom.epot[i]+vij/2.0;
	      }
	      if ((atom.QC==0)||(atom.repatom[j]==1)) {
		atom.fx[j] = atom.fx[j] - ki*drx;
		atom.fy[j] = atom.fy[j] - ki*dry;
		atom.fz[j] = atom.fz[j] - ki*drz;
		if (j != i) {
		  atom.epot[j]=atom.epot[j]+vij/2.0;}
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
		atom.satom[i][0][0] = atom.satom[i][0][0] - kif*drx*drx;
		atom.satom[i][0][1] = atom.satom[i][0][1] - kif*dry*drx;
		atom.satom[i][1][1] = atom.satom[i][1][1] - kif*dry*dry;
		atom.satom[i][0][2] = atom.satom[i][0][2] - kif*drz*drx;
		atom.satom[i][1][2] = atom.satom[i][1][2] - kif*drz*dry;
		atom.satom[i][2][2] = atom.satom[i][2][2] - kif*drz*drz;
		atom.satom[j][0][0] = atom.satom[j][0][0] - kif*drx*drx;
		atom.satom[j][0][1] = atom.satom[j][0][1] - kif*dry*drx;
		atom.satom[j][1][1] = atom.satom[j][1][1] - kif*dry*dry;
		atom.satom[j][0][2] = atom.satom[j][0][2] - kif*drz*drx;
		atom.satom[j][1][2] = atom.satom[j][1][2] - kif*drz*dry;
		atom.satom[j][2][2] = atom.satom[j][2][2] - kif*drz*drz;
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
      atom.satom[i][1][0] = atom.satom[i][0][1];
      atom.satom[i][2][0] = atom.satom[i][0][2];
      atom.satom[i][2][1] = atom.satom[i][1][2];}
    for (int ii=1; ii<=atom.natom; ii++) {
      for (int i=0; i<3; i++) { for (int j=0; j<3; j++) { 
	  cell.dmat[i][j] = cell.dmat[i][j] - atom.satom[ii][i][j];
	} } }
    cell.virx = cell.dmat[0][0]; cell.viry = cell.dmat[1][1]; cell.virz = cell.dmat[2][2];
    cell.volume = cell.Getvolume();
    for (int ii=1; ii<=atom.natom; ii++) {
      for (int i=0; i<3; i++) { for (int j=0; j<3; j++) {
	  atom.satom[ii][i][j] = atom.satom[ii][i][j] * (double)atom.natom / cell.volume;
	} } }
  }
  atom.epotsum=0.0;
  for (int i=1; i<=atom.natom; i++) {
    atom.epotsum += atom.epot[i]; }
}

double geamden(double r, int i)
{
  double den = geam.gfe[i] * exp(-geam.gbet[i]*(r/geam.gre[i]-1.0)) /
    (1.0 + pow(r/geam.gre[i]-geam.glam[i], 20.0));
  return den;
}

double geamdend(double r, int i)
{
  double x = r/geam.gre[i] - geam.glam[i];
  double dend = -geam.gfe[i] * exp(-geam.gbet[i] * (r/geam.gre[i]-1.0))/geam.gre[i]
    * (geam.gbet[i]*(1.0+pow(x,20.0)) + 20.0*pow(x,19.0)) / ((1.0+pow(x,20.0))*(1.0+pow(x,20.0)));
  return dend;
}

double geampair(double r, int i, int j)
{
  double xi, xj, pair;
  xi = geam.gaa[i]*exp(-geam.galp[i]*(r/geam.gre[i]-1.0)) / 
    (1.0+pow(r/geam.gre[i]-geam.gkai[i],20.0))
    -  geam.gbb[i]*exp(-geam.gbet[i]*(r/geam.gre[i]-1.0)) / 
    (1.0+pow(r/geam.gre[i]-geam.glam[i],20.0));
  if (i==j) {
    pair = xi;
  } else {
    xj = geam.gaa[j]*exp(-geam.galp[j]*(r/geam.gre[j]-1.0)) / 
      (1.0+pow(r/geam.gre[j]-geam.gkai[j],20.0))
      -  geam.gbb[j]*exp(-geam.gbet[j]*(r/geam.gre[j]-1.0)) / 
      (1.0+pow(r/geam.gre[j]-geam.glam[j],20.0));
    pair = 0.5*(geamden(r,j)/geamden(r,i)*xi+geamden(r,i)/geamden(r,j)*xj);
  }
  return pair;
}

double geampaird(double r, int i, int j)
{
  double x, y, xi, xj, xdi, xdj, fi, fj, fdi, fdj, paird;
  x = r/geam.gre[i]-geam.gkai[i];
  y = r/geam.gre[i]-geam.glam[i];
  xdi = -geam.gaa[i]*exp(-geam.galp[i]*(r/geam.gre[i]-1.0))/geam.gre[i]
    * (geam.galp[i]*(1.0+pow(x,20.0)) + 20.0*pow(x,19.0)) / ((1.0+pow(x,20.0))*(1.0+pow(x,20.0)))
    +    geam.gbb[i]*exp(-geam.gbet[i]*(r/geam.gre[i]-1.0))/geam.gre[i]
    * (geam.gbet[i]*(1.0+pow(y,20.0)) + 20.0*pow(y,19.0)) / ((1.0+pow(y,20.0))*(1.0+pow(y,20.0)));
  if (i==j) {
    paird = xdi;
  } else {
    x = r/geam.gre[j]-geam.gkai[j];
    y = r/geam.gre[j]-geam.glam[j];
    xdj = -geam.gaa[j]*exp(-geam.galp[j]*(r/geam.gre[j]-1.0))/geam.gre[j]
      * (geam.galp[j]*(1.0+pow(x,20.0)) + 20.0*pow(x,19.0)) / ((1.0+pow(x,20.0))*(1.0+pow(x,20.0)))
      +    geam.gbb[j]*exp(-geam.gbet[j]*(r/geam.gre[j]-1.0))/geam.gre[j]
      * (geam.gbet[j]*(1.0+pow(y,20.0)) + 20.0*pow(y,19.0)) / ((1.0+pow(y,20.0))*(1.0+pow(y,20.0)));
    xi  = geampair(r,i,i); xj  = geampair(r,j,j);
    fi  = geamden(r,i);    fj  = geamden(r,j);
    fdi = geamdend(r,i);   fdj = geamdend(r,j);
    paird = 0.5*( ( fi*fj*xdi+fi*fdj*xi-fdi*fj*xi)/(fi*fi)
		  +(fi*fj*xdj+fdi*fj*xj-fi*fdj*xj)/(fj*fj) );
  }
  return paird;
}

double geamembed(double r, int i)
{
  double x, embed;
  if (r<0.85*geam.grhoe[i]) {
    x = r/(0.85*geam.grhoe[i])-1.0;
    embed = geam.gffn0[i]+geam.gffn1[i]*x+geam.gffn2[i]*x*x+geam.gffn3[i]*x*x*x;
  } else if (r<1.15*geam.grhoe[i]) {
    x = r/geam.grhoe[i]-1.0;
    embed = geam.gff0[i]+geam.gff1[i]*x+geam.gff2[i]*x*x+geam.gff3[i]*x*x*x;
  } else {
    x = pow(r/geam.grhos[i], geam.geta[i]);
    embed = geam.gffe[i]*(1.0-log(x))*x;
  }
  return embed;
}

double geamembedd(double r, int i)
{
  double x, embedd;
  if (r<0.85*geam.grhoe[i]) {
    x = r/(0.85*geam.grhoe[i])-1.0;
    embedd = (geam.gffn1[i]+2.0*geam.gffn2[i]*x+3.0*geam.gffn3[i]*x*x) / 
      (0.85*geam.grhoe[i]);
  } else if (r<1.15*geam.grhoe[i]) {
    x = r/geam.grhoe[i]-1.0;
    embedd = (geam.gff1[i]+2.0*geam.gff2[i]*x+3.0*geam.gff3[i]*x*x) / geam.grhoe[i];
  } else {
    x = pow(r/geam.grhos[i], geam.geta[i]-1.0);
    embedd = -geam.gffe[i]*geam.geta[i]/geam.grhos[i]*x
      *(geam.geta[i]*log(r/geam.grhos[i]));
  }
  return embedd;
}

int geamsp(char *species)
{
  int sp;
  if      (strcmp(species, "Cu") == 0) { sp = 0; }
  else if (strcmp(species, "Ag") == 0) { sp = 1; }
  else if (strcmp(species, "Au") == 0) { sp = 2; }
  else if (strcmp(species, "Ni") == 0) { sp = 3; }
  else if (strcmp(species, "Pd") == 0) { sp = 4; }
  else if (strcmp(species, "Pt") == 0) { sp = 5; }
  else if (strcmp(species, "Al") == 0) { sp = 6; }
  else if (strcmp(species, "Pb") == 0) { sp = 7; }
  else if (strcmp(species, "Fe") == 0) { sp = 8; }
  else if (strcmp(species, "Mo") == 0) { sp = 9; }
  else if (strcmp(species, "Ta") == 0) { sp = 10; }
  else if (strcmp(species, "W ") == 0) { sp = 11; }
  else if (strcmp(species, "Mg") == 0) { sp = 12; }
  else if (strcmp(species, "Co") == 0) { sp = 13; }
  else if (strcmp(species, "Ti") == 0) { sp = 14; }
  else if (strcmp(species, "Zr") == 0) { sp = 15; }
  else { printf("GEAM ERROR: Atom %s is not supported.\n",species); sp = 0; return sp; }
  return sp;
}
int geamsp(int i)
{
  int sp;
  if      (atom.anum[i] == 29) { sp = 0; }
  else if (atom.anum[i] == 47) { sp = 1; }
  else if (atom.anum[i] == 79) { sp = 2; }
  else if (atom.anum[i] == 28) { sp = 3; }
  else if (atom.anum[i] == 46) { sp = 4; }
  else if (atom.anum[i] == 78) { sp = 5; }
  else if (atom.anum[i] == 13) { sp = 6; }
  else if (atom.anum[i] == 82) { sp = 7; }
  else if (atom.anum[i] == 26) { sp = 8; }
  else if (atom.anum[i] == 42) { sp = 9; }
  else if (atom.anum[i] == 73) { sp = 10; }
  else if (atom.anum[i] == 74) { sp = 11; }
  else if (atom.anum[i] == 12) { sp = 12; }
  else if (atom.anum[i] == 27) { sp = 13; }
  else if (atom.anum[i] == 22) { sp = 14; }
  else if (atom.anum[i] == 40) { sp = 15; }
  else { printf("GEAM ERROR: Atom %s is not supported.\n",atom.asp[i]); sp = 0; return sp; }
  return sp;
}
