#include <fstream>
#include <iostream>
#include "myheader.h"

#if !defined __linux__ && !defined __APPLE__
#define snprintf sprintf_s
#endif

void capture();
void writeconfig(const char* fname);
void writecfgs(const char* fname);
void writecfge(const char* fname);
extern int confwrint, autocap;

/*
void writedata()
{
  std::ofstream foutene( "energy.d", std::ios::out | std::ios::app );
  foutene << istep << " " << atom.epotsum << " " << atom.ekinsum << std::endl;
  
}
*/
//void writedata( FILE *fp )
void writedata()
{
  fprintf(energyfile, "%d %22.17e %22.17e %22.17e\n",
	  istep, atom.epotsum/eV/atom.natom, atom.ekinsum/eV/atom.natom,
	  (atom.epotsum+atom.ekinsum)/eV/atom.natom);
  fprintf(stressfile, "%d %10.5e %10.5e %10.5e  %10.5e %10.5e %10.5e\n",
	  istep, cell.sgmmat[0][0]*1e-6, cell.sgmmat[1][1]*1e-6, cell.sgmmat[2][2]*1e-6,
	  cell.sgmmat[0][1]*1e-6, cell.sgmmat[1][2]*1e-6, cell.sgmmat[2][0]*1e-6);
  fprintf(ssnorfile, "%10.5e %10.5e %10.5e  %10.5e %10.5e %10.5e\n",
	  cell.hmat[0][0]*1e10, cell.hmat[1][1]*1e10, cell.hmat[2][2]*1e10,
	  cell.sgmmat[0][0]*1e-6, cell.sgmmat[1][1]*1e-6, cell.sgmmat[2][2]*1e-6);
  fprintf(cellfile, "%d %8.6f %8.6f %8.6f  %8.6f %8.6f %8.6f  %8.6f %8.6f %8.6f\n",
	  istep,
	  cell.hmat[0][0]*1e10, cell.hmat[1][0]*1e10, cell.hmat[2][0]*1e10,
	  cell.hmat[0][1]*1e10, cell.hmat[1][1]*1e10, cell.hmat[2][1]*1e10,
	  cell.hmat[0][2]*1e10, cell.hmat[1][2]*1e10, cell.hmat[2][2]*1e10);
  fprintf(tempfile, "%d %10.4f\n", istep, tempc);
  fflush(energyfile); fflush(stressfile); fflush(cellfile); fflush(ssnorfile);
  fflush(tempfile);
  if (ensemble==9)
    {
      fprintf(msdfile, "%e %e\n",(double)istep*dt,atom.MSD()); fflush(msdfile);
    }

  // Output config data
  if (confwrint > 0) {
    if (istep % confwrint == 0) {
      int num = istep / confwrint; 
      char filepath[80] = "CONFIG.SNAP"; char numc[10];
      if (num<10000) { snprintf(numc, sizeof(numc), "%04d", num);
      } else { snprintf(numc, sizeof(numc), "%d", num); }
      strcat(filepath,numc);
      writeconfig(filepath);
      char filepath_cfg[80] = "config";
      strcat(filepath_cfg,numc); strcat(filepath_cfg,".cfg");
      char filepath_cfge[80] = "config";
      strcat(filepath_cfge,numc); strcat(filepath_cfge,".cfge");
      writecfgs(filepath_cfg); writecfge(filepath_cfge);
      if (autocap) { capture(); }
    } }

  /*
  for (int i=1;i<=atom.natom;i++) {
    printf("%d %22.17e\n",i,atom.epot[i]);
  }
  */
}
/* void writedata_initialize()
{
  std::ofstream foutene( "energy.d", std::ios::trunc );
} */
 //int writedata_initialize(FILE *fp)
void writedata_initialize()
{
  if (energyfile != NULL) { fclose(energyfile); }
  energyfile = fopen("energy.d","w");
  fprintf(energyfile, "# step, Epot (eV, per atom), Ekin, Epot+Ekin\n");
  if (stressfile != NULL) { fclose(stressfile); }
  stressfile = fopen("stress.d","w");
  fprintf(stressfile, "# step, Stress xx (MPa), yy, zz, xy, yz, zx\n");
  if (ssnorfile != NULL) { fclose(ssnorfile); }
  ssnorfile = fopen("ssnormal.d","w");
  fprintf(ssnorfile, "# cell 11 22 33 Stress xx (MPa), yy, zz\n");
  if (cellfile != NULL) { fclose(cellfile); }
  cellfile = fopen("cell.d","w");
  fprintf(cellfile, "# step, Cell matrix 11 (A), 21, 31,  12, 22, 32,  13, 23, 33\n");
  if (ecnorfile != NULL) { fclose(ecnorfile); }
  ecnorfile = fopen("ecnormal.d","w");
  fprintf(ecnorfile, "# cell 11 22 33 Econst xx (GPa), yy,zz\n");
  if (ecallfile != NULL) { fclose(ecallfile); }
  ecallfile = fopen("ecall.d","w");
  fprintf(ecnorfile, "# cell 1x 1y 1z 2x 2y 2z 3x 3y 3z Econst 11 (GPa), 22, 33, 44, 55, 66, 12, 13, 14, 15, 16, 23, 24, 25, 26, 34, 35, 36, 45, 46, 56\n");
  if (msdfile != NULL) { fclose(msdfile); }
  msdfile = fopen("msd.d","w");
  fprintf(msdfile, "# time, MSD\n");
  if (miscfile != NULL) { fclose(miscfile); }
  miscfile = fopen("misc.d","w");
  if (tempfile !=NULL) { fclose(tempfile); }
  tempfile = fopen("temp.d","w");
}

void writedata_elastic_const()
{
  fprintf(ecnorfile, "%12.8e %12.8e %12.8e  %12.8e %12.8e %12.8e\n",
	  cell.hmat[0][0]*1e10, cell.hmat[1][1]*1e10, cell.hmat[2][2]*1e10,
	  cell.ecmat[0][0]/GPa, cell.ecmat[1][1]/GPa, cell.ecmat[2][2]/GPa);
  fflush(ecnorfile);
  fprintf(ecallfile, "%12.8e %12.8e %12.8e  %12.8e %12.8e %12.8e  %12.8e %12.8e %12.8e   %12.8e %12.8e %12.8e  %12.8e %12.8e %12.8e    %12.8e %12.8e %12.8e  %12.8e %12.8e %12.8e  %12.8e %12.8e %12.8e  %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n",
	  cell.hmat[0][0]/ang, cell.hmat[1][0]/ang, cell.hmat[2][0]/ang,
	  cell.hmat[0][1]/ang, cell.hmat[1][1]/ang, cell.hmat[2][1]/ang,
	  cell.hmat[0][2]/ang, cell.hmat[1][2]/ang, cell.hmat[2][2]/ang,
	  cell.ecmat[0][0]/GPa, cell.ecmat[1][1]/GPa, cell.ecmat[2][2]/GPa,
	  cell.ecmat[3][3]/GPa, cell.ecmat[4][4]/GPa, cell.ecmat[5][5]/GPa,
	  cell.ecmat[0][1]/GPa, cell.ecmat[0][2]/GPa, cell.ecmat[0][3]/GPa,
	  cell.ecmat[0][4]/GPa, cell.ecmat[0][5]/GPa,
	  cell.ecmat[1][2]/GPa, cell.ecmat[1][3]/GPa, cell.ecmat[1][4]/GPa,
	  cell.ecmat[1][5]/GPa,
	  cell.ecmat[2][3]/GPa, cell.ecmat[2][4]/GPa, cell.ecmat[2][5]/GPa,
	  cell.ecmat[3][4]/GPa, cell.ecmat[3][5]/GPa,
	  cell.ecmat[4][5]/GPa);
  fflush(ecallfile);
}
	  
