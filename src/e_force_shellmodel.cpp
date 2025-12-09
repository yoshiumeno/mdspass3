//
// Shell model FORCE & ENERGY CALCULATION for PbTiO3
//     
//  Last update 21.02.2007
// re-written for C++ 12.07.2024

#include <iostream>
#include<string.h>
#include<stdlib.h>
#include<fstream>
#include "myheader.h"
#define dsq(a) ((a)*(a))
void shell_delete();
void shell_alloc();
void resetmat(double a[3][3]);
double erfc(double x);
void recipro(double h[3][3], double r[3][3]);
int ksp(char* species);

void e_force_shellmodel()
{
   if (shell.initialize) {
      shell_delete();
      printf("SHELL delete DONE\n");
      shell_alloc();
      printf("SHELL alloc DONE\n");
      shell.initialize = false;
      printf("SHELL INITIALIZE DONE\n");
   }
   //
   //double alarge[4], rho[4], clarge[4], rcsr[4], rtp[4];
   //double q[3][2], springc[3], springc2[3], rec[3][3];
   //int kint[3][3];
   int istat;
   double rec[3][3];
   double pewaldeach1, pewaldeach2, pewaldeach3, pewaldeach4;
   double f, fxx, fyy, fzz;
   double gbx, gby, gbz, gbabs, gbabs2;
   double ew2a, ew2b;

   //Ewald parameter (1/m)
   double ewaldgm = pow(double(atom.natom)*dsq(M_PI/cell.volume), 1.0/6.0);
   //Ewald 1st term
   double acc = pow(1.0e-1, double(shell.naccew));
   double rce1 = sqrt(-log(acc))/ewaldgm; // (m)
   //Ewald 2nd term
   double deth = cell.volume;
   recipro(cell.hmat, rec);
   double rec0 = 2.0*M_PI*sqrt(dsq(rec[0][0])+dsq(rec[1][0])+dsq(rec[2][0]))/deth;
   double rec1 = 2.0*M_PI*sqrt(dsq(rec[0][1])+dsq(rec[1][1])+dsq(rec[2][1]))/deth;
   double rec2 = 2.0*M_PI*sqrt(dsq(rec[0][2])+dsq(rec[1][2])+dsq(rec[2][2]))/deth;
   double gce2 = 2.0*ewaldgm*sqrt(-log(acc));
   int mewald2x = int(gce2/rec0)+1;
   int mewald2y = int(gce2/rec1)+1;
   int mewald2z = int(gce2/rec2)+1;
   //---Cutoff radius
   double rc = 0.0;
   rc = std::max(rce1, rc);
   for (int i = 0; i < 4; i++) {
      rc = std::max(shell.rcsr[i]*1.0e-10, rc);
   }
   double r2;
   double rc2 = dsq(rc);

   //Virial term reset
   cell.virx = 0.0; cell.viry = 0.0; cell.virz = 0.0;
   resetmat(cell.dmat);
   //Reset energy and force
   atom.epotsum = 0.0;
   for (int i = 1; i <= atom.natom; i++) {
      atom.fx[i] = 0.0; atom.fy[i] = 0.0; atom.fz[i] = 0.0;  atom.epot[i] = 0.0;
   }
   double pewald0 = 0.0; double pewald1 = 0.0; double pewald2 = 0.0; double pewald3 = 0.0; double pewald4 = 0.0;
   double psr = 0.0; double pharmonic = 0.0;
   for (int i = 1; i <= atom.natom/2; i++) {
      for (int j = 0; j < 2; j++) {
         shell.tmpfx[i][j] = 0.0; shell.tmpfy[i][j] = 0.0; shell.tmpfz[i][j] = 0.0;
      }
   }
   //Copy r -> tmpr
   for (int i = 1; i <= atom.natom/2; i++) {
      shell.tmprx[i][0] = atom.rx[i];
      shell.tmpry[i][0] = atom.ry[i];
      shell.tmprz[i][0] = atom.rz[i];
      shell.tmprx[i][1] = atom.rx[i+atom.natom/2];
      shell.tmpry[i][1] = atom.ry[i+atom.natom/2];
      shell.tmprz[i][1] = atom.rz[i+atom.natom/2];
   }
   //----Ewald 1st term and Short-range interaction
   int j, ix, iy, iz; double rr2, r, rr;
   for (int i = 1; i <= atom.natom/2; i++) {
      for (int kk = 1; kk <= book.alistnum[i]; kk++) {
         j = book.alist[i][kk][0];
         if ((j>atom.natom/2)||(j<i)) { continue; }
   	   ix = book.alist[i][kk][1]; iy = book.alist[i][kk][2]; iz = book.alist[i][kk][3];
	      //rr2 = atom.Dist2(i,j,ix,iy,iz);
         for (int k = 0; k < 2; k++) {
            for (int l = 0; l < 2; l++) {
               double xx = shell.tmprx[i][k] - shell.tmprx[j][k];
               double yy = shell.tmpry[i][k] - shell.tmpry[j][k];
               double zz = shell.tmprz[i][k] - shell.tmprz[j][k];
               xx += cell.hmat[0][0]*double(ix)+cell.hmat[0][1]*double(iy)+cell.hmat[0][2]*double(iz);
               yy += cell.hmat[1][0]*double(ix)+cell.hmat[1][1]*double(iy)+cell.hmat[1][2]*double(iz);
               zz += cell.hmat[2][0]*double(ix)+cell.hmat[2][1]*double(iy)+cell.hmat[2][2]*double(iz);
               int istat = 0;
               if ((abs(xx)<rc)&&(abs(yy)<rc)&&(abs(yy)<rc)) {
                  if ((abs(xx)>1.0e-15)||(abs(yy)>1.0e-15)||(abs(zz)>1.0e-15)) {
                     istat = 1;
                  }
               }
               if (istat == 1) {
                  r2 = dsq(xx)+dsq(yy)+dsq(zz);
                  r = sqrt(r2);
                  rr = r*1.0e10;
                  //Ewald 1st term (Real space)
                  bool if_first = true;
                  if (r > rce1) { if_first = false; }
                  if ((i == j)&&(rr < 1.0)) { if_first = false; }
                  if (if_first) {
                     //Energy (J)
                     double xxx = ewaldgm*r;
                     pewaldeach1 = shell.q[ksp(atom.asp[i])][k]*shell.q[ksp(atom.asp[j])][l]*shell.efpe*erfc(double(ewaldgm*r))/r;
                     // Force(N)
                     f = shell.q[ksp(atom.asp[i])][k]*shell.q[ksp(atom.asp[j])][l]*shell.efpe*(erfc(ewaldgm*r)/r2
                          +2.0*ewaldgm/sqrt(M_PI)*exp(-dsq(ewaldgm*r))/r);
                     fxx = f/r*xx;
                     fyy = f/r*yy;
                     fzz = f/r*zz;
                     if (i == j) {
                        shell.tmpfx[i][k] += fxx;
                        shell.tmpfy[i][k] += fyy;
                        shell.tmpfz[i][k] += fzz;
                        cell.dmat[0][0] += fxx*xx/2.0;
                        cell.dmat[0][1] += fxx*yy/2.0;
                        cell.dmat[0][2] += fxx*zz/2.0;
                        cell.dmat[1][0]  = cell.dmat[0][1];
                        cell.dmat[1][1] += fyy*yy/2.0;
                        cell.dmat[1][2] += fyy*zz/2.0;
                        cell.dmat[2][0]  = cell.dmat[0][2];
                        cell.dmat[2][1]  = cell.dmat[1][2];
                        cell.dmat[2][2] += fzz*zz/2.0;
                        pewald1 += pewaldeach1/2.0;
                     } else {
                        shell.tmpfx[i][k] += fxx;
                        shell.tmpfy[i][k] += fyy;
                        shell.tmpfz[i][k] += fzz;
                        shell.tmpfx[j][l] -= fxx;
                        shell.tmpfy[j][l] -= fyy;
                        shell.tmpfz[j][l] -= fzz;                       
                        cell.dmat[0][0] += fxx*xx;
                        cell.dmat[0][1] += fxx*yy;
                        cell.dmat[0][2] += fxx*zz;
                        cell.dmat[1][0]  = cell.dmat[0][1];
                        cell.dmat[1][1] += fyy*yy;
                        cell.dmat[1][2] += fyy*zz;
                        cell.dmat[2][0]  = cell.dmat[0][2];
                        cell.dmat[2][1]  = cell.dmat[1][2];
                        cell.dmat[2][2] += fzz*zz;
                        pewald1 += pewaldeach1;
                     }
                  } //if_first end
                  bool if_short = true;
                  if (shell.kint[ksp(atom.asp[i])][ksp(atom.asp[j])]==0) { if_short = false; }
                  if (k != 1) { if_short = false; }
                  if (l != 1) { if_short = false; }
                  int kitr = shell.kint[ksp(atom.asp[i])][ksp(atom.asp[j])];
                  if (rr > shell.rcsr[kitr]) {if_short = false; }
                  if (if_short) {
                     f = (shell.alarge[kitr]/shell.rho[kitr]*exp(-rr/shell.rho[kitr])-6.0*shell.clarge[ kitr]/pow(rr,7))*eV*1.0e10;
                     double psreach = (shell.alarge[kitr]*exp(-rr/shell.rho[kitr])-shell.clarge[ kitr]/pow(rr,6))*eV;
                     // Tapering (rt < rr < rc)
                     if ((rr > shell.rcsr[kitr]-shell.rtp[kitr])&&(rr <= shell.rcsr[kitr])) {
                        psreach = psreach*0.50
                             *(1.0+cos(M_PI*(rr-shell.rcsr[kitr]-shell.rtp[kitr])/(shell.rcsr[kitr]-shell.rcsr[kitr]-shell.rtp[kitr])));
                        f = ((shell.alarge[kitr]/shell.rho[kitr]*exp(-rr/shell.rho[kitr])-6.0*shell.clarge[ kitr]/pow(rr,7))
                             *0.50*(1.0+cos(M_PI*(rr-(shell.rcsr[kitr]-shell.rtp[kitr]))/shell.rtp[kitr]))
                             -(shell.alarge[kitr]*exp(-rr/shell.rho[kitr])-shell.clarge[ kitr]/pow(rr,6))
                             *(-0.50*M_PI/shell.rtp[kitr]*sin(M_PI*(rr-(shell.rcsr[kitr]-shell.rtp[kitr]))/shell.rtp[kitr])))*eV*1.0e10;
                     }
                     fxx=f/r*xx;
                     fyy=f/r*yy;
                     fzz=f/r*zz;
                     if (i == j) {
                        shell.tmpfx[i][1] += fxx;
                        shell.tmpfy[i][1] += fyy;
                        shell.tmpfz[i][1] += fzz;
                        cell.dmat[0][0] += fxx*xx/2.0;
                        cell.dmat[0][1] += fxx*yy/2.0;
                        cell.dmat[0][2] += fxx*zz/2.0;
                        cell.dmat[1][0]  = cell.dmat[0][1];
                        cell.dmat[1][1] += fyy*yy/2.0;
                        cell.dmat[1][2] += fyy*zz/2.0;
                        cell.dmat[2][0]  = cell.dmat[0][2];
                        cell.dmat[2][1]  = cell.dmat[1][2];
                        cell.dmat[2][2] += fzz*zz/2.0;
                        psr += psreach/2.0;
                     } else {
                        shell.tmpfx[i][1] += fxx;
                        shell.tmpfy[i][1] += fyy;
                        shell.tmpfz[i][1] += fzz;
                        shell.tmpfx[j][1] -= fxx;
                        shell.tmpfy[j][1] -= fyy;
                        shell.tmpfz[j][1] -= fzz;
                        cell.dmat[0][0] += fxx*xx;
                        cell.dmat[0][1] += fxx*yy;
                        cell.dmat[0][2] += fxx*zz;
                        cell.dmat[1][0]  = cell.dmat[0][1];
                        cell.dmat[1][1] += fyy*yy;
                        cell.dmat[1][2] += fyy*zz;
                        cell.dmat[2][0]  = cell.dmat[0][2];
                        cell.dmat[2][1]  = cell.dmat[1][2];
                        cell.dmat[2][2] += fzz*zz;
                        psr += psreach;
                     }
                  } //if_short end
               } //if istat == 1 end
            }
         }
      }
   } //Ewald 1st term & short-range end
   //Ewald 2nd term (Reciprocal space)
   for (int mx = -mewald2x; mx <= mewald2x; mx++) {
   for (int my = -mewald2y; my <= mewald2y; my++) {
   for (int mz = -mewald2z; mz <= mewald2z; mz++) {
      if ((mx == 0)&&(my == 0)&&(mz == 0)) { continue; }
         gbx = 2.0*M_PI/cell.volume*(rec[0][0]*double(mx)+rec[0][1]*double(my)+rec[0][2]*double(mz));
         gby = 2.0*M_PI/cell.volume*(rec[1][0]*double(mx)+rec[1][1]*double(my)+rec[1][2]*double(mz));
         gbz = 2.0*M_PI/cell.volume*(rec[2][0]*double(mx)+rec[2][1]*double(my)+rec[2][2]*double(mz));
         gbabs2 = dsq(gbx)+dsq(gby)+dsq(gbz);
         gbabs = sqrt(gbabs2);
         if (gbabs > gce2) { continue; }
         // Make cossum and sinsum
         double cossum = 0.0;
         double sinsum = 0.0;
         for (int i = 1; i <= atom.natom/2; i++) {
            for (int j = 0; j < 2; j++) {
               double gbr = shell.tmprx[i][j]*gbx+shell.tmpry[i][j]*gby+shell.tmprz[i][j]*gbz; 
               cossum += shell.q[ksp(atom.asp[i])][j]*cos(gbr);
               sinsum += shell.q[ksp(atom.asp[i])][j]*sin(gbr);
            }
         }
         // Force (N) 
         for (int i = 1; i <= atom.natom/2; i++) {
            for (int j = 0; j < 2; j++) {
               double gbr = shell.tmprx[i][j]*gbx+shell.tmpry[i][j]*gby+shell.tmprz[i][j]*gbz; 
               f = 4.0*M_PI/cell.volume/gbabs*exp(-gbabs2/4.0/dsq(ewaldgm))
                    *shell.q[ksp(atom.asp[i])][j]*(sin(gbr)*cossum-cos(gbr)*sinsum)*shell.efpe;
               fxx = f/gbabs*gbx;
               fyy = f/gbabs*gby;
               fzz = f/gbabs*gbz;
               shell.tmpfx[i][j] += fxx;
               shell.tmpfy[i][j] += fyy;
               shell.tmpfz[i][j] += fzz;
            }
         }
         // Stress
         ew2a = 2.0*M_PI/cell.volume*shell.efpe*exp(-gbabs2/4.0/dsq(ewaldgm))/gbabs2*(dsq(cossum)+dsq(sinsum));
         ew2b = 2.0/gbabs2*(gbabs2/4.0/dsq(ewaldgm)+1.0);
         cell.dmat[0][0] -= ew2a*(ew2b*gbx*gbx-1.0);
         cell.dmat[0][1] -= ew2a*ew2b*gbx*gby;
         cell.dmat[0][2] -= ew2a*ew2b*gbx*gbz;
         cell.dmat[1][0]  = cell.dmat[0][1];
         cell.dmat[1][1] -= ew2a*(ew2b*gby*gby-1.0);
         cell.dmat[1][2] -= ew2a*ew2b*gby*gbz;
         cell.dmat[2][0]  = cell.dmat[0][2];
         cell.dmat[2][1]  = cell.dmat[1][2];
         cell.dmat[2][2] -= ew2a*(ew2b*gbz*gbz-1.0);
         // Energy (J)
         pewaldeach2 = 2.0*M_PI/cell.volume*shell.efpe*exp(-gbabs2/4.0/dsq(ewaldgm))/gbabs2*(dsq(cossum)+dsq(sinsum));
         pewald2 += pewaldeach2;
   }
   }
   }//Ewald 2nd term (Reciprocal space) end
   //Ewald 3rd term
   for (int i = 1; i <= atom.natom/2; i++) {
      for (int j = 0; j < 2; j++) {
         pewald3 -= dsq(shell.q[ksp(atom.asp[i])][j])*ewaldgm/sqrt(M_PI)*shell.efpe;
      }
   }//Ewald 3rd term end
   //Ewald 4th term(Correction) and Core-Shell interaction(Spring)
   for (int i = 1; i <= atom.natom/2; i++) {
      double xx = shell.tmprx[i][0] - shell.tmprx[i][1];
      double yy = shell.tmpry[i][0] - shell.tmpry[i][1];
      double zz = shell.tmprz[i][0] - shell.tmprz[i][1];
      r2 = dsq(xx)+dsq(yy)+dsq(zz);
      r = sqrt(r2);
      // Ewald 4th term
      // Energy (J)
      if (r < 1.0e-11) {
         pewaldeach4 = -shell.q[ksp(atom.asp[i])][0]*shell.q[ksp(atom.asp[i])][1]*shell.efpe*2.0/sqrt(M_PI)
              *(ewaldgm-r2*dsq(ewaldgm)*ewaldgm/3.0+dsq(r2)*dsq(ewaldgm)*dsq(ewaldgm)*ewaldgm/10.0
              -dsq(r2)*r2*dsq(dsq(ewaldgm)*ewaldgm)*ewaldgm/42.0);
      } else {
         pewaldeach4 = -shell.q[ksp(atom.asp[i])][0]*shell.q[ksp(atom.asp[i])][1]*shell.efpe*(1.0-erfc(ewaldgm*r))/r;
      }
      pewald4 += pewaldeach4;
      // Force (N)
      if (r == 0.0) {
         continue;
      } else if (r < 1.0e-11) {
         f = shell.q[ksp(atom.asp[i])][0]*shell.q[ksp(atom.asp[i])][1]*shell.efpe*2.0/sqrt(M_PI)*(-2.0/3.0*dsq(ewaldgm)*ewaldgm*r
              +2.0/5.0*pow(ewaldgm,5)*r2*r-1.0/7.0*pow(ewaldgm,7)*dsq(r2)*r
              +1.0/27.0*pow(ewaldgm,9)*dsq(r2*r)*r);
      } else {
         f = shell.q[ksp(atom.asp[i])][0]*shell.q[ksp(atom.asp[i])][1]*shell.efpe*(2.0*ewaldgm/sqrt(M_PI)
              *exp(-dsq(ewaldgm*r))/r-(1.0-erfc(ewaldgm*r))/r2);
      }
      fxx = f/r*xx;
      fyy = f/r*yy;
      fzz = f/r*zz;
      shell.tmpfx[i][0] += fxx;
      shell.tmpfy[i][0] += fyy;
      shell.tmpfz[i][0] += fzz;
      shell.tmpfx[i][1] -= fxx;
      shell.tmpfy[i][1] -= fyy;
      shell.tmpfz[i][1] -= fzz;
      cell.dmat[0][0] += fxx*xx;
      cell.dmat[0][1] += fxx*yy;
      cell.dmat[0][2] += fxx*zz;
      cell.dmat[1][0] = cell.dmat[0][1];
      cell.dmat[1][1] += fyy*yy;
      cell.dmat[1][2] += fyy*zz;
      cell.dmat[2][0] = cell.dmat[0][2];
      cell.dmat[2][1] = cell.dmat[1][2];
      cell.dmat[2][2] += fzz*zz;
      // Core-Shell interaction(Spring)
      if (r == 0.0) { continue; }
      double pheach = (0.5*shell.springc[ksp(atom.asp[i])]*(r2*1.0e20)+1.0/24.0*shell.springc2[ksp(atom.asp[i])]*(dsq(r2)*1.0e40))*eV;
      pharmonic += pheach;
      f =-(shell.springc[ksp(atom.asp[i])]*(r*1.0e10)+1.0/6.0*shell.springc2[ksp(atom.asp[i])]*(r2*r*1.0e30))*eV*1.0e10;
      fxx = f/r*xx;
      fyy = f/r*yy;
      fzz = f/r*zz;
      shell.tmpfx[i][0] += fxx;
      shell.tmpfy[i][0] += fyy;
      shell.tmpfz[i][0] += fzz;
      shell.tmpfx[i][1] -= fxx;
      shell.tmpfy[i][1] -= fyy;
      shell.tmpfz[i][1] -= fzz;
      cell.dmat[0][0] += fxx*xx;
      cell.dmat[0][1] += fxx*yy;
      cell.dmat[0][2] += fxx*zz;
      cell.dmat[1][0]  = cell.dmat[0][1];
      cell.dmat[1][1] += fyy*yy;
      cell.dmat[1][2] += fyy*zz;
      cell.dmat[2][0]  = cell.dmat[0][2];
      cell.dmat[2][1]  = cell.dmat[1][2];
      cell.dmat[2][2] += fzz*zz;
   } //Ewald 4th term(Correction) and Core-Shell interaction(Spring) end

   //---Etot---
   atom.epotsum = pewald1+pewald2+pewald3+pewald4+psr+pharmonic;
   cell.virx = cell.dmat[0][0]; cell.viry = cell.dmat[1][1]; cell.virz = cell.dmat[2][2];
   for (int i = 1; i <= atom.natom/2; i++) {
         atom.fx[i] = shell.tmpfx[i][0];
         atom.fy[i] = shell.tmpfy[i][0];
         atom.fz[i] = shell.tmpfz[i][0];
         atom.fx[i+atom.natom/2] = shell.tmpfx[i][1];
         atom.fy[i+atom.natom/2] = shell.tmpfy[i][1];
         atom.fz[i+atom.natom/2] = shell.tmpfz[i][1];
   }
}

int ksp(char* species)
{
   int k;
   if      (strcmp(species, "Pb") == 0) { k = 0; }
   else if (strcmp(species, "O") == 0)  { k = 1; }
   else if (strcmp(species, "Ti") == 0) { k = 2; }
   else {
      printf("Shell Model Error: Atom %s is not supported in Shell Model.\n", species); 
      k = 0; return k; }
   return k;
}

/*
double erfc_(double x)
{
   double y = std::erf(x);
   return 1.0 - y;
}
*/

void shell_delete()
{
   if (shell.tmprx == NULL) {return;}
   for (int i = 0; i <= atom.natom; i++) {
      printf("%d\n",i);
      delete[] shell.tmprx[i];
      delete[] shell.tmpry[i];
      delete[] shell.tmprz[i];
      delete[] shell.tmpfx[i];
      delete[] shell.tmpfy[i];
      delete[] shell.tmpfz[i];
   }
   delete[] shell.tmprx; shell.tmprx = NULL;
   delete[] shell.tmpry; shell.tmpry = NULL;
   delete[] shell.tmprz; shell.tmprz = NULL;
   delete[] shell.tmpfx; shell.tmpfx = NULL;
   delete[] shell.tmpfy; shell.tmpfy = NULL;
   delete[] shell.tmpfz; shell.tmpfz = NULL;
}
void shell_alloc()
{
   shell.tmprx = new double*[atom.natom+1];
   shell.tmpry = new double*[atom.natom+1];
   shell.tmprz = new double*[atom.natom+1];
   shell.tmpfx = new double*[atom.natom+1];
   shell.tmpfy = new double*[atom.natom+1];
   shell.tmpfz = new double*[atom.natom+1];
   for (int i = 0; i <= atom.natom; i++) {
      shell.tmprx[i] = new double[2];
      shell.tmpry[i] = new double[2];
      shell.tmprz[i] = new double[2];
      shell.tmpfx[i] = new double[2];
      shell.tmpfy[i] = new double[2];
      shell.tmpfz[i] = new double[2];
   }
}