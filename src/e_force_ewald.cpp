//
// FORCE & ENERGY CALCULATION using Ewald Method
// 2024.07
//

#include <iostream>
#include<string.h>
#include<stdlib.h>
#include<fstream>
#include "myheader.h"
#define dsq(a) ((a)*(a))
void resetmat(double a[3][3]);
double erfc(double x);
void recipro(double h[3][3], double r[3][3]);
int ksp(char* species);

void e_force_ewald()
{
   if (ew.initialize) {
      ew.initialize = false;
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
   double acc = pow(1.0e-1, double(ew.naccew));
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
   double rc = ew.rcsr * 1.0e-10;
   rc = std::max(rce1, rc);
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
   double pewald0 = 0.0; double pewald1 = 0.0; double pewald2 = 0.0; double pewald3 = 0.0;
   double psr = 0.0; double pharmonic = 0.0;
   //----Ewald 1st term and Short-range interaction
   int j, ix, iy, iz; double rr2, r, rr;
   for (int i = 1; i <= atom.natom; i++) {
      for (int kk = 1; kk <= book.alistnum[i]; kk++) {
         j = book.alist[i][kk][0];
         if (j < i) { continue; }
   	   ix = book.alist[i][kk][1]; iy = book.alist[i][kk][2]; iz = book.alist[i][kk][3];
	      //rr2 = atom.Dist2(i,j,ix,iy,iz);
         //for (int k = 0; k < 2; k++) {
         //   for (int l = 0; l < 2; l++) {
         //double xx = ew.tmprx[i][k] - ew.tmprx[j][k];
         //double yy = ew.tmpry[i][k] - ew.tmpry[j][k];
         //double zz = ew.tmprz[i][k] - ew.tmprz[j][k];
         double xx = atom.rx[i] - atom.rx[j];
         double yy = atom.ry[i] - atom.ry[j];
         double zz = atom.rz[i] - atom.rz[j];
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
               pewaldeach1 = ew.q[ksp(atom.asp[i])]*ew.q[ksp(atom.asp[j])]*ew.efpe*erfc(ewaldgm*r)/r;
               // Force(N)
               f = ew.q[ksp(atom.asp[i])]*ew.q[ksp(atom.asp[j])]*ew.efpe*(erfc(ewaldgm*r)/r2
                          +2.0*ewaldgm/sqrt(M_PI)*exp(-dsq(ewaldgm*r))/r);
               fxx = f/r*xx;
               fyy = f/r*yy;
               fzz = f/r*zz;
               if (i == j) {
                  atom.fx[i] += fxx;
                  atom.fy[i] += fyy;
                  atom.fz[i] += fzz;
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
                  atom.fx[i] += fxx;
                  atom.fy[i] += fyy;
                  atom.fz[i] += fzz;
                  atom.fx[j] -= fxx;
                  atom.fy[j] -= fyy;
                  atom.fz[j] -= fzz;                       
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
            /*
            bool if_short = true;
            if (kint[ksp(atom.asp[i])][ksp(atom.asp[j])]==0) { if_short = false; }
            if (k != 1) { if_short = false; }
            if (l != 1) { if_short = false; }
            int kitr = kint[ksp(atom.asp[i])][ksp(atom.asp[j])];
            if (rr > ew.rcsr[kitr]) {if_short = false; }
            if (if_short) {
               f = (ew.alarge[kitr]/ew.rho[kitr]*exp(-rr/ew.rho[kitr])-6.0*ew.clarge[ kitr]/pow(rr,7))*ev*1.0e10;
               psreach = (ew.alarge[kitr]*exp(-rr/ew.rho[kitr])-ew.clarge[ kitr]/pow(rr,6))*ev;
               // Tapering (rt < rr < rc)
               if ((rr > ew.rcsr[kitr]-ew.rtp[kitr])&&(rr <= ew.rcsr[kitr])) {
                  psreach = psreach*0.50
                       *(1.0+cos(M_PI*(rr-ew.rcsr[kitr]-ew.rtp[kitr])/(ew.rcsr[kitr]-ew.rcsr[kitr]-ew.rtp[kitr])));
                  f = ((ew.alarge[kitr]/ew.rho[kitr]*exp(-rr/ew.rho[kitr])-6.0*ew.clarge[ kitr]/pow(rr,7))
                       *0.50*(1.0+cos(M_PI*(rr-(ew.rcsr[kitr]-ew.rtp[kitr]))/ew.rtp[kitr]))
                       -(ew.alarge[kitr]*exp(-rr/ew.rho[kitr])-ew.clarge[ kitr]/pow(rr,6))
                       *(-0.50*M_PI/ew.rtp[kitr]*sin(M_PI*(rr-(ew.rcsr[kitr]-ew.rtp[kitr]))/ew.rtp[kitr])))*ev*1.0e10;
               }
               fxx=f/r*xx;
               fyy=f/r*yy;
               fzz=f/r*zz;
               if (i == j) {
                  ew.tmpfx[i][1] += fxx;
                  ew.tmpfy[i][1] += fyy;
                  ew.tmpfz[i][1] += fzz;
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
                  ew.tmpfx[i][1] += fxx;
                  ew.tmpfy[i][1] += fyy;
                  ew.tmpfz[i][1] += fzz;
                  ew.tmpfx[j][1] -= fxx;
                  ew.tmpfy[j][1] -= fyy;
                  ew.tmpfz[j][1] -= fzz;
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
            */
         } //if istat == 1 end
            //}
         //}
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
      for (int i = 1; i <= atom.natom; i++) {
         //for (int j = 0; j < 2; j++) {
         double gbr = atom.rx[i]*gbx+atom.ry[i]*gby+atom.rz[i]*gbz; 
         cossum += ew.q[ksp(atom.asp[i])]*cos(gbr);
         sinsum += ew.q[ksp(atom.asp[i])]*sin(gbr);
         //}
      }
      // Force (N) 
      for (int i = 1; i <= atom.natom; i++) {
         //for (int j = 0; j < 2; j++) {
         double gbr = atom.rx[i]*gbx+atom.ry[i]*gby+atom.rz[i]*gbz; 
         f = 4.0*M_PI/cell.volume/gbabs*exp(-gbabs2/4.0/dsq(ewaldgm))
              *ew.q[ksp(atom.asp[i])]*(sin(gbr)*cossum-cos(gbr)*sinsum)*ew.efpe;
         fxx = f/gbabs*gbx;
         fyy = f/gbabs*gby;
         fzz = f/gbabs*gbz;
         atom.fx[i] += fxx;
         atom.fy[i] += fyy;
         atom.fz[i] += fzz;
         //}
      }
      // Stress
      ew2a = 2.0*M_PI/cell.volume*ew.efpe*exp(-gbabs2/4.0/dsq(ewaldgm))/gbabs2*(dsq(cossum)+dsq(sinsum));
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
      pewaldeach2 = 2.0*M_PI/cell.volume*ew.efpe*exp(-gbabs2/4.0/dsq(ewaldgm))/gbabs2*(dsq(cossum)+dsq(sinsum));
      pewald2 += pewaldeach2;
   }
   }
   }//Ewald 2nd term (Reciprocal space) end
   //Ewald 3rd term
   for (int i = 1; i <= atom.natom; i++) {
      //for (int j = 0; j < 2; j++) {
      pewald3 -= dsq(ew.q[ksp(atom.asp[i])])*ewaldgm/sqrt(M_PI)*ew.efpe;
      //}
   }//Ewald 3rd term end
   //Ewald 4th term(Correction) and Core-shell interaction(Spring)
   /*
   for (int i = 1; i <= atom.natom; i++) {
      xx = ew.tmprx[i][0] - ew.tmprx[i][1];
      yy = ew.tmpry[i][0] - ew.tmpry[i][1];
      zz = ew.tmprz[i][0] - ew.tmprz[i][1];
      r2 = dsq(xx)+dsq(yy)+dsq(zz);
      r = sqrt(r2);
      // Ewald 4th term
      // Energy (J)
      if (r < 1.0e-11) {
         pewaldeach4 = -ew.q[ksp(atom.asp[i])][0]*ew.q[ksp(atom.asp[i])][1]*ew.efpe*2.0/sqrt(M_PI)
              *(ewaldgm-r2*dsq(ewaldgm)*ewaldgm/3.0+dsq(r2)*dsq(ewaldgm)*dsq(ewaldgm)*ewaldgm/10.0
              -dsq(r2)*r2*dsq(dsq(ewaldgm)*ewaldgm)*ewaldgm/42.0);
      } else {
         pewaldeach4 = -ew.q[ksp(atom.asp[i])][0]*ew.q[ksp(atom.asp[i])][1]*ew.efpe*(1.0-erfc(ewaldgm*r))/r;
      }
      pewald4 += pewaldeach4;
      // Force (N)
      if (r == 0.0) {
         continue;
      } elseif (r < 1.0e-11) {
         f = ew.q[ksp(atom.asp[i])][0]*ew.q[ksp(atom.asp[i])][1]*ew.efpe*2.0/sqrt(M_PI)*(-2.0/3.0*dsq(ewaldgm)*ewaldgm*r
              +2.0/5.0*pow(ewaldgm,5)*r2*r-1.0/7.0*pow(ewaldgm,7)*dsq(r2)*r
              +1.0/27.0*pow(ewaldgm,9)*dsq(r2*r)*r);
      } else {
         f = ew.q[ksp(atom.asp[i])][0]*ew.q[ksp(atom.asp[i])][1]*efpe*(2.0*ewaldgm/sqrt(M_PI)
              *exp(-dsq(ewaldgm*r))/r-(1.0-erfc(ewaldgm*r))/r2);
      }
      fxx = f/r*xx;
      fyy = f/r*yy;
      fzz = f/r*zz;
      ew.tmpfx[i][0] += fxx;
      ew.tmpfy[i][0] += fyy;
      ew.tmpfz[i][0] += fzz;
      ew.tmpfx[i][1] -= fxx;
      ew.tmpfy[i][1] -= fyy;
      ew.tmpfz[i][1] -= fzz;
      cell.dmat[0][0] += fxx*xx;
      cell.dmat[0][1] += fxx*yy;
      cell.dmat[0][2] += fxx*zz;
      cell.dmat[1][0] = cell.dmat[0][1];
      cell.dmat[1][1] += fyy*yy;
      cell.dmat[1][2] += fyy*zz;
      cell.dmat[2][0] = cell.dmat[0][2];
      cell.dmat[2][1] = cell.dmat[1][2];
      cell.dmat[2][2] += fzz*zz;
      // Core-shell interaction(Spring)
      if (r == 0.0) { continue; }
      pheach = (0.5*ew.springc[ksp(atom.asp[i])]*(r2*1.0e20)+1.0/24.0*ew.springc2[ksp(atom.asp[i])]*(dsq(r2)*1.0e40))*ev;
      pharmonic += pheach;
      f =-(ew.springc[ksp(atom.asp[i])]*(r*1.0e10)+1.0/6.0*ew.springc2[ksp(atom.asp[i])]*(r2*r*1.0e30))*ev*1.0e10;
      fxx = f/r*xx;
      fyy = f/r*yy;
      fzz = f/r*zz;
      ew.tmpfx[i][0] += fxx;
      ew.tmpfy[i][0] += fyy;
      ew.tmpfz[i][0] += fzz;
      ew.tmpfx[i][1] -= fxx;
      ew.tmpfy[i][1] -= fyy;
      ew.tmpfz[i][1] -= fzz;
      cell.dmat[0][0] += fxx*xx;
      cell.dmat[0][1] += fxx*yy;
      cell.dmat[0][2] += fxx*zz;
      cell.dmat[1][0]  = cell.dmat[0][1];
      cell.dmat[1][1] += fyy*yy;
      cell.dmat[1][2] += fyy*zz;
      cell.dmat[2][0]  = cell.dmat[0][2];
      cell.dmat[2][1]  = cell.dmat[1][2];
      cell.dmat[2][2] += fzz*zz;
   } //Ewald 4th term(Correction) and Core-shell interaction(Spring) end
   */

   //---Etot---
   //enpot = pewald1+pewald2+pewald3+pewald4+psr+pharmonic;
   atom.epotsum = pewald1+pewald2+pewald3+psr+pharmonic;
   cell.virx = cell.dmat[0][0]; cell.viry = cell.dmat[1][1]; cell.virz = cell.dmat[2][2];
}

/*
int ksp(char* species)
{
   int k;
   if      (strcmp(species, "Pb") == 0) { k = 0; }
   else if (strcmp(species, "O") == 0) { k = 1; }
   else if (strcmp(species, "Ti") == 0) { k = 2; }
   else {
      printf("Ewakd Error: Atom %s i not supported in Shell Model.\n", species); 
      k = 0; return k; }
   return k;
}
*/

/*
double erfc_(double x)
{
   double y = std::erf(x);
   return 1.0 - y;
}
*/

