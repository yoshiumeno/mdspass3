#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#if defined(_WIN32) || defined(_WIN64)
#include <Windows.h>
#endif
#include <GL/gl.h>
#include <complex>
#include "myheader.h"

using namespace std;

extern float bondlength;
double anint(double x);
void e_force_brenner(void);

void bond_set()
{
  if (book.algo == 3) {
    //std::cout<<"Setting bond using BK algo #3 (Brenner's BK scheme)"<<std::endl;
    //if (!bre.nabors) { e_force_brenner(); } // Call Brenner's BK if not yet called
    //bre.initialize = true;
    //e_force_brenner();
    float bondlength2 = bondlength*bondlength;
    double cube1=cell.hmat[0][0]/ang;
    double cube2=cell.hmat[1][1]/ang;
    double cube3=cell.hmat[2][2]/ang;
    int totalbond = 0;
    for (int i=1; i<=bre.np; i++) {  // bre.np = atom.natom
      int jbegin = bre.nabors[i]; int jend = bre.nabors[i+1]-1;
      //if (jbegin > jend) { continue; }
      if (jbegin > jend) { atom.nneighbor[i] = 0; continue; } // corrected 20191202
      atom.nneighbor[i] = 0;
      for (int j=jbegin; j<=jend; j++) {
	//if (bre.lcheck[j]!=1) { continue; } <<= Not needed!!
	int jn = bre.list[j];
	//if (i >= jn) { continue; }
	float rr1 = bre.r01[jn]-bre.r01[i];
	float rr2 = bre.r02[jn]-bre.r02[i];
	float rr3 = bre.r03[jn]-bre.r03[i];
	//if (cell.pbcx) { rr1 = rr1-cube1*anint(rr1/cube1); }
	//if (cell.pbcy) { rr2 = rr2-cube2*anint(rr2/cube2); }
	//if (cell.pbcz) { rr3 = rr3-cube3*anint(rr3/cube3); }
	double xx=0,yy=0,zz=0;
	/*
	if (cell.pbcx) { xx = -(double)anint(rr1/cube1); }
	if (cell.pbcy) { yy = -(double)anint(rr2/cube2); }
	if (cell.pbcz) { zz = -(double)anint(rr3/cube3); }
	*/
	xx = (double)bre.rep1[j]; yy = (double)bre.rep2[j]; zz = (double)bre.rep3[j];
	rr1 -= xx*cell.hmat[0][0]/ang + yy*cell.hmat[0][1]/ang + zz*cell.hmat[0][2]/ang;
	rr2 -= xx*cell.hmat[1][0]/ang + yy*cell.hmat[1][1]/ang + zz*cell.hmat[1][2]/ang;
	rr3 -= xx*cell.hmat[2][0]/ang + yy*cell.hmat[2][1]/ang + zz*cell.hmat[2][2]/ang;
	float rsq = rr1*rr1+rr2*rr2+rr3*rr3;
	if (rsq < bondlength2) {
	  if (atom.nneighbor[i] == MAXNEIGHBOR-1) {
	    std::cout<<"MAXNEIGHBOR is too small"<<std::endl; break; }
	  atom.neighbor[i][atom.nneighbor[i]] = jn;
	  atom.nneighbor[i] = atom.nneighbor[i] + 1;
	  totalbond++;
	}
      }
    }
    //printf("Number of total bonds = %d\n",totalbond);
  } else if (book.algo == 1) {
    float bondlength2 = bondlength*bondlength*1e-20;
    int totalbond = 0;
    for (int i=1; i<=atom.natom; i++) {
      atom.nneighbor[i] = 0;
      for (int j=1; j<=book.alistnum[i]; j++) {
	int jn = book.alist[i][j][0];
	int ix = book.alist[i][j][1];
	int iy = book.alist[i][j][2];
	int iz = book.alist[i][j][3];
	double rr2 = atom.Dist2(i,jn,ix,iy,iz);
	if (rr2 < bondlength2) {
	  if (atom.nneighbor[i] == MAXNEIGHBOR-1) {
	    std::cout<<"MAXNEIGHBOR is too small"<<std::endl; break; }
	  atom.neighbor[i][atom.nneighbor[i]] = jn;
	  atom.nneighbor[i] = atom.nneighbor[i] + 1;
	  totalbond++;
	}
      }
    }
    //printf("Number of total bonds = %d\n",totalbond);
  }

}
