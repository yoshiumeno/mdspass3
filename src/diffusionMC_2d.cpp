#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "myheader.h"

#if !defined __linux__ && !defined __APPLE__
#define snprintf sprintf_s
#endif
#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

void out_of_cell();

void diffusionMC_2d_initialize()
{
  if (dmc2d.pattern) { 
    for (int i=0; i<dmc2d.pattern_x; i++) { delete[] dmc2d.pattern[i]; }
    delete[] dmc2d.pattern; dmc2d.pattern = NULL; }
  int ix, iy;
  int ii;
  printf("Pattern of phases is read from PATTERN:\n");
  FILE *fp = fopen("PATTERN","r");
  fscanf(fp, "%d %d\n", &ix, &iy);
  dmc2d.pattern_x = ix; dmc2d.pattern_y = iy;
  printf("Pattern size (X x Y) = %d x %d\n",ix,iy);
  dmc2d.pattern = new int*[ix];
  for (int i=0; i<ix; i++) { dmc2d.pattern[i] = new int[iy]; }
  for (int j=0; j<iy; j++) {
    for (int i=0; i<ix; i++) {
      fscanf(fp, "%d ", &dmc2d.pattern[i][j]);
    }
  }
  if (true) {
    for (int j=0; j<iy; j++) {
      for (int i=0; i<ix; i++) {
	printf("%d ", dmc2d.pattern[i][j]);
      } printf("\n");
    }
    printf("\n");
  }

  fclose(fp);
  dmc2d.coef[0] = 1.5;
  dmc2d.coef[1] = 0.5;
  incell = 1;
  for (int i=1; i<=atom.natom; i++) {
    // used for trajectory
    atom.ax[i] = atom.rx[i];
    atom.ay[i] = atom.ry[i];
    atom.az[i] = atom.rz[i];
  }

}

void diffusionMC_2d()
{
  if (istep==0) { diffusionMC_2d_initialize(); }
  double dx, dy, dr, theta;
  double alpha = 1e4;
  double diffusion_coefficient = 2.0;
  for (int i=1; i<=atom.natom; i++) {
    double xx = atom.rx[i]/cell.hmat[0][0];
    double yy = atom.ry[i]/cell.hmat[1][1];
    int ix = (int)((double)dmc2d.pattern_x*xx);
    int iy = (int)((double)dmc2d.pattern_y*yy);
    if (ix>=dmc2d.pattern_x) { ix = dmc2d.pattern_x-1; }
    if (iy>=dmc2d.pattern_x) { iy = dmc2d.pattern_y-1; }
    if (ix<0) { ix = 0; }
    if (iy<0) { iy = 0; }
    diffusion_coefficient = dmc2d.coef[dmc2d.pattern[ix][iy]];
    dr = alpha * dt * sqrt(diffusion_coefficient);
    theta = 360.0*(double)rand()/(double)RAND_MAX;
    dx = dr * cos(M_PI/180.0*theta);
    dy = dr * sin(M_PI/180.0*theta);
    atom.rx[i] += dx; atom.ry[i] += dy;
    atom.ax[i] += dx; atom.ay[i] += dy;
    atom.vx[i] = dx/dt; atom.vy[i] = dy/dt; atom.vz[i] = 0.0;
    if (incell) { out_of_cell(); }
  }

}
