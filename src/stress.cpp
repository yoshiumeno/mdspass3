#include <iostream>
#include <stdio.h>
#include "myheader.h"

void restore_pos(double *hmat_p);
void stretch(double x, double y, double z);
void potential();
void strain(double strmat[3][3]);
void resetmat(double strmat[3][3]);
void calcstresstensor();
void inverse(double a[3][3], double b[3][3]);
void matmul(double a[3][3], double b[3][3], double c[3][3]);

extern float strs_set_xx,strs_set_yy,strs_set_zz,strs_set_xy,strs_set_yz,strs_set_zx;

void calcstresstensor()
{
  double pt[3][3], tm[3][3];
  double xx, yy, zz;
  inverse(cell.hmat,cell.hinmat);
  cell.volume = cell.Getvolume();
  resetmat(pt);
  matmul(cell.hvmat,cell.hinmat,tm);
  for (int i=1; i<=atom.natom; i++) {
    xx = atom.vx[i] - (tm[0][0]*atom.rx[i]+tm[0][1]*atom.ry[i]+tm[0][2]*atom.rz[i]);
    yy = atom.vy[i] - (tm[1][0]*atom.rx[i]+tm[1][1]*atom.ry[i]+tm[1][2]*atom.rz[i]);
    zz = atom.vz[i] - (tm[2][0]*atom.rx[i]+tm[2][1]*atom.ry[i]+tm[2][2]*atom.rz[i]);
    pt[0][0] += xx * xx * atom.wm[i];
    pt[1][0] += yy * xx * atom.wm[i];
    pt[2][0] += zz * xx * atom.wm[i];
    pt[0][1] += xx * yy * atom.wm[i];
    pt[1][1] += yy * yy * atom.wm[i];
    pt[2][1] += zz * yy * atom.wm[i];
    pt[0][2] += xx * zz * atom.wm[i];
    pt[1][2] += yy * zz * atom.wm[i];
    pt[2][2] += zz * zz * atom.wm[i];
  }
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      pt[i][j] += cell.dmat[i][j];
      pt[i][j] /= cell.volume;
      cell.sgmmat[i][j] = -pt[i][j];
    } }
}

void stresscheck(double eps)
{
  double ene0,ene1,sxx,syy,szz,sxy,syz,szx;
  //double eps = 1.0e-3;
  double hmat_p[9];
  double strmat[3][3];
  printf("######### Checking global stress calculation #########\n");
  // Keep original cell and position
  for (int i=0; i<3; i++) { for (int j=0; j<3; j++) {
      hmat_p[3*i+j] = cell.hmat[i][j]; } }
  for (int i=1; i<=atom.natom; i++) {
    atom.rx_p[i] = atom.rx[i];
    atom.ry_p[i] = atom.ry[i];
    atom.rz_p[i] = atom.rz[i]; }
  cell.volume = cell.Getvolume();

  /*
  stretch(1+eps, 1.0, 1.0); potential(); ene0 = atom.epotsum; restore_pos(hmat_p);
  stretch(1-eps, 1.0, 1.0); potential(); ene1 = atom.epotsum; restore_pos(hmat_p);
  sxx=(ene0-ene1)/(eps*2)/cell.volume*1e-6;
  stretch(1.0, 1+eps, 1.0); potential(); ene0 = atom.epotsum; restore_pos(hmat_p);
  stretch(1.0, 1-eps, 1.0); potential(); ene1 = atom.epotsum; restore_pos(hmat_p);
  syy=(ene0-ene1)/(eps*2)/cell.volume*1e-6;
  stretch(1.0, 1.0, 1+eps); potential(); ene0 = atom.epotsum; restore_pos(hmat_p);
  stretch(1.0, 1.0, 1-eps); potential(); ene1 = atom.epotsum; restore_pos(hmat_p);
  szz=(ene0-ene1)/(eps*2)/cell.volume*1e-6;
  printf("%f %f %f (MPa)\n",sxx,syy,szz);
  */
  
  resetmat(strmat);strmat[0][0]= eps;strain(strmat);
  potential(); ene0 = atom.epotsum; restore_pos(hmat_p);
  resetmat(strmat);strmat[0][0]=-eps;strain(strmat);
  potential(); ene1 = atom.epotsum; restore_pos(hmat_p);
  sxx=(ene0-ene1)/(eps*2)/cell.volume*1e-6;
  resetmat(strmat);strmat[1][1]= eps;strain(strmat);
  potential(); ene0 = atom.epotsum; restore_pos(hmat_p);
  resetmat(strmat);strmat[1][1]=-eps;strain(strmat);
  potential(); ene1 = atom.epotsum; restore_pos(hmat_p);
  syy=(ene0-ene1)/(eps*2)/cell.volume*1e-6;
  resetmat(strmat);strmat[2][2]= eps;strain(strmat);
  potential(); ene0 = atom.epotsum; restore_pos(hmat_p);
  resetmat(strmat);strmat[2][2]=-eps;strain(strmat);
  potential(); ene1 = atom.epotsum; restore_pos(hmat_p);
  szz=(ene0-ene1)/(eps*2)/cell.volume*1e-6;
  printf(" xx,yy,zz (numerical)  %10.5f %10.5f %10.5f (MPa)\n",sxx,syy,szz);

  potential();calcstresstensor();
  printf(" xx,yy,zz (analytical) %10.5f %10.5f %10.5f (MPa)\n",
	 cell.sgmmat[0][0]*1e-6, cell.sgmmat[1][1]*1e-6, cell.sgmmat[2][2]*1e-6);

  resetmat(strmat);strmat[0][1]= eps;strain(strmat);
  potential(); ene0 = atom.epotsum; restore_pos(hmat_p);
  resetmat(strmat);strmat[0][1]=-eps;strain(strmat);
  potential(); ene1 = atom.epotsum; restore_pos(hmat_p);
  sxy=(ene0-ene1)/(eps*2)/cell.volume*1e-6;
  resetmat(strmat);strmat[1][2]= eps;strain(strmat);
  potential(); ene0 = atom.epotsum; restore_pos(hmat_p);
  resetmat(strmat);strmat[1][2]=-eps;strain(strmat);
  potential(); ene1 = atom.epotsum; restore_pos(hmat_p);
  syz=(ene0-ene1)/(eps*2)/cell.volume*1e-6;
  resetmat(strmat);strmat[2][0]= eps;strain(strmat);
  potential(); ene0 = atom.epotsum; restore_pos(hmat_p);
  resetmat(strmat);strmat[2][0]=-eps;strain(strmat);
  potential(); ene1 = atom.epotsum; restore_pos(hmat_p);
  szx=(ene0-ene1)/(eps*2)/cell.volume*1e-6;
  printf(" xy,yz,zx (numerical)  %10.5f %10.5f %10.5f (MPa)\n",sxy,syz,szx);

  potential();calcstresstensor();
  printf(" xy,yz,zx (analytical) %10.5f %10.5f %10.5f (MPa)\n",
	 cell.sgmmat[0][1]*1e-6, cell.sgmmat[1][2]*1e-6, cell.sgmmat[2][0]*1e-6);
  
}

void restore_pos(double *hmat_p)
{
  for (int i=0; i<3; i++) { for (int j=0; j<3; j++) {
      cell.hmat[i][j] = hmat_p[3*i+j]; } }
  for (int i=1; i<=atom.natom; i++) {
    atom.rx[i] = atom.rx_p[i];
    atom.ry[i] = atom.ry_p[i];
    atom.rz[i] = atom.rz_p[i]; }
}

void stress_set()
{
  cell.sgmmat_set[0][0] = strs_set_xx*1e6;
  cell.sgmmat_set[1][1] = strs_set_yy*1e6;
  cell.sgmmat_set[2][2] = strs_set_zz*1e6;
  cell.sgmmat_set[0][1] = strs_set_xy*1e6; cell.sgmmat_set[1][0] = strs_set_xy*1e6;
  cell.sgmmat_set[1][2] = strs_set_yz*1e6; cell.sgmmat_set[2][1] = strs_set_yz*1e6;
  cell.sgmmat_set[2][0] = strs_set_zx*1e6; cell.sgmmat_set[0][2] = strs_set_zx*1e6;
}

