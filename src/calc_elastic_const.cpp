#include <iostream>
#include <fstream>
#include "myheader.h"

void matcpy(double a[3][3], double b[3][3]);
void resetmat(double a[3][3]);
void resetmat6(double a[6][6]);
void strain(double strmat[3][3]);
void potential();
void calcstresstensor();
void writedata_elastic_const();

void calc_elastic_const()
{
  double strmat_store[3][3], hmat_store[3][3], strmat[3][3];
  double *rx_store, *ry_store, *rz_store;
  double inf_str;
  double c11, c12, c13, c44, c55, c66;
  double s1, s2, s3, s4, s5, s6;

  //matcpy(strmat_set, strmat_store);
  resetmat(strmat);
  /*
      strmat_store=strmat_set   ! Strain matrix save
                                ! (should be restored before return)
  */
  rx_store = new double[atom.natom+1];
  ry_store = new double[atom.natom+1];
  rz_store = new double[atom.natom+1];
  matcpy(cell.hmat, hmat_store);
  for (int i=1; i<=atom.natom; i++) {
    rx_store[i] = atom.rx[i]; ry_store[i] = atom.ry[i]; rz_store[i] = atom.rz[i]; }
  
  inf_str=1.0e-5;  // Infinitesimal strain
  //inf_str=econstd
      
  //      elasc=0.0d0
  resetmat6(cell.ecmat);

  //  C_11,C_12,C_13
  matcpy(hmat_store, cell.hmat);
  for (int i=1; i<=atom.natom; i++) {
    atom.rx[i] = rx_store[i]; atom.ry[i] = ry_store[i]; atom.rz[i] = rz_store[i]; }
  resetmat(strmat);
  strmat[0][0]=inf_str;
  strain(strmat); potential(); calcstresstensor();
  s1=cell.sgmmat[0][0]; s3=cell.sgmmat[1][1]; s5=cell.sgmmat[2][2];
  matcpy(hmat_store, cell.hmat);
  for (int i=1; i<=atom.natom; i++) {
    atom.rx[i] = rx_store[i]; atom.ry[i] = ry_store[i]; atom.rz[i] = rz_store[i]; }
 
  matcpy(hmat_store, cell.hmat);
  for (int i=1; i<=atom.natom; i++) {
    atom.rx[i] = rx_store[i]; atom.ry[i] = ry_store[i]; atom.rz[i] = rz_store[i]; }
  resetmat(strmat);
  strmat[0][0]=-inf_str;
  strain(strmat); potential(); calcstresstensor();
  s2=cell.sgmmat[0][0]; s4=cell.sgmmat[1][1]; s6=cell.sgmmat[2][2];
  matcpy(hmat_store, cell.hmat);
  for (int i=1; i<=atom.natom; i++) {
    atom.rx[i] = rx_store[i]; atom.ry[i] = ry_store[i]; atom.rz[i] = rz_store[i]; }

  cell.ecmat[0][0] = (s1 - s2)/(inf_str * 2.0);
  cell.ecmat[0][1] = (s3 - s4)/(inf_str * 2.0);
  cell.ecmat[0][2] = (s5 - s6)/(inf_str * 2.0);
  cell.ecmat[1][0] = cell.ecmat[0][1];
  cell.ecmat[2][0] = cell.ecmat[0][2];

  //  C_22,C_23
  matcpy(hmat_store, cell.hmat);
  for (int i=1; i<=atom.natom; i++) {
    atom.rx[i] = rx_store[i]; atom.ry[i] = ry_store[i]; atom.rz[i] = rz_store[i]; }
  resetmat(strmat);
  strmat[1][1]=inf_str;
  strain(strmat); potential(); calcstresstensor();
  s1=cell.sgmmat[1][1]; s3=cell.sgmmat[2][2];
  matcpy(hmat_store, cell.hmat);
  for (int i=1; i<=atom.natom; i++) {
    atom.rx[i] = rx_store[i]; atom.ry[i] = ry_store[i]; atom.rz[i] = rz_store[i]; }
 
  matcpy(hmat_store, cell.hmat);
  for (int i=1; i<=atom.natom; i++) {
    atom.rx[i] = rx_store[i]; atom.ry[i] = ry_store[i]; atom.rz[i] = rz_store[i]; }
  resetmat(strmat);
  strmat[1][1]=-inf_str;
  strain(strmat); potential(); calcstresstensor();
  s2=cell.sgmmat[1][1]; s4=cell.sgmmat[2][2];
  matcpy(hmat_store, cell.hmat);
  for (int i=1; i<=atom.natom; i++) {
    atom.rx[i] = rx_store[i]; atom.ry[i] = ry_store[i]; atom.rz[i] = rz_store[i]; }

  cell.ecmat[1][1] = (s1 - s2)/(inf_str * 2.0);
  cell.ecmat[1][2] = (s3 - s4)/(inf_str * 2.0);
  cell.ecmat[2][1] = cell.ecmat[1][2];

  //  C_33
  matcpy(hmat_store, cell.hmat);
  for (int i=1; i<=atom.natom; i++) {
    atom.rx[i] = rx_store[i]; atom.ry[i] = ry_store[i]; atom.rz[i] = rz_store[i]; }
  resetmat(strmat);
  strmat[2][2]=inf_str;
  strain(strmat); potential(); calcstresstensor();
  s1=cell.sgmmat[2][2];
  matcpy(hmat_store, cell.hmat);
  for (int i=1; i<=atom.natom; i++) {
    atom.rx[i] = rx_store[i]; atom.ry[i] = ry_store[i]; atom.rz[i] = rz_store[i]; }
 
  matcpy(hmat_store, cell.hmat);
  for (int i=1; i<=atom.natom; i++) {
    atom.rx[i] = rx_store[i]; atom.ry[i] = ry_store[i]; atom.rz[i] = rz_store[i]; }
  resetmat(strmat);
  strmat[2][2]=-inf_str;
  strain(strmat); potential(); calcstresstensor();
  s2=cell.sgmmat[2][2];
  matcpy(hmat_store, cell.hmat);
  for (int i=1; i<=atom.natom; i++) {
    atom.rx[i] = rx_store[i]; atom.ry[i] = ry_store[i]; atom.rz[i] = rz_store[i]; }

  cell.ecmat[2][2] = (s1 - s2)/(inf_str * 2.0);

  //  C_44
  matcpy(hmat_store, cell.hmat);
  for (int i=1; i<=atom.natom; i++) {
    atom.rx[i] = rx_store[i]; atom.ry[i] = ry_store[i]; atom.rz[i] = rz_store[i]; }
  resetmat(strmat);
  strmat[0][1]=inf_str/2.0; strmat[1][0]=strmat[0][1];
  strain(strmat); potential(); calcstresstensor();
  s1=cell.sgmmat[0][1];
  matcpy(hmat_store, cell.hmat);
  for (int i=1; i<=atom.natom; i++) {
    atom.rx[i] = rx_store[i]; atom.ry[i] = ry_store[i]; atom.rz[i] = rz_store[i]; }
 
  matcpy(hmat_store, cell.hmat);
  for (int i=1; i<=atom.natom; i++) {
    atom.rx[i] = rx_store[i]; atom.ry[i] = ry_store[i]; atom.rz[i] = rz_store[i]; }
  resetmat(strmat);
  strmat[0][1]=-inf_str/2.0; strmat[1][0]=strmat[0][1];
  strain(strmat); potential(); calcstresstensor();
  s2=cell.sgmmat[0][1];
  matcpy(hmat_store, cell.hmat);
  for (int i=1; i<=atom.natom; i++) {
    atom.rx[i] = rx_store[i]; atom.ry[i] = ry_store[i]; atom.rz[i] = rz_store[i]; }

  cell.ecmat[3][3] = (s1 - s2)/(inf_str * 2.0);

  //  C_55
  matcpy(hmat_store, cell.hmat);
  for (int i=1; i<=atom.natom; i++) {
    atom.rx[i] = rx_store[i]; atom.ry[i] = ry_store[i]; atom.rz[i] = rz_store[i]; }
  resetmat(strmat);
  //  strmat[1][2]=inf_str/2.0; strmat[2][1]=strmat[1][2];
  strain(strmat); potential(); calcstresstensor();
  s1=cell.sgmmat[1][2];
  matcpy(hmat_store, cell.hmat);
  for (int i=1; i<=atom.natom; i++) {
    atom.rx[i] = rx_store[i]; atom.ry[i] = ry_store[i]; atom.rz[i] = rz_store[i]; }
 
  matcpy(hmat_store, cell.hmat);
  for (int i=1; i<=atom.natom; i++) {
    atom.rx[i] = rx_store[i]; atom.ry[i] = ry_store[i]; atom.rz[i] = rz_store[i]; }
  resetmat(strmat);
  strmat[1][2]=-inf_str/2.0; strmat[2][1]=strmat[1][2];
  strain(strmat); potential(); calcstresstensor();
  s2=cell.sgmmat[1][2];
  matcpy(hmat_store, cell.hmat);
  for (int i=1; i<=atom.natom; i++) {
    atom.rx[i] = rx_store[i]; atom.ry[i] = ry_store[i]; atom.rz[i] = rz_store[i]; }

  cell.ecmat[4][4] = (s1 - s2)/(inf_str * 2.0);

  ///  C_66
  matcpy(hmat_store, cell.hmat);
  for (int i=1; i<=atom.natom; i++) {
    atom.rx[i] = rx_store[i]; atom.ry[i] = ry_store[i]; atom.rz[i] = rz_store[i]; }
  resetmat(strmat);
  strmat[2][0]=inf_str/2.0; strmat[0][2]=strmat[2][0];
  strain(strmat); potential(); calcstresstensor();
  s1=cell.sgmmat[2][0];
  matcpy(hmat_store, cell.hmat);
  for (int i=1; i<=atom.natom; i++) {
    atom.rx[i] = rx_store[i]; atom.ry[i] = ry_store[i]; atom.rz[i] = rz_store[i]; }
 
  matcpy(hmat_store, cell.hmat);
  for (int i=1; i<=atom.natom; i++) {
    atom.rx[i] = rx_store[i]; atom.ry[i] = ry_store[i]; atom.rz[i] = rz_store[i]; }
  resetmat(strmat);
  strmat[2][0]=-inf_str/2.0; strmat[0][2]=strmat[2][0];
  strain(strmat); potential(); calcstresstensor();
  s2=cell.sgmmat[2][0];
  matcpy(hmat_store, cell.hmat);
  for (int i=1; i<=atom.natom; i++) {
    atom.rx[i] = rx_store[i]; atom.ry[i] = ry_store[i]; atom.rz[i] = rz_store[i]; }

  cell.ecmat[5][5] = (s1 - s2)/(inf_str * 2.0);

  printf("# Elastic constants (numerical evaluation)\n");
  printf("# Using delta_e = %f\n", inf_str);
  printf("# Results not guaranteed\n");
  printf("%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f \n",
	 cell.ecmat[0][0]/GPa,cell.ecmat[0][1]/GPa,cell.ecmat[0][2]/GPa,0.,0.,0.);
  printf("%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f \n",
	 cell.ecmat[0][1]/GPa,cell.ecmat[1][1]/GPa,cell.ecmat[1][2]/GPa,0.,0.,0.);
  printf("%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f \n",
	 cell.ecmat[0][2]/GPa,cell.ecmat[1][2]/GPa,cell.ecmat[2][2]/GPa,0.,0.,0.);
  printf("%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f \n",
	 0.,0.,0.,cell.ecmat[3][3]/GPa,0.,0.);
  printf("%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f \n",
	 0.,0.,0.,0.,cell.ecmat[4][4]/GPa,0.);
  printf("%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f \n",
	 0.,0.,0.,0.,0.,cell.ecmat[5][5]/GPa);
  /*
      write (6,'(6f10.4)') elasc(1,2)/gpa,elasc(2,2)/gpa,elasc(2,3)/gpa,0.,0.,0.
      write (6,'(6f10.4)') elasc(1,3)/gpa,elasc(2,3)/gpa,elasc(3,3)/gpa,0.,0.,0.
      write (6,'(6f10.4)') 0.,0.,0.,elasc(4,4)/gpa,0.,0.
      write (6,'(6f10.4)') 0.,0.,0.,0.,elasc(5,5)/gpa,0.
      write (6,'(6f10.4)') 0.,0.,0.,0.,0.,elasc(6,6)/gpa
  */

  matcpy(hmat_store, cell.hmat);
  for (int i=1; i<=atom.natom; i++) {
    atom.rx[i] = rx_store[i]; atom.ry[i] = ry_store[i]; atom.rz[i] = rz_store[i]; }
  //     strmat_set = strmat_store
  delete[] rx_store; delete[] ry_store; delete[] rz_store;
  writedata_elastic_const();
}
      
      
