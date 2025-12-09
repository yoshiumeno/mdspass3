#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "myheader.h"

void matcpy(double a[3][3], double b[3][3]);
void resetmat(double a[3][3]);
void resetmat6(double a[6][6]);
void potential();
void e_force_dipole(int mode);

void force_check(double dis)
{
  double threshold=30.0;             // Error exceeding this value is shown
  double rxk, ryk, rzk, fx_cal, fy_cal, fz_cal, enpot0, fx_ene, fy_ene, fz_ene;
  double errx, erry, errz;
  dis *= ang;

  printf("Force check ...\n");
  double xx0,xx1;

  /*
  rxk=atom.rx[10];
  for (int i=-10; i<=10; i++) {
    atom.rx[10]=rxk+dis*(double)i;
    e_force_dipole(0);
    atom.rx[10]=rxk;
    e_force_dipole(1);
    printf("%d %20.10e\n",i,atom.epotsum/eV);
  }
  exit(0);
  */

  for (int i=1; i<=atom.natom; i++) {
    rxk=atom.rx[i]; ryk=atom.ry[i]; rzk=atom.rz[i];
    potential();
    //e_force_dipole(0);
    fx_cal=atom.fx[i]; fy_cal=atom.fy[i]; fz_cal=atom.fz[i];

    atom.rx[i]=rxk-dis;
    potential();
    //e_force_dipole(1);
    enpot0=atom.epotsum;
    atom.rx[i]=rxk+dis;
    potential();
    //e_force_dipole(1);
    //if (i==1) { xx1=dipole.val1; }
    //if (i==1) { printf("%e %e\n",-(xx1-xx0)/(dis*2/ang),dipole.val1); }
    //if (i==1) { printf("%e %e\n",xx1,xx0); }
    fx_ene=-(atom.epotsum-enpot0)/(dis*2.0);
    atom.rx[i]=rxk;
    errx=fabs(fx_cal/fx_ene*100.0-100);
    
    atom.ry[i]=ryk-dis;
    potential();
    //e_force_dipole(1);
    enpot0=atom.epotsum;
    atom.ry[i]=ryk+dis;
    potential();
    //e_force_dipole(1);
    fy_ene=-(atom.epotsum-enpot0)/(dis*2.0);
    atom.ry[i]=ryk;
    erry=fabs(fy_cal/fy_ene*100.0-100);

    atom.rz[i]=rzk-dis;
    potential();
    //e_force_dipole(1);
    enpot0=atom.epotsum;
    atom.rz[i]=rzk+dis;
    potential();
    //e_force_dipole(1);
    fz_ene=-(atom.epotsum-enpot0)/(dis*2.0);
    atom.rz[i]=rzk;
    errz=fabs(fz_cal/fz_ene*100.0-100);

    fx_cal /= eV/ang; fy_cal /= eV/ang; fz_cal /= eV/ang;
    fx_ene /= eV/ang; fy_ene /= eV/ang; fz_ene /= eV/ang;
    if (fabs(fx_ene)<1.e-6) { errx=-1.0; }
    if (fabs(fy_ene)<1.e-6) { erry=-1.0; }
    if (fabs(fz_ene)<1.e-6) { erry=-1.0; }

    if (errx >= 0) {
      printf("Atom = %4d (x) Numer(eV/A): %+10.4e Analy: %+10.4e Err %10.4f %10.4f (%)\n",
	     i,fx_ene,fx_cal,fabs(fx_ene-fx_cal),errx);
    } else {
      printf("Atom = %4d (x) Numer(eV/A): %+10.4e Analy: %+10.4e Err %10.4f        --- (%)\n",
	     i,fx_ene,fx_cal,fabs(fx_ene-fx_cal));
    }
    if (erry >= 0) {
      printf("            (y) Numer(eV/A): %+10.4e Analy: %+10.4e Err %10.4f %10.4f (%)\n",
	     fy_ene,fy_cal,fabs(fy_ene-fy_cal),erry);
    } else {
      printf("            (y) Numer(eV/A): %+10.4e Analy: %+10.4e Err %10.4f        --- (%)\n",
	     fy_ene,fy_cal,fabs(fy_ene-fy_cal),erry);
    }
    if (errz >= 0) {
      printf("            (z) Numer(eV/A): %+10.4e Analy: %+10.4e Err %10.4f %10.4f (%)\n",
	     fz_ene,fz_cal,fabs(fz_ene-fz_cal),errz);
    } else {
      printf("            (z) Numer(eV/A): %+10.4e Analy: %+10.4e Err %10.4f        --- (%)\n",
	     fz_ene,fz_cal,fabs(fz_ene-fz_cal),errz);
    }
  }
}

