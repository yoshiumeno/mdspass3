#include <iostream>
#include "myheader.h"

void inverse(double mat[][3], double imat[][3]);
void stretch_celladjust(float c1x, float c1y, float c1z,
			float c2x, float c2y, float c2z,
			float c3x, float c3y, float c3z);

void create_lattice(int lattice_type, float alat, float clat,
		    float &alplat1,
		    float &c1x, float &c1y, float &c1z,
		    float &c2x, float &c2y, float &c2z,
		    float &c3x, float &c3y, float &c3z)
{
  double xx,yy,zz,xxx,yyy,zzz;
  double cosalpx, alpx, alpd, fac;
  if (lattice_type == 0) {
    printf("Rhombohedral lattice is created:\n");
    printf(" Lattice constant a0 = %5.2f A\n",alat);
    printf(" Angle alpha = %5.2f deg\n",alplat1);
    // when alp<20, xx is searched more carefully
    if (alplat1 < 20) { fac = 0.001; } else { fac = 0.01; }
    // alp is restrained to be equal to or below 90
    if (alplat1 > 90.0) { alplat1 = 90.0; }
    double alpha = alplat1/180.0 * M_PI;
    c1x = alat; c1y = 0.0; c1z = 0.0;
    c2x = alat * cos(alpha); c2y = alat * sin(alpha); c2z = 0.0;
    if (alplat1 >= 90 ) {
      c3x = 0.0; c3y = 0.0; c3z = alat;
    } else {
      xx = alat * 0.5; bool noconv = true; int icount = 0;
      while(noconv) {
	c3x = xx * cos(alpha/2.0); c3y = xx * sin(alpha/2.0);
	c3z = sqrt(alat*alat - c3x*c3x - c3y*c3y);
	yy = sqrt(c1x*c1x + c1y*c1y + c1z*c1z);
	zz = sqrt(c3x*c3x + c3y*c3y + c3z*c3z);
	cosalpx = (c1x*c3x + c1y*c3y + c1z*c3z)/(yy*zz);
	alpx = acos(cosalpx) * 180.0 / M_PI;
	//printf("%d %f %e %f\n",icount,xx,cosalpx,alpx);
	alpd = alpx - alplat1;
	if (fabs(alpd) > 0.1) { // 0.1 degree difference is tolerated
	  xx += alpd * alat * fac;
	} else {
	  noconv = false;
	}
	icount++;
	if (icount > 100) {
	  printf("### ERROR in Rhombohedral Lattice Creation (not converged) ###\n");
	  noconv = false;
	}
      }
    }
  }
  
  int itag = ifix_atoms;
  ifix_atoms = 0;
  stretch_celladjust(cell1x,cell1y,cell1z,cell2x,cell2y,cell2z,cell3x,cell3y,cell3z);
  ifix_atoms = itag;
}

