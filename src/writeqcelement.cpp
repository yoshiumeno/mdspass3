#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include "myheader.h"

void inverse(double mat[][3], double imat[][3]);

void writeqcelement(const char* fname)
{
  printf ("Write qcelement file = %s\n",fname);
  FILE *fp;
  double xx, yy, zz;
  double hinmat[3][3];
  //fp = fopen("qcelement.out", "w");
  fp = fopen(fname, "w");
  fprintf(fp,"%d\n",atom.nelem);
  for (int i=1; i<=atom.nelem; i++) {
    fprintf(fp,"%d %d %d %d %d %d %d %d\n",atom.elem_v[i][1], atom.elem_v[i][2], atom.elem_v[i][3], atom.elem_v[i][4],
	    atom.elem_v_rep[i][1],atom.elem_v_rep[i][2],atom.elem_v_rep[i][3],atom.elem_v_rep[i][4]);
  }
  fclose(fp);
}
