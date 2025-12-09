#include <iostream>
#include "myheader.h"

void add_atom(const char *aasp, float arx, float ary, float arz);

extern char config_atom[];
extern int iatom_from, iatom_to;

void add_polyhed(int type, float size)
{
  float px[20], py[20], pz[20];
  int nver = 0; // number of vertices
  double radc = 0; // circumscribed sphere radius
  double radi = 0; // inscribed sphere radius
  double lenedge = 0; // edge length
  if (type == 0) {
    printf("Tetrahedron is added.\n");
    px[ 0] =  -1.0; py[ 0] = -1.0; pz[ 0] =  -1.0;
    px[ 1] =  -1.0; py[ 1] =  1.0; pz[ 1] =   1.0;
    px[ 2] =   1.0; py[ 2] = -1.0; pz[ 2] =   1.0;
    px[ 3] =   1.0; py[ 3] =  1.0; pz[ 3] =  -1.0;
    nver = 4;
    radc = sqrt(3.0);
    radi = 1.0/sqrt(3.0);
    lenedge = 2.0*sqrt(2.0);
  } else if (type == 1) {
    printf("Hexahedron (cube) is added.\n");
    px[ 0] =  -1.0; py[ 0] = -1.0; pz[ 0] =  -1.0;
    px[ 1] =  -1.0; py[ 1] = -1.0; pz[ 1] =   1.0;
    px[ 2] =  -1.0; py[ 2] =  1.0; pz[ 2] =  -1.0;
    px[ 3] =  -1.0; py[ 3] =  1.0; pz[ 3] =   1.0;
    px[ 4] =   1.0; py[ 4] = -1.0; pz[ 4] =  -1.0;
    px[ 5] =   1.0; py[ 5] = -1.0; pz[ 5] =   1.0;
    px[ 6] =   1.0; py[ 6] =  1.0; pz[ 6] =  -1.0;
    px[ 7] =   1.0; py[ 7] =  1.0; pz[ 7] =   1.0;
    nver = 8;
    radc = sqrt(3.0);
    radi = 1.0;
    lenedge = 2.0;
  } else if (type == 2) {
    printf("Octahedron is added.\n");
    px[ 0] =  -1.0; py[ 0] =  0.0; pz[ 0] =   0.0;
    px[ 1] =   0.0; py[ 1] = -1.0; pz[ 1] =   0.0;
    px[ 2] =   0.0; py[ 2] =  0.0; pz[ 2] =  -1.0;
    px[ 3] =   1.0; py[ 3] =  0.0; pz[ 3] =   0.0;
    px[ 4] =   0.0; py[ 4] =  1.0; pz[ 4] =   0.0;
    px[ 5] =   0.0; py[ 5] =  0.0; pz[ 5] =   1.0;
    nver = 6;
    radc = 1.0;
    radi = 1.0/sqrt(3.0);
    lenedge = sqrt(2.0);
  } else if (type == 3) {
    printf("Dodecahedron is added.\n");
    float phi = (1.0+sqrt(5.0))/2.0;
    float phi2 = phi*phi;
    px[ 0] =    0.0; py[ 0] =  -1.0; pz[ 0] =  -phi2;
    px[ 1] =    0.0; py[ 1] =   1.0; pz[ 1] =  -phi2;
    px[ 2] =    0.0; py[ 2] =  -1.0; pz[ 2] =   phi2;
    px[ 3] =    0.0; py[ 3] =   1.0; pz[ 3] =   phi2;

    px[ 4] =   -1.0; py[ 4] = -phi2; pz[ 4] =   0.0;
    px[ 5] =   -1.0; py[ 5] =  phi2; pz[ 5] =   0.0;
    px[ 6] =    1.0; py[ 6] = -phi2; pz[ 6] =   0.0;
    px[ 7] =    1.0; py[ 7] =  phi2; pz[ 7] =   0.0;

    px[ 8] =  -phi2; py[ 8] =   0.0; pz[ 8] =  -1.0;
    px[ 9] =  -phi2; py[ 9] =   0.0; pz[ 9] =   1.0;
    px[10] =   phi2; py[10] =   0.0; pz[10] =  -1.0;
    px[11] =   phi2; py[11] =   0.0; pz[11] =   1.0;

    px[12] =   -phi; py[12] =  -phi; pz[12] =  -phi;
    px[13] =   -phi; py[13] =  -phi; pz[13] =   phi;
    px[14] =   -phi; py[14] =   phi; pz[14] =  -phi;
    px[15] =   -phi; py[15] =   phi; pz[15] =   phi;
    px[16] =    phi; py[16] =  -phi; pz[16] =  -phi;
    px[17] =    phi; py[17] =  -phi; pz[17] =   phi;
    px[18] =    phi; py[18] =   phi; pz[18] =  -phi;
    px[19] =    phi; py[19] =   phi; pz[19] =   phi;

    nver = 20;
    radc = phi*sqrt(3.0);
    radi = sqrt((25.0+11.0*sqrt(5.0))/10.0);
    lenedge = 2.0;
  } else if (type == 4) {
    printf("Icosahedron is added.\n");
    float phi = (1.0+sqrt(5.0))/2.0;
    px[ 0] =   0.0; py[ 0] = -1.0; pz[ 0] =  -phi;
    px[ 1] =   0.0; py[ 1] =  1.0; pz[ 1] =  -phi;
    px[ 2] =   0.0; py[ 2] = -1.0; pz[ 2] =   phi;
    px[ 3] =   0.0; py[ 3] =  1.0; pz[ 3] =   phi;

    px[ 4] =  -phi; py[ 4] =  0.0; pz[ 4] =  -1.0;
    px[ 5] =  -phi; py[ 5] =  0.0; pz[ 5] =   1.0;
    px[ 6] =   phi; py[ 6] =  0.0; pz[ 6] =  -1.0;
    px[ 7] =   phi; py[ 7] =  0.0; pz[ 7] =   1.0;

    px[ 8] =  -1.0; py[ 8] = -phi; pz[ 8] =   0.0;
    px[ 9] =   1.0; py[ 9] = -phi; pz[ 9] =   0.0;
    px[10] =  -1.0; py[10] =  phi; pz[10] =   0.0;
    px[11] =   1.0; py[11] =  phi; pz[11] =   0.0;
    nver = 12;
    radc = sqrt(sqrt(phi));
    radi = phi*phi/sqrt(3.0);
    lenedge = 2;
  }
  
  // Apply size factor
  radc *= size; radi *= size; lenedge *= size;
  for (int i=0; i<nver; i++) {
    px[i] *= size; py[i] *= size; pz[i] *= size;
  }
  // Create new atoms on vertices
  iatom_from = atom.natom+1;
  for (int i=0; i<nver; i++) {
    add_atom(config_atom, px[i], py[i], pz[i]);
    printf(" Atom %5d : (%5.2f, %5.2f, %5.2f) \n",atom.natom,px[i], py[i], pz[i]);
  }
  iatom_to = atom.natom;
  printf(" %4d vertices\n",nver);
  printf(" circumscribed sphere rad = %5.2f A\n",radc);
  printf(" inscribed sphere rad     = %5.2f A\n",radi);
}
