#ifndef VIEW_CPP
#define VIEW_CPP

#include <fstream>
#include <iostream>
#include "myheader.h"

extern float rotmat[16], obj_pos[3], scl;

void view_save()
{
  FILE *fp = fopen("VIEW","w");
  fprintf(fp, "%f %f %f %f\n",rotmat[ 0],rotmat[ 1],rotmat[ 2],rotmat[ 3]);
  fprintf(fp, "%f %f %f %f\n",rotmat[ 4],rotmat[ 5],rotmat[ 6],rotmat[ 7]);
  fprintf(fp, "%f %f %f %f\n",rotmat[ 8],rotmat[ 9],rotmat[10],rotmat[11]);
  fprintf(fp, "%f %f %f %f\n",rotmat[12],rotmat[13],rotmat[14],rotmat[15]);
  fprintf(fp, "%f %f %f\n", obj_pos[0], obj_pos[1], obj_pos[2]);
  fprintf(fp, "%f\n",scl);
  fclose(fp);
  printf("View data is saved to VIEW\n");
}

void view_read()
{
  FILE *fp;
  if (fp = fopen("VIEW","r")) {
    fscanf(fp, "%f %f %f %f\n",&rotmat[ 0],&rotmat[ 1],&rotmat[ 2],&rotmat[ 3]);
    fscanf(fp, "%f %f %f %f\n",&rotmat[ 4],&rotmat[ 5],&rotmat[ 6],&rotmat[ 7]);
    fscanf(fp, "%f %f %f %f\n",&rotmat[ 8],&rotmat[ 9],&rotmat[10],&rotmat[11]);
    fscanf(fp, "%f %f %f %f\n",&rotmat[12],&rotmat[13],&rotmat[14],&rotmat[15]);
    fscanf(fp, "%f %f %f\n", &obj_pos[0], &obj_pos[1], &obj_pos[2]);
    fscanf(fp, "%f\n",&scl);
    fclose(fp);
    printf("View data is read from VIEW\n");
  }
}

#endif // VIEW_CPP
