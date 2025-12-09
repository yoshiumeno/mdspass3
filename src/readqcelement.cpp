#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "myheader.h"

void get_first_arg(std::string &line, std::string &arg1);
int remove_head_spaces(std::string &line);
int sprit_first_arg(std::string &line, std::string &arg1);
int count_arg_number(std::string line);

void mk_helmat(int iel, double rx[], double ry[], double rz[], double mat[][3]);
void inverse(double mat[][3], double imat[][3]);
double q_x(int i, int iel, double rx[], double ry[], double rz[], double hinelmat[][3]);
double q_y(int i, int iel, double rx[], double ry[], double rz[], double hinelmat[][3]);
double q_z(int i, int iel, double rx[], double ry[], double rz[], double hinelmat[][3]);

void readqcelement(const char* fname)
{
  double helmat[3][3], hinelmat[3][3];
  double qx, qy, qz;

  printf ("Read qcelement file = %s\n", fname);
  strcpy(current_qcelement_name,fname);
  std::ifstream fin( fname );
  double xx, yy, zz;
  int nelem, ielem;
  std::string line;
  std::string arg1, arg2, arg3, arg4, arg5, arg6, arg7;
  char larg1[40], larg2[40], larg3[40];
  bool abs;

  if (atom.elem_v) { for (int i=1; i<=atom.nelem; i++) { if (atom.elem_v[i]) free(atom.elem_v[i]); } free(atom.elem_v); }
  if (atom.elem_v_rep) {
    for (int i=1; i<=atom.nelem; i++) { if (atom.elem_v_rep[i]) free(atom.elem_v_rep[i]); } free(atom.elem_v_rep); }
  //  fin >> nelem;
  //  for (int i=1; i<=nelem; i++) {
  fin >> atom.nelem;
  atom.elem_v     = (int **)malloc(sizeof(int*)*(atom.nelem+1));
  atom.elem_v_rep = (int **)malloc(sizeof(int*)*(atom.nelem+1));
  for (int i=1; i<=atom.nelem; i++) {
    atom.elem_v[i]     = (int *)malloc(sizeof(int)*5);
    atom.elem_v_rep[i] = (int *)malloc(sizeof(int)*5);
  }
  for (int i=1; i<=atom.nelem; i++) {
    //    fin >> ielem;
    fin >> atom.elem_v[i][1]; fin >> atom.elem_v[i][2]; fin >> atom.elem_v[i][3]; fin >> atom.elem_v[i][4];
    fin >> atom.elem_v_rep[i][1]; fin >> atom.elem_v_rep[i][2]; fin >> atom.elem_v_rep[i][3]; fin >> atom.elem_v_rep[i][4];
			      
    //    std::cout<<atom.elem_v[i][1]<<" "<<atom.elem_v[i][2]<<" "<<atom.elem_v[i][3]<<" "<<atom.elem_v[i][4]<<" ";
    //    std::cout<<atom.elem_v_rep[i][1]<<" "<<atom.elem_v_rep[i][2]<<" "<<atom.elem_v_rep[i][3]<<" "<<atom.elem_v_rep[i][4]<<std::endl;
  }

  // set atom.repatom[i] (inside an element -> 0, outside or on a vertex -> 1)
  atom.nrepatom = 0;
  for (int i=1; i<=atom.natom; i++) {
    atom.repatom[i] = 1; atom.elem_id[i] = 0;
    for (int iel=1; iel<=atom.nelem; iel++) {
      mk_helmat(iel, atom.rx, atom.ry, atom.rz, helmat);
      inverse(helmat, hinelmat);
      //      qx = q_x(i, atom.elem_v[iel][1], atom.rx, atom.ry, atom.rz, hinelmat);
      //      qy = q_y(i, atom.elem_v[iel][1], atom.rx, atom.ry, atom.rz, hinelmat);
      //      qz = q_z(i, atom.elem_v[iel][1], atom.rx, atom.ry, atom.rz, hinelmat);
      qx = q_x(i, iel, atom.rx, atom.ry, atom.rz, hinelmat);
      qy = q_y(i, iel, atom.rx, atom.ry, atom.rz, hinelmat);
      qz = q_z(i, iel, atom.rx, atom.ry, atom.rz, hinelmat);
      //      std::cout<<i<<" "<<iel<<" "<<atom.elem_v[iel][1]<<" "<<qx<<" "<<qy<<std::endl;
      //      if ((qx>=0.0)&&(qy>=0.0)&&(qz>=0.0)&&(qx+qy+qz<=1.0)) { atom.repatom[i] = 0; atom.elem_id[i] = iel; }
      float del=5.0e-2;
      if ((qx>=-del)&&(qy>=-del)&&(qz>=-del)&&(qx+qy+qz<=1.0+del)) { atom.repatom[i] = 0; atom.elem_id[i] = iel; }
      if ((fabs(qx)<del)&&(fabs(qy)<del)&&(fabs(qz)<del)) { atom.repatom[i] = 1; atom.elem_id[i] = iel; }
      if ((fabs(qx-1.0)<del)&&(fabs(qy)<del)&&(fabs(qz)<del)) { atom.repatom[i] = 1; atom.elem_id[i] = iel; }
      if ((fabs(qx)<del)&&(fabs(qy-1.0)<del)&&(fabs(qz)<del)) { atom.repatom[i] = 1; atom.elem_id[i] = iel; }
      if ((fabs(qx)<del)&&(fabs(qy)<del)&&(fabs(qz-1.0)<del)) { atom.repatom[i] = 1; atom.elem_id[i] = iel; }
      

    }
    //      std::cout<<i<<" "<<atom.repatom[i]<<" "<<atom.elem_id[i]<<std::endl;
    if (atom.repatom[i] == 1) { atom.nrepatom++; }
  }
  std::cout<<"Atom:"<<atom.natom<<" Repatom:"<<atom.nrepatom<<" Element:"<<atom.nelem<<std::endl;

}

void readqcelement()
{
  double helmat[3][3], hinelmat[3][3];
  double qx, qy, qz;

  // printf ("Read qcelement file = %s\n", fname);
  double xx, yy, zz;
  int nelem, ielem;
  std::string line;
  std::string arg1, arg2, arg3, arg4, arg5, arg6, arg7;
  char larg1[40], larg2[40], larg3[40];
  bool abs;

  // set atom.repatom[i] (inside an element -> 0, outside or on a vertex -> 1)
  atom.nrepatom = 0;
  for (int i=1; i<=atom.natom; i++) {
    atom.repatom[i] = 1; atom.elem_id[i] = 0;
    for (int iel=1; iel<=atom.nelem; iel++) {
      mk_helmat(iel, atom.rx, atom.ry, atom.rz, helmat);
      inverse(helmat, hinelmat);
      qx = q_x(i, iel, atom.rx, atom.ry, atom.rz, hinelmat);
      qy = q_y(i, iel, atom.rx, atom.ry, atom.rz, hinelmat);
      qz = q_z(i, iel, atom.rx, atom.ry, atom.rz, hinelmat);
      float del=5.0e-2;
      if ((qx>=-del)&&(qy>=-del)&&(qz>=-del)&&(qx+qy+qz<=1.0+del)) { atom.repatom[i] = 0; atom.elem_id[i] = iel; }
      if ((fabs(qx)<del)&&(fabs(qy)<del)&&(fabs(qz)<del)) { atom.repatom[i] = 1; atom.elem_id[i] = iel; }
      if ((fabs(qx-1.0)<del)&&(fabs(qy)<del)&&(fabs(qz)<del)) { atom.repatom[i] = 1; atom.elem_id[i] = iel; }
      if ((fabs(qx)<del)&&(fabs(qy-1.0)<del)&&(fabs(qz)<del)) { atom.repatom[i] = 1; atom.elem_id[i] = iel; }
      if ((fabs(qx)<del)&&(fabs(qy)<del)&&(fabs(qz-1.0)<del)) { atom.repatom[i] = 1; atom.elem_id[i] = iel; }
      

    }
    if (atom.repatom[i] == 1) { atom.nrepatom++; }
  }
  std::cout<<"Atom:"<<atom.natom<<" Repatom:"<<atom.nrepatom<<" Element:"<<atom.nelem<<std::endl;

}

/*
void get_first_arg(std::string &line, std::string &arg1)
{
  int istt = line.find_first_not_of(" ");
  if (istt == -1) {
    arg1 = "";
  } else {
    line = line.substr(istt);
    int iend = line.find_first_of(" ");
    if (iend == -1) { iend = line.length(); }
    //  std::cout<<istt<< " " <<iend<<std::endl;
    arg1 = line.substr(0,iend);
    //  std::cout<<line <<" " <<arg1 <<std::endl;
    line = line.substr(iend);
  }
}

int remove_head_spaces(std::string &line)
{
  int istt = line.find_first_not_of(" ");
  if (istt == -1) {
    return -1;
  } else {
    line=line.substr(istt);
    return 0;
  }
}

int sprit_first_arg(std::string &line, std::string &arg1)
{
    int iend = line.find_first_of(" ");
    //    std::cout<<iend<<std::endl;
    if (iend == -1) { // no space, only one arg
      arg1 = line;
      line = "";
      return -1;
    } else  {
     arg1 = line.substr(0,iend);
     line = line.substr(iend);
     if (remove_head_spaces(line) == -1) {
       line = "";
       return -1;
     }
     return 0;
    }
}
  
int count_arg_number(std::string line)
{
  std::string cline, arg;
  cline = line;
  int count = 0;
  if (remove_head_spaces(cline) == -1) {
    return 0;
  } 
  for (int i=1; i<=100; i++) {
    if (sprit_first_arg(cline, arg) == -1) {
      count++;
      return count;
    } else {
      count++;
    }
  }
  std::cout << "Error: a line has more than 100 args" << std::endl;
  return 999;
}
*/
