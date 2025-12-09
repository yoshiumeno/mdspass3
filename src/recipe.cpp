#include <string.h>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "myheader.h"

extern int istep0, inst_mode;
bool check_arg(std::string arg);
void get_first_arg(std::string &line, std::string &arg1);
void writeconfig(const char* fname); void writeconfig_abs(const char* fname);
void readsetdat(const char* fname);
void capture(); void capture(const char* fname);
void stretch(double x, double y, double z);
void stretch_celladjust(float x, float y, float z);
void instability_atomcell();
void instability_atom();
void instability_atom_noremove();
void calc_elastic_const();
void remove_carriage_return(std::string& line);

void recipe()
{
  static int counter;
  char fn[60] = "RECIPE";
  std::string line, arg1, arg2;
  std::ifstream fin(fn); fin.close(); fin.clear();

  fin.open(fn);
  if (!fin) { printf("Recipe file %s does not exist\n",fn); return; }
  printf("recipe file = %s  line = %d\n",fn,counter);
  if (counter>0) { for (int i=0; i<counter; i++) { getline(fin,line); remove_carriage_return(line); } }

 IN:
  getline(fin,line); remove_carriage_return(line);
  get_first_arg(line,arg1); get_first_arg(line, arg2);
  bool repeat = true;
  if (strcmp(arg1.c_str(), "")==0) { //If line is empty or end of the file, reset and exit
    counter = -1; repeat = false;
  } else if (strcmp(arg1.c_str(), "MD")==0) { //MD is turned on
    mdmotion=1; istep0 = istep; repeat = false;
  } else if (strcmp(arg1.c_str(), "WRITE")==0) { //Write config (fra)
    if (check_arg(arg2)) { writeconfig(arg2.c_str()); }
    else { writeconfig("CONFIG.OUT"); }
  } else if (strcmp(arg1.c_str(), "WRITE_ABS")==0) { //Write config (abs)
    if (check_arg(arg2)) { writeconfig_abs(arg2.c_str()); }
    else { writeconfig_abs("CONFIG.OUT.ABS"); }
  } else if (strcmp(arg1.c_str(), "READ_SETDAT")==0) { //Read new SETDAT
    if (check_arg(arg2)) { readsetdat(arg2.c_str()); }
  } else if (strcmp(arg1.c_str(), "CAPTURE")==0) { //Capture
    if (check_arg(arg2)) { capture(arg2.c_str()); }
    else { capture(); }
  } else if (strcmp(arg1.c_str(), "SET_LX")==0) { //Deformation
    if (check_arg(arg2)) { cellx = atof(arg2.c_str()); }
    stretch_celladjust(cellx,celly,cellz);
  } else if (strcmp(arg1.c_str(), "SET_LY")==0) { //Deformation
    if (check_arg(arg2)) { celly = atof(arg2.c_str()); }
    stretch_celladjust(cellx,celly,cellz);
  } else if (strcmp(arg1.c_str(), "SET_LZ")==0) { //Deformation
    if (check_arg(arg2)) { cellz = atof(arg2.c_str()); }
    stretch_celladjust(cellx,celly,cellz);
  } else if (strcmp(arg1.c_str(), "SET_EX")==0) { //Deformation
    if (check_arg(arg2)) { stretch(1.0+atof(arg2.c_str()),1.0,1.0); }
  } else if (strcmp(arg1.c_str(), "SET_EY")==0) { //Deformation
    if (check_arg(arg2)) { stretch(1.0,1.0+atof(arg2.c_str()),1.0); }
  } else if (strcmp(arg1.c_str(), "SET_EZ")==0) { //Deformation
    if (check_arg(arg2)) { stretch(1.0,1.0,1.0+atof(arg2.c_str())); }
  } else if (strcmp(arg1.c_str(), "INSTABILITY")==0) { //Instability
    if (inst_mode == 0) { instability_atomcell(); }
    else if (inst_mode == 1) { instability_atom(); }
    else if (inst_mode == 2) { instability_atom_noremove(); }
  } else if (strcmp(arg1.c_str(), "CALC_EC")==0) { //Calc elastic const
    calc_elastic_const();
  } else if (strcmp(arg1.c_str(), "QUIT")==0) { //Quit
    exit(0);
  }
  counter++;
  if (repeat) goto IN;

  fin.close(); fin.clear();
}

bool check_arg(std::string arg)
{
  if (strcmp(arg.c_str(), "")==0) {
    printf(" Error in RECIPE. Empty argument.\n");
    return false;
  } else {
    return true;
  }
}
