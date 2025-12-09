#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#if defined(_WIN32) || defined(_WIN64)
#include <Windows.h>
#endif
#include <GL/gl.h>
#include "myheader.h"

#if defined __linux__ || defined __APPLE__
#include <unistd.h>
#else
#include <direct.h>
#endif
#if !defined __linux__ && !defined __APPLE__
#define snprintf sprintf_s
#endif

std::string replace_string(std::string string1, std::string string2, std::string string3);
void get_first_arg(std::string &line, std::string &arg1);
void get_first_arg(char* line, char* arg1);
int remove_head_spaces(std::string &line);
int sprit_first_arg(std::string &line, std::string &arg1);
int count_arg_number(std::string line);
void set_atom_weight();
void potential();
void out_of_cell();
void bookkeep();
void pot_initialize_all();
void set_atom_color();
void change_atom_color();
int atom_number(char* at);
void deallocate_arrays();
void allocate_arrays();
void readconfig(const char* fname, int mode);
void set_potfile(const char* fname);
void velocity_random(double target);
void stress_set();
void transpose(double a[3][3]);
void remove_carriage_return(std::string& line);

extern GLuint objects;
//extern GLfloat *color[NOBJECTS];
extern GLfloat **color, **color0;
extern int *iatom, *repidx;
extern GLfloat red[], yellow[];
extern float ex;
extern int incell;
extern char cwdname[80];
extern void mark_display_lists_dirty(); // Dec 2025

void readconfig(const char* fname)
{
  readconfig(fname, 0);
}

void readconfig(const char* fname, int mode)
{ // mode = 0: full, 1: read only (no array reallocation)
  //  std::ofstream foutene( "energy.d", std::ios::out | std::ios::app );
  
  //  foutene << istep << " " << epotsum << " " << ekinsum << std::endl;
  printf ("Read config file = %s\n", fname);
  strcpy(current_config_name,fname);
  std::ifstream finconfig( fname );
  //std::ifstream finconfig; finconfig.open(fname);
  if (!finconfig) { printf("### File [%s] does not exist! ###\n",fname); return; }
  double xx, yy, zz;
  int ix;
  std::string line;
  std::string arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8;
  char larg1[40], larg2[40], larg3[40], larg4[40];
  char larg5[2], larg6[2], larg7[2], larg8[5];
  bool abs = false;
  char species[3];
  double *rxs, *rys, *rzs, hmats[3][3];
  int natoms, narg;

  natoms = 0;
  if ((mode==0)&&(imerge==1)) {
    printf("Read config: Merge mode (Order: Read atoms -> Existing atoms)\n");
    natoms = atom.natom;
    rxs = new double[atom.natom+1]; rys = new double[atom.natom+1]; rzs = new double[atom.natom+1];
    for (int i=1; i<=atom.natom; i++) {
      rxs[i] = atom.rx[i]; rys[i] = atom.ry[i]; rzs[i] = atom.rz[i];
    }
    for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
	hmats[i][j] = cell.hmat[i][j];
      } }
  }
  atom.nwall = 0;

  if (mode==0) {
    deallocate_arrays();

    //printf("glDeletelists");
    //  if (glIsList(objects)) {glDeleteLists(objects, atom.natom*3+1); printf("glDeleteLists\n");}
    //if (glIsList(objects+1)) {glDeleteLists(objects, atom.natom*3+1); printf("glDeleteLists\n"); }

    // Dec. 2025
    //if (glIsList(objects)) {
    //  glDeleteLists(objects, atom.natom*3); printf("glDeleteLists %d %d\n",objects, atom.natom*3);
    //}
    mark_display_lists_dirty();

    //glDeleteLists(objects, atom.natom*3); printf("glDeleteLists %d %d\n",objects, atom.natom*3);
  }

  // 1st line
  getline (finconfig, line); remove_carriage_return(line);
  //  std::cout << count_arg_number(line) << std::endl;
  //  std::cout << "LINE1 = " << "\"" << line << "\"" << std::endl;
  get_first_arg(line, arg1); get_first_arg(line, arg2);
  //  std::cout << "ARG1 = " << "\"" << arg1 << "\"" << std::endl;
  //  std::cout << "ARG2 = " << "\"" << arg2 << "\"" << std::endl;
  //  strcpy(species, arg1.c_str()); printf("Species (1st line) = %s\n",species);
  strncpy(species, arg1.c_str(),sizeof(species)); species[2]='\0';
  if (mode==0) { printf("Species (1st line) = %s\n",species); }

  // 2nd line
  getline (finconfig, line); remove_carriage_return(line); narg = count_arg_number(line);
  //  std::cout << "LINE2 = " << "\"" << line << "\"" << std::endl;
  get_first_arg(line, arg1); get_first_arg(line, arg2);
  //  std::cout << "ARG1 = " << "\"" << arg1 << "\"" << std::endl;
  //  std::cout << "ARG2 = " << "\"" << arg2 << "\"" << std::endl;
  ///strncpy(atom.potential_func, arg1.c_str(),sizeof(atom.potential_func)); //<--???
  strcpy(atom.potential_func, arg1.c_str());
  if (mode==0) { printf("Potential (2nd line) = '%s' ",atom.potential_func); }
  ///strncpy(atom.potential_arg,  arg2.c_str(),sizeof(atom.potential_arg)); //<--???
  strcpy(atom.potential_arg,  arg2.c_str());
  if (mode==0) { printf(": Potential arg (2nd line) = '%s'\n",atom.potential_arg); }
  book.algo = 1;
  //if (strcmp(atom.potential_func,"Tersoff")==0) { book.algo = 2; }
  if (strcmp(atom.potential_func,"Brenner")==0) { book.algo = 3; }
  if (strcmp(atom.potential_func,"AIREBO" )==0) { book.algo = 3; }
  if ((strcmp(atom.potential_func, "Dipole") == 0)||(strcmp(atom.potential_func, "ADP") == 0))
    { if (narg > 1) { set_potfile(atom.potential_arg); } }
  if (strcmp(atom.potential_func, "NiYSZ") == 0) { // Combined mode (Pair+EAM+Dipole)
    if (narg >= 3) { get_first_arg(line, arg3);
      strcpy(atom.potential_func,arg2.c_str()); set_potfile(arg3.c_str()); }
    if (narg >= 5) { get_first_arg(line, arg4); get_first_arg(line, arg5);
      strcpy(atom.potential_func,arg4.c_str()); set_potfile(arg5.c_str()); }
    if (narg >= 7) { get_first_arg(line, arg6); get_first_arg(line, arg7);
      strcpy(atom.potential_func,arg6.c_str()); set_potfile(arg7.c_str()); }
    strcpy(atom.potential_func, "NiYSZ");
  }

  //Correct obsolete tags for potential
  if (strcmp(atom.potential_func,"GEAM")==0) {
    strcpy(atom.potential_func,"EAM"); strcpy(atom.potential_arg,"GEAM"); }
  if (strcmp(atom.potential_func,"SiC_Vashishta")==0) {
    strcpy(atom.potential_func,"SW"); strcpy(atom.potential_arg,"SiC_Vashishta"); }

  ipottype=-99;
  for (int i=0; i<MAXPOTTYPE; i++) {
    //printf("%s %s\n",atom.potential_func, potstring_list[i]);
    std::string line = std::string(atom.potstring_list[i]);
    std::string arg1, arg2; narg = count_arg_number(line);
    if (narg==1) { get_first_arg(line,arg1);
    } else if (narg>=2) { get_first_arg(line,arg1); get_first_arg(line,arg2); }
    //if (strcmp(atom.potential_func, potstring_list[i]) == 0) { ipottype=i; }
    if (strcmp(atom.potential_func, arg1.c_str()) == 0) {
      if (narg==1) { ipottype=i;
      } else {
	if (strcmp(atom.potential_arg, arg2.c_str()) == 0) { ipottype=i; } }
    }
  }

  if (ipottype<0) {
    printf("################################################\n");
    printf("### WARNING: No match in available potential ###\n");
    printf("In CONFIG: %s\n",atom.potential_func);
    printf("Available potentials are:\n");
    for (int i=0; i<MAXPOTTYPE; i++) {
      printf("%d %s\n",i,atom.potstring_list[i]);
    }
    printf("Potential is set to Morse\n");
    printf("################################################\n");
    strcpy(atom.potential_func,"Morse");ipottype=0;
  }


  // number of atoms, lattice constant and cell matrix
  finconfig >> atom.natom; //atom.nrepatom=atom.natom;
  std::cout << "Number of atoms = "<<atom.natom << std::endl;
  atom.natom = atom.natom + natoms;

  if (mode==0) {
    allocate_arrays();

    // Dec. 2025
    // objects = glGenLists(atom.natom*3); printf("objects = %d\n",objects);
    mark_display_lists_dirty();

    //for (int i=1; i<=atom.natom*3; i++) { memcpy(color[i],yellow,sizeof(GLfloat)*4); }
  }

  for (int i=1; i<=atom.natom; i++) { strcpy (atom.asp[i], species); }

  //  std::cout << atom.asp[1] << std::endl;

  finconfig >> cell.alat;
  finconfig >> cell.hmat[0][0];   finconfig >> cell.hmat[1][0];   finconfig >> cell.hmat[2][0]; 
  finconfig >> cell.hmat[0][1];   finconfig >> cell.hmat[1][1];   finconfig >> cell.hmat[2][1]; 
  finconfig >> cell.hmat[0][2];   finconfig >> cell.hmat[1][2];   finconfig >> cell.hmat[2][2]; 
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      cell.hmat[i][j] = cell.hmat[i][j] * cell.alat * 1.0e-10;
      if ((mode==0)&&(imerge==1)) {
	float xab = cell.hmat[i][j]-hmats[i][j]; if (xab<0) { xab=-xab; }
	if (xab>1.0e-10) {
	  printf("##WARNING## hmat[%d][%d] is much different: ",i,j);
	  printf("%8.2f => %8.2f ang\n",hmats[i][j]*1.0e10, cell.hmat[i][j]*1.0e10);
	} } 
    }
  }
  //cellx=cell.hmat[0][0]/ang;celly=cell.hmat[1][1]/ang;cellz=cell.hmat[2][2]/ang;
  cellx=sqrt(cell.hmat[0][0]*cell.hmat[0][0]+cell.hmat[1][0]*cell.hmat[1][0]+cell.hmat[2][0]*cell.hmat[2][0])/ang;
  celly=sqrt(cell.hmat[0][1]*cell.hmat[0][1]+cell.hmat[1][1]*cell.hmat[1][1]+cell.hmat[2][1]*cell.hmat[2][1])/ang;
  cellz=sqrt(cell.hmat[0][2]*cell.hmat[0][2]+cell.hmat[1][2]*cell.hmat[1][2]+cell.hmat[2][2]*cell.hmat[2][2])/ang;
  cell1x=cell.hmat[0][0]/ang;cell1y=cell.hmat[1][0]/ang;cell1z=cell.hmat[2][0]/ang;
  cell2x=cell.hmat[0][1]/ang;cell2y=cell.hmat[1][1]/ang;cell2z=cell.hmat[2][1]/ang;
  cell3x=cell.hmat[0][2]/ang;cell3y=cell.hmat[1][2]/ang;cell3z=cell.hmat[2][2]/ang;
  // fra or abs line
  getline (finconfig, line); remove_carriage_return(line);
  getline (finconfig, line); remove_carriage_return(line);
  //  std::cout << "LINE = " << "\"" << line << "\"" << std::endl;
  get_first_arg(line, arg1); get_first_arg(line, arg2);
  if ((arg1 == "abs")||(arg1 == "ABS")) { abs = true; }
  //  std::cout << "ARG1 = " << "\"" << arg1 << "\"" << std::endl;
  //  std::cout << abs << std::endl;
  //  std::cout << "ARG2 = " << "\"" << arg2 << "\"" << std::endl;

  /*
  if (abs) {
    for (int i=1; i<=atom.natom; i++) {
      //      finconfig >> xx;  finconfig >> yy;   finconfig >> zz;
      if (i<=atom.natom-natoms) {
	getline(finconfig, line);
	get_first_arg(line, arg1); get_first_arg(line, arg2); get_first_arg(line, arg3);
	strcpy (larg1, arg1.c_str()); strcpy (larg2, arg2.c_str()); strcpy (larg3, arg3.c_str()); 
	xx = atof(larg1); yy = atof(larg2); zz = atof(larg3);
	atom.rx[i] = xx*1.0e-10; atom.ry[i] = yy*1.0e-10; atom.rz[i] = zz*1.0e-10;
      } else {
	atom.rx[i] = rxs[i-atom.natom+natoms]/hmats[0][0]*cell.hmat[0][0];
	atom.ry[i] = rys[i-atom.natom+natoms]/hmats[1][1]*cell.hmat[1][1];
	atom.rz[i] = rzs[i-atom.natom+natoms]/hmats[2][2]*cell.hmat[2][2];
      }
      atom.rx_org[i] = atom.rx[i]; atom.ry_org[i] = atom.ry[i]; atom.rz_org[i] = atom.rz[i]; 
      atom.vx[i] = 0.0; atom.vy[i] = 0.0; atom.vz[i] = 0.0;
    }
  } else {
  */
  for (int i=1; i<=atom.natom; i++) {
    if (i<=atom.natom-natoms) {
      getline(finconfig, line); remove_carriage_return(line);
      if (count_arg_number(line)==3) {
	if (strcmp(species,"Mu") == 0) { printf("Error in CONFIG. Specify atom type.\n"); return; }
	get_first_arg(line, arg1); get_first_arg(line, arg2); get_first_arg(line, arg3);
	strcpy (larg1, arg1.c_str()); strcpy (larg2, arg2.c_str()); strcpy (larg3, arg3.c_str()); 
	xx = atof(larg1); yy = atof(larg2); zz = atof(larg3);
	atom.mfx[i] = false; atom.mfy[i] = false; atom.mfz[i] = false; 
      } else if (count_arg_number(line)==4) {
	get_first_arg(line, arg1); get_first_arg(line, arg2); get_first_arg(line, arg3); get_first_arg(line, arg4);
	strcpy (larg1, arg1.c_str()); strcpy (larg2, arg2.c_str()); strcpy (larg3, arg3.c_str());
	strcpy (atom.asp[i], arg4.c_str());
	xx = atof(larg1); yy = atof(larg2); zz = atof(larg3);
	atom.mfx[i] = false; atom.mfy[i] = false; atom.mfz[i] = false; 
      } else if (count_arg_number(line)==5) { // with group #
	get_first_arg(line, arg1); get_first_arg(line, arg2); get_first_arg(line, arg3); get_first_arg(line, arg4);
  get_first_arg(line, arg5);
	strcpy (larg1, arg1.c_str()); strcpy (larg2, arg2.c_str()); strcpy (larg3, arg3.c_str());
	strcpy (atom.asp[i], arg4.c_str()); strcpy (larg8, arg5.c_str());
	xx = atof(larg1); yy = atof(larg2); zz = atof(larg3);
	atom.mfx[i] = false; atom.mfy[i] = false; atom.mfz[i] = false; 
  atom.group[i] = atoi(larg8);
      } else if (count_arg_number(line)==6) {
	if (strcmp(species,"Mu") == 0) { printf("Error in CONFIG. Specify atom type.\n"); return; }
	get_first_arg(line, arg1); get_first_arg(line, arg2); get_first_arg(line, arg3); get_first_arg(line, arg4);
	get_first_arg(line, arg5); get_first_arg(line, arg6);
	strcpy (larg1, arg1.c_str()); strcpy (larg2, arg2.c_str()); strcpy (larg3, arg3.c_str()); 
	strcpy (larg4, arg4.c_str()); strcpy (larg5, arg5.c_str()); strcpy (larg6, arg6.c_str());
	xx = atof(larg1); yy = atof(larg2); zz = atof(larg3);
	atom.mfx[i] = false; atom.mfy[i] = false; atom.mfz[i] = false; 
	if (strcmp(larg4,"F") == 0) { atom.mfx[i] = true; }
	if (strcmp(larg5,"F") == 0) { atom.mfy[i] = true; }
	if (strcmp(larg6,"F") == 0) { atom.mfz[i] = true; }
      } else if (count_arg_number(line)==8) { // with group #
	get_first_arg(line, arg1); get_first_arg(line, arg2); get_first_arg(line, arg3); get_first_arg(line, arg4);
	get_first_arg(line, arg5); get_first_arg(line, arg6); get_first_arg(line, arg7); get_first_arg(line, arg8);
	strcpy (larg1, arg1.c_str()); strcpy (larg2, arg2.c_str()); strcpy (larg3, arg3.c_str()); 
	strcpy (larg5, arg5.c_str()); strcpy (larg6, arg6.c_str()); strcpy (larg7, arg7.c_str()); 
	strcpy (atom.asp[i], arg4.c_str()); strcpy (larg8, arg8.c_str());
	xx = atof(larg1); yy = atof(larg2); zz = atof(larg3);
	atom.mfx[i] = false; atom.mfy[i] = false; atom.mfz[i] = false; 
	if (strcmp(larg5,"F") == 0) { atom.mfx[i] = true; }
	if (strcmp(larg6,"F") == 0) { atom.mfy[i] = true; }
	if (strcmp(larg7,"F") == 0) { atom.mfz[i] = true; }
  atom.group[i] = atoi(larg8);
      } else if (count_arg_number(line)==7) {
	get_first_arg(line, arg1); get_first_arg(line, arg2); get_first_arg(line, arg3); get_first_arg(line, arg4);
	get_first_arg(line, arg5); get_first_arg(line, arg6); get_first_arg(line, arg7);
	strcpy (larg1, arg1.c_str()); strcpy (larg2, arg2.c_str()); strcpy (larg3, arg3.c_str()); 
	strcpy (larg5, arg5.c_str()); strcpy (larg6, arg6.c_str()); strcpy (larg7, arg7.c_str()); 
	strcpy (atom.asp[i], arg4.c_str());
	xx = atof(larg1); yy = atof(larg2); zz = atof(larg3);
	atom.mfx[i] = false; atom.mfy[i] = false; atom.mfz[i] = false; 
	if (strcmp(larg5,"F") == 0) { atom.mfx[i] = true; }
	if (strcmp(larg6,"F") == 0) { atom.mfy[i] = true; }
	if (strcmp(larg7,"F") == 0) { atom.mfz[i] = true; }
      } else {
	std::cout << "Error in reading config: Line has "<<count_arg_number(line)<<" args"<<std::endl;
      }
      if (abs) {
	atom.rx[i] = xx*1.0e-10; atom.ry[i] = yy*1.0e-10; atom.rz[i] = zz*1.0e-10;
      } else {
	atom.rx[i] = cell.hmat[0][0]*xx + cell.hmat[0][1]*yy + cell.hmat[0][2]*zz;
	atom.ry[i] = cell.hmat[1][0]*xx + cell.hmat[1][1]*yy + cell.hmat[1][2]*zz;
	atom.rz[i] = cell.hmat[2][0]*xx + cell.hmat[2][1]*yy + cell.hmat[2][2]*zz;
      }
    } else {
      atom.rx[i] = rxs[i-atom.natom+natoms]/hmats[0][0]*cell.hmat[0][0];
      atom.ry[i] = rys[i-atom.natom+natoms]/hmats[1][1]*cell.hmat[1][1];
      atom.rz[i] = rzs[i-atom.natom+natoms]/hmats[2][2]*cell.hmat[2][2];
      atom.mfx[i] = false; atom.mfy[i] = false; atom.mfz[i] = false; 
    }	
    atom.rx_org[i] = atom.rx[i]; atom.ry_org[i] = atom.ry[i]; atom.rz_org[i] = atom.rz[i]; 
    atom.vx[i] = 0.0; atom.vy[i] = 0.0; atom.vz[i] = 0.0;
    atom.ax[i] = 0.0; atom.ay[i] = 0.0; atom.az[i] = 0.0;
    atom.bx[i] = 0.0; atom.by[i] = 0.0; atom.bz[i] = 0.0;
    atom.cx[i] = 0.0; atom.cy[i] = 0.0; atom.cz[i] = 0.0;
  }
  // Atom weight
  //atom.wm = 1.6726e-27*63.55;
  set_atom_color();
  set_atom_weight();
  for (int i=1; i<=atom.natom; i++) {
    atom.anum[i] = atom_number(atom.asp[i]);
  }
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      cell.hmat_org[i][j]=cell.hmat[i][j];
      cell.hvmat[i][j] = 0; cell.hamat[i][j] = 0; cell.hbmat[i][j] = 0;
      cell.hcmat[i][j] = 0; cell.sgmmat_set[i][j] = 0;
    }
  }
  stress_set();
  //istep = 0;
  ex = 0.0;
  //out_of_cell();
  if (incell) { out_of_cell(); }
  chdir(cwdname);
  if (mode==0) {
    pot_initialize_all();
    bookkeep();
    potential();
    change_atom_color();
    // For output in "Status" area
    f_max=atom.Fmax()/eV*ang;
    epotatom=atom.epotsum/eV/atom.natom;
  }
  if (initial_velocity) { velocity_random(temp_set); }

  finconfig.close(); finconfig.clear();
  if ((mode==0)&&(imerge==1)) { delete[] rxs; delete[] rys; delete[] rzs; }
}

std::string replace_string(std::string string1, std::string string2, std::string string3)
{
  std::string::size_type pos(string1.find(string2));
  while (pos != std::string::npos)
    {
      string1.replace(pos, string2.length(), string3);
      pos = string1.find(string2, pos+string3.length());
    }
  return string1;
}

void get_first_arg(std::string &line, std::string &arg1)
{
  line = replace_string(line, "\t", " ");
  int istt = line.find_first_not_of(" ");
  if (istt == -1) {
    arg1 = "";
  } else {
    line = line.substr(istt);
    int iend = line.find_first_of(" ");
    if (iend == -1) { iend = line.length(); }
    //std::cout<<istt<< " " <<iend<<std::endl;
    arg1 = line.substr(0,iend);
    //std::cout<<line <<" " <<arg1 <<std::endl;
    line = line.substr(iend);
  }
}
void get_first_arg(char* line, char* arg1)
{
  if (strchr(line,' ')==NULL) { strcpy(arg1,line);
  } else {
    int l=strchr(line,' ')-line; int len=strlen(line);
    strncpy(arg1,line,l);
    for (int i=l+1; i<=len; i++) {
      if (line[i]==' ') { l++; }
    }
    for (int i=l+1; i<=len; i++) {
      line[i-l-1]=line[i];
    }
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

int remove_after_sharp(std::string &line)
{
  int shp = line.find("#");
  if (shp>0) {
    line = line.substr(0,shp);
    return 1;
  } else {
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

void set_atom_weight()
{
  double am = 1.6726485e-27;
  for (int i=1; i<=atom.natom; i++) {
    if (strcmp(atom.potential_func,"Shellmodel")==0) {
      if ((strcmp(atom.asp[i],"Pb")==0)||(strcmp(atom.asp[i],"pb")==0)) {
	if (i<=atom.natom/2) { atom.wm[i]=am*207.2; // Core
	} else { atom.wm[i]=am*207.2; } // Shell
      } else if ((strcmp(atom.asp[i],"Ti")==0)||(strcmp(atom.asp[i],"ti")==0)) {
	if (i<=atom.natom/2) { atom.wm[i]=am*47.9; // Core
	} else { atom.wm[i]=am*47.9; } // Shell
      } else if ((strcmp(atom.asp[i],"O")==0)||(strcmp(atom.asp[i],"o")==0)) {
	if (i<=atom.natom/2) { atom.wm[i]=am*15.9994; // Core
	} else { atom.wm[i]=am*15.9994; } // Shell
      }
    } else {
      if ((strcmp(atom.asp[i],"C")==0)||(strcmp(atom.asp[i],"c")==0)) {
	atom.wm[i]=am*12.01;
      } else if ((strcmp(atom.asp[i],"Pb")==0)||(strcmp(atom.asp[i],"pb")==0)) {
	atom.wm[i]=am*207.2;
      } else if ((strcmp(atom.asp[i],"Mo")==0)||(strcmp(atom.asp[i],"mo")==0)) {
	atom.wm[i]=am*95.94;
      } else if ((strcmp(atom.asp[i],"Ta")==0)||(strcmp(atom.asp[i],"ta")==0)) {
	atom.wm[i]=am*180.95;
      } else if ((strcmp(atom.asp[i],"Mg")==0)||(strcmp(atom.asp[i],"mg")==0)) {
	atom.wm[i]=am*24.305;
      } else if ((strcmp(atom.asp[i],"Co")==0)||(strcmp(atom.asp[i],"co")==0)) {
	atom.wm[i]=am*58.93;
      } else if ((strcmp(atom.asp[i],"Ti")==0)||(strcmp(atom.asp[i],"ti")==0)) {
	atom.wm[i]=am*47.9;
      } else if ((strcmp(atom.asp[i],"Ba")==0)||(strcmp(atom.asp[i],"ba")==0)) {
	atom.wm[i]=am*137.34;
      } else if ((strcmp(atom.asp[i],"Zr")==0)||(strcmp(atom.asp[i],"zr")==0)) {
	atom.wm[i]=am*91.22;
      } else if ((strcmp(atom.asp[i],"Cu")==0)||(strcmp(atom.asp[i],"cu")==0)) {
	atom.wm[i]=am*63.55;
      } else if ((strcmp(atom.asp[i],"Al")==0)||(strcmp(atom.asp[i],"al")==0)) {
	atom.wm[i]=am*26.982;
      } else if ((strcmp(atom.asp[i],"Si")==0)||(strcmp(atom.asp[i],"si")==0)) {
	atom.wm[i]=am*28.086;
      } else if ((strcmp(atom.asp[i],"Ge")==0)||(strcmp(atom.asp[i],"ge")==0)) {
	atom.wm[i]=am*72.59;
      } else if ((strcmp(atom.asp[i],"Cr")==0)||(strcmp(atom.asp[i],"cr")==0)) {
	atom.wm[i]=am*52.00;
      } else if ((strcmp(atom.asp[i],"Fe")==0)||(strcmp(atom.asp[i],"fe")==0)) {
	atom.wm[i]=am*55.85;
      } else if ((strcmp(atom.asp[i],"Ni")==0)||(strcmp(atom.asp[i],"ni")==0)) {
	atom.wm[i]=am*58.71;
      } else if ((strcmp(atom.asp[i],"Ag")==0)||(strcmp(atom.asp[i],"ag")==0)) {
	atom.wm[i]=am*107.87;
      } else if ((strcmp(atom.asp[i],"Au")==0)||(strcmp(atom.asp[i],"au")==0)) {
	atom.wm[i]=am*196.97;
      } else if ((strcmp(atom.asp[i],"Pd")==0)||(strcmp(atom.asp[i],"pd")==0)) {
	atom.wm[i]=am*106.40;
      } else if ((strcmp(atom.asp[i],"Sn")==0)||(strcmp(atom.asp[i],"sn")==0)) {
	atom.wm[i]=am*118.69;
      } else if ((strcmp(atom.asp[i],"W")==0)||(strcmp(atom.asp[i],"w")==0)) {
	atom.wm[i]=am*183.85;
      } else if ((strcmp(atom.asp[i],"Pt")==0)||(strcmp(atom.asp[i],"pt")==0)) {
	atom.wm[i]=am*195.09;
      } else if ((strcmp(atom.asp[i],"H")==0)||(strcmp(atom.asp[i],"h")==0)) {
	atom.wm[i]=am*1.0079;
      } else if ((strcmp(atom.asp[i],"O")==0)||(strcmp(atom.asp[i],"o")==0)) {
	atom.wm[i]=am*15.999;
      } else if ((strcmp(atom.asp[i],"Y")==0)||(strcmp(atom.asp[i],"y")==0)) {
	atom.wm[i]=am*88.91;
      } else if ((strcmp(atom.asp[i],"Nd")==0)||(strcmp(atom.asp[i],"nd")==0)) {
	atom.wm[i]=am*144.24;
      } else if ((strcmp(atom.asp[i],"B")==0)||(strcmp(atom.asp[i],"b")==0)) {
	atom.wm[i]=am*10.81;
      } else if ((strcmp(atom.asp[i],"N")==0)||(strcmp(atom.asp[i],"n")==0)) {
	atom.wm[i]=am*14.01;
      } else if ((strcmp(atom.asp[i],"Ga")==0)||(strcmp(atom.asp[i],"ga")==0)) {
	atom.wm[i]=am*69.72;
      } else {
	printf("WARNING: Weight data not found for i= %d, %s\n",i,atom.asp[i]); atom.wm[i] = am*10.0; } 
    }
  }
    /*
      wmave=0.0d0
      do i=1,n
      wmave=wmave+atom.wm[i]
      enddo
      wmave=wmave/dble(n)
    */
}

