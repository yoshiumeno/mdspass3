#include <string>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "myheader.h"

#define NONEINT -99999
#define NONEDOUBLE -1.0e30

extern int radius, segments, show_cnt_wall, show_cnt_wallv, show_cnt_wall_num, show_cnt_ring;
extern int show_axis, show_cell, relax_algo, relax_accel;
extern int draw_bond, draw_force, draw_load, draw_aux;
extern int draw_bond_pbc, bond_display_type, bond_display_color;
extern int mode_cnt_corrugation, read_cntwall;
extern float bondlength, dtm, vscl_force;
extern int mdspeed, itolfor, itolstep, tolstep, irecipe, inogui;
extern int itolstress;
extern float tolstress;
extern int cnt_load_algo;
extern float cnt_pressure;
extern float tolfor;
extern int color_mode, color_mode_auto;
extern float fp_alph_ini, fp_ffinc, fp_ffdec, fp_ffalph, fp_ffdtmax;
extern int fp_nfmin;
extern int cellfix_xx, cellfix_yy, cellfix_zz, cellfix_xy, cellfix_yz, cellfix_zx;
extern float prdamper_val1, prmass_scale, nhmass_scale;
extern float rotmat[16], obj_pos[3], scl;
extern float strs_set_xx, strs_set_yy, strs_set_zz, strs_set_xy, strs_set_yz, strs_set_zx;
extern int cell_relax_rep; extern float cell_relax_tolfor;
extern int enable_atom_click;

void get_first_arg(std::string &line, std::string &arg1);
void stress_set();
void remove_carriage_return(std::string& line);

double find_double(const char* fname, const char* tag)
{
  std::ifstream fin(fname);
  std::string line, arg1, arg2;
  double arg = NONEDOUBLE;
  for (int i=1; i<=100; i++) {
    getline (fin, line); remove_carriage_return(line);
    get_first_arg(line, arg1); get_first_arg(line, arg2);
    if (strcmp(arg1.c_str(), tag)==0) {
      arg = (double)atof(arg2.c_str());
    }
  }
  fin.close(); fin.clear();
  return arg;
}

int find_int(const char* fname, const char* tag)
{
  std::ifstream fin(fname);
  std::string line, arg1, arg2;
  int arg = NONEINT;
  for (int i=1; i<=100; i++) {
    getline (fin, line); remove_carriage_return(line);
    get_first_arg(line, arg1); get_first_arg(line, arg2);
    if (strcmp(arg1.c_str(), tag)==0) {
      arg = atoi(arg2.c_str());
    }
  }
  fin.close(); fin.clear();
  return arg;
}

int find_bool(const char* fname, const char* tag)
{
  std::ifstream fin(fname);
  std::string line, arg1, arg2;
  int arg = -1;
  for (int i=1; i<=100; i++) {
    getline (fin, line); remove_carriage_return(line);
    get_first_arg(line, arg1); get_first_arg(line, arg2);
    if (strcmp(arg1.c_str(), tag)==0) {
      if (strcmp(arg2.c_str(), "Yes")==0) {
	arg = 1;
      } else if (strcmp(arg2.c_str(), "yes")==0) {
	arg = 1;
      } else if (strcmp(arg2.c_str(), "YES")==0) {
	arg = 1;
      } else if (strcmp(arg2.c_str(), "No")==0) {
	arg = 0;
      } else if (strcmp(arg2.c_str(), "no")==0) {
	arg = 0;
      } else if (strcmp(arg2.c_str(), "NO")==0) {
	arg = 0;
      }
    }
  }
  fin.close(); fin.clear();
  return arg;
}

void get_dat(const char* fname, const char* tag, double &val)
{
  if (find_double(fname, tag)!=NONEDOUBLE) {
    val = find_double(fname, tag); printf("# %s = %e (double)\n", tag, val); } 
}
void get_dat(const char* fname, const char* tag, float &val)
{
  if (find_double(fname, tag)!=NONEDOUBLE) {
    val = (float)find_double(fname, tag); printf("# %s = %e (float)\n", tag, val); } 
}
void get_dat(const char* fname, const char* tag, bool &val)
{
  if (find_bool(fname, tag) > -1 ) {
    val = find_bool(fname, tag);
    if (val == true) { printf("# %s = yes\n", tag); }
    else { printf("# %s = no\n", tag); } } 
}
void get_dat(const char* fname, const char* tag, int &val)
{
  if (find_bool(fname, tag) > -1) {
    val = find_bool(fname, tag);
    if (val == true) { printf("# %s = yes\n", tag); }
    else { printf("# %s = no\n", tag); }
  } else if (find_int(fname, tag)!=NONEINT) {
    val = find_int(fname, tag); printf("# %s = %d\n", tag, val);
  }
}

void readsetdat(const char* fname)
{
  printf("### Read %s ###\n",fname);
  get_dat(fname, "dt", dt);
  if (log10(dt)<-13) { printf(" ! You gave dt in [s]. Next time give it in [fs].\n"); }
  else { dt *= 1e-15;  }
  if (find_double(fname, "frc")!=NONEDOUBLE) {
    get_dat(fname, "frc", book.frc);
    if (log10(book.frc)<-8) { printf(" ! You gave frc in [m]. Next time give it in [A].\n"); }
    else { book.frc *= 1e-10; }
    book.frc2 = book.frc*book.frc; }
  get_dat(fname, "frc_margin", frcmar); if (frcmar<0) { frcmar = 0.2; }
  get_dat(fname, "nbk", book.nbk);
  get_dat(fname, "pbcx", cell.pbcx);
  get_dat(fname, "pbcy", cell.pbcy);
  get_dat(fname, "pbcz", cell.pbcz);
  get_dat(fname, "cellfix_xx", cellfix_xx);
  get_dat(fname, "cellfix_yy", cellfix_yy);
  get_dat(fname, "cellfix_zz", cellfix_zz);
  get_dat(fname, "cellfix_xy", cellfix_xy);
  get_dat(fname, "cellfix_yz", cellfix_yz);
  get_dat(fname, "cellfix_zx", cellfix_zx);
  get_dat(fname, "radius", radius);
  get_dat(fname, "segment", segments);
  get_dat(fname, "draw_bond", draw_bond);
  get_dat(fname, "draw_bond_pbc", draw_bond_pbc);
  get_dat(fname, "bond_type", bond_display_type);
  get_dat(fname, "bond_color", bond_display_color);
  get_dat(fname, "draw_force", draw_force);
  get_dat(fname, "draw_load", draw_load);
  get_dat(fname, "draw_aux", draw_aux);
  get_dat(fname, "forcelength", vscl_force);
  get_dat(fname, "bondlength", bondlength);
  //draw_bond = find_bool(fname, "draw_bond");
  //draw_force = find_bool(fname, "draw_force");
  //vscl_force = (float)find_double(fname, "forcelength");
  //bondlength = (float)find_double(fname, "bondlength");
  get_dat(fname, "show_cell", show_cell);
  get_dat(fname, "show_axis", show_axis);
  get_dat(fname, "no_trans", notrans);
  get_dat(fname, "keep_in_cell", incell);
  get_dat(fname, "temp", temp_set);
  get_dat(fname, "show_cnt_wall", show_cnt_wall);
  get_dat(fname, "show_cnt_wallv", show_cnt_wallv);
  get_dat(fname, "show_cnt_wall_num", show_cnt_wall_num);
  get_dat(fname, "show_cnt_ring", show_cnt_ring);
  get_dat(fname, "cnt_corrugation", mode_cnt_corrugation);
  get_dat(fname, "cnt_pressure", cnt_pressure);
  get_dat(fname, "read_cnt_wall", read_cntwall);
  get_dat(fname, "ffdtmax", fire.ffdtmax);
  if (find_double(fname, "ffdtmax")!=NONEDOUBLE) {
    if (log10(fire.ffdtmax)<-13) { printf(" ! You gave ffdtmax in [s]. Next time give it in [fs].\n"); }
    else { fire.ffdtmax *= 1e-15;  } }
  get_dat(fname, "fire_dtmax", fp_ffdtmax);
  //if (find_double(fname, "fire_dtmax")!=NONEDOUBLE) {
  //  if (log10(fp_ffdtmax)<-10) { printf(" ! You gave ffdtmax in [s]. Next time give it in [fs].\n"); }
  //  else { fp_ffdtmax *= 1e-15;  } }
  get_dat(fname, "fire_alph_ini", fp_alph_ini);
  get_dat(fname, "fire_inc", fp_ffinc);
  get_dat(fname, "fire_dec", fp_ffdec);
  get_dat(fname, "fire_alph_rate", fp_ffalph);
  get_dat(fname, "fire_nfmin", fp_nfmin);
  get_dat(fname, "mdmotion", mdmotion);
  get_dat(fname, "ensemble", ensemble);
  get_dat(fname, "recipe", irecipe);
  get_dat(fname, "nogui", inogui);
  get_dat(fname, "tersoff_nocutoff", tersoff.nocutoff);
  get_dat(fname, "initial_step", istep);
  /*
  show_cell = find_bool(fname, "show_cell");
  show_axis = find_bool(fname, "show_axis");
  notrans = find_bool(fname, "no_trans");
  temp_set = find_double(fname, "temp");
  show_cnt_wall = find_bool(fname, "show_cnt_wall");
  show_cnt_wallv = find_bool(fname, "show_cnt_wallv");
  show_cnt_wall_num = find_int(fname, "show_cnt_wall_num");
  mode_cnt_corrugation = find_bool(fname, "cnt_corrugation");
  fire.ffdtmax = find_double(fname, "ffdtmax");
  mdmotion = find_bool(fname, "mdmotion");
  ensemble = find_int(fname, "ensemble");
  */
  get_dat(fname, "dexdt", dexdt);
  get_dat(fname, "deydt", deydt);
  get_dat(fname, "dezdt", dezdt);
  get_dat(fname, "repeat_lz", repeat_lz);
  get_dat(fname, "repeat_lz_min", repeat_lz_min);
  get_dat(fname, "repeat_lz_max", repeat_lz_max);
  get_dat(fname, "yzplane_punch", yzplane_punch);
  get_dat(fname, "yzplane_punch_d", yzplane_punch_d);
  get_dat(fname, "yzplane_punch_dd", yzplane_punch_dd);
  get_dat(fname, "relax_algo", relax_algo);
  get_dat(fname, "relax_gloc_accel", relax_accel);
  get_dat(fname, "redraw_interval", mdspeed);
  get_dat(fname, "cnt_load_algo", cnt_load_algo);
  get_dat(fname, "stop_force", itolfor);
  get_dat(fname, "stop_force_val", tolfor);
  get_dat(fname, "stop_step", itolstep);
  get_dat(fname, "stop_step_num", tolstep);
  get_dat(fname, "stop_stress", itolstress);
  get_dat(fname, "stop_stress_val", tolstress);
  get_dat(fname, "color_mode", color_mode);
  get_dat(fname, "color_mode_auto", color_mode_auto);
  get_dat(fname, "prdamper", prdamper_val1);
  get_dat(fname, "prmass", prmass_scale);
  get_dat(fname, "nhmass", nhmass_scale);
  get_dat(fname, "initial_velocity", initial_velocity);
  get_dat(fname, "stress_set_xx", strs_set_xx);
  get_dat(fname, "stress_set_yy", strs_set_yy);
  get_dat(fname, "stress_set_zz", strs_set_zz);
  get_dat(fname, "stress_set_xy", strs_set_xy);
  get_dat(fname, "stress_set_yz", strs_set_yz);
  get_dat(fname, "stress_set_zx", strs_set_zx);
  get_dat(fname, "cell_relax_times", cell_relax_rep);
  get_dat(fname, "cell_relax_force", cell_relax_tolfor);
  get_dat(fname, "enable_atom_click", enable_atom_click);
  get_dat(fname, "grab_1", ifgrab1);
  get_dat(fname, "grab_1_group", grab1num);
  get_dat(fname, "grab_1_dxdt", grabdxdt1);
  get_dat(fname, "grab_1_dydt", grabdydt1);
  get_dat(fname, "grab_1_dzdt", grabdzdt1);
  get_dat(fname, "grab_2", ifgrab2);
  get_dat(fname, "grab_2_group", grab2num);
  get_dat(fname, "grab_2_dxdt", grabdxdt2);
  get_dat(fname, "grab_2_dydt", grabdydt2);
  get_dat(fname, "grab_2_dzdt", grabdzdt2);
  get_dat(fname, "push_1", ifpush1);
  get_dat(fname, "push_1_group", push1num);
  get_dat(fname, "push_1_fx", pushfx1);
  get_dat(fname, "push_1_fy", pushfy1);
  get_dat(fname, "push_1_fz", pushfz1);
  get_dat(fname, "push_2", ifpush2);
  get_dat(fname, "push_2_group", push2num);
  get_dat(fname, "push_2_fx", pushfx2);
  get_dat(fname, "push_2_fy", pushfy2);
  get_dat(fname, "push_2_fz", pushfz2);
  /*
  repeat_lz = find_bool(fname, "repeat_lz");
  repeat_lz_min = (float)find_double(fname, "repeat_lz_min");
  repeat_lz_max = (float)find_double(fname, "repeat_lz_max");
  relax_algo = find_int(fname, "relax_algo");
  */
  printf("### Read %s end ###\n",fname);

  dtm = (float)dt*1e15;

  if (cellfix_xx) { cell.fix[0][0] = true; }
  if (cellfix_yy) { cell.fix[1][1] = true; }
  if (cellfix_zz) { cell.fix[2][2] = true; }
  if (cellfix_xy) { cell.fix[0][1] = true; cell.fix[1][0] = true; }
  if (cellfix_yz) { cell.fix[1][2] = true; cell.fix[2][1] = true; }
  if (cellfix_zx) { cell.fix[2][0] = true; cell.fix[0][2] = true; }

  stress_set();

  /*
  char string[80];
  FILE *fp = fopen("setdat","r");
  fgets(string, 80, fp);
  printf("%s\n",string);
  */
}
