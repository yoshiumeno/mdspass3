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

void get_first_arg(std::string &line, std::string &arg1);
void stress_set();

void write_dat(const char* fname, const char* tag, double val);
void write_dat(const char* fname, const char* tag, float val);
void write_dat(const char* fname, const char* tag, int val);


void writesetdat(const char* fname)
{
  printf("### Write Settings to %s ###\n",fname);
  std::ofstream fout(fname, std::ios::out);
  fout << "#" << std::endl;
  fout << "# Saved settings" << std::endl;
  fout << "#" << std::endl;
  fout.close(); fout.clear();
  
  write_dat(fname, "dt", dt);
  if (log10(dt)<-13) { printf(" ! You gave dt in [s]. Next time give it in [fs].\n"); }
  else { dt *= 1e-15;  }
  book.frc2 = book.frc*book.frc;
  write_dat(fname, "frc_margin", frcmar); if (frcmar<0) { frcmar = 0.2; }
  write_dat(fname, "nbk", book.nbk);
  write_dat(fname, "pbcx", cell.pbcx);
  write_dat(fname, "pbcy", cell.pbcy);
  write_dat(fname, "pbcz", cell.pbcz);
  write_dat(fname, "cellfix_xx", cellfix_xx);
  write_dat(fname, "cellfix_yy", cellfix_yy);
  write_dat(fname, "cellfix_zz", cellfix_zz);
  write_dat(fname, "cellfix_xy", cellfix_xy);
  write_dat(fname, "cellfix_yz", cellfix_yz);
  write_dat(fname, "cellfix_zx", cellfix_zx);
  write_dat(fname, "radius", radius);
  write_dat(fname, "segment", segments);
  write_dat(fname, "draw_bond", draw_bond);
  write_dat(fname, "draw_bond_pbc", draw_bond_pbc);
  write_dat(fname, "bond_type", bond_display_type);
  write_dat(fname, "bond_color", bond_display_color);
  write_dat(fname, "draw_force", draw_force);
  write_dat(fname, "draw_load", draw_load);
  write_dat(fname, "draw_aux", draw_aux);
  write_dat(fname, "forcelength", vscl_force);
  write_dat(fname, "bondlength", bondlength);
  write_dat(fname, "show_cell", show_cell);
  write_dat(fname, "show_axis", show_axis);
  write_dat(fname, "no_trans", notrans);
  write_dat(fname, "keep_in_cell", incell);
  write_dat(fname, "temp", temp_set);
  write_dat(fname, "show_cnt_wall", show_cnt_wall);
  write_dat(fname, "show_cnt_wallv", show_cnt_wallv);
  write_dat(fname, "show_cnt_wall_num", show_cnt_wall_num);
  write_dat(fname, "show_cnt_ring", show_cnt_ring);
  write_dat(fname, "cnt_corrugation", mode_cnt_corrugation);
  write_dat(fname, "cnt_pressure", cnt_pressure);
  write_dat(fname, "read_cnt_wall", read_cntwall);
  write_dat(fname, "ffdtmax", fire.ffdtmax);
  if (log10(fire.ffdtmax)<-13) { printf(" ! You gave ffdtmax in [s]. Next time give it in [fs].\n"); }
  else { fire.ffdtmax *= 1e-15;  }
  write_dat(fname, "fire_dtmax", fp_ffdtmax);
  write_dat(fname, "fire_alph_ini", fp_alph_ini);
  write_dat(fname, "fire_inc", fp_ffinc);
  write_dat(fname, "fire_dec", fp_ffdec);
  write_dat(fname, "fire_alph_rate", fp_ffalph);
  write_dat(fname, "fire_nfmin", fp_nfmin);
  write_dat(fname, "mdmotion", mdmotion);
  write_dat(fname, "ensemble", ensemble);
  write_dat(fname, "recipe", irecipe);
  write_dat(fname, "nogui", inogui);
  write_dat(fname, "tersoff_nocutoff", tersoff.nocutoff);
  write_dat(fname, "initial_step", istep);
  write_dat(fname, "dexdt", dexdt);
  write_dat(fname, "deydt", deydt);
  write_dat(fname, "dezdt", dezdt);
  write_dat(fname, "repeat_lz", repeat_lz);
  write_dat(fname, "repeat_lz_min", repeat_lz_min);
  write_dat(fname, "repeat_lz_max", repeat_lz_max);
  write_dat(fname, "yzplane_punch", yzplane_punch);
  write_dat(fname, "yzplane_punch_d", yzplane_punch_d);
  write_dat(fname, "yzplane_punch_dd", yzplane_punch_dd);
  write_dat(fname, "relax_algo", relax_algo);
  write_dat(fname, "relax_gloc_accel", relax_accel);
  write_dat(fname, "redraw_interval", mdspeed);
  write_dat(fname, "cnt_load_algo", cnt_load_algo);
  write_dat(fname, "stop_force", itolfor);
  write_dat(fname, "stop_force_val", tolfor);
  write_dat(fname, "stop_step", itolstep);
  write_dat(fname, "stop_step_num", tolstep);
  write_dat(fname, "stop_stress", itolstress);
  write_dat(fname, "stop_stress_val", tolstress);
  write_dat(fname, "color_mode", color_mode);
  write_dat(fname, "color_mode_auto", color_mode_auto);
  write_dat(fname, "prdamper", prdamper_val1);
  write_dat(fname, "prmass", prmass_scale);
  write_dat(fname, "nhmass", nhmass_scale);
  write_dat(fname, "initial_velocity", initial_velocity);
  write_dat(fname, "stress_set_xx", strs_set_xx);
  write_dat(fname, "stress_set_yy", strs_set_yy);
  write_dat(fname, "stress_set_zz", strs_set_zz);
  write_dat(fname, "stress_set_xy", strs_set_xy);
  write_dat(fname, "stress_set_yz", strs_set_yz);
  write_dat(fname, "stress_set_zx", strs_set_zx);
  write_dat(fname, "cell_relax_times", cell_relax_rep);
  write_dat(fname, "cell_relax_force", cell_relax_tolfor);
  write_dat(fname, "enable_atom_click", enable_atom_click);
  write_dat(fname, "grab_1", ifgrab1);
  write_dat(fname, "grab_1_group", grab1num);
  write_dat(fname, "grab_1_dxdt", grabdxdt1);
  write_dat(fname, "grab_1_dydt", grabdydt1);
  write_dat(fname, "grab_1_dzdt", grabdzdt1);
  write_dat(fname, "grab_2", ifgrab2);
  write_dat(fname, "grab_2_group", grab2num);
  write_dat(fname, "grab_2_dxdt", grabdxdt2);
  write_dat(fname, "grab_2_dydt", grabdydt2);
  write_dat(fname, "grab_2_dzdt", grabdzdt2);
  write_dat(fname, "push_1", ifpush1);
  write_dat(fname, "push_1_group", push1num);
  write_dat(fname, "grab_1_dxdt", pushfx1);
  write_dat(fname, "grab_1_dydt", pushfy1);
  write_dat(fname, "grab_1_dzdt", pushfz1);
  write_dat(fname, "push_2", ifpush2);
  write_dat(fname, "push_2_group", push2num);
  write_dat(fname, "grab_2_dxdt", pushfx2);
  write_dat(fname, "grab_2_dydt", pushfy2);
  write_dat(fname, "grab_2_dzdt", pushfz2);

  printf("### Read %s end ###\n",fname);

  dtm = (float)dt*1e15;

  if (cellfix_xx) { cell.fix[0][0] = true; }
  if (cellfix_yy) { cell.fix[1][1] = true; }
  if (cellfix_zz) { cell.fix[2][2] = true; }
  if (cellfix_xy) { cell.fix[0][1] = true; cell.fix[1][0] = true; }
  if (cellfix_yz) { cell.fix[1][2] = true; cell.fix[2][1] = true; }
  if (cellfix_zx) { cell.fix[2][0] = true; cell.fix[0][2] = true; }

  stress_set();

}

void write_dat(const char* fname, const char* tag, double val)
{
  std::ofstream fout(fname, std::ios::out | std::ios::app);
  fout << tag << "   " << val << std::endl;
  fout.close(); fout.clear();
}
void write_dat(const char* fname, const char* tag, float val)
{
  std::ofstream fout(fname, std::ios::out | std::ios::app);
  fout << tag << "   " << val << std::endl;
  fout.close(); fout.clear();
}
void write_dat(const char* fname, const char* tag, int val)
{
  std::ofstream fout(fname, std::ios::out | std::ios::app);
  fout << tag << "   " << val << std::endl;
  fout.close(); fout.clear();
}
