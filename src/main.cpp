#include <string.h>
#include <iostream>
#if defined(_WIN32) || defined(_WIN64)
#include <Windows.h>
#endif
#include <GL/gl.h>
#include <fstream>
#include <math.h>
//#include <unistd.h>
#if defined __linux__ || defined __APPLE__
#include <unistd.h>
#else
#include <direct.h>
#endif
#if !defined __linux__ && !defined __APPLE__
#define snprintf sprintf_s
#endif
//ReaxFF-related
#include <vector>
#include <iomanip>
#include <memory>
// 2024/02 CJS FLTK のインクルードを追加します。
#include <FL/Fl_File_Chooser.H>
#include <FL/Fl_File_Browser.H>
#include <FL/Fl.H>
// 2024/02 CJS ここまで
#include "Configuration.h"
#include "Input.h"
#include "ReaxFF.h"
#include "ReaxPotential.h"
#include "Population.h"
#include "Utils.h"
#include "DifferentialEvolution.h"
#include "PowellMinimization.h"
#include "PFMPI.h"
// added when modifying Chromosome.h (2022.01.15) for compile on Win
const double Chromosome::F_LOWER = 0.1;
const double Chromosome::F_UPPER = 0.9;
// 2024/02 CJS ウィンドウのインクルードを追加します。
#include "MainWindow.h"
// 2024/02 CJS ここまで
#include "myheader.h"

int istep = 0, istep0 = 0;
int mdmotion = 0;
int mdspeed = 10;
double dt = 1.0e-15;
int confwrint = 0, autocap = 0;
int itolfor = 0; float tolfor = 0.02;
int ireheat = 0;
float reheat_temp = 300.0, reheat_dtm = 2.0;
int reheat_step = 10, reheat_count = 5;
int itolstep = 0; int tolstep = 10000;
int itolstress = 0; float tolstress = 10.0;
int irecipe = 0;
int inogui = 0;
double msd = 0;
double rcut, rcut2;
float temp_set = 100.0;
float tempc = 0.0, cellx, celly, cellz, f_max, epotatom, dtm;
float cell1x, cell1y, cell1z;
float cell2x, cell2y, cell2z;
float cell3x, cell3y, cell3z;
float slice_1min = 0, slice_1max = 1, slice_2min = 0, slice_2max = 1, slice_3min = 0, slice_3max = 1;
int trim_mode = 0, trim_cylinder_axis = 1;
float trim_cylinder_diameter = 5.0, trim_cylinder_diameter2 = 5.0;
float rigid_shift_z = 0.5, rigid_shift_dx = 0, rigid_shift_dy = 0, rigid_shift_dz = 0;
int rigid_shift_xtimes = 10, rigid_shift_ytimes = 10, rigid_shift_relax_maxstep=100;
int ifix_atoms;
FILE* energyfile = fopen("energy.d", "w");
FILE* stressfile = fopen("stress.d", "w");
FILE* ssnorfile = fopen("ssnormal.d", "w");
FILE* cellfile = fopen("cell.d", "w");
FILE* ecnorfile = fopen("ecnormal.d", "w");
FILE* ecallfile = fopen("ecall.d", "w");
FILE* msdfile = fopen("msd.d", "w");
FILE* miscfile = fopen("misc.d", "w");
FILE* tempfile = fopen("temp.d", "w");
float strs_xx, strs_yy, strs_zz, strs_xy, strs_yz, strs_zx;
float strs_set_xx, strs_set_yy, strs_set_zz, strs_set_xy, strs_set_yz, strs_set_zx;
int cellfix_xx, cellfix_yy, cellfix_zz, cellfix_xy, cellfix_yz, cellfix_zx;
float rcut_f, frcmar = 0.2, frc_f;
float fp_alph_ini = 0.1, fp_ffinc = 1.1, fp_ffdec = 0.5, fp_ffalph = 0.99, fp_ffdtmax = 10.0;
int fp_nfmin = 5;
float auto_val1, auto_val2;
int group_num = 0, group_range_axis = 0, group_num_show = 0, if_group_num_show = 0;
float group_range_max = 1.0, group_range_min = 0.0;

// for Rotation/Motion by Mouse
int cx, cy; float cs, cpos0, cpos1; double sx, sy;
float cq[4] = { 1.0, 0.0, 0.0, 0.0 }; // Initial quaternion
float tq[4]; // quaternion while mouse drag
bool right_button;

Atom atom;
Cell cell;
Book book;
Tersoff tersoff;
GEAM geam;
Eammis eammis;
Adp adp;
Dipole dipole;
Cond cond;
Dmc2d dmc2d;
//Bre *bre = new Bre;
Bre bre;
Fire fire;
Cg cg;
Combined combined;
Pair* pairPot;
Sw sw; // Kubo 20140618 X0.1
//ReaxFF
MyReax myreax;
Shell shell;
Ew ew;

#ifdef CUDA
int* arrayDint;
int* arrayDrepatom;
double* arrayDrx, * arrayDry, * arrayDrz, * arrayDfx, * arrayDfy, * arrayDfz;
double* arrayDepot, * arrayDhmat, arrayHhmat[9];
int* arrayDalistnum, * arrayDalist;
int iarray[( NMAX + 1 ) * ( NBOOK + 1 ) * 4];
#endif

void qmul(float r[], const float p[], const float q[]);
void qua2rot(float r[], float q[]);
void rot2qua(float q[], float r[]);
void out_of_cell();
void md(), md_set();
//void vscale();
void bookkeep(); void bookkeep_write(); void bookkeep_read();
//void e_force();
void readconfig(const char* fname);
void readconfig(const char* fname, int mode);
void readqcelement(const char* fname);
void readqcelement();
void writedata();
void writedata_initialize();
void writeconfig(const char* fname); void writeconfig_abs(const char* fname);
void writeposcar(const char* fname);
void writecfgs(const char* fname); void writecfge(const char* fname);
void writeqcelement(const char* fname);
void readsetdat(const char* fname);
void writesetdat(const char* fname);
//void instability();
//void instability_dsy();
void instability_atomcell();
void instability_atom();
void instability_atom_noremove();
void instability_QC();
void mk_helmat(int iel, double rx[], double ry[], double rz[], double mat[][3]);
void createconfig();
void cnt_wall_set();
void cnt_wall_read(const char* fname); void cnt_wall_write(const char* fname);
void cnt_wall_discard();
void myGlutDisplay(void);
void bond_set();
void create_lattice(int lattice_type, float alat, float clat,
    float& alplat1,
    float& c1x, float& c1y, float& c1z,
    float& c2x, float& c2y, float& c2z,
    float& c3x, float& c3y, float& c3z);
void stretch_celladjust(float x, float y, float z);
void stretch_celladjust(float c1x, float c1y, float c1z,
    float c2x, float c2y, float c2z,
    float c3x, float c3y, float c3z);
void slice(float r1min, float r1max, float r2min, float r2max, float r3min, float r3max);
void rigid_shift(float, float, float, float);
void rigid_shift_xy(float, float, float, float, int, int);
int  rigid_shift_relax(float, int);
void rotate_cell(float x, float y, float z, int reverse);
void relax_fire_reset();
void relax_cg_reset();
void get_first_arg(std::string& line, std::string& arg1);
int count_arg_number(std::string line);
#ifndef NOPNG
void capture();
#endif
void phonon_calc(); void phonon_calc(int m);
void cntshell_calc_test(double dmu, double dnu, int n);
void cntshell_calc_1(); void cntshell_calc_2(); void cntshell_calc_all();
void neb_calc(int neb_num, const char* fini, const char* fend);
void stresscheck(double eps_strschk);
void stress_set();
void force_check(double fcheck_disp);
void calc_elastic_const();
void potential_set(int control);
void change_atom_color();
double csp(int i, int mi);
void potential();
void loading();
void set_potfile();
void multiply_cell(int ix, int iy, int iz);
void remove_atom(int ia); void move_atom(int ia, float arx, float ary, float arz);
void remove_atom(int ia, int ib);
void move_atom(int ia, const char* aasp, float arx, float ary, float arz);
void move_atom_shift(float arx, float ary, float arz);
void move_atom(int ia, float arx, float ary, float arz);
void move_atom_shift(int ia, float arx, float ary, float arz);
void move_atom_shift(int ia, int ib, float arx, float ary, float arz);
void rotate_atom(float x, float y, float z);
void rotate_atom(int ia, int ib, float x, float y, float z);
void add_polyhed(int polyhedtype, float polyhedsize);
void set_atom_color(); void set_atom_color(int ia);
void add_atom(const char* aasp, float arx, float ary, float arz);
void calc_distance(int select_atom[], int select_atom_repidx[]);
void calc_angle(int select_atom[], int select_atom_repidx[]);
void calc_dihedral(int select_atom[], int select_atom_repidx[]);
void view_save(); void view_read();
void autocalc();
void show_e_force_all();
void set_atom_group(int, int, float, float); void reset_atom_group();
void set_atom_group_visibility();
void trim(int, int, float, float); void trim(int, int, float);
void ( *integral_a[NINTEGRAL] )( ), ( *integral_b[NINTEGRAL] )( );
int integral_type;

float org_x = 0.0, org_y = 0.0, org_z = 0.0;
int size_w, size_h;
float xy_aspect;
int   last_x, last_y;
float rotationX = 0.0, rotationY = 0.0;
float rotmat[16] = {
  1,0,0,0,
  0,1,0,0,
  0,0,1,0,
  0,0,0,1
};
float obj_pos[] = { 0.0, 0.0, 0.0 };
float scl = 1.0;


float vscl = 0.2, vscl_force = 1.0;


GLfloat red[] = { 1.0, 0.0, 0.0, 1.0 };
GLfloat yellow[] = { 1.0, 1.0, 0.0, 1.0 };
GLfloat white[] = { 2.0, 2.0, 2.0, 1.0 };
GLfloat gray[] = { 0.7, 0.7, 0.5, 1.0 };
GLfloat blue[] = { 0.1, 0.1, 0.9 };
GLfloat blue2[] = { 0.0, 0.0, 1.0 };
GLfloat green[] = { 0.0, 1.0, 0.0 };
GLfloat green2[] = { 0.0, 0.7, 0.0 };
GLfloat black[] = { 0.1, 0.1, 0.1, 1.0 };
GLfloat purple[] = { 1.0, 0.0, 1.0, 1.0 };

//GLfloat *color[NOBJECTS];
GLfloat** color, ** color0;
//int iatom[NREPLICA]; int repidx[NREPLICA]; int ibase, icnt;
int* iatom;
int* repidx;
int ibase;
int icnt;


/** These are the live variables passed into GLUI ***/
int   ortho_state = 1;
// 原子のヒット機能のオンオフを指定する変数です。0: ヒット機能オフ、0 以外: ヒット機能オン 2024/11 CJS
int   enable_atom_click = 1;
int   draw_bond = 0, draw_bond_pbc = 0, draw_force = 0, draw_load = 0, draw_aux;
int   bond_display_type = 1; // 0:line, 1:cylinder
int   bond_display_color = 0;
int   ensemble = 0; int ensemble_p = ensemble;
int   relax_algo = 0, relax_accel = 1, relax_accel_interval = 10;
float relax_accel_threshold = 0.2, relax_damper_value = 0.0;
int   cell_relax_rep = 1; float cell_relax_tolfor = 0.01;
int   config_type = 0;
int   lattice_type = 0;
char  config_atom[3] = "Al";
char  config_atom2[3] = "Al";
int   irepx = 1, irepy = 1, irepz = 1;
int   icntm = 8, icntn = 0;
float cscnt = 50.0, rotz = 0.0, shiftz = 0.0;
float alat = 4.0, clat = 4.0, alplat1 = 90;
int   imerge = 0;
char  current_config_name[100], current_qcelement_name[100];
char  cwdname[80] = "aaa";
//const char *potstring_list[] = {"Morse", "Pair", "GEAM", "Tersoff Si(B)", "Tersoff Si(C)", "Tersoff Si(B*)", "Tersoff C_Si_Ge", "Tersoff B_N_C", "Brenner", "EAM Mishin", "Dipole", "ADP", "AIREBO", "NiYSZ", "Vashishta", "AIREBO_BN", "SW GaN" };
const char* ensemble_list[] = { "NVE", "NVT (scaling)", "NVT (thermostat)", "Relaxation (atom)", "NPH", "NPT", "NPH + PRdamper", "NPT + PRdamper", "Full relaxation (PRdampaer)", "Full relaxation (static)", "Diffusion MC" };
const char* kpstring_list[] = { "fcc/dia", "bcc", "sc", "hcp", "chain" };
const char* color_list[] = { "Red", "Blue", "Green", "Yellow", "Purple", "Gray", "Black", "White" };
const char* sw_arg_list[] = { "GaN", "SiN" };
const char* tersoff_arg_list[] = { "B", "C" };
int   ipottype = 0; int ipotarg = 0;
int   notrans = 0;
int   incell = 0;
int   segments = 15;
int   radius = 8;
int   show_only_elem = 0;
int   show_axis = 0;
int   show_cell = 1;
int   show_cnt_wall = 0, show_cnt_wallv = 0, show_cnt_ring = 0;
int   show_cnt_wall_num = 1;
int   mode_cnt_corrugation = 0, read_cntwall = 0;
int   show_cnt_all = 1;
int   show_cnt_layer[6];
int   cnt_load_algo = 0;
float dexdt = 0.0, ex = 0.0;
float deydt = 0.0, ey = 0.0;
float dezdt = 0.0, ez = 0.0;
int repeat_lz; float repeat_lz_min, repeat_lz_max;
int ifgrab1 = 0, ifgrab2 = 0;
int grab1num = 0, grab2num = 0;
float grabdxdt1 = 0.0, grabdydt1 = 0.0, grabdzdt1 = 0.0;
float grabdxdt2 = 0.0, grabdydt2 = 0.0, grabdzdt2 = 0.0;
int ifpush1 = 0, ifpush2 = 0;
int push1num = 0, push2num = 0;
float pushfx1 = 0.0, pushfy1 = 0.0, pushfz1 = 0.0;
float pushfx2 = 0.0, pushfy2 = 0.0, pushfz2 = 0.0;
int   yzplane_punch = 0;
float yzplane_punch_d = 10.0, yzplane_punch_dd = 0.0, yzplane_punch_ftot = 0.0;
float cnt_pressure = 0.1, cnt_pressure_ftot = 0, cnt_pressure_gpa, cnt_pressure_gpa2;
float outermost_radius_f, cylinder_side_area_f, cylinder_side_area0_f;
float cnt_ring_radius = 15.0, cnt_ring_fmax = 10.0, cnt_ring_sharpness = 2.0;
float cntshell_dmu = 0.0, cntshell_dnu = 0.0, cntshell_eps = 0.0001;
int cntshell_n = 2;
int   ievec = 0;
int   ievec_num = 1;
float evec_len = 10.0, eigval = 0;
double mat[3][3];
int b_state = 1;
float xx1, yy1, zz1, xx2, yy2, zz2, xx3, yy3, zz3, xx4, yy4, zz4;
GLuint objects;
int edit_elem_mode = 0, select_atom[10], select_atom_repidx[10];
int measure_mode = 0;
int draw_replica = 0;
float replica_range = 0.0;
int hoge = 0;
float bondlength = 1.52; // in ang
int bondthickness = 5;
int capture_count = 0;
float eps_strschk = 0.000001;
float fcheck_disp = 0.0001;
int color_mode = 0, color_mode_auto = 0;
float color_mode_vmin, color_mode_vmax;
float prmass_scale = 0, nhmass_scale = 0;
bool iprdamper = false, ivscale = false, irelax = false;
float prdamper_val1 = 1, prdamper_val2 = 0, prlimit = 0;
int hessian_read = 0, hessian_write = 0, inst_mode;
int phonon_rep = 2, phonon_kp = 0, phonon_knum = 50;
float dos_gauss_width = 0.05; int dos_kmesh = 10;
int neb_num = 10, neb_node = 0, neb_ite = 100, neb_init_read = 0;
float neb_tol_fac = -11;
int neb_kconst = 0;
float neb_kfactor = 1.0;
int iatom_pick = 0, createconfig_mode = 0, editatom_mode = 0;
int iatom_from = 1, iatom_to = 1;
float atomrx, atomry, atomrz;
int milX1 = 1, milX2 = 0, milX3 = 0;
int milY1 = 0, milY2 = 1, milY3 = 0;
int milZ1 = 0, milZ2 = 0, milZ3 = 1;
float cellrot_x = 0, cellrot_y = 0, cellrot_z = 0;
int marked_atom[3] = { 0,0,0 }, marked_atom_color[3] = { 0,0,0 };
int initial_velocity = 0;
float rotatomx = 0, rotatomy = 0, rotatomz = 0;
int polyhedtype = 0; float polyhedsize = 1;

int capture_width = 0;
int capture_height = 0;


// Using a std::string as a live variable is safe.
std::string text = "Test string..";
std::string potfile;


// 2024/02 CJS pointer_cb の引数にかかわる処理を FLTK に合わせて変更しました。
void pointer_cb(Fl_Widget* widget, long process_id)
{
    int id = process_id;

    // id が -1 の場合は何もしません。
    if (id == -1)
    {
        return;
    }
    else if (id == CB_QCELM_FB)
    {
        Fl_File_Chooser* chooser = (Fl_File_Chooser*)( widget );

#if defined _WIN32 || defined _WIN64
        // Windows 環境の場合、文字コードを変換します。
        std::string original_text = chooser->value();
        std::string text = utf8_to_shift_jis(original_text);
#else
        std::string text = chooser->value();
#endif
        char fname[300] = "aaa";
        strcpy(fname, text.c_str());
        readqcelement(fname);
    }
    else if (id == CB_CNTWALL_FB)
    {
        Fl_File_Chooser* chooser = (Fl_File_Chooser*)( widget );

#if defined _WIN32 || defined _WIN64
        // Windows 環境の場合、文字コードを変換します。
        std::string original_text = chooser->value();
        std::string text = utf8_to_shift_jis(original_text);
#else
        std::string text = chooser->value();
#endif
        char fname[300] = "aaa";
        strcpy(fname, text.c_str());
        cnt_wall_read(fname);
    }
    else if (id == CB_SETDAT_FB)
    {
        Fl_File_Chooser* chooser = (Fl_File_Chooser*)( widget );
#if defined _WIN32 || defined _WIN64
        // Windows 環境の場合、文字コードを変換します。
        std::string original_text = chooser->value();
        std::string text = utf8_to_shift_jis(original_text);
#else
        std::string text = chooser->value();
#endif
        char fname[300] = "aaa";
        strcpy(fname, text.c_str());
        readsetdat(fname);
    }
    else if (id == WRITE_SETDAT_ID)
    {
        chdir(cwdname);
        writesetdat("SETDAT.OUT");
    }
    else if (id == RESET_ID)
    {
        //    org_x=0.0; org_y=0.0; org_z=0.0;
        //    rotationX=0.0; rotationY=0.0;
        //    for (int i=0; i<=15; i++) { rotmat[i]=0; }
        //    rotmat[0]=1;rotmat[5]=1;rotmat[10]=1;rotmat[15]=1;
        relax_fire_reset(); relax_cg_reset(); md_set();
        for (int i = 1; i <= atom.natom; i++)
        {
            atom.rx[i] = atom.rx_org[i];
            atom.ry[i] = atom.ry_org[i];
            atom.rz[i] = atom.rz_org[i];
            atom.vx[i] = 0.0; atom.vy[i] = 0.0; atom.vz[i] = 0.0;
            atom.fx[i] = 0.0; atom.fy[i] = 0.0; atom.fz[i] = 0.0;
            atom.ax[i] = 0.0; atom.ay[i] = 0.0; atom.az[i] = 0.0;
            atom.bx[i] = 0.0; atom.by[i] = 0.0; atom.bz[i] = 0.0;
            atom.cx[i] = 0.0; atom.cy[i] = 0.0; atom.cz[i] = 0.0;
        }
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                cell.hmat[i][j] = cell.hmat_org[i][j];
            }
        }
        //cellx=cell.hmat[0][0]/ang;celly=cell.hmat[1][1]/ang;cellz=cell.hmat[2][2]/ang;
        cellx = sqrt(cell.hmat[0][0] * cell.hmat[0][0] + cell.hmat[1][0] * cell.hmat[1][0] + cell.hmat[2][0] * cell.hmat[2][0]) / ang;
        celly = sqrt(cell.hmat[0][1] * cell.hmat[0][1] + cell.hmat[1][1] * cell.hmat[1][1] + cell.hmat[2][1] * cell.hmat[2][1]) / ang;
        cellz = sqrt(cell.hmat[0][2] * cell.hmat[0][2] + cell.hmat[1][2] * cell.hmat[1][2] + cell.hmat[2][2] * cell.hmat[2][2]) / ang;
        cell1x = cell.hmat[0][0] / ang; cell1y = cell.hmat[1][0] / ang; cell1z = cell.hmat[2][0] / ang;
        cell2x = cell.hmat[0][1] / ang; cell2y = cell.hmat[1][1] / ang; cell2z = cell.hmat[2][1] / ang;
        cell3x = cell.hmat[0][2] / ang; cell3y = cell.hmat[1][2] / ang; cell3z = cell.hmat[2][2] / ang;
        cell.Reset(); cell.volume = cell.Getvolume();
        // give 'dummy' value for istep and then reset it (to sync display)
        istep = -1;
        istep = 0;
        ex = 0.0; ey = 0.0; ez = 0.0; //dexdt = 0.0; deydt = 0.0; dezdt = 0.0;
        bond_set();
        writedata_initialize();
    }
    else if (id == INST_ID)
    {
        if (atom.QC)
        {
            instability_QC();
        }
        else
        {
            //if (inst_mode == 0) {
            //instability();
            //} else if (inst_mode == 1) {
            //instability_dsy();
            //} else if (inst_mode == 2) {
            if (inst_mode == 0)
            {
                instability_atomcell();
            }
            else if (inst_mode == 1)
            {
                instability_atom();
            }
            else if (inst_mode == 2)
            {
                instability_atom_noremove();
            }
        }
    }
    else if (id == WRITECONFIG_ID)
    {
        chdir(cwdname);
        writeconfig("CONFIG.OUT"); writeconfig_abs("CONFIG.OUT.ABS");
        writeposcar("POSCAR");
        writecfgs("config.cfg"); writecfge("config.cfge");
    }
    else if (id == WRITEQCELEMENT_ID)
    {
        chdir(cwdname);
        writeqcelement("QCELEMENT.OUT");
    }
    else if (id == MDSWITCH_ID)
    {
        if (mdmotion == 0)
        {
            mdmotion = 1;
            istep0 = istep;
        }
        else
        {
            mdmotion = 0;
        }
    }
    else if (id == STRSCHK_ID)
    {
        stresscheck((double)eps_strschk);
    }
    else if (id == PHONON_CALC_ID)
    {
        phonon_calc();
    }
    else if (id == DOS_CALC_ID)
    {
        phonon_calc(1);
    }
    else if (id == NEB_CALC_ID)
    {
        neb_calc(neb_num, "CONFIG.NEB.INI", "CONFIG.NEB.END");
    }
    else if (id == CNTSHELL_TEST_ID)
    {
        cntshell_calc_test((double)cntshell_dmu, (double)cntshell_dnu, cntshell_n);
    }
    else if (id == CNTSHELL_CALC_ID)
    {
        cntshell_calc_1();
        //cntshell_calc_2();
    }
    else if (id == EXTRA_AUTO_ID)
    {
        autocalc();
    }
    else if (id == EXTRA_BKWRITE_ID)
    {
        bookkeep_write();
    }
    else if (id == EXTRA_BKREAD_ID)
    {
        bookkeep_read();
    }
    else if (id == MEASURE_CLEAR_ID)
    {
        for (int i = 0; i < 3; i++) { select_atom[i] = 0; select_atom_repidx[i] = 0; }
        set_atom_color();
    }
    else if (id == MEASURE_CALC2_ID)
    {
        calc_distance(select_atom, select_atom_repidx);
    }
    else if (id == MEASURE_CALC3_ID)
    {
        calc_angle(select_atom, select_atom_repidx);
    }
    // 2024/02 CJS 使用されていない処理のため、コメント アウトしました。
    // else if (id == MEASURE_CALC4_ID)
    // {
    //     calc_dihedral(select_atom, select_atom_repidx);
    // }
    else if (id == EDIT_ELEM_TAKE_ID)
    {
        if (( select_atom[0] != 0 ) && ( select_atom[1] != 0 ) && ( select_atom[2] != 0 ) && ( select_atom[3] != 0 ))
        {
            printf("New element: %d %d %d %d  %d %d %d %d\n",
                select_atom[0], select_atom[1], select_atom[2], select_atom[3],
                select_atom_repidx[0], select_atom_repidx[1], select_atom_repidx[2], select_atom_repidx[3]);
            int i = atom.nelem + 1;
            atom.elem_v = (int**)realloc(atom.elem_v, sizeof(int*) * ( i + 1 ));
            atom.elem_v_rep = (int**)realloc(atom.elem_v_rep, sizeof(int*) * ( i + 1 ));
            atom.elem_v[i] = (int*)malloc(sizeof(int) * 5);
            atom.elem_v_rep[i] = (int*)malloc(sizeof(int) * 5);
            atom.elem_v[i][1] = select_atom[0]; atom.elem_v_rep[i][1] = select_atom_repidx[0];
            atom.elem_v[i][2] = select_atom[1]; atom.elem_v_rep[i][2] = select_atom_repidx[1];
            atom.elem_v[i][3] = select_atom[2]; atom.elem_v_rep[i][3] = select_atom_repidx[2];
            atom.elem_v[i][4] = select_atom[3]; atom.elem_v_rep[i][4] = select_atom_repidx[3];
            atom.nelem = i;
            readqcelement();
            for (int ii = 0; ii < 10; ii++) { select_atom[ii] = 0; select_atom_repidx[ii] = 0; }
            for (int ii = 1; ii <= atom.natom + icnt; ii++) { memcpy(color[ii], yellow, sizeof(GLfloat) * 4); }
        }
    }
    else if (id == EDIT_ELEM_XZ_TAKE_ID)
    {
        if (( select_atom[0] != 0 ) && ( select_atom[1] != 0 ) && ( select_atom[2] != 0 ))
        {
            //      printf("New element: %d %d %d %d  %d %d %d %d\n",
            //	     select_atom[0],select_atom[1],select_atom[2],select_atom[3],
            //	     select_atom_repidx[0],select_atom_repidx[1],select_atom_repidx[2],select_atom_repidx[3]);
            int i = atom.nelem + 3;
            atom.elem_v = (int**)realloc(atom.elem_v, sizeof(int*) * ( i + 1 ));
            atom.elem_v_rep = (int**)realloc(atom.elem_v_rep, sizeof(int*) * ( i + 1 ));
            for (int ii = atom.nelem + 1; ii <= atom.nelem + 3; ii++)
            {
                atom.elem_v[ii] = (int*)malloc(sizeof(int) * 5);
                atom.elem_v_rep[ii] = (int*)malloc(sizeof(int) * 5);
            }
            i = atom.nelem + 1;
            atom.elem_v[i][1] = select_atom[0]; atom.elem_v_rep[i][1] = select_atom_repidx[0];
            atom.elem_v[i][2] = select_atom[1]; atom.elem_v_rep[i][2] = select_atom_repidx[1];
            atom.elem_v[i][3] = select_atom[2]; atom.elem_v_rep[i][3] = select_atom_repidx[2];
            atom.elem_v[i][4] = select_atom[2]; atom.elem_v_rep[i][4] = select_atom_repidx[2] + 2;
            i = atom.nelem + 2;
            atom.elem_v[i][1] = select_atom[0]; atom.elem_v_rep[i][1] = select_atom_repidx[0];
            atom.elem_v[i][2] = select_atom[1]; atom.elem_v_rep[i][2] = select_atom_repidx[1];
            atom.elem_v[i][3] = select_atom[1]; atom.elem_v_rep[i][3] = select_atom_repidx[1] + 2;
            atom.elem_v[i][4] = select_atom[2]; atom.elem_v_rep[i][4] = select_atom_repidx[2] + 2;
            i = atom.nelem + 3;
            atom.elem_v[i][1] = select_atom[0]; atom.elem_v_rep[i][1] = select_atom_repidx[0];
            atom.elem_v[i][2] = select_atom[0]; atom.elem_v_rep[i][2] = select_atom_repidx[0] + 2;
            atom.elem_v[i][3] = select_atom[1]; atom.elem_v_rep[i][3] = select_atom_repidx[1] + 2;
            atom.elem_v[i][4] = select_atom[2]; atom.elem_v_rep[i][4] = select_atom_repidx[2] + 2;
            atom.nelem = i;
            readqcelement();
            for (int ii = 0; ii < 10; ii++) { select_atom[ii] = 0; select_atom_repidx[ii] = 0; }
            for (int ii = 1; ii <= atom.natom + icnt; ii++) { memcpy(color[ii], yellow, sizeof(GLfloat) * 4); }
        }
    }
    else if (id == CALC_ID)
    {
        bookkeep(); potential(); loading(); change_atom_color();
        // For output in "Status" area
        //cellx=cell.hmat[0][0]/ang;celly=cell.hmat[1][1]/ang;cellz=cell.hmat[2][2]/ang;
        cellx = sqrt(cell.hmat[0][0] * cell.hmat[0][0] + cell.hmat[1][0] * cell.hmat[1][0] + cell.hmat[2][0] * cell.hmat[2][0]) / ang;
        celly = sqrt(cell.hmat[0][1] * cell.hmat[0][1] + cell.hmat[1][1] * cell.hmat[1][1] + cell.hmat[2][1] * cell.hmat[2][1]) / ang;
        cellz = sqrt(cell.hmat[0][2] * cell.hmat[0][2] + cell.hmat[1][2] * cell.hmat[1][2] + cell.hmat[2][2] * cell.hmat[2][2]) / ang;
        cell1x = cell.hmat[0][0] / ang; cell1y = cell.hmat[1][0] / ang; cell1z = cell.hmat[2][0] / ang;
        cell2x = cell.hmat[0][1] / ang; cell2y = cell.hmat[1][1] / ang; cell2z = cell.hmat[2][1] / ang;
        cell3x = cell.hmat[0][2] / ang; cell3y = cell.hmat[1][2] / ang; cell3z = cell.hmat[2][2] / ang;
        f_max = atom.Fmax() / eV * ang;
        epotatom = atom.epotsum / eV / atom.natom;
        show_e_force_all();
    }
    // 2024/02 CJS 使用されていない処理のため、コメント アウトしました。
    // else if (id == CELLSIZE_APPLY_ID)
    // {
    //     stretch_celladjust(cellx, celly, cellz);
    // }
    // メイン ウィンドウ -> Status パネルの Cell(A) ~ の入力があった際に呼び出される処理です。
    else if (( id == STATUS_LX ) || ( id == STATUS_LY ) || ( id == STATUS_LZ ))
    {
        stretch_celladjust(cellx, celly, cellz);
    }
    // メイン ウィンドウ -> Status パネルの dt (fs) の入力があった際に呼び出される処理です。
    else if (id == STATUS_DT)
    {
        dt = (double)dtm * 1e-15;
    }
    else if (id == STRSSET_ID)
    {
        stress_set();
    }
    else if (id == REMOVEATOM_DO_ID)
    {
        remove_atom(iatom_pick);
    }
    else if (id == REMOVEPART_DO_ID)
    {
        remove_atom(iatom_from, iatom_to);
    }
    else if (id == ADDATOM_DO_ID)
    {
        add_atom(config_atom, atomrx, atomry, atomrz);
        iatom_pick = atom.natom;
        atomrx = atom.rx[iatom_pick] / ang; atomry = atom.ry[iatom_pick] / ang; atomrz = atom.rz[iatom_pick] / ang;
    }
    else if (id == MOVEATOM_DO_ID)
    {
        //move_atom(iatom_pick, atomrx, atomry, atomrz);
        move_atom(iatom_pick, config_atom, atomrx, atomry, atomrz);
    }
    else if (id == SHIFTALL_DO_ID)
    {
        move_atom_shift(atomrx, atomry, atomrz);
    }
    else if (id == SHIFTONE_DO_ID)
    {
        move_atom_shift(iatom_pick, atomrx, atomry, atomrz);
    }
    else if (id == SHIFTPART_DO_ID)
    {
        move_atom_shift(iatom_from, iatom_to, atomrx, atomry, atomrz);
    }
    else if (id == ROTATEALL_DO_ID)
    {
        rotate_atom(rotatomx, rotatomy, rotatomz);
    }
    else if (id == ROTATEPART_DO_ID)
    {
        rotate_atom(iatom_from, iatom_to, rotatomx, rotatomy, rotatomz);
    }
    else if (id == PACKALL_DO_ID)
    {
        out_of_cell();
    }
    else if (id == ADDPOLYHED_DO_ID)
    {
        add_polyhed(polyhedtype, polyhedsize);
    }
    else if (id == CREATECONFIG_DO_ID)
    {
        createconfig();
    }
    else if (id == MULTIPLYCELL_DO_ID)
    {
        multiply_cell(irepx, irepy, irepz);
    }
    else if (id == CNTWALL_DO_ID)
    {
        cnt_wall_set();
    }
    else if (id == CNTWALL_WRITE_ID)
    {
        cnt_wall_write("CNTWALL.SAVE");
    }
    else if (id == SETGROUP_DO_ID)
    {
        set_atom_group(group_num, group_range_axis, group_range_max, group_range_min);
        set_atom_group_visibility();
    }
    else if (id == RESETGROUP_DO_ID)
    {
        reset_atom_group();
        set_atom_group_visibility();
    }
}

// 2024/02 CJS control_cb の引数にかかわる処理を FLTK に合わせて変更しました。
void control_cb(Fl_Widget* widget, long process_id)
{
    int id = process_id;

    // id が -1 の場合は何もしません。
    if (id == -1)
    {
        return;
    }
    if (show_cnt_wall_num > atom.nwall)
    {
        show_cnt_wall_num = atom.nwall;
    }
    if (id == CB_INST)
    {
        if (ievec_num < 1) { ievec_num = 1; }
        if (ievec_num > MAXMODE) { ievec_num = MAXMODE; }
        eigval = atom.eigval[ievec_num - 1];
    }
    if (id == CB_PHONON)
    {
        //if (phonon_rep > 4) { phonon_rep = 4; }
        if (phonon_rep < 1) { phonon_rep = 1; }
    }
    if (id == CB_NEB)
    {
        if (neb_num < 3) { neb_num = 3; }
        if (neb_num > 999) { neb_num = 999; }
    }
    if (id == CB_NEB_SHOW)
    {
        if (neb_node >= neb_num) { neb_node = neb_num - 1; }
        if (neb_node < 0) { neb_node = 0; }
        char filepath[80] = "CONFIG.NEB."; char numc[10];
        if (neb_node == neb_num - 1)
        {
            strcpy(numc, "END");
        }
        else if (neb_node == 0)
        {
            strcpy(numc, "INI");
        }
        else { snprintf(numc, sizeof(numc), "%02d", neb_node); }

        // スピナーの値を更新します。
        Fl_Spinner* spinner = (Fl_Spinner*)( widget );
        spinner->value(neb_node);

        strcat(filepath, numc);
        readconfig(filepath, 1);
        bookkeep(); potential(); loading(); change_atom_color(); //added 20210718
    }
    if (id == CB_EDITATOM)
    {
        if (iatom_pick > atom.natom) { iatom_pick = atom.natom; }
        if (iatom_pick < 0) { iatom_pick = 0; }
        if (( iatom_pick > 0 ) && ( iatom_pick <= atom.natom ))
        {
            strcpy(config_atom, atom.asp[iatom_pick]);
            atomrx = atom.rx[iatom_pick] / ang; atomry = atom.ry[iatom_pick] / ang; atomrz = atom.rz[iatom_pick] / ang;
        }
        if (iatom_from <= 0) { iatom_from = 1; }
        if (iatom_from > atom.natom) { iatom_from = atom.natom; }
        if (iatom_to <= 0) { iatom_to = 1; }
        if (iatom_to > atom.natom) { iatom_to = atom.natom; }
    }
    if (id == CB_BOND) { bond_set(); }
    if (id == CB_COLOR_MODE) { change_atom_color(); }
    if (id == CB_MARKED_ATOM)
    {
        for (int i = 0; i < 3; i++)
        {
            if (marked_atom[i] < 0) { marked_atom[i] = 0; }
            else if (marked_atom[i] > atom.natom) { marked_atom[i] = atom.natom; }
        }
        change_atom_color();
    }
    if (id == CB_ENSEMBLE) { md_set(); }//if (ensemble == 1) { velocity_random(temp_set); } }
    if (id == CB_CELLFIX_ID)
    {
        if (cellfix_xx) { cell.fix[0][0] = true; }
        else { cell.fix[0][0] = false; }
        if (cellfix_yy) { cell.fix[1][1] = true; }
        else { cell.fix[1][1] = false; }
        if (cellfix_zz) { cell.fix[2][2] = true; }
        else { cell.fix[2][2] = false; }
        if (cellfix_xy) { cell.fix[0][1] = true; cell.fix[1][0] = true; }
        else { cell.fix[0][1] = false; cell.fix[1][0] = false; }
        if (cellfix_yz) { cell.fix[1][2] = true; cell.fix[2][1] = true; }
        else { cell.fix[1][2] = false; cell.fix[2][1] = false; }
        if (cellfix_zx) { cell.fix[2][0] = true; cell.fix[0][2] = true; }
        else { cell.fix[2][0] = false; cell.fix[0][2] = false; }
    }
    if (id == CB_POTFILE_FB)
    {
        // ウィジェットをファイル ブラウザにキャストします。
        Fl_File_Browser* file_browser = (Fl_File_Browser*)( widget );

        potfile = cwdname;
#if defined __linux__ || defined __APPLE__
        potfile += "/pot/";
#else
        potfile += "\\pot\\";
#endif
        // ファイル ブラウザの選択されたファイル名を取得します。
        potfile += file_browser->text(file_browser->value());
        std::string argtext = file_browser->text(file_browser->value());
        strcpy(atom.potential_arg, argtext.c_str());

        set_potfile();
    }
    if (id == CB_ROTATE) { rotate_cell(cellrot_x, cellrot_y, cellrot_z, 0); }
    if (id == CB_ROTATE_REV) { rotate_cell(cellrot_x, cellrot_y, cellrot_z, 1); }
    if (id == CB_RESET_VIEW)
    {
        org_x = 0.0; org_y = 0.0; org_z = 0.0;
        rotationX = 0.0; rotationY = 0.0;
        obj_pos[0] = 0; obj_pos[1] = 0; obj_pos[2] = 0; scl = 1.0;
        for (int i = 0; i <= 15; i++) { rotmat[i] = 0; }
        rotmat[0] = 1; rotmat[6] = -1; rotmat[9] = 1; rotmat[15] = 1;
    }
    if (id == CB_CELL_DIM)
    {
        stretch_celladjust(cell1x, cell1y, cell1z, cell2x, cell2y, cell2z, cell3x, cell3y, cell3z);
    }
    if (id == CB_LATTICE)
    {
        create_lattice(lattice_type, alat, clat, alplat1,
            cell1x, cell1y, cell1z, cell2x, cell2y, cell2z, cell3x, cell3y, cell3z);
    }
    if (id == CB_SLICE)
    {
        slice(slice_1min, slice_1max, slice_2min, slice_2max, slice_3min, slice_3max);
    }
    if (id == CB_TRIM)
    {
        if (trim_mode == 0) {
            trim(trim_mode, trim_cylinder_axis, trim_cylinder_diameter);
        } else if (trim_mode == 1) {
            trim(trim_mode, trim_cylinder_axis, trim_cylinder_diameter, trim_cylinder_diameter2);
        } else if (trim_mode == 2) {
            trim(trim_mode, trim_cylinder_axis, trim_cylinder_diameter);
        } else { printf("trim mode %d not implemented!!!\n",trim_mode); }
    }
    if (id == CB_RIGID_SHIFT)
    {
        rigid_shift(rigid_shift_z, rigid_shift_dx, rigid_shift_dy, rigid_shift_dz);
    }
    if (id == CB_RIGID_SHIFT_XY)
    {
        printf("Rigid shift output in gsf.d\n");
        rigid_shift_xy(rigid_shift_z, rigid_shift_dx, rigid_shift_dy, rigid_shift_dz,
            rigid_shift_xtimes, rigid_shift_ytimes);
    }
    if (id == CB_RIGID_SHIFT_RELAX) {
        int irigrlx = rigid_shift_relax(rigid_shift_z, rigid_shift_relax_maxstep);
        printf("Relaxation along z (above sep-z) took %d steps.\n",irigrlx); 
    }
    if (( id == CB_FIRE ) || ( id == CB_REL ))
    {
        relax_fire_reset(); relax_cg_reset();
    }
    if (id == CB_ELASC)
    {
        calc_elastic_const();
    }
    if (id == CB_FCHECK)
    {
        force_check((double)fcheck_disp);
    }
    if (id == CB_SHOW_GROUP)
    {
        set_atom_group_visibility();
    }
    if (id == CB_QUIT)
    {
        exit(0);
    }
}

// 2024/02 CJS idle 関数の実装を FLTK に合わせて変更しました。
void idle(void*)
{
    if (mdmotion == 1)
    {
        for (int i = 1; i <= mdspeed; i++)
        {
            if (istep % book.nbk == 0)
            {
                if (incell) { out_of_cell(); }
                bookkeep();
            }
            md(); if (mdmotion == 0) { break; }
            writedata();
            istep++;
        }
    }
}


/**************************************** main() ********************/
int main(int argc, char* argv[])
{
    getcwd(cwdname, 80);

    // ReaxFF initialization (tantatively written here)
    PFMPI::Init(&argc, &argv);
    myreax.myreaxinput.ReadFile("pot/ReaxFF.input");
    myreax.myreaxpot.ReadFile("pot/ReaxFF.pot");
    myreax.myreaxconf.ReadFile("pot/ReaxFF.config");
    ReaxFF::Initialize(myreax.myreaxpot, myreax.myreaxconf, myreax.myreaxinput);
    printf("## ReaxFF initialized ##\n");

#if defined __linux__ || defined __APPLE__
    if (argc > 1)
    {
        printf("Change directory to %s\n", argv[1]); chdir(argv[1]);
    }
    // default initial configuration
    strcpy(atom.potential_func, "GEAM"); book.algo = 1;
    if (argc <= 2)
    {
        readsetdat("SETDAT"); readconfig("CONFIG"); //readqcelement("QCELEMENT");
    }
    else if (argc == 3)
    {
        readsetdat("SETDAT"); readconfig(argv[2]); //readqcelement("QCELEMENT");
    }
    else if (argc == 4)
    {
        readsetdat(argv[3]); readconfig(argv[2]); //readqcelement("QCELEMENT");
    }
#else
    readsetdat("SETDAT"); readconfig("CONFIG");
#endif
    /*
    //Reax
    PFMPI::Init(&argc, &argv);
    input.ReadFile("pot/ReaxFF.input");
    reaxpot.ReadFile("pot/ReaxFF.pot");
    reaxconf.ReadFile("pot/ReaxFF.config");
    ReaxFF::Initialize(reaxpot,reaxconf,input);
    printf("Number of Atoms (ReaxFF): %d\n",ReaxFF::GetNumberOfAtoms());
    ReaxFF::CalcAll();
    double xxx,yyy,zzz;
    ReaxFF::GetForce(1,xxx,yyy,zzz);printf("%f %f %f\n",xxx,yyy,zzz);
    printf("ReaxFF OK\n");
    //Reax-end
    */
    if (( mode_cnt_corrugation ) && ( read_cntwall ) && ( cnt_load_algo == 0 ))
    {
        cnt_wall_read("CNTWALL");
    }
    //getcwd(cwdname,80); // Moved to above Oct2015
    printf("Current Working Directory = %s\n", cwdname);
    writedata_initialize(); md_set();

    if (!initial_velocity)
    {
        for (int i = 1; i <= atom.natom; i++)
        {
            atom.vx[i] = 0.0e3;   atom.vy[i] = 0.0e3; atom.vz[i] = 0.0e3;
        }
    }
    for (int i = 1; i <= atom.natom; i++)
    {
        atom.rx_org[i] = atom.rx[i]; atom.ry_org[i] = atom.ry[i]; atom.rz_org[i] = atom.rz[i];
    }

    // if SETDAT000 has "nogui yes"
    if (inogui > 0)
    {
        while (mdmotion == 1) { md(); istep++; }
        exit(0);
    }

    qua2rot(rotmat, cq);

    // 2024/02 CJS 追加
    // コントロール ウィンドウを作成します。
    std::unique_ptr<MainWindow> main_window = std::make_unique<MainWindow>();
    main_window->show();

    // アイドル時に呼び出される関数を登録します。
    Fl::add_idle(idle);

    // イベント ループを開始します。
    return Fl::run();

    // 2024/02 CJS 追加 ここまで
}
