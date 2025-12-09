#ifndef _MYGLUTDISPLAY_H
#define _MYGLUTDISPLAY_H

#include <string.h>
#include <iostream>
#include<fstream>
#include<math.h>
#if defined(_WIN32) || defined(_WIN64)
#include <Windows.h>
#endif
#include <GL/gl.h>
#include <FL/Fl_Gl_Window.H>
#include <FL/glu.h>
#include <FL/glut.H>

#include "myheader.h"

extern float obj_pos[3];
extern float rotmat[16];
// 原子のヒット機能のオンオフを指定する変数です。2024/11 CJS
extern int enable_atom_click;
extern int mdspeed, radius, show_axis, show_only_elem, show_cell;
extern int segments, ibase, icnt, draw_replica, draw_bond, draw_bond_pbc, draw_force;
extern int bondthickness, bond_display_type, bond_display_color;
extern int draw_load, draw_aux;
extern int *iatom, *repidx, ievec, ievec_num, show_cnt_wall, show_cnt_wall_num, show_cnt_ring;
extern int show_cnt_wallv;
extern float scl, evec_len, vscl, vscl_force, cnt_ring_radius;
extern double mat[3][3];
extern GLuint objects;
extern int edit_elem_mode, select_atom[10], select_atom_repidx[10];
extern GLfloat red[],blue[],green[],yellow[],purple[],gray[],black[],white[];
extern GLfloat blue2[],green2[];
extern GLfloat **color, **color0;
extern int createconfig_mode, iatom_pick;
extern float replica_range;
extern int mode_cnt_corrugation, show_cnt_all, show_cnt_layer[6];
extern float bondlength;

// NOTE: myglutdisplay.cpp 外で実態が宣言されています。
void md();
void bookkeep();
void stretch(double x, double y, double z);
void writedata();
void mk_helmat(int iel, double rx[], double ry[], double rz[], double mat[][3]);
void out_of_cell();
void change_atom_color();
double dsquare(double x);
void calc_center_xy(double& x, double& y);


void myGlutDisplay(void);

#endif // _MYGLUTDISPLAY_H
