#include <string>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#if defined(_WIN32) || defined(_WIN64)
#include <Windows.h>
#endif
#include <GL/gl.h>
#include <algorithm>
#include "myheader.h"

extern GLfloat** color, ** color0;
extern int color_mode, color_mode_auto;
extern float color_mode_vmin, color_mode_vmax;
extern int marked_atom[], marked_atom_color[];
extern GLfloat red[], blue[], green[], yellow[], purple[], gray[], black[], white[];

void hsv2rgb(float h, float s, float v, float& r, float& g, float& b);
void set_atom_color();
double csp(int i, int mi);
void sort(int n, double* val, int* nod);
double dsquare(double d);

void change_atom_color()
{
	if (color_mode == 0) { // color by species
		set_atom_color();
	}
	else if (color_mode >= 1) {
		// 1: color by energy, 2:CSP (Central Symmetric Parameter)
		// 3: Mises stress
		GLfloat hsv[3], rgb[3];
		double vmin = color_mode_vmin; double vmax = color_mode_vmax;
		if (color_mode_auto) { // calculate max and min values
			vmin = 1.0e20; vmax = -1.0e20;
			for (int i = 1; i <= atom.natom; i++) {
				double val;
				if (color_mode == 1) {
					val = atom.epot[i] / eV;
				}
				else if (color_mode == 2) {
					int ni = atom.nneighbor[i]; int mi = (std::min)(12, ni); mi = floor((double)mi / 2) * 2;
					val = csp(i, mi);
				}
				else if (color_mode == 3) {
					val = atom.Mises(i) / MPa;
				}
				if (vmin > val) { vmin = val; }	if (vmax < val) { vmax = val; }
			}
			color_mode_vmin = vmin; color_mode_vmax = vmax;
		}
		for (int i = 1; i <= atom.natom; i++) { // calculate values and set colors
			double val;
			if (color_mode == 1) {
				val = atom.epot[i] / eV;
			}
			else if (color_mode == 2) {
				int ni = atom.nneighbor[i]; int mi = (std::min)(12, ni); mi = floor((double)mi / 2) * 2;
				val = csp(i, mi);
			}
			else if (color_mode == 3) {
				val = atom.Mises(i) / MPa;
			}
			if (val > vmax) { val = vmax; } if (val < vmin) { val = vmin; }
			hsv[0] = 240 - 240 * (val - vmin) / (vmax - vmin);
			hsv[1] = 1.0;
			hsv[2] = 1.0;
			hsv2rgb(hsv[0], hsv[1], hsv[2], rgb[0], rgb[1], rgb[2]);
			memcpy(color[i], rgb, sizeof(GLfloat) * 4);
		}
	}
	for (int i = 0; i < 3; i++) {
		if ((marked_atom[i] > 0) && (marked_atom[i] <= atom.natom)) {
			if (marked_atom_color[i] == 0) {
				memcpy(color[marked_atom[i]], red, sizeof(GLfloat) * 4);
			}
			else if (marked_atom_color[i] == 1) {
				memcpy(color[marked_atom[i]], blue, sizeof(GLfloat) * 4);
			}
			else if (marked_atom_color[i] == 2) {
				memcpy(color[marked_atom[i]], green, sizeof(GLfloat) * 4);
			}
			else if (marked_atom_color[i] == 3) {
				memcpy(color[marked_atom[i]], yellow, sizeof(GLfloat) * 4);
			}
			else if (marked_atom_color[i] == 4) {
				memcpy(color[marked_atom[i]], purple, sizeof(GLfloat) * 4);
			}
			else if (marked_atom_color[i] == 5) {
				memcpy(color[marked_atom[i]], gray, sizeof(GLfloat) * 4);
			}
			else if (marked_atom_color[i] == 6) {
				memcpy(color[marked_atom[i]], black, sizeof(GLfloat) * 4);
			}
			else if (marked_atom_color[i] == 7) {
				memcpy(color[marked_atom[i]], white, sizeof(GLfloat) * 4);
			}
		}
	}
	for (int i = 1; i <= atom.natom; i++) { memcpy(color0[i], color[i], sizeof(GLfloat) * 4); }
}


void hsv2rgb(float h, float s, float v, float& r, float& g, float& b)
{
	float hi = floor(h / 60.0);
	float f = (h / 60.0) - hi;
	float p = v * (1.0 - s);
	float q = v * (1.0 - s * f);
	float t = v * (1.0 - (s * (1.0 - f)));
	if (hi == 0) { r = v; g = t; b = p; }
	else if (hi == 1) { r = q; g = v; b = p; }
	else if (hi == 2) { r = p; g = v; b = t; }
	else if (hi == 3) { r = p; g = q; b = v; }
	else if (hi == 4) { r = t; g = p; b = v; }
	else { r = v; g = p; b = q; }
}

#define NNMAX 30
double csp(int i, int mi)
{
	//mi=floor(mi/2)*2;
	int ix, iy, iz;
	bool flg[NNMAX];
	int nod[NNMAX]; int id[NNMAX];
	double dd[NNMAX]; double ddt[NNMAX];
	double d[NNMAX]; double d0[NNMAX];
	double dx[NNMAX]; double dy[NNMAX]; double dz[NNMAX];
	if (mi != 12) { return 0.02; }
	int ni = atom.nneighbor[i];
	for (int j = 0; j < ni; j++) {
		int jn = atom.neighbor[i][j];
		d[j] = atom.Dist2Closest(i, jn, ix, iy, iz);
		d[j] = sqrt(d[j]);
		d0[j] = d[j];
		dx[j] = atom.Dx(i, jn, ix, iy, iz);
		dy[j] = atom.Dy(i, jn, ix, iy, iz);
		dz[j] = atom.Dz(i, jn, ix, iy, iz);
	}
	sort(ni, d0, nod);
	//for (int ii=1; ii<=ni; ii++) {
	  //printf("%d %e %e %e %e\n",ii,d[ii-1],dx[ii-1],dy[ii-1],dz[ii-1]);
	  //printf("%d %e   %e %e %e\n",ii,d[nod[ii-1]],dx[nod[ii-1]],dy[nod[ii-1]],dz[nod[ii-1]]);
	//}
	//printf("\n");
	for (int j = 1; j <= mi; j++) { flg[j] = true; }
	int imin = 1;
	for (int icnt = 1; icnt <= mi / 2; icnt++) {
		int jj = 0;
		for (int j = imin + 1; j <= mi; j++) { // calc ~D_j
			if (flg[j]) {
				//jj++; ddt[jj] = dsquare(d[nod[imin]]+d[nod[j]]);
				jj++; ddt[jj] = dsquare(dx[nod[imin - 1]] + dx[nod[j - 1]])
					+ dsquare(dy[nod[imin - 1]] + dy[nod[j - 1]]) + dsquare(dz[nod[imin - 1]] + dz[nod[j - 1]]);
				id[jj] = j;
				//jj++; ddt[jj] = dsquare(dx[imin-1]+dx[j-1])
				//	+ dsquare(dy[imin-1]+dy[j-1]) + dsquare(dz[imin-1]+dz[j-1]);
				//printf("--> %d %d   %d/%d  %e  %e %e %e\n",icnt,jj, j,mi, ddt[jj],dx[imin-1]+dx[j-1],dy[imin-1]+dy[j-1],dz[imin-1]+dz[j-1]);
			}
		}
		double dmin = 1e20; int jpr = 0; // find min(~D_j) and j'
		for (int jjj = 1; jjj <= jj; jjj++) {
			if (dmin > ddt[jjj]) {
				dmin = ddt[jjj]; jpr = id[jjj];//jpr = jjj+imin;
			}
		}
		dd[icnt] = dmin;
		//printf("dd[%d]= %e imin = %d  jpr = %d\n",icnt,dmin,imin,jpr);
		//printf("  imin %e %e %e\n",dx[nod[imin-1]],dy[nod[imin-1]],dz[nod[imin-1]]);
		//printf("  jpr  %e %e %e\n",dx[nod[jpr-1]], dy[nod[jpr-1]], dz[nod[jpr-1]]);
		flg[imin] = false; flg[jpr] = false;
		if (icnt < mi / 2) {
			for (int jjj = imin; jjj < mi; jjj++) {
				imin++;
				if (flg[imin]) { break; }
			}
		}
	}
	// calc ci
	double nume = 0, deno = 0;
	for (int k = 1; k <= mi / 2; k++) { nume += dd[k]; } //printf(" %d %e\n",k,dd[k]); }
	for (int j = 1; j <= mi; j++) { deno += dsquare(d[nod[j - 1]]); }
	//for (int j=1; j<=mi;   j++) { deno += dsquare(d[j-1]); }
	//if (i==1) { printf("%d %e %e %f\n",i,nume,deno,nume/deno/2.0);}
	return nume / deno / 2.0;

}
