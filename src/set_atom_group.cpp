#include <string>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#if defined(_WIN32) || defined(_WIN64)
#include <Windows.h>
#endif
#include <GL/gl.h>
#include "myheader.h"

void set_atom_group(int group, int axis, float range_max, float range_min)
{
    printf("Setting atoms into group:\n");
    if (axis == 0) printf("qx from %f to %f --> group %d\n",range_min, range_max, group);
    if (axis == 1) printf("qy from %f to %f --> group %d\n",range_min, range_max, group);
    if (axis == 2) printf("qz from %f to %f --> group %d\n",range_min, range_max, group);
    atom.Calcq();
    for (int i = 1; i <= atom.natom; i++) {
        if (axis == 0) { // x-axis
            if ((atom.qx[i] >= range_min)&&(atom.qx[i] <= range_max)) {
                atom.group[i] = group;
            }
        } else if (axis == 1) { // y-axis
            if ((atom.qy[i] >= range_min)&&(atom.qy[i] <= range_max)) {
                atom.group[i] = group;
            }
        } else if (axis == 2) { // z-axis
            if ((atom.qz[i] >= range_min)&&(atom.qz[i] <= range_max)) {
                atom.group[i] = group;
            }
        }
    }
}

void reset_atom_group()
{
    printf("Atom group has been reset.\n");
    for (int i = 1; i <= atom.natom; i++) {
        atom.group[i] = 0;
    }
}

void set_atom_group_visibility()
{
    printf("Changing visuality for atom grouping\n");
    if (if_group_num_show) {
        for (int i = 1; i <= atom.natom; i++) {
            if (atom.group[i] == group_num_show) {
                atom.visible[i] = true; } else { atom.visible[i] = false; }
        }
    } else {
        for (int i = 1; i <= atom.natom; i++) {
            atom.visible[i] = true;
        }
    }
}