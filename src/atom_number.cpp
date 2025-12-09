#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include "myheader.h"

int atom_number(char* at)
{
  int type = 0;
  if (strcmp(at,"H")==0) { type = 1; }
  else if (strcmp(at,"He")==0) { type = 2; }
  else if (strcmp(at,"Li")==0) { type = 3; }
  else if (strcmp(at,"Be")==0) { type = 4; }
  else if (strcmp(at,"B")==0) { type = 5; }
  else if (strcmp(at,"C")==0) { type = 6; }
  else if (strcmp(at,"N")==0) { type = 7; }
  else if (strcmp(at,"O")==0) { type = 8; }
  else if (strcmp(at,"F")==0) { type = 9; }
  else if (strcmp(at,"Ne")==0) { type = 10; }
  else if (strcmp(at,"Na")==0) { type = 11; }
  else if (strcmp(at,"Mg")==0) { type = 12; }
  else if (strcmp(at,"Al")==0) { type = 13; }
  else if (strcmp(at,"Si")==0) { type = 14; }
  else if (strcmp(at,"P")==0) { type = 15; }
  else if (strcmp(at,"S")==0) { type = 16; }
  else if (strcmp(at,"Cl")==0) { type = 17; }
  else if (strcmp(at,"Ar")==0) { type = 18; }
  else if (strcmp(at,"K")==0) { type = 19; }
  else if (strcmp(at,"Ca")==0) { type = 20; }
  else if (strcmp(at,"Sc")==0) { type = 21; }
  else if (strcmp(at,"Ti")==0) { type = 22; }
  else if (strcmp(at,"V")==0) { type = 23; }
  else if (strcmp(at,"Cr")==0) { type = 24; }
  else if (strcmp(at,"Mn")==0) { type = 25; }
  else if (strcmp(at,"Fe")==0) { type = 26; }
  else if (strcmp(at,"Co")==0) { type = 27; }
  else if (strcmp(at,"Ni")==0) { type = 28; }
  else if (strcmp(at,"Cu")==0) { type = 29; }
  else if (strcmp(at,"Zn")==0) { type = 30; }
  else if (strcmp(at,"Ga")==0) { type = 31; }
  else if (strcmp(at,"Ge")==0) { type = 32; }
  else if (strcmp(at,"As")==0) { type = 33; }
  else if (strcmp(at,"Se")==0) { type = 34; }
  else if (strcmp(at,"Br")==0) { type = 35; }
  else if (strcmp(at,"Kr")==0) { type = 36; }
  else if (strcmp(at,"Rb")==0) { type = 37; }
  else if (strcmp(at,"Sr")==0) { type = 38; }
  else if (strcmp(at,"Y")==0) { type = 39; }
  else if (strcmp(at,"Zr")==0) { type = 40; }
  else if (strcmp(at,"Nb")==0) { type = 41; }
  else if (strcmp(at,"Mo")==0) { type = 42; }
  else if (strcmp(at,"Tc")==0) { type = 43; }
  else if (strcmp(at,"Ru")==0) { type = 44; }
  else if (strcmp(at,"Rh")==0) { type = 45; }
  else if (strcmp(at,"Pd")==0) { type = 46; }
  else if (strcmp(at,"Ag")==0) { type = 47; }
  else if (strcmp(at,"Cd")==0) { type = 48; }
  else if (strcmp(at,"In")==0) { type = 49; }
  else if (strcmp(at,"Sn")==0) { type = 50; }
  else if (strcmp(at,"Sb")==0) { type = 51; }
  else if (strcmp(at,"Te")==0) { type = 52; }
  else if (strcmp(at,"I")==0) { type = 53; }
  else if (strcmp(at,"Xe")==0) { type = 54; }
  else if (strcmp(at,"Cs")==0) { type = 55; }
  else if (strcmp(at,"Ba")==0) { type = 56; }
  else if (strcmp(at,"La")==0) { type = 57; }
  else if (strcmp(at,"Ce")==0) { type = 58; }
  else if (strcmp(at,"Pr")==0) { type = 59; }
  else if (strcmp(at,"Nd")==0) { type = 60; }
  else if (strcmp(at,"Pm")==0) { type = 61; }
  else if (strcmp(at,"Sm")==0) { type = 62; }
  else if (strcmp(at,"Eu")==0) { type = 63; }
  else if (strcmp(at,"Gd")==0) { type = 64; }
  else if (strcmp(at,"Tb")==0) { type = 65; }
  else if (strcmp(at,"Dy")==0) { type = 66; }
  else if (strcmp(at,"Ho")==0) { type = 67; }
  else if (strcmp(at,"Er")==0) { type = 68; }
  else if (strcmp(at,"Tm")==0) { type = 69; }
  else if (strcmp(at,"Yb")==0) { type = 70; }
  else if (strcmp(at,"Lu")==0) { type = 71; }
  else if (strcmp(at,"Hf")==0) { type = 72; }
  else if (strcmp(at,"Ta")==0) { type = 73; }
  else if (strcmp(at,"W")==0) { type = 74; }
  else if (strcmp(at,"Re")==0) { type = 75; }
  else if (strcmp(at,"Os")==0) { type = 76; }
  else if (strcmp(at,"Ir")==0) { type = 77; }
  else if (strcmp(at,"Pt")==0) { type = 78; }
  else if (strcmp(at,"Au")==0) { type = 79; }
  else if (strcmp(at,"Hg")==0) { type = 80; }
  else if (strcmp(at,"Tl")==0) { type = 81; }
  else if (strcmp(at,"Pb")==0) { type = 82; }
  else if (strcmp(at,"Bi")==0) { type = 83; }
  else if (strcmp(at,"Po")==0) { type = 84; }
  else if (strcmp(at,"At")==0) { type = 85; }
  else if (strcmp(at,"Rn")==0) { type = 86; }
  else if (strcmp(at,"Fr")==0) { type = 87; }
  else if (strcmp(at,"Ra")==0) { type = 88; }
  else if (strcmp(at,"Ac")==0) { type = 89; }
  else if (strcmp(at,"Th")==0) { type = 90; }
  else if (strcmp(at,"Pa")==0) { type = 91; }
  else if (strcmp(at,"U")==0) { type = 92; }
  else if (strcmp(at,"Np")==0) { type = 93; }
  else if (strcmp(at,"Pu")==0) { type = 94; }
  else if (strcmp(at,"Am")==0) { type = 95; }
  else if (strcmp(at,"Cm")==0) { type = 96; }
  else if (strcmp(at,"Bk")==0) { type = 97; }
  else if (strcmp(at,"Cf")==0) { type = 98; }
  else if (strcmp(at,"Es")==0) { type = 99; }
  else if (strcmp(at,"Fm")==0) { type = 100; }
  return type;
}
