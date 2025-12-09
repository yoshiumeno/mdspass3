#include "myheader.h"
#define dsq(a) ((a)*(a))

double Cell::Getvolume()
{
  double vol;
  vol = hmat[0][0]*(hmat[1][1]*hmat[2][2]-hmat[2][1]*hmat[1][2])
    + hmat[1][0]*(hmat[2][1]*hmat[0][2]-hmat[0][1]*hmat[2][2])
    + hmat[2][0]*(hmat[0][1]*hmat[1][2]-hmat[1][1]*hmat[0][2]);
  return vol;
}
void Cell::Setlen()
{
  len[0] = sqrt(hmat[0][0]*hmat[0][0]+hmat[1][0]*hmat[1][0]+hmat[2][0]*hmat[2][0]);
  len[1] = sqrt(hmat[0][1]*hmat[0][1]+hmat[1][1]*hmat[1][1]+hmat[2][1]*hmat[2][1]);
  len[2] = sqrt(hmat[0][2]*hmat[0][2]+hmat[1][2]*hmat[1][2]+hmat[2][2]*hmat[2][2]);
}
double Cell::Dstress()
{
  double dstr = 0;
  double x;
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      if (!cell.fix[i][j]) {
	x = fabs(sgmmat[i][j]-sgmmat_set[i][j]);
	if (x > dstr) { dstr = x; }
      }
    } }
  return dstr;
}

double Atom::MSD()
{
  double msd = 0;
  for (int i=1; i<=natom; i++) {
    msd += dsq(atom.ax[i]-atom.rx_org[i])
      +dsq(atom.ay[i]-atom.ry_org[i])
      +dsq(atom.az[i]-atom.rz_org[i]);
  }
  return msd;
}

double Atom::Enkin()
{
  double ekin;
  ekin=0.0;
  for (int i=1; i<=natom; i++) {
        ekin = ekin + 0.5 * wm[i] * ( vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i] );
  }
  return ekin;
}

double Atom::Temp()
{
  return Enkin()*2.0/3.0/(double)natom/1.380662e-23;
}

double Atom::Mises(int i)
{
  double xx
    =(satom[i][0][0]-satom[i][1][1])*(satom[i][0][0]-satom[i][1][1])
    +(satom[i][1][1]-satom[i][2][2])*(satom[i][1][1]-satom[i][2][2])
    +(satom[i][2][2]-satom[i][0][0])*(satom[i][2][2]-satom[i][0][0])
    +3.0*( satom[i][0][1]*satom[i][0][1]*2.0
	   +satom[i][1][2]*satom[i][1][2]*2.0
	   +satom[i][2][0]*satom[i][2][0]*2.0);
  return sqrt(xx/2.0);
}

double Atom::Dist2(int i, int j) // calc dr_ij
{
  double xx, yy, zz;
  xx = rx[j] - rx[i]; yy = ry[j] - ry[i]; zz = rz[j] - rz[i]; 
  dist2 = ( xx*xx + yy*yy + zz*zz );
  return dist2;
}

double Atom::Dist2(int i, int j, int ix, int iy, int iz) // calc dr_ij, j is translated by ix,iy,iz
{
  double xx, yy, zz;
  xx = rx[j]+cell.hmat[0][0]*ix+cell.hmat[0][1]*iy+cell.hmat[0][2]*iz - rx[i];
  yy = ry[j]+cell.hmat[1][0]*ix+cell.hmat[1][1]*iy+cell.hmat[1][2]*iz - ry[i];
  zz = rz[j]+cell.hmat[2][0]*ix+cell.hmat[2][1]*iy+cell.hmat[2][2]*iz - rz[i];
  dist2 = ( xx*xx + yy*yy + zz*zz );
  return dist2;
}
double Atom::Angle(int i, int j, int k, int ix0, int iy0, int iz0, int ix1, int iy1, int iz1)
{
  double rij2 = Dist2(i, j, ix0, iy0, iz0);
  double rik2 = Dist2(i, k, ix1, iy1, iz1);
  double rjk2 = Dist2(j, k, ix1-ix0, iy1-iy0, iz1-iz0);
  double costheta = (rij2+rjk2-rik2)/(2.0*sqrt(rij2*rjk2));
  double angle = acos(costheta) * 180 / M_PI;
  return angle;
}
double Atom::Dx(int i, int j, int ix, int iy, int iz)
{
  return rx[j]+cell.hmat[0][0]*ix+cell.hmat[0][1]*iy+cell.hmat[0][2]*iz - rx[i];
}
double Atom::Dy(int i, int j, int ix, int iy, int iz)
{
  return ry[j]+cell.hmat[1][0]*ix+cell.hmat[1][1]*iy+cell.hmat[1][2]*iz - ry[i];
}
double Atom::Dz(int i, int j, int ix, int iy, int iz)
{
  return rz[j]+cell.hmat[2][0]*ix+cell.hmat[2][1]*iy+cell.hmat[2][2]*iz - rz[i];
}
double Atom::Fmax()
{
  double F=0.0, F2=0.0;
  for (int i=1; i<=natom; i++) {
    double x = fx[i]*fx[i] + fy[i]*fy[i] + fz[i]*fz[i];
    if ( x > F2 ) { F2 = x; }
  }
  F = sqrt( F2 );
  return F;
}

// Functions for orthogonal cell
double Atom::Dist2Closest_ortho(int i, int j, int &ix, int &iy, int &iz)
{
  double xx, yy, zz, xxx, yyy, zzz;
  xx = rx[j] - rx[i]; yy = ry[j] - ry[i]; zz = rz[j] - rz[i];
  xxx = xx; yyy = yy; zzz = zz;
  ix = 0; iy = 0; iz = 0;
  if (cell.pbcx) {
    if (xxx < -book.frc) {
      xx += cell.hmat[0][0]; ix =  1;
    } else if (xxx > book.frc) {
      xx -= cell.hmat[0][0]; ix = -1;
    }
  }
  if (cell.pbcy) {
    if (yyy < -book.frc) {
      yy += cell.hmat[1][1]; iy =  1;
    } else if (yyy > book.frc) {
      yy -= cell.hmat[1][1]; iy = -1;
    }
  }
  if (cell.pbcz) {
    if (zzz < -book.frc) {
      zz += cell.hmat[2][2]; iz =  1;
    } else if (zzz > book.frc) {
      zz -= cell.hmat[2][2]; iz = -1;
    }
  }
  dist2 = ( xx*xx + yy*yy + zz*zz );
  return dist2;
}
double Atom::Dist2Closest_yz_ortho(int i, int j, int ix, int &iy, int &iz)
{
  double xx, yy, zz, xxx, yyy, zzz;
  xx = rx[j] - rx[i]; yy = ry[j] - ry[i]; zz = rz[j] - rz[i];
  xxx = xx; yyy = yy; zzz = zz;
  iy = 0; iz = 0;
  if (cell.pbcx) {
    xx += cell.hmat[0][0]*(double)ix;
  }
  if (cell.pbcy) {
    if (yyy < -book.frc) {
      yy += cell.hmat[1][1]; iy =  1;
    } else if (yyy > book.frc) {
      yy -= cell.hmat[1][1]; iy = -1;
    }
  }
  if (cell.pbcz) {
    if (zzz < -book.frc) {
      zz += cell.hmat[2][2]; iz =  1;
    } else if (zzz > book.frc) {
      zz -= cell.hmat[2][2]; iz = -1;
    }
  }
  dist2 = ( xx*xx + yy*yy + zz*zz );
  return dist2;
}
double Atom::Dist2Closest_xz_ortho(int i, int j, int &ix, int iy, int &iz)
{
  double xx, yy, zz, xxx, yyy, zzz;
  xx = rx[j] - rx[i]; yy = ry[j] - ry[i]; zz = rz[j] - rz[i];
  xxx = xx; yyy = yy; zzz = zz;
  ix = 0; iz = 0;
  if (cell.pbcx) {
    if (xxx < -book.frc) {
      xx += cell.hmat[0][0]; ix =  1;
    } else if (xxx > book.frc) {
      xx -= cell.hmat[0][0]; ix = -1;
    }
  }
  if (cell.pbcy) {
    xx += cell.hmat[0][1]*(double)iy;
    yy += cell.hmat[1][1]*(double)iy;
    zz += cell.hmat[2][1]*(double)iy;
  }
  if (cell.pbcz) {
    if (zzz < -book.frc) {
      zz += cell.hmat[2][2]; iz =  1;
    } else if (zzz > book.frc) {
      zz -= cell.hmat[2][2]; iz = -1;
    }
  }
  dist2 = ( xx*xx + yy*yy + zz*zz );
  return dist2;
}
double Atom::Dist2Closest_xy_ortho(int i, int j, int &ix, int &iy, int iz)
{
  double xx, yy, zz, xxx, yyy, zzz;
  xx = rx[j] - rx[i]; yy = ry[j] - ry[i]; zz = rz[j] - rz[i];
  xxx = xx; yyy = yy; zzz = zz;
  ix = 0; iy = 0;
  if (cell.pbcx) {
    if (xxx < -book.frc) {
      xx = xx+cell.hmat[0][0]; ix =  1;
    } else if (xxx > book.frc) {
      xx = xx-cell.hmat[0][0]; ix = -1;
    }
  }
  if (cell.pbcy) {
    if (yyy < -book.frc) {
      yy += cell.hmat[1][1]; iy =  1;
    } else if (yyy > book.frc) {
      yy -= cell.hmat[1][1]; iy = -1;
    }
  }
  if (cell.pbcz) {
    xx += cell.hmat[0][2]*(double)iz;
    yy += cell.hmat[1][2]*(double)iz;
    zz += cell.hmat[2][2]*(double)iz;
  }
  dist2 = ( xx*xx + yy*yy + zz*zz );
  return dist2;
}

// Functions for non-orthogonal cell
double Atom::Dist2Closest(int i, int j, int &ix, int &iy, int &iz)
{
  double xx, yy, zz, xxx, yyy, zzz;
  xx  = rx[j] - rx[i]; yy  = ry[j] - ry[i]; zz  = rz[j] - rz[i];
  xxx = qx[j] - qx[i]; yyy = qy[j] - qy[i]; zzz = qz[j] - qz[i];
  ix = 0; iy = 0; iz = 0;
  if (cell.pbcx) {
    if (xxx < -book.frc/cell.len[0]) {
      xx += cell.hmat[0][0]; yy += cell.hmat[1][0]; zz += cell.hmat[2][0]; ix =  1;
    } else if (xxx > book.frc/cell.len[0]) {
      xx -= cell.hmat[0][0]; yy -= cell.hmat[1][0]; zz -= cell.hmat[2][0]; ix = -1;
    }
  }
  if (cell.pbcy) {
    if (yyy < -book.frc/cell.len[1]) {
      xx += cell.hmat[0][1]; yy += cell.hmat[1][1]; zz += cell.hmat[2][1]; iy =  1;
    } else if (yyy > book.frc/cell.len[1]) {
      xx -= cell.hmat[0][1]; yy -= cell.hmat[1][1]; zz -= cell.hmat[2][1]; iy = -1;
    }
  }
  if (cell.pbcz) {
    if (zzz < -book.frc/cell.len[2]) {
      xx += cell.hmat[0][2]; yy += cell.hmat[1][2]; zz += cell.hmat[2][2]; iz =  1;
    } else if (zzz > book.frc/cell.len[2]) {
      xx -= cell.hmat[0][2]; yy -= cell.hmat[1][2]; zz -= cell.hmat[2][2]; iz = -1;
    }
  }
  dist2 = ( xx*xx + yy*yy + zz*zz );
  return dist2;
}
double Atom::Dist2Closest_yz(int i, int j, int ix, int &iy, int &iz)
{
  double xx, yy, zz, xxx, yyy, zzz;
  xx  = rx[j] - rx[i]; yy  = ry[j] - ry[i]; zz  = rz[j] - rz[i];
  xxx = qx[j] - qx[i]; yyy = qy[j] - qy[i]; zzz = qz[j] - qz[i];
  iy = 0; iz = 0;
  if (cell.pbcx) {
    xx += cell.hmat[0][0]*(double)ix;
    yy += cell.hmat[1][0]*(double)ix;
    zz += cell.hmat[2][0]*(double)ix;
  }
  if (cell.pbcy) {
    if (yyy < -book.frc/cell.len[1]) {
      xx += cell.hmat[0][1]; yy += cell.hmat[1][1]; zz += cell.hmat[2][1]; iy =  1;
    } else if (yyy > book.frc/cell.len[1]) {
      xx -= cell.hmat[0][1]; yy -= cell.hmat[1][1]; zz -= cell.hmat[2][1]; iy = -1;
    }
  }
  if (cell.pbcz) {
    if (zzz < -book.frc/cell.len[2]) {
      xx += cell.hmat[0][2]; yy += cell.hmat[1][2]; zz += cell.hmat[2][2]; iz =  1;
    } else if (zzz > book.frc/cell.len[2]) {
      xx -= cell.hmat[0][2]; yy -= cell.hmat[1][2]; zz -= cell.hmat[2][2]; iz = -1;
    }
  }
  dist2 = ( xx*xx + yy*yy + zz*zz );
  return dist2;
}
double Atom::Dist2Closest_xz(int i, int j, int &ix, int iy, int &iz)
{
  double xx, yy, zz, xxx, yyy, zzz;
  xx  = rx[j] - rx[i]; yy  = ry[j] - ry[i]; zz  = rz[j] - rz[i];
  xxx = qx[j] - qx[i]; yyy = qy[j] - qy[i]; zzz = qz[j] - qz[i];
  ix = 0; iz = 0;
  if (cell.pbcx) {
    if (xxx < -book.frc/cell.len[0]) {
      xx += cell.hmat[0][0]; yy += cell.hmat[1][0]; zz += cell.hmat[2][0]; ix =  1;
    } else if (xxx > book.frc/cell.len[0]) {
      xx -= cell.hmat[0][0]; yy -= cell.hmat[1][0]; zz -= cell.hmat[2][0]; ix = -1;
    }
  }
  if (cell.pbcy) {
    xx += cell.hmat[0][1]*(double)iy;
    yy += cell.hmat[1][1]*(double)iy;
    zz += cell.hmat[2][1]*(double)iy;
  }
  if (cell.pbcz) {
    if (zzz < -book.frc/cell.len[2]) {
      xx += cell.hmat[0][2]; yy += cell.hmat[1][2]; zz += cell.hmat[2][2]; iz =  1;
    } else if (zzz > book.frc/cell.len[2]) {
      xx -= cell.hmat[0][2]; yy -= cell.hmat[1][2]; zz -= cell.hmat[2][2]; iz = -1;
    }
  }
  dist2 = ( xx*xx + yy*yy + zz*zz );
  return dist2;
}
double Atom::Dist2Closest_xy(int i, int j, int &ix, int &iy, int iz)
{
  double xx, yy, zz, xxx, yyy, zzz;
  xx  = rx[j] - rx[i]; yy  = ry[j] - ry[i]; zz  = rz[j] - rz[i];
  xxx = qx[j] - qx[i]; yyy = qy[j] - qy[i]; zzz = qz[j] - qz[i];
  ix = 0; iy = 0;
  if (cell.pbcx) {
    if (xxx < -book.frc/cell.len[0]) {
      xx += cell.hmat[0][0]; yy += cell.hmat[1][0]; zz += cell.hmat[2][0]; ix =  1;
    } else if (xxx > book.frc/cell.len[0]) {
      xx -= cell.hmat[0][0]; yy -= cell.hmat[1][0]; zz -= cell.hmat[2][0]; ix = -1;
    }
  }
  if (cell.pbcy) {
    if (yyy < -book.frc/cell.len[1]) {
      xx += cell.hmat[0][1]; yy += cell.hmat[1][1]; zz += cell.hmat[2][1]; iy =  1;
    } else if (yyy > book.frc/cell.len[1]) {
      xx -= cell.hmat[0][1]; yy -= cell.hmat[1][1]; zz -= cell.hmat[2][1]; iy = -1;
    }
  }
  if (cell.pbcz) {
    xx += cell.hmat[0][2]*(double)iz;
    yy += cell.hmat[1][2]*(double)iz;
    zz += cell.hmat[2][2]*(double)iz;
  }
  dist2 = ( xx*xx + yy*yy + zz*zz );
  return dist2;
}

void Atom::Calcq()
{
  for (int i=1; i<=natom; i++) {
    qx[i] = cell.hinmat[0][0]*rx[i] + cell.hinmat[0][1]*ry[i] + cell.hinmat[0][2]*rz[i];
    qy[i] = cell.hinmat[1][0]*rx[i] + cell.hinmat[1][1]*ry[i] + cell.hinmat[1][2]*rz[i];
    qz[i] = cell.hinmat[2][0]*rx[i] + cell.hinmat[2][1]*ry[i] + cell.hinmat[2][2]*rz[i];
  }
}
void Atom::Calcq(int i)
{
  qx[i] = cell.hinmat[0][0]*rx[i] + cell.hinmat[0][1]*ry[i] + cell.hinmat[0][2]*rz[i];
  qy[i] = cell.hinmat[1][0]*rx[i] + cell.hinmat[1][1]*ry[i] + cell.hinmat[1][2]*rz[i];
  qz[i] = cell.hinmat[2][0]*rx[i] + cell.hinmat[2][1]*ry[i] + cell.hinmat[2][2]*rz[i];
}

Atom::Atom()
{
  QC = false;
  instcenter = 1;
  natom = -1;
  for (int i=0; i<MAXPOTTYPE; i++) { for (int j=0; j<MAXPOTARGTYPE; j++) {
      strcpy(potargstring_list[i][j], ""); } }
  strcpy(potstring_list[ 0], "Morse");     potarg_number[ 0] = 0; potarg_readable[ 0] = false;
  strcpy(potstring_list[ 1], "Buckingham");potarg_number[ 1] = 0; potarg_readable[ 1] = false;
  strcpy(potstring_list[ 2], "EAM");       potarg_number[ 2] = 2; potarg_readable[ 2] = false;
  strcpy(potargstring_list[ 2][0], "GEAM");
  strcpy(potargstring_list[ 2][1], "Mishin");
  strcpy(potstring_list[ 3], "Tersoff");   potarg_number[ 3] = 6; potarg_readable[ 3] = false;
  strcpy(potargstring_list[ 3][0], "Si(B)");
  strcpy(potargstring_list[ 3][1], "Si(C)");
  strcpy(potargstring_list[ 3][2], "Si(B*)");
  strcpy(potargstring_list[ 3][3], "C_Si_Ge");
  strcpy(potargstring_list[ 3][4], "B_N_C");
  strcpy(potargstring_list[ 3][5], "SiC_Erhart");
  strcpy(potstring_list[ 4], "Brenner");   potarg_number[ 4] = 0; potarg_readable[ 4] = false;
  strcpy(potstring_list[ 5], "AIREBO");    potarg_number[ 5] = 0; potarg_readable[ 5] = false;
  strcpy(potstring_list[ 6], "AIREBO_BN"); potarg_number[ 6] = 0; potarg_readable[ 6] = false;
  strcpy(potstring_list[ 7], "Dipole");    potarg_number[ 7] = 0; potarg_readable[ 7] = true;
  strcpy(potstring_list[ 8], "ADP");       potarg_number[ 8] = 0; potarg_readable[ 8] = true;
  strcpy(potstring_list[ 9], "NiYSZ");     potarg_number[ 9] = 0; potarg_readable[ 9] = false;
  strcpy(potstring_list[10], "SW");        potarg_number[10] = 4; potarg_readable[10] = false;
  strcpy(potargstring_list[10][0], "GaN_Bere");
  strcpy(potargstring_list[10][1], "GaN_Aichoune");
  strcpy(potargstring_list[10][2], "SiC_Vashishta");
  strcpy(potargstring_list[10][3], "Si_Stillinger");
  strcpy(potstring_list[11], "ReaxFF");    potarg_number[11] = 0; potarg_readable[11] = true;
  strcpy(potstring_list[12], "Sinusoidal");potarg_number[12] = 0; potarg_readable[12] = false;
  strcpy(potstring_list[13], "ShellModel");potarg_number[13] = 0; potarg_readable[13] = false;
}

Book::Book()
{
  frc  = 10.0e-10;
  frc2 = frc*frc;
  nbk  = 50;
  nbook = 1000;
  nbook_new = nbook;
  natom = 0;
  alloc = true;
}

Cell::Cell()
{
  //  pbcx=true; pbcy=true; pbcz=true;
  pbcx=false; pbcy=false; pbcz=false;
  prmass = 1.0e-24;
  damper_param = 1.0e13;
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      hvmat[i][j] = 0; hamat[i][j] = 0; hbmat[i][j] = 0;
      hcmat[i][j] = 0; sgmmat_set[i][j] = 0;
    } }
  s = 1.0; sv = 0; sa = 0; sb = 0; sc = 0;
  nhmass = 1.0e-45;
}
void Cell::Reset()
{
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      hvmat[i][j] = 0; hamat[i][j] = 0; hbmat[i][j] = 0;
      hcmat[i][j] = 0;
    } }
  s = 1.0; sv = 0; sa = 0; sb = 0; sc = 0;
}

Cond::Cond()
{
}

Fire::Fire()
{
  sp = 0.0; fnorm = 0.0; vnorm = 0.0;
  fire_alph = 0.10;
  fire_alph_ini = 0.10;
  nfmin = 5;
  ffinc = 1.1; ffdec = 0.5; ffalph = 0.99;
  ifire = 0;
  ffdtmax = 1.0e-14;
}

Cg::Cg()
{
  lam_init = 1.0e-3;
}

GEAM::GEAM()
{
  gre[0] = 2.556162; gre[1] = 2.891814; gre[2] = 2.885034; gre[3] = 2.488746; gre[4] = 2.750897;
  gre[5] = 2.771916; gre[6] = 2.886166; gre[7] = 3.499723; gre[8] = 2.481987; gre[9] = 2.728100;
  gre[10]= 2.860082; gre[11]= 2.740840; gre[12]= 3.196291; gre[13]= 2.505979; gre[14] =2.933872;
  gre[15]= 3.199978;

  gfe[0] = 1.554485; gfe[1] = 1.106232; gfe[2] = 1.529021; gfe[3] = 2.007018; gfe[4] = 1.595417;
  gfe[5] = 2.336509; gfe[6] = 1.392302; gfe[7] = 0.647872; gfe[8] = 1.885957; gfe[9] = 2.723710;
  gfe[10]= 3.086341; gfe[11]= 3.487340; gfe[12]= 0.544323; gfe[13]= 1.975299; gfe[14]= 1.863200;
  gfe[15]= 2.230909;

  grhoe[0] = 22.150141; grhoe[1] = 15.539255; grhoe[2] = 21.319637; grhoe[3] = 27.984706;
  grhoe[4] = 22.770550; grhoe[5] = 34.108882; grhoe[6] = 20.226537; grhoe[7] =  8.906840;
  grhoe[8] = 20.041463; grhoe[9] = 29.354065; grhoe[10]= 33.787168; grhoe[11]= 37.234847;
  grhoe[12]=  7.132600; grhoe[13]= 27.206789; grhoe[14]= 25.565138; grhoe[15]= 30.879991;

  grhos[0] = 22.150141; grhos[1] = 15.539255; grhos[2] = 21.319637; grhos[3] = 27.984706;
  grhos[4] = 22.770550; grhos[5] = 34.108882; grhos[6] = 20.226537; grhos[7] =  8.906840;
  grhos[8] = 20.041463; grhos[9] = 29.354065; grhos[10]= 33.787168; grhos[11]= 37.234847;
  grhos[12]=  7.132600; grhos[13]= 27.206789; grhos[14]= 25.565138; grhos[15]= 30.879991;

  galp[0] = 7.669911; galp[1] = 7.944536; galp[2] = 8.086176; galp[3] = 8.029633;
  galp[4] = 7.605017; galp[5] = 7.079952; galp[6] = 6.942419; galp[7] = 8.468412;
  galp[8] = 9.818270; galp[9] = 8.393531; galp[10]= 8.489528; galp[11]= 8.900114;
  galp[12]=10.228708; galp[13]= 8.679625; galp[14]= 8.775431; galp[15]= 8.559190;

  gbet[0] = 4.090619; gbet[1] = 4.237086; gbet[2] = 4.312627; gbet[3] = 4.282471;
  gbet[4] = 4.056009; gbet[5] = 3.775974; gbet[6] = 3.702623; gbet[7] = 4.516486;
  gbet[8] = 5.236411; gbet[9] = 4.476550; gbet[10]= 4.527748; gbet[11]= 4.746728;
  gbet[12]= 5.455311; gbet[13]= 4.629134; gbet[14]= 4.680230; gbet[15]= 4.564902;

  gaa[0] = 0.327584; gaa[1] = 0.266074; gaa[2] = 0.230728; gaa[3] = 0.439664;
  gaa[4] = 0.385412; gaa[5] = 0.449644; gaa[6] = 0.251519; gaa[7] = 0.134878;
  gaa[8] = 0.392811; gaa[9] = 0.708787; gaa[10]= 0.611679; gaa[11]= 0.882435;
  gaa[12]= 0.137518; gaa[13]= 0.421378; gaa[14]= 0.373601; gaa[15]= 0.424667;

  gbb[0] = 0.468735; gbb[1] = 0.386272; gbb[2] = 0.336695; gbb[3] = 0.632771;
  gbb[4] = 0.545121; gbb[5] = 0.593713; gbb[6] = 0.313394; gbb[7] = 0.203093;
  gbb[8] = 0.646243; gbb[9] = 1.120373; gbb[10]= 1.032101; gbb[11]= 1.394592;
  gbb[12]= 0.225930; gbb[13]= 0.640107; gbb[14]= 0.570968; gbb[15]= 0.640054;
  
  gkai[0] = 0.431307; gkai[1] = 0.425351; gkai[2] = 0.420755; gkai[3] = 0.413436;
  gkai[4] = 0.425578; gkai[5] = 0.413484; gkai[6] = 0.395132; gkai[7] = 0.425877;
  gkai[8] = 0.170306; gkai[9] = 0.137640; gkai[10]= 0.176977; gkai[11]= 0.139209;
  gkai[12]= 0.5; gkai[13]= 0.5; gkai[14]= 0.5; gkai[15]= 0.5;

  glam[0] = 0.862140; glam[1] = 0.850703; glam[2] = 0.841511; glam[3] = 0.826873;
  glam[4] = 0.851156; glam[5] = 0.826967; glam[6] = 0.790264; glam[7] = 0.851753;
  glam[8] = 0.340613; glam[9] = 0.275280; glam[10]= 0.353954; glam[11]= 0.278417;
  glam[12]= 1.0; glam[13]= 1.0; glam[14]= 1.0; glam[15] = 1.0;

  gffn0[0] = -2.176490; gffn0[1] = -1.729619; gffn0[2] = -2.930281; gffn0[3] = -2.693996;
  gffn0[4] = -2.320473; gffn0[5] = -4.099542; gffn0[6] = -2.806783; gffn0[7] = -1.419644;
  gffn0[8] = -2.534992; gffn0[9] = -3.692913; gffn0[10]= -5.103845; gffn0[11]= -4.946281;
  gffn0[12]= -0.896473; gffn0[13]= -2.541799; gffn0[14]= -3.203773; gffn0[15]= -4.485793;

  gffn1[0] = -0.140035; gffn1[1] = -0.221025; gffn1[2] = -0.554034; gffn1[3] = -0.066073;
  gffn1[4] = -0.421263; gffn1[5] = -0.754764; gffn1[6] = -0.276173; gffn1[7] = -0.228622;
  gffn1[8] = -0.059605; gffn1[9] = -0.178812; gffn1[10]= -0.405524; gffn1[11]= -0.148818;
  gffn1[12]= -0.044291; gffn1[13]= -0.219415; gffn1[14]= -0.198262; gffn1[15]= -0.293129;

  gffn2[0] = 0.285621; gffn2[1] = 0.541558; gffn2[2] = 1.489437; gffn2[3] = 0.170482;
  gffn2[4] = 0.966525; gffn2[5] = 1.766503; gffn2[6] = 0.893409; gffn2[7] = 0.630069;
  gffn2[8] = 0.193065; gffn2[9] = 0.380450; gffn2[10]= 1.112997; gffn2[11]= 0.365057;
  gffn2[12]= 0.162232; gffn2[13]= 0.733381; gffn2[14]= 0.683779; gffn2[15]= 0.990148;

  gffn3[0] = -1.750834; gffn3[1] = -0.967036; gffn3[2] = -0.886809; gffn3[3] = -2.457442;
  gffn3[4] = -0.932685; gffn3[5] = -1.578274; gffn3[6] = -1.637201; gffn3[7] = -0.560952;
  gffn3[8] = -2.282322; gffn3[9] = -3.133650; gffn3[10]= -3.585325; gffn3[11]= -4.432406;
  gffn3[12]= -0.689950; gffn3[13]= -1.589003; gffn3[14]= -2.321732; gffn3[15]= -3.202516;

  gff0[0] = -2.19; gff0[1] = -1.75; gff0[2] = -2.98; gff0[3] = -2.70;
  gff0[4] = -2.36; gff0[5] = -4.17; gff0[6] = -2.83; gff0[7] = -1.44;
  gff0[8] = -2.54; gff0[9] = -3.71; gff0[10]= -5.14; gff0[11]= -4.96;
  gff0[12]= -0.90; gff0[13]= -2.56; gff0[14]= -3.22; gff0[15]= -4.51;

  gff1[0] = 0.0; gff1[1] = 0.0; gff1[2] = 0.0; gff1[3] = 0.0;
  gff1[4] = 0.0; gff1[5] = 0.0; gff1[6] = 0.0; gff1[7] = 0.0;
  gff1[8] = 0.0; gff1[9] = 0.0; gff1[10]= 0.0; gff1[11]= 0.0;
  gff1[12]= 0.0; gff1[13]= 0.0; gff1[14]= 0.0; gff1[15] = 0.0;

  gff2[0] = 0.702991; gff2[1] = 0.983967; gff2[2] = 2.283863; gff2[3] = 0.282257;
  gff2[4] = 1.966273; gff2[5] = 3.474733; gff2[6] = 0.929508; gff2[7] = 0.921049;
  gff2[8] = 0.200269; gff2[9] = 0.875874; gff2[10]= 1.640098; gff2[11]= 0.661935;
  gff2[12]= 0.122838; gff2[13]= 0.705845; gff2[14]= 0.608587; gff2[15]= 0.928602;
  
  gff3[0] = 0.683705; gff3[1] =  0.520904; gff3[2] =  0.494127; gff3[3] =  0.102879;
  gff3[4] = 1.396717; gff3[5] =  2.288323; gff3[6] = -0.682320; gff3[7] =  0.108847;
  gff3[8] =-0.148770; gff3[9] =  0.776222; gff3[10]=  0.221375; gff3[11]=  0.348147;
  gff3[12]=-0.226010; gff3[13]= -0.687140; gff3[14]= -0.750710; gff3[15]= -0.981870;
  
  geta[0] = 0.921150; geta[1] = 1.149461; geta[2] = 1.286960; geta[3] =  0.509860;
  geta[4] = 1.399758; geta[5] = 1.393490; geta[6] = 0.779208; geta[7] =  1.172361;
  geta[8] = 0.391750; geta[9] = 0.790879; geta[10]= 0.848843; geta[11]= -0.582714;
  geta[12]= 0.431425; geta[13]= 0.694608; geta[14]= 0.558572; geta[15]=  0.597133;
  
  gffe[0] = -2.191675; gffe[1] = -1.751274; gffe[2] = -2.981365; gffe[3] = -2.700493;
  gffe[4] = -2.362609; gffe[5] = -4.174332; gffe[6] = -2.829437; gffe[7] = -1.440494;
  gffe[8] = -2.539945; gffe[9] = -3.712093; gffe[10]= -5.141526; gffe[11]= -4.961306;
  gffe[12]= -0.899702; gffe[13]= -2.559307; gffe[14]= -3.219176; gffe[15]= -4.509025;

  initialize = true;
}

Tersoff::Tersoff()
{
  initialize = true;
  maxnei = 50;
  large = false;
  debug_large = false;
}

Sw::Sw(void)
{ // Constructor
  initialize = true;
  for(int i=0;i<200;i++){AtomNum2Index[i] = 0;}
  nSpec = 0;
  nPair = 0;
  nTrio = 0;
  Phi = NULL;
  R   = NULL;
  P   = NULL;
  PairType = NULL;
  TrioType = NULL;
}
void Sw::DeleteArray(void)
{
  if(nSpec>0){
    DeletePairType();
    DeleteTrioType();
  } // if: nSpec
  
  if(nPair>0){
    DeletePhi();
  } // if: nPair
  
  if(nTrio>0){
    DeleteR();
    DeleteP();
  } // if: nTrio
}
void Sw::AllocateArray(void)
{
  AllocatePairType();
  AllocateTrioType();
  AllocatePhi();
  AllocateR();
  AllocateP();
}
void Sw::DeletePhi(void)
{
  if(Phi!=NULL){
    for(int i=0;i<nPair+1;i++){
      if(Phi[i]!=NULL){
	delete Phi[i];
	Phi[i] = NULL;
      } // if: Phi[i]
    } // for: i
    delete [] Phi;
    Phi = NULL;
  } //if: Phi
}
void Sw::DeleteR(void)
{
  if(R!=NULL){
    for(int i=0;i<nTrio+1;i++){
      if(R[i]!=NULL){
	delete R[i];
	R[i] = NULL;
      } // if: R[i]
    } // for: i
    delete [] R;
    R = NULL;
  } // if: R
}
void Sw::DeleteP(void)
{
  if(P!=NULL){
    for(int i=0;i<nTrio+1;i++){
      if(P[i]!=NULL){
	delete P[i];
	P[i] = NULL;
      } // if: P[i]
    } // for: i
    delete [] P;
    P = NULL;
  } // if: P
}
void Sw::DeletePairType(void)
{
  if(PairType!=NULL){
    for(int i=0;i<nSpec+1;i++){
      if(PairType[i]!=NULL){
	delete [] PairType[i];
	PairType[i] = NULL;
      } // if: PairType[i]
    } // for: i
    delete [] PairType;
    PairType = NULL;
  } //if: PairType
}
void Sw::DeleteTrioType(void)
{
  if(TrioType!=NULL){
    for(int i=0;i<nSpec+1;i++){
      if(TrioType[i]!=NULL){
	for(int j=0;j<nSpec+1;j++){
	  if(TrioType[i][j]!=NULL){
	    delete [] TrioType[i][j];
	    TrioType[i][j] = NULL;
	  } // if: TrioType[i][j]
	} // for: j
	delete [] TrioType[i];
	TrioType[i] = NULL;
      } // if: TrioType[i]
    } // for: i
    delete [] TrioType;
    TrioType = NULL;
  } //if: PairType
}
void Sw::AllocatePairType(void)
{
  if(nSpec==0){return;}
  PairType = new int*[nSpec+1];
  for(int i=0;i<nSpec+1;i++){
    PairType[i] = new int[nSpec+1];
    for(int j=0;j<nSpec+1;j++){
      PairType[i][j] = 0;
    } // for: j
  } // for: i
}
void Sw::AllocateTrioType(void)
{
  if(nSpec==0){return;}
  TrioType = new int**[nSpec+1];
  for(int i=0;i<nSpec+1;i++){
    TrioType[i] = new int*[nSpec+1];
    for(int j=0;j<nSpec+1;j++){
      TrioType[i][j] = new int[nSpec+1];
      for(int k=0;k<nSpec+1;k++){
	TrioType[i][j][k] = 0;
      } // for: k
    } // for: j
  } // for: i
}
void Sw::AllocatePhi(void)
{
  if(nPair==0){return;}
  Phi = new PairFunction*[nPair+1];
  for(int i=0;i<nPair+1;i++){
    Phi[i] = NULL; // Instances are created later...
  } // for: i
}
void Sw::AllocateR(void)
{
  if(nTrio==0){return;}
  R = new PairFunction*[nTrio+1];
  for(int i=0;i<nTrio+1;i++){
    R[i] = NULL; // Instances are created later...
  } // for: i
}
void Sw::AllocateP(void)
{
  if(nTrio==0){return;}
  P = new PairFunction*[nTrio+1];
  for(int i=0;i<nTrio+1;i++){
    P[i] = NULL; // Instances are created later...
  } // for: i
}

