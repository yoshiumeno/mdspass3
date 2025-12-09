#include <stdio.h>
#include <math.h>
#include "myheader.h"

void outproduct (double a[3], double b[3], double x[3]);
double innerproduct(double a[3], double b[3]);

void inverse(double h[3][3], double hi[3][3])
{
  double deth
    = h[0][0]*h[1][1]*h[2][2] + h[0][1]*h[1][2]*h[2][0] 
    + h[0][2]*h[1][0]*h[2][1] - h[0][0]*h[1][2]*h[2][1] 
    - h[0][1]*h[1][0]*h[2][2] - h[0][2]*h[1][1]*h[2][0];
    
  hi[0][0]=(h[1][1]*h[2][2]-h[1][2]*h[2][1])/deth;
  hi[0][1]=(h[0][2]*h[2][1]-h[0][1]*h[2][2])/deth;
  hi[0][2]=(h[0][1]*h[1][2]-h[0][2]*h[1][1])/deth;
  hi[1][0]=(h[1][2]*h[2][0]-h[1][0]*h[2][2])/deth;
  hi[1][1]=(h[0][0]*h[2][2]-h[0][2]*h[2][0])/deth;
  hi[1][2]=(h[0][2]*h[1][0]-h[0][0]*h[1][2])/deth;
  hi[2][0]=(h[1][0]*h[2][1]-h[1][1]*h[2][0])/deth;
  hi[2][1]=(h[0][1]*h[2][0]-h[0][0]*h[2][1])/deth;
  hi[2][2]=(h[0][0]*h[1][1]-h[0][1]*h[1][0])/deth;
}

void resetvec(double a[3])
{
  for (int i=0; i<3; i++) {
    a[i]=0;
  }
}

void resetmat(double a[3][3])
{
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      a[i][j]=0; } }
}

void resetmat6(double a[6][6])
{
  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      a[i][j]=0; } }
}

void transpose(double a[3][3])
{
  double x;
  for (int i=0; i<3; i++) {
    for (int j=i+1; j<3; j++) {
      x=a[i][j]; a[i][j]=a[j][i]; a[j][i]=x; } }
}
void transpose(double a[3][3], double b[3][3])
{
  double x;
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      b[i][j]=a[j][i]; } }
}

//c = product(a,b)
void matmul(double a[3][3], double b[3][3], double c[3][3])
{
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      c[i][j]=0.0;
      for (int k=0; k<3; k++) {
	c[i][j] += a[i][k]*b[k][j];
      } } }
}
void matadd(double a[3][3], double b[3][3], double c[3][3])
{
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      c[i][j]=a[i][j]+b[i][j];
    } }
}
void matsub(double a[3][3], double b[3][3], double c[3][3])
{
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      c[i][j]=a[i][j]-b[i][j];
    } }
}
void matcpy(double a[3][3], double b[3][3])
{
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      b[i][j] = a[i][j];
    } }
}

void recipro(double h[3][3], double r[3][3])
{
  double a[3], b[3], c[3], x[3];
  for (int i=0; i<3; i++) {
    a[i] = h[i][0];
    b[i] = h[i][1];
    c[i] = h[i][2];
  }
  outproduct(b,c,x);
  for (int i=0; i<3; i++) {
    r[i][0] = x[i];
  }
  outproduct(c,a,x);
  for (int i=0; i<3; i++) {
    r[i][1] = x[i];
  }
  outproduct(a,b,x);
  for (int i=0; i<3; i++) {
    r[i][2] = x[i];
  }
}

void outproduct (double a[3], double b[3], double x[3])
{
  x[0] = a[1]*b[2]-a[2]*b[1];
  x[1] = a[2]*b[0]-a[0]*b[2];
  x[2] = a[0]*b[1]-a[1]*b[0];
}

double innerproduct(double a[3], double b[3])
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

double veclength(double x[3])
{
  return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}

void tensor_rotation(double s[3][3],  double x1[3],  double x2[3],  double x3[3],
		     double sp[3][3], double xp1[3], double xp2[3], double xp3[3])
{
  double q[3][3];
  q[0][0] = innerproduct(xp1,x1);
  q[0][1] = innerproduct(xp1,x2);
  q[0][2] = innerproduct(xp1,x3);
  q[1][0] = innerproduct(xp2,x1);
  q[1][1] = innerproduct(xp2,x2);
  q[1][2] = innerproduct(xp2,x3);
  q[2][0] = innerproduct(xp3,x1);
  q[2][1] = innerproduct(xp3,x2);
  q[2][2] = innerproduct(xp3,x3);
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      sp[i][j] = 0;
      for (int m=0; m<3; m++) {
	for (int n=0; n<3; n++) {
	  sp[i][j] += s[m][n]*q[i][m]*q[j][n];
	}
      }
    }
  }
}

double theta(double x, double y)
{
  double t;
  if (x == 0) {
    if (y > 0) { t = M_PI/2.0; }
    if (y < 0) { t = M_PI/2.0*3.0; }
  } else {
    t = atan(y/x);
    if ((y>=0)&&(x<0)) { t += M_PI; }
    if ((y<0)&&(x<0)) { t += M_PI; }
    if ((y<0)&&(x>0)) { t += M_PI*2; }
  }
  return t;
}

// product of quaternions  r <- p x q
void qmul(float r[], const float p[], const float q[])
{
  r[0] = p[0] * q[0] - p[1] * q[1] - p[2] * q[2] - p[3] * q[3];
  r[1] = p[0] * q[1] + p[1] * q[0] + p[2] * q[3] - p[3] * q[2];
  r[2] = p[0] * q[2] - p[1] * q[3] + p[2] * q[0] + p[3] * q[1];
  r[3] = p[0] * q[3] + p[1] * q[2] - p[2] * q[1] + p[3] * q[0];
}
// generate rotation matrix r from quaternion q
void qua2rot(float r[], float q[])
{
  double x2 = q[1] * q[1] * 2.0;
  double y2 = q[2] * q[2] * 2.0;
  double z2 = q[3] * q[3] * 2.0;
  double xy = q[1] * q[2] * 2.0;
  double yz = q[2] * q[3] * 2.0;
  double zx = q[3] * q[1] * 2.0;
  double xw = q[1] * q[0] * 2.0;
  double yw = q[2] * q[0] * 2.0;
  double zw = q[3] * q[0] * 2.0;
  r[ 0] = 1.0 - y2 - z2;
  r[ 1] = xy + zw;
  r[ 2] = zx - yw;
  r[ 4] = xy - zw;
  r[ 5] = 1.0 - z2 - x2;
  r[ 6] = yz + xw;
  r[ 8] = zx + yw;
  r[ 9] = yz - xw;
  r[10] = 1.0 - x2 - y2;
  r[15] = 1.0;
  r[ 3] = 0.0;
  r[ 7] = 0.0;
  r[11] = 0.0;
  r[12] = 0.0;
  r[13] = 0.0;
  r[14] = 0.0;
}
// generate quaternion q from rotation matrix r
void rot2qua(float q[], float r[])
{
  float elem[4];
  elem[0] =  r[0] - r[5] - r[10] + 1.0;
  elem[1] = -r[0] + r[5] - r[10] + 1.0;
  elem[2] = -r[0] - r[5] + r[10] + 1.0;
  elem[3] =  r[0] + r[5] + r[10] + 1.0;
  unsigned biggestIndex = 0;
  for (int i=1; i<4; i++) {
    if (elem[i] > elem[biggestIndex])
      biggestIndex = i;
  }
  if (elem[biggestIndex] < 0.0)
    printf("rot2qua error\n");

  //float *qq[4] = {&q[0],&q[1],&q[2],&q[3]};
  float *qq[4] = {&q[1],&q[2],&q[3],&q[0]};
  float v = sqrtf(elem[biggestIndex])*0.5;
  *qq[biggestIndex] = v;
  float mult = 0.25 / v;
  
  switch (biggestIndex) {
  case 0:
    *qq[1] = (r[1] + r[4]) * mult;
    *qq[2] = (r[8] + r[2]) * mult;
    *qq[3] = (r[6] - r[9]) * mult;
    break;
  case 1:
    *qq[0] = (r[1] + r[4]) * mult;
    *qq[2] = (r[6] + r[9]) * mult;
    *qq[3] = (r[8] - r[2]) * mult;
    break;
  case 2:
    *qq[0] = (r[8] + r[2]) * mult;
    *qq[1] = (r[6] + r[9]) * mult;
    *qq[3] = (r[1] - r[4]) * mult;
    break;
  case 3:
    *qq[0] = (r[6] - r[9]) * mult;
    *qq[1] = (r[8] - r[2]) * mult;
    *qq[2] = (r[1] - r[4]) * mult;
    break;
  }
}

void remove_carriage_return(std::string& line)
{
    if (line.size() == 0) return;
    if (*line.rbegin() == '\r')
    {
        line.erase(line.length() - 1);
    }
}
