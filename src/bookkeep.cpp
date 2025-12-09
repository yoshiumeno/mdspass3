#include <iostream>
#include <stdlib.h>
#include <fstream>
#include "myheader.h"

void bond_set();
void bookkeep_alloc(int size);
void outproduct (double a[3], double b[3], double x[3]);
double veclength(double x[3]);
void bookkeep_write();
void bookkeep_read();
void inverse(double mat[3][3], double imat[3][3]);

void bookkeep()
{
 START:
  if (book.alloc) {
    bookkeep_alloc(book.nbook_new);
  }

  int ix, iy, iz, nrpx, nrpy, nrpz;
  double r2;
  if (rcut>0) { book.frc = ((double)rcut_f+frcmar)*1e-10; }
  //if (rcut>0) { book.frc = (rcut+frcmar)*1e-10; }
  book.frc2 = book.frc * book.frc;
  frc_f = book.frc*1e10;
  //printf("##Bookkeep is called: step = %d##\n",istep);
  //printf(" Cutoff radius = %f A,",rcut_f);
  //printf(" Bookkeep radius = %f A\n",frc_f);
  if (book.algo == 3) { return; }
  nrpx = 0; nrpy = 0; nrpz = 0;

  double vol = cell.hmat[0][0]*cell.hmat[1][1]*cell.hmat[2][2] + cell.hmat[0][1]*cell.hmat[1][2]*cell.hmat[2][0] 
    + cell.hmat[0][2]*cell.hmat[1][0]*cell.hmat[2][1] - cell.hmat[0][0]*cell.hmat[1][2]*cell.hmat[2][1] 
    - cell.hmat[0][1]*cell.hmat[1][0]*cell.hmat[2][2] - cell.hmat[0][2]*cell.hmat[1][1]*cell.hmat[2][0];

  double area;

  if (cell.pbcx)
  { 
	double b[3] = {cell.hmat[1][0], cell.hmat[1][1], cell.hmat[1][2]};
	double c[3] = {cell.hmat[2][0], cell.hmat[2][1], cell.hmat[2][2]};
	double a[3]; outproduct(b, c, a);
	area = veclength(a);
	nrpx=(int)(book.frc*2/(vol/area)+1);
  }

  if (cell.pbcy)
  { 
	double c[3] = {cell.hmat[2][0], cell.hmat[2][1], cell.hmat[2][2]};
	double a[3] = {cell.hmat[0][0], cell.hmat[0][1], cell.hmat[0][2]};
	double b[3]; outproduct(c, a, b);
	area = veclength(b);
	nrpy=(int)(book.frc*2/(vol/area)+1);
  }

  if (cell.pbcz) 
  { 
	double a[3] = {cell.hmat[0][0], cell.hmat[0][1], cell.hmat[0][2]};
	double b[3] = {cell.hmat[1][0], cell.hmat[1][1], cell.hmat[1][2]};
	double c[3]; outproduct(a, b, c);
	area = veclength(c);
	nrpz=(int)(book.frc*2/(vol/area)+1); 
  }

  //std::cout<<"Bookkeep: Range= "<<nrpx<<" x "<<nrpy<<" x "<<nrpz<<std::endl;
  cell.Setlen(); // calculate length of cell edges
  inverse(cell.hmat, cell.hinmat); // calculate inverse cell matrix
  atom.Calcq();
  // Check if cell is orthogonal
  bool ortho = true;
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      if ((i != j)&&(cell.hmat[i][j] != 0.0)) { ortho = false; } 
    }
  }

  //  std::cout << ortho << nrpx << nrpy << nrpz;
  if (ortho&&(nrpx<=1)&&(nrpy<=1)&&(nrpz<=1)) { // large & orthogonal cell
    for (int i=1; i<=atom.natom; i++) {
      book.alistnum[i] = 0;
      for (int j=1; j<=atom.natom; j++) {
	double d2 = atom.Dist2Closest_ortho(i,j,ix,iy,iz);
	if ((d2<book.frc2)&&(d2>1.0e-30)) {
	  book.alistnum[i]++;
	  book.alist[i][book.alistnum[i]][0] = j;
	  book.alist[i][book.alistnum[i]][1] = ix;
	  book.alist[i][book.alistnum[i]][2] = iy;
	  book.alist[i][book.alistnum[i]][3] = iz;
	}
      }
    }
  } else if ((!ortho)&&(nrpx<=1)&&(nrpy<=1)&&(nrpz<=1)) { // large & non-orthogonal cell
    for (int i=1; i<=atom.natom; i++) {
      book.alistnum[i] = 0;
      for (int j=1; j<=atom.natom; j++) {
	double d2 = atom.Dist2Closest(i,j,ix,iy,iz);
	if ((d2<book.frc2)&&(d2>1.0e-30)) {
	  book.alistnum[i]++;
	  book.alist[i][book.alistnum[i]][0] = j;
	  book.alist[i][book.alistnum[i]][1] = ix;
	  book.alist[i][book.alistnum[i]][2] = iy;
	  book.alist[i][book.alistnum[i]][3] = iz;
	}
      }
    }
  } else if (ortho&&(nrpx>1)&&(nrpy<=1)&&(nrpz<=1)) { // thin in x and orthogonal
    for (int i=1; i<=atom.natom; i++) {
      book.alistnum[i] = 0;
      for (int j=1; j<=atom.natom; j++) {
	//if (atom.Dist2Closest_yz_ortho(i,j,0,iy,iz)<book.frc2) { // <-- more boost, but risky?
	for (int ixx=-nrpx; ixx<=nrpx; ixx++) {
	  double val = atom.Dist2Closest_yz_ortho(i,j,ixx,iy,iz);
	  if ((val<book.frc2)&&(val>1.0e-30)) {
	    book.alistnum[i]++;
	    book.alist[i][book.alistnum[i]][0] = j;
	    book.alist[i][book.alistnum[i]][1] = ixx;
	    book.alist[i][book.alistnum[i]][2] = iy;
	    book.alist[i][book.alistnum[i]][3] = iz;
	  }
	}
	//}
      }
    }
  } else if ((!ortho)&&(nrpx>1)&&(nrpy<=1)&&(nrpz<=1)) { // thin in x and non-orthogonal
    for (int i=1; i<=atom.natom; i++) {
      book.alistnum[i] = 0;
      for (int j=1; j<=atom.natom; j++) {
	//if (atom.Dist2Closest_yz(i,j,0,iy,iz)<book.frc2) { // <-- more boost, but risky?
	for (int ixx=-nrpx; ixx<=nrpx; ixx++) {
	  double val = atom.Dist2Closest_yz(i,j,ixx,iy,iz);
	  if ((val<book.frc2)&&(val>1.0e-30)) {
	    book.alistnum[i]++;
	    book.alist[i][book.alistnum[i]][0] = j;
	    book.alist[i][book.alistnum[i]][1] = ixx;
	    book.alist[i][book.alistnum[i]][2] = iy;
	    book.alist[i][book.alistnum[i]][3] = iz;
	  }
	}
	//}
      }
    }
  } else if (ortho&&(nrpx<=1)&&(nrpy>1)&&(nrpz<=1)) { // thin in y and orthogonal
    for (int i=1; i<=atom.natom; i++) {
      book.alistnum[i] = 0;
      for (int j=1; j<=atom.natom; j++) {
	//if (atom.Dist2Closest_xz_ortho(i,j,ix,0,iz)<book.frc2) { // <-- more boost, but risky?
	for (int iyy=-nrpy; iyy<=nrpy; iyy++) {
	  double val = atom.Dist2Closest_xz_ortho(i,j,ix,iyy,iz);
	  if ((val<book.frc2)&&(val>1.0e-30)) {
	    book.alistnum[i]++;
	    book.alist[i][book.alistnum[i]][0] = j;
	    book.alist[i][book.alistnum[i]][1] = ix;
	    book.alist[i][book.alistnum[i]][2] = iyy;
	    book.alist[i][book.alistnum[i]][3] = iz;
	  }
	}
	//}
      }
    }
  } else if ((!ortho)&&(nrpx<=1)&&(nrpy>1)&&(nrpz<=1)) { // thin in y and non-orthogonal
    for (int i=1; i<=atom.natom; i++) {
      book.alistnum[i] = 0;
      for (int j=1; j<=atom.natom; j++) {
	//if (atom.Dist2Closest_xz(i,j,ix,0,iz)<book.frc2) { // <-- more boost, but risky?
	for (int iyy=-nrpy; iyy<=nrpy; iyy++) {
	  double val = atom.Dist2Closest_xz(i,j,ix,iyy,iz);
	  if ((val<book.frc2)&&(val>1.0e-30)) {
	    book.alistnum[i]++;
	    book.alist[i][book.alistnum[i]][0] = j;
	    book.alist[i][book.alistnum[i]][1] = ix;
	    book.alist[i][book.alistnum[i]][2] = iyy;
	    book.alist[i][book.alistnum[i]][3] = iz;
	  }
	}
	//}
      }
    }
  } else if (ortho&&(nrpx<=1)&&(nrpy<=1)&&(nrpz>1)) { // thin in z and orthogonal
    for (int i=1; i<=atom.natom; i++) {
      book.alistnum[i] = 0;
      for (int j=1; j<=atom.natom; j++) {
	//if (atom.Dist2Closest_xy_ortho(i,j,ix,iy,0)<book.frc2) { // <-- more boost, but risky?
	for (int izz=-nrpz; izz<=nrpz; izz++) {
	  double val = atom.Dist2Closest_xy_ortho(i,j,ix,iy,izz);
	  if ((val<book.frc2)&&(val>1.0e-30)) {
	    book.alistnum[i]++;
	    book.alist[i][book.alistnum[i]][0] = j;
	    book.alist[i][book.alistnum[i]][1] = ix;
	    book.alist[i][book.alistnum[i]][2] = iy;
	    book.alist[i][book.alistnum[i]][3] = izz;
	  }
	}
	//}
      }
    }
  } else if ((!ortho)&&(nrpx<=1)&&(nrpy<=1)&&(nrpz>1)) { // thin in z and non-orthogonal
    for (int i=1; i<=atom.natom; i++) {
      book.alistnum[i] = 0;
      for (int j=1; j<=atom.natom; j++) {
	//if (atom.Dist2Closest_xy(i,j,ix,iy,0)<book.frc2) { // <-- more boost, but risky?
	for (int izz=-nrpz; izz<=nrpz; izz++) {
	  double val = atom.Dist2Closest_xy(i,j,ix,iy,izz);
	  if ((val<book.frc2)&&(val>1.0e-30)) {
	    book.alistnum[i]++;
	    book.alist[i][book.alistnum[i]][0] = j;
	    book.alist[i][book.alistnum[i]][1] = ix;
	    book.alist[i][book.alistnum[i]][2] = iy;
	    book.alist[i][book.alistnum[i]][3] = izz;
	  }
	}
	//}
      }
    }
    // When cell is small block (short in 3 directions)
  } else {
    for (int i=1; i<=atom.natom; i++) {
      book.alistnum[i] = 0;
      for (int j=1; j<=atom.natom; j++) {
	for (int ixx=-nrpx; ixx<=nrpx; ixx++) {
	  for (int iyy=-nrpy; iyy<=nrpy; iyy++) {
	    for (int izz=-nrpz; izz<=nrpz; izz++) {
	      if (atom.Dist2(i,j,ixx,iyy,izz)<book.frc2) {
		if (atom.Dist2(i,j,ixx,iyy,izz)>1.0e-30) {
		  book.alistnum[i]++;
		  if (book.alistnum[i]>=book.nbook) { goto OUT; }
		  book.alist[i][book.alistnum[i]][0] = j;
		  book.alist[i][book.alistnum[i]][1] = ixx;
		  book.alist[i][book.alistnum[i]][2] = iyy;
		  book.alist[i][book.alistnum[i]][3] = izz;
		}
	      }
	    }
	  }
	}
      }
    }
  }
   
 OUT:
  //if (book.alistnum[i]>=NBOOK) {
  for (int i=1; i<=atom.natom; i++) {
    if (book.alistnum[i]>=book.nbook) {
      //std::cout<<"Bookkeep Error: NBOOK too small"<<std::endl;
      book.nbook_new = book.nbook*1.5;
      book.alloc = true;
      std::cout<<"Bookkeep: Reallocate arrays (from ";
      std::cout<<book.nbook<<" to "<<book.nbook_new<<")"<<std::endl;
      goto START;
    }
    if (book.alloc) goto START;
  }
  bond_set();
}

void bookkeep_alloc(int size)
{
  //alist[NMAX+1][NBOOK+1][4];
  if (book.alist) {
    for (int i=0; i<=book.natom; i++) {
      for (int j=0; j<=book.nbook; j++) {
	delete[] book.alist[i][j]; }
      delete[] book.alist[i]; }
    delete[] book.alist; book.alist = NULL;
  }
  book.natom = atom.natom;
  book.alist = new int**[atom.natom+1];
  for (int i=0; i<=atom.natom; i++) {
    book.alist[i] = new int*[size+1];
    for (int j=0; j<=size; j++) {
      book.alist[i][j] = new int[4];
    }
  }
  book.nbook = size;
  book.alloc = false;
}

void bookkeep_write()
{
  printf("## book-keeping table is written to BOOKKEEP_TABLE\n");
  FILE* fp = fopen("BOOKKEEP_TABLE", "w");
  for (int i=1; i<=atom.natom; i++) {
    fprintf(fp, "%d %d\n", i, book.alistnum[i]);
    for (int j=1; j<=book.alistnum[i]; j++) {
      int alist0 = book.alist[i][j][0];
      int alist1 = book.alist[i][j][1];
      int alist2 = book.alist[i][j][2];
      int alist3 = book.alist[i][j][3];
      fprintf(fp, "%d %d %d %d\n", alist0, alist1, alist2, alist3);
    }
  }
  fclose(fp);
}
void bookkeep_read()
{
  // return value:
  // 0 = success
  // 1 = data structure is wrong
  // 2 = book.nbook is too small
  printf("## book-keeping table is read from BOOKKEEP_TABLE\n");
  int i0, i1;
  int alist0, alist1, alist2, alist3;
  FILE* fp = fopen("BOOKKEEP_TABLE", "r");
 READ:
  for (int i=1; i<=atom.natom; i++) {
    fscanf(fp, "%d %d\n", &i0, &i1);
    if (i0 != i) {
      printf("### ERROR! ###\n");
      printf("  error in bookkeep_read error\n");
      return;
    }
    book.alistnum[i] = i1;
    if (book.alistnum[i]>=book.nbook) {
      printf(" bookkeep_read: book.nbook is expanded\n");
      book.nbook_new = book.nbook*1.5;
      book.alloc = true; bookkeep_alloc(book.nbook_new);
      goto READ;
    }
    for (int j=1; j<=book.alistnum[i]; j++) {
      fscanf(fp, "%d %d %d %d\n", &alist0, &alist1, &alist2, &alist3);
      book.alist[i][j][0] = alist0;
      book.alist[i][j][1] = alist1;
      book.alist[i][j][2] = alist2;
      book.alist[i][j][3] = alist3;
    }
  }
  fclose(fp);
  return;
}

/*    
int search_neighbor(i, j, &ix, &iy, &iz)
{
  double sdx, sdy, sdz;
  sdx = atom.rx[i] - atom.rx[j];
  sdy = atom.ry[i] - atom.ry[j];
  sdz = atom.rz[i] - atom.rz[j];
}
*/

