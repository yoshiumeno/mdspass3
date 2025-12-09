#ifndef _MYHEADER_H
#define _MYHEADER_H


#include<math.h>
#include<string.h>

#include <FL/Fl.H>
#include <FL/Fl_Widget.H>

#define NMAX 2500
#define NBOOK 1000
#define ang 1.0e-10
#define eV  1.6021892e-19
#define GPa 1.0e9
#define MPa 1.0e6
#define SELECTIONS 1000000
#define MAXPOTTYPE 14
#define MAXPOTARGTYPE 20
#define MAXENSTYPE 11
#define MAXNEIGHBOR 100
#define NINTEGRAL 10
#define MAXMODE 30
#define MAXKP 5
#define MAXCOLOR 8
#define MAXGROUP 10

#define FLOAT double
#define FLOAT2 double

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif
#include "output.h"
#include "stringutility.h"

#include "pair_function.h" // Kubo 20140224
#include "pair_function2.h"

#include "ReaxPotential.h"
#include "Configuration.h"
#include "ReaxFF.h"

// 2024/02 CJS 追加
// コールバック関数のプロトタイプ宣言です。
/**
* 主に、GUI のボタンが押された時に呼び出されるコールバック関数です。
*
* @param widget ボタンのインスタンス
* @param process_id プロセス ID
**/
//void pointer_cb(Fl_Widget* widget, void* process_id);
void pointer_cb(Fl_Widget* widget, long process_id);

/**
* 主に、OpenGL の描画が行われる時に呼び出されるコールバック関数です。
*
* @param widget ボタンのインスタンス
* @param process_id プロセス ID
**/
//void control_cb(Fl_Widget* widget, void* process_id);
void control_cb(Fl_Widget* widget, long process_id);
// 2024/02 CJS 追加 ここまで

class Atom
{
public:
    Atom();
    int natom, nrepatom;
    char** asp; //Atom species
    int* anum; //Atomic number
    char potential_func[30];
    char potential_arg[30];
    char potstring_list[MAXPOTTYPE][30];
    char potargstring_list[MAXPOTTYPE][MAXPOTARGTYPE][30];
    int potarg_number[MAXPOTTYPE];
    bool potarg_readable[MAXPOTTYPE];
    double* wm;
    double* rx, * ry, * rz, * fx, * fy, * fz;
    double* rx_org, * ry_org, * rz_org, * rx_p, * ry_p, * rz_p;
    double* vx, * vy, * vz;
    double* qx, * qy, * qz;
    double* ax, * ay, * az, * bx, * by, * bz, * cx, * cy, * cz; // for p-c integral
    double* fx_l, * fy_l, * fz_l; //loading
    int* group; //group of atoms (YU2025)
    double epotsum, ekinsum;
    double* epot, * epot_p;
    FLOAT* rx_float, * ry_float, * rz_float, * fx_float, * fy_float, * fz_float;
    FLOAT2* epot_float;
    bool* mfx, * mfy, * mfz; // To fix atom motion
    bool* lock; // To lock atom (for gloc-accel)
    int nelem; // number of elements
    int* repatom; // whether an atom is repatom or not (1=yes, 0=no)
    int* elem_id; // element number to which an atom belongs
    int** elem_v, ** elem_v_rep;
    int nwall; // number of wall-triangle (for CNT)
    int** wall_v; // vertices of wall (for CNT)  ==> wall_v[nwall][3]
    double* wall_area; // area of wall (for CNT) ==> wall_area[nwall]
    int* wall_vrev; // reverse normal vector or not
    double** wall_nvec; // normal vector of wall (for CNT) ==> wall_nvec[nwall][3]
    int** neighbor; // list of neiboring atoms
    int* nneighbor; // number of neiboring atoms
    bool* visible; // whether atom is drawn or not
    double Enkin();
    double Temp();
    double MSD();
    double Mises(int i);
    double Dist2(int i, int j);
    double Dist2(int i, int j, int ix, int iy, int iz);
    double Dist2Closest(int i, int j, int& ix, int& iy, int& iz);
    double Dist2Closest_yz(int i, int j, int ix, int& iy, int& iz);
    double Dist2Closest_xz(int i, int j, int& ix, int iy, int& iz);
    double Dist2Closest_xy(int i, int j, int& ix, int& iy, int iz);
    //double Dist2Closest_yz(int i, int j, int &iy, int &iz);
    //double Dist2Closest_xz(int i, int j, int &ix, int &iz);
    //double Dist2Closest_xy(int i, int j, int &ix, int &iy);
    double Dist2Closest_ortho(int i, int j, int& ix, int& iy, int& iz);
    double Dist2Closest_yz_ortho(int i, int j, int ix, int& iy, int& iz);
    double Dist2Closest_xz_ortho(int i, int j, int& ix, int iy, int& iz);
    double Dist2Closest_xy_ortho(int i, int j, int& ix, int& iy, int iz);
    double Dx(int i, int j, int ix, int iy, int iz);
    double Dy(int i, int j, int ix, int iy, int iz);
    double Dz(int i, int j, int ix, int iy, int iz);
    double Angle(int i, int j, int k, int ix0, int iy0, int iz0, int ix1, int iy1, int iz1);
    double Fmax();
    int QC;
    int instcenter;
    double** evecx, ** evecy, ** evecz;
    double eigval[MAXMODE];
    double*** satom;
    int getAtomNumber(char* atom)
    {
        if (( strcmp(atom, "C") == 0 ) || ( strcmp(atom, "c") == 0 ))
        {
            return 6;
        }
        else if (( strcmp(atom, "Cu") == 0 ) || ( strcmp(atom, "cu") == 0 ))
        {
            return 29;
        }
        else if (( strcmp(atom, "Al") == 0 ) || ( strcmp(atom, "al") == 0 ))
        {
            return 13;
        }
        else if (( strcmp(atom, "Si") == 0 ) || ( strcmp(atom, "si") == 0 ))
        {
            return 14;
        }
        else if (( strcmp(atom, "Ni") == 0 ) || ( strcmp(atom, "ni") == 0 ))
        {
            return 22;
        }
        else if (( strcmp(atom, "Fe") == 0 ) || ( strcmp(atom, "fe") == 0 ))
        {
            return 26;
        }
        else if (( strcmp(atom, "Ti") == 0 ) || ( strcmp(atom, "ti") == 0 ))
        {
            return 28;
        }
        else if (( strcmp(atom, "Ge") == 0 ) || ( strcmp(atom, "ge") == 0 ))
        {
            return 32;
        }
        else if (( strcmp(atom, "Ag") == 0 ) || ( strcmp(atom, "ag") == 0 ))
        {
            return 47;
        }
        else if (( strcmp(atom, "Sn") == 0 ) || ( strcmp(atom, "sn") == 0 ))
        {
            return 50;
        }
        else if (( strcmp(atom, "Au") == 0 ) || ( strcmp(atom, "au") == 0 ))
        {
            return 79;
        }
        else if (( strcmp(atom, "Pd") == 0 ) || ( strcmp(atom, "pd") == 0 ))
        {
            return 46;
        }
        else if (( strcmp(atom, "Pt") == 0 ) || ( strcmp(atom, "pt") == 0 ))
        {
            return 78;
        }
        else { return 1; }
    }
    void Calcq();
    void Calcq(int i);
private:
    double dist2;
};

class Book
{
public:
    Book();
    int* alistnum;
    int*** alist;
    double frc, frc2;
    int nbk;
    int algo;
    int natom, nbook, nbook_new;
    bool alloc;
};

class Cell
{
public:
    Cell();
    double hmat[3][3], hinmat[3][3], hvmat[3][3];
    double hamat[3][3], hbmat[3][3], hcmat[3][3];
    double hmat_org[3][3];
    double len[3];
    int pbcx, pbcy, pbcz;
    double alat, clat;
    double volume;
    double virx, viry, virz;
    double sgmmat[3][3], dmat[3][3], sgmmat_set[3][3];
    double sgmmat_p[3][3], hmat_p[3][3];
    bool relax_static_initial;
    double ecmat[6][6];
    double Getvolume();
    void Reset();
    void Setlen();
    double Dstress();
    double ww, prmass, damper_param;
    double wq, nhmass, s, sv, sa, sb, sc, sf;
    bool fix[3][3];
};

class Tersoff
{
public:
    Tersoff();
    //Tersoff() {
    //  initialize = true;
    //  nocutoff = false;
    //};
    double** zmat, ** b;
    double** z_ij, ** z_ji, ** b_ij, ** b_ji; int maxnei;
    double* terrr, * terss, * teraa, * terbb, * terlambda, * termu, * terbeta;
    double* tern, * terc2, * terd2, * terh, * terw, * termum, * terchi;
    int* term;
    //double terc2, terd2, termum;
    bool initialize;
    int ntype, nptype, * type, ** ptype, * typen;
    int natom, functype;
    bool nocutoff; // no-cutoff scheme (removal of cutoff after 2nd step)
    bool large; // large system
    bool debug_large; // debug for large system mode
};

class GEAM
{
public:
    GEAM();
    double gre[16], gfe[16], grhoe[16], grhos[16], galp[16], gbet[16], gaa[16], gbb[16], gkai[16],
        glam[16], gffn0[16], gffn1[16], gffn2[16], gffn3[16], gff0[16], gff1[16], gff2[16], gff3[16],
        geta[16], gffe[16];
    int* sp;
    double* rhob;
    bool initialize;
};

class Eammis
{
public:
    Eammis()
    {
        initialize = true;
        mesh = 10000;
        //mesh = 1000000;
    };
    double* pair, * paird, * den, * dend, * embed, * embedd;
    double pairmin, pairmax, denmin, denmax, embedmin, embedmax;
    double* pairx, * pairy, ** pairc, * pairxn;
    double* denx, * deny, ** denc, * denxn;
    double* embedx, * embedy, ** embedc, * embedxn;
    double* rhob;
    double* aden, * bden, * cden, * dden;
    double* apair, * bpair, * cpair, * dpair;
    double* aembed, * bembed, * cembed, * dembed;
    int mesh;
    bool initialize;
    char species[3];
    int natoms; // number of atoms, that eam is applied to (needed for combined mode)
    int atom_number; // number of chemical element eam is applied to (needed to avoid string comparison)
};

class Adp
{
public:
    Adp()
    { // Constructor
        initialize = true; cut = 7.0;
#if defined __linux__ || defined __APPLE__
        strcpy(fname, "pot/ADP_NdFeB.pot");
#else
        strcpy(fname, "pot\\ADP_NdFeB.pot");
#endif
    };
    // Potential parameters should come here
    double cut;
    bool initialize;
    double* rho, * mu_x, * mu_y, * mu_z, * nu;
    double* lambda_xx, * lambda_yy, * lambda_zz;
    double* lambda_xy, * lambda_yz, * lambda_zx;
    //double *param_pair, *param_dp, *param_qp;
    //double *param_den, *param_eb, *gradF;
    double** param_pair, ** param_dp, ** param_qp;
    double** param_den, ** param_eb, * gradF;
    int ntype, nptype;
    int* type, ** ptype, * typen;
    char fname[120];
    // Kubo 20140224 -------------//
    PairFunction** Phi;
    PairFunction** Rho;
    PairFunction** F;
    PairFunction** U;
    PairFunction** W;
    // End -----------------------//
};

class Dipole
{
public:
    Dipole()
    {
        dp_eps = 14.40; dp_cut = 8.0; dp_tol = 1.e-7; dp_mix = 0.2;
        // ntype = 3; nptype = 3;
        initialize = true;
        shift = true;
        //    shift = false; // Kubo
#if defined __linux__ || defined __APPLE__
        strcpy(fname, "pot/Dipole_YSZ.pot");
#else
        strcpy(fname, "pot\\Dipole_YSZ.pot");
#endif
    };
    double* dp_alpha, * dp_b, * dp_c, * r_cut;
    double* E_statx, * E_staty, * E_statz;
    double* E_indx, * E_indy, * E_indz;
    double* E_oldx, * E_oldy, * E_oldz;
    double* E_totx, * E_toty, * E_totz;
    double* p_srx, * p_sry, * p_srz;
    double* p_indx, * p_indy, * p_indz;
    double dp_eps, dp_cut, dp_tol, dp_mix, ew_rcut;
    double* ratio, * charge, last_charge, dp_kappa;
    double* buck_a, * buck_s, * buck_c;
    double* buck_vrc, * buck_vprc;
    int sw_kappa;
    bool initialize, shift;
    int ntype, nptype;
    int* type, ** ptype, * typen;
    char fname[120];
    double val1, val2, val3;
    int natoms; // for combined mode
    //int mode; // for combined mode
};

class Combined
{
public:
    Combined()
    {
        initialize = true;
    }

    bool initialize;
};

class Bre
{
public:
    Bre()
    {
        //    std::cout<< "Bre constructor was called" << std::endl;
        initialize = true;
        npmax = 7000; nlmax = 1000000, ntab = 10000;
        npm1 = npmax + 1; ntypes = 8; nnma = 5000 * 40;
        ikl = 2; nmabig = nnma * 5; ver = 3.00;
        xqm = 3.70; att = 3.20;
        ljmin = 2.0; ljmax = 10.0; ljmesh = 10000;
        //    std::cout<<nmabig<<std::endl;
    };
    int npmax, nlmax, ntab, npm1, ntypes, nnma, ikl,
        nmabig;
    double ver;
    int igh[26], in3[65][4], ndihed;
    double* bww, * dww, * rcor, * cor1, * cor2, * cor3;
    int* rep1, * rep2, * rep3;
    double xh[3][11][11], xh1[3][11][11], xh2[3][11][11];
    int* list, * nabors, * lcheck;
    double ad[5][5], axl[5][5], bd[5][5], bxl[5][5],
        cd[5][5], cxl[5][5], dd[5][5], dxl[5][5], ed[5][5],
        xn1[3], rb1[5][5], rb2[5][5], pid[5][5], rmax[5][5],
        rlist[5][5], spgc[7][6], spgh[7][4], xq, att, xqm, pq,
        xtn2[5], xtn1[5], adb[5], cdb[5], cd2[5], ddb[5],
        ddb2[5], hdb[5], chi[5][5], xtm[5], xtl[5];
    double epslj[5][5], siglj[5][5], epsts[5][5], bmin[5][5], bmax[5][5]; //AIREBO
    int igc[26];
    double xxdb, xdb[5][5][5], reg[5][5][5];
    double clmn[4][11][11][11][65], clm[3][11][11][17],
        tlmn[11][11][11][65], pidt;
    int in2[17][3];
    double* exx1, * dexx1, * exx2;
    int* ivct2b, * jvct2b;
    double* xhc1, * xhc2;
    double*** atable, *** datable, *** rtable, *** drtable,
        ddtab[5][5], *** tabfc, *** tabdfc;
    bool initialize;

    double** eps, ** sig, * xmass;
    int* noa;
    double rhh, rch, sigma, epsi, rll;
    double rc, rlis;
    int* ktype, kt[101], kt2[101], np;
    double* r01, * r02, * r03, * rnp1, * rnp2, * rnp3;
    double* rpp1, * rpp2, * rpp3;
    double tote;
    int kend;
    //AIREBO
    double* lj, * ljd, ljmin, ljmax;
    int ljmesh;
};

//// Kubo 20140618 ---------------------------------------------
class Sw
{
public:
    Sw(void); // Constructor
    int nSpec; // Num of Species
    int nPair; // Num of Pair
    int nTrio; // Num of Trio (Triplet)
    int AtomNum2Index[200]; // Mapping: AtomicNumber -> SpecIndex (0~nSpec).
    // 0 means undefined (= unused) interaction.
    int** PairType; // nSpec x nSpec -> PairIndex (0~nPair).
    int*** TrioType; // nSpec x nSpec x nSpec -> TrioIndex (0~nTrio).
    PairFunction** Phi; // Pair Interaction
    PairFunction** R;   // Pair Interaction in Trio Term
    PairFunction** P;   // Angular-Dependent Function
    bool initialize;
    void DeleteArray(void);
    void AllocateArray(void);
    void DeletePhi(void);
    void DeleteR(void);
    void DeleteP(void);
    void DeletePairType(void);
    void DeleteTrioType(void);
    void AllocatePairType(void);
    void AllocateTrioType(void);
    void AllocatePhi(void);
    void AllocateR(void);
    void AllocateP(void);
};
//// End -------------------------------------------------------

/*
class Sw : public Vashishta{
 public:
  Sw(void); // Constructor
};
*/

class Cond
{
public:
    Cond();
};

class Fire
{
public:
    Fire();
    double sp, fnorm, vnorm, ffinc, ffdec, ffalph, fire_alph, fire_alph_ini, ffdtmax;
    int ifire, nfmin;
};

class Cg
{
public:
    Cg();
    double lamtr, lam_init;
    double** sgx, ** sgy, ** sgz, ** shx, ** shy, ** shz;
    double* rstx, * rsty, * rstz;
    double pot0, pot1, grad0, grad1, normsh;
    int step;
};

class Dmc2d
{
public:
    int pattern_x, pattern_y;
    int** pattern;
    double coef[2];
};

class MyReax
{
public:
    MyReax()
    {
#if defined __linux__ || defined __APPLE__
        strcpy(potfname, "pot/ReaxFF.pot");
#else
        strcpy(potfname, "pot\\ReaxFF.pot");
#endif
    };
    bool initialize;
    ReaxPotential myreaxpot;
    Input myreaxinput;
    ConfigurationSet myreaxconf;
    char potfname[120];
};

class Shell
{
public:
    double alarge[4], rho[4], clarge[4], rcsr[4], rtp[4];
    double q[3][2], springc[3], springc2[3];
    double efpe;
    int kint[3][3];
    double epsil0 = 8.854187816e-12; // Permittivity in vacuum (F/m)
    double elc = 1.6021892e-19;  // elemetary charge e (C)
    int naccew = 10; // Accuracy parameter for Ewald method
    bool initialize;
    double** tmpfx, ** tmpfy, ** tmpfz;
    double** tmprx, ** tmpry, ** tmprz;
    Shell()
    {
        initialize = true;
        //-----Shell model parameters for PbTiO3
        // Pb-0, O-1, Ti-2  (ksp) 
        // Interaction number
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                kint[i][j] = 0;
            }            
        }
        kint[0][1] = 0; kint[1][0] = 0; //Pb-O
        kint[1][1] = 1; //O-O
        kint[2][1] = 2; kint[1][2] = 2; //Ti-O
        kint[0][2] = 3; kint[2][0] = 3; //Pb-Ti      
        // Charge (e)
        q[0][0] =  5.495849e0;        // Pb Core
        q[0][1] = -3.633216e0;        // Pb Shell
        q[1][0] =  2.548431e0;        // O  Core
        q[1][1] = -4.199168e0;        // O  Shell
        q[2][0] = 19.369090e0;        // Ti Core
        q[2][1] =-16.279518e0;        // Ti Shell 
        // Spring constant k2, k4 (eV/Angst**2,eV/Angst**4)
        springc[0] =  154.87128e0;  // Pb
        springc[1] =  180.91341e0;  // O
        springc[2] = 8829.4096e0;   // Ti
        springc2[0]= 22416.670e0;   // Pb
        springc2[1]=  6945.7814e0;  // O
        springc2[2]=1928581.7e0;    // Ti
        // Short-range potential parameters
        // Aij (eV)
        alarge[0]=2538.4110e0;     // Pb-O
        alarge[1]=1698.6653e0;     // O -O
        alarge[2]=2555.2075e0;     // Ti-O
        alarge[3]=387.73157e0;     // Pb-Ti
        // rhoij (Angst)
        rho[0]=0.300698e0;         // Pb-O
        rho[1]=0.271756e0;         // O -O
        rho[2]=0.278391e0;         // Ti-O
        rho[3]=0.394957e0;         // Pb-Ti
        // Cij (eV*Angst**6)
        clarge[0]=2.6167608e0;     // Pb-O
        clarge[1]=61.843537e0;     // O -O
        clarge[2]=2.2555672e0;     // Ti-O
        clarge[3]=223.24409e0;     // Pb-Ti
        // Cutoff radius (Angst)
        rcsr[0]=8.0e0;             // Pb-O
        rcsr[1]=8.0e0;             // O -O
        rcsr[2]=8.0e0;             // Ti-O
        rcsr[3]=8.0e0;             // Pb-Ti
        // Taper range (Angst)
        rtp[0]=2.0e0;              // Pb-O
        rtp[1]=2.0e0;              // O -O
        rtp[2]=2.0e0;              // Ti-O
        rtp[3]=2.0e0;              // Pb-Ti
        //-Parameters for Ewald method  -Molecular Simulation, vol.1, pp.207-224 (1988)-
        efpe = elc*elc/(4.0*M_PI*epsil0);
    };
};

//Ewald method
class Ew
{
public:
    double q[3];
    double efpe;
    double epsil0 = 8.854187816e-12; // Permittivity in vacuum (F/m)
    double elc = 1.6021892e-19;  // elemetary charge e (C)
    int naccew = 10; // Accuracy parameter for Ewald method
    bool initialize;
    double rcsr = 8.0; // Cutoff radius (Angst)
    Ew()
    {
        initialize = true;
        // Charge (e)
        q[0] =  2.0; //Pb
        q[1] = -2.0; //O
        q[2] =  3.0; //Ti
        //-Parameters for Ewald method  -Molecular Simulation, vol.1, pp.207-224 (1988)-
        efpe = elc*elc/(4.0*M_PI*epsil0);
    };
};

extern int istep;
extern int mdmotion;
extern double rcut, rcut2;
extern float rcut_f, frcmar, frc_f;
extern float temp_set;
extern float tempc, cellx, celly, cellz, f_max, epotatom, dtm;
extern float cell1x, cell1y, cell1z, cell2x, cell2y, cell2z, cell3x, cell3y, cell3z;
extern float slice_1min, slice_1max, slice_2min, slice_2max, slice_3min, slice_3max;
extern int trim_mode, trim_cylinder_axis;
extern float trim_cylinder_diameter, trim_cylinder_diameter2;
extern int ifix_atoms;
extern float strs_xx, strs_yy, strs_zz, strs_xy, strs_yz, strs_zx;
extern int ensemble, ensemble_p;
extern int config_type;
extern int irepx, irepy, irepz, icntm, icntn;
extern float alat, clat, cscnt, rotz, shiftz;
extern int imerge;
extern char current_config_name[100], current_qcelement_name[100];
extern const char* potstring_list[];
extern int ipottype, ipotarg;
extern int notrans, incell;
extern FILE* energyfile, * stressfile, * cellfile, * ssnorfile, * ecnorfile, * ecallfile;
extern FILE* msdfile, * miscfile, * tempfile;
extern double msd;
extern double dt;
extern int repeat_lz;
extern float dexdt, deydt, dezdt, repeat_lz_min, repeat_lz_max;
extern int ifgrab1, ifgrab2, grab1num, grab2num;
extern int ifpush1, ifpush2, push1num, push2num;
extern float grabdxdt1, grabdydt1, grabdzdt1, grabdxdt2, grabdydt2, grabdzdt2;
extern float pushfx1, pushfy1, pushfz1, pushfx2, pushfy2, pushfz2;
extern int yzplane_punch;
extern float yzplane_punch_d, yzplane_punch_dd, yzplane_punch_ftot;
//extern float ex;
//extern float dexdt;
extern void ( *integral_a[NINTEGRAL] )( ), ( *integral_b[NINTEGRAL] )( );
extern int integral_type;
extern float prmass_scale, nhmass_scale;
extern int initial_velocity;

extern Atom atom;
extern Cell cell;
extern Book book;
extern Tersoff tersoff;
extern GEAM geam;
extern Eammis eammis;
extern Dipole dipole;
extern Cond cond;
extern Bre bre;
extern Fire fire;
extern Adp adp;
extern Cg cg;
extern Combined combined;
extern Pair* pairPot;
extern Dmc2d dmc2d;
//extern Vashishta vashishta;
extern Sw sw;
extern MyReax myreax;
extern Shell shell;
extern Ew ew;


extern int istep0;
extern int mdspeed;
extern int confwrint, autocap;
extern int itolfor; extern float tolfor;
extern int ireheat;
extern float reheat_temp, reheat_dtm;
extern int reheat_step, reheat_count;
extern int itolstep; extern int tolstep;
extern int itolstress; extern float tolstress;
extern int irecipe;
extern int inogui;
extern float rigid_shift_z, rigid_shift_dx, rigid_shift_dy, rigid_shift_dz;
extern int rigid_shift_xtimes, rigid_shift_ytimes;
extern float strs_set_xx, strs_set_yy, strs_set_zz, strs_set_xy, strs_set_yz, strs_set_zx;
extern int cellfix_xx, cellfix_yy, cellfix_zz, cellfix_xy, cellfix_yz, cellfix_zx;
extern float fp_alph_ini, fp_ffinc, fp_ffdec, fp_ffalph, fp_ffdtmax;
extern int fp_nfmin;
extern float auto_val1, auto_val2;
extern int group_num, group_range_axis, group_num_show, if_group_num_show;
extern float group_range_max, group_range_min;

// for Rotation/Motion by Mouse
extern int cx, cy;
extern float cs, cpos0, cpos1;
extern double sx, sy;
extern float cq[4]; // Initial quaternion
extern float tq[4]; // quaternion while mouse drag
extern bool right_button;

#define WRITE_SETDAT_ID 1102
#define RESET_ID 102
#define WRITECONFIG_ID 1061
#define CREATECONFIG_ID 1062
#define CREATECONFIG_DO_ID 1064
#define MULTIPLYCELL_DO_ID 1069
#define REMOVEATOM_DO_ID 1070
#define REMOVEPART_DO_ID 1082
#define SETGROUP_DO_ID 2001 //YU2025
#define RESETGROUP_DO_ID 2002 //YU2025
#define ADDATOM_DO_ID 1071
#define MOVEATOM_DO_ID 1072
#define SHIFTALL_DO_ID 1073
#define SHIFTONE_DO_ID 1081
#define SHIFTPART_DO_ID 1078
#define ROTATEALL_DO_ID 1074
#define ROTATEPART_DO_ID 1079
#define ADDPOLYHED_DO_ID 1080
#define PACKALL_DO_ID 1083
#define CNTWALL_DO_ID 1065
#define CNTWALL_READ_ID 1066
#define CNTWALL_WRITE_ID 1067
#define WRITEQCELEMENT_ID 116
#define INST_ID 105
#define MDSWITCH_ID 109
#define EDIT_ELEM_TAKE_ID 112
#define EDIT_ELEM_XZ_TAKE_ID 115
#define MEASURE_CALC2_ID 303
#define MEASURE_CALC3_ID 304
#define MEASURE_CALC4_ID 305
#define MEASURE_CLEAR_ID 306
#define SETPARAM_READ_ID 204
#define CELLSIZE_APPLY_ID 202
#define CAPTURE_ID 203
#define STRSCHK_ID 207
#define STRSSET_ID 211
#define EXTRA_AUTO_ID 219
#define EXTRA_BKWRITE_ID 220
#define EXTRA_BKREAD_ID 221
#define PHONON_CALC_ID 214
#define DOS_CALC_ID 218
#define NEB_CALC_ID 215
#define CNTSHELL_TEST_ID 216
#define CNTSHELL_CALC_ID 217
#define CALC_ID 208
#define CB_BOND 1001
#define CB_ENSEMBLE 1004
#define CB_CFG_FB 1005
#define CB_SETDAT_FB 1030
#define CB_QCELM_FB 1006
#define CB_CNTWALL_FB 1007
#define CB_COLOR_MODE 1008
#define CB_POTFILE_FB 1009
#define CB_INST 1010
#define CB_PHONON 1011
#define CB_EDITATOM 1012
#define CB_NEB 1013
#define CB_NEB_SHOW 1014
#define CB_ROTATE 1015
#define CB_ROTATE_REV 1016
#define CB_RESET_VIEW 1017
#define CB_CELL_DIM 1018
#define CB_LATTICE 1029
#define CB_SLICE 1019
#define CB_TRIM 1033
#define CB_RIGID_SHIFT 1027
#define CB_RIGID_SHIFT_XY 1028
#define CB_RIGID_SHIFT_RELAX 1031
#define CB_FIRE 1020
#define CB_CG 1023
#define CB_REL 1024
#define CB_ELASC 1021
#define CB_FCHECK 1022
#define CB_MARKED_ATOM 1025
#define CB_CELLFIX_ID 1026
#define CB_SHOW_GROUP 1032
#define CB_QUIT 9999

// 2024/02 CJS 追加
#define STATUS_LX 2020
#define STATUS_LY 2021
#define STATUS_LZ 2022
#define STATUS_DT 2023
// ウィジェットと変数の値の同期を実行する間隔 (秒) です。
#define POLLING_INTERVAL 0.1
// 2024/02 CJS 追加 ここまで

extern float xy_aspect;
extern int   last_x, last_y;
extern float rotationX, rotationY;
extern float rotmat[16];
extern float obj_pos[];
extern float scl;
extern float vscl;
extern float vscl_force;
extern int   main_window;


/** These are the live variables passed into GLUI ***/
// 原子のヒット機能のオンオフを指定する変数です。2024/11 CJS
extern int enable_atom_click;
extern int ortho_state;
extern int draw_bond, draw_bond_pbc, draw_force, draw_load, draw_aux;
extern int bond_display_type, bond_display_color;
extern int relax_algo, relax_accel, relax_accel_interval;
extern float relax_accel_threshold, relax_damper_value;
extern int cell_relax_rep;
extern float cell_relax_tolfor;
extern int lattice_type;
extern char config_atom[3];
extern char config_atom2[3];
extern float alplat1;
extern char cwdname[80];
extern const char* ensemble_list[];
extern const char* kpstring_list[];
extern const char* color_list[];
extern const char* sw_arg_list[];
extern const char* tersoff_arg_list[];
extern int segments, radius;
extern int show_only_elem, show_axis, show_cell;
extern int show_cnt_wall, show_cnt_wallv, show_cnt_ring;
extern int show_cnt_wall_num;
extern int mode_cnt_corrugation, read_cntwall;
extern int show_cnt_all;
extern int show_cnt_layer[6];
extern int cnt_load_algo;
extern float yzplane_punch_d, yzplane_punch_dd, yzplane_punch_ftot;
extern float cnt_pressure, cnt_pressure_ftot, cnt_pressure_gpa, cnt_pressure_gpa2;
extern float outermost_radius_f, cylinder_side_area_f, cylinder_side_area0_f;
extern float cnt_ring_radius, cnt_ring_fmax, cnt_ring_sharpness;
extern float cntshell_dmu, cntshell_dnu, cntshell_eps;
extern int cntshell_n;
extern int ievec, ievec_num;
extern float evec_len, eigval;
extern double mat[3][3];
extern float org_x, org_y, org_z;
extern int size_w, size_h;
extern int b_state;
extern float xx1, yy1, zz1, xx2, yy2, zz2, xx3, yy3, zz3, xx4, yy4, zz4;
//extern GLuint objects;
extern int edit_elem_mode, select_atom[10], select_atom_repidx[10];
extern int measure_mode;
extern int draw_replica;
extern float replica_range;
extern int hoge;
extern float bondlength;
extern int bondthickness;
extern int capture_count;
extern float eps_strschk;
extern float fcheck_disp;
extern int color_mode, color_mode_auto;
extern float color_mode_vmin, color_mode_vmax;
extern bool iprdamper, ivscale, irelax;
extern float prdamper_val1, prdamper_val2, prlimit;
extern int hessian_read, hessian_write, inst_mode;
extern int phonon_rep, phonon_kp, phonon_knum;
extern float dos_gauss_width;
extern int dos_kmesh;
extern int neb_num, neb_node, neb_ite, neb_init_read;
extern float neb_tol_fac;
extern int neb_kconst;
extern float neb_kfactor;
extern int iatom_pick, createconfig_mode, editatom_mode;
extern int iatom_from, iatom_to;
extern float atomrx, atomry, atomrz;
extern int milX1, milX2, milX3;
extern int milY1, milY2, milY3;
extern int milZ1, milZ2, milZ3;
extern float cellrot_x, cellrot_y, cellrot_z;
extern int marked_atom[3], marked_atom_color[3];
extern float rotatomx, rotatomy, rotatomz;
extern int polyhedtype;
extern float polyhedsize;
extern int icnt;


// 2024/02 CJS 追加
// キャプチャ用の変数です。
extern int capture_width;
extern int capture_height;
// 2024/02 CJS ここまで


#endif // _MYHEADER_H



#ifdef CUDA
extern  int* arrayDint;
extern  int* arrayDrepatom;
extern  double* arrayDrx, * arrayDry, * arrayDrz, * arrayDfx, * arrayDfy, * arrayDfz;
extern  double* arrayDepot, * arrayDhmat, arrayHhmat[9];
extern  int* arrayDalistnum, * arrayDalist;
extern  int iarray[( NMAX + 1 ) * ( NBOOK + 1 ) * 4];
#endif

