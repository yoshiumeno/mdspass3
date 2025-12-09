#include "myheader.h"
//#include "pair_function.h"

void PairFunction::Initialize(void)
{
  NumOfParam = 0;
  CutoffRadius = 0.0;
  ParamTable = NULL;
  Func = NULL;
  ID ="";
  
  FlgWarningForNumericalDerivative = true;
}

int PairFunction::GetNumOfParam(string funcform){
    if(strcmp(funcform.c_str(),"const")==0)         {return(1);}
    if(strcmp(funcform.c_str(),"csw2")==0)          {return(4);}
    if(strcmp(funcform.c_str(),"csw2_sc")==0)       {return(5);}
    if(strcmp(funcform.c_str(),"eopp_exp")==0)      {return(6);}
    if(strcmp(funcform.c_str(),"eopp_exp_sc")==0)   {return(7);}
    if(strcmp(funcform.c_str(),"exp_plus")==0)      {return(3);}
    if(strcmp(funcform.c_str(),"exp_plus_sc")==0)   {return(4);}
    if(strcmp(funcform.c_str(),"meopp")==0)         {return(7);}
    if(strcmp(funcform.c_str(),"meopp_sc")==0)      {return(8);}
    if(strcmp(funcform.c_str(),"mishin")==0)        {return(6);}
    if(strcmp(funcform.c_str(),"mishin_sc")==0)     {return(7);}
    if(strcmp(funcform.c_str(),"ms")==0)            {return(3);}
    if(strcmp(funcform.c_str(),"poly_5")==0)        {return(5);}
    if(strcmp(funcform.c_str(),"vashpair")==0)      {return(7);} // Kubo 20140619
    if(strcmp(funcform.c_str(),"vashpair_shift")==0){return(7);} //
    if(strcmp(funcform.c_str(),"vashpair2")==0)     {return(8);} // Kubo 20140718
    if(strcmp(funcform.c_str(),"vashpair2_sc")==0)  {return(9);} // Kubo 20140718
    if(strcmp(funcform.c_str(),"vashtrior")==0)     {return(2);} // Kubo 20140619
    if(strcmp(funcform.c_str(),"vashtriop")==0)     {return(3);} //
    if(strcmp(funcform.c_str(),"swpair")==0)        {return(5);} //
    if(strcmp(funcform.c_str(),"swtrior")==0)       {return(3);} //
    if(strcmp(funcform.c_str(),"swtriop")==0)       {return(2);} //
    return(0);
  }

// PSI Function for Cutoff -------------------------------------
void PairFunction::Psi(
        double rr,      //  IN: Distance (Angst.)
        double h,       //  IN: Smoothing Parameter (Angst.)
        double &val,    // OUT: Value (-)
        double &grad    // OUT: Derivative (1/Angst.)
  ){
    double x  = (rr-CutoffRadius)/h;
    double x2 = x*x;
    double x3 = x2*x;
    double x4 = x2*x2;
    if(rr<CutoffRadius){
      val  = x4/(1.0+x4);
      grad = 4.0*x3/(1.0+x4)/(1.0+x4)/h;
    }else{
      val  = 0.0;
      grad = 0.0;
    }
  }

// Template ----------------------------------------------------
void PairFunction::Func_Template(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
  ){
    val  = pf->GetParam(0)*0.0;
    grad = pf->GetParam(0)*0.0;
  }
void PairFunction::Func_Template_sc(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
  ){
    double val_tmp, grad_tmp;
    double sc_tmp,  scp_tmp;
    double h;

    h = pf->GetParam(0);
    Func_Template(rr,pf,val_tmp,grad_tmp);
    (pf->Psi)(rr,h,sc_tmp,scp_tmp);

    val  = val_tmp*sc_tmp;
    grad = val_tmp*scp_tmp + grad_tmp*sc_tmp;
  }


// Zero --------------------------------------------------------
void PairFunction::Zero(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
  ){
    val  = 0.0;
    grad = 0.0;
  }


// Const -------------------------------------------------------
void PairFunction::Const(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
  ){
    val  = pf->GetParam(0);
    grad = 0.0;
  }

// CSW2 --------------------------------------------------------
void PairFunction::Csw2(
	  double rr,        // IN : Distance (Angst.)
	  PairFunction *pf, // IN : Pointer to PairFunction
	  double &val,      // OUT: Energy (eV) or Others
	  double &grad      // OUT: Force (eV/Angst.) or Others
	  )
{
  double p = pow(rr, pf->GetParam(3));
  val  = (1.0+pf->GetParam(0)*cos(pf->GetParam(1)*rr+pf->GetParam(2)))/p;
  grad = (-pf->GetParam(0)*pf->GetParam(1)*sin(pf->GetParam(1)*rr+pf->GetParam(2)))/p
    - pf->GetParam(3)*val/rr;
}
void PairFunction::Csw2_sc(
	     double rr,        // IN : Distance (Angst.)
	     PairFunction *pf, // IN : Pointer to PairFunction
	     double &val,      // OUT: Energy (eV) or Others
	     double &grad      // OUT: Force (eV/Angst.) or Others
	     )
{
  double val_tmp, grad_tmp;
  double sc_tmp,  scp_tmp;
  double h;
  
  h = pf->GetParam(pf->GetNumOfParam()-1);
  Csw2(rr,pf,val_tmp,grad_tmp);
  (pf->Psi)(rr,h,sc_tmp,scp_tmp);
  
  val  = val_tmp*sc_tmp;
  grad = val_tmp*scp_tmp + grad_tmp*sc_tmp;
}

// Eopp_Exp ----------------------------------------------------
void PairFunction::Eopp_Exp(
	      double rr,        // IN : Distance (Angst.)
	      PairFunction *pf, // IN : Pointer to PairFunction
	      double &val,      // OUT: Energy (eV) or Others
	      double &grad      // OUT: Force (eV/Angst.) or Others
	      )
{
  double tmp[5];
  tmp[0] = pf->GetParam(0)*exp(-pf->GetParam(1)*rr);
  tmp[1] = pf->GetParam(2)/pow(rr,pf->GetParam(3));
  tmp[2] = pf->GetParam(4)*rr + pf->GetParam(5);
  tmp[3] = cos(tmp[2]);
  tmp[4] = sin(tmp[2]);
  
  val  = tmp[0] + tmp[1]*tmp[3];
  grad = - pf->GetParam(1)*tmp[0]
    - pf->GetParam(3)*tmp[1]*tmp[3]/rr
    - pf->GetParam(4)*tmp[1]*tmp[4];
}
void PairFunction::Eopp_Exp_sc(
		 double rr,        // IN : Distance (Angst.)
		 PairFunction *pf, // IN : Pointer to PairFunction
		 double &val,      // OUT: Energy (eV) or Others
		 double &grad      // OUT: Force (eV/Angst.) or Others
		 )
{
  double val_tmp, grad_tmp;
  double sc_tmp,  scp_tmp;
  double h;
  
  h = pf->GetParam(pf->GetNumOfParam()-1);
  Eopp_Exp(rr,pf,val_tmp,grad_tmp);
  (pf->Psi)(rr,h,sc_tmp,scp_tmp);
  
  val  = val_tmp*sc_tmp;
  grad = val_tmp*scp_tmp + grad_tmp*sc_tmp;
}

// Exp_Plus ----------------------------------------------------
void PairFunction::Exp_Plus(
	      double rr,        // IN : Distance (Angst.)
	      PairFunction *pf, // IN : Pointer to PairFunction
	      double &val,      // OUT: Energy (eV) or Others
	      double &grad      // OUT: Force (eV/Angst.) or Others
	      )
{
  double tmp = pf->GetParam(0)*exp(-pf->GetParam(1)*rr);
  val  = tmp + pf->GetParam(2);
  grad = -pf->GetParam(1) * tmp;
}
void PairFunction::Exp_Plus_sc(
		 double rr,        // IN : Distance (Angst.)
		 PairFunction *pf, // IN : Pointer to PairFunction
		 double &val,      // OUT: Energy (eV) or Others
		 double &grad      // OUT: Force (eV/Angst.) or Others
		 )
{
  double val_tmp, grad_tmp;
  double sc_tmp,  scp_tmp;
  double h;
  
  h = pf->GetParam(pf->GetNumOfParam()-1);
  Exp_Plus(rr,pf,val_tmp,grad_tmp);
  (pf->Psi)(rr,h,sc_tmp,scp_tmp);
  
  val  = val_tmp*sc_tmp;
  grad = val_tmp*scp_tmp + grad_tmp*sc_tmp;
}

// Meopp -------------------------------------------------------
void PairFunction::Meopp(
	   double rr,        // IN : Distance (Angst.)
	   PairFunction *pf, // IN : Pointer to PairFunction
	   double &val,      // OUT: Energy (eV) or Others
	   double &grad      // OUT: Force (eV/Angst.) or Others
	   )
{
  double tmp[6];
  
  tmp[0] = rr - pf->GetParam(6);
  tmp[1] = pf->GetParam(0)/pow(tmp[0],pf->GetParam(1));
  tmp[2] = pf->GetParam(2)/pow(rr,pf->GetParam(3));
  tmp[3] = pf->GetParam(4)*rr + pf->GetParam(5);
  tmp[4] = cos(tmp[3]);
  tmp[5] = sin(tmp[3]);
  
  val  = tmp[1] + tmp[2]*tmp[4];
  //    grad = - pf->GetParam(2)*tmp[1]/tmp[0]
  grad = - pf->GetParam(1)*tmp[1]/tmp[0] // Kubo 20140318
    - pf->GetParam(3)*tmp[2]*tmp[4]/rr
    - pf->GetParam(4)*tmp[2]*tmp[5];
  
}
void PairFunction::Meopp_sc(
	      double rr,        // IN : Distance (Angst.)
	      PairFunction *pf, // IN : Pointer to PairFunction
	      double &val,      // OUT: Energy (eV) or Others
	      double &grad      // OUT: Force (eV/Angst.) or Others
	      )
{
  double val_tmp, grad_tmp;
  double sc_tmp,  scp_tmp;
  double h;
  
  h = pf->GetParam(pf->GetNumOfParam()-1);
  Meopp(rr,pf,val_tmp,grad_tmp);
  (pf->Psi)(rr,h,sc_tmp,scp_tmp);
  
  val  = val_tmp*sc_tmp;
  grad = val_tmp*scp_tmp + grad_tmp*sc_tmp;
}

// Mishin ------------------------------------------------------
void PairFunction::Mishin(
	    double rr,        // IN : Distance (Angst.)
	    PairFunction *pf, // IN : Pointer to PairFunction
	    double &val,      // OUT: Energy (eV) or Others
	    double &grad      // OUT: Force (eV/Angst.) or Others
	    )
{
  double z  = rr - pf->GetParam(3);
  //    double e  = exp(-pf->GetParam(5)*z);
  double e  = exp(-pf->GetParam(5)*rr); // Defition in POTFIT
  double p  = pow(z,pf->GetParam(4));
  double pp = pf->GetParam(4)*p/z;
  double ep = -pf->GetParam(5)*e;
  val  = pf->GetParam(0)*  p*e*(1.0+pf->GetParam(1)*e) + pf->GetParam(2);
  grad = pf->GetParam(0)*(pp*e*(1.0+pf->GetParam(1)*e)
			  + p*ep*(1.0+pf->GetParam(1)* e)
			  + p* e*(    pf->GetParam(1)*ep));
}
void PairFunction::Mishin_sc(
	       double rr,        // IN : Distance (Angst.)
	       PairFunction *pf, // IN : Pointer to PairFunction
	       double &val,      // OUT: Energy (eV) or Others
	       double &grad      // OUT: Force (eV/Angst.) or Others
	       )
{
  double val_tmp, grad_tmp;
  double sc_tmp,  scp_tmp;
  double h;
  
  h = pf->GetParam(pf->GetNumOfParam()-1);
  Mishin(rr,pf,val_tmp,grad_tmp);
  (pf->Psi)(rr,h,sc_tmp,scp_tmp);
  
  val  = val_tmp*sc_tmp;
  grad = val_tmp*scp_tmp + grad_tmp*sc_tmp;
}

// MS (Morse Stretch) ------------------------------------------
void PairFunction::MS(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
	)
{
  double de,a,r0;
  de = pf->GetParam(0);
  a  = pf->GetParam(1);
  r0 = pf->GetParam(2);
  
  double e1 = exp((1.0-rr/r0)*a/2.0);
  double e2 = e1*e1;
  
  val  = de*(e2-2.0*e1);
  grad = -de*a/r0*(e2-e1);
}

// Poly_5 ------------------------------------------------------
void PairFunction::Poly_5(
	    double rr,        // IN : Distance (Angst.)
	    PairFunction *pf, // IN : Pointer to PairFunction
	    double &val,      // OUT: Energy (eV) or Others
	    double &grad      // OUT: Force (eV/Angst.) or Others
	    )
{
  double  r1 = rr - 1.0;
  double dr1 = r1*r1;
  val =     pf->GetParam(0)
    + 0.5*pf->GetParam(1)*dr1
    +     pf->GetParam(2)*r1*dr1
    +     pf->GetParam(3)*dr1*dr1
    +     pf->GetParam(4)*r1*dr1*dr1;
  grad =    pf->GetParam(1)*r1
    + 3.0*pf->GetParam(2)*dr1
    + 4.0*pf->GetParam(3)*r1*dr1
    + 5.0*pf->GetParam(4)*dr1*dr1;
}

// Stillinger-Weber
//void PairFunction::PairFunction::SwPair( // under construction!!
void PairFunction::SwPair( // under construction!!
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
  ){
    // param[0]: epsilon
    // param[1]: A
    // param[2]: B
    // param[3]: a
    // param[4]: sigma
    double term[4], decay[2], rij,r2,r4;
    double br4, rija;
    rij  = rr / pf->GetParam(4);
    rija = rij - pf->GetParam(3); //r_ij-a
    if (rija >= -1e-10) {
      val = 0; grad = 0;
    } else {
      r2   = rij*rij;
      r4   = r2*r2;
      br4  = pf->GetParam(2) / r4; //B*r_ij^-4
      term[0] = pf->GetParam(0)*pf->GetParam(1);
      term[1] = br4 - 1.0;
      term[2] = exp( 1.0 / rija );
      val  = term[0] * term[1] * term[2];
      grad = - term[0] / pf->GetParam(4) * (4. * br4 / rij + (br4 - 1.0)/rija/rija) * term[2];
    }
  }
//void PairFunction::PairFunction::SwTrioR(
void PairFunction::SwTrioR(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
  ){
    // param[0]: gamma
    // param[1]: a
    // param[2]: sigma
    double rij, dr, x;
    rij = rr / pf->GetParam(2);
    dr  = rij - pf->GetParam(1);
    if(dr>-1e-10){val = 0; grad = 0; return;}
    x  = pf->GetParam(0)/dr;
    val  = exp(x);
    grad = -val*x/dr/pf->GetParam(2);
  }
//void PairFunction::PairFunction::SwTrioP(
void PairFunction::SwTrioP(
        double cosjik,    // IN : Angle j-i-k
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
  ){
    // param[0]: epsilon
    // param[1]: lambda
    double dcos, dcos2;
    dcos  = cosjik + 1.0 / 3.0;
    dcos2 = dcos*dcos;
    if(fabs(dcos)<1e-10){val = 0; grad = 0; return;}

    val  = pf->GetParam(0) * pf->GetParam(1) * dcos2;
    grad = 2.0 * pf->GetParam(0) * pf->GetParam(1) * dcos;
  }
  // End of Stillinger-Weber


// VashPair ----------------------------------------------------
//void PairFunction::PairFunction::VashPair(
void PairFunction::VashPair(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
  ){
    double term[4], decay[2], r2,r4,r6;
    r2 = rr*rr;
    r4 = r2*r2;
    r6 = r2*r4;
    decay[0] = exp(-rr/pf->GetParam(5));
    decay[1] = exp(-rr/pf->GetParam(6));
    term[0] = pf->GetParam(0)/pow(rr,pf->GetParam(4));
    term[1] = pf->GetParam(1)*decay[0]/rr;
    term[2] = pf->GetParam(2)*decay[1]/r4/2.0;
    term[3] = pf->GetParam(3)/r6;
    val  = + term[0]
           + term[1]
           - term[2]
           - term[3];
    grad = - term[0]*pf->GetParam(4)/rr
           - term[1]*(1.0/rr+1.0/pf->GetParam(5))
           + term[2]*(4.0/rr+1.0/pf->GetParam(6))
           + term[3]*6/rr;
  }
//void PairFunction::PairFunction::VashPair_shift(
void PairFunction::VashPair_shift(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
  ){
//// This routine can be improved by storing val_rc and grad_rc.
    double val_rr, grad_rr;
    double val_rc, grad_rc;
    double rc = pf->GetCutoffRadius();

    if(rr>=rc){val = 0.0; grad = 0.0; return;}

    VashPair(rr,pf,val_rr,grad_rr);
    VashPair(rc,pf,val_rc,grad_rc);

    val  = val_rr  - val_rc - (rr-rc)*grad_rc;
    grad = grad_rr - grad_rc;
  }

// VashPair2 (Another notation of VashPair) --------------------
//void PairFunction::PairFunction::VashPair2(
void PairFunction::VashPair2(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
  ){
    double term[4], decay[2], r2,r4,r6;
    double pot[7];
    r2 = rr*rr;
    r4 = r2*r2;
    r6 = r2*r4;
    pot[0] =       pf->GetParam(0); //---------------- H_ij
    pot[1] = 14.40*pf->GetParam(1)*pf->GetParam(2); // const*Z_i*Z_j
    pot[2] = 14.40*pf->GetParam(3); //---------------- const*D_ij
    pot[3] =       pf->GetParam(4); //---------------- W_ij
    pot[4] =       pf->GetParam(5); //---------------- eta_ij
    pot[5] =       pf->GetParam(6); //---------------- lambda
    pot[6] =       pf->GetParam(7); //---------------- xi
    decay[0] = exp(-rr/pot[5]);
    decay[1] = exp(-rr/pot[6]);
    term[0] = pot[0]/pow(rr,pot[4]);
    term[1] = pot[1]*decay[0]/rr;
    term[2] = pot[2]*decay[1]/r4/2.0;
    term[3] = pot[3]/r6;
    val  = + term[0]
           + term[1]
           - term[2]
           - term[3];
    grad = - term[0]*pot[4]/rr
           - term[1]*(1.0/rr+1.0/pot[5])
           + term[2]*(4.0/rr+1.0/pot[6])
           + term[3]*6/rr;
  }
//void PairFunction::PairFunction::VashPair2_sc(
void PairFunction::VashPair2_sc(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
  ){
    double val_tmp, grad_tmp;
    double sc_tmp,  scp_tmp;
    double h;

    h = pf->GetParam(pf->GetNumOfParam()-1);
    VashPair2(rr,pf,val_tmp,grad_tmp);
    (pf->Psi)(rr,h,sc_tmp,scp_tmp);

    val  = val_tmp*sc_tmp;
    grad = val_tmp*scp_tmp + grad_tmp*sc_tmp;
  }


// VashTrioR (Vashishta Angular Part) --------------------------
//void PairFunction::PairFunction::VashTrioR(
void PairFunction::VashTrioR(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
  ){
    double dr, x;
    dr = rr - pf->GetParam(1);
    if(dr>-1e-10){val = 0; grad = 0; return;}
    x  = pf->GetParam(0)/dr;

    val  = exp(x);
    grad = -val*x/dr;
  }

// VashTrioP (Vashishta Angular Part) --------------------------
//void PairFunction::PairFunction::VashTrioP(
void PairFunction::VashTrioP(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
  ){
    double dr, dr2;
    dr = rr - pf->GetParam(2);
    dr2 = dr*dr;
    if(fabs(dr)<1e-10){val = 0; grad = 0; return;}

    val  = pf->GetParam(0)*dr2/(1.0+pf->GetParam(1)*dr2);
    grad = 2.0*val*(1.0-pf->GetParam(1)/pf->GetParam(0)*val)/dr;
  }
