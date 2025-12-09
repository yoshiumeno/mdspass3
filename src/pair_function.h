#ifndef _PAIR_FUNCTION_H
#define _PAIR_FUNCTION_H

#include<iostream>
#include<iomanip>
#include<fstream>

using namespace std;

class PairFunction{ //========================================//


 protected: //==================================================

//////// Member Variables //////////////////////////////////////
  int     NumOfParam;   // Number of Parameters in Function.
  double  CutoffRadius; // Cutoff Radius for Interaction. Nonpositive Value kills Cutoff.
  double *ParamTable;   // Potential Paramater Table
  void  (*Func)(double, PairFunction*, double&, double&);
  string  ID;           // Identifier for Error Message (Unimplemented...).

  bool FlgWarningForNumericalDerivative;

 private: //====================================================

  // Initializer ---------------
  // Initialize Members by naught.
  void Initialize(void);

  // Set NumOfParam, Func, ParamTable.
  void Initialize(string funcform, double *param){
    Initialize();

    Func = GetFunc(funcform);
    NumOfParam = GetNumOfParam(funcform);
    if(NumOfParam>0){
      ParamTable = new double[NumOfParam];
      for(int i=0;i<NumOfParam;i++){ParamTable[i] = param[i];}
    }
  }

  // Set Above + CutoffRadius.
  void Initialize(string funcform, double *param, double rc){
    Initialize(funcform,param);
    CutoffRadius = rc;
  }

  // Set Above + ID.
  void Initialize(string funcform, double *param,double rc,string tag){
    Initialize(funcform,param,rc);
    ID = tag;
  }

  // Numerical Derivative
  double NumericalDerivative(double rr,double dr){
    if(FlgWarningForNumericalDerivative){
      cout<<"Numerical Derivative Applied";
      if(strcmp(ID.c_str(),"")==0) {cout << "." << endl;}
      else{cout<<" for "<< ID.c_str() << "."<<endl;}
      FlgWarningForNumericalDerivative = false;
    }

    double fx0,fx1,fx2;
    double dfdx0,dfdx1,dfdx2;

    (*Func)(rr,   this,fx0,dfdx0);
    (*Func)(rr+dr,this,fx1,dfdx1);
    (*Func)(rr-dr,this,fx2,dfdx2);

    dfdx0 = (fx1-fx2)/(2.0*dr);
    return(dfdx0);
  }

 public: //=====================================================
  // Constructor ---------------
  PairFunction(void){
    Initialize();
  }

  PairFunction(string funcform, double *param){
    Initialize(funcform,param);
  }

  PairFunction(string funcform, double *param, double rc){
    Initialize(funcform,param,rc);
  }

  PairFunction(string funcform, double *param, double rc,string s){
    Initialize(funcform,param,rc,s);
  }


  // Destructor ----------------
 ~PairFunction(){
    if(ParamTable!=NULL){
      delete [] ParamTable;
      ParamTable = NULL;
    }
  }

//////// Common Functions //////////////////////////////////////
  // Calculate Function
  void Calc(double rr, double &val, double &grad){
    if((CutoffRadius>0.0)   // Cutoff Working?
     &&(rr>CutoffRadius)){  // Out of Cutoff?
      val  = 0.0;
      grad = 0.0;
    }else{
      (*Func)(rr,this,val,grad);
//      grad = NumericalDerivative(rr,0.001); // For Debug;
    }
  }

  // Plotter for Debugging
  void Plot(
    double min,  // Minmum  Value to Plot
    double max,  // Maximum Value to Plot
    int nmesh,   // Number of Mesh Data
    string fname // File Name to Output
  ){
    if(Func==NULL){return;} // Kubo 20140619

    double x,fx,dfdx,dfdx_nmr;
    double dx = fabs(max-min)/nmesh;
    ofstream ofs(fname.c_str());

    if(min>max){min = max;}

    ofs << "# r, f(r), df/dr(r)_analytic, df/dr(r)_numerical" << endl;

    bool flgtmp = FlgWarningForNumericalDerivative;
    FlgWarningForNumericalDerivative = false;

    for(int i=0;i<=nmesh;i++){
      x = min + dx*i;

      Calc(x,fx,dfdx);
      dfdx_nmr = NumericalDerivative(x,0.001);

      ofs << setiosflags(ios::scientific)
          << setw(16) << setprecision(6)  << x;
      ofs << setiosflags(ios::scientific)
          << setw(22) << setprecision(10) << fx;
      ofs << setiosflags(ios::scientific)
          << setw(22) << setprecision(10) << dfdx;
      ofs << setiosflags(ios::scientific)
          << setw(22) << setprecision(10) << dfdx_nmr;
      ofs << endl;
    }

    FlgWarningForNumericalDerivative = flgtmp;
  }

//////// Reference to Members //////////////////////////////////
  int GetNumOfParam(void){
    return(NumOfParam);
  }

  double GetCutoffRadius(void){ //// Kubo 20140619
    return(CutoffRadius);
  }

  double GetParam(int id){
    if(id<NumOfParam){
      return(ParamTable[id]);
    }else{
      cout << "WARNING: Mismatch in Number of Parameters." << endl;
      return(0.0);
    }
  }





////////////////////////////////////////////////////////////////X
////                                                        ////X
////          EDIT BELOW TO ADD NEW FUNCTION FORMS          ////X
////                                                        ////X
////////////////////////////////////////////////////////////////X

  int GetNumOfParam(string funcform);

  // protected: //==================================================
  void (*GetFunc(string funcform))(double, PairFunction*, double&, double&){
    if(strcmp(funcform.c_str(),"const")==0)         {return(Const);}
    if(strcmp(funcform.c_str(),"csw2")==0)          {return(Csw2);}
    if(strcmp(funcform.c_str(),"csw2_sc")==0)       {return(Csw2_sc);}
    if(strcmp(funcform.c_str(),"eopp_exp")==0)      {return(Eopp_Exp);}
    if(strcmp(funcform.c_str(),"eopp_exp_sc")==0)   {return(Eopp_Exp_sc);}
    if(strcmp(funcform.c_str(),"exp_plus")==0)      {return(Exp_Plus);}
    if(strcmp(funcform.c_str(),"exp_plus_sc")==0)   {return(Exp_Plus_sc);}
    if(strcmp(funcform.c_str(),"meopp")==0)         {return(Meopp);}
    if(strcmp(funcform.c_str(),"meopp_sc")==0)      {return(Meopp_sc);}
    if(strcmp(funcform.c_str(),"mishin")==0)        {return(Mishin);}
    if(strcmp(funcform.c_str(),"mishin_sc")==0)     {return(Mishin_sc);}
    if(strcmp(funcform.c_str(),"ms")==0)            {return(MS);}
    if(strcmp(funcform.c_str(),"poly_5")==0)        {return(Poly_5);}
    if(strcmp(funcform.c_str(),"vashpair")==0)      {return(VashPair);}       // Kubo 20140619
    if(strcmp(funcform.c_str(),"vashpair_shift")==0){return(VashPair_shift);} //
    if(strcmp(funcform.c_str(),"vashpair2")==0)     {return(VashPair2);}      // Kubo 20140718
    if(strcmp(funcform.c_str(),"vashpair2_sc")==0)  {return(VashPair2_sc);}   //
    if(strcmp(funcform.c_str(),"vashtrior")==0)     {return(VashTrioR);}      // Kubo 20140619
    if(strcmp(funcform.c_str(),"vashtriop")==0)     {return(VashTrioP);}      //
    if(strcmp(funcform.c_str(),"swpair")==0)        {return(SwPair);}
    if(strcmp(funcform.c_str(),"swtrior")==0)       {return(SwTrioR);}
    if(strcmp(funcform.c_str(),"swtriop")==0)       {return(SwTrioP);}
    cout << "Unknown Function ("<< funcform.c_str() <<") Type Found." << endl;
    cout << "Zero Function, f(x)=0, Applied for First Aid." << endl;
    return(Zero);
  }

// PSI Function for Cutoff -------------------------------------
  void Psi(
        double rr,      //  IN: Distance (Angst.)
        double h,       //  IN: Smoothing Parameter (Angst.)
        double &val,    // OUT: Value (-)
        double &grad    // OUT: Derivative (1/Angst.)
	   );

// Template ----------------------------------------------------
  static void Func_Template(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
			    );
  static void Func_Template_sc(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
			       );
// Zero --------------------------------------------------------
  static void Zero(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
		   );
// Const -------------------------------------------------------
  static void Const(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
		    );
// CSW2 --------------------------------------------------------
  static void Csw2(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
		   );
  static void Csw2_sc(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
		      );

// Eopp_Exp ----------------------------------------------------
  static void Eopp_Exp(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
		       );
  static void Eopp_Exp_sc(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
			  );

// Exp_Plus ----------------------------------------------------
  static void Exp_Plus(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
		       );
  static void Exp_Plus_sc(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
			  );

// Meopp -------------------------------------------------------
  static void Meopp(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
		    );
  static void Meopp_sc(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
		       );

// Mishin ------------------------------------------------------
  static void Mishin(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
		     );
  static void Mishin_sc(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
			);

// MS (Morse Stretch) ------------------------------------------
  static void MS(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
		 );

// Poly_5 ------------------------------------------------------
  static void Poly_5(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
		     );

  // Stillinger-Weber
  static void SwPair( // under construction!!
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
		      );
  static void SwTrioR(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
		      );
  static void SwTrioP(
        double cosjik,    // IN : Angle j-i-k
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
		      );
  // End of Stillinger-Weber

// VashPair ----------------------------------------------------
  static void VashPair(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
		       );
  static void VashPair_shift(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
			     );

// VashPair2 (Another notation of VashPair) --------------------
  static void VashPair2(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
			);
  static void VashPair2_sc(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
			   );

// VashTrioR (Vashishta Angular Part) --------------------------
  static void VashTrioR(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
			);

// VashTrioP (Vashishta Angular Part) --------------------------
  static void VashTrioP(
        double rr,        // IN : Distance (Angst.)
        PairFunction *pf, // IN : Pointer to PairFunction
        double &val,      // OUT: Energy (eV) or Others
        double &grad      // OUT: Force (eV/Angst.) or Others
			);

}; // End of Class ===========================================//


#endif // _PAIR_FUNCTION_H
