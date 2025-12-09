#ifndef _PAIR_FUNCTION2_H
#define _PAIR_FUNCTION2_H

enum PairType
{
	BUCKINGHAMPOTTYPE = 1, // Buckingham potential type
	MORSEPOTTYPE = 2 // Morse potential type
};

enum PotentialMode
{
	NORMALPOTMODE = 1, // Normal potential type
	COMBINEDPOTMODE = 2 // Combined potential type (initial values of forces, energies and stresses are not zeroed)
};


/* Base class for all pair interactions,
so far it is supposed that pairwise interactions between all atom types must be specified in the potential file,
even if they are zeroes
 */
class Pair
{
	public:
		/* class constructor */
		Pair(int, double, PotentialMode);
		
		/* class destructor */
		~Pair();
		
		/* read potential parameters from a POTFIT-format file*/
		int ReadParams();

		/* calculate forces, energies and stresses */
		void Calculate();

		/* reset atomic forces, energies and stresses */
		void Reset();

		/* virtual function for energy calculation of pair interaction */
		virtual double pairEnergy(double, int) = 0;

		/* virtual function for force calculation of pair interaction */
		virtual double pairForce(double, int) = 0;	

		// name of the file with potential parameters
		char fname[60];

		// is successfully initialized
		bool initialized; 

	protected:

		//number of parameters of each pair interation
		int nparams;

		// total number of atomic species
		int ntype;
		int* types;

		//total number of atomic pairs to describe
		int npairs;

		//array with parameters of pair interations (first index - pair number, second - parameter number)
		double** params;
		int** pairnums; // two dimensional array (if 0 and 1 are indices of Zr and Y respictively in types array, then pairnums[0][1] gives index of Zr-Y interaction in params array)
		int* atomtypes; // list with indices of atomic species of all atoms (if ith atom is Zr and Zr is first element in types, then atomtypes[i] = 0)
	
		// type of the potential
		PairType type;

		// is successfully initialized
		// bool initialized; // moved to public

		// name of the file with potential parameters
		// char fname[60]; // moved to public

		// double cut-off radius for pair interactions and its square
		double rcut;
		double r2;

		// mode of calculation (0 - normal, 1 - combined). For combined mode initial values of forces, energies and stresses are not zeroed
		int mode;

};

/* Double morse pairwise potential */
class Morse : public Pair 
{
	public:
	Morse(int, double, PotentialMode);

	protected:
	virtual double pairEnergy(double, int);
	virtual double pairForce(double, int);

	/* help functions */
	double MorseExp(double, double, double);
	double MorseExpDeriv(double, double, double);	
};

/* Buckongham pairwise potential */
class Buckingham : public Pair 
{
	public:
	Buckingham(int, double, PotentialMode);

	protected:
	virtual double pairEnergy(double, int);
	virtual double pairForce(double, int);	
};

#endif // _PAIR_FUNCTION2_H
