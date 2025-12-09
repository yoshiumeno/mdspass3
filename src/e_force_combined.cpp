#include <iostream>
#include<string>
#include<stdlib.h>
#include<fstream>
#include "myheader.h"

void e_force_dipole(PotentialMode mode);
void e_force_mishin(PotentialMode mode);

void e_force_combined()
{	
	if (combined.initialize) {
		std::cout << "Initialize of the Combined potential.." << std::endl;
		combined.initialize = false;
	}

	/*! IMPORTANT: ONLY THIS ORDER OF FUNCTION CALLS HAS BEEN TESTED (theoretically any order would work properly) */
	pairPot->Reset();

	/* calculate dipole contribution (!!!combined mode must be used!!!) */
	e_force_dipole(COMBINEDPOTMODE);

	/* calculate EAM contribution (!!!combined mode must be used!!!) */
	e_force_mishin(COMBINEDPOTMODE);

	/* calculate Buckingham interaction */
	pairPot->Calculate();

	/* need to recalculate total potential energy */
	atom.epotsum = 0;
	for (int i=1; i<=atom.natom; i++) { atom.epotsum += atom.epot[i]; }
}	
