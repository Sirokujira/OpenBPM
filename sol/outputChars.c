/*
outputChars.c

calculate and output to oth.log
*/

#include "obpm.h"
#include "user_define_complex.h"
#include "obpm_prototype.h"


void outputChars(FILE *fp)
{
	// setup far field
	alloc_farfield();
	setup_farfield();

	// input imepedanece
	if (NFeed && NFreq1) {
		calcZin();
		outputZin(fp);
		calcPin();  // for post
	}

	// S-parameters
	if (NPoint && NFreq1) {
		calcSpara();
		outputSpara(fp);
	}

	// coupling
	if (NFeed && NPoint && NFreq1) {
		outputCoupling(fp);
	}

	// cross section
	if (IPlanewave && NFreq2) {
		outputCross(fp);
	}
}
